subroutine HFB(dotest,maxSCF,thresh,max_diis,level_shift,nNuc,ZNuc,rNuc,ENuc, & 
               nBas,nOrb,nO,S,T,V,Hc,ERI,dipole_int,X,EHFB,eHF,c,P,Panom,F,Delta)

! Perform Hartree-Fock Bogoliubov calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh
  double precision,intent(in)   :: level_shift

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO
  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas) 
  double precision,intent(in)   :: X(nBas,nOrb)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  integer                       :: maxSCF_chem_pot
  integer                       :: iorb
  integer                       :: nSCF
  integer                       :: nSCF_chem_pot
  integer                       :: nOrb2
  integer                       :: nBas2
  integer                       :: nBas_Sq
  integer                       :: n_diis
  double precision              :: ET
  double precision              :: EV
  double precision              :: EJ
  double precision              :: EK
  double precision              :: EL
  double precision              :: chem_pot
  double precision              :: dipole(ncart)

  double precision              :: Conv
  double precision              :: rcond
  double precision              :: delta_chem_pot
  double precision              :: trace_1rdm
  double precision              :: thrs_N
  double precision,external     :: trace_matrix
  double precision,allocatable  :: eHFB_(:)
  double precision,allocatable  :: err(:,:)
  double precision,allocatable  :: err_diis(:,:)
  double precision,allocatable  :: F_diis(:,:)
  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: K(:,:)
  double precision,allocatable  :: cp(:,:)
  double precision,allocatable  :: H_hfb(:,:)
  double precision,allocatable  :: R(:,:)

! Output variables

  double precision,intent(out)  :: EHFB
  double precision,intent(inout):: eHF(nOrb)
  double precision,intent(inout):: c(nBas,nOrb)
  double precision,intent(out)  :: P(nBas,nBas)
  double precision,intent(out)  :: Panom(nBas,nBas)
  double precision,intent(out)  :: F(nBas,nBas)
  double precision,intent(out)  :: Delta(nBas,nBas)

! Hello world

  write(*,*)
  write(*,*)'*****************************'
  write(*,*)'* HF Bogoliubov Calculation *'
  write(*,*)'*****************************'
  write(*,*)

! Useful quantities

  nBas_Sq = nBas*nBas
  nBas2 = nBas+nBas
  nOrb2 = nOrb+nOrb

! Memory allocation

  allocate(J(nBas,nBas))
  allocate(K(nBas,nBas))

  allocate(err(nBas,nBas))

  allocate(cp(nOrb2,nOrb2))
  allocate(H_hfb(nOrb2,nOrb2))
  allocate(R(nOrb2,nOrb2))
  allocate(eHFB_(nOrb2))

  allocate(err_diis(nBas_Sq,max_diis))
  allocate(F_diis(nBas_Sq,max_diis))

! Guess coefficients chem. pot

  chem_pot = 0.5d0*(eHF(nO)+eHF(nO+1))

! Initialization

  thrs_N          = 1d-6
  maxSCF_chem_pot = 1000
  n_diis          = 0
  F_diis(:,:)     = 0d0
  err_diis(:,:)   = 0d0
  rcond           = 0d0

  P(:,:)     = 2d0 * matmul(c(:,1:nO), transpose(c(:,1:nO)))
  Panom(:,:) = 0d0 ! Do sth TODO

  Conv = 1d0
  nSCF = 0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------


  do while(Conv > thresh .and. nSCF < maxSCF)

    ! Increment 

    nSCF = nSCF + 1

    write(*,*)
    write(*,*)'-------------------------------------'
    write(*,'(1X,A1,1X,A15,1X,A1,1X,A15,1X,A1)') &
            '|','Tr[1D]','|','Chem. Pot.','|'
    write(*,*)'-------------------------------------'

    ! Loop to adjust chem_pot
    delta_chem_pot=1d0
    nSCF_chem_pot=0
    do

     ! Build Fock matrix
     
     call Hartree_matrix_AO_basis(nBas,P,ERI,J)
     call exchange_matrix_AO_basis(nBas,P,ERI,K)
     call anomalous_matrix_AO_basis(nBas,Panom,ERI,Delta)
     
     F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:) - chem_pot*S(:,:)
     
     ! Diagonalize H_HFB matrix
     
     H_hfb(:,:) = 0d0
     H_hfb(1:nOrb      ,1:nOrb      ) = matmul(transpose(X),matmul(F,X))
     H_hfb(nOrb+1:nOrb2,nOrb+1:nOrb2) = -H_hfb(1:nOrb,1:nOrb)
     H_hfb(1:nOrb      ,nOrb+1:nOrb2) = matmul(transpose(X),matmul(Delta,X))
     H_hfb(nOrb+1:nOrb2,1:nOrb      ) = transpose(H_hfb(1:nOrb,nOrb+1:nOrb2))
     
     cp(:,:) = H_hfb(:,:)
     call diagonalize_matrix(nOrb2,cp,eHFB_)
     
     ! Build R and extract P and Panom
       
     trace_1rdm = 0d0 
     R(:,:)     = 0d0
     do iorb=1,nOrb
      R(:,:) = R(:,:) + matmul(cp(:,iorb:iorb),transpose(cp(:,iorb:iorb))) 
     enddo
     do iorb=1,nOrb
      trace_1rdm = trace_1rdm + R(iorb,iorb) 
     enddo
     write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1)') &
      '|',trace_1rdm,'|',chem_pot,'|'

     nSCF_chem_pot = nSCF_chem_pot + 1
     if( abs(trace_1rdm-nO) < thrs_N .or. nSCF_chem_pot > maxSCF_chem_pot ) exit

     
     
    enddo

    ! Extract P and Panom from R
    
    P(:,:)     = 0d0
    Panom(:,:) = 0d0
    P(:,:)     = 2d0*matmul(X,matmul(R(1:nOrb,1:nOrb),transpose(X)))
    Panom(:,:) = matmul(X,matmul(R(1:nOrb,nOrb+1:nOrb2),transpose(X)))


    ! Kinetic energy

    ET = trace_matrix(nBas,matmul(P,T))

    ! Potential energy

    EV = trace_matrix(nBas,matmul(P,V))

    ! Hartree energy

    EJ = 0.5d0*trace_matrix(nBas,matmul(P,J))

    ! Exchange energy

    EK = 0.25d0*trace_matrix(nBas,matmul(P,K))

    ! Anomalous energy

    EL = -0.25d0*trace_matrix(nBas,matmul(Panom,Delta))

    ! Total energy

    EHFB = ET + EV + EJ + EK !+ EL

    ! Dump results
    write(*,*)'-----------------------------------------------------------------------------------------------'
    write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A16,1X,A1,1X,A16,1X,A1A16,1X,A1,1X,A10,2X,A1,1X)') &
            '|','#','|','E(HFB)','|','EJ(HFB)','|','EK(HFB)','|','EL(HFB)','|','Conv','|'
    write(*,*)'-----------------------------------------------------------------------------------------------'

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1XF16.10,1X,A1,1X,E10.2,1X,A1,1X)') &
      '|',nSCF,'|',EHFB + ENuc,'|',EJ,'|',EK,'|',EL,'|',Conv,'|'

    write(*,*)'-----------------------------------------------------------------------------------------------'
    write(*,*)
    ! Check convergence TODO HFB does not fulfill [F,1D]=0 ... 

    err = matmul(F,matmul(P,S)) - matmul(matmul(S,P),F)
    if(nSCF > 1) Conv = maxval(abs(err))

  end do
!------------------------------------------------------------------------
! End of SCF loop
!------------------------------------------------------------------------

! Did it actually converge?

  if(nSCF == maxSCF) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'                 Convergence failed                 '
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)

    deallocate(J,K,err,cp,H_hfb,R,eHFB_,err_diis,F_diis)

    stop

  end if

! Compute dipole moments

eHF(:)=eHF(:)-chem_pot!TODO remove

  call dipole_moment(nBas,P,nNuc,ZNuc,rNuc,dipole_int,dipole)
  call print_HFB(nBas,nOrb,nO,eHF,c,ENuc,ET,EV,EJ,EK,EL,EHFB,chem_pot,dipole)

! Testing zone

  if(dotest) then
 
    call dump_test_value('R','HFB energy',EHFB)
    call dump_test_value('R','HFB HOMO energy',eHF(nO))
    call dump_test_value('R','HFB LUMO energy',eHF(nO+1))
    call dump_test_value('R','HFB dipole moment',norm2(dipole))

  end if

! Memory deallocation

  deallocate(J,K,err,cp,H_hfb,R,eHFB_,err_diis,F_diis)

end subroutine 



   ! ! DIIS extrapolation TODO check and adapt
   !
   ! if(max_diis > 1) then
   !
   !   n_diis = min(n_diis+1,max_diis)
   !   call DIIS_extrapolation(rcond,nBas_Sq,nBas_Sq,n_diis,err_diis,F_diis,err,F)
   !
   ! end if

   ! ! Level shift TODO check and adapt
   !
   ! if(level_shift > 0d0 .and. Conv > thresh) then
   !   call level_shifting(level_shift,nBas,nOrb,nO,S,c,F)
   ! endif

