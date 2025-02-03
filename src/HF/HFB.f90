subroutine HFB(dotest,maxSCF,thresh,max_diis,level_shift,nNuc,ZNuc,rNuc,ENuc,     & 
               nBas,nOrb,nO,S,T,V,Hc,ERI,dipole_int,X,EHFB,eHF,c,P,Panom,F,Delta, &
               temperature)

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
  double precision,intent(in)   :: temperature
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas) 
  double precision,intent(in)   :: X(nBas,nOrb)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  integer                       :: iorb
  integer                       :: nSCF
  integer                       :: nOrb2
  integer                       :: nOrb_Sq
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
  double precision              :: trace_1rdm
  double precision              :: thrs_N
  double precision              :: norm_anom
  double precision,external     :: trace_matrix
  double precision,allocatable  :: eigVAL(:)
  double precision,allocatable  :: Occ(:)
  double precision,allocatable  :: err(:,:)
  double precision,allocatable  :: err_diis(:,:)
  double precision,allocatable  :: H_HFB_diis(:,:)
  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: K(:,:)
  double precision,allocatable  :: eigVEC(:,:)
  double precision,allocatable  :: H_HFB(:,:)
  double precision,allocatable  :: R(:,:)
  double precision,allocatable  :: R_old(:,:)

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

  nOrb_Sq = nOrb*nOrb
  nOrb2 = nOrb+nOrb

! Memory allocation

  allocate(J(nBas,nBas))
  allocate(K(nBas,nBas))

  allocate(err(nOrb2,nOrb2))

  allocate(eigVEC(nOrb2,nOrb2))
  allocate(H_HFB(nOrb2,nOrb2))
  allocate(R(nOrb2,nOrb2))
  allocate(R_old(nOrb2,nOrb2))
  allocate(eigVAL(nOrb2))

  allocate(err_diis(nOrb_Sq,max_diis))
  allocate(H_HFB_diis(nOrb_Sq,max_diis))

! Guess chem. pot.

  chem_pot = 0.5d0*(eHF(nO)+eHF(nO+1))

! Initialization

  thrs_N          = 1d-8
  n_diis          = 0
  H_HFB_diis(:,:)     = 0d0
  err_diis(:,:)   = 0d0
  rcond           = 0d0

  P(:,:)         = matmul(c(:,1:nO), transpose(c(:,1:nO)))
  Panom(:,:)     = 0d0

  ! Use Fermi-Dirac occupancies to compute P, Panom, and chem_pot
  
  if(abs(temperature)>1d-4) then
   allocate(Occ(nOrb))
   Occ(:)     = 0d0
   Occ(1:nO)  = 1d0
   call fermi_dirac_occ(nO,nOrb,thrs_N,temperature,chem_pot,Occ,eHF)
   P(:,:)      = 0d0
   Panom(:,:)  = 0d0
   do iorb=1,nOrb
    P(:,:)     = P(:,:)     + Occ(iorb)                        * &
                matmul(c(:,iorb:iorb),transpose(c(:,iorb:iorb))) 
    Panom(:,:) = Panom(:,:) + sqrt(Occ(iorb)*(1d0-Occ(iorb)))  * &
                matmul(c(:,iorb:iorb),transpose(c(:,iorb:iorb))) 
   enddo
   deallocate(Occ)
  endif

  P(:,:)       = 2d0*P(:,:)

  Conv = 1d0
  nSCF = 0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------

  write(*,*)
  write(*,*) 'Enterning HFB SCF procedure'  
  write(*,*)
  do while(Conv > thresh .and. nSCF < maxSCF)
   
 
    ! Increment 

    nSCF = nSCF + 1

    ! Build Fock and Delta matrices
    
    call Hartree_matrix_AO_basis(nBas,P,ERI,J)
    call exchange_matrix_AO_basis(nBas,P,ERI,K)
    call anomalous_matrix_AO_basis(nBas,Panom,ERI,Delta)
    
    F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:) - chem_pot*S(:,:)

    ! Diagonalize H_HFB matrix
    
    H_HFB(:,:) = 0d0
    H_HFB(1:nOrb      ,1:nOrb      ) = matmul(transpose(X),matmul(F,X))
    H_HFB(nOrb+1:nOrb2,nOrb+1:nOrb2) = -H_HFB(1:nOrb,1:nOrb)
    H_HFB(1:nOrb      ,nOrb+1:nOrb2) = matmul(transpose(X),matmul(Delta,X))
    H_HFB(nOrb+1:nOrb2,1:nOrb      ) = transpose(H_HFB(1:nOrb,nOrb+1:nOrb2))
    
    eigVEC(:,:) = H_HFB(:,:)
    call diagonalize_matrix(nOrb2,eigVEC,eigVAL)
    
    ! Build R 
      
    R(:,:)     = 0d0
    do iorb=1,nOrb
     R(:,:) = R(:,:) + matmul(eigVEC(:,iorb:iorb),transpose(eigVEC(:,iorb:iorb))) 
    enddo
    trace_1rdm = 0d0
    do iorb=1,nOrb
     trace_1rdm = trace_1rdm+R(iorb,iorb) 
    enddo

    ! Adjust the chemical potential 

    if( abs(trace_1rdm-nO) > thrs_N ) & 
     call fix_chem_pot(nO,nOrb,nOrb2,nSCF,thrs_N,trace_1rdm,chem_pot,H_HFB,eigVEC,R,eigVAL)

   ! DIIS extrapolation TODO

    !if(max_diis > 1 .and. .false.) then
    if(max_diis > 1 .and. nSCF>1) then

     write(*,*) ' Doing DIIS'

     err   = matmul(H_HFB,R_old) - matmul(R_old,H_HFB)
     n_diis = min(n_diis+1,max_diis)
     call DIIS_extrapolation(rcond,nOrb_Sq,nOrb_Sq,n_diis,err_diis,H_HFB_diis,err,H_HFB)

     eigVEC(:,:) = H_HFB(:,:)
     call diagonalize_matrix(nOrb2,eigVEC,eigVAL)
     
     ! Build R and check trace
       
     trace_1rdm = 0d0 
     R(:,:)     = 0d0
     do iorb=1,nOrb
      R(:,:) = R(:,:) + matmul(eigVEC(:,iorb:iorb),transpose(eigVEC(:,iorb:iorb))) 
     enddo
     do iorb=1,nOrb
      trace_1rdm = trace_1rdm + R(iorb,iorb) 
     enddo

     ! Adjust the chemical potential 
     
     if( abs(trace_1rdm-nO) > thrs_N ) & 
      call fix_chem_pot(nO,nOrb,nOrb2,nSCF,thrs_N,trace_1rdm,chem_pot,H_HFB,eigVEC,R,eigVAL)
   
    end if

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
    EL = trace_matrix(nBas,matmul(Panom,Delta))
    ! Total energy
    EHFB = ET + EV + EJ + EK + EL

    ! Check convergence

    if(nSCF > 1) then

     err   = matmul(H_HFB,R_old) - matmul(R_old,H_HFB)
     Conv  = maxval(abs(err))

    endif

    ! Update R_old

    R_old(:,:) = R(:,:)

    ! Dump results
    write(*,*)'-----------------------------------------------------------------------------------------------'
    write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A16,1X,A1,1X,A16,1X,A1A16,1X,A1,1X,A10,2X,A1,1X)') &
            '|','#','|','E(HFB)','|','EJ(HFB)','|','EK(HFB)','|','EL(HFB)','|','Conv','|'
    write(*,*)'-----------------------------------------------------------------------------------------------'

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1XF16.10,1X,A1,1X,E10.2,1X,A1,1X)') &
      '|',nSCF,'|',EHFB + ENuc,'|',EJ,'|',EK,'|',EL,'|',Conv,'|'

    write(*,*)'-----------------------------------------------------------------------------------------------'
    write(*,*)

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

    deallocate(J,K,err,eigVEC,H_HFB,R,R_old,eigVAL,err_diis,H_HFB_diis)

    stop

  end if

! Compute dipole moments, occupation numbers and || Anomalous density||

  eigVEC(:,:) = 0d0
  eigVEC(1:nOrb,1:nOrb) = R(1:nOrb,1:nOrb)
  call diagonalize_matrix(nOrb2,eigVEC,eigVAL)
  eigVAL(:)   = 2d0*eigVAL(:)
  norm_anom = trace_matrix(nOrb,matmul(transpose(R(1:nOrb,nOrb+1:nOrb2)),R(1:nOrb,nOrb+1:nOrb2)))
  call dipole_moment(nBas,P,nNuc,ZNuc,rNuc,dipole_int,dipole)
  call print_HFB(nBas,nOrb,nO,norm_anom,eigVAL,ENuc,ET,EV,EJ,EK,EL,EHFB,chem_pot,dipole)

! Testing zone

  if(dotest) then
 
    call dump_test_value('R','HFB energy',EHFB)
    call dump_test_value('R','HFB HOMO energy',eHF(nO))
    call dump_test_value('R','HFB LUMO energy',eHF(nO+1))
    call dump_test_value('R','HFB dipole moment',norm2(dipole))

  end if

! Memory deallocation

  deallocate(J,K,err,eigVEC,H_HFB,R,R_old,eigVAL,err_diis,H_HFB_diis)

end subroutine 

