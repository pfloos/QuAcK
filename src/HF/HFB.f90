subroutine HFB(dotest,maxSCF,thresh,max_diis,level_shift,nNuc,ZNuc,rNuc,ENuc,                & 
               nBas,nOrb,nOrb_twice,nO,S,T,V,Hc,ERI,dipole_int,X,EHFB,eHF,c,P,Panom,F,Delta, &
               temperature,sigma,chem_pot_hf,chem_pot,restart_hfb,U_QP,eHFB_state)

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
  integer,intent(in)            :: nOrb_twice
  integer,intent(in)            :: nO
  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: temperature,sigma
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas) 
  double precision,intent(in)   :: X(nBas,nOrb)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  logical                       :: chem_pot_hf
  logical                       :: restart_hfb
  integer                       :: nBas2
  integer                       :: iorb
  integer                       :: nSCF
  integer                       :: nBas2_Sq
  integer                       :: n_diis
  double precision              :: ET
  double precision              :: EV
  double precision              :: EJ
  double precision              :: EK
  double precision              :: EL
  double precision              :: Delta_HL
  double precision              :: dipole(ncart)

  double precision              :: Conv
  double precision              :: rcond
  double precision              :: trace_1rdm
  double precision              :: thrs_N
  double precision              :: norm_anom
  double precision,external     :: trace_matrix
  double precision,allocatable  :: eigVAL(:)
  double precision,allocatable  :: Occ(:)
  double precision,allocatable  :: err_diis(:,:)
  double precision,allocatable  :: H_HFB_diis(:,:)
  double precision,allocatable  :: cHF(:,:)
  double precision,allocatable  :: c_tmp(:,:)
  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: K(:,:)
  double precision,allocatable  :: eigVEC(:,:)
  double precision,allocatable  :: H_HFB(:,:)
  double precision,allocatable  :: R(:,:)

  double precision,allocatable  :: err_ao(:,:)
  double precision,allocatable  :: S_ao(:,:)
  double precision,allocatable  :: X_ao(:,:)
  double precision,allocatable  :: c_ao(:,:)
  double precision,allocatable  :: R_ao_old(:,:)
  double precision,allocatable  :: H_HFB_ao(:,:)

! Output variables

  double precision,intent(out)  :: EHFB,chem_pot
  double precision,intent(inout):: c(nBas,nOrb)
  double precision,intent(out)  :: P(nBas,nBas)
  double precision,intent(out)  :: Panom(nBas,nBas)
  double precision,intent(out)  :: F(nBas,nBas)
  double precision,intent(out)  :: Delta(nBas,nBas)
  double precision,intent(out)  :: U_QP(nOrb_twice,nOrb_twice)
  double precision,intent(out)  :: eHFB_state(nOrb_twice)

! Hello world

  write(*,*)
  write(*,*)'*****************************'
  write(*,*)'* HF Bogoliubov Calculation *'
  write(*,*)'*****************************'
  write(*,*)

! Useful quantities

  nBas2 = nBas+nBas
  nBas2_Sq = nBas2*nBas2

! Memory allocation

  allocate(Occ(nOrb))

  allocate(J(nBas,nBas))
  allocate(K(nBas,nBas))

  allocate(cHF(nBas,nOrb))

  allocate(eigVEC(nOrb_twice,nOrb_twice))
  allocate(H_HFB(nOrb_twice,nOrb_twice))
  allocate(R(nOrb_twice,nOrb_twice))
  allocate(eigVAL(nOrb_twice))

  allocate(err_ao(nBas2,nBas2))
  allocate(S_ao(nBas2,nBas2))
  allocate(X_ao(nBas2,nOrb_twice))
  allocate(R_ao_old(nBas2,nBas2))
  allocate(H_HFB_ao(nBas2,nBas2))

  allocate(err_diis(nBas2_Sq,max_diis))
  allocate(H_HFB_diis(nBas2_Sq,max_diis))

! Guess chem. pot.

  chem_pot = 0.5d0*(eHF(nO)+eHF(nO+1))

! Initialization

  thrs_N          = 1d-8
  n_diis          = 0
  H_HFB_diis(:,:) = 0d0
  err_diis(:,:)   = 0d0
  rcond           = 0d0

  P(:,:)         = matmul(c(:,1:nO), transpose(c(:,1:nO)))
  Panom(:,:)     = 0d0
  cHF(:,:)       = c(:,:)

  ! Use Fermi-Dirac occupancies to compute P, Panom, and chem_pot
  
  if(abs(temperature)>1d-4) then
   Occ(:)     = 0d0
   Occ(1:nO)  = 1d0
   call fermi_dirac_occ(nO,nOrb,thrs_N,temperature,chem_pot,Occ,eHF)
   if(chem_pot_hf) chem_pot = 0.5d0*(eHF(nO)+eHF(nO+1))
   P(:,:)      = 0d0
   Panom(:,:)  = 0d0
   do iorb=1,nOrb
    P(:,:)     = P(:,:)     + Occ(iorb)                        * &
                matmul(c(:,iorb:iorb),transpose(c(:,iorb:iorb))) 
    Panom(:,:) = Panom(:,:) + sqrt(abs(Occ(iorb)*(1d0-Occ(iorb))))  * &
                matmul(c(:,iorb:iorb),transpose(c(:,iorb:iorb))) 
   enddo
  endif

  ! Read restart file

  if(restart_hfb) then
   call read_restart_HFB(nBas, nOrb, Occ, c, S, chem_pot)
   P(:,:)      = 0d0
   Panom(:,:)  = 0d0
   do iorb=1,nOrb
    P(:,:)     = P(:,:)     + Occ(iorb)                        * &
                matmul(c(:,iorb:iorb),transpose(c(:,iorb:iorb))) 
    Panom(:,:) = Panom(:,:) + sqrt(abs(Occ(iorb)*(1d0-Occ(iorb))))  * &
                matmul(c(:,iorb:iorb),transpose(c(:,iorb:iorb))) 
   enddo
  endif

  P(:,:)       = 2d0*P(:,:)
  S_ao(:,:)    = 0d0
  S_ao(1:nBas      ,1:nBas      ) = S(1:nBas,1:nBas)
  S_ao(nBas+1:nBas2,nBas+1:nBas2) = S(1:nBas,1:nBas)
  X_ao(:,:)    = 0d0
  X_ao(1:nBas      ,1:nOrb      )      = X(1:nBas,1:nOrb)
  X_ao(nBas+1:nBas2,nOrb+1:nOrb_twice) = X(1:nBas,1:nOrb)

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
    call anomalous_matrix_AO_basis(nBas,sigma,Panom,ERI,Delta)
    
    F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:) - chem_pot*S(:,:)

    ! Diagonalize H_HFB matrix
    
    H_HFB(:,:) = 0d0
    H_HFB(1:nOrb      ,1:nOrb      )           = matmul(transpose(X),matmul(F,X))
    H_HFB(nOrb+1:nOrb_twice,nOrb+1:nOrb_twice) = -H_HFB(1:nOrb,1:nOrb)
    H_HFB(1:nOrb      ,nOrb+1:nOrb_twice) = matmul(transpose(X),matmul(Delta,X))
    H_HFB(nOrb+1:nOrb_twice,1:nOrb      ) = H_HFB(1:nOrb,nOrb+1:nOrb_twice)
    
    eigVEC(:,:) = H_HFB(:,:)
    call diagonalize_matrix(nOrb_twice,eigVEC,eigVAL)
    
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
     call fix_chem_pot(nO,nOrb,nOrb_twice,nSCF,thrs_N,trace_1rdm,chem_pot,H_HFB,eigVEC,R,eigVAL)

    ! DIIS extrapolation

    if(max_diis > 1 .and. nSCF>1) then

     write(*,*) ' Doing DIIS'

     F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:) - chem_pot*S(:,:)
     H_HFB_ao(:,:)    = 0d0
     H_HFB_ao(1:nBas      ,1:nBas      ) =  F(1:nBas,1:nBas)
     H_HFB_ao(nBas+1:nBas2,nBas+1:nBas2) = -F(1:nBas,1:nBas)
     H_HFB_ao(1:nBas      ,nBas+1:nBas2) = Delta(1:nBas,1:nBas)
     H_HFB_ao(nBas+1:nBas2,1:nBas      ) = Delta(1:nBas,1:nBas)
     err_ao = matmul(H_HFB_ao,matmul(R_ao_old,S_ao)) - matmul(matmul(S_ao,R_ao_old),H_HFB_ao)

     n_diis = min(n_diis+1,max_diis)
     call DIIS_extrapolation(rcond,nBas2_Sq,nBas2_Sq,n_diis,err_diis,H_HFB_diis,err_ao,H_HFB_ao)

     H_HFB = matmul(transpose(X_ao),matmul(H_HFB_ao,X_ao))
     eigVEC(:,:) = H_HFB(:,:)
     call diagonalize_matrix(nOrb_twice,eigVEC,eigVAL)
     
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
      call fix_chem_pot(nO,nOrb,nOrb_twice,nSCF,thrs_N,trace_1rdm,chem_pot,H_HFB,eigVEC,R,eigVAL)
   
    end if

    ! Extract P and Panom from R
   
    P(:,:)     = 0d0
    Panom(:,:) = 0d0
    P(:,:)     = 2d0*matmul(X,matmul(R(1:nOrb,1:nOrb),transpose(X)))
    Panom(:,:) = matmul(X,matmul(R(1:nOrb,nOrb+1:nOrb_twice),transpose(X)))

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

     F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:) - chem_pot*S(:,:)
     H_HFB_ao(:,:)    = 0d0
     H_HFB_ao(1:nBas      ,1:nBas      ) =  F(1:nBas,1:nBas)
     H_HFB_ao(nBas+1:nBas2,nBas+1:nBas2) = -F(1:nBas,1:nBas)
     H_HFB_ao(1:nBas      ,nBas+1:nBas2) = Delta(1:nBas,1:nBas)
     H_HFB_ao(nBas+1:nBas2,1:nBas      ) = Delta(1:nBas,1:nBas)
     err_ao = matmul(H_HFB_ao,matmul(R_ao_old,S_ao)) - matmul(matmul(S_ao,R_ao_old),H_HFB_ao)
     Conv  = maxval(abs(err_ao))

    endif

    ! Update R_old

    R_ao_old(:,:)    = 0d0
    R_ao_old(1:nBas      ,1:nBas      ) = 0.5d0*P(1:nBas,1:nBas)
    R_ao_old(nBas+1:nBas2,nBas+1:nBas2) = matmul(X(1:nBas,1:nOrb), transpose(X(1:nBas,1:nOrb)))-0.5d0*P(1:nBas,1:nBas)
    R_ao_old(1:nBas      ,nBas+1:nBas2) = Panom(1:nBas,1:nBas)
    R_ao_old(nBas+1:nBas2,1:nBas      ) = Panom(1:nBas,1:nBas)

    ! Dump results
    write(*,*)'-------------------------------------------------------------------------------------------------'
    write(*,'(1X,A1,1X,A4,1X,A1,1X,A16,1X,A1,1X,A16,1X,A1,1X,A16,1X,A1A16,1X,A1,1X,A10,2X,A1,1X)') &
            '|','#','|','E(HFB)','|','EJ(HFB)','|','EK(HFB)','|','EL(HFB)','|','Conv','|'
    write(*,*)'-------------------------------------------------------------------------------------------------'

    write(*,'(1X,A1,1X,I4,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1XF16.10,1X,A1,1X,E10.2,1X,A1,1X)') &
      '|',nSCF,'|',EHFB + ENuc,'|',EJ,'|',EK,'|',EL,'|',Conv,'|'

    write(*,*)'-------------------------------------------------------------------------------------------------'
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

    deallocate(J,K,eigVEC,H_HFB,R,eigVAL,err_diis,H_HFB_diis,Occ)
    deallocate(err_ao,S_ao,X_ao,R_ao_old,H_HFB_ao,cHF)

    stop

  end if

! Compute final energy before printing summary

  F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:) - chem_pot*S(:,:)

  ! Diagonalize H_HFB matrix
  
  H_HFB(:,:) = 0d0
  H_HFB(1:nOrb      ,1:nOrb      )           = matmul(transpose(X),matmul(F,X))
  H_HFB(nOrb+1:nOrb_twice,nOrb+1:nOrb_twice) = -H_HFB(1:nOrb,1:nOrb)
  H_HFB(1:nOrb      ,nOrb+1:nOrb_twice) = matmul(transpose(X),matmul(Delta,X))
  H_HFB(nOrb+1:nOrb_twice,1:nOrb      ) = H_HFB(1:nOrb,nOrb+1:nOrb_twice)
  
  eigVEC(:,:) = H_HFB(:,:)
  call diagonalize_matrix(nOrb_twice,eigVEC,eigVAL)

  ! Build R and check trace
    
  trace_1rdm = 0d0 
  R(:,:)     = 0d0
  do iorb=1,nOrb
   R(:,:) = R(:,:) + matmul(eigVEC(:,iorb:iorb),transpose(eigVEC(:,iorb:iorb))) 
  enddo

  ! Extract P and Panom from R
 
  P(:,:)     = 0d0
  Panom(:,:) = 0d0
  P(:,:)     = 2d0*matmul(X,matmul(R(1:nOrb,1:nOrb),transpose(X)))
  Panom(:,:) = matmul(X,matmul(R(1:nOrb,nOrb+1:nOrb_twice),transpose(X)))

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

! Compute dipole moments, occupation numbers, || Anomalous density||,
! organize the coefs c with natural orbitals (descending occ numbers), and
! also print the restart file

  call dipole_moment(nBas,P,nNuc,ZNuc,rNuc,dipole_int,dipole)
  eHFB_state(:) = eigVAL(:)
  Delta_HL=eHFB_state(nOrb+1)-eHFB_state(nOrb)
  norm_anom = trace_matrix(nOrb,matmul(transpose(R(1:nOrb,nOrb+1:nOrb_twice)),R(1:nOrb,nOrb+1:nOrb_twice)))
  deallocate(eigVEC,eigVAL)
  allocate(eigVEC(nOrb,nOrb),eigVAL(nOrb))
  eigVEC(:,:) = 0d0
  eigVEC(1:nOrb,1:nOrb) = R(1:nOrb,1:nOrb)
  call diagonalize_matrix(nOrb,eigVEC,eigVAL)
  Occ(1:nOrb)   = eigVAL(1:nOrb)
  c = matmul(X,eigVEC)
  call write_restart_HFB(nBas,nOrb,Occ,c,chem_pot) ! Warning: orders Occ and their c in descending order w.r.t. occupation numbers.
  call print_HFB(nBas,nOrb,nOrb_twice,nO,norm_anom,Occ,eHFB_state,ENuc,ET,EV,EJ,EK,EL,EHFB,chem_pot, &
                 dipole,Delta_HL)

! Choose the NO representation where the 1-RDM is diag.

  allocate(c_ao(nBas2,nOrb_twice))

  if(.true.) then ! NO basis

! Compute W_no and V_no (i.e. diag[H_HFB^no] built in NO basis and get W and V)

   deallocate(eigVEC,eigVAL)
   allocate(eigVEC(nOrb_twice,nOrb_twice),eigVAL(nOrb_twice))

   c_ao(:,:) = 0d0
   c_ao(1:nBas      ,1:nOrb      )      = c(1:nBas,1:nOrb)
   c_ao(nBas+1:nBas2,nOrb+1:nOrb_twice) = c(1:nBas,1:nOrb)
   H_HFB = matmul(transpose(c_ao),matmul(H_HFB_ao,c_ao)) ! H_HFB is in the NO basis
   eigVEC(:,:) = H_HFB(:,:)

   call diagonalize_matrix(nOrb_twice,eigVEC,eigVAL)

   ! Build R (as R^no) and save the eigenvectors
     
   trace_1rdm = 0d0 
   R(:,:)     = 0d0
   do iorb=1,nOrb
    R(:,:) = R(:,:) + matmul(eigVEC(:,iorb:iorb),transpose(eigVEC(:,iorb:iorb))) 
   enddo
   U_QP(:,:) = eigVEC(:,:)

   ! Check trace of R
   do iorb=1,nOrb
    trace_1rdm = trace_1rdm + R(iorb,iorb)
   enddo
   trace_1rdm = 2d0*trace_1rdm
   write(*,*)
   write(*,'(A33,1X,F16.10,A3)') ' Trace [ 1D^NO ]     = ',trace_1rdm,'   '
   write(*,*)
 
  endif 

  deallocate(c_ao)

! Testing zone

  if(dotest) then
 
    call dump_test_value('R','HFB energy',EHFB)
    call dump_test_value('R','Trace 1D',trace_1rdm)
    call dump_test_value('R','HFB dipole moment',norm2(dipole))

  end if

! Memory deallocation

  deallocate(J,K,eigVEC,H_HFB,R,eigVAL,err_diis,H_HFB_diis,Occ)
  deallocate(err_ao,S_ao,X_ao,R_ao_old,H_HFB_ao,cHF)

end subroutine 

