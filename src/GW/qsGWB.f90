subroutine qsGWB(dotest,maxSCF,thresh,max_diis,level_shift,nNuc,ZNuc,rNuc,ENuc,eta,shift,        & 
               nBas,nOrb,nOrb_twice,nO,S,T,V,Hc,ERI,dipole_int,X,EqsGWB,eqsGW,c,P,Panom,F,Delta, &
               sigma,chem_pot,restart_hfb,U_QP,eqsGWB_state,nfreqs,ntimes,wcoord,wweight)

! Perform qsGW Bogoliubov calculation

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
  integer,intent(in)            :: nfreqs
  integer,intent(in)            :: ntimes
  double precision,intent(in)   :: wcoord(nfreqs),wweight(nfreqs)
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: shift
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: sigma
  double precision,intent(in)   :: eqsGW(nOrb)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas) 
  double precision,intent(in)   :: X(nBas,nOrb)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  logical                       :: restart_hfb
  integer                       :: nBas2
  integer                       :: iorb,jorb,korb,lorb
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
  double precision,allocatable  :: H_qsGWB_diis(:,:)
  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: K(:,:)
  double precision,allocatable  :: Sigc_he(:,:)
  double precision,allocatable  :: Sigc_hh(:,:)
  double precision,allocatable  :: eigVEC(:,:)
  double precision,allocatable  :: H_qsGWB(:,:)
  double precision,allocatable  :: R(:,:)
  double precision,allocatable  :: vMAT(:,:)
  double precision,allocatable  :: ERI_MO(:,:,:,:)

  double precision,allocatable  :: err_ao(:,:)
  double precision,allocatable  :: S_ao(:,:)
  double precision,allocatable  :: X_ao(:,:)
  double precision,allocatable  :: c_ao(:,:)
  double precision,allocatable  :: R_ao_old(:,:)
  double precision,allocatable  :: H_qsGWB_ao(:,:)

! Output variables

  double precision,intent(inout) :: chem_pot
  double precision,intent(out)   :: EqsGWB
  double precision,intent(inout) :: c(nBas,nOrb)
  double precision,intent(inout) :: P(nBas,nBas)
  double precision,intent(inout) :: Panom(nBas,nBas)
  double precision,intent(inout) :: F(nBas,nBas)
  double precision,intent(inout) :: Delta(nBas,nBas)
  double precision,intent(inout) :: U_QP(nOrb_twice,nOrb_twice)
  double precision,intent(inout) :: eqsGWB_state(nOrb_twice)

! Hello world

  write(*,*)
  write(*,*)'*******************************'
  write(*,*)'* qsGW Bogoliubov Calculation *'
  write(*,*)'*******************************'
  write(*,*)

! Useful quantities

  nBas2 = nBas+nBas
  nBas2_Sq = nBas2*nBas2

! Memory allocation

  allocate(Occ(nOrb))

  allocate(J(nBas,nBas))
  allocate(K(nBas,nBas))
  allocate(Sigc_he(nBas,nBas))
  allocate(Sigc_hh(nBas,nBas))

  allocate(eigVEC(nOrb_twice,nOrb_twice))
  allocate(H_qsGWB(nOrb_twice,nOrb_twice))
  allocate(R(nOrb_twice,nOrb_twice))
  allocate(eigVAL(nOrb_twice))
  allocate(vMAT(nOrb*nOrb,nOrb*nOrb))

  allocate(err_ao(nBas2,nBas2))
  allocate(S_ao(nBas2,nBas2))
  allocate(X_ao(nBas2,nOrb_twice))
  allocate(R_ao_old(nBas2,nBas2))
  allocate(H_qsGWB_ao(nBas2,nBas2))

  allocate(err_diis(nBas2_Sq,max_diis))
  allocate(H_qsGWB_diis(nBas2_Sq,max_diis))

  allocate(ERI_MO(nOrb,nOrb,nOrb,nOrb))
  call AOtoMO_ERI_RHF(nBas,nOrb,X,ERI,ERI_MO)
  do iorb=1,nOrb
   do jorb=1,nOrb
    do korb=1,nOrb
     do lorb=1,nOrb
      vMAT(1+(korb-1)+(iorb-1)*nOrb,1+(lorb-1)+(jorb-1)*nOrb)=ERI_MO(iorb,jorb,korb,lorb)
     enddo
    enddo
   enddo
  enddo
  deallocate(ERI_MO)

! Initialization

  thrs_N          = 1d-8
  n_diis          = 0
  H_qsGWB_diis(:,:) = 0d0
  err_diis(:,:)   = 0d0
  rcond           = 0d0

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
  write(*,*) 'Enterning qsGWB SCF procedure'  
  write(*,*)
  do while(Conv > thresh .and. nSCF < maxSCF)
   
 
    ! Increment 

    nSCF = nSCF + 1

    ! Build Fock and Delta matrices
    
    call Hartree_matrix_AO_basis(nBas,P,ERI,J)
    call exchange_matrix_AO_basis(nBas,P,ERI,K)
    call anomalous_matrix_AO_basis(nBas,sigma,Panom,ERI,Delta)
    if(Conv<1d-4 .or. nSCF==1) then
     call sigc_AO_basis_HFB(nBas,nOrb,nOrb_twice,eta,shift,X,U_QP,eqsGWB_state, & 
                            S,vMAT,nfreqs,ntimes,wcoord,wweight,Sigc_he,Sigc_hh)
    endif
    
    F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:) + Sigc_he(:,:) - chem_pot*S(:,:)

    ! Diagonalize H_qsGWB matrix
    
    H_qsGWB(:,:) = 0d0
    H_qsGWB(1:nOrb      ,1:nOrb      )           = matmul(transpose(X),matmul(F,X))
    H_qsGWB(nOrb+1:nOrb_twice,nOrb+1:nOrb_twice) = -H_qsGWB(1:nOrb,1:nOrb)
    H_qsGWB(1:nOrb      ,nOrb+1:nOrb_twice) = matmul(transpose(X),matmul(Delta+Sigc_hh,X))
    H_qsGWB(nOrb+1:nOrb_twice,1:nOrb      ) = H_qsGWB(1:nOrb,nOrb+1:nOrb_twice)
    
    eigVEC(:,:) = H_qsGWB(:,:)
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
     call fix_chem_pot(nO,nOrb,nOrb_twice,nSCF,thrs_N,trace_1rdm,chem_pot,H_qsGWB,eigVEC,R,eigVAL)

    ! DIIS extrapolation

    if(max_diis > 1 .and. nSCF>1) then

     write(*,*) ' Doing DIIS'

     F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:) + Sigc_he(:,:) - chem_pot*S(:,:)
     H_qsGWB_ao(:,:)    = 0d0
     H_qsGWB_ao(1:nBas      ,1:nBas      ) =  F(1:nBas,1:nBas)
     H_qsGWB_ao(nBas+1:nBas2,nBas+1:nBas2) = -F(1:nBas,1:nBas)
     H_qsGWB_ao(1:nBas      ,nBas+1:nBas2) = Delta(1:nBas,1:nBas) + Sigc_hh(1:nBas,1:nBas)
     H_qsGWB_ao(nBas+1:nBas2,1:nBas      ) = Delta(1:nBas,1:nBas) + Sigc_hh(1:nBas,1:nBas)
     err_ao = matmul(H_qsGWB_ao,matmul(R_ao_old,S_ao)) - matmul(matmul(S_ao,R_ao_old),H_qsGWB_ao)

     n_diis = min(n_diis+1,max_diis)
     call DIIS_extrapolation(rcond,nBas2_Sq,nBas2_Sq,n_diis,err_diis,H_qsGWB_diis,err_ao,H_qsGWB_ao)

     H_qsGWB = matmul(transpose(X_ao),matmul(H_qsGWB_ao,X_ao))
     eigVEC(:,:) = H_qsGWB(:,:)
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
      call fix_chem_pot(nO,nOrb,nOrb_twice,nSCF,thrs_N,trace_1rdm,chem_pot,H_qsGWB,eigVEC,R,eigVAL)
   
    end if

    ! Extract P and Panom from R
    
    eqsGWB_state(:) = eigVAL(:)
    U_QP(:,:)  = eigVEC(:,:)
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
    EqsGWB = ET + EV + EJ + EK + EL

    ! Check convergence

    if(nSCF > 1) then

     F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:) + Sigc_he(:,:) - chem_pot*S(:,:)
     H_qsGWB_ao(:,:)    = 0d0
     H_qsGWB_ao(1:nBas      ,1:nBas      ) =  F(1:nBas,1:nBas)
     H_qsGWB_ao(nBas+1:nBas2,nBas+1:nBas2) = -F(1:nBas,1:nBas)
     H_qsGWB_ao(1:nBas      ,nBas+1:nBas2) = Delta(1:nBas,1:nBas) + Sigc_hh(1:nBas,1:nBas)
     H_qsGWB_ao(nBas+1:nBas2,1:nBas      ) = Delta(1:nBas,1:nBas) + Sigc_hh(1:nBas,1:nBas)
     err_ao = matmul(H_qsGWB_ao,matmul(R_ao_old,S_ao)) - matmul(matmul(S_ao,R_ao_old),H_qsGWB_ao)
     Conv  = maxval(abs(err_ao))

    endif

    ! Update R_old

    R_ao_old(:,:)    = 0d0
    R_ao_old(1:nBas      ,1:nBas      ) = 0.5d0*P(1:nBas,1:nBas)
    R_ao_old(nBas+1:nBas2,nBas+1:nBas2) = matmul(X(1:nBas,1:nOrb), transpose(X(1:nBas,1:nOrb)))-0.5d0*P(1:nBas,1:nBas)
    R_ao_old(1:nBas      ,nBas+1:nBas2) = Panom(1:nBas,1:nBas)
    R_ao_old(nBas+1:nBas2,1:nBas      ) = Panom(1:nBas,1:nBas)

    ! Dump results
    write(*,*)'------------------------------------------------------------------------------------------&
    &-------'
    write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A16,1X,A1,1X,A16,1X,A1,1X,A16,1X,A1,1X,A10,2X,A1,1X)') &
            '|','#','|','E(qsGWB)','|','EJ(qsGWB)','|','EK(qsGWB)','|','EL(qsGWB)','|','Conv','|'
    write(*,*)'------------------------------------------------------------------------------------------&
    &-------'

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,&
     &  1X,E10.2,1X,A2,1X)') &
      '|',nSCF,'|',EqsGWB + ENuc,'|',EJ,'|',EK,'|',EL,'|',Conv,' |'

    write(*,*)'------------------------------------------------------------------------------------------&
    &-------'
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

    deallocate(J,K,Sigc_he,Sigc_hh,eigVEC,H_qsGWB,R,eigVAL,err_diis,H_qsGWB_diis,Occ)
    deallocate(err_ao,S_ao,X_ao,R_ao_old,H_qsGWB_ao)
    deallocate(vMAT)

    stop

  end if

! Compute dipole moments, occupation numbers, || Anomalous density||,
! organize the coefs c with natural orbitals (descending occ numbers), and
! also print the restart file

  eqsGWB_state(:) = eigVAL(:)
  Delta_HL=eqsGWB_state(nOrb+1)-eqsGWB_state(nOrb)
  deallocate(eigVEC,eigVAL)
  allocate(eigVEC(nOrb,nOrb),eigVAL(nOrb))
  eigVEC(1:nOrb,1:nOrb) = 0d0
  eigVEC(1:nOrb,1:nOrb) = R(1:nOrb,1:nOrb)
  call diagonalize_matrix(nOrb,eigVEC,eigVAL)
  Occ(1:nOrb)   = eigVAL(1:nOrb)
  c = matmul(X,eigVEC)
  norm_anom = trace_matrix(nOrb,matmul(transpose(R(1:nOrb,nOrb+1:nOrb_twice)),R(1:nOrb,nOrb+1:nOrb_twice)))
  call dipole_moment(nBas,P,nNuc,ZNuc,rNuc,dipole_int,dipole)
  call write_restart_qsGWB(nBas,nOrb,Occ,c,chem_pot) ! orders Occ and their c in descending order w.r.t. occupation numbers.
  call print_qsGWB(nBas,nOrb,nOrb_twice,nO,norm_anom,Occ,eqsGWB_state,ENuc,ET,EV,EJ,EK,EL,EqsGWB,chem_pot,&
                   dipole,Delta_HL)

! Compute W_no and V_no (i.e. diag[H_qsGWB^no] built in NO basis and get W and V).

  deallocate(eigVEC,eigVAL)
  allocate(eigVEC(nOrb_twice,nOrb_twice),eigVAL(nOrb_twice),c_ao(nBas2,nOrb_twice))
  c_ao(:,:)    = 0d0
  c_ao(1:nBas      ,1:nOrb      )      = c(1:nBas,1:nOrb)
  c_ao(nBas+1:nBas2,nOrb+1:nOrb_twice) = c(1:nBas,1:nOrb)
  H_qsGWB = matmul(transpose(c_ao),matmul(H_qsGWB_ao,c_ao)) ! H_qsGWB is in the NO basis
  eigVEC(:,:) = H_qsGWB(:,:)
  call diagonalize_matrix(nOrb_twice,eigVEC,eigVAL)
  deallocate(c_ao)
  
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

! Testing zone

  if(dotest) then
 
    call dump_test_value('R','qsGWB energy',EqsGWB)
    call dump_test_value('R','Trace 1D',trace_1rdm)
    call dump_test_value('R','qsGWB dipole moment',norm2(dipole))

  end if

! Memory deallocation

  deallocate(J,K,Sigc_he,Sigc_hh,eigVEC,H_qsGWB,R,eigVAL,err_diis,H_qsGWB_diis,Occ)
  deallocate(err_ao,S_ao,X_ao,R_ao_old,H_qsGWB_ao)
  deallocate(vMAT)

end subroutine 

