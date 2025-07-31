subroutine qsRGWBim(dotest,maxSCF,thresh,max_diis,level_shift,nNuc,ZNuc,rNuc,ENuc,eta,shift, & 
               nBas,nOrb,nOrb_twice,nO,verbose,S,T,V,Hc,ERI,dipole_int,X,EqsGWB,c,P,Panom,   &
               F,Delta,sigma,chem_pot,U_QP,eqsGWB_state,nfreqs,ntimes,wcoord,wweight)

! Perform qsGW Bogoliubov calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh
  double precision,intent(in)   :: level_shift

  integer,intent(in)            :: verbose
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
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas) 
  double precision,intent(in)   :: X(nBas,nOrb)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  logical                       :: offdiag0

  integer                       :: nBas_twice
  integer                       :: iorb,jorb,korb,lorb
  integer                       :: nSCF
  integer                       :: nBas_twice_Sq
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
  double precision              :: pivot_Eqp
  double precision              :: pivot_Occ
  double precision              :: trace_1rdm
  double precision              :: thrs_N
  double precision              :: N_anom
  double precision              :: EcRPA,EcGM
  double precision,external     :: trace_matrix
  double precision,allocatable  :: eigVAL(:)
  double precision,allocatable  :: Occ(:)
  double precision,allocatable  :: eFock(:)
  double precision,allocatable  :: pivot_U_QP(:)
  double precision,allocatable  :: err_diis(:,:)
  double precision,allocatable  :: H_qsGWB_diis(:,:)
  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: K(:,:)
  double precision,allocatable  :: Sigc_ao_he(:,:)
  double precision,allocatable  :: Sigc_ao_hh(:,:)
  double precision,allocatable  :: Sigc_mo_he(:,:)
  double precision,allocatable  :: Sigc_mo_hh(:,:)
  double precision,allocatable  :: H_qsGWB(:,:)
  double precision,allocatable  :: U(:,:)
  double precision,allocatable  :: R(:,:)
  double precision,allocatable  :: vMAT(:,:)
  double precision,allocatable  :: ERI_MO(:,:,:,:)

  double precision,allocatable  :: err_ao(:,:)
  double precision,allocatable  :: S_ao(:,:)
  double precision,allocatable  :: X_ao(:,:)
  double precision,allocatable  :: c_no(:,:)
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

  nBas_twice    = nBas+nBas
  nBas_twice_Sq = nBas_twice*nBas_twice
  offdiag0=.false. ! Set it to true if you want to try qsGWB version 2 where all the off-diag. elements of Sigc_mo are eval at Fermi level

! Memory allocation

  allocate(Occ(nOrb))
  allocate(eFock(nOrb))

  allocate(Sigc_mo_he(nOrb,nOrb))
  allocate(Sigc_mo_hh(nOrb,nOrb))

  allocate(J(nBas,nBas))
  allocate(K(nBas,nBas))
  allocate(Sigc_ao_he(nBas,nBas))
  allocate(Sigc_ao_hh(nBas,nBas))

  allocate(pivot_U_QP(nOrb_twice))
  allocate(H_qsGWB(nOrb_twice,nOrb_twice))
  allocate(R(nOrb_twice,nOrb_twice))
  allocate(eigVAL(nOrb_twice))
  allocate(vMAT(nOrb*nOrb,nOrb*nOrb))

  allocate(err_ao(nBas_twice,nBas_twice))
  allocate(S_ao(nBas_twice,nBas_twice))
  allocate(X_ao(nBas_twice,nOrb_twice))
  allocate(R_ao_old(nBas_twice,nBas_twice))
  allocate(H_qsGWB_ao(nBas_twice,nBas_twice))

  allocate(err_diis(nBas_twice_Sq,max_diis))
  allocate(H_qsGWB_diis(nBas_twice_Sq,max_diis))


! Initialization

  thrs_N            = 1d-8
  n_diis            = 0
  H_qsGWB_diis(:,:) = 0d0
  err_diis(:,:)     = 0d0
  rcond             = 0d0
  Conv              = 1d0
  nSCF              = 0

  S_ao(:,:)       = 0d0
  X_ao(:,:)       = 0d0
  S_ao(1:nBas      ,1:nBas      )           = S(1:nBas,1:nBas)
  S_ao(nBas+1:nBas_twice,nBas+1:nBas_twice) = S(1:nBas,1:nBas)
  X_ao(1:nBas      ,1:nOrb      )           = X(1:nBas,1:nOrb)
  X_ao(nBas+1:nBas_twice,nOrb+1:nOrb_twice) = X(1:nBas,1:nOrb)

! Compute initial U_QP, HFB energies, and electronic major Occ numbers 

  call Hartree_matrix_AO_basis(nBas,P,ERI,J)
  call exchange_matrix_AO_basis(nBas,P,ERI,K)
  call anomalous_matrix_AO_basis(nBas,sigma,Panom,ERI,Delta)
  F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:) - chem_pot*S(:,:)
  allocate(c_ao(nBas_twice,nOrb_twice),U(nOrb,nOrb))
  U = matmul(transpose(c),matmul(F,c))
  call diagonalize_matrix(nOrb,U,eFock)
  c=matmul(c,U)
  c_ao(:,:)       = 0d0
  c_ao(1:nBas           ,1:nOrb           ) = c(1:nBas,1:nOrb)
  c_ao(nBas+1:nBas_twice,nOrb+1:nOrb_twice) = c(1:nBas,1:nOrb)
  H_qsGWB_ao(:,:)  = 0d0
  H_qsGWB_ao(1:nBas           ,1:nBas           ) =  F(1:nBas,1:nBas)
  H_qsGWB_ao(nBas+1:nBas_twice,nBas+1:nBas_twice) = -F(1:nBas,1:nBas)
  H_qsGWB_ao(1:nBas           ,nBas+1:nBas_twice) = Delta(1:nBas,1:nBas)
  H_qsGWB_ao(nBas+1:nBas_twice,1:nBas           ) = Delta(1:nBas,1:nBas)
  U_QP = matmul(transpose(c_ao),matmul(H_qsGWB_ao,c_ao))
  call diagonalize_matrix(nOrb_twice,U_QP,eigVAL)
  Occ=0d0
  do iorb=1,nOrb
   do jorb=1,nOrb
    Occ(iorb)=Occ(iorb)+U_QP(jorb,iorb)**2d0
   enddo
  enddo
  do iorb=1,nOrb
   do jorb=iorb,nOrb
    if(Occ(jorb)>Occ(iorb)) then
     pivot_U_QP(:)= U_QP(:,iorb)
     U_QP(:,iorb) = U_QP(:,jorb)
     U_QP(:,jorb) = pivot_U_QP(:)
     pivot_Eqp    = eigVAL(iorb)
     eigVAL(iorb) = eigVAL(jorb)
     eigVAL(jorb) = pivot_Eqp  
     pivot_Occ    = Occ(iorb)
     Occ(iorb)    = Occ(jorb)
     Occ(jorb)    = pivot_Occ
    endif
   enddo
  enddo
  do iorb=1,nOrb
   do jorb=iorb,nOrb
    if(abs(Occ(iorb))>0.5d0 .and. abs(Occ(jorb))>0.5d0 .and. eigVAL(jorb)<eigVAL(iorb)) then
     pivot_U_QP(:)= U_QP(:,iorb)
     U_QP(:,iorb) = U_QP(:,jorb)
     U_QP(:,jorb) = pivot_U_QP(:)
     pivot_Eqp    = eigVAL(iorb)
     eigVAL(iorb) = eigVAL(jorb)
     eigVAL(jorb) = pivot_Eqp  
     pivot_Occ    = Occ(iorb)
     Occ(iorb)    = Occ(jorb)
     Occ(jorb)    = pivot_Occ
    endif
   enddo
  enddo
  do iorb=1,nOrb
   do jorb=iorb,nOrb
    if(abs(Occ(iorb))<=0.5d0 .and. abs(Occ(jorb))<=0.5d0 .and. eigVAL(jorb)>eigVAL(iorb)) then
     pivot_U_QP(:)= U_QP(:,iorb)
     U_QP(:,iorb) = U_QP(:,jorb)
     U_QP(:,jorb) = pivot_U_QP(:)
     pivot_Eqp    = eigVAL(iorb)
     eigVAL(iorb) = eigVAL(jorb)
     eigVAL(jorb) = pivot_Eqp  
     pivot_Occ    = Occ(iorb)
     Occ(iorb)    = Occ(jorb)
     Occ(jorb)    = pivot_Occ
    endif
   enddo
  enddo
  deallocate(U,c_ao)
  eqsGWB_state(:) = eigVAL(:)


! Compute vMAT in the c basis 

  allocate(ERI_MO(nOrb,nOrb,nOrb,nOrb))
  call AOtoMO_ERI_RHF(nBas,nOrb,c,ERI,ERI_MO)
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
  write(*,*)
  write(*,*) 'Computing phRPA@HFB and GM@HFB energies'
  write(*,*) 'assuming a previous HFB calcution was performed'
  write(*,*)
  call EcRPA_EcGM_w_RHFB(nOrb,nOrb_twice,1,eqsGWB_state,nfreqs,ntimes,wweight,wcoord,vMAT, &
                         U_QP,EqsGWB+ENuc,EcRPA,EcGM)

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
    allocate(ERI_MO(nOrb,nOrb,nOrb,nOrb))
    call AOtoMO_ERI_RHF(nBas,nOrb,c,ERI,ERI_MO)
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
    call sigc_MO_basis_RHFB(nOrb,nOrb_twice,offdiag0,eta,shift,Occ,U_QP,eqsGWB_state, & 
                            vMAT,nfreqs,ntimes,wcoord,wweight,Sigc_mo_he,Sigc_mo_hh)
    Sigc_mo_he = 0.5d0 * (Sigc_mo_he + transpose(Sigc_mo_he))
    Sigc_mo_hh = 0.5d0 * (Sigc_mo_hh + transpose(Sigc_mo_hh))
    if(verbose/=0) then
     write(*,*) 'Sigma_c he MO'
     do iorb=1,nOrb
      write(*,'(*(f10.5))') Sigc_mo_he(iorb,:)
     enddo
     write(*,*) 'Sigma_c hh MO'
     do iorb=1,nOrb
      write(*,'(*(f10.5))') Sigc_mo_hh(iorb,:)
     enddo
    endif
    call MOtoAO(nBas,nOrb,S,c,Sigc_mo_he,Sigc_ao_he)
    call MOtoAO(nBas,nOrb,S,c,Sigc_mo_hh,Sigc_ao_hh)
    if(verbose/=0) then
     write(*,*) 'Sigma_c he AO'
     do iorb=1,nBas
      write(*,'(*(f10.5))') Sigc_ao_he(iorb,:)
     enddo
     write(*,*) 'Sigma_c hh AO'
     do iorb=1,nBas
      write(*,'(*(f10.5))') Sigc_ao_hh(iorb,:)
     enddo
    endif


    F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:) + Sigc_ao_he(:,:) - chem_pot*S(:,:)

    ! Diagonalize H_qsGWB matrix
    
    H_qsGWB(:,:) = 0d0
    H_qsGWB(1:nOrb           ,1:nOrb           ) = matmul(transpose(X),matmul(F,X))
    H_qsGWB(nOrb+1:nOrb_twice,nOrb+1:nOrb_twice) = -H_qsGWB(1:nOrb,1:nOrb)
    H_qsGWB(1:nOrb           ,nOrb+1:nOrb_twice) = matmul(transpose(X),matmul(Delta+Sigc_ao_hh,X))
    H_qsGWB(nOrb+1:nOrb_twice,1:nOrb           ) = H_qsGWB(1:nOrb,nOrb+1:nOrb_twice)
    
    U_QP(:,:) = H_qsGWB(:,:)
    call diagonalize_matrix(nOrb_twice,U_QP,eigVAL)

    ! Build R 
      
    R(:,:)     = 0d0
    do iorb=1,nOrb
     R(:,:) = R(:,:) + matmul(U_QP(:,iorb:iorb),transpose(U_QP(:,iorb:iorb))) 
    enddo
    trace_1rdm = 0d0
    do iorb=1,nOrb
     trace_1rdm = trace_1rdm+R(iorb,iorb) 
    enddo

    ! Adjust the chemical potential 

    if( abs(trace_1rdm-nO) > thrs_N ) & 
     call fix_chem_pot(nO,nOrb,nOrb_twice,nSCF,thrs_N,trace_1rdm,chem_pot,H_qsGWB,U_QP,R,eigVAL)

    ! DIIS extrapolation

    if(max_diis > 1 .and. nSCF>1) then

     write(*,*) ' Doing DIIS'

     F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:) + Sigc_ao_he(:,:) - chem_pot*S(:,:)
     H_qsGWB_ao(:,:)    = 0d0
     H_qsGWB_ao(1:nBas           ,1:nBas           ) =  F(1:nBas,1:nBas)
     H_qsGWB_ao(nBas+1:nBas_twice,nBas+1:nBas_twice) = -F(1:nBas,1:nBas)
     H_qsGWB_ao(1:nBas           ,nBas+1:nBas_twice) = Delta(1:nBas,1:nBas) + Sigc_ao_hh(1:nBas,1:nBas)
     H_qsGWB_ao(nBas+1:nBas_twice,1:nBas           ) = Delta(1:nBas,1:nBas) + Sigc_ao_hh(1:nBas,1:nBas)
     err_ao = matmul(H_qsGWB_ao,matmul(R_ao_old,S_ao)) - matmul(matmul(S_ao,R_ao_old),H_qsGWB_ao)

     n_diis = min(n_diis+1,max_diis)
     call DIIS_extrapolation(rcond,nBas_twice_Sq,nBas_twice_Sq,n_diis,err_diis,H_qsGWB_diis,err_ao,H_qsGWB_ao)

     H_qsGWB = matmul(transpose(X_ao),matmul(H_qsGWB_ao,X_ao))
     U_QP(:,:) = H_qsGWB(:,:)
     call diagonalize_matrix(nOrb_twice,U_QP,eigVAL)

     ! Build R and check trace
       
     trace_1rdm = 0d0 
     R(:,:)     = 0d0
     do iorb=1,nOrb
      R(:,:) = R(:,:) + matmul(U_QP(:,iorb:iorb),transpose(U_QP(:,iorb:iorb))) 
     enddo
     do iorb=1,nOrb
      trace_1rdm = trace_1rdm + R(iorb,iorb) 
     enddo

     ! Adjust the chemical potential 
     
     if( abs(trace_1rdm-nO) > thrs_N ) & 
      call fix_chem_pot(nO,nOrb,nOrb_twice,nSCF,thrs_N,trace_1rdm,chem_pot,H_qsGWB,U_QP,R,eigVAL)
   
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
    EqsGWB = ET + EV + EJ + EK + EL

    ! Check convergence

    F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:) + Sigc_ao_he(:,:) - chem_pot*S(:,:)
    H_qsGWB_ao(:,:)    = 0d0
    H_qsGWB_ao(1:nBas           ,1:nBas           ) =  F(1:nBas,1:nBas)
    H_qsGWB_ao(nBas+1:nBas_twice,nBas+1:nBas_twice) = -F(1:nBas,1:nBas)
    H_qsGWB_ao(1:nBas           ,nBas+1:nBas_twice) = Delta(1:nBas,1:nBas) + Sigc_ao_hh(1:nBas,1:nBas)
    H_qsGWB_ao(nBas+1:nBas_twice,1:nBas           ) = Delta(1:nBas,1:nBas) + Sigc_ao_hh(1:nBas,1:nBas)
    if(nSCF > 1) then
     err_ao = matmul(H_qsGWB_ao,matmul(R_ao_old,S_ao)) - matmul(matmul(S_ao,R_ao_old),H_qsGWB_ao)
     Conv  = maxval(abs(err_ao))
    endif

    ! Update R_old

    R_ao_old(:,:)    = 0d0
    R_ao_old(1:nBas           ,1:nBas           ) = 0.5d0*P(1:nBas,1:nBas)
    R_ao_old(nBas+1:nBas_twice,nBas+1:nBas_twice) = matmul(X(1:nBas,1:nOrb), transpose(X(1:nBas,1:nOrb))) &
                                                  -0.5d0*P(1:nBas,1:nBas)
    R_ao_old(1:nBas           ,nBas+1:nBas_twice) = Panom(1:nBas,1:nBas)
    R_ao_old(nBas+1:nBas_twice,1:nBas           ) = Panom(1:nBas,1:nBas)

    ! Compute new U_QP and electronic major Occ numbers

    allocate(c_ao(nBas_twice,nOrb_twice),U(nOrb,nOrb))
    U = matmul(transpose(c),matmul(F,c))
    call diagonalize_matrix(nOrb,U,eFock)
    c=matmul(c,U)
    c_ao(:,:)       = 0d0
    c_ao(1:nBas           ,1:nOrb           ) = c(1:nBas,1:nOrb)
    c_ao(nBas+1:nBas_twice,nOrb+1:nOrb_twice) = c(1:nBas,1:nOrb)
    U_QP = matmul(transpose(c_ao),matmul(H_qsGWB_ao,c_ao))
    call diagonalize_matrix(nOrb_twice,U_QP,eigVAL)
    Occ=0d0
    do iorb=1,nOrb
     do jorb=1,nOrb
      Occ(iorb)=Occ(iorb)+U_QP(jorb,iorb)**2d0
     enddo
    enddo
    do iorb=1,nOrb-1
     do jorb=iorb+1,nOrb
      if(abs(Occ(jorb))>abs(Occ(iorb))) then
       pivot_U_QP(:)= U_QP(:,iorb)
       U_QP(:,iorb) = U_QP(:,jorb)
       U_QP(:,jorb) = pivot_U_QP(:)
       pivot_Eqp    = eigVAL(iorb)
       eigVAL(iorb) = eigVAL(jorb)
       eigVAL(jorb) = pivot_Eqp
       pivot_Occ    = Occ(iorb)
       Occ(iorb)    = Occ(jorb)
       Occ(jorb)    = pivot_Occ
      endif
     enddo
    enddo
    do iorb=1,nOrb-1
     do jorb=iorb+1,nOrb
      if(abs(Occ(iorb))>0.5d0 .and. abs(Occ(jorb))>0.5d0 .and. eigVAL(jorb)<eigVAL(iorb)) then
       pivot_U_QP(:)= U_QP(:,iorb)
       U_QP(:,iorb) = U_QP(:,jorb)
       U_QP(:,jorb) = pivot_U_QP(:)
       pivot_Eqp    = eigVAL(iorb)
       eigVAL(iorb) = eigVAL(jorb)
       eigVAL(jorb) = pivot_Eqp
       pivot_Occ    = Occ(iorb)
       Occ(iorb)    = Occ(jorb)
       Occ(jorb)    = pivot_Occ
      endif
     enddo
    enddo
    do iorb=1,nOrb-1
     do jorb=iorb+1,nOrb
      if(abs(Occ(iorb))<=0.5d0 .and. abs(Occ(jorb))<=0.5d0 .and. eigVAL(jorb)>eigVAL(iorb)) then
       pivot_U_QP(:)= U_QP(:,iorb)
       U_QP(:,iorb) = U_QP(:,jorb)
       U_QP(:,jorb) = pivot_U_QP(:)
       pivot_Eqp    = eigVAL(iorb)
       eigVAL(iorb) = eigVAL(jorb)
       eigVAL(jorb) = pivot_Eqp
       pivot_Occ    = Occ(iorb)
       Occ(iorb)    = Occ(jorb)
       Occ(jorb)    = pivot_Occ
      endif
     enddo
    enddo
    deallocate(U,c_ao)
    eqsGWB_state(:) = eigVAL(:)


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

    deallocate(J,K,Sigc_ao_he,Sigc_ao_hh,H_qsGWB,R,eigVAL,err_diis,H_qsGWB_diis,eFock,Occ)
    deallocate(err_ao,S_ao,X_ao,R_ao_old,H_qsGWB_ao)
    deallocate(pivot_U_QP,Sigc_mo_he,Sigc_mo_hh)
    deallocate(vMAT)

    stop

  end if

! Compute dipole moments, occupation numbers, || Anomalous density||,
! organize the coefs c with natural orbitals (descending occ numbers), and
! also print the restart file

  allocate(U(nOrb,nOrb),c_no(nBas,nOrb))
  R(:,:)  = 0d0
  do iorb=1,nOrb
   R(:,:) = R(:,:) + matmul(U_QP(:,iorb:iorb),transpose(U_QP(:,iorb:iorb))) 
  enddo
  U(1:nOrb,1:nOrb) = R(1:nOrb,1:nOrb)
  call diagonalize_matrix(nOrb,U,Occ)
  c_no=matmul(c,U)
  Delta_HL=eqsGWB_state(nOrb+1)-eqsGWB_state(nOrb)
  N_anom = trace_matrix(nOrb,matmul(transpose(2d0*R(1:nOrb,nOrb+1:nOrb_twice)), &
              2d0*R(1:nOrb,nOrb+1:nOrb_twice)))
  call dipole_moment(nBas,P,nNuc,ZNuc,rNuc,dipole_int,dipole)
  call write_restart_qsGWB(nBas,nOrb,Occ,c_no,chem_pot) ! Warning: orders Occ and their c in descending order w.r.t. occupation numbers.
  call print_qsGWB(nBas,nOrb,nOrb_twice,nO,N_anom,Occ,eqsGWB_state,ENuc,ET,EV,EJ,EK,EL,EqsGWB,chem_pot,&
                   dipole,Delta_HL)
  deallocate(U,c_no)

! Testing zone

  if(dotest) then
 
    call dump_test_value('R','qsGWB energy',EqsGWB)
    call dump_test_value('R','Trace 1D',trace_1rdm)
    call dump_test_value('R','qsGWB dipole moment',norm2(dipole))

  end if

! Memory deallocation

  deallocate(J,K,Sigc_ao_he,Sigc_ao_hh,H_qsGWB,R,eigVAL,err_diis,H_qsGWB_diis,eFock,Occ)
  deallocate(err_ao,S_ao,X_ao,R_ao_old,H_qsGWB_ao)
  deallocate(pivot_U_QP,Sigc_mo_he,Sigc_mo_hh)
  deallocate(vMAT)

end subroutine 

