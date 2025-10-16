subroutine OOGG0W0(dotest,doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,TDA_W,TDA,dBSE,dTDA,doppBSE,& 
                linearize,eta,doSRG,do_linDM,nBas,nBas2,nC,nO,nV,nR,nS,mu,ENuc,EGHF,ERI_AO,ERI_MO,                             &
                 dipole_int,eHF,cHF,Sovl,XHF,Tkin,Vpot,Hc,PHF,FHF,eGW_out)

! Perform optimized orbital G0W0 calculation (optimal for excitation mu)

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: doXBS
  logical,intent(in)            :: dophBSE
  logical,intent(in)            :: dophBSE2
  logical,intent(in)            :: doppBSE
  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta
  logical,intent(in)            :: doSRG
  logical,intent(in)            :: do_linDM

  integer,intent(in)            :: nBas,nBas2
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  integer,intent(in)            :: mu
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EGHF
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(inout):: ERI_MO(nBas2,nBas2,nBas2,nBas2)
  double precision,intent(in)   :: dipole_int(nBas2,nBas2,ncart)
  double precision,intent(inout):: eHF(nBas2)
  double precision,intent(inout):: cHF(nBas2,nBas2)
  double precision,intent(inout):: PHF(nBas2,nBas2)
  double precision,intent(inout):: FHF(nBas2,nBas2)
  double precision,intent(in)   :: Sovl(nBas2,nBas)
  double precision,intent(in)   :: Tkin(nBas2,nBas2)
  double precision,intent(in)   :: Vpot(nBas2,nBas2)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: XHF(nBas2,nBas2)

! Local variables

  logical                       :: print_W   = .false.
  logical                       :: plot_self = .false.
  logical                       :: dRPA_W
  integer                       :: isp_W
  double precision              :: flow
  double precision              :: EcRPA
  double precision              :: EcBSE(nspin)
  double precision              :: EcGM
  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: SigC(:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)
  double precision,allocatable  :: X(:,:)
  double precision,allocatable  :: X_inv(:,:)
  double precision,allocatable  :: Xbar(:,:)
  double precision,allocatable  :: Xbar_inv(:,:)
  double precision,allocatable  :: Y(:,:)
  double precision,allocatable  :: lambda(:,:)
  double precision,allocatable  :: t(:,:)
  double precision,allocatable  :: rho(:,:,:)
  double precision,allocatable  :: rampl(:,:)
  double precision,allocatable  :: lampl(:,:)
  double precision,allocatable  :: rp(:)
  double precision,allocatable  :: lp(:)
  double precision,allocatable  :: Ca(:,:),Cb(:,:)
  double precision,allocatable  :: ERI_tmp(:,:,:,:)

  double precision,allocatable  :: eGWlin(:)
  double precision,allocatable  :: eGW(:)
 
  double precision              :: OOConv
  double precision              :: thresh = 1.0e-8
  integer                       :: OOi
  double precision,allocatable  :: h(:,:)
  double precision,allocatable  :: Kap(:,:)
  double precision,allocatable  :: ExpKap(:,:)
  double precision,allocatable  :: hess(:,:)
  double precision,allocatable  :: hessInv(:,:)
  double precision,allocatable  :: rdm1(:,:)
  double precision,allocatable  :: rdm1_hf(:,:)
  double precision,allocatable  :: rdm1_rpa(:,:)
  double precision,allocatable  :: rdm2(:,:,:,:)
  double precision,allocatable  :: rdm2_rpa(:,:,:,:)
  double precision,allocatable  :: rdm2_hf(:,:,:,:)
  integer                       :: r,s,rs,p,q,pq
  integer                       :: N,O,V,Nsq
  double precision,allocatable  :: c(:,:)
  double precision,allocatable  :: Fp(:,:)
  double precision,allocatable  :: J(:,:),K(:,:)
  double precision              :: Emu, EOld
  
  double precision              :: rdm2_trace
  double precision,external     :: trace_matrix
  double precision              :: Eblub1
  double precision              :: Eblub2
  double precision              :: EHF_rdm
  double precision              :: ERPA_rdm

! Output variables

  double precision,intent(out)  :: eGW_out(nBas2)
  
  V = nV - nR
  O = nO - nC
  N = V + O
  Nsq = N*N
  write(*,*) "cHF in the begin of OOGG0W0"
  call matout(nBas2,nBas2,cHF)

! Output variables

  write(*,*) "Under developpement, not working yet !!! QuAcK QuAcK !!!"

! Hello world

  write(*,*)
  write(*,*)'*********************************************************'
  write(*,*)'* Restricted G0W0 Calculation with orbital optimization *'
  write(*,*)'*********************************************************'
  write(*,*)
  write(*,*) "Targeted exited state mu: ", mu

! Spin manifold and TDA for dynamical screening

  isp_W = 1
  dRPA_W = .true.
 
   if(TDA_W) then 
     write(*,*) 'Tamm-Dancoff approximation for dynamical screening!'
     write(*,*)
   end if
 
 ! SRG regularization
  flow = 500d0

  if(doSRG) then

    write(*,*) '*** SRG regularized G0W0 scheme ***'
    write(*,*)

  end if
  
! Memory allocation

  allocate(Aph(nS,nS),Bph(nS,nS),SigC(nBas2),Z(nBas2),Om(nS),XpY(nS,nS),XmY(nS,nS),rho(nBas2,nBas2,nS), & 
           eGW(nBas2),eGWlin(nBas2),X(nS,nS),X_inv(nS,nS),Y(nS,nS),Xbar(nS,nS),Xbar_inv(nS,nS),lambda(nS,nS),t(nS,nS),&
           rampl(nS,N),lampl(nS,N),rp(N),lp(N),h(N,N),c(nBas2,nBas2),&
           rdm1(N,N),rdm2(N,N,N,N),rdm1_hf(N,N),rdm2_hf(N,N,N,N),rdm1_rpa(N,N),rdm2_rpa(N,N,N,N),&
           J(nBas2,nBas2),K(nBas2,nBas2))

! Initialize variables for OO  
  OOi               = 1d0
  OOConv            = 1d0
  c(:,:)            = cHF(:,:)
  rdm1(:,:)         = 0d0 
  rdm1_hf(:,:)      = 0d0 
  rdm1_rpa(:,:)     = 0d0 
  rdm2(:,:,:,:)     = 0d0
  rdm2_hf(:,:,:,:)  = 0d0
  rdm2_rpa(:,:,:,:) = 0d0
  rampl(:,:)        = 0d0
  lampl(:,:)        = 0d0
  rp(:)             = 0d0
  lp(:)             = 0d0
  t(:,:)            = 0d0
  lambda(:,:)       = 0d0
  Emu               = 0d0
  eGW(:)            = eHF(:)
  h(:,:)            = 0d0

  write(*,*) "Start orbital optimization loop..."

  do while (OOConv > thresh)

    write(*,*) "Orbital optimiation Iteration: ", OOi 
  
   !----------------------------------!
   ! AO to MO integral transformation !
   !----------------------------------!
   
    write(*,*)
    write(*,*) 'AO to MO transformation... Please be patient'
    write(*,*)

    h(:,:)                       = 0d0
    h(     1:nBas ,     1:nBas ) = Hc(1:nBas,1:nBas)
    h(nBas+1:nBas2,nBas+1:nBas2) = Hc(1:nBas,1:nBas)
    h(:,:)                       = matmul(transpose(c),matmul(h,c))
    
    write(*,*) "Core Hamiltonian"
    call matout(nBas2,nBas2,h) 
    
    allocate(Ca(nBas,nBas2),Cb(nBas,nBas2),ERI_tmp(nBas2,nBas2,nBas2,nBas2))
   
    Ca(:,:) = c(1:nBas,1:nBas2)
    Cb(:,:) = c(nBas+1:nBas2,1:nBas2)
   
    ! 4-index transform 
   
    call AOtoMO_ERI_GHF(nBas,nBas2,Ca,Ca,ERI_AO,ERI_tmp)
    ERI_MO(:,:,:,:) = ERI_tmp(:,:,:,:)
   
    call AOtoMO_ERI_GHF(nBas,nBas2,Ca,Cb,ERI_AO,ERI_tmp)
    ERI_MO(:,:,:,:) = ERI_MO(:,:,:,:) + ERI_tmp(:,:,:,:)
   
    call AOtoMO_ERI_GHF(nBas,nBas2,Cb,Ca,ERI_AO,ERI_tmp)
    ERI_MO(:,:,:,:) = ERI_MO(:,:,:,:) + ERI_tmp(:,:,:,:)
   
    call AOtoMO_ERI_GHF(nBas,nBas2,Cb,Cb,ERI_AO,ERI_tmp)
    ERI_MO(:,:,:,:) = ERI_MO(:,:,:,:) + ERI_tmp(:,:,:,:)
    deallocate(Ca,Cb,ERI_tmp)

  !-------------------!
  ! Compute screening !
  !-------------------!
                   call phGLR_A(dRPA_W,nBas2,nC,nO,nV,nR,nS,1d0,eHF,ERI_MO,Aph)
    if(.not.TDA_W) call phGLR_B(dRPA_W,nBas2,nC,nO,nV,nR,nS,1d0,ERI_MO,Bph)
    call phGLR(TDA_W,nS,Aph,Bph,EcRPA,Om,XpY,XmY)
    
    write(*,*) "EcRPA = ", EcRPA
    if(print_W) call print_excitation_energies('phRPA@GHF','generalized',nS,Om)
  
  !--------------------------!
  ! Compute spectral weights !
  !--------------------------!
  
    call GGW_excitation_density(nBas2,nC,nO,nR,nS,ERI_MO,XpY,rho)
  
  !------------------------!
  ! Compute GW self-energy !
  !------------------------!
  
    if(doSRG) then 
      call GGW_SRG_self_energy_diag(flow,nBas2,nC,nO,nV,nR,nS,eHF,Om,rho,EcGM,SigC,Z)
    else
      call GGW_self_energy_diag(eta,nBas2,nC,nO,nV,nR,nS,eHF,Om,rho,EcGM,SigC,Z,ERI_MO)
    end if
  
  !-----------------------------------!
  ! Solve the quasi-particle equation !
  !-----------------------------------!
  
    ! Linearized or graphical solution?
    eGWlin(:) = eHF(:) + Z(:)*SigC(:)
  
    if(linearize) then 
   
      write(*,*) ' *** Quasiparticle energies obtained by linearization *** '
      write(*,*)
  
      eGW(:) = eGWlin(:)
  
    else 
  
      write(*,*) ' *** Quasiparticle energies obtained by root search *** '
      write(*,*)
    
      call GGW_QP_graph(doSRG,eta,flow,nBas2,nC,nO,nV,nR,nS,eHF,Om,rho,eGWlin,eHF,eGW,Z)
  
    end if
  
  ! Compute the RPA correlation energy

                 call phGLR_A(dRPA_W,nBas2,nC,nO,nV,nR,nS,1d0,eGW,ERI_MO,Aph)
  if(.not.TDA_W) call phGLR_B(dRPA_W,nBas2,nC,nO,nV,nR,nS,1d0,ERI_MO,Bph)
  call phGLR(TDA_W,nS,Aph,Bph,EcRPA,Om,XpY,XmY)
  write(*,*) "EcRPA = ", EcRPA
  
  !--------------!
  ! Dump results !
  !--------------!
  
    call print_GG0W0(nBas2,nC,nO,nV,nR,eHF,ENuc,EGHF,SigC,Z,eGW,EcRPA,EcGM)
  
    eGW_out(:) = eGW(:)
    
    ! Useful quantities to calculate rdms

    X = transpose(0.5*(XpY + XmY))
    Y = transpose(0.5*(XpY - XmY))
    call inverse_matrix(nS,X,X_inv)
    t = matmul(Y,X_inv)
    Xbar = - matmul(t,Y) + X
    call inverse_matrix(nS,Xbar,Xbar_inv)
    lambda = 0.5*matmul(Y,Xbar_inv)
    
    ! Calculate rdm1
    call GG0W0_rdm1_hf(O,V,N,nS,lampl,rampl,lp,rp,lambda,t,rdm1_hf)
    call GG0W0_rdm1_rpa(O,V,N,nS,lampl,rampl,lp,rp,lambda,t,rdm1_rpa)
    rdm1 = rdm1_hf  + rdm1_rpa
    write(*,*) "Trace rdm1: ", trace_matrix(N,rdm1)
    call matout(N,N,rdm1)
    ! Calculate rdm2
    call GG0W0_rdm2_hf(O,V,N,nS,lampl,rampl,lp,rp,lambda,t,rdm2_hf)
    call GG0W0_rdm2_rpa(O,V,N,nS,lampl,rampl,lp,rp,lambda,t,rdm2_rpa)
    rdm2 = rdm2_hf  + rdm2_rpa
    do p=1,nBas2
      do q=1,nBas2
        rdm2_trace = rdm2_trace + rdm2(p,q,p,q)
      end do
    end do
    write(*,*) "Trace rdm2: ", rdm2_trace
    call matout(Nsq,Nsq,rdm2)
    EOld = Emu 
    call energy_from_rdm(N,h,ERI_MO,rdm1,rdm2,Emu)
    call energy_from_rdm(N,h,ERI_MO,rdm1_hf,rdm2_hf,EHF_rdm)
    call energy_from_rdm(N,h,ERI_MO,rdm1_rpa,rdm2_rpa,ERPA_rdm)
    write(*,*) "EGHF = ", EGHF
    write(*,*) "EHF from rdm = " , EHF_rdm
    write(*,*) "EcRPA from rdm = " , ERPA_rdm
    write(*,*) "EcRPA = ", EcRPA
    write(*,*) "E elec", Emu
    write(*,*) "ENuc", ENuc
    
    call G_optimize_orbitals(nBas,nBas2,nV,nR,nC,nO,N,Nsq,O,V,ERI_AO,ERI_MO,h,rdm1,rdm2,c,OOConv)
   
    write(*,*) '----------------------------------------------------------'
    write(*,'(A10,I4,A30)') ' Iteration', OOi ,'for GG0W0 orbital optimization'
    write(*,*) '----------------------------------------------------------'
    write(*,'(A40,F16.10,A3)') ' Convergence of orbital gradient = ',OOConv,' au'
    write(*,'(A40,F16.10,A3)') ' Energy difference = ',Emu-EOld,' au'
    write(*,*) '----------------------------------------------------------'
    write(*,*)
    
    if (OOi==1) then
      OOConv = 0d0 ! remove only for debugging
    end if

    OOi = OOi + 1 
  end do
  cHF(:,:) = c(:,:)
  
  deallocate(rdm1,rdm2,c,Aph,Bph,SigC,Z,Om,XpY,XmY,rho,eGW,&
             eGWlin,X,X_inv,Y,Xbar,Xbar_inv,lambda,t,rampl,lampl,rp,lp,h)

end subroutine
