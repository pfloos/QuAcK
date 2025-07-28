subroutine OORG0W0(dotest,doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,TDA_W,TDA,dBSE,dTDA,doppBSE,singlet,triplet, & 
                 linearize,eta,doSRG,nBas,nOrb,nC,nO,nV,nR,nS,mu,ENuc,ERHF,ERI_AO,ERI_MO,                             &
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
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta
  logical,intent(in)            :: doSRG

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  integer,intent(in)            :: mu
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(inout):: ERI_MO(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: dipole_int(nOrb,nOrb,ncart)
  double precision,intent(inout):: eHF(nOrb)
  double precision,intent(inout):: cHF(nBas,nOrb)
  double precision,intent(inout):: PHF(nBas,nBas)
  double precision,intent(inout):: FHF(nBas,nBas)
  double precision,intent(in)   :: Sovl(nBas,nBas)
  double precision,intent(in)   :: Tkin(nBas,nBas)
  double precision,intent(in)   :: Vpot(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: XHF(nBas,nOrb)

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

  double precision,allocatable  :: eGWlin(:)
  double precision,allocatable  :: eGW(:)
 
  double precision              :: OOConv
  integer                       :: OOi
  double precision,allocatable  :: h(:,:)
  double precision,allocatable  :: Kap(:,:)
  double precision,allocatable  :: ExpKap(:,:)
  double precision,allocatable  :: hess(:,:)
  double precision,allocatable  :: hessInv(:,:)
  double precision,allocatable  :: grad(:)
  double precision,allocatable  :: rdm1(:,:)
  double precision,allocatable  :: rdm2(:,:,:,:)
  integer                       :: r,s,rs,p,q,pq
  integer                       :: N,O,V,Nsq
  double precision,allocatable  :: c(:,:)
  double precision              :: Emu
  double precision              :: buff

  double precision,external     :: trace_matrix
! Output variables

  double precision,intent(out)  :: eGW_out(nOrb)
  
  V = nV - nR
  O = nO - nC
  N = V + O
  Nsq = N*N

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

  allocate(Aph(nS,nS),Bph(nS,nS),SigC(nOrb),Z(nOrb),Om(nS),XpY(nS,nS),XmY(nS,nS),rho(nOrb,nOrb,nS), & 
           eGW(nOrb),eGWlin(nOrb),X(nS,nS),X_inv(nS,nS),Y(nS,nS),Xbar(nS,nS),Xbar_inv(nS,nS),lambda(nS,nS),t(nS,nS),&
           rampl(nS,N),lampl(nS,N),rp(N),lp(N),h(N,N),c(nBas,nOrb),&
           rdm1(N,N),rdm2(N,N,N,N))

! Initialize variables for OO  
  OOi           = 1d0
  OOConv        = 1d0
  c(:,:)        = cHF(:,:)
  h(:,:)        = 0d0 
  rdm1(:,:)     = 0d0 
  rdm2(:,:,:,:) = 0d0
  rampl(:,:)    = 0d0
  lampl(:,:)    = 0d0
  rp(:)         = 0d0
  lp(:)         = 0d0
  t(:,:)        = 0d0
  lambda(:,:)   = 0d0
  eGW(:)        = eHF(:)

  write(*,*) "Start orbital optimization loop..."

  do while (OOConv > 1e-3)
  
    write(*,*) "Orbital optimiation Iteration: ", OOi
    
    h = matmul(transpose(c),matmul(Hc,c))
    call AOtoMO_ERI_RHF(nBas,nOrb,c,ERI_AO,ERI_MO)

  !-------------------!
  ! Compute screening !
  !-------------------!
  
                   call phRLR_A(isp_W,dRPA_W,nOrb,nC,nO,nV,nR,nS,1d0,eHF,ERI_MO,Aph)
    if(.not.TDA_W) call phRLR_B(isp_W,dRPA_W,nOrb,nC,nO,nV,nR,nS,1d0,ERI_MO,Bph)
  
    call phRLR(TDA_W,nS,Aph,Bph,EcRPA,Om,XpY,XmY)
    
    write(*,*) "EcRPA = ", EcRPA
    if(print_W) call print_excitation_energies('phRPA@RHF','singlet',nS,Om)
  
  !--------------------------!
  ! Compute spectral weights !
  !--------------------------!
  
    call RGW_excitation_density(nOrb,nC,nO,nR,nS,ERI_MO,XpY,rho)
  
  !------------------------!
  ! Compute GW self-energy !
  !------------------------!
  
    if(doSRG) then 
      call RGW_SRG_self_energy_diag(flow,nBas,nOrb,nC,nO,nV,nR,nS,eHF,Om,rho,EcGM,SigC,Z)
    else
      call RGW_self_energy_diag(eta,nBas,nOrb,nC,nO,nV,nR,nS,eHF,Om,rho,EcGM,SigC,Z)
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
    
      call RGW_QP_graph(doSRG,eta,flow,nOrb,nC,nO,nV,nR,nS,eHF,Om,rho,eGWlin,eHF,eGW,Z)
  
    end if
  
  ! Plot self-energy, renormalization factor, and spectral function
  
    if(plot_self) call RGW_plot_self_energy(nOrb,eta,nC,nO,nV,nR,nS,eHF,eGW,Om,rho)
  
  ! Cumulant expansion 
  
  ! call RGWC(dotest,eta,nOrb,nC,nO,nV,nR,nS,Om,rho,eHF,eHF,eGW,Z)
   
  !--------------!
  ! Dump results !
  !--------------!
  
    call print_RG0W0(nOrb,nC,nO,nV,nR,eHF,ENuc,ERHF,SigC,Z,eGW,EcRPA,EcGM)
  
    eGW_out(:) = eGW(:)
    
    ! Useful quantities to calculate rdms

    X = 0.5*(XpY + XmY)
    Y = 0.5*(XpY - XmY)
    call inverse_matrix(nS,X,X_inv)
    t = matmul(Y,X_inv)
    Xbar = - matmul(t,Y) + X
    call inverse_matrix(nS,Xbar,Xbar_inv)
    lambda = 0.5*matmul(Y,Xbar_inv)
    
    ! Calculate rdm1
    call RG0W0_rdm1(O,V,N,nS,lampl,rampl,lp,rp,lambda,t,rdm1)
    write(*,*) "Trace rdm1: ", trace_matrix(N,rdm1)

    ! Calculate rdm2
    call RG0W0_rdm2(O,V,N,nS,lampl,rampl,lp,rp,lambda,t,rdm2)
    write(*,*) "Trace rdm2: ", trace_matrix(Nsq,rdm2)
    
    call RGW_energy_from_rdm(N,h,ERI_MO,rdm1,rdm2,Emu)
    write(*,*) "Emu = ", Emu
    
    !--------------------------!
    ! Compute orbital gradient !
    !--------------------------!
 
    allocate(grad(Nsq))
    grad(:) = 0d0
    call pCCD_orbital_gradient(O,V,N,Nsq,h,ERI_MO,rdm1,rdm2,grad)
   
    ! Check convergence of orbital optimization
    OOConv = maxval(grad)
    write(*,*) '-----------------------------------------------------------'
    write(*,'(A10,I4,A30)') ' Iteration',OOi,'for G0W0 orbital optimization'
    write(*,*) '-----------------------------------------------------------'
    write(*,'(A40,F16.10,A3)') ' Convergence of orbital gradient = ',OOConv,' au'
    write(*,*) '----------------------------------------------------------'
    write(*,*)

    !-------------------------!
    ! Compute orbital Hessian !
    !-------------------------!
 
    allocate(hess(Nsq,Nsq))
    hess(:,:) = 0d0 
    call pCCD_orbital_hessian(O,V,N,Nsq,h,ERI_MO,rdm1,rdm2,hess)
    
!    write(*,*) "Hessian"
!    call matout(Nsq,Nsq,hess)
 
    allocate(hessInv(Nsq,Nsq))
 
    call inverse_matrix(Nsq,hess,hessInv)
    
!    write(*,*) "Inv Hessian"
!    call matout(Nsq,Nsq,hessInv)
    
    deallocate(hess)
 
    allocate(Kap(N,N))
 
    Kap(:,:) = 0d0
 
    pq = 0
    do p=1,nOrb
      do q=1,nOrb
 
        pq = pq + 1
 
        rs = 0
        do r=1,nOrb
          do s=1,nOrb
 
            rs = rs + 1
 
              Kap(p,q) = Kap(p,q) - hessInv(pq,rs)*grad(rs)
 
          end do
        end do
 
      end do
    end do
 
    deallocate(hessInv,grad)

!    write(*,*) 'kappa'
!    call matout(nOrb,nOrb,Kap)
!    write(*,*)
! 
    allocate(ExpKap(N,N))
    call matrix_exponential(nOrb,Kap,ExpKap)
    deallocate(Kap)
! 
!    write(*,*) 'e^kappa'
!    call matout(nOrb,nOrb,ExpKap)
!    write(*,*)
! 
!    write(*,*) 'Old orbitals'
!    call matout(nBas,nOrb,c)
!    write(*,*)
! 
    c = matmul(c,ExpKap)
    deallocate(ExpKap)
! 
!    write(*,*) 'Rotated orbitals'
!    call matout(nBas,nOrb,c)
!    write(*,*)

    ! Compute all quantities new for rotated basis


 !   if (OOi==3) then
 !     OOConv = 0d0 ! remove only for debugging
 !   end if

    OOi = OOi + 1 
  end do
  cHF(:,:) = c(:,:)
  
  deallocate(rdm1,rdm2,c,Aph,Bph,SigC,Z,Om,XpY,XmY,rho,eGW,&
             eGWlin,X,X_inv,Y,Xbar,Xbar_inv,lambda,t,rampl,lampl,rp,lp,h)

end subroutine
