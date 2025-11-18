subroutine OORG0W0(dotest,doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,TDA_W,TDA,dBSE,dTDA,doppBSE,singlet,triplet, & 
                 linearize,eta,doSRG,nBas,nOrb,nC,nO,nV,nR,nS,                                                        &
                 maxIter,thresh,dRPA,mu,diagHess,                                                                     &
                 ENuc,ERHF,ERI_AO,ERI_MO,                                                                             &
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
  integer,intent(in)            :: mu,maxIter
  double precision,intent(in)   :: thresh
  logical,intent(in)            :: dRPA,diagHess
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
  logical                       :: OVRotOnly = .true.
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
  
  double precision,allocatable  :: ERI_AO_AS(:,:,:,:)
  double precision,allocatable  :: ERI_MO_AS(:,:,:,:)

  double precision,allocatable  :: eGWlin(:)
  double precision,allocatable  :: eGW(:)
 
  double precision              :: OOConv
  integer                       :: OOi
  double precision,allocatable  :: h(:,:)
  double precision,allocatable  :: rdm1(:,:)
  double precision,allocatable  :: rdm1_hf(:,:)
  double precision,allocatable  :: rdm1_rpa(:,:)
  double precision,allocatable  :: rdm2(:,:,:,:)
  double precision,allocatable  :: rdm2_hf(:,:,:,:)
  double precision,allocatable  :: rdm2_rpa(:,:,:,:)
  integer                       :: r,s,rs,p,q,pq
  integer                       :: i,jind,a,b,ia,jb,l
  integer                       :: ind
  integer                       :: N,O,V,Nsq
  double precision,allocatable  :: c(:,:)
  double precision,allocatable  :: F(:,:)
  double precision,allocatable  :: J(:,:),K(:,:)
  double precision              :: Emu, EOld
  double precision              :: EHF_rdm,EcRPA_rdm,EcRPA_HF

  double precision,external     :: trace_matrix
  double precision,external     :: Kronecker_delta
  double precision              :: trace_rdm2

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
  
   if(.not. dRPA) then
      write(*,*) 'xRPA instead of RPA'
      write(*,*)
      ! Antisymmetrize 2-electron integrals.
      allocate(ERI_AO_AS(nBas,nBas,nBas,nBas),ERI_MO_AS(nOrb,nOrb,nOrb,nOrb))
      ERI_AO_AS(:,:,:,:) = 0d0
      ERI_MO_AS(:,:,:,:) = 0d0
      do p=1,N
        do q=1,N
          do r=1,N
            do s=1,N
              ERI_AO_AS(p,q,r,s) = ERI_AO(p,q,r,s) - 0.5d0*ERI_AO(p,q,s,r)
            enddo
          enddo
        enddo
      enddo
   else
      allocate(ERI_AO_AS(1,1,1,1),ERI_MO_AS(1,1,1,1)) 
   endif

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
           rdm1(N,N),rdm2(N,N,N,N),rdm1_hf(N,N),rdm2_hf(N,N,N,N),rdm1_rpa(N,N),rdm2_rpa(N,N,N,N),&
           J(nBas,nBas),K(nBas,nBas),F(nOrb,nOrb))

! Initialize variables for OO  
  OOi           = 1
  OOConv        = 1d0
  c(:,:)        = cHF(:,:)
  rdm1(:,:)         = 0d0 
  rdm1_hf(:,:)      = 0d0 
  rdm1_rpa(:,:)     = 0d0 
  rdm2(:,:,:,:)     = 0d0
  rdm2_hf(:,:,:,:)  = 0d0
  rdm2_rpa(:,:,:,:) = 0d0
  rampl(:,:)    = 0d0
  lampl(:,:)    = 0d0
  rp(:)         = 0d0
  lp(:)         = 0d0
  t(:,:)        = 0d0
  lambda(:,:)   = 0d0
  Emu           = ERHF
  EcRPA_HF      = 0d0
  eGW(:)        = eHF(:)
  h(:,:)        = 0d0
  
  ! Transform integrals (afterwards this is done in orbital optimization)
  call AOtoMO(nBas,nOrb,c,Hc,h)
  call AOtoMO_ERI_RHF(nBas,nOrb,c,ERI_AO,ERI_MO) 
  if(.not. dRPA) call AOtoMO_ERI_RHF(nBas,nOrb,c,ERI_AO_AS,ERI_MO_AS)

  write(*,*) "Start orbital optimization loop..."

  do while (OOConv > thresh)
    
    EOld = Emu
    ! Build Fock
    PHF(:,:) = 2d0 * matmul(c(:,1:nO), transpose(c(:,1:nO))) 
    J(:,:) = 0d0
    call Hartree_matrix_AO_basis(nBas,PHF,ERI_AO,J)
    call exchange_matrix_AO_basis(nBas,PHF,ERI_AO,K)
    FHF(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:)
    call AOtoMO(nBas,nOrb,C,FHF,F)
    write(*,*) "Orbital optimization Iteration: ", OOi 
   

  !-------------------!
  ! Compute screening !
  !-------------------!
  
                   call OO_phRLR_A(isp_W,dRPA,nOrb,nC,nO,nV,nR,nS,1d0,F,ERI_MO,Aph)
    if(.not.TDA_W) call phRLR_B(isp_W,dRPA,nOrb,nC,nO,nV,nR,nS,1d0,ERI_MO,Bph)
    
    call phRLR(TDA_W,nS,Aph,Bph,EcRPA,Om,XpY,XmY)
    
    if(OOi == 1) EcRPA_HF = EcRPA

    if(print_W) call print_excitation_energies('phRPA@RHF','singlet',nS,Om)
  
    call RG0W0_rdm2_hf(O,V,N,nS,rdm2_hf)
    call RG0W0_rdm1_hf(O,V,N,nS,rdm1_hf)
    
    ! Useful quantities
    X = transpose(0.5*(XpY + XmY))
    Y = transpose(0.5*(XpY - XmY))
    call inverse_matrix(nS,X,X_inv)
    t = matmul(Y,X_inv)
    Xbar = - matmul(t,Y) + X
    call inverse_matrix(nS,Xbar,Xbar_inv)
    lambda = matmul(Y,Xbar_inv)

    ! Calculate rdm1
    call RG0W0_rdm1_rpa(O,V,N,nS,lampl,rampl,lp,rp,lambda,t,rdm1_rpa)
    rdm1 = rdm1_hf + rdm1_rpa
    ! Calculate rdm2
    call RG0W0_rdm2_rpa(O,V,N,nS,lampl,rampl,lp,rp,lambda,t,rdm1_hf,rdm1_rpa,rdm2_rpa)
    rdm2 = rdm2_hf + rdm2_rpa
    call energy_from_rdm(N,h,ERI_MO,rdm1_hf,rdm2_hf,EHF_rdm,.false.)
    if(dRPA) then
      call energy_from_rdm(N,F,ERI_MO,rdm1_rpa,rdm2_rpa,EcRPA_rdm,.false.)
    else
      call energy_from_rdm(N,F,ERI_MO_AS,rdm1_rpa,rdm2_rpa,EcRPA_rdm,.false.)
    endif
    Emu = EcRPA + EHF_rdm
    write(*,*) "ERHF = ", ERHF
    write(*,*) "ERHF + EcRPA@HF = ", ERHF + EcRPA_HF
    write(*,*) "EcRPA = ", EcRPA_rdm
    write(*,*) "EcRPA (no rdm) = ", EcRPA
    write(*,*) "ERPA = ", EcRPA_rdm + EHF_rdm
   
    call R_optimize_orbitals(diagHess,OVRotOnly,dRPA,nBas,nOrb,nV,nR,nC,nO,N,Nsq,O,V,ERI_AO,ERI_AO_AS,ERI_MO,ERI_MO_AS,Hc,h,F,&
                             rdm1_hf,rdm1_rpa,rdm2_hf,rdm2_rpa,c,OOConv)
    
    write(*,*) '----------------------------------------------------------'
    write(*,'(A10,I4,A30)') ' Iteration', OOi ,'for RG0W0 orbital optimization'
    write(*,*) '----------------------------------------------------------'
    write(*,'(A40,F16.10,A3)') ' Convergence of orbital gradient = ',OOConv,' au'
    write(*,'(A40,F16.10,A3)') ' Energy difference = ',Emu-EOld,' au'
    write(*,*) '----------------------------------------------------------'
    write(*,*)
    
 !   if (OOi==1) then
 !     OOConv = 0d0 ! remove only for debugging
 !   end if
    OOi = OOi + 1 
    if(OOi > maxIter) then
      write(*,*) "Orbital optimiaztion failed !!!"
      OOConv = 0d0
    endif
  end do
  cHF(:,:) = c(:,:)
  
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
 
!--------------!
! Dump results !
!--------------!

  call print_RG0W0(nOrb,nC,nO,nV,nR,eHF,ENuc,ERHF,SigC,Z,eGW,EcRPA,EcGM)

    eGW_out(:) = eGW(:)
  deallocate(rdm1,rdm2,c,Aph,Bph,SigC,Z,Om,XpY,XmY,rho,eGW,&
             eGWlin,X,X_inv,Y,Xbar,Xbar_inv,lambda,t,rampl,lampl,rp,lp,h,&
             F)

end subroutine
