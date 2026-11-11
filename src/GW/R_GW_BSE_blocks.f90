subroutine R_GW_BSE_screening(dotest,TDA_W,lin_GW,eta_GW,reg_GW,reg_PA,eta_2b, &
                              nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF,eGW,OmW,rhoW)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  include 'parameters.h'
  logical,intent(in) :: dotest,TDA_W,lin_GW,reg_GW,reg_PA
  integer,intent(in) :: nOrb,nC,nO,nV,nR,nS
  double precision,intent(in) :: eta_GW,eta_2b,ENuc,ERHF,ERI(nOrb,nOrb,nOrb,nOrb),dipole_int(nOrb,nOrb,ncart),eHF(nOrb)
  double precision,intent(out) :: eGW(nOrb),OmW(nS),rhoW(nOrb,nOrb,nS)
  double precision,allocatable :: A(:,:),B(:,:),XpY(:,:),XmY(:,:)
  double precision :: EcRPA,weight,screen_factor
  integer :: m

  if(nO<=nC.or.nOrb-nR<=nO.or.nS/=(nO-nC)*(nOrb-nR-nO)) error stop 'GW/BSE requires nonempty active spaces'
  if(eta_GW<0d0) error stop 'GW/BSE requires nonnegative broadening'
  if(reg_PA.and.eta_2b<=0d0) error stop 'GW/BSE requires a positive static screening flow'
! Build G0W0 and the HF screening modes shared by the static BSE kernels.

  call RG0W0(dotest,.false.,.false.,.false.,.false.,.false.,TDA_W,.false., &
              .false.,.false.,.false.,.true.,.true.,lin_GW,eta_GW,reg_GW,.false., &
              nOrb,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF,eGW)
  if(.not.all(ieee_is_finite(eGW))) error stop 'GW/BSE: nonfinite G0W0 energies'

! Reconstruct the same HF W0 used by G0W0, not a new W[eGW] iteration.

  allocate(A(nS,nS),B(nS,nS),XpY(nS,nS),XmY(nS,nS))
  A=0d0
  B=0d0
  call phRLR_A(1,.true.,nOrb,nC,nO,nV,nR,nS,1d0,eHF,ERI,A)
  if(.not.TDA_W) call phRLR_B(1,.true.,nOrb,nC,nO,nV,nR,nS,1d0,ERI,B)
  call phRLR(TDA_W,nS,A,B,EcRPA,OmW,XpY,XmY)
  if(.not.all(ieee_is_finite(OmW)) .or. any(OmW <= 0d0)) &
    error stop 'GW/BSE: unstable or zero-frequency W0 response'
  call RGW_excitation_density(nOrb,nC,nO,nR,nS,ERI,XpY,rhoW)
  deallocate(A,B,XpY,XmY)

! An SRG factor on each screening mode regularizes the common static kernel.
! Using sqrt(factor) on both half vertices keeps the BSE matrices and the
! external channel couplings consistent. reg_PA=F gives the native GW kernel.

  if(reg_PA) then
    do m=1,nS
      weight=2d0*eta_2b*OmW(m)**2
      if(weight < 1d-5) then
        screen_factor=weight*(1d0-weight/2d0+weight**2/6d0)
      else
        screen_factor=1d0-exp(-weight)
      end if
      rhoW(:,:,m)=sqrt(screen_factor)*rhoW(:,:,m)
    end do
  end if

end subroutine

subroutine R_GW_BSE_eh_blocks(ispin,TDAeh,nOrb,nC,nO,nV,nR,nS,eta_GW,eGW,ERI,OmW,rhoW,aa,bb)
  implicit none
  integer,intent(in) :: ispin,nOrb,nC,nO,nV,nR,nS
  logical,intent(in) :: TDAeh
  double precision,intent(in) :: eta_GW,eGW(nOrb),ERI(nOrb,nOrb,nOrb,nOrb),OmW(nS),rhoW(nOrb,nOrb,nS)
  double precision,intent(out) :: aa(nS,nS),bb(nS,nS)
  double precision,allocatable :: ka(:,:),kb(:,:)
  allocate(ka(nS,nS),kb(nS,nS))
    bb=0d0
    kb=0d0
    call phRLR_A(ispin,.false.,nOrb,nC,nO,nV,nR,nS,1d0,eGW,ERI,aa)
    call RGW_phBSE_static_kernel_A(eta_GW,nOrb,nC,nO,nV,nR,nS,1d0,ERI,OmW,rhoW,ka)
    if(.not.TDAeh) then
      call phRLR_B(ispin,.false.,nOrb,nC,nO,nV,nR,nS,1d0,ERI,bb)
      call RGW_phBSE_static_kernel_B(eta_GW,nOrb,nC,nO,nV,nR,nS,1d0,ERI,OmW,rhoW,kb)
    end if
    aa=aa+ka
    bb=bb+kb
end subroutine

subroutine R_GW_BSE_pp_blocks(ispin,TDApp,nOrb,nC,nO,nV,nR,nS,nOO,nVV,eta_GW,eGW,ERI,OmW,rhoW,bb,cc,dd)
  implicit none
  integer,intent(in) :: ispin,nOrb,nC,nO,nV,nR,nS,nOO,nVV
  logical,intent(in) :: TDApp
  double precision,intent(in) :: eta_GW,eGW(nOrb),ERI(nOrb,nOrb,nOrb,nOrb),OmW(nS),rhoW(nOrb,nOrb,nS)
  double precision,intent(out) :: bb(nVV,nOO),cc(nVV,nVV),dd(nOO,nOO)
  double precision,allocatable :: kb(:,:),kc(:,:),kd(:,:)
  allocate(kb(nVV,nOO),kc(nVV,nVV),kd(nOO,nOO))
    bb=0d0
    kb=0d0
    kc=0d0
    kd=0d0
    call ppRLR_C(ispin,nOrb,nC,nO,nV,nR,nVV,1d0,eGW,ERI,cc)
    call ppRLR_D(ispin,nOrb,nC,nO,nV,nR,nOO,1d0,eGW,ERI,dd)
    if(nVV > 0) call RGW_ppBSE_static_kernel_C(ispin,eta_GW,nOrb,nC,nO,nV,nR,nS,nVV, &
                                             1d0,ERI,OmW,rhoW,kc)
    if(nOO > 0) call RGW_ppBSE_static_kernel_D(ispin,eta_GW,nOrb,nC,nO,nV,nR,nS,nOO, &
                                             1d0,ERI,OmW,rhoW,kd)
    if(.not.TDApp .and. nOO*nVV > 0) then
      call ppRLR_B(ispin,nOrb,nC,nO,nV,nR,nOO,nVV,1d0,ERI,bb)
      call RGW_ppBSE_static_kernel_B(ispin,eta_GW,nOrb,nC,nO,nV,nR,nS,nOO,nVV, &
                                     1d0,ERI,OmW,rhoW,kb)
    end if
    cc=cc+kc
    dd=dd+kd
    bb=bb+kb
end subroutine
