subroutine R_GW_Parquet(dotest,TDA_W,lin_GW,eta_GW,reg_GW,TDAeh,TDApp,linearize, &
                        eta_1b,eta_2b,reg_PA,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF, &
                        ERI,dipole_int,eHF,eGW,ePA,ZPA)

! Frozen-input, native parquet calculation from G0W0 and static BSE@G0W0.
! W0 is built from the HF dRPA response, as in RG0W0/RGW_phBSE/RGW_ppBSE.
! Both BSE gaps and every internal self-energy energy use the frozen eGW.
! There are no parquet one-body/two-body iterations and no ADC completion.

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  include 'parameters.h'

  logical,intent(in)           :: dotest,TDA_W,lin_GW,reg_GW,TDAeh,TDApp,linearize,reg_PA
  integer,intent(in)           :: nBas,nOrb,nC,nO,nV,nR,nS
  double precision,intent(in)  :: eta_GW,eta_1b,eta_2b,ENuc,ERHF
  double precision,intent(in)  :: ERI(nOrb,nOrb,nOrb,nOrb),dipole_int(nOrb,nOrb,ncart),eHF(nOrb)
  double precision,intent(out) :: eGW(nOrb),ePA(nOrb),ZPA(nOrb)

  integer :: nOcc,nVir,nOOs,nVVs,nOOt,nVVt,p,q,r,s,m,it,ls
  integer,parameter :: max_root_it=128,max_backtrack=24
  double precision,parameter :: root_tol=1d-9
  double precision :: EcRPA,EcGM,Ec_eh(nspin),Ec_pp(nspin),screen_factor,weight
  double precision :: sigma,dsigma,w,f,step,trial,trial_sigma,trial_deriv,trial_f
  double precision,allocatable :: A(:,:),B(:,:),OmW(:),XpY(:,:),XmY(:,:),rhoW(:,:,:),W0(:,:,:,:)
  double precision,allocatable :: OmES(:),OmET(:),rhoES(:,:,:),rhoET(:,:,:)
  double precision,allocatable :: OmPS(:),OmHS(:),OmPT(:),OmHT(:)
  double precision,allocatable :: rhoPS(:,:,:),rhoHS(:,:,:),rhoPT(:,:,:),rhoHT(:,:,:)
  double precision,allocatable :: SigC(:),Z(:),eLin(:),omega(:),sig_work(:),z_work(:)

  nOcc=nO-nC
  nVir=nOrb-nR-nO
  if(nOcc <= 0 .or. nVir <= 0 .or. nS /= nOcc*nVir) &
    error stop 'R_GW_Parquet requires nonempty active occupied and virtual spaces'
  if(eta_GW < 0d0) error stop 'R_GW_Parquet: GW broadening must be nonnegative'
  if(reg_PA .and. (eta_1b <= 0d0 .or. eta_2b <= 0d0)) &
    error stop 'R_GW_Parquet: positive reg_1b and reg_2b are required for SRG'
  nOOs=nOcc*(nOcc+1)/2
  nVVs=nVir*(nVir+1)/2
  nOOt=nOcc*(nOcc-1)/2
  nVVt=nVir*(nVir-1)/2

  write(*,*)
  write(*,*) '***********************************************'
  write(*,*) '* Restricted single-shot GW-Parquet calculation *'
  write(*,*) '***********************************************'
  write(*,*) 'Frozen RHF orbitals; G0W0 energies; static eh/pp BSE@G0W0'
  write(*,*) 'No one-body or two-body parquet self-consistency'
  write(*,'(1X,A,L2)') 'TDA for W0 screening: ',TDA_W
  write(*,'(1X,A,L2)') 'TDA for ehBSE: ',TDAeh
  write(*,'(1X,A,L2)') 'TDA for ppBSE: ',TDApp
  write(*,'(1X,A,L2)') 'Linearize final Dyson equation about e_GW: ',linearize
  write(*,'(1X,A,L2)') 'Parquet SRG regularization: ',reg_PA
  if(reg_PA) write(*,'(1X,A,2ES14.5)') 'SRG flows (one-body, static screening): ',eta_1b,eta_2b

! G0W0 only. Its optional BSE/density-matrix/ACFDT work is deliberately off.
! The two BSE calculations below retain their response amplitudes for PA.
! All arrays passed here are MO-sized, including in a rank-reduced AO basis.

  call RG0W0(dotest,.false.,.false.,.false.,.false.,.false.,TDA_W,.false., &
              .false.,.false.,.false.,.true.,.true.,lin_GW,eta_GW,reg_GW,.false., &
              nOrb,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF,eGW)
  if(.not.all(ieee_is_finite(eGW))) error stop 'R_GW_Parquet: nonfinite G0W0 energies'

! Reconstruct the same HF W0 used by G0W0, not a new W[eGW] iteration.

  allocate(A(nS,nS),B(nS,nS),OmW(nS),XpY(nS,nS),XmY(nS,nS),rhoW(nOrb,nOrb,nS))
  A=0d0
  B=0d0
  call phRLR_A(1,.true.,nOrb,nC,nO,nV,nR,nS,1d0,eHF,ERI,A)
  if(.not.TDA_W) call phRLR_B(1,.true.,nOrb,nC,nO,nV,nR,nS,1d0,ERI,B)
  call phRLR(TDA_W,nS,A,B,EcRPA,OmW,XpY,XmY)
  if(.not.all(ieee_is_finite(OmW)) .or. any(OmW <= 0d0)) &
    error stop 'R_GW_Parquet: unstable or zero-frequency W0 response'
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

  allocate(W0(nOrb,nOrb,nOrb,nOrb))
  W0=ERI
  do s=nC+1,nOrb-nR
    do r=nC+1,nOrb-nR
      do q=nC+1,nOrb-nR
        do p=nC+1,nOrb-nR
          do m=1,nS
            W0(p,q,r,s)=W0(p,q,r,s)-4d0*rhoW(p,r,m)*rhoW(q,s,m) &
                         *OmW(m)/(OmW(m)**2+eta_GW**2)
          end do
        end do
      end do
    end do
  end do

  allocate(OmES(nS),OmET(nS),rhoES(nOrb,nOrb,nS),rhoET(nOrb,nOrb,nS))
  call eh_bse(1,OmES,rhoES,Ec_eh(1))
  call eh_bse(2,OmET,rhoET,Ec_eh(2))

  allocate(OmPS(nVVs),OmHS(nOOs),rhoPS(nOrb,nOrb,nVVs),rhoHS(nOrb,nOrb,nOOs))
  allocate(OmPT(nVVt),OmHT(nOOt),rhoPT(nOrb,nOrb,nVVt),rhoHT(nOrb,nOrb,nOOt))
  call pp_bse(1,nOOs,nVVs,OmPS,rhoPS,OmHS,rhoHS,Ec_pp(1))
  call pp_bse(2,nOOt,nVVt,OmPT,rhoPT,OmHT,rhoHT,Ec_pp(2))
  deallocate(W0,rhoW,OmW)

! The BSE poles and Gamma_BSE times eigenvectors are the spectral vertices
! needed by the native PA contraction. No four-index Phi feedback is made.

  allocate(SigC(nOrb),Z(nOrb),eLin(nOrb),omega(nOrb),sig_work(nOrb),z_work(nOrb))
  omega=eGW
  call R_Parquet_self_energy_alt_omega(reg_PA,.true.,0,eta_1b,omega, &
         nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt,eGW,ERI, &
         rhoES,OmES,rhoET,OmET,rhoPS,OmPS,rhoPT,OmPT,rhoHS,OmHS,rhoHT,OmHT,EcGM,SigC,Z)
  if(.not.all(ieee_is_finite(SigC))) error stop 'R_GW_Parquet: nonfinite native self-energy'
  if(.not.all(ieee_is_finite(Z))) error stop 'R_GW_Parquet: singular native self-energy derivative'

  eLin=eHF
  eLin(nC+1:nOrb-nR)=eGW(nC+1:nOrb-nR)+Z(nC+1:nOrb-nR) &
       *(eHF(nC+1:nOrb-nR)+SigC(nC+1:nOrb-nR)-eGW(nC+1:nOrb-nR))
  ePA=eLin
  ZPA=Z

  if(.not.linearize) then
    write(*,*) 'Root search for the frozen-input Dyson equation (no feedback)'
    do p=nC+1,nOrb-nR
      w=eLin(p)
      do it=1,max_root_it
        call evaluate(p,w,sigma,dsigma)
        f=w-eHF(p)-sigma
        if(abs(f) < root_tol) exit
        if(abs(1d0-dsigma) < 1d-12) error stop 'R_GW_Parquet: singular Dyson root derivative'
        step=f/(1d0-dsigma)
        step=sign(min(abs(step),0.5d0),step)
        do ls=1,max_backtrack
          trial=w-step
          call evaluate(p,trial,trial_sigma,trial_deriv)
          trial_f=trial-eHF(p)-trial_sigma
          if(abs(trial_f) < abs(f)) exit
          step=0.5d0*step
        end do
        if(ls > max_backtrack) error stop 'R_GW_Parquet: Dyson root search stalled'
        w=trial
      end do
      if(it > max_root_it) error stop 'R_GW_Parquet: frozen-input Dyson root did not converge'
      ePA(p)=w
      SigC(p)=sigma
      ZPA(p)=1d0/(1d0-dsigma)
    end do
  end if

  write(*,*)
  write(*,*) ' Native GW-Parquet self-energy and quasiparticle energies'
  if(linearize) then
    write(*,*) ' Sigma_PA and Z_PA are evaluated at e_GW; e_PA uses the HF reference.'
  else
    write(*,*) ' Sigma_PA and Z_PA are evaluated at the converged frozen-input Dyson roots.'
  end if
  write(*,'(1X,A5,5A17)') 'Orb.','e_HF (eV)','e_GW (eV)','Sigma_PA (eV)','Z_PA','e_PA (eV)'
  do p=nC+1,nOrb-nR
    write(*,'(1X,I5,5F17.9)') p,eHF(p)*HaToeV,eGW(p)*HaToeV,SigC(p)*HaToeV,ZPA(p),ePA(p)*HaToeV
  end do
  write(*,'(1X,A,F17.9,A)') 'GW-Parquet HOMO energy = ',maxval(ePA(nC+1:nO))*HaToeV,' eV'
  write(*,'(1X,A,F17.9,A)') 'GW-Parquet LUMO energy = ',minval(ePA(nO+1:nOrb-nR))*HaToeV,' eV'
  if(dotest) then
    call dump_test_value('R','GW-Parquet HOMO energy',ePA(nO))
    call dump_test_value('R','GW-Parquet LUMO energy',ePA(nO+1))
  end if

contains

  subroutine eh_bse(ispin,Om,rho,Ec)
    integer,intent(in) :: ispin
    double precision,intent(out) :: Om(nS),rho(nOrb,nOrb,nS),Ec
    double precision,allocatable :: aa(:,:),bb(:,:),ka(:,:),kb(:,:),xp(:,:),xm(:,:)
    allocate(aa(nS,nS),bb(nS,nS),ka(nS,nS),kb(nS,nS),xp(nS,nS),xm(nS,nS))
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
    call phRLR(TDAeh,nS,aa,bb,Ec,Om,xp,xm)
    if(.not.all(ieee_is_finite(Om)) .or. any(Om <= 0d0)) &
      error stop 'R_GW_Parquet: unstable ehBSE@G0W0 response'
    call R_GW_Parquet_eh_integrals(ispin,nOrb,nC,nO,nR,nS,ERI,W0,xp,xm,rho)
    if(ispin == 1) call print_excitation_energies('ehBSE@G0W0','singlet',nS,Om)
    if(ispin == 2) call print_excitation_energies('ehBSE@G0W0','triplet',nS,Om)
  end subroutine eh_bse

  subroutine pp_bse(ispin,nOO,nVV,Om1,rho1,Om2,rho2,Ec)
    integer,intent(in) :: ispin,nOO,nVV
    double precision,intent(out) :: Om1(nVV),rho1(nOrb,nOrb,nVV),Om2(nOO),rho2(nOrb,nOrb,nOO),Ec
    logical :: shifted_pp
    double precision,allocatable :: bb(:,:),cc(:,:),dd(:,:),kb(:,:),kc(:,:),kd(:,:)
    double precision,allocatable :: x1(:,:),y1(:,:),x2(:,:),y2(:,:)
    rho1=0d0
    rho2=0d0
    Ec=0d0
    if(nOO+nVV == 0) return
    allocate(bb(nVV,nOO),cc(nVV,nVV),dd(nOO,nOO),kb(nVV,nOO),kc(nVV,nVV),kd(nOO,nOO))
    allocate(x1(nVV,nVV),y1(nOO,nVV),x2(nVV,nOO),y2(nOO,nOO))
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
    call ppRLR(TDApp,nOO,nVV,bb,cc,dd,Om1,x1,y1,Om2,x2,y2,Ec)
    ! ppRLR optionally uses poles measured from twice the chemical potential.
    ! Native PA denominators require absolute pair energies, so undo that shift.
    inquire(file='use_chem_pot_in_ppRPA',exist=shifted_pp)
    if(shifted_pp) then
      Om1=Om1+eGW(nO)+eGW(nO+1)
      Om2=Om2+eGW(nO)+eGW(nO+1)
    end if
    if(.not.all(ieee_is_finite(Om1)) .or. .not.all(ieee_is_finite(Om2))) &
      error stop 'R_GW_Parquet: nonfinite ppBSE@G0W0 response'
    call R_GW_Parquet_pp_integrals(ispin,nOrb,nC,nO,nR,nOO,nVV,W0,x1,y1,x2,y2,rho1,rho2)
    if(ispin == 1) then
      if(nVV > 0) call print_excitation_energies('ppBSE@G0W0','2p singlet',nVV,Om1)
      if(nOO > 0) call print_excitation_energies('ppBSE@G0W0','2h singlet',nOO,Om2)
    else
      if(nVV > 0) call print_excitation_energies('ppBSE@G0W0','2p triplet',nVV,Om1)
      if(nOO > 0) call print_excitation_energies('ppBSE@G0W0','2h triplet',nOO,Om2)
    end if
  end subroutine pp_bse

  subroutine evaluate(orb,freq,sig,deriv)
    integer,intent(in) :: orb
    double precision,intent(in) :: freq
    double precision,intent(out) :: sig,deriv
    omega=eGW
    omega(orb)=freq
    call R_Parquet_self_energy_alt_omega(reg_PA,.true.,orb,eta_1b,omega, &
           nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt,eGW,ERI, &
           rhoES,OmES,rhoET,OmET,rhoPS,OmPS,rhoPT,OmPT,rhoHS,OmHS,rhoHT,OmHT,EcGM,sig_work,z_work)
    sig=sig_work(orb)
    deriv=1d0-1d0/z_work(orb)
    if(.not.ieee_is_finite(sig) .or. .not.ieee_is_finite(deriv)) &
      error stop 'R_GW_Parquet: nonfinite frozen-input Dyson evaluation'
  end subroutine evaluate

end subroutine R_GW_Parquet


subroutine R_GW_Parquet_eh_integrals(ispin,nOrb,nC,nO,nR,nS,ERI,W0,XpY,XmY,rho)

! Gamma_eh^S(pq;rs)=2 v(pq;rs)-W0(pq;sr), Gamma_eh^T=-W0(pq;sr).
! Contract the same static GW kernel used by the BSE, including its Y part.

  implicit none
  integer,intent(in) :: ispin,nOrb,nC,nO,nR,nS
  double precision,intent(in) :: ERI(nOrb,nOrb,nOrb,nOrb),W0(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in) :: XpY(nS,nS),XmY(nS,nS)
  double precision,intent(out) :: rho(nOrb,nOrb,nS)
  integer :: p,q,j,b,jb
  double precision :: direct
  double precision,allocatable :: gx(:,:,:),gy(:,:,:),x(:,:),y(:,:)
  allocate(gx(nOrb,nOrb,nS),gy(nOrb,nOrb,nS),x(nS,nS),y(nS,nS))
  gx=0d0
  gy=0d0
  direct=0d0
  if(ispin == 1) direct=2d0
  x=0.5d0*(XpY+XmY)
  y=0.5d0*(XpY-XmY)
  do q=nC+1,nOrb-nR
    do p=nC+1,nOrb-nR
      jb=0
      do j=nC+1,nO
        do b=nO+1,nOrb-nR
          jb=jb+1
          gx(p,q,jb)=direct*ERI(p,b,q,j)-W0(p,b,j,q)
          gy(p,q,jb)=direct*ERI(p,j,q,b)-W0(p,j,b,q)
        end do
      end do
    end do
  end do
  call dgemm('N','T',nOrb*nOrb,nS,nS,1d0,gx,nOrb*nOrb,x,nS,0d0,rho,nOrb*nOrb)
  call dgemm('N','T',nOrb*nOrb,nS,nS,1d0,gy,nOrb*nOrb,y,nS,1d0,rho,nOrb*nOrb)
end subroutine R_GW_Parquet_eh_integrals


subroutine R_GW_Parquet_pp_integrals(ispin,nOrb,nC,nO,nR,nOO,nVV,W0,X1,Y1,X2,Y2,rho1,rho2)

! Singlet/triplet symmetrized W0, with unique-pair normalization. This is
! the full-index counterpart of RGW_ppBSE_static_kernel_{B,C,D} plus bare v.

  implicit none
  integer,intent(in) :: ispin,nOrb,nC,nO,nR,nOO,nVV
  double precision,intent(in) :: W0(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in) :: X1(nVV,nVV),Y1(nOO,nVV),X2(nVV,nOO),Y2(nOO,nOO)
  double precision,intent(out) :: rho1(nOrb,nOrb,nVV),rho2(nOrb,nOrb,nOO)
  integer :: p,q,c,d,k,l,cd,kl,offset
  double precision :: exchange,val
  double precision,allocatable :: gv(:,:,:),go(:,:,:)
  allocate(gv(nOrb,nOrb,nVV),go(nOrb,nOrb,nOO))
  gv=0d0
  go=0d0
  rho1=0d0
  rho2=0d0
  offset=0
  exchange=1d0
  if(ispin == 2) then
    offset=1
    exchange=-1d0
  end if
  do q=nC+1,nOrb-nR
    do p=nC+1,nOrb-nR
      cd=0
      do c=nO+1,nOrb-nR
        do d=c+offset,nOrb-nR
          cd=cd+1
          val=W0(p,q,c,d)+exchange*W0(p,q,d,c)
          if(c == d) val=val/sqrt(2d0)
          gv(p,q,cd)=val
        end do
      end do
      kl=0
      do k=nC+1,nO
        do l=k+offset,nO
          kl=kl+1
          val=W0(p,q,k,l)+exchange*W0(p,q,l,k)
          if(k == l) val=val/sqrt(2d0)
          go(p,q,kl)=val
        end do
      end do
    end do
  end do
  if(nVV > 0) then
    call dgemm('N','N',nOrb*nOrb,nVV,nVV,1d0,gv,nOrb*nOrb,X1,max(1,nVV),0d0,rho1,nOrb*nOrb)
    if(nOO > 0) call dgemm('N','N',nOrb*nOrb,nVV,nOO,1d0,go,nOrb*nOrb,Y1,nOO,1d0,rho1,nOrb*nOrb)
  end if
  if(nOO > 0) then
    if(nVV > 0) call dgemm('N','N',nOrb*nOrb,nOO,nVV,1d0,gv,nOrb*nOrb,X2,nVV,0d0,rho2,nOrb*nOrb)
    call dgemm('N','N',nOrb*nOrb,nOO,nOO,1d0,go,nOrb*nOrb,Y2,nOO,1d0,rho2,nOrb*nOrb)
  end if
end subroutine R_GW_Parquet_pp_integrals
