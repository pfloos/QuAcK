subroutine G_TDA_Faddeev_self_energy(adc2,reg_PA,eta,nOrb,nC,nO,nR,nS,nOO,nVV, &
                                    eHF,eQP,ERI,Weh,Whh,Wee,ePA,ZPA, &
    dressed,Teh,Thh,Tee,static3)
! Complete-TDA physical reduction of faddeev_parquet_notes.tex.
! Weh=Aeh-diag(ea-ei), Whh=Dpp+diag(ei+ej), Wee=Cpp-diag(ea+eb).
! Endpoints: ADC(2), ADC(3), or static one-sided dressed sources.
  implicit none
  logical,intent(in) :: adc2,reg_PA
  integer,intent(in) :: nOrb,nC,nO,nR,nS,nOO,nVV
  double precision,intent(in) :: eta,eHF(nOrb),eQP(nOrb),ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in) :: Weh(nS,nS),Whh(nOO,nOO),Wee(nVV,nVV)
  double precision,intent(out) :: ePA(nOrb),ZPA(nOrb)
  integer :: nH
  logical,intent(in) :: dressed,static3
  double precision,intent(in) :: Teh(nS,nS),Thh(nOO,nVV),Tee(nOO,nVV)
  double precision,allocatable :: H(:,:),sigma(:,:)
  if(nC<0.or.nR<0.or.nO<=nC.or.nOrb-nR<=nO) error stop 'G TDA Faddeev: invalid active space'
  nH=nOrb+nOO*(nOrb-nR-nO)+nVV*(nO-nC)
  allocate(H(nH,nH))
  call G_TDA_Faddeev_build(adc2,reg_PA,eta,0,nOrb,nC,nO,nR,nS,nOO,nVV,nH,eHF,eQP,ERI,Weh,Whh,Wee,H, &
    dressed,Teh,Thh,Tee)
  if(static3) then
    allocate(sigma(nOrb,nOrb))
    call TDA_Faddeev_static3(nOrb,nOrb,nOrb,nC,nO,nR,eHF,ERI,sigma)
    H(1:nOrb,1:nOrb)=H(1:nOrb,1:nOrb)+sigma
  end if
  call G_TDA_Faddeev_spectrum(dressed,adc2,nOrb,nC,nR,nH,eHF,H,ePA,ZPA)
end subroutine

subroutine TDA_Faddeev_strip_pair(nOrb,first,last,equal_pair,removal,nPair,e,shift,K)
  implicit none
  integer,intent(in) :: nOrb,first,last,nPair
  logical,intent(in) :: equal_pair,removal
  double precision,intent(in) :: e(nOrb),shift
  double precision,intent(inout) :: K(nPair,nPair)
  integer :: i,j,r,offset
  double precision :: sign_free
  sign_free=1d0
  if(removal) sign_free=-1d0
  offset=1
  if(equal_pair) offset=0
  r=0
  do i=first,last
    do j=i+offset,last
      r=r+1
      K(r,r)=K(r,r)-sign_free*(e(i)+e(j)-shift)
    end do
  end do
  if(r/=nPair) error stop 'TDA Faddeev: inconsistent pair stripping dimension'
end subroutine

subroutine G_TDA_Faddeev_build(adc2,reg_PA,eta,spin_block,nOrb,nC,nO,nR,nS,nOO,nVV,nH, &
                              eHF,eQP,ERI,Weh,Whh,Wee,H, &
    dressed,Teh,Thh,Tee)
! spin_block=0: general spin orbitals; 1: alpha charged block, with
! interleaved alpha/beta orbitals (used by the restricted wrapper).
! Kernel signs follow the actual QuAcK TDA pair Hamiltonians, avoiding
! four-index crossing ambiguities. Pair bases must be complete.
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  logical,intent(in) :: adc2,reg_PA
  integer,intent(in) :: spin_block,nOrb,nC,nO,nR,nS,nOO,nVV,nH
  double precision,intent(in) :: eta,eHF(nOrb),eQP(nOrb),ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in) :: Weh(nS,nS),Whh(nOO,nOO),Wee(nVV,nVV)
  double precision,intent(out) :: H(nH,nH)
  integer :: noact,nvact,nphys,nm,np,p,q,i,j,k,l,a,b,c,d,r,s,ij,kl,ab,cd
  logical,intent(in) :: dressed
  double precision,intent(in) :: Teh(nS,nS),Thh(nOO,nVV),Tee(nOO,nVV)
  integer,allocatable :: ph(:,:),hh(:,:),ee(:,:),minus(:,:),plus(:,:),phys(:)
  double precision :: x,asym
  double precision,allocatable :: Um(:,:),Up(:,:)

  noact=nO-nC
  nvact=nOrb-nR-nO
  if(noact<1.or.nvact<1) error stop 'TDA Faddeev: empty occupied or virtual space'
  if(nS/=noact*nvact.or.nOO/=noact*(noact-1)/2.or.nVV/=nvact*(nvact-1)/2) &
    error stop 'TDA Faddeev: inconsistent active pair dimensions'
  if(spin_block/=0.and.spin_block/=1) error stop 'TDA Faddeev: invalid spin block'
  if(spin_block==1) then
    if(mod(nOrb,2)/=0.or.mod(nC,2)/=0.or.mod(nO,2)/=0.or.mod(nR,2)/=0) &
      error stop 'TDA Faddeev alpha block requires paired spin orbitals'
  end if
  if(reg_PA.and.(eta<=0d0.or..not.ieee_is_finite(eta))) error stop 'TDA Faddeev: invalid SRG flow'
  nphys=nOrb
  if(spin_block==1) nphys=nOrb/2
  allocate(ph(nOrb,nOrb),hh(nOrb,nOrb),ee(nOrb,nOrb), &
           minus(3,nOO*nvact),plus(3,nVV*noact),phys(nphys))
  ph=0; hh=0; ee=0
  r=0
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      r=r+1; ph(i,a)=r
    end do
  end do
  r=0; nm=0
  do i=nC+1,nO
    do j=i+1,nO
      r=r+1; hh(i,j)=r
      do a=nO+1,nOrb-nR
        if(spin_block==1) then
          if(spin(i)+spin(j)-spin(a)/=1) cycle
        end if
        nm=nm+1; minus(:,nm)=[i,j,a]
      end do
    end do
  end do
  r=0
  do a=nO+1,nOrb-nR
    do b=a+1,nOrb-nR
      r=r+1; ee(a,b)=r
    end do
  end do
  np=0
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do b=a+1,nOrb-nR
        if(spin_block==1) then
          if(spin(a)+spin(b)-spin(i)/=1) cycle
        end if
        np=np+1; plus(:,np)=[i,a,b]
      end do
    end do
  end do
  if(nH/=nphys+nm+np) error stop 'TDA Faddeev: incorrect charged dimension'
  do p=1,nphys
    phys(p)=p
    if(spin_block==1) phys(p)=2*p-1
  end do
  allocate(Um(nm,nphys),Up(np,nphys))
  call G_TDA_Faddeev_vertices(adc2,reg_PA,eta,nOrb,nOrb,nC,nO,nR,nphys,nm,np, &
                              phys,minus,plus,eQP,ERI,Um,Up, &
    dressed,Teh,Thh,Tee)
  H=0d0
  do p=1,nphys
    H(p,p)=eHF(phys(p))
  end do

! Removal block with the selected external couplings.
  do r=1,nm
    i=minus(1,r); j=minus(2,r); a=minus(3,r); ij=hh(i,j)
    H(nphys+r,nphys+r)=eQP(i)+eQP(j)-eQP(a)
    do s=1,nm
      k=minus(1,s); l=minus(2,s); c=minus(3,s); kl=hh(k,l)
      x=0d0
      if(a==c) x=x-Whh(ij,kl)
      if(i==k) x=x-Weh(ph(j,a),ph(l,c))
      if(j==l) x=x-Weh(ph(i,a),ph(k,c))
      if(i==l) x=x+Weh(ph(j,a),ph(k,c))
      if(j==k) x=x+Weh(ph(i,a),ph(l,c))
      H(nphys+r,nphys+s)=H(nphys+r,nphys+s)+x
    end do
    do q=1,nphys
      H(q,nphys+r)=Um(r,q)
      H(nphys+r,q)=Um(r,q)
    end do
  end do

! Addition block with the selected external couplings.
  do r=1,np
    i=plus(1,r); a=plus(2,r); b=plus(3,r); ab=ee(a,b)
    H(nphys+nm+r,nphys+nm+r)=eQP(a)+eQP(b)-eQP(i)
    do s=1,np
      k=plus(1,s); c=plus(2,s); d=plus(3,s); cd=ee(c,d)
      x=0d0
      if(i==k) x=x+Wee(ab,cd)
      if(a==c) x=x+Weh(ph(i,b),ph(k,d))
      if(b==d) x=x+Weh(ph(i,a),ph(k,c))
      if(a==d) x=x-Weh(ph(i,b),ph(k,c))
      if(b==c) x=x-Weh(ph(i,a),ph(k,d))
      H(nphys+nm+r,nphys+nm+s)=H(nphys+nm+r,nphys+nm+s)+x
    end do
    do q=1,nphys
      H(q,nphys+nm+r)=Up(r,q)
      H(nphys+nm+r,q)=Up(r,q)
    end do
  end do
  if(.not.all(ieee_is_finite(H))) error stop 'TDA Faddeev: nonfinite matrix'
  asym=maxval(abs(H-transpose(H)))
  if(asym>1d-9*max(1d0,maxval(abs(H)))) then
    write(*,*) 'TDA Faddeev: asymmetric pair kernels, max |H-H^T| = ',asym
    error stop 'TDA Faddeev requires symmetric static kernels'
  end if
contains
  integer function spin(t)
    integer,intent(in) :: t
    spin=1-2*mod(t-1,2)
  end function

end subroutine

subroutine G_TDA_Faddeev_spectrum(dressed,adc2,nOrb,nC,nR,nH,eHF,H,ePA,ZPA)
  implicit none
  logical,intent(in) :: dressed,adc2
  integer,intent(in) :: nOrb,nC,nR,nH
  double precision,intent(in) :: eHF(nOrb)
  double precision,intent(inout) :: H(nH,nH)
  double precision,intent(out) :: ePA(nOrb),ZPA(nOrb)
  if(dressed) then
    write(*,*) 'TDA Faddeev-Parquet: static dressed source approximation'
  else if(adc2) then
    write(*,*) 'TDA Faddeev-Parquet with ADC(2) external couplings'
  else
    write(*,*) 'TDA Faddeev-Parquet with ADC(3) external couplings'
  end if
  call Faddeev_spectrum('TDA Faddeev-Parquet',nOrb,nC,nR,nH,eHF,H,ePA,ZPA)
end subroutine

subroutine Faddeev_spectrum(label,nOrb,nC,nR,nH,eHF,H,ePA,ZPA)
! The upfolded matrix realizes Sigma=U^T (omega I-Mbar)^(-1) U.
! Dominant poles use orbital-specific weights summed over degeneracies.
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  include 'parameters.h'
  character(*),intent(in) :: label
  integer,intent(in) :: nOrb,nC,nR,nH
  double precision,intent(in) :: eHF(nOrb)
  double precision,intent(inout) :: H(nH,nH)
  double precision,intent(out) :: ePA(nOrb),ZPA(nOrb)
  double precision,allocatable :: poles(:),weight(:)
  double precision :: best,w
  integer :: p,s,t,idx
  allocate(poles(nH),weight(nH))
  write(*,*) 'Physical upfolded dimension = ',nH
  call diagonalize_matrix(nH,H,poles)
  if(.not.all(ieee_is_finite(poles))) error stop 'Faddeev: nonfinite spectrum'
  weight=sum(H(nC+1:nOrb-nR,:)**2,dim=1)
  write(*,'(A,F18.10)') ' Sum of active 1p weights = ',sum(weight)
  if(abs(sum(weight)-dble(nOrb-nC-nR))>1d-7) error stop 'Faddeev: spectral sum rule failed'
  write(*,*) label// ' poles: index, energy (eV), active 1p weight'
  do s=1,nH
    if(weight(s)>1d-3) write(*,'(I7,2F20.10)') s,poles(s)*HaToeV,weight(s)
  end do
  ePA=eHF; ZPA=1d0
  do p=nC+1,nOrb-nR
    best=-1d0; idx=1; s=1
    do while(s<=nH)
      t=s
      do while(t<nH)
        if(abs(poles(t+1)-poles(s))>1d-8) exit
        t=t+1
      end do
      w=sum(H(p,s:t)**2)
      if(w>best) then
        best=w; idx=s
      end if
      s=t+1
    end do
    ePA(p)=poles(idx); ZPA(p)=best
  end do
  write(*,*) label// ' dominant poles: orbital, HF (eV), QP (eV), Z'
  do p=nC+1,nOrb-nR
    write(*,'(I7,3F20.10)') p,eHF(p)*HaToeV,ePA(p)*HaToeV,ZPA(p)
  end do
end subroutine

subroutine G_TDA_Faddeev_vertices(adc2,reg_PA,eta,nERI,nOrb,nC,nO,nR,nphys,nm,np, &
                                 phys,minus,plus,eQP,ERI,Um,Up, &
    dressed,Teh,Thh,Tee)
! ADC(2), ADC(3), or dressed endpoints shared by dense and matrix-free solvers.
! nERI=nOrb: general spin-orbital ERIs; 2*nERI=nOrb: spatial Coulomb
! integrals with interleaved alpha/beta spin indices. Cache MP2 amplitudes
! once instead of evaluating the same denominators for every external p.
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  logical,intent(in) :: adc2,reg_PA
  integer,intent(in) :: nERI,nOrb,nC,nO,nR,nphys,nm,np
  integer,intent(in) :: phys(nphys),minus(3,nm),plus(3,np)
  double precision,intent(in) :: eta,eQP(nOrb),ERI(nERI,nERI,nERI,nERI)
  double precision,intent(out) :: Um(nm,nphys),Up(np,nphys)
  logical,intent(in) :: dressed
  double precision,intent(in) :: Teh((nO-nC)*(nOrb-nR-nO),(nO-nC)*(nOrb-nR-nO))
  double precision,intent(in) :: Thh((nO-nC)*(nO-nC-1)/2,(nOrb-nR-nO)*(nOrb-nR-nO-1)/2)
  double precision,intent(in) :: Tee((nO-nC)*(nO-nC-1)/2,(nOrb-nR-nO)*(nOrb-nR-nO-1)/2)
  double precision,allocatable :: t2(:,:,:,:)
  double precision :: x
  integer :: i,j,k,l,a,b,c,d,p,q,r,ns,nh,ne
  integer,allocatable :: ph(:,:),hh(:,:),ee(:,:)
  if(adc2.and.dressed) error stop 'Faddeev vertices: incompatible coupling selections'
  if(nERI/=nOrb.and.2*nERI/=nOrb) error stop 'Faddeev vertices: invalid ERI dimension'
  Um=0d0; Up=0d0
  if(adc2.and..not.dressed) then
    do q=1,nphys
      p=phys(q)
      if(p<=nC.or.p>nOrb-nR) cycle
      do r=1,nm
        Um(r,q)=v(minus(1,r),minus(2,r),p,minus(3,r))
      end do
      do r=1,np
        Up(r,q)=v(plus(2,r),plus(3,r),p,plus(1,r))
      end do
    end do
    return
  end if
  if(reg_PA.and.(eta<=0d0.or..not.ieee_is_finite(eta))) error stop 'Faddeev vertices: invalid SRG flow'
  if(dressed) then
    ns=(nO-nC)*(nOrb-nR-nO); nh=(nO-nC)*(nO-nC-1)/2
    ne=(nOrb-nR-nO)*(nOrb-nR-nO-1)/2
    allocate(ph(nOrb,nOrb),hh(nOrb,nOrb),ee(nOrb,nOrb))
    ph=0; hh=0; ee=0; r=0
    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        r=r+1; ph(i,a)=r
      end do
    end do
    r=0
    do i=nC+1,nO
      do j=i+1,nO
        r=r+1; hh(i,j)=r
      end do
    end do
    r=0
    do a=nO+1,nOrb-nR
      do b=a+1,nOrb-nR
        r=r+1; ee(a,b)=r
      end do
    end do
  else
  allocate(t2(nC+1:nO,nC+1:nO,nO+1:nOrb-nR,nO+1:nOrb-nR))
  do b=nO+1,nOrb-nR
    do a=nO+1,nOrb-nR
      do j=nC+1,nO
        do i=nC+1,nO
          t2(i,j,a,b)=v(i,j,a,b)*inv(eQP(i)+eQP(j)-eQP(a)-eQP(b))
        end do
      end do
    end do
  end do
  end if
  Um=0d0; Up=0d0
  do q=1,nphys
    p=phys(q)
    if(p<=nC.or.p>nOrb-nR) cycle
    do r=1,nm
      i=minus(1,r); j=minus(2,r); a=minus(3,r)
      x=v(i,j,p,a)
      do c=nO+1,nOrb-nR
        do d=c+1,nOrb-nR
          if(dressed) then
            x=x+v(c,d,p,a)*Thh(hh(i,j),ee(c,d))
          else
            x=x+v(c,d,p,a)*t2(i,j,c,d)
          end if
        end do
      end do
      do k=nC+1,nO
        do c=nO+1,nOrb-nR
          if(dressed) then
            x=x+v(c,j,p,k)*Teh(ph(k,c),ph(i,a))-v(c,i,p,k)*Teh(ph(k,c),ph(j,a))
          else
            x=x-v(c,j,p,k)*t2(i,k,c,a)+v(c,i,p,k)*t2(j,k,c,a)
          end if
        end do
      end do
      Um(r,q)=x
    end do
    do r=1,np
      i=plus(1,r); a=plus(2,r); b=plus(3,r)
      x=v(a,b,p,i)
      do k=nC+1,nO
        do l=k+1,nO
          if(dressed) then
            x=x+v(k,l,p,i)*Tee(hh(k,l),ee(a,b))
          else
            x=x+v(k,l,p,i)*t2(k,l,a,b)
          end if
        end do
      end do
      do k=nC+1,nO
        do c=nO+1,nOrb-nR
          if(dressed) then
            x=x+v(k,b,p,c)*Teh(ph(k,c),ph(i,a))-v(k,a,p,c)*Teh(ph(k,c),ph(i,b))
          else
            x=x-v(k,b,p,c)*t2(k,i,a,c)+v(k,a,p,c)*t2(k,i,b,c)
          end if
        end do
      end do
      Up(r,q)=x
    end do
  end do
contains
  double precision function v(t,u,w,z)
    integer,intent(in) :: t,u,w,z
    if(nERI==nOrb) then
      v=ERI(t,u,w,z)-ERI(t,u,z,w)
    else
      v=0d0
      if(mod(t-w,2)==0.and.mod(u-z,2)==0) v=ERI((t+1)/2,(u+1)/2,(w+1)/2,(z+1)/2)
      if(mod(t-z,2)==0.and.mod(u-w,2)==0) v=v-ERI((t+1)/2,(u+1)/2,(z+1)/2,(w+1)/2)
    end if
  end function
  double precision function inv(den)
    double precision,intent(in) :: den
    double precision :: y
    if(reg_PA) then
      y=2d0*eta*den**2
      if(y<1d-5) then
        inv=2d0*eta*den*(1d0-y/2d0+y*y/6d0)
      else
        inv=(1d0-exp(-y))/den
      end if
    else
      if(abs(den)<1d-12) error stop 'Faddeev vertices: singular ADC(3) denominator'
      inv=1d0/den
    end if
  end function
end subroutine

subroutine TDA_Faddeev_sources(reg,eta,nOrb,nC,nO,nR,ns,nh,ne,equal_pair,e,Wph,Wh,We,Bph,Bpp,Teh,Thh,Tee)
! Static one-sided source approximation. B matrices are sources only:
! no backward amplitudes are included in the TDA charged propagator.
! Dph T + T Aph = -Bph; Kh Th - Th De = Bpp^T;
! Dh Te - Te Ke = Bpp^T. With SRG, regularize each spectral denominator.
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  logical,intent(in) :: reg,equal_pair
  integer,intent(in) :: nOrb,nC,nO,nR,ns,nh,ne
  double precision,intent(in) :: eta,e(nOrb),Wph(ns,ns),Wh(nh,nh),We(ne,ne),Bph(ns,ns),Bpp(ne,nh)
  double precision,intent(out) :: Teh(ns,ns),Thh(nh,ne),Tee(nh,ne)
  double precision :: X(ns,ns),Y(nh,nh),Z(ne,ne),om(ns),oh(nh),oe(ne),dp(ns),dh(nh),de(ne)
  double precision :: tmp(ns,ns),tmph(nh,ne),tmpe(nh,ne)
  integer :: i,j,a,b,r,s,offset
  if(reg.and.(eta<=0d0.or..not.ieee_is_finite(eta))) error stop 'Faddeev source: invalid SRG flow'
  call check_kernel(ns,Wph)
  call check_kernel(nh,Wh)
  call check_kernel(ne,We)
  if(.not.all(ieee_is_finite(Bph)).or..not.all(ieee_is_finite(Bpp))) &
    error stop 'Faddeev source: nonfinite source kernel'
  r=0
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      r=r+1; dp(r)=e(a)-e(i)
    end do
  end do
  if(r/=ns) error stop 'Faddeev source: incomplete eh space'
  offset=1
  if(equal_pair) offset=0
  r=0
  do i=nC+1,nO
    do j=i+offset,nO
      r=r+1; dh(r)=e(i)+e(j)
    end do
  end do
  if(r/=nh) error stop 'Faddeev source: incomplete hh space'
  r=0
  do a=nO+1,nOrb-nR
    do b=a+offset,nOrb-nR
      r=r+1; de(r)=e(a)+e(b)
    end do
  end do
  if(r/=ne) error stop 'Faddeev source: incomplete ee space'
  X=Wph; Y=-Wh; Z=We
  do r=1,ns
    X(r,r)=X(r,r)+dp(r)
  end do
  do r=1,nh
    Y(r,r)=Y(r,r)+dh(r)
  end do
  do r=1,ne
    Z(r,r)=Z(r,r)+de(r)
  end do
  call diagonalize_matrix(ns,X,om)
  if(nh>0) call diagonalize_matrix(nh,Y,oh)
  if(ne>0) call diagonalize_matrix(ne,Z,oe)
  tmp=matmul(Bph,X)
  do s=1,ns
    do r=1,ns
      tmp(r,s)=-tmp(r,s)*inv(dp(r)+om(s))
    end do
  end do
  Teh=matmul(tmp,transpose(X))
  tmph=matmul(transpose(Y),transpose(Bpp))
  tmpe=matmul(transpose(Bpp),Z)
  do s=1,ne
    do r=1,nh
      tmph(r,s)=tmph(r,s)*inv(oh(r)-de(s))
      tmpe(r,s)=tmpe(r,s)*inv(dh(r)-oe(s))
    end do
  end do
  Thh=matmul(Y,tmph); Tee=matmul(tmpe,transpose(Z))
  if(.not.all(ieee_is_finite(Teh)).or..not.all(ieee_is_finite(Thh)).or. &
     .not.all(ieee_is_finite(Tee))) error stop 'Faddeev source: nonfinite amplitude'
contains
  subroutine check_kernel(n,W)
    integer,intent(in) :: n
    double precision,intent(in) :: W(n,n)
    if(n==0) return
    if(.not.all(ieee_is_finite(W))) error stop 'Faddeev source: nonfinite pair kernel'
    if(maxval(abs(W-transpose(W)))>1d-9*max(1d0,maxval(abs(W)))) &
      error stop 'Faddeev source requires symmetric static pair kernels'
  end subroutine
  double precision function inv(den)
    double precision,intent(in) :: den
    double precision :: y
    if(.not.ieee_is_finite(den)) error stop 'Faddeev source: nonfinite denominator'
    if(reg) then
      y=2d0*eta*den**2
      if(y<1d-5) then
        inv=2d0*eta*den*(1d0-y/2d0+y*y/6d0)
      else
        inv=(1d0-exp(-y))/den
      end if
    else
      if(abs(den)<1d-12) error stop 'Faddeev source: singular static denominator'
      inv=1d0/den
    end if
  end function
end subroutine
