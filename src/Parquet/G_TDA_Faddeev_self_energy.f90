subroutine G_TDA_Faddeev_self_energy(reg_PA,eta,nOrb,nC,nO,nR,nS,nOO,nVV, &
                                    eHF,eQP,ERI,Weh,Whh,Wee,ePA,ZPA)
! Complete-TDA physical reduction of faddeev_parquet_notes.tex.
! Weh=Aeh-diag(ea-ei), Whh=Dpp+diag(ei+ej), Wee=Cpp-diag(ea+eb).
! External vertices are U^(1)+U^(2) of ADC(3), with bare ERIs.
  implicit none
  logical,intent(in) :: reg_PA
  integer,intent(in) :: nOrb,nC,nO,nR,nS,nOO,nVV
  double precision,intent(in) :: eta,eHF(nOrb),eQP(nOrb),ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in) :: Weh(nS,nS),Whh(nOO,nOO),Wee(nVV,nVV)
  double precision,intent(out) :: ePA(nOrb),ZPA(nOrb)
  integer :: nH
  double precision,allocatable :: H(:,:)
  if(nC<0.or.nR<0.or.nO<=nC.or.nOrb-nR<=nO) error stop 'G TDA Faddeev: invalid active space'
  nH=nOrb+nOO*(nOrb-nR-nO)+nVV*(nO-nC)
  allocate(H(nH,nH))
  call G_TDA_Faddeev_build(reg_PA,eta,0,nOrb,nC,nO,nR,nS,nOO,nVV,nH,eHF,eQP,ERI,Weh,Whh,Wee,H)
  call G_TDA_Faddeev_spectrum(nOrb,nC,nR,nH,eHF,H,ePA,ZPA)
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

subroutine G_TDA_Faddeev_build(reg_PA,eta,spin_block,nOrb,nC,nO,nR,nS,nOO,nVV,nH, &
                              eHF,eQP,ERI,Weh,Whh,Wee,H)
! spin_block=0: general spin orbitals; 1: alpha charged block, with
! interleaved alpha/beta orbitals (used by the restricted wrapper).
! Kernel signs follow the actual QuAcK TDA pair Hamiltonians, avoiding
! four-index crossing ambiguities. Pair bases must be complete.
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  logical,intent(in) :: reg_PA
  integer,intent(in) :: spin_block,nOrb,nC,nO,nR,nS,nOO,nVV,nH
  double precision,intent(in) :: eta,eHF(nOrb),eQP(nOrb),ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in) :: Weh(nS,nS),Whh(nOO,nOO),Wee(nVV,nVV)
  double precision,intent(out) :: H(nH,nH)
  integer :: noact,nvact,nphys,nm,np,p,q,i,j,k,l,a,b,c,d,r,s,ij,kl,ab,cd
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
  call G_TDA_Faddeev_vertices(reg_PA,eta,nOrb,nOrb,nC,nO,nR,nphys,nm,np, &
                              phys,minus,plus,eQP,ERI,Um,Up)
  H=0d0
  do p=1,nphys
    H(p,p)=eHF(phys(p))
  end do

! Removal block and the three ADC(3) endpoint terms, notes Eq. UR2_minus.
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

! Addition block and Eq. UR2_plus. Ordered pair sums absorb the 1/2.
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

subroutine G_TDA_Faddeev_spectrum(nOrb,nC,nR,nH,eHF,H,ePA,ZPA)
! The upfolded matrix realizes Sigma=U^T (omega I-Mbar)^(-1) U.
! Dominant poles use orbital-specific weights summed over degeneracies.
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  include 'parameters.h'
  integer,intent(in) :: nOrb,nC,nR,nH
  double precision,intent(in) :: eHF(nOrb)
  double precision,intent(inout) :: H(nH,nH)
  double precision,intent(out) :: ePA(nOrb),ZPA(nOrb)
  double precision,allocatable :: poles(:),weight(:)
  double precision :: best,w
  integer :: p,s,t,idx
  allocate(poles(nH),weight(nH))
  write(*,*) 'TDA Faddeev-Parquet with ADC(3) external couplings'
  write(*,*) 'Physical upfolded dimension = ',nH
  call diagonalize_matrix(nH,H,poles)
  if(.not.all(ieee_is_finite(poles))) error stop 'TDA Faddeev: nonfinite spectrum'
  weight=sum(H(nC+1:nOrb-nR,:)**2,dim=1)
  write(*,'(A,F18.10)') ' Sum of active 1p weights = ',sum(weight)
  if(abs(sum(weight)-dble(nOrb-nC-nR))>1d-7) error stop 'TDA Faddeev: spectral sum rule failed'
  write(*,*) 'TDA Faddeev-Parquet poles: index, energy (eV), active 1p weight'
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
  write(*,*) 'TDA Faddeev-Parquet dominant poles: orbital, HF (eV), QP (eV), Z'
  do p=nC+1,nOrb-nR
    write(*,'(I7,3F20.10)') p,eHF(p)*HaToeV,ePA(p)*HaToeV,ZPA(p)
  end do
end subroutine

subroutine G_TDA_Faddeev_vertices(reg_PA,eta,nERI,nOrb,nC,nO,nR,nphys,nm,np, &
                                 phys,minus,plus,eQP,ERI,Um,Up)
! ADC(3) endpoints shared by the dense and matrix-free implementations.
! nERI=nOrb: general spin-orbital ERIs; 2*nERI=nOrb: spatial Coulomb
! integrals with interleaved alpha/beta spin indices. Cache MP2 amplitudes
! once instead of evaluating the same denominators for every external p.
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  logical,intent(in) :: reg_PA
  integer,intent(in) :: nERI,nOrb,nC,nO,nR,nphys,nm,np
  integer,intent(in) :: phys(nphys),minus(3,nm),plus(3,np)
  double precision,intent(in) :: eta,eQP(nOrb),ERI(nERI,nERI,nERI,nERI)
  double precision,intent(out) :: Um(nm,nphys),Up(np,nphys)
  double precision,allocatable :: t2(:,:,:,:)
  double precision :: x
  integer :: i,j,k,l,a,b,c,d,p,q,r
  if(nERI/=nOrb.and.2*nERI/=nOrb) error stop 'Faddeev vertices: invalid ERI dimension'
  if(reg_PA.and.(eta<=0d0.or..not.ieee_is_finite(eta))) error stop 'Faddeev vertices: invalid SRG flow'
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
  Um=0d0; Up=0d0
  do q=1,nphys
    p=phys(q)
    if(p<=nC.or.p>nOrb-nR) cycle
    do r=1,nm
      i=minus(1,r); j=minus(2,r); a=minus(3,r)
      x=v(i,j,p,a)
      do c=nO+1,nOrb-nR
        do d=c+1,nOrb-nR
          x=x+v(c,d,p,a)*t2(i,j,c,d)
        end do
      end do
      do k=nC+1,nO
        do c=nO+1,nOrb-nR
          x=x-v(c,j,p,k)*t2(i,k,c,a)+v(c,i,p,k)*t2(j,k,c,a)
        end do
      end do
      Um(r,q)=x
    end do
    do r=1,np
      i=plus(1,r); a=plus(2,r); b=plus(3,r)
      x=v(a,b,p,i)
      do k=nC+1,nO
        do l=k+1,nO
          x=x+v(k,l,p,i)*t2(k,l,a,b)
        end do
      end do
      do k=nC+1,nO
        do c=nO+1,nOrb-nR
          x=x-v(k,b,p,c)*t2(k,i,a,c)+v(k,a,p,c)*t2(k,i,b,c)
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
