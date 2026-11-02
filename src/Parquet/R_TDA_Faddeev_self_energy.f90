subroutine R_TDA_Faddeev_self_energy(adc2,reg_PA,eta,nOrb,nC,nO,nR,nS,nOOs,nOOt,nVVs,nVVt, &
                                    eHF,eQP,ERI,Wehs,Weht,Whhs,Whht,Wees,Weet,ePA,ZPA)
! Restricted complete-TDA Faddeev-Parquet, ADC(2) or ADC(3) endpoints.
! Lift singlet/triplet static pair interactions to spin orbitals and retain
! only the alpha charged determinant block (Ms=1/2). This is an exact spin
! representation, including dark quartet states; no numerical projection
! or duplicate physical 1p states. The external orbitals remain spatial.
! Inputs contain interactions only, with all free pair energies removed.
  implicit none
  logical,intent(in) :: adc2,reg_PA
  integer,intent(in) :: nOrb,nC,nO,nR,nS,nOOs,nOOt,nVVs,nVVt
  double precision,intent(in) :: eta,eHF(nOrb),eQP(nOrb),ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in) :: Wehs(nS,nS),Weht(nS,nS),Whhs(nOOs,nOOs),Whht(nOOt,nOOt)
  double precision,intent(in) :: Wees(nVVs,nVVs),Weet(nVVt,nVVt)
  double precision,intent(out) :: ePA(nOrb),ZPA(nOrb)
  integer :: ng,nsG,nooG,nvvG,noact,nv,nH,p,q,r,s
  double precision,allocatable :: eg(:),fg(:),vg(:,:,:,:),Weh(:,:),Whh(:,:),Wee(:,:),H(:,:)
  noact=nO-nC; nv=nOrb-nR-nO
  if(nC<0.or.nR<0.or.noact<1.or.nv<1) error stop 'R TDA Faddeev: invalid active space'
  if(nS/=noact*nv.or.nOOs/=noact*(noact+1)/2.or.nOOt/=noact*(noact-1)/2.or. &
     nVVs/=nv*(nv+1)/2.or.nVVt/=nv*(nv-1)/2) error stop 'R TDA Faddeev: inconsistent pair spaces'
  ng=2*nOrb; nsG=4*nS; nooG=noact*(2*noact-1); nvvG=nv*(2*nv-1)
  nH=nOrb+nv*(noact*noact+noact*(noact-1)/2)+noact*(nv*nv+nv*(nv-1)/2)
  allocate(eg(ng),fg(ng),vg(ng,ng,ng,ng),Weh(nsG,nsG),Whh(nooG,nooG),Wee(nvvG,nvvG),H(nH,nH))
  do p=1,ng
    eg(p)=eQP((p+1)/2); fg(p)=eHF((p+1)/2)
  end do
! Physicist-order Coulomb integrals: spin conservation on each electron.
  vg=0d0
  do s=1,ng
    do r=1,ng
      do q=1,ng
        if(mod(q-s,2)/=0) cycle
        do p=1,ng
          if(mod(p-r,2)/=0) cycle
          vg(p,q,r,s)=ERI((p+1)/2,(q+1)/2,(r+1)/2,(s+1)/2)
        end do
      end do
    end do
  end do
  call R_TDA_Faddeev_lift(nOrb,nC,nO,nR,nS,nOOs,nOOt,nVVs,nVVt, &
                         Wehs,Weht,Whhs,Whht,Wees,Weet,Weh,Whh,Wee)
  call G_TDA_Faddeev_build(adc2,reg_PA,eta,1,ng,2*nC,2*nO,2*nR,nsG,nooG,nvvG,nH,fg,eg,vg,Weh,Whh,Wee,H)
  write(*,*) 'Restricted TDA Faddeev-Parquet: alpha charged determinant block'
  call G_TDA_Faddeev_spectrum(adc2,nOrb,nC,nR,nH,eHF,H,ePA,ZPA)
end subroutine

subroutine R_TDA_Faddeev_IP(adc2,reg,eta,nOrb,nC,nO,nR,nS,nOOs,nOOt,nVVs,nVVt, &
                           eHF,eQP,ERI,Wehs,Weht,Whhs,Whht,Wees,Weet,omega,Z,homo_weight,residual,info)
! Matrix-free restricted principal IP. Spatial ERIs are passed directly;
! the shared endpoint builder applies spin deltas without expanding ERIs.
  implicit none
  logical,intent(in) :: adc2,reg
  integer,intent(in) :: nOrb,nC,nO,nR,nS,nOOs,nOOt,nVVs,nVVt
  double precision,intent(in) :: eta,eHF(nOrb),eQP(nOrb),ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in) :: Wehs(nS,nS),Weht(nS,nS),Whhs(nOOs,nOOs),Whht(nOOt,nOOt)
  double precision,intent(in) :: Wees(nVVs,nVVs),Weet(nVVt,nVVt)
  double precision,intent(out) :: omega,Z,homo_weight,residual
  integer,intent(out) :: info
  integer :: ng,nsG,nooG,nvvG,noact,nv,p
  double precision,allocatable :: eg(:),fg(:),Weh(:,:),Whh(:,:),Wee(:,:)
  noact=nO-nC; nv=nOrb-nR-nO
  if(nC<0.or.nR<0.or.noact<1.or.nv<1) error stop 'R Faddeev IP: invalid active space'
  if(nS/=noact*nv.or.nOOs/=noact*(noact+1)/2.or.nOOt/=noact*(noact-1)/2.or. &
     nVVs/=nv*(nv+1)/2.or.nVVt/=nv*(nv-1)/2) error stop 'R Faddeev IP: incomplete pair spaces'
  ng=2*nOrb; nsG=4*nS; nooG=noact*(2*noact-1); nvvG=nv*(2*nv-1)
  allocate(eg(ng),fg(ng),Weh(nsG,nsG),Whh(nooG,nooG),Wee(nvvG,nvvG))
  do p=1,ng
    eg(p)=eQP((p+1)/2); fg(p)=eHF((p+1)/2)
  end do
  call R_TDA_Faddeev_lift(nOrb,nC,nO,nR,nS,nOOs,nOOt,nVVs,nVVt, &
                         Wehs,Weht,Whhs,Whht,Wees,Weet,Weh,Whh,Wee)
  call G_TDA_Faddeev_IP(adc2,reg,eta,1,nOrb,ng,2*nC,2*nO,2*nR,nsG,nooG,nvvG, &
                       fg,eg,ERI,Weh,Whh,Wee,omega,Z,homo_weight,residual,info)
end subroutine

subroutine R_TDA_Faddeev_IP_roots(nroots,adc2,reg,eta,nOrb,nC,nO,nR,nS,nOOs,nOOt,nVVs,nVVt, &
                           eHF,eQP,ERI,Wehs,Weht,Whhs,Whht,Wees,Weet,omega,Z,residual,info)
! Matrix-free restricted IP roots. Spatial ERIs are passed directly;
! the shared endpoint builder applies spin deltas without expanding ERIs.
  implicit none
  logical,intent(in) :: adc2,reg
  integer,intent(in) :: nroots,nOrb,nC,nO,nR,nS,nOOs,nOOt,nVVs,nVVt
  double precision,intent(in) :: eta,eHF(nOrb),eQP(nOrb),ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in) :: Wehs(nS,nS),Weht(nS,nS),Whhs(nOOs,nOOs),Whht(nOOt,nOOt)
  double precision,intent(in) :: Wees(nVVs,nVVs),Weet(nVVt,nVVt)
  double precision,intent(out) :: omega(nroots),Z(nroots),residual(nroots)
  integer,intent(out) :: info
  integer :: ng,nsG,nooG,nvvG,noact,nv,p
  double precision,allocatable :: eg(:),fg(:),Weh(:,:),Whh(:,:),Wee(:,:)
  noact=nO-nC; nv=nOrb-nR-nO
  if(nC<0.or.nR<0.or.noact<1.or.nv<1) error stop 'R Faddeev IP: invalid active space'
  if(nS/=noact*nv.or.nOOs/=noact*(noact+1)/2.or.nOOt/=noact*(noact-1)/2.or. &
     nVVs/=nv*(nv+1)/2.or.nVVt/=nv*(nv-1)/2) error stop 'R Faddeev IP: incomplete pair spaces'
  ng=2*nOrb; nsG=4*nS; nooG=noact*(2*noact-1); nvvG=nv*(2*nv-1)
  allocate(eg(ng),fg(ng),Weh(nsG,nsG),Whh(nooG,nooG),Wee(nvvG,nvvG))
  do p=1,ng
    eg(p)=eQP((p+1)/2); fg(p)=eHF((p+1)/2)
  end do
  call R_TDA_Faddeev_lift(nOrb,nC,nO,nR,nS,nOOs,nOOt,nVVs,nVVt, &
                         Wehs,Weht,Whhs,Whht,Wees,Weet,Weh,Whh,Wee)
  call G_TDA_Faddeev_IP_roots(nroots,adc2,reg,eta,1,nOrb,ng,2*nC,2*nO,2*nR,nsG,nooG,nvvG, &
                       fg,eg,ERI,Weh,Whh,Wee,omega,Z,residual,info)
end subroutine

subroutine R_TDA_Faddeev_lift(nOrb,nC,nO,nR,nS,nOOs,nOOt,nVVs,nVVt, &
                             Wehs,Weht,Whhs,Whht,Wees,Weet,Weh,Whh,Wee)
! Orthonormal Clebsch-Gordan transforms. The eh singlet is (aa+bb)/sqrt2,
! the m=0 triplet (aa-bb)/sqrt2; spin-flip blocks are triplets as well.
  implicit none
  integer,intent(in) :: nOrb,nC,nO,nR,nS,nOOs,nOOt,nVVs,nVVt
  double precision,intent(in) :: Wehs(nS,nS),Weht(nS,nS),Whhs(nOOs,nOOs),Whht(nOOt,nOOt)
  double precision,intent(in) :: Wees(nVVs,nVVs),Weet(nVVt,nVVt)
  double precision,intent(out) :: Weh(4*nS,4*nS)
  double precision,intent(out) :: Whh((nO-nC)*(2*(nO-nC)-1),(nO-nC)*(2*(nO-nC)-1))
  double precision,intent(out) :: Wee((nOrb-nR-nO)*(2*(nOrb-nR-nO)-1),(nOrb-nR-nO)*(2*(nOrb-nR-nO)-1))
  integer :: i,a,j,b,ia,jb,u,v,nv,si,sa,sj,sb
  nv=nOrb-nR-nO
  Weh=0d0
  ia=0
  do i=2*nC+1,2*nO
    do a=2*nO+1,2*(nOrb-nR)
      ia=ia+1
      u=(((i+1)/2)-nC-1)*nv+(a+1)/2-nO
      si=mod(i-1,2); sa=mod(a-1,2)
      jb=0
      do j=2*nC+1,2*nO
        do b=2*nO+1,2*(nOrb-nR)
          jb=jb+1
          v=(((j+1)/2)-nC-1)*nv+(b+1)/2-nO
          sj=mod(j-1,2); sb=mod(b-1,2)
          if(si==sa.and.sj==sb) then
            if(si==sj) then
              Weh(ia,jb)=0.5d0*(Wehs(u,v)+Weht(u,v))
            else
              Weh(ia,jb)=0.5d0*(Wehs(u,v)-Weht(u,v))
            end if
          else if(si/=sa.and.si==sj.and.sa==sb) then
            Weh(ia,jb)=Weht(u,v)
          end if
        end do
      end do
    end do
  end do
  call lift_pair(nC+1,nO,nOOs,nOOt,size(Whh,1),Whhs,Whht,Whh)
  call lift_pair(nO+1,nOrb-nR,nVVs,nVVt,size(Wee,1),Wees,Weet,Wee)
contains
  subroutine lift_pair(first,last,ns,nt,ng,Ws,Wt,Wg)
    integer,intent(in) :: first,last,ns,nt,ng
    double precision,intent(in) :: Ws(ns,ns),Wt(nt,nt)
    double precision,intent(out) :: Wg(ng,ng)
    integer :: smap(last,last),tmap(last,last),sid(ng),tid(ng),m(ng)
    double precision :: cs(ng),ct(ng)
    integer :: p,q,u,v,r,s,n
    smap=0; tmap=0; r=0; s=0
    do u=first,last
      do v=u,last
        r=r+1; smap(u,v)=r
        if(u/=v) then
          s=s+1; tmap(u,v)=s
        end if
      end do
    end do
    n=0; cs=0d0; ct=0d0; sid=0; tid=0
    do p=2*first-1,2*last
      do q=p+1,2*last
        n=n+1; u=(p+1)/2; v=(q+1)/2
        m(n)=1-mod(p-1,2)-mod(q-1,2)
        if(m(n)==0) then
          sid(n)=smap(u,v)
          if(u==v) then
            cs(n)=1d0
          else
            cs(n)=dble(1-2*mod(p-1,2))/sqrt(2d0)
            ct(n)=1d0/sqrt(2d0); tid(n)=tmap(u,v)
          end if
        else
          ct(n)=1d0; tid(n)=tmap(u,v)
        end if
      end do
    end do
    if(n/=ng) error stop 'R TDA Faddeev: wrong spin pair dimension'
    Wg=0d0
    do r=1,ng
      do s=1,ng
        if(sid(r)>0.and.sid(s)>0) Wg(r,s)=cs(r)*cs(s)*Ws(sid(r),sid(s))
        if(tid(r)>0.and.tid(s)>0.and.m(r)==m(s)) Wg(r,s)=Wg(r,s)+ct(r)*ct(s)*Wt(tid(r),tid(s))
      end do
    end do
  end subroutine
end subroutine
