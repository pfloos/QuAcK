module R_FRPA_operator
  use TDA_Faddeev_IP_operator, only: Faddeev_operator, initialize, apply, FRPA_channel, davidson_roots
  implicit none
  private
  public :: build_FRPA, run_faddeev_pairs, ccd_sources, FRPA_channel
contains

  subroutine build_FRPA(nOrb,nC,nO,nR,eHF,ERI,backward,H,op,qp,kernels,coupling,static3,external_gw, &
                        iterative,reg,eta,Teh,Thh,Tee,external_energies,ccd_ph,ccd_pp)
! Restricted alpha charged determinant block. Absent options reproduce the legacy FRPA calculation.
    integer,intent(in) :: nOrb,nC,nO,nR
    double precision,intent(in) :: eHF(nOrb),ERI(nOrb,nOrb,nOrb,nOrb)
    logical,intent(in) :: backward
    double precision,allocatable,intent(out) :: H(:,:)
    type(Faddeev_operator),intent(out) :: op
    double precision,optional,intent(in) :: qp(nOrb),eta,Teh(:,:),Thh(:,:),Tee(:,:),external_energies(nOrb)
    type(FRPA_channel),optional,intent(in) :: kernels(3)
    double precision,optional,intent(in) :: ccd_ph(:,:),ccd_pp(:,:)
    integer,optional,intent(in) :: coupling
    logical,optional,intent(in) :: static3,external_gw,iterative,reg
    logical :: use_static,use_gw,use_iterative,adc2,dressed,regularize
    double precision :: flow
    double precision,allocatable :: hf(:),trial(:),source_ph(:,:),source_hh(:,:),source_ee(:,:),correction(:,:)
    integer :: nocc_spin,nv,ns,nh,ne,i,j,a,b,k,l,c,d,r,s,p
    double precision :: mu
    integer,allocatable :: holes(:,:),particles(:,:)
    double precision,allocatable :: e(:),Aph(:,:),Bph(:,:),Cee(:,:),Dhh(:,:),Bpp(:,:)
    double precision,allocatable :: X(:,:),Y(:,:),w(:),v(:),Gph(:,:),Kph(:,:),Ghh(:,:),Khh(:,:),Gee(:,:),Kee(:,:)
    double precision,allocatable :: em(:),ep(:),Um(:,:),Up(:,:),sigma(:,:)

    nocc_spin=2*(nO-nC); nv=2*(nOrb-nR-nO)
    if(nocc_spin<2.or.nv<2) error stop 'FRPA requires active occupied and virtual orbitals'
    ns=nocc_spin*nv; nh=nocc_spin*(nocc_spin-1)/2; ne=nv*(nv-1)/2
    allocate(e(2*nOrb),hf(2*nOrb),holes(2,nh),particles(2,ne),Aph(ns,ns),Bph(ns,ns),Cee(ne,ne),Dhh(nh,nh),Bpp(ne,nh))
    do p=1,2*nOrb
      hf(p)=eHF((p+1)/2)
      e(p)=hf(p)
      if(present(qp)) e(p)=qp((p+1)/2)
    end do
    mu=0.5d0*(e(2*nO)+e(2*nO+1))
    use_static=.true.; use_gw=.false.; use_iterative=.false.; adc2=.false.; dressed=.false.
    regularize=.false.; flow=0d0
    if(present(static3)) use_static=static3
    if(present(external_gw)) use_gw=external_gw
    if(present(iterative)) use_iterative=iterative
    if(present(reg)) regularize=reg
    if(present(eta)) flow=eta
    if(present(coupling)) then
      adc2=coupling==2
      dressed=coupling>=4
    end if

    r=0
    do i=2*nC+1,2*nO
      do j=i+1,2*nO
        r=r+1; holes(:,r)=[i,j]
      end do
    end do
    r=0
    do a=2*nO+1,2*(nOrb-nR)
      do b=a+1,2*(nOrb-nR)
        r=r+1; particles(:,r)=[a,b]
      end do
    end do

    if(present(kernels)) then
      Aph=kernels(1)%hamiltonian; Bph=kernels(1)%metric
      Dhh=kernels(2)%hamiltonian
      Cee=kernels(3)%hamiltonian; Bpp=kernels(3)%metric
    else
    do r=1,ns
      i=2*nC+1+(r-1)/nv; a=2*nO+1+mod(r-1,nv)
      do s=1,ns
        j=2*nC+1+(s-1)/nv; b=2*nO+1+mod(s-1,nv)
        Aph(r,s)=antisym(a,j,i,b)
        Bph(r,s)=antisym(a,b,i,j)
      end do
    end do
    do r=1,ne
      a=particles(1,r); b=particles(2,r)
      do s=1,ne
        c=particles(1,s); d=particles(2,s)
        Cee(r,s)=antisym(a,b,c,d)
      end do
      do s=1,nh
        i=holes(1,s); j=holes(2,s)
        Bpp(r,s)=antisym(a,b,i,j)
      end do
    end do
    do r=1,nh
      i=holes(1,r); j=holes(2,r)
      do s=1,nh
        k=holes(1,s); l=holes(2,s)
        Dhh(r,s)=antisym(i,j,k,l)
      end do
    end do

    end if

    if(.not.backward.and.present(ccd_ph).and.present(ccd_pp)) then
! Static CCD closure: Hermitian parts of the three forward Riccati blocks
! A+B*Tph, C+B*Tpp, and D+B^T*Tpp^T. The free energies remain unchanged.
      correction=matmul(Bph,ccd_ph)
      Aph=Aph+0.5d0*(correction+transpose(correction))
      correction=matmul(Bpp,ccd_pp)
      Cee=Cee+0.5d0*(correction+transpose(correction))
      correction=matmul(transpose(Bpp),transpose(ccd_pp))
      Dhh=Dhh+0.5d0*(correction+transpose(correction))
      deallocate(correction)
    end if

    allocate(source_ph(ns,ns),source_hh(nh,ne),source_ee(nh,ne))
    source_ph=0d0; source_hh=0d0; source_ee=0d0
    if(dressed) then
      if(.not.present(Teh).or..not.present(Thh).or..not.present(Tee)) error stop 'Missing Faddeev coupling amplitudes'
      source_ph=Teh; source_hh=Thh; source_ee=Tee
    end if
    call initialize(op,adc2,regularize,flow,1,nOrb,2*nOrb,2*nC,2*nO,2*nR,ns,nh,ne,hf,e,ERI, &
                    Aph,Dhh,Cee,dressed,source_ph,source_hh,source_ee)
    deallocate(source_ph,source_hh,source_ee)
    allocate(op%sigma_static(nOrb,nOrb))
    op%sigma_static=0d0
    if(use_static) call TDA_Faddeev_static3(nOrb,2*nOrb,nOrb,2*nC,2*nO,2*nR,hf,ERI,op%sigma_static)
    do p=1,nOrb
      if(use_gw) op%sigma_static(p,p)=op%sigma_static(p,p)+e(2*p)-eHF(p)
      if(present(external_energies)) op%sigma_static(p,p)=op%sigma_static(p,p)+external_energies(p)-eHF(p)
      op%diag(p)=op%diag(p)+op%sigma_static(p,p)
    end do
    if(.not.backward) then
      if(.not.use_iterative) then
        allocate(H(op%n,op%n),trial(op%n))
        trial=0d0
        do p=1,op%n
          trial(p)=1d0
          call apply(op,trial,H(:,p))
          trial(p)=0d0
        end do
      end if
      return
    end if
    do r=1,ns
      i=2*nC+1+(r-1)/nv; a=2*nO+1+mod(r-1,nv)
      Aph(r,r)=Aph(r,r)+e(a)-e(i)
    end do
    do r=1,ne
      a=particles(1,r); b=particles(2,r)
      Cee(r,r)=Cee(r,r)+e(a)+e(b)-2d0*mu
    end do
    do r=1,nh
      i=holes(1,r); j=holes(2,r)
      Dhh(r,r)=Dhh(r,r)-e(i)-e(j)+2d0*mu
    end do
    if(present(ccd_ph).and.present(ccd_pp)) then
      call ccd_pair_pencil('ph',Aph,Bph,Aph,ccd_ph,Gph,Kph)
      call ccd_pair_pencil('pp addition',Cee,Bpp,Dhh,ccd_pp,Gee,Kee)
      call ccd_pair_pencil('pp removal',Dhh,transpose(Bpp),Cee,transpose(ccd_pp),Ghh,Khh)
      ! Restore the absolute pair energies after projection.
      Kee=Kee+2d0*mu*Gee
      Khh=-Khh+2d0*mu*Ghh
    else
      call solve_RPA('ph',Aph,Bph,Aph,X,Y,w,v)
      call pair_pencil(X,w,Gph,Kph)
! Undo the pp chemical-potential shift before adding the spectator energies.
      call solve_RPA('pp',Cee,Bpp,Dhh,X,Y,w,v)
      call pair_pencil(X,w+2d0*mu,Gee,Kee)
      call pair_pencil(Y,v+2d0*mu,Ghh,Khh)
    end if

    allocate(op%removal(2),op%addition(2))
    call make_channels(op%im,op%hh,e(2*nO+1:2*(nOrb-nR)),e(2*nC+1:2*nO), &
                       Gph,-Kph,Ghh,Khh,.true.,op%removal)
    call make_channels(reshape(op%ip,[nv,nv,nocc_spin],order=[3,1,2]),op%ee, &
                       e(2*nC+1:2*nO),e(2*nO+1:2*(nOrb-nR)),Gph,Kph,Gee,Kee,.false.,op%addition)
    if(use_iterative) then
      allocate(op%metric_diag(op%n))
      op%metric_diag(:nOrb)=1d0
      call channel_diagonal(op%removal,op%free_m,op%metric_diag(nOrb+1:nOrb+op%nm),op%diag(nOrb+1:nOrb+op%nm))
      call channel_diagonal(op%addition,op%free_p,op%metric_diag(nOrb+op%nm+1:),op%diag(nOrb+op%nm+1:))
      return
    end if

    call sector_poles(op%removal,op%free_m,op%Um,em,Um)
    call sector_poles(op%addition,op%free_p,op%Up,ep,Up)
    allocate(H(op%n,op%n))
    H=0d0
    H(:nOrb,:nOrb)=op%sigma_static
    do p=1,nOrb
      H(p,p)=H(p,p)+eHF(p)
    end do
    H(nOrb+1:nOrb+op%nm,:nOrb)=Um
    H(nOrb+op%nm+1:,:nOrb)=Up
    H(:nOrb,nOrb+1:)=transpose(H(nOrb+1:,:nOrb))
    do r=1,op%nm
      H(nOrb+r,nOrb+r)=em(r)
    end do
    do r=1,op%np
      H(nOrb+op%nm+r,nOrb+op%nm+r)=ep(r)
    end do

  contains

    double precision function antisym(p,q,r,s) result(val)
      integer,intent(in) :: p,q,r,s
      val=0d0
      if(mod(p-r,2)==0.and.mod(q-s,2)==0) val=ERI((p+1)/2,(q+1)/2,(r+1)/2,(s+1)/2)
      if(mod(p-s,2)==0.and.mod(q-r,2)==0) val=val-ERI((p+1)/2,(q+1)/2,(s+1)/2,(r+1)/2)
    end function

  end subroutine

  subroutine solve_RPA(channel,A,B,C,X,Y,w,v)
! For S=L L^T, J S is similar to the symmetric L^T J L.
! Normalize the positive and negative modes in the RPA metric J.
    character(*),intent(in) :: channel
    double precision,intent(in) :: A(:,:),B(:,:),C(:,:)
    double precision,allocatable,intent(out) :: X(:,:),Y(:,:),w(:),v(:)
    double precision,allocatable :: L(:,:),Z(:,:),F(:,:),omega(:)
    integer :: na,nc,n,i,j,info
    na=size(A,1); nc=size(C,1); n=na+nc
    allocate(L(n,n),Z(n,n),F(n,n),omega(n))
    L(:na,:na)=A; L(:na,na+1:)=B; L(na+1:,:na)=transpose(B); L(na+1:,na+1:)=C
    call dpotrf('L',n,L,n,info)
    if(info/=0) then
      write(*,*) 'FRPA: unstable ',channel,' RPA Hessian'
      error stop 'FRPA requires stable HF RPA channels'
    end if
    do j=1,n
      L(:j-1,j)=0d0
    end do
    F=L
    F(na+1:,:)=-F(na+1:,:)
    call dgemm('T','N',n,n,n,1d0,L,n,F,n,0d0,Z,n)
    call diagonalize_matrix(n,Z,omega)
    if(omega(nc)>=0d0.or.omega(nc+1)<=0d0) error stop 'FRPA: singular RPA spectrum'
    do i=1,n
      Z(:,i)=Z(:,i)*sqrt(abs(omega(i)))
    end do
    call dtrsm('L','L','T','N',n,n,1d0,L,n,Z,n)
    X=Z(:na,nc+1:); w=omega(nc+1:)
    Y=Z(na+1:,:nc); v=omega(:nc)
  end subroutine

  subroutine ccd_pair_pencil(channel,A,B,C,T,G,K)
! Project the pair Hessian and signature onto the columns of [I; T].
! When T is an exact RPA amplitude Y*X^(-1), these are the same pair
! matrices as in pair_pencil. With CCD T this is a projected approximation.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    character(*),intent(in) :: channel
    double precision,intent(in) :: A(:,:),B(:,:),C(:,:),T(:,:)
    double precision,allocatable,intent(out) :: G(:,:),K(:,:)
    double precision,allocatable :: check(:,:),BT(:,:)
    integer :: n,i,info

    n=size(A,1)
    if(.not.all(ieee_is_finite(T))) error stop 'FRPA: nonfinite CCD pair amplitudes'
    G=-matmul(transpose(T),T)
    do i=1,n
      G(i,i)=G(i,i)+1d0
    end do
    check=G
    call dpotrf('U',n,check,n,info)
    if(info/=0) then
      write(*,*) 'FRPA: non-positive CCD pair metric in ',channel
      error stop 'FRPA: CCD pair amplitudes require a positive metric'
    end if
    BT=matmul(B,T)
    K=A+BT+transpose(BT)+matmul(transpose(T),matmul(C,T))
    K=0.5d0*(K+transpose(K))
    check=K
    call dpotrf('U',n,check,n,info)
    if(info/=0) then
      write(*,*) 'FRPA: non-positive projected CCD pair Hessian in ',channel
      error stop 'FRPA: CCD pair projection requires a stable pair spectrum'
    end if
  end subroutine

  subroutine pair_pencil(X,w,G,K)
    double precision,intent(in) :: X(:,:),w(:)
    double precision,allocatable,intent(out) :: G(:,:),K(:,:)
    double precision,allocatable :: inverse(:,:),work(:,:),scaled(:,:)
    integer :: n,i,info
    integer,allocatable :: pivot(:)
    n=size(w)
    allocate(inverse(n,n),work(n,n),scaled(n,n),pivot(n),G(n,n),K(n,n))
    work=X; inverse=0d0
    do i=1,n
      inverse(i,i)=1d0
    end do
    call dgesv(n,n,work,n,pivot,inverse,n,info)
    if(info/=0) error stop 'FRPA: singular forward RPA amplitudes'
    do i=1,n
      scaled(i,:)=w(i)*inverse(i,:)
    end do
    call dgemm('T','N',n,n,n,1d0,inverse,n,inverse,n,0d0,G,n)
    call dgemm('T','N',n,n,n,1d0,inverse,n,scaled,n,0d0,K,n)
  end subroutine

  subroutine make_channels(map,pair,spectator,like_energy,Gph,Kph,Gpp,Kpp,removal,channels)
    integer,intent(in) :: map(:,:,:),pair(:,:)
    double precision,intent(in) :: spectator(:),like_energy(:),Gph(:,:),Kph(:,:),Gpp(:,:),Kpp(:,:)
    logical,intent(in) :: removal
    type(FRPA_channel),intent(out) :: channels(2)
    integer :: nlike,nother,i,j,a,b,m
    nlike=size(map,1); nother=size(map,3)
    channels(1)%metric=Gph; channels(1)%hamiltonian=Kph; channels(1)%shift=like_energy
    channels(2)%metric=Gpp; channels(2)%hamiltonian=Kpp; channels(2)%shift=-spectator
    allocate(channels(1)%configs(size(Gph,1),nlike),channels(2)%configs(size(Gpp,1),nother))
    do j=1,nlike
      do i=1,nlike
        do a=1,nother
          m=(a-1)*nlike+i
          if(removal) m=(i-1)*nother+a
          channels(1)%configs(m,j)=map(i,j,a)
        end do
      end do
    end do
    do a=1,nother
      do i=1,nlike
        do b=i+1,nlike
          channels(2)%configs(pair(i,b),a)=map(i,b,a)
        end do
      end do
    end do
  end subroutine

  subroutine channel_diagonal(channels,free,metric,matrix)
    type(FRPA_channel),intent(in) :: channels(:)
    double precision,intent(in) :: free(:)
    double precision,intent(out) :: metric(:),matrix(:)
    integer :: c,i,j,r
    metric=-2d0; matrix=-2d0*free
    do c=1,size(channels)
      associate(ch=>channels(c))
        do j=1,size(ch%configs,2)
          do i=1,size(ch%configs,1)
            r=abs(ch%configs(i,j))
            if(r==0) cycle
            metric(r)=metric(r)+ch%metric(i,i)
            matrix(r)=matrix(r)+ch%hamiltonian(i,i)+ch%shift(j)*ch%metric(i,i)
          end do
        end do
      end associate
    end do
  end subroutine

  subroutine sector_poles(channels,free,U,omega,vertex)
! The projected physical propagator is (omega*metric-matrix)^(-1).
! Generalized eigenvectors are metric-normalized for the residues.
    type(FRPA_channel),intent(in) :: channels(:)
    double precision,intent(in) :: free(:),U(:,:)
    double precision,allocatable,intent(out) :: omega(:),vertex(:,:)
    double precision,allocatable :: metric(:,:),matrix(:,:),work(:)
    double precision :: phase
    integer :: n,i,j,r,s,c,k,info,lwork

    n=size(free)
    allocate(metric(n,n),matrix(n,n),omega(n))
    metric=0d0; matrix=0d0
    do r=1,n
      metric(r,r)=-2d0
      matrix(r,r)=-2d0*free(r)
    end do
    do c=1,size(channels)
      associate(ch=>channels(c))
        do k=1,size(ch%configs,2)
          do j=1,size(ch%configs,1)
            s=ch%configs(j,k)
            if(s==0) cycle
            do i=1,size(ch%configs,1)
              r=ch%configs(i,k)
              if(r==0) cycle
              phase=dble(sign(1,r)*sign(1,s))
              metric(abs(r),abs(s))=metric(abs(r),abs(s))+phase*ch%metric(i,j)
              matrix(abs(r),abs(s))=matrix(abs(r),abs(s))+phase*(ch%hamiltonian(i,j)+ch%shift(k)*ch%metric(i,j))
            end do
          end do
        end do
      end associate
    end do
    allocate(work(1))
    call dsygv(1,'V','U',n,matrix,n,metric,n,omega,work,-1,info)
    if(info/=0) error stop 'FRPA: failed to query generalized eigensolver'
    lwork=max(1,int(work(1))); deallocate(work); allocate(work(lwork))
    call dsygv(1,'V','U',n,matrix,n,metric,n,omega,work,lwork,info)
    if(info>n) error stop 'FRPA: physical Faddeev metric is not positive definite'
    if(info/=0) error stop 'FRPA: physical Faddeev eigensolver failed'
    vertex=matmul(transpose(matrix),U)

  end subroutine

  subroutine run_faddeev_pairs(backward,nOrb,nC,nO,nR,eHF,qp,ERI,coupling,static3,external_gw,iterative,nroots, &
                               reg,eta,kernels,Teh,Thh,Tee,ePA,ZPA,external_energies,ccd_ph,ccd_pp)
    include 'parameters.h'
    logical,intent(in) :: backward,static3,external_gw,iterative,reg
    integer,intent(in) :: nOrb,nC,nO,nR,coupling,nroots
    double precision,intent(in) :: eHF(nOrb),qp(nOrb),ERI(nOrb,nOrb,nOrb,nOrb),eta
    type(FRPA_channel),optional,intent(in) :: kernels(3)
    double precision,optional,intent(in) :: ccd_ph(:,:),ccd_pp(:,:)
    double precision,optional,intent(in) :: Teh(:,:),Thh(:,:),Tee(:,:),external_energies(nOrb)
    double precision,intent(out) :: ePA(nOrb),ZPA(nOrb)
    type(Faddeev_operator) :: op
    double precision,allocatable :: H(:,:),omega(:),Z(:),residual(:)
    integer :: iterations,info,j
    character(4) :: label

    label='FTDA'
    if(backward) label='FRPA'
    write(*,*) 'Restricted ',label
    select case(coupling)
    case(2,3)
      write(*,'(A,I1,A)') ' Faddeev with ADC(',coupling,') couplings'
    case(4)
      write(*,*) 'Faddeev: static dressed source approximation'
    case(5)
      write(*,*) 'Faddeev: CCD-amplitude couplings'
    end select
    if(static3) write(*,*) label//': third-order static self-energy enabled'
    call build_FRPA(nOrb,nC,nO,nR,eHF,ERI,backward,H,op,qp,kernels,coupling,static3,external_gw,iterative,reg,eta,Teh,Thh,Tee, &
                    external_energies,ccd_ph,ccd_pp)
    if(iterative) then
      allocate(omega(nroots),Z(nroots),residual(nroots))
      call davidson_roots(op,nO,nroots,omega,Z,residual,iterations,info)
      if(info/=0) error stop 'Faddeev Davidson did not converge'
      write(*,*) 'Faddeev-Parquet IP roots: root, IP (eV), pole (eV), 1p weight Z, residual (Ha)'
      do j=1,nroots
        write(*,'(I7,3F20.10,ES14.4)') j,-omega(j)*HaToeV,omega(j)*HaToeV,Z(j),residual(j)
      end do
      ePA=eHF; ZPA=0d0
    else
      call Faddeev_spectrum(label,nOrb,nC,nR,op%n,eHF,H,ePA,ZPA)
    end if
  end subroutine

  subroutine ccd_sources(t,Teh,Thh)
    double precision,intent(in) :: t(:,:,:,:)
    double precision,allocatable,intent(out) :: Teh(:,:),Thh(:,:)
    integer :: no,nv,i,j,a,b,r,s
    no=2*size(t,1); nv=2*size(t,3)
    allocate(Teh(no*nv,no*nv),Thh(no*(no-1)/2,nv*(nv-1)/2))
    do j=1,no
      do b=1,nv
        do i=1,no
          do a=1,nv
            Teh((j-1)*nv+b,(i-1)*nv+a)=amplitude(i,j,a,b)
          end do
        end do
      end do
    end do
    r=0
    do i=1,no
      do j=i+1,no
        r=r+1; s=0
        do a=1,nv
          do b=a+1,nv
            s=s+1
            Thh(r,s)=amplitude(i,j,a,b)
          end do
        end do
      end do
    end do
  contains
    double precision function amplitude(i,j,a,b) result(value)
      integer,intent(in) :: i,j,a,b
      value=0d0
      if(mod(i-a,2)==0.and.mod(j-b,2)==0) value=t((i+1)/2,(j+1)/2,(a+1)/2,(b+1)/2)
      if(mod(i-b,2)==0.and.mod(j-a,2)==0) value=value-t((i+1)/2,(j+1)/2,(b+1)/2,(a+1)/2)
    end function
  end subroutine

end module

subroutine R_FRPA(nOrb,nC,nO,nR,eHF,ERI,ePA,ZPA)
  use R_FRPA_operator, only: build_FRPA
  use TDA_Faddeev_IP_operator, only: Faddeev_operator
  implicit none
  integer,intent(in) :: nOrb,nC,nO,nR
  double precision,intent(in) :: eHF(nOrb),ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(out) :: ePA(nOrb),ZPA(nOrb)
  double precision,allocatable :: H(:,:)
  type(Faddeev_operator) :: op
  call build_FRPA(nOrb,nC,nO,nR,eHF,ERI,.true.,H,op)
  write(*,*) 'Restricted FRPA with ADC(3) couplings and third-order static self-energy'
  call Faddeev_spectrum('FRPA',nOrb,nC,nR,size(H,1),eHF,H,ePA,ZPA)
end subroutine

subroutine R_FRPA_BSE(nOrb,nC,nO,nR,nS,nhs,nht,nes,net,eHF,eGW,ERI,Wehs,Weht,Whhs,Whht,Wees,Weet, &
                      Behs,Beht,Bpps,Bppt,Teh,Thh,Tee,coupling,static3,external_gw,iterative,nroots,reg,eta,ePA,ZPA)
  use R_FRPA_operator, only: FRPA_channel, run_faddeev_pairs
  implicit none
  integer,intent(in) :: nOrb,nC,nO,nR,nS,nhs,nht,nes,net,coupling,nroots
  double precision,intent(in) :: eHF(nOrb),eGW(nOrb),ERI(nOrb,nOrb,nOrb,nOrb),eta
  double precision,intent(in) :: Wehs(nS,nS),Weht(nS,nS),Whhs(nhs,nhs),Whht(nht,nht),Wees(nes,nes),Weet(net,net)
  double precision,intent(in) :: Behs(nS,nS),Beht(nS,nS),Bpps(nes,nhs),Bppt(net,nht)
  double precision,intent(in) :: Teh(4*nS,4*nS),Thh(nhs+3*nht,nes+3*net),Tee(nhs+3*nht,nes+3*net)
  logical,intent(in) :: static3,external_gw,iterative,reg
  double precision,intent(out) :: ePA(nOrb),ZPA(nOrb)
  type(FRPA_channel) :: kernels(3)
  double precision,allocatable :: Bpp(:,:),unused(:,:)

  allocate(kernels(1)%hamiltonian(4*nS,4*nS),kernels(1)%metric(4*nS,4*nS))
  allocate(kernels(2)%hamiltonian(nhs+3*nht,nhs+3*nht),kernels(3)%hamiltonian(nes+3*net,nes+3*net))
  allocate(Bpp(nhs+3*nht,nes+3*net),unused(nhs+3*nht,nes+3*net))
  call R_TDA_Faddeev_lift(nOrb,nC,nO,nR,nS,nhs,nht,nes,net,Wehs,Weht,Whhs,Whht,Wees,Weet, &
                         kernels(1)%hamiltonian,kernels(2)%hamiltonian,kernels(3)%hamiltonian)
  call R_TDA_Faddeev_lift_sources(nOrb,nC,nO,nR,nS,nhs,nht,nes,net,Behs,Beht, &
                                 transpose(Bpps),transpose(Bppt),transpose(Bpps),transpose(Bppt), &
                                 kernels(1)%metric,Bpp,unused)
  kernels(3)%metric=transpose(Bpp)
  call run_faddeev_pairs(.true.,nOrb,nC,nO,nR,eHF,eGW,ERI,merge(4,coupling,coupling==0),static3,external_gw, &
                         iterative,nroots,reg,eta,kernels,Teh,Thh,Tee,ePA,ZPA)
end subroutine

subroutine R_Faddeev_ccd_sources(no,nv,t,Teh,Thh)
  use R_FRPA_operator, only: ccd_sources
  implicit none
  integer,intent(in) :: no,nv
  double precision,intent(in) :: t(no,no,nv,nv)
  double precision,intent(out) :: Teh(4*no*nv,4*no*nv),Thh(no*(2*no-1),nv*(2*nv-1))
  double precision,allocatable :: ph(:,:),pp(:,:)
  call ccd_sources(t,ph,pp)
  Teh=ph; Thh=pp
end subroutine
