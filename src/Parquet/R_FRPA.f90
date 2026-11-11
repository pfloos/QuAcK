module R_FRPA_operator
  use TDA_Faddeev_IP_operator, only: Faddeev_operator, initialize
  implicit none
  private
  public :: build_FRPA
contains

  subroutine build_FRPA(nOrb,nC,nO,nR,eHF,ERI,backward,H)
! Frozen-HF FRPA in the alpha charged determinant block.
! The backward flag is internal and permits testing the FTDA limit.
    integer,intent(in) :: nOrb,nC,nO,nR
    double precision,intent(in) :: eHF(nOrb),ERI(nOrb,nOrb,nOrb,nOrb)
    logical,intent(in) :: backward
    double precision,allocatable,intent(out) :: H(:,:)
    type(Faddeev_operator) :: op
    integer :: nocc_spin,nv,ns,nh,ne,i,j,a,b,k,l,c,d,r,s,p
    double precision :: mu
    integer,allocatable :: holes(:,:),particles(:,:)
    double precision,allocatable :: e(:),Aph(:,:),Bph(:,:),Cee(:,:),Dhh(:,:),Bpp(:,:)
    double precision,allocatable :: X(:,:),Y(:,:),w(:),v(:),Gph(:,:),Kph(:,:),Ghh(:,:),Khh(:,:),Gee(:,:),Kee(:,:)
    double precision,allocatable :: em(:),ep(:),Um(:,:),Up(:,:),sigma(:,:)

    nocc_spin=2*(nO-nC); nv=2*(nOrb-nR-nO)
    if(nocc_spin<2.or.nv<2) error stop 'FRPA requires active occupied and virtual orbitals'
    ns=nocc_spin*nv; nh=nocc_spin*(nocc_spin-1)/2; ne=nv*(nv-1)/2
    allocate(e(2*nOrb),holes(2,nh),particles(2,ne),Aph(ns,ns),Bph(ns,ns),Cee(ne,ne),Dhh(nh,nh),Bpp(ne,nh))
    do p=1,2*nOrb
      e(p)=eHF((p+1)/2)
    end do
    mu=0.5d0*(eHF(nO)+eHF(nO+1))

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

! Initialize the shared configuration maps and ADC(3) vertices before adding
! the free pair energies. The dressed source arguments are unused.
    call initialize(op,.false.,.false.,0d0,1,nOrb,2*nOrb,2*nC,2*nO,2*nR,ns,nh,ne,e,e,ERI, &
                    Aph,Dhh,Cee,.false.,Aph,transpose(Bpp),transpose(Bpp))
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
    if(.not.backward) then
      Bph=0d0
      Bpp=0d0
    end if

    call solve_RPA('ph',Aph,Bph,Aph,X,Y,w,v)
    call pair_pencil(X,w,Gph,Kph)
! Undo the pp chemical-potential shift before adding the spectator energies.
    call solve_RPA('pp',Cee,Bpp,Dhh,X,Y,w,v)
    call pair_pencil(X,w+2d0*mu,Gee,Kee)
    call pair_pencil(Y,v+2d0*mu,Ghh,Khh)

    call sector_poles(op%im,op%hh,op%free_m,e(2*nO+1:2*(nOrb-nR)), &
                      e(2*nC+1:2*nO),Gph,-Kph,Ghh,Khh,.true.,op%Um,em,Um)
! Put the identical-particle indices first in the addition map too.
    call sector_poles(reshape(op%ip,[nv,nv,nocc_spin],order=[3,1,2]),op%ee,op%free_p,e(2*nC+1:2*nO), &
                      e(2*nO+1:2*(nOrb-nR)),Gph,Kph,Gee,Kee,.false.,op%Up,ep,Up)

    allocate(H(op%n,op%n),sigma(nOrb,nOrb))
    H=0d0
    call TDA_Faddeev_static3(nOrb,2*nOrb,nOrb,2*nC,2*nO,2*nR,e,ERI,sigma)
    H(:nOrb,:nOrb)=sigma
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

  subroutine sector_poles(map,pair,free,spectator,like_energy,Gph,Kph,Gpp,Kpp,removal,U,omega,vertex)
! With R_i the three channel factors of Degroote et al. Eq. (29),
! (I+R_i)^(-1)=(omega-E_free)^(-1) X_i^(-T)(omega-D_i)X_i^(-1).
! Sum the three equations for the physical amplitude and antisymmetrize.
! The metric is the projected sum of the inverse channel norms minus 2I.
! The Hamiltonian is the sum of the matching pole-energy terms minus 2E_free.
! The resulting propagator is (omega*metric-matrix)^(-1), so the
! generalized eigenvectors must be normalized in metric for the residues.
    integer,intent(in) :: map(:,:,:),pair(:,:)
    double precision,intent(in) :: free(:),spectator(:),like_energy(:),Gph(:,:),Kph(:,:),Gpp(:,:),Kpp(:,:),U(:,:)
    logical,intent(in) :: removal
    double precision,allocatable,intent(out) :: omega(:),vertex(:,:)
    double precision,allocatable :: metric(:,:),matrix(:,:),work(:)
    integer,allocatable :: configs(:)
    integer :: n,nlike,nother,nph,npp,i,j,a,b,r,m,info,lwork

    n=size(free); nlike=size(map,1); nother=size(map,3)
    nph=size(Gph,1); npp=size(Gpp,1)
    allocate(metric(n,n),matrix(n,n),omega(n),configs(max(nph,npp)))
    metric=0d0; matrix=0d0
    do r=1,n
      metric(r,r)=-2d0
      matrix(r,r)=-2d0*free(r)
    end do
! The two ph embeddings are related by exchange. In an ordered
! determinant basis their sum is obtained by using both orientations.
    do j=1,nlike
      do i=1,nlike
        do a=1,nother
          m=(a-1)*nlike+i
          if(removal) m=(i-1)*nother+a
          configs(m)=map(i,j,a)
        end do
      end do
      call add_channel(nph,configs,Gph,Kph,like_energy(j))
    end do
    do a=1,nother
      do i=1,nlike
        do b=i+1,nlike
          configs(pair(i,b))=map(i,b,a)
        end do
      end do
      call add_channel(npp,configs,Gpp,Kpp,-spectator(a))
    end do

    allocate(work(1))
    call dsygv(1,'V','U',n,matrix,n,metric,n,omega,work,-1,info)
    if(info/=0) error stop 'FRPA: failed to query generalized eigensolver'
    lwork=max(1,int(work(1))); deallocate(work); allocate(work(lwork))
    call dsygv(1,'V','U',n,matrix,n,metric,n,omega,work,lwork,info)
    if(info>n) error stop 'FRPA: physical Faddeev metric is not positive definite'
    if(info/=0) error stop 'FRPA: physical Faddeev eigensolver failed'
    vertex=matmul(transpose(matrix),U)

  contains

    subroutine add_channel(nc,idx,G,K,shift)
      integer,intent(in) :: nc,idx(nc)
      double precision,intent(in) :: G(nc,nc),K(nc,nc),shift
      integer :: p,q,r,s
      double precision :: phase
      do q=1,nc
        s=abs(idx(q))
        if(s==0) cycle
        do p=1,nc
          r=abs(idx(p))
          if(r==0) cycle
          phase=dble(sign(1,idx(p))*sign(1,idx(q)))
          metric(r,s)=metric(r,s)+phase*G(p,q)
          matrix(r,s)=matrix(r,s)+phase*(K(p,q)+shift*G(p,q))
        end do
      end do
    end subroutine

  end subroutine

end module

subroutine R_FRPA(nOrb,nC,nO,nR,eHF,ERI,ePA,ZPA)
  use R_FRPA_operator, only: build_FRPA
  implicit none
  integer,intent(in) :: nOrb,nC,nO,nR
  double precision,intent(in) :: eHF(nOrb),ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(out) :: ePA(nOrb),ZPA(nOrb)
  double precision,allocatable :: H(:,:)
  call build_FRPA(nOrb,nC,nO,nR,eHF,ERI,.true.,H)
  write(*,*) 'Restricted FRPA with ADC(3) couplings and third-order static self-energy'
  call Faddeev_spectrum('FRPA',nOrb,nC,nR,size(H,1),eHF,H,ePA,ZPA)
end subroutine
