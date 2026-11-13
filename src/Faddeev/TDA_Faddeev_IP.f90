module TDA_Faddeev_IP_operator
! Matrix-free action of the complete physical TDA Faddeev matrix.
! Only pair kernels, endpoint vertices and O(N) configuration maps are stored.
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: Faddeev_operator,initialize,apply,davidson,davidson_roots,FRPA_channel,apply_channels
  type FRPA_channel
    integer,allocatable :: configs(:,:)
    double precision,allocatable :: metric(:,:),hamiltonian(:,:),shift(:)
  end type
  type Faddeev_operator
    type(FRPA_channel),allocatable :: removal(:),addition(:)
    double precision,allocatable :: metric_diag(:)
    integer :: n,nphys,nm,np,no,nv,ns,noo,nvv,nC,nOccEnd
    integer,allocatable :: minus(:,:),plus(:,:),im(:,:,:),ip(:,:,:),hh(:,:),ee(:,:),phys(:)
    double precision,allocatable :: f(:),free_m(:),free_p(:),diag(:),Um(:,:),Up(:,:),Weh(:,:),Whh(:,:),Wee(:,:)
    double precision,allocatable :: sigma_static(:,:)
    double precision,allocatable :: phx(:,:),phy(:,:),ppx(:,:),ppy(:,:)
  end type
contains
  subroutine initialize(op,adc2,reg,eta,block,nERI,nOrb,nC,nO,nR,ns,noo,nvv,eHF,eQP,ERI,Weh,Whh,Wee, &
    dressed,Teh,Thh,Tee)
    type(Faddeev_operator),intent(out) :: op
    logical,intent(in) :: adc2,reg
    integer,intent(in) :: block,nERI,nOrb,nC,nO,nR,ns,noo,nvv
    double precision,intent(in) :: eta,eHF(nOrb),eQP(nOrb),ERI(nERI,nERI,nERI,nERI)
    double precision,intent(in) :: Weh(ns,ns),Whh(noo,noo),Wee(nvv,nvv)
  logical,intent(in) :: dressed
  double precision,intent(in) :: Teh(nS,nS),Thh(nOO,nVV),Tee(nOO,nVV)
    integer :: i,j,a,b,r,s,q,p,ii,jj,aa,bb
    if(nC<0.or.nR<0.or.nO<=nC.or.nOrb-nR<=nO) error stop 'Faddeev IP: invalid active space'
    op%no=nO-nC; op%nv=nOrb-nR-nO
    if(ns/=op%no*op%nv.or.noo/=op%no*(op%no-1)/2.or.nvv/=op%nv*(op%nv-1)/2) &
      error stop 'Faddeev IP: incomplete pair spaces'
    if(block/=0.and.block/=1) error stop 'Faddeev IP: invalid spin block'
    if(block==1) then
      if(mod(nOrb,2)/=0.or.mod(nC,2)/=0.or.mod(nO,2)/=0.or.mod(nR,2)/=0) &
        error stop 'Faddeev IP: alpha block requires paired spin orbitals'
    end if
    call check_kernel(ns,Weh)
    call check_kernel(noo,Whh)
    call check_kernel(nvv,Wee)
    op%nphys=nOrb
    if(block==1) op%nphys=nOrb/2
    op%ns=ns; op%noo=noo; op%nvv=nvv; op%nC=nC; op%nOccEnd=nO
    allocate(op%minus(3,noo*op%nv),op%plus(3,nvv*op%no),op%phys(op%nphys))
    allocate(op%im(op%no,op%no,op%nv),op%ip(op%no,op%nv,op%nv), &
             op%hh(op%no,op%no),op%ee(op%nv,op%nv))
    op%im=0; op%ip=0; op%hh=0; op%ee=0; op%nm=0; op%np=0; s=0
    do i=nC+1,nO
      ii=i-nC
      do j=i+1,nO
        jj=j-nC; s=s+1; op%hh(ii,jj)=s
        do a=nO+1,nOrb-nR
          if(block==1) then
            if(spin(i)+spin(j)-spin(a)/=1) cycle
          end if
          aa=a-nO; op%nm=op%nm+1; r=op%nm
          op%minus(:,r)=[i,j,a]; op%im(ii,jj,aa)=r; op%im(jj,ii,aa)=-r
        end do
      end do
    end do
    s=0
    do a=1,op%nv
      do b=a+1,op%nv
        s=s+1; op%ee(a,b)=s
      end do
    end do
    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        do b=a+1,nOrb-nR
          if(block==1) then
            if(spin(a)+spin(b)-spin(i)/=1) cycle
          end if
          ii=i-nC; aa=a-nO; bb=b-nO; op%np=op%np+1; r=op%np
          op%plus(:,r)=[i,a,b]; op%ip(ii,aa,bb)=r; op%ip(ii,bb,aa)=-r
        end do
      end do
    end do
    op%n=op%nphys+op%nm+op%np
    allocate(op%f(op%nphys),op%free_m(op%nm),op%free_p(op%np),op%diag(op%n), &
             op%Um(op%nm,op%nphys),op%Up(op%np,op%nphys))
    op%Weh=Weh; op%Whh=Whh; op%Wee=Wee
    do q=1,op%nphys
      p=q
      if(block==1) p=2*q-1
      op%phys(q)=p; op%f(q)=eHF(p)
    end do
    call G_TDA_Faddeev_vertices(adc2,reg,eta,nERI,nOrb,nC,nO,nR,op%nphys,op%nm,op%np, &
                                op%phys,op%minus,op%plus,eQP,ERI,op%Um,op%Up, &
    dressed,Teh,Thh,Tee)
    op%diag(1:op%nphys)=op%f
    do r=1,op%nm
      i=op%minus(1,r); j=op%minus(2,r); a=op%minus(3,r)
      ii=i-nC; jj=j-nC; aa=a-nO
      op%free_m(r)=eQP(i)+eQP(j)-eQP(a)
      s=op%hh(ii,jj)
      op%diag(op%nphys+r)=op%free_m(r)-Whh(s,s) &
                           -Weh((ii-1)*op%nv+aa,(ii-1)*op%nv+aa) &
                           -Weh((jj-1)*op%nv+aa,(jj-1)*op%nv+aa)
    end do
    do r=1,op%np
      i=op%plus(1,r); a=op%plus(2,r); b=op%plus(3,r)
      ii=i-nC; aa=a-nO; bb=b-nO
      op%free_p(r)=eQP(a)+eQP(b)-eQP(i)
      s=op%ee(aa,bb)
      op%diag(op%nphys+op%nm+r)=op%free_p(r)+Wee(s,s) &
                           +Weh((ii-1)*op%nv+aa,(ii-1)*op%nv+aa) &
                           +Weh((ii-1)*op%nv+bb,(ii-1)*op%nv+bb)
    end do
    if(.not.all(ieee_is_finite(op%diag)).or..not.all(ieee_is_finite(op%Um)).or. &
       .not.all(ieee_is_finite(op%Up))) error stop 'Faddeev IP: nonfinite operator'
    allocate(op%phx(ns,max(op%no,op%nv)),op%phy(ns,max(op%no,op%nv)), &
             op%ppx(max(1,noo,nvv),max(op%no,op%nv)),op%ppy(max(1,noo,nvv),max(op%no,op%nv)))
  end subroutine

  subroutine check_kernel(n,W)
    integer,intent(in) :: n
    double precision,intent(in) :: W(n,n)
    if(n==0) return
    if(.not.all(ieee_is_finite(W))) error stop 'Faddeev IP: nonfinite kernel'
    if(maxval(abs(W-transpose(W)))>1d-9*max(1d0,maxval(abs(W)))) &
      error stop 'Faddeev IP requires symmetric static pair kernels'
  end subroutine

  integer function spin(i)
    integer,intent(in) :: i
    spin=1-2*mod(i-1,2)
  end function

  subroutine apply(op,x,y)
    type(Faddeev_operator),intent(inout) :: op
    double precision,intent(in) :: x(op%n)
    double precision,intent(out) :: y(op%n)
    integer :: m,p,i,j,a,b,r,s,idx
    m=op%nphys; p=m+op%nm
    y(1:m)=op%f*x(1:m)
    if(allocated(op%sigma_static)) call dgemv('N',m,m,1d0,op%sigma_static,m,x,1,1d0,y,1)
    if(op%nm>0) then
      call dgemv('T',op%nm,m,1d0,op%Um,op%nm,x(m+1),1,1d0,y,1)
      y(m+1:p)=op%free_m*x(m+1:p)
      call dgemv('N',op%nm,m,1d0,op%Um,op%nm,x,1,1d0,y(m+1),1)
    end if
    if(op%np>0) then
      call dgemv('T',op%np,m,1d0,op%Up,op%np,x(p+1),1,1d0,y,1)
      y(p+1:)=op%free_p*x(p+1:)
      call dgemv('N',op%np,m,1d0,op%Up,op%np,x,1,1d0,y(p+1),1)
    end if
    if(allocated(op%removal)) then
      y(m+1:p)=y(m+1:p)-3d0*op%free_m*x(m+1:p)
      y(p+1:)=y(p+1:)-3d0*op%free_p*x(p+1:)
      call apply_channels(op%removal,x(m+1:p),y(m+1:p),.false.)
      call apply_channels(op%addition,x(p+1:),y(p+1:),.false.)
      return
    end if
! Removal: -Whh on the ordered pair, -Weh on each of the two eh pairs.
    op%phx=0d0; op%ppx=0d0
    do j=1,op%no
      do i=1,op%no
        do a=1,op%nv
          idx=op%im(i,j,a)
          if(idx/=0) op%phx((i-1)*op%nv+a,j)=dble(sign(1,idx))*x(m+abs(idx))
        end do
      end do
    end do
    call dgemm('N','N',op%ns,op%no,op%ns,1d0,op%Weh,op%ns,op%phx,op%ns,0d0,op%phy,op%ns)
    do r=1,op%nm
      i=op%minus(1,r)-op%nC; j=op%minus(2,r)-op%nC; a=op%minus(3,r)-op%nOccEnd
      op%ppx(op%hh(i,j),a)=x(m+r)
      y(m+r)=y(m+r)-op%phy((i-1)*op%nv+a,j)+op%phy((j-1)*op%nv+a,i)
    end do
    if(op%noo>0) then
      call dgemm('N','N',op%noo,op%nv,op%noo,1d0,op%Whh,op%noo,op%ppx,size(op%ppx,1), &
                  0d0,op%ppy,size(op%ppy,1))
      do r=1,op%nm
        i=op%minus(1,r)-op%nC; j=op%minus(2,r)-op%nC; a=op%minus(3,r)-op%nOccEnd
        y(m+r)=y(m+r)-op%ppy(op%hh(i,j),a)
      end do
    end if
! Addition: +Wee and the antisymmetrized sum of the two eh interactions.
    op%phx=0d0; op%ppx=0d0
    do b=1,op%nv
      do i=1,op%no
        do a=1,op%nv
          idx=op%ip(i,a,b)
          if(idx/=0) op%phx((i-1)*op%nv+a,b)=dble(sign(1,idx))*x(p+abs(idx))
        end do
      end do
    end do
    call dgemm('N','N',op%ns,op%nv,op%ns,1d0,op%Weh,op%ns,op%phx,op%ns,0d0,op%phy,op%ns)
    do r=1,op%np
      i=op%plus(1,r)-op%nC; a=op%plus(2,r)-op%nOccEnd; b=op%plus(3,r)-op%nOccEnd
      op%ppx(op%ee(a,b),i)=x(p+r)
      y(p+r)=y(p+r)+op%phy((i-1)*op%nv+a,b)-op%phy((i-1)*op%nv+b,a)
    end do
    if(op%nvv>0) then
      call dgemm('N','N',op%nvv,op%no,op%nvv,1d0,op%Wee,op%nvv,op%ppx,size(op%ppx,1), &
                  0d0,op%ppy,size(op%ppy,1))
      do r=1,op%np
        i=op%plus(1,r)-op%nC; a=op%plus(2,r)-op%nOccEnd; b=op%plus(3,r)-op%nOccEnd
        y(p+r)=y(p+r)+op%ppy(op%ee(a,b),i)
      end do
    end if
  end subroutine

  subroutine apply_channels(channels,x,y,metric_only)
    type(FRPA_channel),intent(in) :: channels(:)
    double precision,intent(in) :: x(:)
    double precision,intent(inout) :: y(:)
    logical,intent(in) :: metric_only
    integer :: c,j,k,r,n
    double precision,allocatable :: u(:),v(:),g(:)
    do c=1,size(channels)
      associate(ch=>channels(c))
        n=size(ch%configs,1)
        allocate(u(n),v(n),g(n))
        do j=1,size(ch%configs,2)
          u=0d0
          do k=1,n
            r=ch%configs(k,j)
            if(r/=0) u(k)=sign(1,r)*x(abs(r))
          end do
          call dgemv('N',n,n,1d0,ch%metric,n,u,1,0d0,g,1)
          v=g
          if(.not.metric_only) then
            v=ch%shift(j)*g
            call dgemv('N',n,n,1d0,ch%hamiltonian,n,u,1,1d0,v,1)
          end if
          do k=1,n
            r=ch%configs(k,j)
            if(r/=0) y(abs(r))=y(abs(r))+sign(1,r)*v(k)
          end do
        end do
        deallocate(u,v,g)
      end associate
    end do
  end subroutine

  subroutine apply_metric(op,x,y)
    type(Faddeev_operator),intent(in) :: op
    double precision,intent(in) :: x(op%n)
    double precision,intent(out) :: y(op%n)
    integer :: m,p
    y=x
    if(.not.allocated(op%removal)) return
    m=op%nphys; p=m+op%nm
    y(m+1:)=-2d0*x(m+1:)
    call apply_channels(op%removal,x(m+1:p),y(m+1:p),.true.)
    call apply_channels(op%addition,x(p+1:),y(p+1:),.true.)
  end subroutine

  subroutine harmonic(m,V,AV,shift,coeff)
! Harmonic extraction: B y = mu C y, with B=V^T(A-shift)V,
! C=((A-shift)V)^T((A-shift)V), theta=shift+1/mu.
! Whiten the positive part of C rather than factor an indefinite B.
! Null directions of C are retained: they can contain a converged root.
    integer,intent(in) :: m
    double precision,intent(in) :: V(:,:),AV(:,:),shift
    double precision,intent(out) :: coeff(:,:)
    double precision,allocatable :: W(:,:),B(:,:),C(:,:),Z(:,:),D(:,:),eval(:),mu(:)
    double precision :: cutoff,norm
    integer :: i,k,rank,nn
    allocate(W(size(AV,1),size(AV,2)),B(m,m),C(m,m))
    W=AV-shift*V
    call dgemm('T','N',m,m,size(V,1),1d0,V,size(V,1),W,size(W,1),0d0,B,m)
    B=0.5d0*(B+transpose(B))
    call dgemm('T','N',m,m,size(W,1),1d0,W,size(W,1),W,size(W,1),0d0,C,m)
    allocate(eval(m))
    call diagonalize_matrix(m,C,eval)
    cutoff=max(1d-28,1d-14*maxval(abs(eval)))
    rank=count(eval>cutoff); nn=m-rank
    coeff=0d0
    if(nn>0) coeff(:,1:nn)=C(:,1:nn)
    if(rank>0) then
      allocate(Z(m,rank),D(rank,rank),mu(rank))
      do i=1,rank
        Z(:,i)=C(:,nn+i)/sqrt(eval(nn+i))
      end do
      D=matmul(transpose(Z),matmul(B,Z)); D=0.5d0*(D+transpose(D))
      call diagonalize_matrix(rank,D,mu)
      coeff(:,nn+1:m)=matmul(Z,D)
    end if
    do k=1,m
      norm=sqrt(sum(coeff(:,k)**2))
      if(norm<1d-30) error stop 'Faddeev IP: singular harmonic vector'
      coeff(:,k)=coeff(:,k)/norm
    end do
  end subroutine

  subroutine davidson(op,nocc,omega,Z,homo_weight,residual,iterations,info,max_space,max_iterations)
    type(Faddeev_operator),intent(inout) :: op
    integer,intent(in) :: nocc
    double precision,intent(out) :: omega,Z,homo_weight,residual
    integer,intent(out) :: iterations,info
    integer,intent(in),optional :: max_space,max_iterations
    double precision,parameter :: tolerance=1d-8
    double precision,allocatable :: V(:,:),AV(:,:),coef(:,:),x(:,:),ax(:,:), &
                                    theta(:),hw(:),rn(:),score(:),t(:),previous(:,:),keep(:,:)
    integer,allocatable :: homo(:),order(:)
    integer :: n,limit,m,nh,nt,k,j,i,added,oldm,nkeep,chosen,nmv,maxiter,first_occ
    double precision :: shift,fermi,norm,den,overlap
    logical :: converged
    n=op%n; limit=min(n,96)
    maxiter=256
    if(present(max_space)) limit=min(n,max_space)
    if(present(max_iterations)) maxiter=max_iterations
    if(maxiter<1) error stop 'Faddeev IP: invalid iteration limit'
    first_occ=count(op%phys<=op%nC)+1
    if(nocc<first_occ.or.nocc>=op%nphys) error stop 'Faddeev IP: invalid HOMO index'
    nh=count(abs(op%f(first_occ:nocc)-op%f(nocc))<1d-7)
    if(nh<1.or.nh+2>limit) error stop 'Faddeev IP: HOMO degeneracy exceeds Davidson subspace'
    allocate(homo(nh)); j=0
    do i=first_occ,nocc
      if(abs(op%f(i)-op%f(nocc))<1d-7) then
        j=j+1; homo(j)=i
      end if
    end do
    nt=min(limit,max(4,nh+2))
    allocate(V(n,limit),AV(n,limit),coef(limit,limit),x(n,limit),ax(n,limit), &
             theta(limit),hw(limit),rn(limit),score(limit),t(n),previous(n,nh), &
             order(limit),keep(n,nt))
    V=0d0; previous=0d0; m=nh; nmv=0
    do j=1,nh
      V(homo(j),j)=1d0
      previous(:,j)=V(:,j)
    end do
    ! Include nearby occupied one-particle states as additional guesses.
    do i=nocc,first_occ,-1
      if(any(homo==i)) cycle
      if(m>=min(limit,nh+2)) exit
      m=m+1; V(i,m)=1d0
    end do
    do j=1,m
      call apply(op,V(:,j),AV(:,j)); nmv=nmv+1
    end do
    shift=op%f(nocc); fermi=0.5d0*(op%f(nocc)+op%f(nocc+1))
    info=1; omega=shift; Z=0d0; homo_weight=0d0; residual=huge(1d0)
    do iterations=1,maxiter
      call harmonic(m,V(:,1:m),AV(:,1:m),shift,coef(1:m,1:m))
      call dgemm('N','N',n,m,m,1d0,V,n,coef,limit,0d0,x,n)
      call dgemm('N','N',n,m,m,1d0,AV,n,coef,limit,0d0,ax,n)
      do k=1,m
        norm=sqrt(sum(x(:,k)**2)); x(:,k)=x(:,k)/norm; ax(:,k)=ax(:,k)/norm
        theta(k)=dot_product(x(:,k),ax(:,k))
        hw(k)=sum(x(homo,k)**2)
        rn(k)=sqrt(sum((ax(:,k)-theta(k)*x(:,k))**2))
        overlap=sum(matmul(transpose(previous),x(:,k))**2)
        ! HOMO spectral character is primary; overlap breaks near ties.
        score(k)=hw(k)+0.02d0*overlap
        if(theta(k)>=fermi) score(k)=-1d0-abs(theta(k)-shift)
        order(k)=k
      end do
      do i=1,m-1
        do j=i+1,m
          if(score(order(j))>score(order(i))) then
            k=order(i); order(i)=order(j); order(j)=k
          end if
        end do
      end do
      converged=.true.
      do j=1,nh
        k=order(j)
        if(hw(k)<1d-10.or.theta(k)>=fermi.or.rn(k)>tolerance) converged=.false.
        previous(:,j)=x(:,k)
      end do
      if(converged) then
        ! Among a split HOMO manifold, report the smallest principal IP.
        chosen=order(1)
        do j=2,nh
          k=order(j)
          if(theta(k)>theta(chosen)) chosen=k
        end do
        omega=theta(chosen); Z=sum(x(1:op%nphys,chosen)**2)
        homo_weight=hw(chosen); residual=rn(chosen); info=0
        call apply(op,x(:,chosen),t); nmv=nmv+1
        residual=sqrt(sum((t-omega*x(:,chosen))**2))
        if(residual>tolerance) info=1
        exit
      end if
      oldm=m; added=0
      if(m+nt>limit) then
        ! Thick restart retains several candidate poles and their character.
        nkeep=min(nt,m)
        do j=1,nkeep
          keep(:,j)=x(:,order(j))
        end do
        m=0
        do j=1,nkeep
          t=keep(:,j)
          call orthogonalize(n,m,V,t,norm)
          if(norm<1d-10) cycle
          m=m+1; V(:,m)=t/norm
          call apply(op,V(:,m),AV(:,m)); nmv=nmv+1
        end do
      end if
      do j=1,min(nt,oldm)
        k=order(j)
        if(rn(k)<=tolerance) cycle
        t=ax(:,k)-theta(k)*x(:,k)
        do i=1,n
          den=theta(k)-op%diag(i)
          if(abs(den)<0.05d0) den=sign(0.05d0,den)
          t(i)=t(i)/den
        end do
        call orthogonalize(n,m,V,t,norm)
        if(norm<1d-10) cycle
        m=m+1; V(:,m)=t/norm
        call apply(op,V(:,m),AV(:,m)); nmv=nmv+1; added=added+1
        if(m==limit) exit
      end do
      if(added==0) then
        ! An unpreconditioned residual can escape a nearly dependent correction.
        k=order(1); t=ax(:,k)-theta(k)*x(:,k)
        call orthogonalize(n,m,V,t,norm)
        if(norm<1d-12.or.m==limit) exit
        m=m+1; V(:,m)=t/norm
        call apply(op,V(:,m),AV(:,m)); nmv=nmv+1
      end if
      shift=theta(order(1))
    end do
    write(*,'(A,I5,A,I6)') ' Faddeev IP Davidson iterations: ',min(iterations,maxiter),'; matrix-vector products: ',nmv
    if(info/=0) then
      k=order(1); omega=theta(k); Z=sum(x(1:op%nphys,k)**2); homo_weight=hw(k); residual=rn(k)
    end if
    write(*,'(A,ES12.4)') ' Faddeev IP residual (Ha): ',residual
    if(.not.ieee_is_finite(omega).or..not.ieee_is_finite(residual)) info=2
  end subroutine

  subroutine davidson_roots(op,nocc,nroots,omega,Z,residual,iterations,info,max_space,max_iterations)
! Follow the dominant poles of the highest nroots active occupied HF orbitals.
! HF hole vectors seed the search; a one-to-one maximum-weight assignment
! keeps distinct Ritz roots attached to their HF orbital character.
    type(Faddeev_operator),intent(inout) :: op
    integer,intent(in) :: nocc,nroots
    double precision,intent(out) :: omega(nroots),Z(nroots),residual(nroots)
    integer,intent(out) :: iterations,info
    integer,intent(in),optional :: max_space,max_iterations
    double precision,parameter :: tolerance=1d-8
    double precision,allocatable :: V(:,:),AV(:,:),coef(:,:),x(:,:),ax(:,:), &
                                    theta(:),rn(:),weight(:,:),character(:),t(:),keep(:,:)
    double precision,allocatable :: BV(:,:),bx(:,:),projected_metric(:,:),work(:)
    integer,allocatable :: order(:),assigned(:),target(:),reported(:)
    integer :: n,limit,maxiter,nt,m,i,j,k,first_occ,first_target,nh,nkeep,nmv,added,lapack_info,metric_space
    double precision :: fermi,norm,den
    logical :: converged,generalized
    n=op%n; maxiter=256
    first_occ=count(op%phys<=op%nC)+1
    if(nocc<first_occ.or.nocc>=op%nphys) error stop 'Faddeev IP roots: invalid HOMO index'
    if(nroots<1.or.nroots>nocc-first_occ+1) &
      error stop 'Faddeev IP roots: root count exceeds active occupied HF orbitals'
    first_target=nocc-nroots+1
! Retain a whole HF degeneracy at the lower boundary; report nroots poles.
    do while(first_target>first_occ)
      if(abs(op%f(first_target-1)-op%f(nocc-nroots+1))>1d-7) exit
      first_target=first_target-1
    end do
    nh=nocc-first_target+1
    limit=min(n,max(96,4*nh+8))
    if(present(max_space)) limit=min(n,max_space)
    if(present(max_iterations)) maxiter=max_iterations
    if(limit<min(n,nh+2).or.maxiter<1) error stop 'Faddeev IP roots: invalid solver limits'
    fermi=0.5d0*(op%f(nocc)+op%f(nocc+1))
    nt=min(limit-1,max(nh+4,2*nh))
    allocate(V(n,limit),AV(n,limit),coef(limit,limit),x(n,limit),ax(n,limit), &
             theta(limit),rn(limit),weight(nh,limit),character(limit),t(n),keep(n,nt), &
             order(limit),assigned(nh),target(nh),reported(nh))
    generalized=allocated(op%metric_diag)
! Zero-sized metric workspace in the ordinary FTDA eigenproblem.
    metric_space=merge(limit,0,generalized)
    allocate(BV(n,metric_space),bx(n,metric_space), &
             projected_metric(metric_space,metric_space),work(3*metric_space))
    V=0d0; m=nh; nmv=0
    do j=1,nh
      target(j)=nocc-j+1
      V(target(j),j)=1d0
    end do
! Only occupied HF holes are used as initial guesses, never charged satellites.
    do i=first_target-1,first_occ,-1
      if(m>=min(limit,nh+2)) exit
      m=m+1; V(i,m)=1d0
    end do
    do j=1,m
      call apply(op,V(:,j),AV(:,j)); nmv=nmv+1
      if(generalized) call apply_metric(op,V(:,j),BV(:,j))
    end do
    info=1; omega=op%f(nocc); Z=0d0; residual=huge(1d0)
    do iterations=1,maxiter
! Ordinary Ritz extraction allows each HF hole to follow its own energy,
! including deeper ionizations far from the HOMO or the gap midpoint.
      call dgemm('T','N',m,m,n,1d0,V,n,AV,n,0d0,coef,limit)
      coef(1:m,1:m)=0.5d0*(coef(1:m,1:m)+transpose(coef(1:m,1:m)))
      if(generalized) then
        call dgemm('T','N',m,m,n,1d0,V,n,BV,n,0d0,projected_metric,limit)
        call dsygv(1,'V','U',m,coef,limit,projected_metric,limit,theta,work,size(work),lapack_info)
        if(lapack_info/=0) error stop 'FRPA: projected Davidson metric or eigensolver failed'
        call dgemm('N','N',n,m,m,1d0,BV,n,coef,limit,0d0,bx,n)
      else
        call diagonalize_matrix(m,coef(1:m,1:m),theta(1:m))
      end if
      call dgemm('N','N',n,m,m,1d0,V,n,coef,limit,0d0,x,n)
      call dgemm('N','N',n,m,m,1d0,AV,n,coef,limit,0d0,ax,n)
      do k=1,m
        rn(k)=sqrt(sum((ax(:,k)-theta(k)*x(:,k))**2))
        if(generalized) rn(k)=norm2(ax(:,k)-theta(k)*bx(:,k))
        weight(:,k)=x(target,k)**2
        character(k)=sum(weight(:,k))
        if(theta(k)>=fermi) weight(:,k)=-1d0
      end do
      call assign_holes(nh,m,weight(:,1:m),assigned)
      order(1:nh)=assigned
! Keep additional roots with HF target character through thick restarts.
      nkeep=nh
      do k=1,m
        if(any(assigned==k)) cycle
        nkeep=nkeep+1; order(nkeep)=k
      end do
      do i=nh+1,m-1
        do j=i+1,m
          if(character(order(j))>character(order(i))) then
            k=order(i); order(i)=order(j); order(j)=k
          end if
        end do
      end do
      converged=.true.
      do j=1,nh
        k=assigned(j)
        if(theta(k)>=fermi.or.character(k)<1d-10.or.rn(k)>tolerance) converged=.false.
      end do
      reported=assigned
      do i=1,nh-1
        do j=i+1,nh
          if(theta(reported(j))>theta(reported(i))) then
            k=reported(i); reported(i)=reported(j); reported(j)=k
          end if
        end do
      end do
      do j=1,nroots
        k=reported(j); omega(j)=theta(k); Z(j)=sum(x(1:op%nphys,k)**2); residual(j)=rn(k)
      end do
      if(converged) then
        do j=1,nroots
          k=reported(j)
          call apply(op,x(:,k),t); nmv=nmv+1
          if(generalized) then
            call apply_metric(op,x(:,k),bx(:,k))
            residual(j)=norm2(t-omega(j)*bx(:,k))
          else
            residual(j)=norm2(t-omega(j)*x(:,k))
          end if
        end do
        if(all(residual<=tolerance)) then
          info=0
          write(*,*) 'Faddeev HF assignments: root, HF orbital, orbital weight'
          do j=1,nroots
            k=reported(j)
            i=findloc(assigned,k,dim=1)
            write(*,'(2I7,F20.10)') j,target(i),weight(i,k)
          end do
        end if
        exit
      end if
      if(m+nh>limit) then
        nkeep=min(nt,m)
        do j=1,nkeep
          keep(:,j)=x(:,order(j))
        end do
        m=0
        do j=1,nkeep
          t=keep(:,j)
          call orthogonalize(n,m,V,t,norm)
          if(norm<1d-10) cycle
          m=m+1; V(:,m)=t/norm
          call apply(op,V(:,m),AV(:,m)); nmv=nmv+1
          if(generalized) call apply_metric(op,V(:,m),BV(:,m))
        end do
      end if
      added=0
      do j=1,nh
        k=assigned(j)
        if(rn(k)<=tolerance) cycle
        t=ax(:,k)-theta(k)*x(:,k)
        if(generalized) t=ax(:,k)-theta(k)*bx(:,k)
        do i=1,n
          den=theta(k)-op%diag(i)
          if(generalized) den=theta(k)*op%metric_diag(i)-op%diag(i)
          if(abs(den)<0.05d0) den=sign(0.05d0,den)
          t(i)=t(i)/den
        end do
        call orthogonalize(n,m,V,t,norm)
        if(norm<1d-10) cycle
        m=m+1; V(:,m)=t/norm
        call apply(op,V(:,m),AV(:,m)); nmv=nmv+1; added=added+1
        if(generalized) call apply_metric(op,V(:,m),BV(:,m))
        if(m==limit) exit
      end do
      if(added==0) then
        do j=1,nh
          k=assigned(j); t=ax(:,k)-theta(k)*x(:,k)
          if(generalized) t=ax(:,k)-theta(k)*bx(:,k)
          call orthogonalize(n,m,V,t,norm)
          if(norm<1d-12.or.m==limit) cycle
          m=m+1; V(:,m)=t/norm
          call apply(op,V(:,m),AV(:,m)); nmv=nmv+1; added=1
          if(generalized) call apply_metric(op,V(:,m),BV(:,m))
          exit
        end do
        if(added==0) exit
      end if
    end do
    write(*,'(A,I5,A,I6)') ' Faddeev IP Davidson iterations: ',min(iterations,maxiter),'; matrix-vector products: ',nmv
    write(*,'(A,ES12.4)') ' Largest Faddeev IP residual (Ha): ',maxval(residual)
    if(.not.all(ieee_is_finite(omega)).or..not.all(ieee_is_finite(residual))) info=2
  end subroutine

  subroutine assign_holes(nh,m,weight,assigned)
! Rectangular Hungarian assignment maximizes the total HF-hole weight.
! Each HF orbital gets a distinct root, even when several orbitals mix.
    integer,intent(in) :: nh,m
    double precision,intent(in) :: weight(nh,m)
    integer,intent(out) :: assigned(nh)
    double precision :: u(0:nh),v(0:m),best(m),delta,cost
    integer :: owner(0:m),way(m),i,j,j0,j1,i0
    logical :: used(0:m)
    if(nh>m) error stop 'Faddeev IP roots: fewer Ritz vectors than targets'
    u=0d0; v=0d0; owner=0
    do i=1,nh
      owner(0)=i; j0=0; best=huge(1d0); used=.false.
      do
        used(j0)=.true.; i0=owner(j0); delta=huge(1d0); j1=0
        do j=1,m
          if(used(j)) cycle
          cost=-weight(i0,j)-u(i0)-v(j)
          if(cost<best(j)) then
            best(j)=cost; way(j)=j0
          end if
          if(best(j)<delta) then
            delta=best(j); j1=j
          end if
        end do
        do j=0,m
          if(used(j)) then
            u(owner(j))=u(owner(j))+delta; v(j)=v(j)-delta
          else if(j>0) then
            best(j)=best(j)-delta
          end if
        end do
        j0=j1
        if(owner(j0)==0) exit
      end do
      do
        j1=way(j0); owner(j0)=owner(j1); j0=j1
        if(j0==0) exit
      end do
    end do
    assigned=0
    do j=1,m
      if(owner(j)>0) assigned(owner(j))=j
    end do
  end subroutine

  subroutine orthogonalize(n,m,V,t,norm)
    integer,intent(in) :: n,m
    double precision,intent(in) :: V(:,:)
    double precision,intent(inout) :: t(n)
    double precision,intent(out) :: norm
    integer :: pass,j
    do pass=1,2
      do j=1,m
        t=t-dot_product(V(:,j),t)*V(:,j)
      end do
    end do
    norm=sqrt(sum(t**2))
  end subroutine
end module

subroutine G_TDA_Faddeev_IP(adc2,reg,eta,block,nERI,nOrb,nC,nO,nR,ns,noo,nvv, &
                           eHF,eQP,ERI,Weh,Whh,Wee,omega,Z,homo_weight,residual,info, &
    dressed,Teh,Thh,Tee)
  use TDA_Faddeev_IP_operator
  implicit none
  include 'parameters.h'
  logical,intent(in) :: adc2,reg
  integer,intent(in) :: block,nERI,nOrb,nC,nO,nR,ns,noo,nvv
  double precision,intent(in) :: eta,eHF(nOrb),eQP(nOrb),ERI(nERI,nERI,nERI,nERI)
  double precision,intent(in) :: Weh(ns,ns),Whh(noo,noo),Wee(nvv,nvv)
  double precision,intent(out) :: omega,Z,homo_weight,residual
  integer,intent(out) :: info
  logical,intent(in) :: dressed
  double precision,intent(in) :: Teh(nS,nS),Thh(nOO,nVV),Tee(nOO,nVV)
  type(Faddeev_operator) :: op
  integer :: nocc,iterations
  double precision :: start,finish
  call wall_time(start)
  call initialize(op,adc2,reg,eta,block,nERI,nOrb,nC,nO,nR,ns,noo,nvv,eHF,eQP,ERI,Weh,Whh,Wee, &
    dressed,Teh,Thh,Tee)
  call wall_time(finish)
  if(dressed) then
    write(*,*) 'TDA Faddeev-Parquet: static dressed source approximation'
  else if(adc2) then
    write(*,*) 'Matrix-free TDA Faddeev-Parquet principal IP with ADC(2) couplings'
  else
    write(*,*) 'Matrix-free TDA Faddeev-Parquet principal IP with ADC(3) couplings'
  end if
  write(*,*) 'Physical upfolded dimension = ',op%n
  write(*,'(A,F12.3,A)') ' Faddeev operator setup: ',finish-start,' seconds'
  nocc=nO
  if(block==1) nocc=nO/2
  call wall_time(start)
  call davidson(op,nocc,omega,Z,homo_weight,residual,iterations,info)
  call wall_time(finish)
  write(*,'(A,F12.3,A)') ' Faddeev Davidson solve: ',finish-start,' seconds'
  if(info/=0) return
  write(*,'(A,F20.10)') ' Principal IP (eV): ', -omega*HaToeV
  write(*,'(A,F20.10)') ' Principal pole (eV): ', omega*HaToeV
  write(*,'(A,F20.10)') ' Principal pole 1p weight Z: ',Z
  write(*,'(A,F20.10)') ' Principal pole HOMO-subspace weight: ',homo_weight
  if(homo_weight<0.5d0) write(*,*) 'HOMO character is fragmented; compare nearby poles with the full spectrum.'
end subroutine

subroutine G_TDA_Faddeev_IP_roots(nroots,adc2,reg,eta,block,nERI,nOrb,nC,nO,nR,ns,noo,nvv, &
                           eHF,eQP,ERI,Weh,Whh,Wee,omega,Z,residual,info, &
    dressed,Teh,Thh,Tee,static3,gw_external)
  use TDA_Faddeev_IP_operator
  implicit none
  include 'parameters.h'
  logical,intent(in) :: adc2,reg
  integer,intent(in) :: nroots,block,nERI,nOrb,nC,nO,nR,ns,noo,nvv
  double precision,intent(in) :: eta,eHF(nOrb),eQP(nOrb),ERI(nERI,nERI,nERI,nERI)
  double precision,intent(in) :: Weh(ns,ns),Whh(noo,noo),Wee(nvv,nvv)
  double precision,intent(out) :: omega(nroots),Z(nroots),residual(nroots)
  integer,intent(out) :: info
  logical,intent(in) :: dressed,static3,gw_external
  double precision,intent(in) :: Teh(nS,nS),Thh(nOO,nVV),Tee(nOO,nVV)
  type(Faddeev_operator) :: op
  integer :: nocc,iterations,j
  double precision :: start,finish,homo_weight
  call wall_time(start)
  call initialize(op,adc2,reg,eta,block,nERI,nOrb,nC,nO,nR,ns,noo,nvv,eHF,eQP,ERI,Weh,Whh,Wee, &
    dressed,Teh,Thh,Tee)
  if(static3.or.gw_external) then
    allocate(op%sigma_static(op%nphys,op%nphys))
    op%sigma_static=0d0
    if(static3) call TDA_Faddeev_static3(nERI,nOrb,op%nphys,nC,nO,nR,eHF,ERI,op%sigma_static)
    ! Add the diagonal shift to the operator while retaining HF orbital targets.
    do j=1,op%nphys
      if(gw_external) op%sigma_static(j,j)=op%sigma_static(j,j)+eQP(op%phys(j))-eHF(op%phys(j))
      op%diag(j)=op%diag(j)+op%sigma_static(j,j)
    end do
  end if
  call wall_time(finish)
  if(dressed) then
    write(*,*) 'TDA Faddeev-Parquet: static dressed source approximation'
  else if(adc2) then
    write(*,*) 'Matrix-free TDA Faddeev-Parquet IP roots with ADC(2) couplings'
  else
    write(*,*) 'Matrix-free TDA Faddeev-Parquet IP roots with ADC(3) couplings'
  end if
  write(*,*) 'Physical upfolded dimension = ',op%n
  write(*,'(A,F12.3,A)') ' Faddeev operator setup: ',finish-start,' seconds'
  nocc=nO
  if(block==1) nocc=nO/2
  call wall_time(start)
  if(nroots==1) then
    call davidson(op,nocc,omega(1),Z(1),homo_weight,residual(1),iterations,info)
  else
    call davidson_roots(op,nocc,nroots,omega,Z,residual,iterations,info)
  end if
  call wall_time(finish)
  write(*,'(A,F12.3,A)') ' Faddeev Davidson solve: ',finish-start,' seconds'
  if(info/=0) return
  if(nroots==1) then
    write(*,*) 'Faddeev HF assignments: root, HF orbital, orbital weight'
    write(*,'(2I7,F20.10)') 1,nocc,homo_weight
  end if
  write(*,*) 'Faddeev-Parquet IP roots: root, IP (eV), pole (eV), 1p weight Z, residual (Ha)'
  do j=1,nroots
    write(*,'(I7,3F20.10,ES14.4)') j,-omega(j)*HaToeV,omega(j)*HaToeV,Z(j),residual(j)
  end do
end subroutine

subroutine G_TDA_Faddeev_IP_action(adc2,reg,eta,block,nERI,nOrb,nC,nO,nR,ns,noo,nvv,n,nvec, &
                                  eHF,eQP,ERI,Weh,Whh,Wee,x,y,diag, &
    dressed,Teh,Thh,Tee)
! Regression interface: compare the operator with the dense construction.
  use TDA_Faddeev_IP_operator
  implicit none
  logical,intent(in) :: adc2,reg
  integer,intent(in) :: block,nERI,nOrb,nC,nO,nR,ns,noo,nvv,n,nvec
  double precision,intent(in) :: eta,eHF(nOrb),eQP(nOrb),ERI(nERI,nERI,nERI,nERI)
  double precision,intent(in) :: Weh(ns,ns),Whh(noo,noo),Wee(nvv,nvv),x(n,nvec)
  double precision,intent(out) :: y(n,nvec),diag(n)
  logical,intent(in) :: dressed
  double precision,intent(in) :: Teh(nS,nS),Thh(nOO,nVV),Tee(nOO,nVV)
  type(Faddeev_operator) :: op
  integer :: j
  call initialize(op,adc2,reg,eta,block,nERI,nOrb,nC,nO,nR,ns,noo,nvv,eHF,eQP,ERI,Weh,Whh,Wee, &
    dressed,Teh,Thh,Tee)
  if(n/=op%n) error stop 'Faddeev IP action: incorrect dimension'
  diag=op%diag
  do j=1,nvec
    call apply(op,x(:,j),y(:,j))
  end do
end subroutine
