subroutine STEOM_eB_davidson(TDA,symmetric,n,nRoots,tolerance,matrix,factor,diagonal,roots,iterations,residual,info)

! Lowest real roots of the singles matrix or L^T (A+B) L, where A-B = L L^T.
! Cache images of the orthonormal Davidson basis and restart with Ritz vectors.

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  logical,intent(in) :: TDA,symmetric
  integer,intent(in) :: n,nRoots
  double precision,intent(in) :: tolerance,matrix(n,n),factor(n,merge(0,n,TDA)),diagonal(n)
  double precision,intent(out) :: roots(nRoots),residual
  integer,intent(out) :: iterations,info

  integer,parameter :: maxIterations=200
  integer :: limit,m,previous,k,p,seed,lapackInfo
  logical :: available(n)
  double precision :: floor,guessMix,dummy(1,1)
  integer,allocatable :: order(:)
  double precision,allocatable :: basis(:,:),image(:,:),projected(:,:),values(:),work(:)
  double precision,allocatable :: imaginary(:),projectedRight(:,:)
  double precision,allocatable :: vectors(:,:),ritzImage(:,:),errors(:,:),norms(:),trial(:),denominator(:),buffer(:)

  limit=min(n,max(32,6*nRoots+8))
  allocate(basis(n,limit),image(n,limit),projected(limit,limit),values(limit),work(max(1,4*limit)))
  p=merge(0,limit,symmetric)
  allocate(imaginary(p),projectedRight(p,p),order(p))
  allocate(vectors(n,nRoots),ritzImage(n,nRoots),errors(n,nRoots),norms(nRoots),trial(n),denominator(n),buffer(n))
  roots=0d0
  residual=huge(1d0)
  info=1
  floor=sqrt(epsilon(1d0))*max(1d0,maxval(abs(diagonal)))

  available=.true.
  m=0
  guessMix=0.05d0/sqrt(dble(n))
  do k=1,min(n,nRoots+2)
    seed=minloc(diagonal,dim=1,mask=available)
    available(seed)=.false.
    ! Give each starting vector components outside its diagonal seed block.
    do p=1,n
      trial(p)=guessMix*sin(dble(p)*dble(k)*sqrt(2d0))
    end do
    trial(seed)=trial(seed)+1d0
    call append_trial()
  end do

  iterations=0
  if(m < nRoots) return

  do iterations=1,maxIterations
    call dgemm('T','N',m,m,n,1d0,basis,n,image,n,0d0,projected,limit)
    if(symmetric) then
      call dsyev('V','U',m,projected,limit,values,work,size(work),lapackInfo)
    else
      call dgeev('N','V',m,projected,limit,values,imaginary,dummy,1,projectedRight,limit,work,size(work),lapackInfo)
    end if
    if(lapackInfo /= 0 .or. .not.all(ieee_is_finite(values(1:m)))) then
      info=2
      return
    end if

    if(.not.symmetric) then
      order(1:m)=[(k,k=1,m)]
      call quick_sort(values,order,m)
      ! Complex Ritz vectors require a complex subspace; use the dense
      ! nonsymmetric solver when any targeted Ritz pair is complex.
      if(any(imaginary(order(1:nRoots)) /= 0d0)) then
        info=3
        return
      end if
      projected(1:m,1:nRoots)=projectedRight(1:m,order(1:nRoots))
    end if

    roots=values(1:nRoots)
    call dgemm('N','N',n,nRoots,m,1d0,basis,n,projected,limit,0d0,vectors,n)
    call dgemm('N','N',n,nRoots,m,1d0,image,n,projected,limit,0d0,ritzImage,n)
    do k=1,nRoots
      errors(:,k)=ritzImage(:,k)-roots(k)*vectors(:,k)
      norms(k)=norm2(errors(:,k))
    end do
    residual=maxval(norms)
    if(residual <= tolerance) then
      info=0
      return
    end if

    if(m > limit-count(norms > tolerance) .and. limit < n) then
      ! Right Ritz vectors of a general matrix need not be orthogonal.
      if(symmetric) then
        basis(:,1:nRoots)=vectors
        image(:,1:nRoots)=ritzImage
        m=nRoots
      else
        m=0
        do k=1,nRoots
          trial=vectors(:,k)
          call append_trial()
        end do
      end if
    end if

    previous=m
    do k=1,nRoots
      if(m == limit) exit
      if(norms(k) <= tolerance) cycle
      denominator=roots(k)-diagonal
      where(abs(denominator) < floor) denominator=sign(floor,denominator)
      trial=errors(:,k)/denominator
      p=m
      call append_trial()
      if(m == p) then
        trial=errors(:,k)
        call append_trial()
      end if
    end do
    if(m == previous) return
  end do
  iterations=maxIterations

contains

  subroutine append_trial()
    integer :: pass
    double precision :: length,initialLength

    initialLength=norm2(trial)
    do pass=1,2
      if(m > 0) trial=trial-matmul(basis(:,1:m),matmul(trial,basis(:,1:m)))
    end do
    length=norm2(trial)
    if(length <= 100d0*epsilon(1d0)*initialLength) return
    m=m+1
    basis(:,m)=trial/length

    buffer=basis(:,m)
    if(.not.TDA) call dtrmv('L','N','N',n,factor,n,buffer,1)
    if(symmetric) then
      call dsymv('U',n,1d0,matrix,n,buffer,1,0d0,image(:,m),1)
    else
      call dgemv('N',n,n,1d0,matrix,n,buffer,1,0d0,image(:,m),1)
    end if
    if(.not.TDA) call dtrmv('L','T','N',n,factor,n,image(:,m),1)
  end subroutine append_trial

end subroutine STEOM_eB_davidson
