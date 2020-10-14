double precision function dfac(n)
  implicit none
  integer :: n
  double precision, external :: fact
  dfac = fact(n)
end

!-------------------
! The following functions were taken with from Quantum Package
! (https://github.com/QuantumPackage/qp2 , AGPL license)

double precision function fact(n)
  implicit none
  integer                        :: n
  double precision, save         :: memo(1:100)
  integer, save                  :: memomax = 1
  integer                        :: i
  double precision :: logfact

  if (n<=memomax) then
    if (n<2) then
      fact = 1.d0
    else
      fact = memo(n)
    endif
    return
  endif

  memo(1) = 1.d0
  do i=memomax+1,min(n,100)
    memo(i) = memo(i-1)*dble(i)
  enddo
  memomax = min(n,100)
  fact = dexp(logfact(n))
end function



double precision function logfact(n)
  implicit none
  integer                        :: n
  double precision, save         :: memo(1:100)
  integer, save                  :: memomax = 1
  integer                        :: i

  if (n<=memomax) then
    if (n<2) then
      logfact = 0.d0
    else
      logfact = memo(n)
    endif
    return
  endif

  memo(1) = 0.d0
  do i=memomax+1,min(n,100)
    memo(i) = memo(i-1)+dlog(dble(i))
  enddo
  memomax = min(n,100)
  logfact = memo(memomax)
  do i=101,n
    logfact = logfact + dlog(dble(i))
  enddo
end function
            

