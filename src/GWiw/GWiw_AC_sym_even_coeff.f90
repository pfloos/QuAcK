subroutine GWiw_AC_sym_even_coeff(nZ,Zw,F,nblkmax,nblk,coef)

! Symmetric-even analytic continuation of an even, real-analytic function f(w) = f(z), z = w^2,
! from reference values F(i) = f(Zw(i)) given at arbitrary complex frequencies Zw(i).
!
! The function is represented by the continued fraction
!
!   f = a1 / ( b1 + c1 (z-z1) + ( (z-z1)^2 + d1 ) a2 / ( b2 + c2 (z-z2) + ( (z-z2)^2 + d2 ) ... ) )
!
! with real coefficients, so that f(z*) = f(z)* holds exactly. Each block interpolates either a pair of
! real-valued reference points (w real, or w purely imaginary) or one general complex point together with
! its complex conjugate. Blocks are chosen greedily among the stamps with the largest residual.
!
! The algorithm follows I. Duchemin (beDeft, GPLv3, AnalyticContinuation.h,
! SymmetricEvenAnalyticContinuation_coeff) and was checked against the Python port in TITOU/SRG.
!
! coef(1:5,k) = a_k, b_k, c_k, d_k, z_k

  implicit none

! Input variables

  integer,intent(in)            :: nZ
  integer,intent(in)            :: nblkmax
  complex(kind=8),intent(in)    :: Zw(nZ)
  complex(kind=8),intent(in)    :: F(nZ)

! Output variables

  integer,intent(out)           :: nblk
  double precision,intent(out)  :: coef(5,nblkmax)

! Local variables

  integer,parameter             :: RSTACK = 8
  integer,parameter             :: ISTACK = 8
  integer,parameter             :: CSTACK = 32
  double precision,parameter    :: eps = 1d-16
  double precision,parameter    :: rel = 1d-12

  integer                       :: i,j,k,p,m,ncl,nrem(3),ic
  integer                       :: cls(nZ)
  integer                       :: top(CSTACK,3),ntop(3),stk(3)
  logical                       :: used(nZ)
  double precision              :: res(nZ),resmax(3),thresh,best,rtot,fk
  double precision              :: a,b,c,d,z
  double precision              :: abest,bbest,cbest,dbest,zbest
  integer                       :: kbest,pbest,cbest_cls
  complex(kind=8)               :: zz(nZ),FF(nZ),GG(nZ),A0(nZ),A1(nZ),B0(nZ),B1(nZ),Sf(nZ)
  complex(kind=8)               :: an,bn,A2,B2,Sn,r
  double precision              :: nrm
  logical                       :: found

  stk(1) = RSTACK
  stk(2) = ISTACK
  stk(3) = CSTACK

! Classify stamps: 1 = real axis, 2 = imaginary axis, 3 = general complex

  do i=1,nZ
    if(abs(aimag(Zw(i))) < eps) then
      cls(i) = 1
    else if(abs(real(Zw(i))) < eps) then
      cls(i) = 2
    else
      cls(i) = 3
    end if
    zz(i) = Zw(i)*Zw(i)
    if(cls(i) == 3) then
      FF(i) = F(i)
    else
      FF(i) = cmplx(real(F(i)),0d0,kind=8)
      zz(i) = cmplx(real(zz(i)),0d0,kind=8)
    end if
  end do

  GG   = FF
  A0   = (1d0,0d0)
  A1   = (0d0,0d0)
  B0   = (0d0,0d0)
  B1   = (1d0,0d0)
  Sf   = (1d0,0d0)
  used = .false.

! Convergence threshold from the initial residuals

  resmax = 0d0
  do i=1,nZ
    resmax(cls(i)) = max(resmax(cls(i)),abs(FF(i))**2)
  end do
  thresh = maxval(resmax)*rel

  nblk = 0

  do while(nblk < nblkmax)

!   Residuals of the current convergent for the stamps not yet used

    nrem   = 0
    resmax = 0d0
    do i=1,nZ
      if(used(i)) then
        res(i) = -1d0
      else
        nrem(cls(i)) = nrem(cls(i)) + 1
        res(i) = abs(A1(i)/B1(i)-FF(i))**2
        if(.not.(res(i) < huge(1d0))) res(i) = -1d0
        resmax(cls(i)) = max(resmax(cls(i)),res(i))
      end if
    end do

    if(.not.( (resmax(1) > thresh .and. nrem(1) >= 2) .or. &
              (resmax(2) > thresh .and. nrem(2) >= 2) .or. &
              (resmax(3) > thresh .and. nrem(3) >= 1) )) exit

!   Collect the stamps with the largest residual in each class

    ntop = 0
    do ic=1,3
      do m=1,stk(ic)
        found = .false.
        k = 0
        rtot = 0d0
        do i=1,nZ
          if(cls(i) /= ic .or. used(i) .or. res(i) <= 0d0) cycle
          if(any(top(1:ntop(ic),ic) == i)) cycle
          if(res(i) > rtot) then
            rtot  = res(i)
            k     = i
            found = .true.
          end if
        end do
        if(.not.found) exit
        ntop(ic) = ntop(ic) + 1
        top(ntop(ic),ic) = k
      end do
    end do

!   Try all candidate blocks

    best = 1d301
    kbest = 0
    pbest = 0
    cbest_cls = 0
    abest = 0d0; bbest = 0d0; cbest = 0d0; dbest = 0d0; zbest = 0d0

    do ic=1,2
      if(nrem(ic) < 2) cycle
      do i=1,ntop(ic)
        do j=i+1,ntop(ic)
          k = top(i,ic)
          p = top(j,ic)
          a = real(GG(k))*real(GG(p))
          b = 0.5d0*(real(GG(k))+real(GG(p)))
          c = -(real(GG(k))-real(GG(p)))/(real(zz(k))-real(zz(p)))
          d = -(real(zz(k))-real(zz(p)))**2/4d0
          z = 0.5d0*(real(zz(k))+real(zz(p)))
          if(abs(a) <= 0d0) cycle
          fk = 1d0/sqrt(abs(a))
          a = a*fk
          b = b*fk
          c = c*fk
          rtot = 0d0
          do m=1,nZ
            if(used(m) .or. m == k .or. m == p) cycle
            an = a*Sf(m)
            bn = b + c*(zz(m)-z)
            A2 = bn*A1(m) + an*A0(m)
            B2 = bn*B1(m) + an*B0(m)
            r  = A2/B2 - FF(m)
            if(abs(r)**2 < huge(1d0)) then
              rtot = rtot + abs(r)**2
            else
              rtot = rtot + 1d300
            end if
          end do
          if(rtot < best) then
            best = rtot
            kbest = k
            pbest = p
            cbest_cls = ic
            abest = a; bbest = b; cbest = c; dbest = d; zbest = z
          end if
        end do
      end do
    end do

    if(nrem(3) >= 1) then
      do i=1,ntop(3)
        k = top(i,3)
        a = abs(GG(k))**2
        if(a <= 0d0) cycle
        b = real(GG(k))
        c = -aimag(GG(k))/aimag(zz(k))
        d = aimag(zz(k))**2
        z = real(zz(k))
        fk = 1d0/sqrt(abs(a))
        a = a*fk
        b = b*fk
        c = c*fk
        rtot = 0d0
        do m=1,nZ
          if(used(m) .or. m == k) cycle
          an = a*Sf(m)
          bn = b + c*(zz(m)-z)
          A2 = bn*A1(m) + an*A0(m)
          B2 = bn*B1(m) + an*B0(m)
          r  = A2/B2 - FF(m)
          if(abs(r)**2 < huge(1d0)) then
            rtot = rtot + abs(r)**2
          else
            rtot = rtot + 1d300
          end if
        end do
        if(rtot < best) then
          best = rtot
          kbest = k
          pbest = 0
          cbest_cls = 3
          abest = a; bbest = b; cbest = c; dbest = d; zbest = z
        end if
      end do
    end if

    if(kbest == 0) exit

!   Accept the best block

    nblk = nblk + 1
    coef(1,nblk) = abest
    coef(2,nblk) = bbest
    coef(3,nblk) = cbest
    coef(4,nblk) = dbest
    coef(5,nblk) = zbest

    used(kbest) = .true.
    if(pbest > 0) used(pbest) = .true.

!   Update the working arrays of all remaining stamps

    do m=1,nZ
      if(used(m)) cycle
      an = abest*Sf(m)
      bn = bbest + cbest*(zz(m)-zbest)
      Sn = (zz(m)-zbest)**2 + dbest
      GG(m) = (abest/GG(m)-bn)/Sn
      A2 = bn*A1(m) + an*A0(m)
      B2 = bn*B1(m) + an*B0(m)
      Sf(m) = Sn
      A0(m) = A1(m)
      A1(m) = A2
      B0(m) = B1(m)
      B1(m) = B2
      nrm = 1d0/max(abs(A1(m)),abs(B1(m)))
      A0(m) = A0(m)*nrm
      A1(m) = A1(m)*nrm
      B0(m) = B0(m)*nrm
      B1(m) = B1(m)*nrm
    end do

  end do

end subroutine
