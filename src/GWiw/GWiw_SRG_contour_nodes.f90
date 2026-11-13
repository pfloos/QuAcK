subroutine GWiw_SRG_contour_nodes(nOrb,e,flow,delta,Omax,nnodes_max,nnodes,xnode,wnode)

! Quadrature nodes x_j and weights w_j for the contour integral along the straight line z = x + i delta
! that evaluates the SRG static self-energy,
!
!   St_pq = -(1/pi) Im  int dx  K_s(e_p - z, e_q - z) Sigma_c,pq(z)
!
! Panels of Gauss-Legendre nodes:
!   - near the kernel centres (e_p+e_q)/2, |x-x0| < 3/sqrt(s): panels of width 0.5 delta (order 8)
!     (the Gaussian part of K_s oscillates as exp(2 i s x delta) there),
!   - elsewhere inside the pole region [e_min-Omax,e_max+Omax]: panels of width 3 delta (order 10),
!   - outside the pole region the panel width grows linearly with the distance to the region (no poles
!     of Sigma_c there) until the contour reaches Lext from the region edges.

  implicit none

! Input variables

  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nnodes_max
  double precision,intent(in)   :: e(nOrb)
  double precision,intent(in)   :: flow,delta,Omax

! Output variables

  integer,intent(out)           :: nnodes
  double precision,intent(out)  :: xnode(nnodes_max),wnode(nnodes_max)

! Local variables

  integer,parameter             :: ordN = 8, ordF = 10
  double precision,parameter    :: thN = 0.5d0, thF = 3d0, gfar = 0.3d0, Lext = 500d0
  integer                       :: p,q,npair,i,k,nz,ord
  double precision              :: Rn,elo,ehi,x,w,xa,xb,dist,xs,xe
  double precision              :: tN(ordN),wN(ordN),tF(ordF),wF(ordF)
  double precision,allocatable  :: zlo(:),zhi(:),tmp(:)
  logical                       :: innear

  Rn = 3d0/sqrt(flow)

! Gauss-Legendre rules on (-1,1)

  call cdgqf(ordN,1,0d0,0d0,tN,wN)
  call cdgqf(ordF,1,0d0,0d0,tF,wF)

! Near zones: merged intervals around all kernel centres

  npair = nOrb*(nOrb+1)/2
  allocate(zlo(npair),zhi(npair),tmp(npair))
  k = 0
  do p=1,nOrb
    do q=p,nOrb
      k = k + 1
      zlo(k) = 0.5d0*(e(p)+e(q)) - Rn
      zhi(k) = 0.5d0*(e(p)+e(q)) + Rn
    end do
  end do
! sort by lower edge (insertion sort)
  do i=2,npair
    xa = zlo(i); xb = zhi(i); k = i-1
    do while(k >= 1)
      if(zlo(k) <= xa) exit
      zlo(k+1) = zlo(k); zhi(k+1) = zhi(k); k = k-1
    end do
    zlo(k+1) = xa; zhi(k+1) = xb
  end do
! merge
  nz = 1
  do i=2,npair
    if(zlo(i) <= zhi(nz)) then
      zhi(nz) = max(zhi(nz),zhi(i))
    else
      nz = nz + 1
      zlo(nz) = zlo(i); zhi(nz) = zhi(i)
    end if
  end do

  elo = minval(e) - Omax
  ehi = maxval(e) + Omax
  xs  = elo - Lext
  xe  = ehi + Lext

! Build the panels

  nnodes = 0
  x = xs
  do while(x < xe - 1d-12)

    ! inside a near zone?
    innear = .false.
    xb = x
    do i=1,nz
      if(x >= zlo(i)-1d-12 .and. x < zhi(i)-1d-12) then
        innear = .true.
        xb = zhi(i)
        exit
      end if
    end do

    if(innear) then
      w   = thN*delta
      ord = ordN
    else
      ! distance to the pole region
      dist = 0d0
      if(x < elo) dist = elo - x
      if(x > ehi) dist = x - ehi
      w   = max(thF*delta,gfar*dist)
      ord = ordF
      ! do not step over the start of a near zone
      do i=1,nz
        if(zlo(i) > x+1d-12 .and. zlo(i) < x+w) w = zlo(i) - x
      end do
    end if

    xa = x
    xb = min(x+w,xe)
    if(innear) xb = min(xb,zhi(i))
    if(nnodes + ord > nnodes_max) then
      write(*,*) ' Error: too many contour nodes, increase nnodes_max'
      error stop
    end if
    if(ord == ordN) then
      do k=1,ord
        nnodes = nnodes + 1
        xnode(nnodes) = 0.5d0*(xa+xb) + 0.5d0*(xb-xa)*tN(k)
        wnode(nnodes) = 0.5d0*(xb-xa)*wN(k)
      end do
    else
      do k=1,ord
        nnodes = nnodes + 1
        xnode(nnodes) = 0.5d0*(xa+xb) + 0.5d0*(xb-xa)*tF(k)
        wnode(nnodes) = 0.5d0*(xb-xa)*wF(k)
      end do
    end if
    x = xb

  end do

  deallocate(zlo,zhi,tmp)

end subroutine
