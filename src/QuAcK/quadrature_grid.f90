subroutine quadrature_grid(nRad,nAng,nGrid,root,weight)

! Build roots and weights of quadrature grid

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                 :: nRad,nAng,nGrid

! Local variables

  integer                            :: i,j,k
  double precision                   :: scale
  double precision,allocatable       :: Radius(:)
  double precision,allocatable       :: RadWeight(:)
  double precision,allocatable       :: XYZ(:,:)
  double precision,allocatable       :: XYZWeight(:)

! Output variables

  double precision,intent(out)       :: root(3,nGrid)
  double precision,intent(out)       :: weight(nGrid)

! Memory allocation

  allocate(Radius(nRad),RadWeight(nRad),XYZ(3,nAng),XYZWeight(nAng))

! Findthe radial grid

  scale = 1d0
  call EulMac(Radius,RadWeight,nRad,scale)

  write(*,20)
  write(*,30)
  write(*,20)
  do i = 1,nRad
     write(*,40) i,Radius(i),RadWeight(i)
  end do
  write(*,20)
  write(*,*)

! Find the angular grid

  call Lebdev(XYZ,XYZWeight,nAng)

  write(*,20)
  write(*,50)
  write(*,20)
  do j = 1,nAng
     write(*,60) j,(XYZ(k,j),k=1,3), XYZWeight(j)
  end do
  write(*,20)

! Form the roots and weights

  k = 0
  do i=1,nRad
    do j=1,nAng
       k = k + 1
       root(:,k) = Radius(i)*XYZ(:,j)
       weight(k) = RadWeight(i)*XYZWeight(j)
    enddo
  enddo
 
! Compute values of the basis functions (and the its gradient if required) at each grid point

20 format(T2,58('-'))
30 format(T20,'Radial Quadrature',/, &
          T6,'I',T26,'Radius',T50,'Weight')
40 format(T3,I4,T18,F17.10,T35,F25.10)
50 format(T20,'Angular Quadrature',/, &
          T6,'I',T19,'X',T29,'Y',T39,'Z',T54,'Weight')
60 format(T3,I4,T13,3F10.5,T50,F10.5)

end subroutine quadrature_grid
