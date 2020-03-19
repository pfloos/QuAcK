subroutine G0F2(linearize,nBas,nC,nO,nV,nR,V,e0)

! Perform a one-shot second-order Green function calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: linearize
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nO
  integer,intent(in)            :: nC
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: e0(nBas)
  double precision,intent(in)   :: V(nBas,nBas,nBas,nBas)

! Local variables

  double precision              :: eps
  double precision              :: VV
  double precision,allocatable  :: eGF2(:)
  double precision,allocatable  :: Bpp(:,:)
  double precision,allocatable  :: Z(:)

  integer                       :: i,j,a,b,p

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|     One-shot second-order Green function     |'
  write(*,*)'************************************************'
  write(*,*)

! Memory allocation

  allocate(Bpp(nBas,2),Z(nBas),eGF2(nBas))


  ! Frequency-dependent second-order contribution

  Bpp(:,:) = 0d0
  Z(:)     = 0d0

  do p=nC+1,nBas-nR
    do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nBas-nR

          eps = e0(p) + e0(a) - e0(i) - e0(j)
          VV  = (2d0*V(p,a,i,j) - V(p,a,j,i))*V(p,a,i,j)
          Bpp(p,1) = Bpp(p,1) + VV/eps
          Z(p)     = Z(p)     + VV/eps**2

        end do
      end do
    end do
  end do

  do p=nC+1,nBas-nR
    do i=nC+1,nO
      do a=nO+1,nBas-nR
        do b=nO+1,nBas-nR

          eps = e0(p) + e0(i) - e0(a) - e0(b)
          VV  = (2d0*V(p,i,a,b) - V(p,i,b,a))*V(p,i,a,b)
          Bpp(p,2) = Bpp(p,2) + VV/eps
          Z(p)     = Z(p)     + VV/eps**2

        end do
      end do
    end do
  end do

  Z(:) = 1d0/(1d0 + Z(:))

  if(linearize) then

    eGF2(:) = e0(:) + Z(:)*(Bpp(:,1) + Bpp(:,2))

  else

    eGF2(:) = e0(:) + Bpp(:,1) + Bpp(:,2)

  end if
  ! Print results

  call print_G0F2(nBas,nO,e0,eGF2,Z)

end subroutine G0F2
