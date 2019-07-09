subroutine fc_srlda(nBas,nGrid,weight,MO,rho,mu,eG0W0)

! Compute sr-lda ec

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: MO(nBas,nGrid)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: mu(nGrid)
  double precision,intent(in)   :: eG0W0(nBas)

! Local variables

  integer                       :: iG,p
  double precision              :: r
  double precision              :: rs
  double precision              :: ec,vcup,vcdw
  double precision,allocatable  :: de(:)

! Memory allocation

  allocate(de(nBas))

! Initialization

  de(:) = 0d0

  do iG=1,ngrid

    r = max(0d0,rho(iG))

    if(r > threshold) then

      rs = (4d0*pi*r/3d0)**(-1d0/3d0)

      call lsdsr(rs,0d0,mu(iG),ec,vcup,vcdw)

      do p=1,nBas

        de(p)= de(p) + weight(iG)*vcup*MO(p,iG)**2

      end do
 
    end if

  end do

  print*, 'Eigenvalues correction from srDFT (in eV)'
  call matout(nBas,1,de(:)*HaToeV)

  print*, 'Corrected G0W0 eigenvalues (in eV)'
  call matout(nBas,1,(eG0W0(:) + de(:))*HaToeV)

  

end subroutine fc_srlda
