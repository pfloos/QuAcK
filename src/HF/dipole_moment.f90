subroutine dipole_moment(nBas,P,nNuc,ZNuc,rNuc,dipole_int,dipole)

! Compute density matrix based on the occupation numbers

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: P(nBas,nBas)
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  integer                       :: ixyz
  integer                       :: iNuc
  integer                       :: mu,nu

! Output variables

  double precision,intent(out)  :: dipole(ncart)

! Initialization

  dipole(:) = 0d0

! Loop over cartesian components

  do ixyz=1,ncart

    ! Nuclear part

    do iNuc=1,nNuc

      dipole(ixyz) = dipole(ixyz) + ZNuc(iNuc)*rNuc(iNuc,ixyz)

    end do

    ! Electronic part

    do nu=1,nBas
      do mu=1,nBas
        dipole(ixyz) = dipole(ixyz) - P(mu,nu)*dipole_int(mu,nu,ixyz)
      end do
    end do

  end do

end subroutine 
