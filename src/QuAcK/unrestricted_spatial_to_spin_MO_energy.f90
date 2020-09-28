subroutine unrestricted_spatial_to_spin_MO_energy(nBas,e,nBas2,se)

! Convert MO energies from unrestricted spatial to spin orbitals

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nBas2
  double precision,intent(in)   :: e(nBas,nspin)

! Local variables

  integer                       :: p

! Output variables

  double precision,intent(out)  :: se(nBas2)

  do p=1,nBas2,2
    se(p) = e(p,1)
  enddo

  do p=2,nBas2,2
    se(p) = e(p,2)
  enddo

end subroutine unrestricted_spatial_to_spin_MO_energy
