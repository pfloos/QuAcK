subroutine spatial_to_spin_MO_energy(nBas,e,nBas2,se)

! Convert MO energies from spatial to spin orbitals

  implicit none

! Input variables

  integer,intent(in)            :: nBas,nBas2
  double precision,intent(in)   :: e(nBas)

! Local variables

  integer                       :: p

! Output variables

  double precision,intent(out)  :: se(nBas2)

    do p=1,nBas2

      se(p) = e((p+1)/2)

    enddo

! print*,'MO energies in spinorbital basis'
! call matout(nBas2,1,se)

end subroutine spatial_to_spin_MO_energy
