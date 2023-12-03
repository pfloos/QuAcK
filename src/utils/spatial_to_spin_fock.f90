subroutine spatial_to_spin_fock(nBas,F,nBas2,sF)

! Convert Fock matrix from spatial to spin orbitals

  implicit none

! Input variables

  integer,intent(in)            :: nBas,nBas2
  double precision,intent(in)   :: F(nBas,nBas)

! Local variables

  integer                       :: p,q
  double precision,external     :: Kronecker_delta

! Output variables

  double precision,intent(out)  :: sF(nBas2,nBas2)

    do p=1,nBas2
      do q=1,nBas2

        sF(p,q) = Kronecker_delta(mod(p,2),mod(q,2))*F((p+1)/2,(q+1)/2)

      end do
    end do

end subroutine spatial_to_spin_fock
