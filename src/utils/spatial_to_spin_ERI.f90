subroutine spatial_to_spin_ERI(nBas,ERI,nBas2,sERI)

! Convert ERIs from spatial to spin orbitals

  implicit none

! Input variables

  integer,intent(in)            :: nBas,nBas2
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: p,q,r,s
  double precision,external     :: Kronecker_delta

! Output variables

  double precision,intent(out)  :: sERI(nBas2,nBas2,nBas2,nBas2)

    do p=1,nBas2
      do q=1,nBas2
        do r=1,nBas2
          do s=1,nBas2

            sERI(p,q,r,s) = Kronecker_delta(mod(p,2),mod(r,2)) &
                          * Kronecker_delta(mod(q,2),mod(s,2)) &
                          * ERI((p+1)/2,(q+1)/2,(r+1)/2,(s+1)/2)

          enddo
        enddo
      enddo
    enddo

end subroutine spatial_to_spin_ERI
