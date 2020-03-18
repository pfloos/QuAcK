subroutine hartree_coulomb(nBas,P,ERI,J)

! Compute Coulomb matrix

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: P(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: mu,nu,la,si

! Output variables

  double precision,intent(out)  :: J(nBas,nBas)

  J = 0d0
  do si=1,nBas
    do la=1,nBas
      do nu=1,nBas
        do mu=1,nBas
          J(mu,nu) = J(mu,nu) + P(la,si)*ERI(mu,la,nu,si)
        enddo
      enddo
    enddo
  enddo


end subroutine hartree_coulomb
