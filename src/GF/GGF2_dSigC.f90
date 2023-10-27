double precision function GGF2_dSigC(p,w,eta,nBas,nC,nO,nV,nR,eHF,ERI)

! Compute diagonal of the correlation part of the self-energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: p
  double precision,intent(in)   :: w
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas,nC,nO,nV,nR
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: i,j,a,b
  double precision              :: eps

! Initialize 

  GGF2_dSigC = 0d0

  ! Occupied part of the correlation self-energy

  do i=nC+1,nO
    do j=nC+1,nO
      do a=nO+1,nBas-nR

        eps = w + eHF(a) - eHF(i) - eHF(j)
        GGF2_dSigC = GGF2_dSigC - 0.5d0*(ERI(p,a,i,j) - ERI(p,a,j,i))**2*(eps**2 - eta**2)/(eps**2 + eta**2)**2

      end do
    end do
  end do

  ! Virtual part of the correlation self-energy

  do i=nC+1,nO
    do a=nO+1,nBas-nR
      do b=nO+1,nBas-nR

        eps = w + eHF(i) - eHF(a) - eHF(b)
        GGF2_dSigC = GGF2_dSigC - 0.5d0*(ERI(p,i,a,b) - ERI(p,i,b,a))**2*(eps**2 - eta**2)/(eps**2 + eta**2)**2

      end do
    end do
  end do

end function 
