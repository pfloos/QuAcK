double precision function RGF2_SigC(p,w,eta,nBas,nC,nO,nV,nR,eHF,ERI)

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

  RGF2_SigC = 0d0

  ! Occupied part of the correlation self-energy

  do i=nC+1,nO
    do j=nC+1,nO
      do a=nO+1,nBas-nR

        eps = w + eHF(a) - eHF(i) - eHF(j)
        RGF2_SigC = RGF2_SigC + (2d0*ERI(p,a,i,j) - ERI(p,a,j,i))*ERI(p,a,i,j)*eps/(eps**2 + eta**2)

        end do
     end do
  end do

  ! Virtual part of the correlation self-energy

  do i=nC+1,nO
    do a=nO+1,nBas-nR
      do b=nO+1,nBas-nR

        eps = w + eHF(i) - eHF(a) - eHF(b)
        RGF2_SigC = RGF2_SigC + (2d0*ERI(p,i,a,b) - ERI(p,i,b,a))*ERI(p,i,a,b)*eps/(eps**2 + eta**2)

      end do
    end do
  end do

end function 
