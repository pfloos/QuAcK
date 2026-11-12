double precision function RGF2_SRG_SigC(p,w,flow,nBas,nC,nO,nV,nR,eHF,ERI)

! Compute diagonal of the correlation part of the self-energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: p
  double precision,intent(in)   :: w
  double precision,intent(in)   :: flow
  integer,intent(in)            :: nBas,nC,nO,nV,nR
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: i,j,a,b
  double precision              :: eps,s,ratio
  double precision,external     :: RGF2_SRG_ratio

  RGF2_SRG_SigC = 0d0
  s = flow

  ! Occupied part of the correlation self-energy

  do i=nC+1,nO
    do j=nC+1,nO
      do a=nO+1,nBas-nR

        eps = w + eHF(a) - eHF(i) - eHF(j)
        ratio = RGF2_SRG_ratio(2d0*s,eps**2)
        RGF2_SRG_SigC = RGF2_SRG_SigC + &
          (eps*ratio)*(2d0*ERI(p,a,i,j) - ERI(p,a,j,i))*ERI(p,a,i,j)

        end do
     end do
  end do

  ! Virtual part of the correlation self-energy

  do i=nC+1,nO
    do a=nO+1,nBas-nR
      do b=nO+1,nBas-nR

        eps = w + eHF(i) - eHF(a) - eHF(b)
        ratio = RGF2_SRG_ratio(2d0*s,eps**2)
        RGF2_SRG_SigC = RGF2_SRG_SigC + &
          (eps*ratio)*(2d0*ERI(p,i,a,b) - ERI(p,i,b,a))*ERI(p,i,a,b)

      end do
    end do
  end do

end function 


double precision function RGF2_SRG_ratio(scale,delta_sq)

! Evaluate (1 - exp(-scale*delta_sq))/delta_sq, including its limit at zero.

  implicit none
  double precision,intent(in) :: scale,delta_sq
  double precision            :: x

  x = scale*delta_sq
  if(x < 1d-5) then
    RGF2_SRG_ratio = scale*(1d0 - x*(0.5d0 - x*(1d0/6d0 - x/24d0)))
  else
    RGF2_SRG_ratio = (1d0 - exp(-x))/delta_sq
  end if

end function
