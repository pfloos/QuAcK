double precision function cRGF2_Re_SigC(p,Re_w,Im_w,eta,nBas,nC,nO,nV,nR,eHF,e_cap,ERI)

! Compute diagonal of the correlation part of the self-energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: p
  double precision,intent(in)   :: Re_w
  double precision,intent(in)   :: Im_w
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas,nC,nO,nV,nR
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: e_cap(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: i,j,a,b
  double precision              :: eps
  double precision              :: eta_tilde
  double precision              :: num

! Initialize 

  cRGF2_Re_SigC = 0d0

  ! Occupied part of the correlation self-energy

  do i=nC+1,nO
    do j=nC+1,nO
      do a=nO+1,nBas-nR

        eps = Re_w + eHF(a) - eHF(i) - eHF(j)
        eta_tilde = eta - Im_w + e_cap(i) + e_cap(a) - e_cap(j)
        num = (2d0*ERI(p,a,i,j) - ERI(p,a,j,i))*ERI(p,a,i,j)
        cRGF2_Re_SigC = cRGF2_Re_SigC + num*eps/(eps**2 + eta_tilde**2)
      end do
    end do
  end do

  ! Virtual part of the correlation self-energy

  do i=nC+1,nO
    do a=nO+1,nBas-nR
      do b=nO+1,nBas-nR

        eps = Re_w + eHF(i) - eHF(a) - eHF(b)
        num = (2d0*ERI(p,i,a,b) - ERI(p,i,b,a))*ERI(p,i,a,b)
        eta_tilde = eta + Im_w  - e_cap(a)  - e_cap(b) + e_cap(i)
        cRGF2_Re_SigC = cRGF2_Re_SigC + num*eps/(eps**2 + eta_tilde**2)
      end do
    end do
  end do

end function 
