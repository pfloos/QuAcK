double precision function UGW_SigC(p,w,eta,nBas,nC,nO,nV,nR,nS,e,Om,rho)

! Compute diagonal of the correlation part of the self-energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: p
  double precision,intent(in)   :: w
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)

! Local variables

  integer                       :: i,a,jb
  double precision              :: num,eps

! Initialize 

  UGW_SigC = 0d0

! Occupied part of the correlation self-energy

  do i=nC+1,nO
    do jb=1,nS
      eps = w - e(i) + Om(jb)
      num = rho(p,i,jb)**2
      UGW_SigC = UGW_SigC + num*eps/(eps**2 + eta**2)  
    end do
  end do

! Virtual part of the correlation self-energy

  do a=nO+1,nBas-nR
    do jb=1,nS
      eps = w - e(a) - Om(jb)
      num = rho(p,a,jb)**2
      UGW_SigC = UGW_SigC + num*eps/(eps**2 + eta**2)
    end do
  end do

end function 