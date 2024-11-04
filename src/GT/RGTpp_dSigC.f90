double precision function RGTpp_dSigC(p,w,eta,nBas,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,e,Om1s,rho1s,Om2s,rho2s, &
                                     Om1t,rho1t,Om2t,rho2t)

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
  integer,intent(in)            :: nOOs,nOOt
  integer,intent(in)            :: nVVs,nVVt

  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: Om1s(nVVs),Om1t(nVVt)
  double precision,intent(in)   :: rho1s(nBas,nBas,nVVs),rho1t(nBas,nBas,nVVt)
  double precision,intent(in)   :: Om2s(nOOs),Om2t(nOOt)
  double precision,intent(in)   :: rho2s(nBas,nBas,nOOs),rho2t(nBas,nBas,nOOt)

! Local variables

  integer                       :: i,a,cd,kl
  double precision              :: eps

! Initialize 

  RGTpp_dSigC = 0d0

!----------------------------------------------
! Occupied part of the T-matrix self-energy 
!----------------------------------------------

  do i=nC+1,nO

    do cd=1,nVVs
      eps = w + e(i) - Om1s(cd)
      RGTpp_dSigC = RGTpp_dSigC - (1d0/2d0)*rho1s(p,i,cd)**2*(eps**2 - eta**2)/(eps**2 + eta**2)**2
    end do

    do cd=1,nVVt
      eps = w + e(i) - Om1t(cd)
      RGTpp_dSigC = RGTpp_dSigC - (3d0/2d0)*rho1t(p,i,cd)**2*(eps**2 - eta**2)/(eps**2 + eta**2)**2
    end do

  end do

!----------------------------------------------
! Virtual part of the T-matrix self-energy
!----------------------------------------------

  do a=nO+1,nBas-nR

    do kl=1,nOOs
      eps = w + e(a) - Om2s(kl)
      RGTpp_dSigC = RGTpp_dSigC - (1d0/2d0)*rho2s(p,a,kl)**2*(eps**2 - eta**2)/(eps**2 + eta**2)**2
    end do

    do kl=1,nOOt
      eps = w + e(a) - Om2t(kl)
      RGTpp_dSigC = RGTpp_dSigC - (3d0/2d0)*rho2t(p,a,kl)**2*(eps**2 - eta**2)/(eps**2 + eta**2)**2
    end do

  end do
     
end function 
