subroutine GGTpp_self_energy_diag(eta,nBas,nC,nO,nV,nR,nOO,nVV,e,Om1,rho1,Om2,rho2,EcGM,Sig,Z)

! Compute diagonal of the correlation part of the T-matrix self-energy

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: Om1(nVV)
  double precision,intent(in)   :: rho1(nBas,nBas,nVV)
  double precision,intent(in)   :: Om2(nOO)
  double precision,intent(in)   :: rho2(nBas,nBas,nOO)

! Local variables

  integer                       :: i,j,a,b,p,cd,kl
  double precision              :: num,eps

! Output variables

  double precision,intent(inout)  :: EcGM
  double precision,intent(inout)  :: Sig(nBas)
  double precision,intent(inout)  :: Z(nBas)

! Initialization

  Sig(:) = 0d0
  Z(:)   = 0d0
  EcGM   = 0d0

!--------------------------------------!
! Occupied part of the Tpp self-energy !
!--------------------------------------!

  do p=nC+1,nBas-nR
    do i=nC+1,nO

      do cd=1,nVV
        eps = e(p) + e(i) - Om1(cd)
        num = rho1(p,i,cd)**2
        Sig(p) = Sig(p) + num*eps/(eps**2 + eta**2)
        Z(p)   = Z(p)   - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
      end do

    end do
  end do

!------------------------------------------!
! Virtual part of the T-matrix self-energy !
!------------------------------------------!

  do p=nC+1,nBas-nR
    do a=nO+1,nBas-nR

      do kl=1,nOO
        eps = e(p) + e(a) - Om2(kl)
        num = rho2(p,a,kl)**2
        Sig(p) = Sig(p) + num*eps/(eps**2 + eta**2)
        Z(p)   = Z(p)   - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
      end do

    end do
  end do

!-------------------------------------!
! Galitskii-Migdal correlation energy !
!-------------------------------------!

  do i=nC+1,nO
    do j=nC+1,nO

      do cd=1,nVV
        eps = e(i) + e(j) - Om1(cd)
        num = rho1(i,j,cd)**2
        EcGM = EcGM + num*eps/(eps**2 + eta**2)
      end do

    end do
  end do

  do a=nO+1,nBas-nR
    do b=nO+1,nBas-nR

      do kl=1,nOO
        eps = e(a) + e(b) - Om2(kl)
        num = rho2(a,b,kl)**2
        EcGM = EcGM - num*eps/(eps**2 + eta**2)
      end do

    end do
  end do

  Z(:) = 1d0/(1d0 - Z(:))

end subroutine 
