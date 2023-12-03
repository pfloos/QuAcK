subroutine GTpp_self_energy_diag(eta,nBas,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,e,Om1s,rho1s,Om2s,rho2s,Om1t,rho1t,Om2t,rho2t, &
                                 EcGM,Sig,Z)

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
  integer,intent(in)            :: nOOs,nOOt
  integer,intent(in)            :: nVVs,nVVt
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: Om1s(nVVs),Om1t(nVVt)
  double precision,intent(in)   :: rho1s(nBas,nBas,nVVs),rho1t(nBas,nBas,nVVt)
  double precision,intent(in)   :: Om2s(nOOs),Om2t(nOOt)
  double precision,intent(in)   :: rho2s(nBas,nBas,nOOs),rho2t(nBas,nBas,nOOt)

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

      do cd=1,nVVs
        eps = e(p) + e(i) - Om1s(cd)
        num = rho1s(p,i,cd)**2
        Sig(p) = Sig(p) + num*eps/(eps**2 + eta**2)
        Z(p)   = Z(p)   - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
      end do

      do cd=1,nVVt
        eps = e(p) + e(i) - Om1t(cd)
        num = rho1t(p,i,cd)**2
        Sig(p) = Sig(p) + num*eps/(eps**2 + eta**2)
        Z(p)   = Z(p)   - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
      end do

    end do
  end do

!----------------------------------------------
! Virtual part of the T-matrix self-energy
!----------------------------------------------

  do p=nC+1,nBas-nR
    do a=nO+1,nBas-nR

      do kl=1,nOOs
        eps = e(p) + e(a) - Om2s(kl)
        num = rho2s(p,a,kl)**2
        Sig(p) = Sig(p) + num*eps/(eps**2 + eta**2)
        Z(p)   = Z(p)   - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
      end do

      do kl=1,nOOt
        eps = e(p) + e(a) - Om2t(kl)
        num = rho2t(p,a,kl)**2
        Sig(p) = Sig(p) + num*eps/(eps**2 + eta**2)
        Z(p)   = Z(p)   - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
      end do

    end do
  end do

!----------------------------------------------
! Galitskii-Migdal correlation energy
!----------------------------------------------

  do i=nC+1,nO
    do j=nC+1,nO

      do cd=1,nVVs
        eps = e(i) + e(j) - Om1s(cd)
        num = rho1s(i,j,cd)**2
        EcGM = EcGM + num*eps/(eps**2 + eta**2)
      end do

      do cd=1,nVVt
        eps = e(i) + e(j) - Om1t(cd)
        num = rho1t(i,j,cd)**2
        EcGM = EcGM + num*eps/(eps**2 + eta**2)
      end do

    end do
  end do

  do a=nO+1,nBas-nR
    do b=nO+1,nBas-nR

      do kl=1,nOOs
        eps = e(a) + e(b) - Om2s(kl)
        num = rho2s(a,b,kl)**2
        EcGM = EcGM - num*eps/(eps**2 + eta**2)
      end do

      do kl=1,nOOt
        eps = e(a) + e(b) - Om2t(kl)
        num = rho2t(a,b,kl)**2
        EcGM = EcGM - num*eps/(eps**2 + eta**2)
      end do

    end do
  end do

  Z(:) = 1d0/(1d0 - Z(:))

end subroutine 
