subroutine RGTpp_self_energy_iomega(eta,wcoord,nBas,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,e,Om1s,rho1s,Om2s,rho2s,Om1t,rho1t,Om2t,rho2t,EcGM,Sig)

! Compute the correlation part of the T-matrix self-energy and the renormalization factor

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
  complex*16,intent(in)         :: wcoord

! Local variables

  integer                       :: i,j,a,b,p,q,cd,kl
  complex*16                    :: num,eps

! Output variables

  double precision,intent(inout):: EcGM
  complex*16,intent(inout)      :: Sig(nBas,nBas)

! Initialization

  Sig(:,:) = czero
  EcGM     = 0d0

!----------------------------------------------
! Occupied part of the T-matrix self-energy 
!----------------------------------------------

  do p=nC+1,nBas-nR
    do q=nC+1,nBas-nR
      do i=nC+1,nO

        do cd=1,nVVs
          eps = wcoord + e(i) - Om1s(cd)
          num = (1d0/2d0)*rho1s(p,i,cd)*rho1s(q,i,cd)
          Sig(p,q) = Sig(p,q) + num/eps
        end do

        do cd=1,nVVt
          eps = wcoord + e(i) - Om1t(cd)
          num = (3d0/2d0)*rho1t(p,i,cd)*rho1t(q,i,cd)
          Sig(p,q) = Sig(p,q) + num/eps
        end do

      end do
    end do
  end do

!----------------------------------------------
! Virtual part of the T-matrix self-energy
!----------------------------------------------

  do p=nC+1,nBas-nR
    do q=nC+1,nBas-nR
      do a=nO+1,nBas-nR

        do kl=1,nOOs
          eps = wcoord + e(a) - Om2s(kl)
          num = (1d0/2d0)*rho2s(p,a,kl)*rho2s(q,a,kl)
          Sig(p,q) = Sig(p,q) + num/eps
        end do

        do kl=1,nOOt
          eps = wcoord + e(a) - Om2t(kl)
          num = (3d0/2d0)*rho2t(p,a,kl)*rho2t(q,a,kl)
          Sig(p,q) = Sig(p,q) + num/eps
        end do

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
        num = (1d0/2d0)*rho1s(i,j,cd)*rho1s(i,j,cd)
        EcGM = EcGM + real(num*eps/(eps**2 + eta**2))
      end do

      do cd=1,nVVt
        eps = e(i) + e(j) - Om1t(cd)
        num = (3d0/2d0)*rho1t(i,j,cd)*rho1t(i,j,cd)
        EcGM = EcGM + real(num*eps/(eps**2 + eta**2))
      end do

    end do
  end do

  do a=nO+1,nBas-nR
    do b=nO+1,nBas-nR

      do kl=1,nOOs
        eps = e(a) + e(b) - Om2s(kl)
        num = (1d0/2d0)*rho2s(a,b,kl)*rho2s(a,b,kl)
        EcGM = EcGM - real(num*eps/(eps**2 + eta**2))
      end do

      do kl=1,nOOt
        eps = e(a) + e(b) - Om2t(kl)
        num = (3d0/2d0)*rho2t(a,b,kl)*rho2t(a,b,kl)
        EcGM = EcGM - real(num*eps/(eps**2 + eta**2))
      end do

    end do
  end do


end subroutine 
