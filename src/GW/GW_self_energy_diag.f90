subroutine GW_self_energy_diag(eta,nBas,nC,nO,nV,nR,nS,e,Omega,rho,EcGM,Sig,Z)

! Compute diagonal of the correlation part of the self-energy and the renormalization factor

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: Omega(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)

! Local variables

  integer                       :: i,a,p,q,jb
  double precision              :: num,eps

! Output variables

  double precision,intent(out)  :: Sig(nBas)
  double precision,intent(out)  :: Z(nBas)
  double precision,intent(out)  :: EcGM

! Initialize 

  Sig(:) = 0d0
  Z(:)   = 0d0

!----------------!
! GW self-energy !
!----------------!
 
! Occupied part of the correlation self-energy

  do p=nC+1,nBas-nR
    do i=nC+1,nO
      do jb=1,nS

        eps = e(p) - e(i) + Omega(jb)
        num = 2d0*rho(p,i,jb)**2
        Sig(p) = Sig(p) + num*eps/(eps**2 + eta**2)
        Z(p)   = Z(p)   - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

      end do
    end do
  end do

! Virtual part of the correlation self-energy

  do p=nC+1,nBas-nR
    do a=nO+1,nBas-nR
      do jb=1,nS

        eps = e(p) - e(a) - Omega(jb)
        num = 2d0*rho(p,a,jb)**2
        Sig(p) = Sig(p) + num*eps/(eps**2 + eta**2)
        Z(p)   = Z(p)   - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

      end do
    end do
  end do

! Galitskii-Migdal correlation energy

  EcGM = 0d0
  do i=nC+1,nO
    do a=nO+1,nBas-nR
      do jb=1,nS

        eps = e(a) - e(i) + Omega(jb)
        num = 4d0*rho(a,i,jb)**2
        EcGM = EcGM - num*eps/(eps**2 + eta**2)

      end do
    end do
  end do

! Compute renormalization factor from derivative 

  Z(:) = 1d0/(1d0 - Z(:))

end subroutine 
