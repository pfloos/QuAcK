subroutine unrestricted_renormalization_factor(eta,nBas,nC,nO,nV,nR,nSa,nSb,nSt,e,Omega,rho,Z)

! Compute the renormalization factor in the unrestricted formalism

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nSa
  integer,intent(in)            :: nSb
  integer,intent(in)            :: nSt
  double precision,intent(in)   :: e(nBas,nspin)
  double precision,intent(in)   :: Omega(nSt)
  double precision,intent(in)   :: rho(nBas,nBas,nSt,nspin)

! Local variables

  integer                       :: i,j,a,b,p,q,jb
  double precision              :: eps

! Output variables

  double precision,intent(out)  :: Z(nBas,nspin)

! Initialize 

  Z(:,:) = 0d0

!--------------!
! Spin-up part !
!--------------!

  ! Occupied part of the renormalization factor

  do p=nC(1)+1,nBas-nR(1)
    do i=nC(1)+1,nO(1)
      do jb=1,nSt
        eps = e(p,1) - e(i,1) + Omega(jb)
        Z(p,1) = Z(p,1) + rho(p,i,jb,1)**2*(eps/(eps**2 + eta**2))**2
      end do
    end do
  end do

  ! Virtual part of the correlation self-energy

  do p=nC(1)+1,nBas-nR(1)
    do a=nO(1)+1,nBas-nR(1)
      do jb=1,nSt
        eps = e(p,1) - e(a,1) - Omega(jb)
        Z(p,1) = Z(p,1) + rho(p,a,jb,1)**2*(eps/(eps**2 + eta**2))**2
      end do
    end do
  end do

!----------------!
! Spin-down part !
!----------------!

  ! Occupied part of the correlation self-energy

  do p=nC(2)+1,nBas-nR(2)
    do i=nC(2)+1,nO(2)
      do jb=1,nSt
        eps = e(p,2) - e(i,2) + Omega(jb)
        Z(p,2) = Z(p,2) + rho(p,i,jb,2)**2*(eps/(eps**2 + eta**2))**2
      end do
    end do
  end do

  ! Virtual part of the correlation self-energy

  do p=nC(2)+1,nBas-nR(2)
    do a=nO(2)+1,nBas-nR(2)
      do jb=1,nSt
        eps = e(p,2) - e(a,2) - Omega(jb)
        Z(p,2) = Z(p,2) + rho(p,a,jb,2)**2*(eps/(eps**2 + eta**2))**2
      end do
    end do
  end do

! Final rescaling

  Z(:,:) = 1d0/(1d0 + Z(:,:))


end subroutine unrestricted_renormalization_factor
