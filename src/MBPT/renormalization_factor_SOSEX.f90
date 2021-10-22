subroutine renormalization_factor_SOSEX(eta,nBas,nC,nO,nV,nR,nS,e,Omega,rho,Z)

! Compute renormalization factor for the SOSEX version of GW

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
  double precision,intent(in)   :: Omega(nS,nspin)
  double precision,intent(in)   :: rho(nBas,nBas,nS,nspin)

! Local variables

  integer                       :: ispin
  integer                       :: i,j,a,b,p,jb
  double precision              :: eps

! Output variables

  double precision,intent(out)  :: Z(nBas)

! Initialize

  Z(:)  = 0d0

  ! Occupied part of the correlation self-energy

  do p=nC+1,nBas-nR
    do i=nC+1,nO
      jb = 0
      do j=nC+1,nO
        do b=nO+1,nBas-nR
          jb = jb + 1
          do ispin=1,nspin
            eps = e(p) - e(i) + Omega(jb,ispin) 
            Z(p) = Z(p)  - 2d0*rho(p,i,jb,ispin)**2*(eps/(eps**2 + eta**2))**2
          end do
        end do
      end do
    end do
  end do

  ! Virtual part of the correlation self-energy

  do p=nC+1,nBas-nR
    do a=nO+1,nBas-nR
      jb = 0
      do j=nC+1,nO
        do b=nO+1,nBas-nR
          jb = jb + 1
          do ispin=1,nspin
            eps = e(p) - e(a) - Omega(jb,ispin) 
            Z(p) = Z(p)  - 2d0*rho(p,a,jb,ispin)**2*(eps/(eps**2 + eta**2))**2
          end do
        end do
      end do
    end do
  end do

! Compute renormalization factor from derivative of SigC
 
  Z(:) = 1d0/(1d0 - Z(:))

end subroutine renormalization_factor_SOSEX
