subroutine excitation_density_SOSEX(nBas,nC,nO,nR,nS,ERI,XpY,rho)

! Compute excitation densities for SOSEX

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: XpY(nS,nS,nspin)

! Local variables

  integer                       :: ispin
  integer                       :: p,q
  integer                       :: i,a,j,b
  integer                       :: ia,jb

! Output variables

  double precision,intent(out)  :: rho(nBas,nBas,nS,nspin)

  rho(:,:,:,:) = 0d0

! Singlet part

  ispin = 1

  do ia=1,nS
    do p=nC+1,nBas-nR
      do q=nC+1,nBas-nR
        jb = 0
        do j=nC+1,nO
          do b=nO+1,nBas-nR
            jb = jb + 1
            rho(p,q,ia,ispin) = rho(p,q,ia,ispin) + ERI(p,j,q,b)*XpY(ia,jb,ispin)
          enddo
        enddo
      enddo
    enddo
  enddo

! Triplet part

  ispin = 2

  do ia=1,nS
    do p=nC+1,nBas-nR
      do q=nC+1,nBas-nR
        jb = 0
        do j=nC+1,nO
          do b=nO+1,nBas-nR
            jb = jb + 1
            rho(p,q,ia,ispin) = rho(p,q,ia,ispin) + (ERI(p,j,q,b) - ERI(p,j,b,q))*XpY(ia,jb,ispin)
          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine excitation_density_SOSEX
