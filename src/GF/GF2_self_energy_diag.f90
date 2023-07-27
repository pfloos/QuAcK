subroutine GF2_self_energy_diag(eta,nBas,nC,nO,nV,nR,eHF,eGF2,ERI,SigC,Z)

! Compute diagonal part of the GF2 self-energy and its renormalization factor

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: eGF2(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: i,j,a,b
  integer                       :: p
  double precision              :: eps
  double precision              :: num

! Output variables

  double precision,intent(out)  :: SigC(nBas)
  double precision,intent(out)  :: Z(nBas)

! Initialize 

  SigC(:) = 0d0
  Z(:)    = 0d0

! Compute GF2 self-energy

  do p=nC+1,nBas-nR
    do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nBas-nR

          eps = eGF2(p) + eHF(a) - eHF(i) - eHF(j)
          num = (2d0*ERI(p,a,i,j) - ERI(p,a,j,i))*ERI(p,a,i,j)

          SigC(p) = SigC(p) + num*eps/(eps**2 + eta**2)
          Z(p)    = Z(p)   - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

        end do
      end do
    end do
  end do

  do p=nC+1,nBas-nR
    do i=nC+1,nO
      do a=nO+1,nBas-nR
        do b=nO+1,nBas-nR

          eps = eGF2(p) + eHF(i) - eHF(a) - eHF(b)
          num = (2d0*ERI(p,i,a,b) - ERI(p,i,b,a))*ERI(p,i,a,b)

          SigC(p) = SigC(p) + num*eps/(eps**2 + eta**2)
          Z(p)    = Z(p)   - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

        end do
      end do
    end do
  end do

  Z(:) = 1d0/(1d0 - Z(:))

end subroutine 
