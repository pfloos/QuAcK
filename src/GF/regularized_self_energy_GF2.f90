subroutine regularized_self_energy_GF2(eta,nBas,nC,nO,nV,nR,nS,eHF,eGF2,ERI,SigC,Z)

! Compute GF2 self-energy and its renormalization factor

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: eGF2(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: i,j,a,b
  integer                       :: p,q
  double precision              :: eps
  double precision              :: num

  double precision              :: kappa
  double precision              :: fk,dfk

! Output variables

  double precision,intent(out)  :: SigC(nBas,nBas)
  double precision,intent(out)  :: Z(nBas)

! Initialize 

  SigC(:,:) = 0d0
  Z(:)      = 0d0

!-----------------------------------------!
! Parameters for regularized calculations !
!-----------------------------------------!

  kappa = 1d0

!----------------------------------------------------!
! Compute GF2 self-energy and renormalization factor !
!----------------------------------------------------!

  do p=nC+1,nBas-nR
    do q=nC+1,nBas-nR
      do i=nC+1,nO
        do j=nC+1,nO
          do a=nO+1,nBas-nR

            eps = eGF2(p) + eHF(a) - eHF(i) - eHF(j)
            num = (2d0*ERI(p,a,i,j) - ERI(p,a,j,i))*ERI(q,a,i,j)

            fk  = (1d0 - exp(-2d0*eps**2/kappa**2))/eps
            dfk = - fk/eps + 4d0*kappa**2*exp(-2d0*eps**2/kappa**2)

            SigC(p,q) = SigC(p,q) + num*fk
            if(p == q) Z(p) = Z(p) - num*dfk

          end do
        end do
      end do
    end do
  end do

  do p=nC+1,nBas-nR
    do q=nC+1,nBas-nR
      do i=nC+1,nO
        do a=nO+1,nBas-nR
          do b=nO+1,nBas-nR

            eps = eGF2(p) + eHF(i) - eHF(a) - eHF(b)
            num = (2d0*ERI(p,i,a,b) - ERI(p,i,b,a))*ERI(q,i,a,b)

             fk  = (1d0 - exp(-2d0*eps**2/kappa**2))/eps
             dfk = - fk/eps + 4d0*kappa**2*exp(-2d0*eps**2/kappa**2)

            SigC(p,q) = SigC(p,q) + num*fk
            if(p == q) Z(p) = Z(p) - num*dfk

          end do
        end do
      end do
    end do
  end do

  Z(:) = 1d0/(1d0 - Z(:))

end subroutine
