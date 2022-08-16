subroutine excitation_density_Tmatrix_so(nBas,nC,nO,nV,nR,nOO,nVV,ERI,X1,Y1,rho1,X2,Y2,rho2)

! Compute excitation densities for T-matrix self-energy

  implicit none

! Input variables

  integer,intent(in)            :: nBas,nC,nO,nV,nR,nOO,nVV
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: X1(nVV,nVV)
  double precision,intent(in)   :: Y1(nOO,nVV)
  double precision,intent(in)   :: X2(nVV,nOO)
  double precision,intent(in)   :: Y2(nOO,nOO)

! Local variables

  integer                       :: j,k,l
  integer                       :: b,c,d
  integer                       :: p,q
  integer                       :: ab,cd,ij,kl
  double precision,external     :: Kronecker_delta

! Output variables

  double precision,intent(out)  :: rho1(nBas,nBas,nVV)
  double precision,intent(out)  :: rho2(nBas,nBAs,nOO)

! Initialization

  rho1(:,:,:) = 0d0
  rho2(:,:,:) = 0d0

  do p=nC+1,nBas-nR
    do q=nC+1,nBas-nR

      do ab=1,nVV

        cd = 0
        do c=nO+1,nBas-nR
          do d=c+1,nBas-nR
            cd = cd + 1
            rho1(p,q,ab) = rho1(p,q,ab) + (ERI(p,q,c,d) - ERI(p,q,d,c))*X1(cd,ab)
          end do
        end do

        kl = 0
        do k=nC+1,nO
          do l=k+1,nO
            kl = kl + 1
            rho1(p,q,ab) = rho1(p,q,ab) + (ERI(p,q,k,l) - ERI(p,q,l,k))*Y1(kl,ab)
          end do
        end do

      end do

      do ij=1,nOO

        cd = 0
        do c=nO+1,nBas-nR
          do d=c+1,nBas-nR
            cd = cd + 1
            rho2(p,q,ij) = rho2(p,q,ij) + (ERI(p,q,c,d) - ERI(p,q,d,c))*X2(cd,ij)
          end do
        end do

        kl = 0
        do k=nC+1,nO
          do l=k+1,nO
            kl = kl + 1
            rho2(p,q,ij) = rho2(p,q,ij) + (ERI(p,q,k,l) - ERI(p,q,l,k))*Y2(kl,ij)
          end do
        end do

      end do
    end do

  end do

end subroutine excitation_density_Tmatrix_so