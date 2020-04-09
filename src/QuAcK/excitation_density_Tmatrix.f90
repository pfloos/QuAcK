subroutine excitation_density_Tmatrix(ispin,dERI,xERI,nBas,nC,nO,nV,nR,nOO,nVV,ERI,X1,Y1,rho1,X2,Y2,rho2)

! Compute excitation densities for T-matrix self-energy

  implicit none

! Input variables

  integer,intent(in)            :: ispin
  double precision,intent(in)   :: dERI
  double precision,intent(in)   :: xERI
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV 
  double precision,intent(in)   :: X1(nVV,nVV)
  double precision,intent(in)   :: Y1(nOO,nVV)
  double precision,intent(in)   :: X2(nVV,nOO)
  double precision,intent(in)   :: Y2(nOO,nOO)

! Local variables

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: p
  integer                       :: ab,cd,ij,kl
  double precision,external     :: Kronecker_delta

! Output variables

  double precision,intent(out)  :: rho1(nBas,nO,nVV)
  double precision,intent(out)  :: rho2(nBas,nV,nOO)

! Initialization

  rho1(:,:,:) = 0d0
  rho2(:,:,:) = 0d0

!----------------------------------------------
! Singlet manifold
!----------------------------------------------

  if(ispin == 1) then

    do p=nC+1,nBas-nR
 
      do i=nC+1,nO
        do ab=1,nVV
 
          cd = 0
          do c=nO+1,nBas-nR
            do d=c,nBas-nR
              cd = cd + 1
              rho1(p,i,ab) = rho1(p,i,ab) & 
                           + (dERI*ERI(p,i,c,d) + xERI*ERI(p,i,d,c))*X1(cd,ab) 
            end do
          end do
 
          kl = 0
          do k=nC+1,nO
            do l=k,nO
              kl = kl + 1
              rho1(p,i,ab) = rho1(p,i,ab) & 
                           + (dERI*ERI(p,i,k,l) + xERI*ERI(p,i,l,k))*Y1(kl,ab) 
            end do
          end do
 
        end do
      end do
 
      do a=1,nV-nR
        do ij=1,nOO
 
          cd = 0
          do c=nO+1,nBas-nR
            do d=c,nBas-nR
              cd = cd + 1
              rho2(p,a,ij) = rho2(p,a,ij) & 
                           + (dERI*ERI(p,nO+a,c,d) + xERI*ERI(p,nO+a,d,c))*X2(cd,ij)
            end do
          end do
 
          kl = 0
          do k=nC+1,nO
            do l=k,nO
              kl = kl + 1
              rho2(p,a,ij) = rho2(p,a,ij) & 
                           + (dERI*ERI(p,nO+a,k,l) + xERI*ERI(p,nO+a,l,k))*Y2(kl,ij) 
            end do
          end do
 
        end do
      end do
 
    end do

  end if

!----------------------------------------------
! Triplet manifold
!----------------------------------------------

  if(ispin == 2) then

    do p=nC+1,nBas-nR
 
      do i=nC+1,nO
        do ab=1,nVV
 
          cd = 0
          do c=nO+1,nBas-nR
            do d=c+1,nBas-nR
              cd = cd + 1
              rho1(p,i,ab) = rho1(p,i,ab) & 
                           + (dERI*ERI(p,i,c,d) + xERI*ERI(p,i,d,c))*X1(cd,ab) 
            end do
          end do
 
          kl = 0
          do k=nC+1,nO
            do l=k+1,nO
              kl = kl + 1
              rho1(p,i,ab) = rho1(p,i,ab) & 
                           + (dERI*ERI(p,i,k,l) + xERI*ERI(p,i,l,k))*Y1(kl,ab) 
            end do
          end do
 
        end do
      end do
 
      do a=1,nV-nR
        do ij=1,nOO
 
          cd = 0
          do c=nO+1,nBas-nR
            do d=c+1,nBas-nR
              cd = cd + 1
              rho2(p,a,ij) = rho2(p,a,ij) & 
                           + (dERI*ERI(p,nO+a,c,d) + xERI*ERI(p,nO+a,d,c))*X2(cd,ij) 
            end do
          end do
 
          kl = 0
          do k=nC+1,nO
            do l=k+1,nO
              kl = kl + 1
              rho2(p,a,ij) = rho2(p,a,ij) & 
                           + (dERI*ERI(p,nO+a,k,l) + xERI*ERI(p,nO+a,l,k))*Y2(kl,ij) 
            end do
          end do
 
        end do
      end do
 
    end do

  end if

end subroutine excitation_density_Tmatrix
