subroutine excitation_density_Tmatrix(ispin,nBas,nC,nO,nV,nR,nOO,nVV,ERI,X1,Y1,rho1,X2,Y2,rho2)

! Compute excitation densities for T-matrix self-energy

  implicit none

! Input variables

  integer,intent(in)            :: ispin
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
  integer                       :: p,q
  integer                       :: ab,cd,ij,kl
  double precision,external     :: Kronecker_delta

! Output variables

  double precision,intent(out)  :: rho1(nBas,nBas,nVV)
  double precision,intent(out)  :: rho2(nBas,nBas,nOO)

! Initialization

  rho1(:,:,:) = 0d0
  rho2(:,:,:) = 0d0

!----------------------------------------------
! Singlet manifold
!----------------------------------------------

  if(ispin == 1) then

    do p=nC+1,nBas-nR
      do q=nC+1,nBas-nR

!       do ab=1,nVV
        ab = 0
        do a=nO+1,nBas-nR
          do b=a,nBas-nR
            ab = ab + 1
 
          cd = 0
          do c=nO+1,nBas-nR
!           do d=nO+1,c
            do d=c,nBas-nR
              cd = cd + 1
              rho1(p,q,ab) = rho1(p,q,ab) & 
              + (1d0*ERI(p,q,c,d) + 0d0*ERI(p,q,d,c))*X1(cd,ab)
!                          + ERI(p,q,c,d)*X1(cd,ab)/sqrt((1d0 + Kronecker_delta(c,d)))
!                          + (ERI(p,q,c,d) + ERI(p,q,d,c))*X1(cd,ab)/sqrt((1d0 + Kronecker_delta(a,b))*(1d0 + Kronecker_delta(c,d)))
            end do
          end do
 
          kl = 0
          do k=nC+1,nO
!           do l=nC+1,k
            do l=k,nO
              kl = kl + 1
              rho1(p,q,ab) = rho1(p,q,ab) & 
              + (1d0*ERI(p,q,k,l) + 0d0*ERI(p,q,l,k))*Y1(kl,ab)
!                          + ERI(p,q,k,l)*Y1(kl,ab)/sqrt((1d0 + Kronecker_delta(k,l)))
!                          + (ERI(p,q,k,l) + ERI(p,q,l,k))*Y1(kl,ab)/sqrt((1d0 + Kronecker_delta(a,b))*(1d0 + Kronecker_delta(k,l)))
            end do
          end do
 
        end do
        end do
 
!       do ij=1,nOO
        ij = 0
        do i=nC+1,nO
          do j=i,nO
            ij = ij + 1
 
          cd = 0
          do c=nO+1,nBas-nR
!           do d=nO+1,c
            do d=c,nBas-nR
              cd = cd + 1
              rho2(p,q,ij) = rho2(p,q,ij) &
              + (1d0*ERI(p,q,c,d) + 0d0*ERI(p,q,d,c))*X2(cd,ij)
!                          + ERI(p,q,c,d)*X2(cd,ij)/sqrt((1d0 + Kronecker_delta(c,d)))
!                          + (ERI(p,q,c,d) + ERI(p,q,d,c))*X2(cd,ij)/sqrt((1d0 + Kronecker_delta(i,j))*(1d0 + Kronecker_delta(c,d)))
            end do
          end do
 
          kl = 0
          do k=nC+1,nO
!           do l=nC+1,k
            do l=k,nO
              kl = kl + 1
              rho2(p,q,ij) = rho2(p,q,ij) &
              + (1d0*ERI(p,q,k,l) + 0d0*ERI(p,q,l,k))*Y2(kl,ij)
!                          + ERI(p,q,k,l)*Y2(kl,ij)/sqrt((1d0 + Kronecker_delta(k,l)))
!                          + (ERI(p,q,k,l) + ERI(p,q,l,k))*Y2(kl,ij)/sqrt((1d0 + Kronecker_delta(i,j))*(1d0 + Kronecker_delta(k,l)))
            end do
          end do
 
        end do
        end do

      end do
    end do

  end if

!----------------------------------------------
! Triplet manifold
!----------------------------------------------

  if(ispin == 2 .or. ispin == 4) then

    do p=nC+1,nBas-nR
      do q=nC+1,nBas-nR
 
        do ab=1,nVV
 
          cd = 0
          do c=nO+1,nBas-nR
            do d=c+1,nBas-nR
              cd = cd + 1
              rho1(p,q,ab) = rho1(p,q,ab) & 
                           + (ERI(p,q,c,d) - ERI(p,q,d,c))*X1(cd,ab) 
            end do
          end do
 
          kl = 0
          do k=nC+1,nO
            do l=k+1,nO
              kl = kl + 1
              rho1(p,q,ab) = rho1(p,q,ab) & 
                           + (ERI(p,q,k,l) - ERI(p,q,l,k))*Y1(kl,ab) 
            end do
          end do
 
        end do
 
        do ij=1,nOO
 
          cd = 0
          do c=nO+1,nBas-nR
            do d=c+1,nBas-nR
              cd = cd + 1
              rho2(p,q,ij) = rho2(p,q,ij) & 
                           + (ERI(p,q,c,d) - ERI(p,q,d,c))*X2(cd,ij) 
            end do
          end do
 
          kl = 0
          do k=nC+1,nO
            do l=k+1,nO
              kl = kl + 1
              rho2(p,q,ij) = rho2(p,q,ij) & 
                           + (ERI(p,q,k,l) - ERI(p,q,l,k))*Y2(kl,ij) 
            end do
          end do
 
        end do

      end do
    end do

  end if

!----------------------------------------------
! alpha-beta block
!----------------------------------------------

  if(ispin == 3) then

    do p=nC+1,nBas-nR
      do q=nC+1,nBas-nR
 
        do ab=1,nVV
 
          cd = 0
          do c=nO+1,nBas-nR
            do d=nO+1,nBas-nR
              cd = cd + 1
              rho1(p,q,ab) = rho1(p,q,ab) + ERI(p,q,c,d)*X1(cd,ab) 
            end do
          end do
 
          kl = 0
          do k=nC+1,nO
            do l=nC+1,nO
              kl = kl + 1
              rho1(p,q,ab) = rho1(p,q,ab) + ERI(p,q,k,l)*Y1(kl,ab) 
            end do
          end do
 
        end do
 
        do ij=1,nOO
 
          cd = 0
          do c=nO+1,nBas-nR
            do d=nO+1,nBas-nR
              cd = cd + 1
              rho2(p,q,ij) = rho2(p,q,ij) + ERI(p,q,c,d)*X2(cd,ij) 
            end do
          end do
 
          kl = 0
          do k=nC+1,nO
            do l=nC+1,nO
              kl = kl + 1
              rho2(p,q,ij) = rho2(p,q,ij) + ERI(p,q,k,l)*Y2(kl,ij) 
            end do
          end do
 
        end do

      end do
    end do

  end if

end subroutine excitation_density_Tmatrix
