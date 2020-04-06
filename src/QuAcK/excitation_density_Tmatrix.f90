subroutine excitation_density_Tmatrix(nBas,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,ERI,rho1st,rho2st, & 
                                      X1s,Y1s,rho1s,X2s,Y2s,rho2s,X1t,Y1t,rho1t,X2t,Y2t,rho2t)

! Compute excitation densities for T-matrix self-energy

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  integer,intent(in)            :: nOOs
  integer,intent(in)            :: nOOt
  double precision,intent(in)   :: X1s(nVVs,nVVs)
  double precision,intent(in)   :: Y1s(nOOs,nVVs)
  double precision,intent(in)   :: X2s(nVVs,nOOs)
  double precision,intent(in)   :: Y2s(nOOs,nOOs)
  integer,intent(in)            :: nVVs
  integer,intent(in)            :: nVVt
  double precision,intent(in)   :: X1t(nVVt,nVVt)
  double precision,intent(in)   :: Y1t(nOOt,nVVt)
  double precision,intent(in)   :: X2t(nVVt,nOOt)
  double precision,intent(in)   :: Y2t(nOOt,nOOt)

! Local variables

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: p
  integer                       :: ab,cd,ij,kl
  double precision,external     :: Kronecker_delta

! Output variables

  double precision,intent(out)  :: rho1s(nBas,nO,nVVs)
  double precision,intent(out)  :: rho2s(nBas,nV,nOOs)
  double precision,intent(out)  :: rho1t(nBas,nO,nVVt)
  double precision,intent(out)  :: rho2t(nBas,nV,nOOt)
  double precision,intent(out)  :: rho1st(nBas,nO,nVVt)
  double precision,intent(out)  :: rho2st(nBas,nV,nOOt)

! Initialization

  rho1s(:,:,:) = 0d0
  rho2s(:,:,:) = 0d0
  rho1t(:,:,:) = 0d0
  rho2t(:,:,:) = 0d0
  rho1st(:,:,:) = 0d0
  rho2st(:,:,:) = 0d0

!----------------------------------------------
! Singlet manifold
!----------------------------------------------

  do p=nC+1,nBas-nR

    do i=nC+1,nO
      do ab=1,nVVs

        cd = 0
        do c=nO+1,nBas-nR
          do d=c,nBas-nR
            cd = cd + 1
            rho1s(p,i,ab) = rho1s(p,i,ab) + 2.0d0*ERI(p,i,c,d)*X1s(cd,ab) 
          end do
        end do

        kl = 0
        do k=nC+1,nO
          do l=k,nO
            kl = kl + 1
            rho1s(p,i,ab) = rho1s(p,i,ab) + 2.0d0*ERI(p,i,k,l)*Y1s(kl,ab) 
          end do
        end do

      end do
    end do

    do a=1,nV-nR
      do ij=1,nOOs

        cd = 0
        do c=nO+1,nBas-nR
          do d=c,nBas-nR
            cd = cd + 1
            rho2s(p,a,ij) = rho2s(p,a,ij) + 2.0d0*ERI(p,nO+a,c,d)*X2s(cd,ij)
          end do
        end do

        kl = 0
        do k=nC+1,nO
          do l=k,nO
            kl = kl + 1
            rho2s(p,a,ij) = rho2s(p,a,ij) + 2.0d0*ERI(p,nO+a,k,l)*Y2s(kl,ij) 
          end do
        end do

      end do
    end do

  end do

!----------------------------------------------
! Triplet manifold
!----------------------------------------------

  do p=nC+1,nBas-nR

    do i=nC+1,nO
      do ab=1,nVVt

        cd = 0
        do c=nO+1,nBas-nR
          do d=c+1,nBas-nR
            cd = cd + 1
            rho1t(p,i,ab) = rho1t(p,i,ab) + 1d0*(ERI(p,i,c,d) - ERI(p,i,d,c))*X1t(cd,ab) 
          end do
        end do

        kl = 0
        do k=nC+1,nO
          do l=k+1,nO
            kl = kl + 1
            rho1t(p,i,ab) = rho1t(p,i,ab) + 1d0*(ERI(p,i,k,l) - ERI(p,i,l,k))*Y1t(kl,ab) 
          end do
        end do

      end do
    end do

    do a=1,nV-nR
      do ij=1,nOOt

        cd = 0
        do c=nO+1,nBas-nR
          do d=c+1,nBas-nR
            cd = cd + 1
            rho2t(p,a,ij) = rho2t(p,a,ij) + 1d0*(ERI(p,nO+a,c,d) - ERI(p,nO+a,d,c))*X2t(cd,ij) 
          end do
        end do

        kl = 0
        do k=nC+1,nO
          do l=k+1,nO
            kl = kl + 1
            rho2t(p,a,ij) = rho2t(p,a,ij) + 1d0*(ERI(p,nO+a,k,l) - ERI(p,nO+a,l,k))*Y2t(kl,ij) 
          end do
        end do

      end do
    end do

  end do

!----------------------------------------------
! Singlet-triplet crossed term
!----------------------------------------------

  do p=nC+1,nBas-nR

    do i=nC+1,nO
      do ab=1,nVVt

        cd = 0
        do c=nO+1,nBas-nR
          do d=c+1,nBas-nR
            cd = cd + 1
            rho1st(p,i,ab) = rho1st(p,i,ab) + 2d0*ERI(p,i,c,d)*X1t(cd,ab) 
          end do
        end do

        kl = 0
        do k=nC+1,nO
          do l=k+1,nO
            kl = kl + 1
            rho1st(p,i,ab) = rho1st(p,i,ab) + 2d0*ERI(p,i,k,l)*Y1t(kl,ab) 
          end do
        end do

      end do
    end do

    do a=1,nV-nR
      do ij=1,nOOt

        cd = 0
        do c=nO+1,nBas-nR
          do d=c+1,nBas-nR
            cd = cd + 1
            rho2st(p,a,ij) = rho2st(p,a,ij) + 2d0*ERI(p,nO+a,c,d)*X2t(cd,ij) 
          end do
        end do

        kl = 0
        do k=nC+1,nO
          do l=k+1,nO
            kl = kl + 1
            rho2st(p,a,ij) = rho2st(p,a,ij) + 2d0*ERI(p,nO+a,k,l)*Y2t(kl,ij) 
          end do
        end do

      end do
    end do

  end do

end subroutine excitation_density_Tmatrix
