subroutine unrestricted_excitation_density_Tmatrix(ispin,nBas,nC,nO,nV,nR,nH,nP,ERI_aaaa,ERI_aabb,ERI_bbbb,X1,Y1,rho1,X2,Y2,rho2)

! Compute excitation densities for T-matrix self-energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)
  integer,intent(in)            :: nH
  integer,intent(in)            :: nP 
  double precision,intent(in)   :: X1(nP,nP)
  double precision,intent(in)   :: Y1(nH,nP)
  double precision,intent(in)   :: X2(nP,nH)
  double precision,intent(in)   :: Y2(nH,nH)

! Local variables

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: p,q
  integer                       :: ab,cd,ij,kl
  double precision,external     :: Kronecker_delta

! Output variables

  double precision,intent(out)  :: rho1(nBas,nBas,nP)
  double precision,intent(out)  :: rho2(nBas,nBas,nH)

! Initialization

  rho1(:,:,:) = 0d0
  rho2(:,:,:) = 0d0

!----------------------------------------------
! alpha-beta block
!----------------------------------------------

  if(ispin == 3) then

    do p=nC(1)+1,nBas-nR(1)
      do q=nC(2)+1,nBas-nR(2)
        do ab=1,nP 
        cd = 0
          do c=nO(1)+1,nBas-nR(1)
            do d=nO(2)+1,nBas-nR(1)
              cd = cd + 1

              rho1(p,q,ab) = rho1(p,q,ab) & 
              + (1d0*ERI_aabb(p,q,c,d) + 0d0*ERI_aabb(p,q,d,c))*X1(cd,ab)

            end do
          end do
 
          kl = 0
          do k=nC(1)+1,nO(1)
            do l=nC(1)+1,nO(1)
              kl = kl + 1

              rho1(p,q,ab) = rho1(p,q,ab) & 
              + (1d0*ERI_aabb(p,q,k,l) + 0d0*ERI_aabb(p,q,l,k))*Y1(kl,ab)

            end do
          end do
        end do
        
 
        ij = 0
        do i=nC(1)+1,nO(1)
          do j=nC(2)+1,nO(2)
            ij = ij + 1
 
          cd = 0
          do c=nO(1)+1,nBas-nR(1)
            do d=nO(2)+1,nBas-nR(2)
              cd = cd + 1
              rho2(p,q,ij) = rho2(p,q,ij) &
              + (1d0*ERI_aabb(p,q,c,d) + 0d0*ERI_aabb(p,q,d,c))*X2(cd,ij)

            end do
          end do
 
          kl = 0
          do k=nC(1)+1,nO(1)
            do l=nC(1)+1,nO(1)
              kl = kl + 1
              rho2(p,q,ij) = rho2(p,q,ij) &
              + (1d0*ERI_aabb(p,q,k,l) + 0d0*ERI_aabb(p,q,l,k))*Y2(kl,ij)

            end do
          end do
 
        end do
        end do

      end do
    end do

  end if

!----------------------------------------------
! alpha-alpha block 
!----------------------------------------------

  if(ispin == 4) then

    do p=nC(1)+1,nBas-nR(1)
      do q=nC(1)+1,nBas-nR(1)
 
        do ab=1,nP
 
          cd = 0
          do c=nO(1)+1,nBas-nR(1)
            do d=c+1,nBas-nR(1)
              cd = cd + 1
              rho1(p,q,ab) = rho1(p,q,ab) & 
                           + (ERI_aaaa(p,q,c,d) - ERI_aaaa(p,q,d,c))*X1(cd,ab) 
            end do
          end do
 
          kl = 0
          do k=nC(1)+1,nO(1)
            do l=k+1,nO(1)
              kl = kl + 1
              rho1(p,q,ab) = rho1(p,q,ab) & 
                           + (ERI_aaaa(p,q,k,l) - ERI_aaaa(p,q,l,k))*Y1(kl,ab) 
            end do
          end do
 
        end do
 
        do ij=1,nH
 
          cd = 0
          do c=nO(1)+1,nBas-nR(1)
            do d=c+1,nBas-nR(1)
              cd = cd + 1
              rho2(p,q,ij) = rho2(p,q,ij) & 
                           + (ERI_aaaa(p,q,c,d) - ERI_aaaa(p,q,d,c))*X2(cd,ij) 
            end do
          end do
 
          kl = 0
          do k=nC(1)+1,nO(1)
            do l=k+1,nO(1)
              kl = kl + 1
              rho2(p,q,ij) = rho2(p,q,ij) & 
                           + (ERI_aaaa(p,q,k,l) - ERI_aaaa(p,q,l,k))*Y2(kl,ij) 
            end do
          end do
 
        end do

      end do
    end do

  end if

!----------------------------------------------
! beta-beta block
!----------------------------------------------

  if(ispin == 7) then

    do p=nC(2)+1,nBas-nR(2)
      do q=nC(2)+1,nBas-nR(2)
 
        do ab=1,nP
 
          cd = 0
          do c=nO(2)+1,nBas-nR(2)
            do d=c+1,nBas-nR(2)
              cd = cd + 1
              rho1(p,q,ab) = rho1(p,q,ab) + (ERI_bbbb(p,q,c,d)-ERI_bbbb(p,q,d,c))*X1(cd,ab) 
            end do
          end do
 
          kl = 0
          do k=nC(2)+1,nO(2)
            do l=k+1,nO(2)
              kl = kl + 1
              rho1(p,q,ab) = rho1(p,q,ab) + (ERI_bbbb(p,q,k,l)-ERI_bbbb(p,q,l,k))*Y1(kl,ab) 
            end do
          end do
 
        end do
 
        do ij=1,nH
 
          cd = 0
          do c=nO(2)+1,nBas-nR(2)
            do d=c+1,nBas-nR(2)
              cd = cd + 1
              rho2(p,q,ij) = rho2(p,q,ij) + (ERI_bbbb(p,q,c,d)-ERI_bbbb(p,q,d,c))*X2(cd,ij) 
            end do
          end do
 
          kl = 0
          do k=nC(2)+1,nO(2)
            do l=k+1,nO(2)
              kl = kl + 1
              rho2(p,q,ij) = rho2(p,q,ij) + (ERI_bbbb(p,q,k,l)-ERI_bbbb(p,q,l,k))*Y2(kl,ij) 
            end do
          end do
 
        end do

      end do
    end do

  end if

end subroutine unrestricted_excitation_density_Tmatrix
