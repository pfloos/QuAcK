subroutine CVS_UGW_excitation_density(nBas,nC,nO,nR,nSa,nSb,nSt,nCVS,nFC,occupations,virtuals,ERI_aaaa,ERI_aabb,ERI_bbbb,XpY,rho)

! Compute excitation densities for unrestricted reference

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nSa
  integer,intent(in)            :: nSb
  integer,intent(in)            :: nSt
  integer,intent(in)            :: nCVS(nspin),nFC(nspin)
  integer,intent(in)            :: occupations(maxval(nO-nFC),nspin)
  integer,intent(in)            :: virtuals(nBas-minval(nO),nspin)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: XpY(nSt,nSt)

! Local variables

  integer                       :: ia,jb,p,q,j,b

! Output variables

  double precision,intent(out)  :: rho(nBas,nBas,nSt,nspin)

! Initialization

  rho(:,:,:,:) = 0d0   

!-------------!
! alpha block !
!-------------!

  do p=1,nBas
    do q=1,nBas

      ! Same-spin contribution
      do ia=1,nSt
        jb = 0
        do j=1,nO(1) -nFC(1)
          do b=nCVS(1)+1,nBas -nO(1)
            jb = jb + 1

            rho(p,q,ia,1) = rho(p,q,ia,1) + ERI_aaaa(p,occupations(j,1),q,virtuals(b,1))*XpY(ia,jb)

          end do
        end do
      end do

      ! Opposite-spin contribution
      do ia=1,nSt
        jb = nSa
        do j=1,nO(2) - nFC(2)
          do b=1+nCVS(2),nBas-nO(2)
            jb = jb + 1

            rho(p,q,ia,1) = rho(p,q,ia,1) + ERI_aabb(p,occupations(j,2),q,virtuals(b,2))*XpY(ia,jb)

          end do
        end do
      end do

    end do
  end do

!------------!
! Beta block !
!------------!

  do p=1,nBas
    do q=1,nBas

      ! Opposite-spin contribution
      do ia=1,nSt
        jb = 0
        do j=1,nO(1) - nFC(1)
          do b=nCVS(1)+1,nBas-nO(1)
            jb = jb + 1

            rho(p,q,ia,2) = rho(p,q,ia,2) + ERI_aabb(occupations(j,1),p,virtuals(b,1),q)*XpY(ia,jb)

          end do
        end do
      end do

      ! Same-spin contribution
      do ia=1,nSt
        jb = nSa
        do j=1,nO(2) - nFC(2)
          do b=nCVS(2)+1,nBas-nO(2)
            jb = jb + 1

            rho(p,q,ia,2) = rho(p,q,ia,2) + ERI_bbbb(p,occupations(j,2),q,virtuals(b,2))*XpY(ia,jb)

          end do
        end do
      end do

    end do
  end do

end subroutine 
