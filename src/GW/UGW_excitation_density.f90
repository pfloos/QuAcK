subroutine UGW_excitation_density(nBas,nC,nO,nR,nSa,nSb,nSt,ERI_aaaa,ERI_aabb,ERI_bbbb,XpY,rho)

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

  do p=nC(1)+1,nBas-nR(1)
    do q=nC(1)+1,nBas-nR(1)

      ! Same-spin contribution
      do ia=1,nSt
        jb = 0
        do j=nC(1)+1,nO(1)
          do b=nO(1)+1,nBas-nR(1)
            jb = jb + 1

            rho(p,q,ia,1) = rho(p,q,ia,1) + ERI_aaaa(p,j,q,b)*XpY(ia,jb)

          enddo
        enddo
      enddo

      ! Opposite-spin contribution
      do ia=1,nSt
        jb = nSa
        do j=nC(2)+1,nO(2)
          do b=nO(2)+1,nBas-nR(2)
            jb = jb + 1

            rho(p,q,ia,1) = rho(p,q,ia,1) + ERI_aabb(p,j,q,b)*XpY(ia,jb)

          enddo
        enddo
      enddo

    enddo
  enddo

!------------!
! Beta block !
!------------!

  do p=nC(2)+1,nBas-nR(2)
    do q=nC(2)+1,nBas-nR(2)

      ! Opposite-spin contribution
      do ia=1,nSt
        jb = 0
        do j=nC(1)+1,nO(1)
          do b=nO(1)+1,nBas-nR(1)
            jb = jb + 1

            rho(p,q,ia,2) = rho(p,q,ia,2) + ERI_aabb(j,p,b,q)*XpY(ia,jb)

          enddo
        enddo
      enddo

      ! Same-spin contribution
      do ia=1,nSt
        jb = nSa
        do j=nC(2)+1,nO(2)
          do b=nO(2)+1,nBas-nR(2)
            jb = jb + 1

            rho(p,q,ia,2) = rho(p,q,ia,2) + ERI_bbbb(p,j,q,b)*XpY(ia,jb)

          enddo
        enddo
      enddo

    enddo
  enddo

end subroutine 
