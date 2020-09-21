subroutine unrestricted_excitation_density(nBas,nC,nO,nR,nSa,nSb,nSt,ERI_aa,ERI_ab,ERI_bb,XpY_a,XpY_b,rho)

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
  double precision,intent(in)   :: ERI_aa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_ab(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: XpY_a(nSa,nSa)
  double precision,intent(in)   :: XpY_b(nSb,nSb)

! Local variables

  integer                       :: ia,jb,p,q,j,b

! Output variables

  double precision,intent(out)  :: rho(nBas,nBas,nSt,nspin)

! Initialization

  rho(:,:,:,:) = 0d0   

!-------------
! alpha block
!-------------

  do p=nC(1)+1,nBas-nR(1)
    do q=nC(1)+1,nBas-nR(1)

      ! Same-spin contribution
      do ia=1,nSa
        jb = 0
        do j=nC(1)+1,nO(1)
          do b=nO(1)+1,nBas-nR(1)
            jb = jb + 1

            rho(p,q,ia,1) = rho(p,q,ia,1) + ERI_aa(p,j,q,b)*XpY_a(ia,jb)

          enddo
        enddo
      enddo

      ! Opposite-spin contribution
      do ia=1,nSb
        jb = 0
        do j=nC(2)+1,nO(2)
          do b=nO(2)+1,nBas-nR(2)
            jb = jb + 1

            rho(p,q,nSa+ia,1) = rho(p,q,nSa+ia,1) + ERI_ab(p,j,q,b)*XpY_b(ia,jb)

          enddo
        enddo
      enddo

    enddo
  enddo

!------------
! Beta block
!------------

  do p=nC(2)+1,nBas-nR(2)
    do q=nC(2)+1,nBas-nR(2)

      ! Same-spin contribution
      do ia=1,nSb
        jb = 0
        do j=nC(2)+1,nO(2)
          do b=nO(2)+1,nBas-nR(2)
            jb = jb + 1

            rho(p,q,ia,2) = rho(p,q,ia,2) + ERI_bb(p,j,q,b)*XpY_b(ia,jb)

          enddo
        enddo
      enddo

      ! Opposite-spin contribution
      do ia=1,nSa
        jb = 0
        do j=nC(1)+1,nO(1)
          do b=nO(1)+1,nBas-nR(1)
            jb = jb + 1

            rho(p,q,nSb+ia,2) = rho(p,q,nSb+ia,2) + ERI_ab(j,p,b,q)*XpY_a(ia,jb)

          enddo
        enddo
      enddo

    enddo
  enddo

end subroutine unrestricted_excitation_density
