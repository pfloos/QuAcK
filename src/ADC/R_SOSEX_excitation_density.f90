subroutine R_SOSEX_excitation_density(nOrb,nC,nO,nR,nS,e,Om,ERI,XpY,rhoL,rhoR)

! Compute excitation densities for SOSEX

  implicit none

! Input variables

  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: e(nOrb)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: XpY(nS,nS)

! Local variables

  integer                       :: ia,jb,p,q,j,a,b,k,c,m
  double precision, allocatable :: tmp(:,:,:)

! Output variables

  double precision,intent(out)  :: rhoL(nOrb,nOrb,nS)
  double precision,intent(out)  :: rhoR(nOrb,nOrb,nS)

!--------------------------!
! Left effective integrals !
!--------------------------!

  rhoR(:,:,:) = 0d0   
  !$OMP PARALLEL &
  !$OMP SHARED(nC,nOrb,nR,nO,nS,rhoR,ERI,XpY) &
  !$OMP PRIVATE(q,p,jb,ia) &
  !$OMP DEFAULT(NONE)
  !$OMP DO
  do q=nC+1,nOrb-nR
     do p=nC+1,nOrb-nR
        jb = 0
        do j=nC+1,nO
           do b=nO+1,nOrb-nR
              jb = jb + 1
              do ia=1,nS
                 rhoR(p,q,ia) = rhoR(p,q,ia) + ERI(p,j,q,b)*XpY(ia,jb)
              end do
           end do
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

!---------------------------!
! Right effective integrals !
!---------------------------!

  rhoL(:,:,:) = rhoR(:,:,:)

  do p=nC+1,nOrb-nR
    do k=nC+1,nO
      do m=1,nS

        do j=nC+1,nO
          do a=nO+1,nOrb-nR

            rhoL(p,k,m) = rhoL(p,k,m) + rhoR(a,j,m)*(ERI(p,j,a,k)/(e(a) - e(j) + Om(m)) + ERI(p,a,j,k)/(e(a) - e(j) - Om(m)))

          end do
        end do

      end do
    end do
  end do

  do p=nC+1,nOrb-nR
    do c=nO+1,nOrb-nR
      do m=1,nS

        do j=nC+1,nO
          do a=nO+1,nOrb-nR

            rhoL(p,c,m) = rhoL(p,c,m) + rhoR(a,j,m)*(ERI(p,j,a,c)/(e(a) - e(j) - Om(m)) + ERI(p,a,j,c)/(e(a) - e(j) + Om(m)))

          end do
        end do

      end do
    end do
  end do

end subroutine 
