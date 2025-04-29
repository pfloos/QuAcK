subroutine complex_RGW_excitation_density(nOrb,nC,nO,nR,nS,ERI,XpY,rho)

! Compute excitation densities

  implicit none

! Input variables

  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  complex*16,intent(in)         :: ERI(nOrb,nOrb,nOrb,nOrb)
  complex*16,intent(in)         :: XpY(nS,nS)

! Local variables

  integer                       :: ia,jb,p,q,j,b

! Output variables

  complex*16,intent(out)        :: rho(nOrb,nOrb,nS)
  
  rho(:,:,:) = 0d0   
  !$OMP PARALLEL &
  !$OMP SHARED(nC,nOrb,nR,nO,nS,rho,ERI,XpY) &
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
                 rho(p,q,ia) = rho(p,q,ia) + ERI(p,j,q,b)*XpY(ia,jb)
              end do
           end do
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
end subroutine 
