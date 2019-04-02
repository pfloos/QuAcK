subroutine excitation_density_SOSEX(nBas,nC,nO,nR,nS,G,XpY,rho)

! Compute excitation densities

  implicit none

! Input variables

  integer,intent(in)            :: nBas,nC,nO,nR,nS
  double precision,intent(in)   :: G(nBas,nBas,nBas,nBas),XpY(nS,nS)

! Local variables

  integer                       :: ia,jb,x,y,j,b

! Output variables

  double precision,intent(out)  :: rho(nBas,nBas,nS)

  rho(:,:,:) = 0d0
  do ia=1,nS
    do x=nC+1,nBas-nR
      do y=nC+1,nBas-nR
        jb = 0
        do j=nC+1,nO
          do b=nO+1,nBas-nR
            jb = jb + 1
            rho(x,y,ia) = rho(x,y,ia) + G(x,y,b,j)*XpY(ia,jb)
          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine excitation_density_SOSEX
