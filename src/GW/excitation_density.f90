subroutine excitation_density(nBas,nC,nO,nR,nS,ERI,XpY,rho)

! Compute excitation densities

  implicit none

! Input variables

  integer,intent(in)            :: nBas,nC,nO,nR,nS
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: XpY(nS,nS)

! Local variables

  integer                       :: ia,jb,p,q,j,b

! Output variables

  double precision,intent(out)  :: rho(nBas,nBas,nS)

  rho(:,:,:) = 0d0   

  do p=nC+1,nBas-nR
    do q=nC+1,nBas-nR
      do ia=1,nS
        jb = 0
        do j=nC+1,nO
          do b=nO+1,nBas-nR
            jb = jb + 1
            rho(p,q,ia) = rho(p,q,ia) + ERI(p,j,q,b)*XpY(ia,jb)
          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine excitation_density
