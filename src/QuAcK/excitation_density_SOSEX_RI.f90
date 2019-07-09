subroutine excitation_density_SOSEX_RI(nBas,nC,nO,nR,nS,c,G,XpY,rho)

! Compute excitation densities

  implicit none

! Input variables

  integer,intent(in)            :: nBas,nC,nO,nR,nS
  double precision,intent(in)   :: c(nBas,nBas),G(nBas,nBas,nBas,nBas),XpY(nS,nS)

! Local variables

  double precision,allocatable  :: scr(:,:,:)
  integer                       :: mu,nu,la,si,ia,jb,x,y,j,b

! Output variables

  double precision,intent(out)  :: rho(nBas,nBas,nS)

! Memory allocation
  allocate(scr(nBas,nBas,nS))

  rho(:,:,:) = 0d0
  do nu=1,nBas
    do si=1,nBas
      do ia=1,nS
        jb = 0
        do j=nC+1,nO
          do b=nO+1,nBas-nR
            jb = jb + 1
            rho(nu,si,ia) = rho(nu,si,ia) + c(nu,j)*XpY(ia,jb)*c(si,b)
          enddo
        enddo
      enddo
    enddo
  enddo
 
  scr(:,:,:) = 0d0   
  do mu=1,nBas
    do la=1,nBas
      do ia=1,nS
        do nu=1,nBas
          do si=1,nBas
            scr(mu,la,ia) = scr(mu,la,ia) + G(mu,nu,la,si)*rho(nu,si,ia)
          enddo
        enddo
      enddo
    enddo
  enddo

  rho(:,:,:) = 0d0   
  do ia=1,nS
    do x=nC+1,nBas-nR
      do y=nC+1,nBas-nR
        do mu=1,nBas
          do la=1,nBas
            rho(x,y,ia) = rho(x,y,ia) + c(mu,x)*scr(mu,la,ia)*c(la,y)
          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine excitation_density_SOSEX_RI
