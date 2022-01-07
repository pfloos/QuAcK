subroutine xc_potential_grid(nBas,nGrid,AO,rho,Fx,Vxgrid)


! Compute the exchange-correlation potential on the grid
                                                                                                                       
  implicit none
  include 'parameters.h'                                                                                               

! Input variables

  integer,intent(in)            :: nBas, nGrid
  double precision,intent(in)   :: rho(nGrid,nspin)
  double precision,intent(in)   :: Fx(nBas,nBas,nspin)
  double precision,intent(in)   :: AO(nBas,nGrid)

! Local variables

  integer                       :: mu,nu
  integer                       :: ispin,iG
  double precision              :: r 
  double precision              :: Fxgrid(nGrid,nspin)

! Output variables

  double precision,intent(out)  :: Vxgrid(nGrid)

! Compute Vx 

  Vxgrid(:) = 0d0
  Fxgrid(:,:) = 0d0
 
  do iG=1,nGrid 
    do ispin=1,nspin 
      do mu=1,nBas
        do nu=1,nBas
          r = max(0d0,rho(iG,ispin))
          if(r > threshold) then
          Fxgrid(iG,ispin) = Fxgrid(iG,ispin) + AO(mu,iG)*AO(nu,iG)*4d0/3d0*CxLSDA*r**(1d0/3d0) 
          endif
        enddo
      enddo
    enddo
  enddo

  Vxgrid(:)=Fxgrid(:,1)+Fxgrid(:,2)
  open(411, file = 'Vxgrid.dat', status = 'new') 
  do iG=1,nGrid
    write(411,*) iG, Vxgrid(iG) 
  end do 
  close(411) 


end subroutine xc_potential_grid

