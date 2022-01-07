subroutine S51_lda_exchange_potential(nGrid,weight,nBas,AO,rho,Fx)

! Compute Slater's LDA exchange potential

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: rho(nGrid)

! Local variables

  integer                       :: mu,nu,iG
  double precision              :: r,vAO

! Output variables

  double precision,intent(out)  :: Fx(nBas,nBas)

! Compute LDA exchange matrix in the AO basis

  Fx(:,:) = 0d0
  do mu=1,nBas
    do nu=1,nBas
      do iG=1,nGrid

        r = max(0d0,rho(iG))

        if(r > threshold) then

          vAO = weight(iG)*AO(mu,iG)*AO(nu,iG)
          Fx(mu,nu) = Fx(mu,nu) + vAO*4d0/3d0*CxLSDA*r**(1d0/3d0)

        endif

      enddo
    enddo
  enddo

end subroutine S51_lda_exchange_potential
