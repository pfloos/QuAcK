subroutine UW38_lda_correlation_potential(nGrid,weight,nBas,AO,rho,Fc)

! Compute the unrestricted version of the Wigner's LDA correlation potential

  implicit none
include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: rho(nGrid,nspin)

! Local variables

  integer                       :: mu,nu,iG
  double precision              :: ra,rb,r
  double precision              :: a,d,ec
  double precision              :: dFcdr

! Output variables

  double precision,intent(out)  :: Fc(nBas,nBas,nspin)

! Coefficients for Wigner's LDA correlation

  a = 0.04918d0
  d = 0.349d0

! Compute LDA correlation matrix in the AO basis

  Fc(:,:,:) = 0d0

  do mu=1,nBas
    do nu=1,nBas
      do iG=1,nGrid

        ra = max(0d0,rho(iG,1))
        rb = max(0d0,rho(iG,2))
      
!       Spin-up part contribution

        if(ra > threshold) then

          r = ra + rb

          ec    = ra*rb/(r + d*r**(2d0/3d0))
          dFcdr = ec*(d/(3d0*r**(4d0/3d0)*(1d0 + d*r**(-1d0/3d0))) - 1d0/r + 1d0/ra)
  
          Fc(mu,nu,1) = Fc(mu,nu,1) + weight(iG)*AO(mu,iG)*AO(nu,iG)*dFcdr

        endif

!       Spin-down part contribution

        if(rb > threshold) then

          r = ra + rb

          ec    = ra*rb/(r + d*r**(2d0/3d0))
          dFcdr = ec*(d/(3d0*r**(4d0/3d0)*(1d0 + d*r**(-1d0/3d0))) - 1d0/r + 1d0/rb)
  
          Fc(mu,nu,2) = Fc(mu,nu,2) + weight(iG)*AO(mu,iG)*AO(nu,iG)*dFcdr

        endif

      enddo

    enddo
  enddo

  Fc(:,:,:) = -4d0*a*Fc(:,:,:)

end subroutine UW38_lda_correlation_potential
