subroutine US51_lda_exchange_individual_energy(nEns,nGrid,weight,rhow,rho,LZx,Ex)

! Compute the restricted version of Slater's LDA exchange individual energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nEns
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid,nspin)
  double precision,intent(in)   :: rho(nGrid,nspin,nEns)

! Local variables

  integer                       :: iG
  integer                       :: iEns
  integer                       :: ispin
  double precision              :: r
  double precision              :: rI
  double precision              :: e
  double precision              :: dedr

! Output variables

  double precision,intent(out)  :: LZx(nspin)
  double precision,intent(out)  :: Ex(nspin,nEns)

  LZx(:)  = 0d0
  Ex(:,:) = 0d0

  do ispin=1,nspin

    do iG=1,nGrid
 
      r = max(0d0,rhow(iG,ispin))

      if(r > threshold) then
 
        e    =         CxLSDA*r**(+1d0/3d0)
        dedr = 1d0/3d0*CxLSDA*r**(-2d0/3d0)
 
        LZx(ispin) = LZx(ispin) - weight(iG)*dedr*r*r      
 
        do iEns=1,nEns
 
          rI = max(0d0,rho(iG,ispin,iEns))
 
          if(rI > threshold) Ex(ispin,iEns) = Ex(ispin,iEns) + weight(iG)*(e+dedr*r)*rI 
 
        end do
 
      endif

    enddo

  enddo

end subroutine US51_lda_exchange_individual_energy
