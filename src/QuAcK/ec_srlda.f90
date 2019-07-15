subroutine ec_srlda(nGrid,weight,rho,mu)

! Compute sr-lda ec

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: mu(nGrid)

! Local variables

  integer                       :: iG
  double precision              :: r,ra,rb
  double precision              :: rs
  double precision              :: ec_lda,ecmd_lda
  double precision              :: ec,ecmd,vcup,vcdw

! Initialization

  ec_lda   = 0d0
  ecmd_lda = 0d0

  do iG=1,ngrid

    r = max(0d0,rho(iG))
    ra = 0.5d0*r
    rb = 0.5d0*r

    if(r > threshold) then

      rs = (4d0*pi*r/3d0)**(-1d0/3d0)

      call lsdsr(rs,0d0,mu(iG),ec,vcup,vcdw)
      call ESRC_MD_LDAERF(mu(iG),ra,rb,.true.,ecmd)
      ec_lda   = ec_lda   + weight(iG)*ec*r
      ecmd_lda = ecmd_lda + weight(iG)*ecmd*r

    end if

  end do

  write(*,*) 
  write(*,'(A32,1X,F16.10)') 'Ec(sr-LDA)   = ',ec_lda
  write(*,'(A32,1X,F16.10)') 'Ecmd(sr-LDA) = ',ec_lda + ecmd_lda
  write(*,*) 

end subroutine ec_srlda
