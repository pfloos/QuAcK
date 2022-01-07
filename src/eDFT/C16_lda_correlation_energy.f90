subroutine C16_lda_correlation_energy(nGrid,weight,rho,Ec)

! Compute unrestricted Chachiyo's LDA correlation energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid,nspin)

! Local variables

  integer                       :: iG
  double precision              :: ra,rb,r,rs
  double precision              :: a_p,b_p,ec_p
  double precision              :: a_f,b_f,ec_f
  double precision              :: z,fz,ec_z

! Output variables

  double precision              :: Ec(nsp)

! Coefficients for Chachiyo's LDA correlation

  a_p = (log(2d0) - 1d0)/(2d0*pi**2)
  b_p = 20.4562557d0

  a_f = (log(2d0) - 1d0)/(4d0*pi**2)
  b_f = 27.4203609d0

! Compute LDA correlation energy

  Ec(:) = 0d0

  do iG=1,nGrid

!   Spin-up and spin-down densities

    ra = max(0d0,rho(iG,1))
    rb = max(0d0,rho(iG,2))

!   Total density

    r = ra + rb

!   Spin-up part contribution

    if(ra > threshold) then

      rs = (4d0*pi*ra/3d0)**(-1d0/3d0)
      ec_f = a_f*log(1d0 + b_f/rs + b_f/rs**2)

      Ec(1) = Ec(1) + weight(iG)*ec_f*ra

    endif

!   Opposite-spin contribution

    if(r > threshold) then

      rs  = (4d0*pi*r/3d0)**(-1d0/3d0)
      ec_p = a_p*log(1d0 + b_p/rs + b_p/rs**2)
      ec_f = a_f*log(1d0 + b_f/rs + b_f/rs**2)

      z = (ra-rb)/r
      fz = (1d0 + z)**(4d0/3d0) + (1d0 - z)**(4d0/3d0) - 2d0
      fz = fz/(2d0*(2d0**(1d0/3d0) - 1d0))

      ec_z = ec_p + (ec_f - ec_p)*fz

      Ec(2) = Ec(2) + weight(iG)*ec_z*r

      endif

!   Spin-down contribution

    if(rb > threshold) then

      rs = (4d0*pi*rb/3d0)**(-1d0/3d0)
      ec_f = a_f*log(1d0 + b_f/rs + b_f/rs**2)

      Ec(3) = Ec(3) + weight(iG)*ec_f*rb

    endif

  enddo

  Ec(2) = Ec(2) - Ec(1) - Ec(3)

end subroutine C16_lda_correlation_energy
