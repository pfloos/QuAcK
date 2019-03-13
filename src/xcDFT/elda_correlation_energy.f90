subroutine elda_correlation_energy(nEns,aLF,nGrid,weight,rho,Ec)

! Compute LDA correlation energy of 2-glomium for various states

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nEns
  double precision,intent(in)   :: aLF(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid,nspin)

! Local variables

  integer                       :: iG
  double precision              :: ra,rb,r,ec_p

! Output variables

  double precision,intent(out)  :: Ec(nsp)


! Compute eLDA correlation energy

  Ec(:) = 0d0

  do iG=1,nGrid

    ra = max(0d0,rho(iG,1))
    rb = max(0d0,rho(iG,2))
    r  = ra + rb

!   Spin-up contribution

    if(ra > threshold) then

      ec_p = aLF(1)/(1d0 + aLF(2)*ra**(-1d0/6d0) + aLF(3)*ra**(-1d0/3d0))

      Ec(1) = Ec(1) + weight(iG)*ec_p*ra

    end if

!   Opposite-spin contribution

    if(r > threshold) then

      ec_p = aLF(1)/(1d0 + aLF(2)*r**(-1d0/6d0) + aLF(3)*r**(-1d0/3d0))

      Ec(2) = Ec(2) + weight(iG)*ec_p*r

    end if

!   Spin-down contribution

    if(rb > threshold) then

      ec_p = aLF(1)/(1d0 + aLF(2)*rb**(-1d0/6d0) + aLF(3)*rb**(-1d0/3d0))

      Ec(3) = Ec(3) + weight(iG)*ec_p*rb

    end if

  end do

  Ec(2) = Ec(2) - Ec(1) - Ec(3)

end subroutine elda_correlation_energy
