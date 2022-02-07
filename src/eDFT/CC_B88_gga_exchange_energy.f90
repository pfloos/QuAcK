subroutine CC_B88_gga_exchange_energy(nEns,wEns,nCC,aCC,nGrid,weight,&
                                     rho,drho,Cx_choice,Ex)

! Compute the unrestricted version of the curvature-corrected exchange functional

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nCC
  double precision,intent(in)   :: aCC(nCC,nEns-1)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: drho(3,nGrid)
  integer,intent(in)            :: Cx_choice

! Local variables

  integer                       :: iG
  double precision              :: b
  double precision              :: r,g,x

  double precision              :: a1,b1,c1,d1,w1
  double precision              :: a2,b2,c2,d2,w2
  double precision              :: Fx1,Fx2,Cx

! Output variables

  double precision              :: Ex

! Coefficients for B88 GGA exchange functional

  b = 0.0042d0

! Defining enhancements factor for weight-dependent functionals

! Parameters for first state

     a1 = aCC(1,1)
     b1 = aCC(2,1)
     c1 = aCC(3,1)
     d1 = aCC(4,1)

! Parameters for second state

     a2 = aCC(1,2)
     b2 = aCC(2,2)
     c2 = aCC(3,2)
     d2 = aCC(4,2)


     w1 = wEns(2)
     Fx1 = 1d0 + a1*w1 + b1*w1**2 + c1*w1**3 + d1*w1**4

     w2 = wEns(3)
     Fx2 = 1d0 + a2*w2 + b2*w2**2 + c2*w2**3 + d2*w2**4


  select case (Cx_choice)

    case(1)
      Cx = Fx1

    case(2)
      Cx = Fx2

    case(3)
      Cx = Fx2*Fx1

    case default
      Cx = 1.d0

  end select


! Compute GIC-GGA exchange energy

  Ex = 0d0

  do iG=1,nGrid

    r = max(0d0,rho(iG))

    if(r > threshold) then
      g = drho(1,iG)**2 + drho(2,iG)**2 + drho(3,iG)**2
      x = sqrt(g)/r**(4d0/3d0)

      Ex = Ex + weight(iG)*r**(4d0/3d0)*(CxLSDA - b*x**2/(1d0 + 6d0*b*x*asinh(x)))

    end if

  end do

  Ex = Cx*Ex

end subroutine CC_B88_gga_exchange_energy
