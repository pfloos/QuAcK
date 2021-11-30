subroutine unrestricted_correlation_derivative_discontinuity(rung,DFA,nEns,wEns,nGrid,weight,rhow,drhow,kappa,Ec)

! Compute the correlation part of the derivative discontinuity

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: rung
  integer,intent(in)            :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid,nspin)
  double precision,intent(in)   :: drhow(ncart,nGrid,nspin)
  double precision,intent(in)   :: kappa(nEns)

! Local variables

! Output variables

  double precision,intent(out)  :: Ec(nsp,nEns)

  select case (rung)

!   Hartree calculation

    case(0) 

      Ec(:,:) = 0d0

!   LDA functionals

    case(1) 

      call unrestricted_lda_correlation_derivative_discontinuity(DFA,nEns,wEns,nGrid,weight,rhow,Ec)

!   GGA functionals

    case(2) 

      call unrestricted_gga_correlation_derivative_discontinuity(DFA,nEns,wEns,nGrid,weight,rhow,Ec)

!   MGGA functionals

    case(3) 

      call unrestricted_mgga_correlation_derivative_discontinuity(DFA,nEns,wEns,nGrid,weight,rhow,Ec)

!   Hybrid functionals

    case(4) 

      call unrestricted_hybrid_correlation_derivative_discontinuity(DFA,nEns,wEns,nGrid,weight,rhow,Ec)

  end select
 
end subroutine unrestricted_correlation_derivative_discontinuity
