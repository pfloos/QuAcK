subroutine unrestricted_correlation_energy(rung,DFA,nEns,wEns,nGrid,weight,rho,drho,Ec)

! Compute the unrestricted version of the correlation energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: rung
  character(len=12),intent(in)  :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid,nspin)
  double precision,intent(in)   :: drho(ncart,nGrid,nspin)

! Local variables

! Output variables

  double precision,intent(out)  :: Ec(nsp)

  select case (rung)

!   Hartree calculation

    case(0) 

      Ec(:) = 0d0

!   LDA functionals

    case(1) 

      call unrestricted_lda_correlation_energy(DFA,nEns,wEns,nGrid,weight,rho,Ec)

!   GGA functionals

    case(2) 

      call unrestricted_gga_correlation_energy(DFA,nEns,wEns,nGrid,weight,rho,drho,Ec)

!   MGGA functionals

    case(3) 

      call unrestricted_mgga_correlation_energy(DFA,nEns,wEns,nGrid,weight,rho,drho,Ec)

!   Hybrid functionals

    case(4) 

      call unrestricted_hybrid_correlation_energy(DFA,nEns,wEns,nGrid,weight,rho,drho,Ec)

  end select
 
end subroutine unrestricted_correlation_energy
