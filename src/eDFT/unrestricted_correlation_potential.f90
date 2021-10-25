subroutine unrestricted_correlation_potential(rung,DFA,nEns,wEns,nGrid,weight,nBas,AO,dAO,rho,drho,Fc)

! Compute the correlation potential 

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: rung
  integer,intent(in)            :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: dAO(ncart,nBas,nGrid)
  double precision,intent(in)   :: rho(nGrid,nspin)
  double precision,intent(in)   :: drho(ncart,nGrid,nspin)

! Local variables

  double precision,allocatable  :: FcLDA(:,:,:)
  double precision,allocatable  :: FcGGA(:,:,:)
  double precision              :: aC

! Output variables

  double precision,intent(out)  :: Fc(nBas,nBas,nspin)

! Memory allocation

  select case (rung)

!   Hartree calculation

    case(0)

      Fc(:,:,:) = 0d0

!   LDA functionals

    case(1)

      call unrestricted_lda_correlation_potential(DFA,nEns,wEns,nGrid,weight,nBas,AO,rho,Fc)

!   GGA functionals

    case(2)

      call unrestricted_gga_correlation_potential(DFA,nEns,wEns,nGrid,weight,nBas,AO,dAO,rho,drho,Fc)

!   MGGA functionals

    case(3)

      call unrestricted_mgga_correlation_potential(DFA,nEns,wEns,nGrid,weight,nBas,AO,dAO,rho,drho,Fc)

!   Hybrid functionals

    case(4)

      call unrestricted_hybrid_correlation_potential(DFA,nEns,wEns,nGrid,weight,nBas,AO,dAO,rho,drho,Fc)

  end select

end subroutine unrestricted_correlation_potential
