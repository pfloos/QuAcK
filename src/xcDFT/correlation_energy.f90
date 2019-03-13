subroutine correlation_energy(rung,DFA,nEns,wEns,nGrid,weight,rho,drho,Ec)

! Compute the correlation energy

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
  double precision,intent(in)   :: drho(3,nGrid,nspin)

! Local variables

  double precision              :: EcLDA(nsp)
  double precision              :: EcGGA(nsp)
  double precision              :: aC

! Output variables

  double precision,intent(out)  :: Ec(nsp)

  select case (rung)

!   Hartree calculation

    case(0) 

      Ec(:) = 0d0

!   LDA functionals

    case(1) 

      call lda_correlation_energy(DFA,nEns,wEns(:),nGrid,weight(:),rho(:,:),Ec(:))

!   GGA functionals

    case(2) 

      call gga_correlation_energy(DFA,nEns,wEns(:),nGrid,weight(:),rho(:,:),drho(:,:,:),Ec(:))

!   Hybrid functionals

    case(4) 

      aC = 0.81d0

      call lda_correlation_energy(DFA,nEns,wEns(:),nGrid,weight(:),rho(:,:),EcLDA(:))
      call gga_correlation_energy(DFA,nEns,wEns(:),nGrid,weight(:),rho(:,:),drho(:,:,:),EcGGA(:))

      Ec(:) = EcLDA(:) + aC*(EcGGA(:) - EcLDA(:)) 

!   Hartree-Fock calculation

    case(666) 

      Ec(:) = 0d0

  end select
 
end subroutine correlation_energy
