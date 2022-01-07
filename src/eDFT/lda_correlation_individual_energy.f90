subroutine lda_correlation_individual_energy(DFA,LDA_centered,nEns,wEns,nGrid,weight,rhow,rho,LZc,Ec)

! Compute LDA correlation energy for individual states

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: LDA_centered
  integer,intent(in)            :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid,nspin)
  double precision,intent(in)   :: rho(nGrid,nspin,nEns)

! Output variables

  double precision              :: LZc(nsp)
  double precision              :: Ec(nsp,nEns)

! Select correlation functional

  select case (DFA)

    case (1)

!     call W38_lda_correlation_individual_energy(nGrid,weight,rhow,rho,LZc,Ec)

    case (2)

!     call PW92_lda_correlation_individual_energy(nGrid,weight,rhow,rho,LZc,Ec)

    case (3)

!     call VWN3_lda_correlation_individual_energy(nEns,nGrid,weight,rhow,rho,LZc,Ec)

    case (4)

      call VWN5_lda_correlation_individual_energy(nEns,nGrid,weight,rhow,rho,LZc,Ec)

    case default

      call print_warning('!!! LDA correlation functional not available !!!')
      stop

  end select

end subroutine lda_correlation_individual_energy
