subroutine restricted_lda_correlation_energy(DFA,nEns,wEns,nGrid,weight,rho,Ec)

! Select LDA correlation functional

  implicit none
  include 'parameters.h'

! Input variables

  character(len=12),intent(in)  :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)

! Output variables

  double precision,intent(out)  :: Ec

! Select correlation functional

  select case (DFA)

!   Hartree-Fock

    case ('HF')

      Ec = 0d0

    case ('RVWN5')

      call RVWN5_lda_correlation_energy(nGrid,weight(:),rho(:),Ec)

    case ('RMFL20')

      call RMFL20_lda_correlation_energy(nEns,wEns(:),nGrid,weight(:),rho(:),Ec)

    case default

      call print_warning('!!! LDA correlation functional not available !!!')
      stop

  end select

end subroutine restricted_lda_correlation_energy
