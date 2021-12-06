subroutine unrestricted_lda_exchange_individual_energy(DFA,LDA_centered,nEns,wEns,nCC,aCC,nGrid,weight,rhow,&
                                                       rho,Cx_choice,doNcentered,LZx,Ex)

! Compute LDA exchange energy for individual states

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: LDA_centered
  integer,intent(in)            :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nCC
  double precision,intent(in)   :: aCC(nCC,nEns-1)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid,nspin)
  double precision,intent(in)   :: rho(nGrid,nspin,nEns)
  integer,intent(in)            :: Cx_choice
  logical,intent(in)            :: doNcentered

! Output variables

  double precision              :: LZx(nspin)
  double precision              :: Ex(nspin,nEns)

! Select correlation functional

  select case (DFA)

    case (1)

      call US51_lda_exchange_individual_energy(nEns,nGrid,weight,rhow,rho,LZx,Ex)

    case (2)

      call UCC_lda_exchange_individual_energy(nEns,wEns,nCC,aCC,nGrid,weight,rhow,rho, &
                                              Cx_choice,doNcentered,LZx,Ex)

    case default

      call print_warning('!!! LDA exchange individual energy not available !!!')
      stop

  end select

end subroutine unrestricted_lda_exchange_individual_energy
