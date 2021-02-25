subroutine unrestricted_hybrid_correlation_potential(DFA,nEns,wEns,nGrid,weight,nBas,AO,dAO,rho,drho,Fc)

! Compute the correlation potential for hybrid functionals

  implicit none
  include 'parameters.h'

! Input variables

  character(len=12),intent(in)  :: DFA
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

  select case (DFA)

    case('HF')

      Fc(:,:,:) = 0d0

    case('B3LYP')

      allocate(FcLDA(nBas,nBas,nspin),FcGGA(nBas,nBas,nspin))

      aC = 0.81d0

      call unrestricted_lda_correlation_potential('VWN3        ',nEns,wEns,nGrid,weight,nBas,AO,rho,FcLDA)
      call unrestricted_gga_correlation_potential('LYP         ',nEns,wEns,nGrid,weight,nBas,AO,dAO,rho,drho,FcGGA)

      Fc(:,:,:) = FcLDA(:,:,:) + aC*(FcGGA(:,:,:) - FcLDA(:,:,:))

    case('BHHLYP')

      allocate(FcGGA(nBas,nBas,nspin))

      call unrestricted_gga_correlation_potential('LYP         ',nEns,wEns,nGrid,weight,nBas,AO,dAO,rho,drho,Fc)

    case('PBE')

      allocate(FcGGA(nBas,nBas,nspin))

      call unrestricted_gga_correlation_potential('PBE         ',nEns,wEns,nGrid,weight,nBas,AO,dAO,rho,drho,Fc)

    case default

      call print_warning('!!! Hybrid correlation potential not available !!!')
      stop

  end select

end subroutine unrestricted_hybrid_correlation_potential
