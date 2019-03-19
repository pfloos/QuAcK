subroutine SPTDHF(singlet_manifold,triplet_manifold,nBas,nC,nO,nV,nR,nS,ERI,e)

! Perform random phase approximation calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: singlet_manifold,triplet_manifold
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas),e(nBas)

! Local variables

  logical                       :: dRPA,TDA,BSE
  integer                       :: ispin
  double precision,allocatable  :: Omega(:,:),XpY(:,:,:)

  double precision              :: rho
  double precision              :: EcRPA

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|  Time-dependent Hartree-Fock calculation     |'
  write(*,*)'************************************************'
  write(*,*)

! Switch on exchange for TDHF

  dRPA = .false.
 
! Switch off Tamm-Dancoff approximation for TDHF

  TDA = .false.
 
! Switch off Bethe-Salpeter equation for TDHF

  BSE = .false. 

! Memory allocation

  allocate(Omega(nS,nspin),XpY(nS,nS,nspin))

! Singlet manifold

  if(singlet_manifold) then 

    ispin = 1

    call SP_linear_response(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,e,ERI, & 
                         rho,EcRPA,Omega(:,ispin),XpY(:,:,ispin))
    call print_excitation('TDHF ',ispin,nS,Omega(:,ispin))

  endif

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A27,F15.6)') 'RPA correlation energy    =',EcRPA
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)


! Triplet manifold 

  if(triplet_manifold) then 

    ispin = 2

    call SP_linear_response(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,e,ERI, &
                         rho,EcRPA,Omega(:,ispin),XpY(:,:,ispin))
    call print_excitation('TDHF ',ispin,nS,Omega(:,ispin))

  endif

end subroutine SPTDHF
