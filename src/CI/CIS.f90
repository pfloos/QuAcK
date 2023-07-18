subroutine CIS(singlet,triplet,doCIS_D,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eHF)

! Perform configuration interaction single calculation`

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  logical,intent(in)            :: doCIS_D
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  logical                       :: dump_matrix = .false.
  logical                       :: dump_trans = .false.
  logical                       :: dRPA = .false.
  integer                       :: ispin
  integer                       :: maxS = 10
  double precision              :: lambda
  double precision,allocatable  :: A(:,:),Om(:)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|      Configuration Interaction Singles       |'
  write(*,*)'************************************************'
  write(*,*)

! Adiabatic connection scaling

  lambda = 1d0

! Memory allocation

  allocate(A(nS,nS),Om(nS))

! Compute CIS matrix

  if(singlet) then

    ispin = 1
    call phLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,lambda,eHF,ERI,A)
 
    if(dump_matrix) then
      print*,'CIS matrix (singlet state)'
      call matout(nS,nS,A)
      write(*,*)
    endif

    call diagonalize_matrix(nS,A,Om)
    call print_excitation('CIS         ',ispin,nS,Om)
    call print_transition_vectors_ph(.true.,nBas,nC,nO,nV,nR,nS,dipole_int,Om,transpose(A),transpose(A))
 
    if(dump_trans) then
      print*,'Singlet CIS transition vectors'
      call matout(nS,nS,A)
      write(*,*)
    endif

    ! Compute CIS(D) correction 

    maxS = min(maxS,nS)
    if(doCIS_D) call D_correction(ispin,nBas,nC,nO,nV,nR,nS,maxS,eHF,ERI,Om(1:maxS),A(:,1:maxS))

  endif

  if(triplet) then

    ispin = 2
    call phLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,lambda,eHF,ERI,A)
 
    if(dump_matrix) then
      print*,'CIS matrix (triplet state)'
      call matout(nS,nS,A)
      write(*,*)
    endif
 
    call diagonalize_matrix(nS,A,Om)
    call print_excitation('CIS          ',ispin,nS,Om)
    call print_transition_vectors_ph(.false.,nBas,nC,nO,nV,nR,nS,dipole_int,Om,transpose(A),transpose(A))

    if(dump_trans) then
      print*,'Triplet CIS transition vectors'
      call matout(nS,nS,A)
      write(*,*)
    endif

    ! Compute CIS(D) correction 

    maxS = min(maxS,nS)
    if(doCIS_D) call D_correction(ispin,nBas,nC,nO,nV,nR,nS,maxS,eHF,ERI,Om(1:maxS),A(:,1:maxS))

  endif

end subroutine 
