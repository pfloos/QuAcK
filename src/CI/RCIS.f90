subroutine RCIS(dotest,singlet,triplet,doCIS_D,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eHF)

! Perform configuration interaction single calculation`

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

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
  write(*,*)'******************************'
  write(*,*)'* Restricted CIS Calculation *'
  write(*,*)'******************************'
  write(*,*)

! Adiabatic connection scaling

  lambda = 1d0

! Memory allocation

  allocate(A(nS,nS), Om(nS))

! Compute CIS matrix

  if(singlet) then

    ispin = 1
    call phRLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,lambda,eHF,ERI,A)
 
    if(dump_matrix) then
      print*,'CIS matrix (singlet state)'
      call matout(nS,nS,A)
      write(*,*)
    end if

    call diagonalize_matrix(nS,A,Om)
    call print_excitation_energies('CIS@RHF','singlet',nS,Om)
    call phLR_transition_vectors(.true.,nBas,nC,nO,nV,nR,nS,dipole_int,Om,transpose(A),transpose(A))
 
    if(dump_trans) then
      print*,'Singlet CIS transition vectors'
      call matout(nS,nS,A)
      write(*,*)
    end if

    ! Compute CIS(D) correction 

    maxS = min(maxS,nS)
    if(doCIS_D) call CIS_D(ispin,nBas,nC,nO,nV,nR,nS,maxS,eHF,ERI,Om(1:maxS),A(:,1:maxS))

    ! Testing zone
  
    if(dotest) then

      call dump_test_value('R','CIS singlet excitation energy',Om(1))

    end if

  end if

  if(triplet) then

    ispin = 2
    call phRLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,lambda,eHF,ERI,A)
 
    if(dump_matrix) then
      print*,'CIS matrix (triplet state)'
      call matout(nS,nS,A)
      write(*,*)
    end if
 
    call diagonalize_matrix(nS,A,Om)
    call print_excitation_energies('CIS@RHF','triplet',nS,Om)
    call phLR_transition_vectors(.false.,nBas,nC,nO,nV,nR,nS,dipole_int,Om,transpose(A),transpose(A))

    if(dump_trans) then
      print*,'Triplet CIS transition vectors'
      call matout(nS,nS,A)
      write(*,*)
    end if

    ! Compute CIS(D) correction 

    maxS = min(maxS,nS)
    if(doCIS_D) call CIS_D(ispin,nBas,nC,nO,nV,nR,nS,maxS,eHF,ERI,Om(1:maxS),A(:,1:maxS))

    ! Testing zone
 
    if(dotest) then
 
      call dump_test_value('R','CIS triplet excitation energy',Om(1))
 
    end if

  end if

  deallocate(A,Om)

end subroutine 
