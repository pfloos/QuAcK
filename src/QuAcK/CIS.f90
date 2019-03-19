subroutine CIS(singlet_manifold,triplet_manifold, & 
               nBas,nC,nO,nV,nR,nS,ERI,eHF)

! Perform configuration interaction single calculation`

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: singlet_manifold,triplet_manifold
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas),eHF(nBas)

! Local variables

  logical                       :: dRPA
  logical                       :: dump_matrix = .false.
  logical                       :: dump_trans = .false.
  integer                       :: ispin
  double precision,allocatable  :: A(:,:),Omega(:)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|      Configuration Interaction Singles       |'
  write(*,*)'************************************************'
  write(*,*)

! Switch on exchange for CIS

  dRPA = .false.

! Memory allocation

  allocate(A(nS,nS),Omega(nS))

! Compute CIS matrix

  if(singlet_manifold) then

    ispin = 1
    call linear_response_A_matrix(ispin,dRPA,nBas,nC,nO,nV,nR,nS,eHF,ERI,A)
 
    if(dump_matrix) then
      print*,'CIS matrix (singlet state)'
      call matout(nS,nS,A)
      write(*,*)
    endif

    call diagonalize_matrix(nS,A,Omega)
    call print_excitation('CIS  ',ispin,nS,Omega)
 
    if(dump_trans) then
      print*,'Singlet CIS transition vectors'
      call matout(nS,nS,A)
      write(*,*)
    endif

  endif

  if(triplet_manifold) then

    ispin = 2
    call linear_response_A_matrix(ispin,dRPA,nBas,nC,nO,nV,nR,nS,eHF,ERI,A)
 
    if(dump_matrix) then
      print*,'CIS matrix (triplet state)'
      call matout(nS,nS,A)
      write(*,*)
    endif
 
    call diagonalize_matrix(nS,A,Omega)
    call print_excitation('CIS  ',ispin,nS,Omega)

    if(dump_trans) then
      print*,'Triplet CIS transition vectors'
      call matout(nS,nS,A)
      write(*,*)
    endif

  endif

end subroutine CIS
