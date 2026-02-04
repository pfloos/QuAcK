subroutine print_cUHF(nBas,nO,eHF,c,ENuc,ET,EV,EJ,Ex,EW,EUHF)

! Print one- and two-electron energies and other stuff for UHF calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nO(nspin)
  complex*16,intent(in)              :: eHF(nBas,nspin)
  complex*16,intent(in)              :: c(nBas,nBas,nspin)
  double precision,intent(in)        :: ENuc
  complex*16,intent(in)              :: ET(nspin)
  complex*16,intent(in)              :: EV(nspin)
  complex*16,intent(in)              :: EJ(nsp)
  complex*16,intent(in)              :: Ex(nspin)
  complex*16,intent(in)              :: EW(nspin)
  complex*16,intent(in)              :: EUHF

! Local variables

  logical                            :: dump_orb = .false.

! Dump results

  write(*,*)
  write(*,'(A79)')              '---------------------------------------------------------------------------------'
  write(*,'(A33)')              ' Summary                  '
  write(*,'(A79)')              '---------------------------------------------------------------------------------'
  write(*,'(A40,1X,F16.10,1X,A1,F16.10,A1,A3)') ' One-electron    energy = ',real(sum(ET + EV +EW)),'+',aimag(sum(ET+EV+EW)),'i',' au'
  write(*,'(A40,1X,F16.10,1X,A1,F16.10,A1,A3)') ' One-electron  a energy = ',real(ET(1) + EV(1) +EW(1)),'+',aimag(ET(1)+EV(1)+EW(1)),'i',' au'
  write(*,'(A40,1X,F16.10,1X,A1,F16.10,A1,A3)') ' One-electron  b energy = ',real(ET(2) + EV(2) +EW(2)),'+',aimag(ET(2)+EV(2)+EW(2)),'i',' au'
  write(*,'(A40,1X,F16.10,1X,A1,F16.10,A1,A3)') ' Kinetic         energy = ',real(sum(ET)),'+',aimag(sum(ET)),'i',' au'
  write(*,'(A40,1X,F16.10,1X,A1,F16.10,A1,A3)') ' Kinetic       a energy = ',real(ET(1)),'+',aimag(ET(1)),'i',' au'
  write(*,'(A40,1X,F16.10,1X,A1,F16.10,A1,A3)') ' Kinetic       b energy = ',real(ET(2)),'+',aimag(ET(2)),'i',' au'
  write(*,'(A40,1X,F16.10,1X,A1,F16.10,A1,A3)') ' Potential       energy = ',real(sum(EV)),'+',aimag(sum(EV)),'i',' au'
  write(*,'(A40,1X,F16.10,1X,A1,F16.10,A1,A3)') ' Potential     a energy = ',real(EV(1)),'+',aimag(EV(1)),'i',' au'
  write(*,'(A40,1X,F16.10,1X,A1,F16.10,A1,A3)') ' Potential     b energy = ',real(EV(2)),'+',aimag(EV(2)),'i',' au'
  write(*,'(A40,1X,F16.10,1X,A1,F16.10,A1,A3)') ' CAP             energy = ',real(sum(EW)),'+',aimag(sum(EW)),'i',' au'
  write(*,'(A40,1X,F16.10,1X,A1,F16.10,A1,A3)') ' CAP           a energy = ',real(EW(1)),'+',aimag(EW(1)),'i',' au'
  write(*,'(A40,1X,F16.10,1X,A1,F16.10,A1,A3)') ' CAP           b energy = ',real(EW(2)),'+',aimag(EW(2)),'i',' au'
  write(*,'(A79)')              '---------------------------------------------------------------------------------'
  write(*,'(A40,1X,F16.10,1X,A1,F16.10,A1,A3)') ' Two-electron    energy = ',real(sum(EJ) + sum(Ex)),'+',aimag(sum(EJ) + sum(Ex)),'i',' au'
  write(*,'(A40,1X,F16.10,1X,A1,F16.10,A1,A3)') ' Two-electron aa energy = ',real(EJ(1) + Ex(1)),'+',aimag(EJ(1)+Ex(1)),'i',' au'
  write(*,'(A40,1X,F16.10,1X,A1,F16.10,A1,A3)') ' Two-electron ab energy = ',real(EJ(2)),'+',aimag(EJ(2)),'i',' au'
  write(*,'(A40,1X,F16.10,1X,A1,F16.10,A1,A3)') ' Two-electron bb energy = ',real(EJ(2) + Ex(2)),'+',aimag(EJ(2)+Ex(2)),'i',' au'
  write(*,'(A40,1X,F16.10,1X,A1,F16.10,A1,A3)') ' Hartree         energy = ',real(sum(EJ)),'+',aimag(sum(EJ)),'i',' au'
  write(*,'(A40,1X,F16.10,1X,A1,F16.10,A1,A3)') ' Hartree      aa energy = ',real(EJ(1)),'+',aimag(EJ(1)),'i',' au'
  write(*,'(A40,1X,F16.10,1X,A1,F16.10,A1,A3)') ' Hartree      ab energy = ',real(EJ(2)),'+',aimag(EJ(2)),'i',' au'
  write(*,'(A40,1X,F16.10,1X,A1,F16.10,A1,A3)') ' Hartree      bb energy = ',real(EJ(3)),'+',aimag(EJ(3)),'i',' au'
  write(*,'(A40,1X,F16.10,1X,A1,F16.10,A1,A3)') ' Exchange        energy = ',real(sum(Ex)),'+',aimag(sum(Ex)),'i',' au'
  write(*,'(A40,1X,F16.10,1X,A1,F16.10,A1,A3)') ' Exchange      a energy = ',real(Ex(1)),'+',aimag(Ex(1)),'i',' au'
  write(*,'(A40,1X,F16.10,1X,A1,F16.10,A1,A3)') ' Exchange      b energy = ',real(Ex(2)),'+',aimag(Ex(2)),'i',' au'
  write(*,'(A79)')              '---------------------------------------------------------------------------------'
  write(*,'(A40,1X,F16.10,1X,A1,F16.10,A1,A3)') ' Electronic   energy = ',real(EUHF),'+',aimag(EUHF),'i',' au'
  write(*,'(A40,1X,F16.10,19X,A3)') ' Nuclear   repulsion = ',ENuc,' au'
  write(*,'(A40,1X,F16.10,1X,A1,F16.10,A1,A3)') ' cUHF         energy = ',real(EUHF + ENuc),'+',aimag(EUHF+ENuc),'i',' au'
  write(*,'(A79)')              '---------------------------------------------------------------------------------'
  write(*,*)

! Print results

  if(dump_orb) then
    write(*,'(A69)') '-----------------------------------------'
    write(*,'(A69)') 'UHF spin-up   orbital coefficients '
    write(*,'(A69)') '-----------------------------------------'
    call complex_matout(nBas,nBas,c(:,:,1))
    write(*,*)
    write(*,'(A69)') '-----------------------------------------'
    write(*,'(A69)') 'UHF spin-down orbital coefficients '
    write(*,'(A69)') '-----------------------------------------'
    call complex_matout(nBas,nBas,c(:,:,2))
    write(*,*)
  end if

  write(*,'(A79)') '---------------------------------------------------------------------------------'
  write(*,'(A33)') ' UHF spin-up   orbital energies  '
  write(*,'(A79)') '---------------------------------------------------------------------------------'
  call complex_vecout(nBas,eHF(:,1))
  write(*,*)
  write(*,'(A79)') '---------------------------------------------------------------------------------'
  write(*,'(A33)') ' UHF spin-down orbital energies  '
  write(*,'(A79)') '---------------------------------------------------------------------------------'
  call complex_vecout(nBas,eHF(:,2))
  write(*,*)

end subroutine 
