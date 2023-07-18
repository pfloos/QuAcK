subroutine print_RHF(nBas,nO,eHF,cHF,ENuc,ET,EV,EJ,EK,ERHF,dipole)

! Print one-electron energies and other stuff for G0W0

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nO
  double precision,intent(in)        :: eHF(nBas)
  double precision,intent(in)        :: cHF(nBas,nBas)
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ET
  double precision,intent(in)        :: EV
  double precision,intent(in)        :: EJ
  double precision,intent(in)        :: EK
  double precision,intent(in)        :: ERHF
  double precision,intent(in)        :: dipole(ncart)

! Local variables

  integer                            :: ixyz
  integer                            :: HOMO
  integer                            :: LUMO
  double precision                   :: Gap

! HOMO and LUMO

  HOMO = nO
  LUMO = HOMO + 1
  Gap = eHF(LUMO)-eHF(HOMO)

! Dump results


  write(*,*)
  write(*,'(A50)')           '-----------------------------------------'
  write(*,'(A32)')           ' Summary              '
  write(*,'(A50)')           '-----------------------------------------'
  write(*,'(A32,1X,F16.10,A3)') ' One-electron energy: ',ET + EV,' au'
  write(*,'(A32,1X,F16.10,A3)') ' Kinetic      energy: ',ET,' au'
  write(*,'(A32,1X,F16.10,A3)') ' Potential    energy: ',EV,' au'
  write(*,'(A50)')           '-----------------------------------------'
  write(*,'(A32,1X,F16.10,A3)') ' Two-electron energy: ',EJ + EK,' au'
  write(*,'(A32,1X,F16.10,A3)') ' Hartree      energy: ',EJ,' au'
  write(*,'(A32,1X,F16.10,A3)') ' Exchange     energy: ',EK,' au'
  write(*,'(A50)')           '-----------------------------------------'
  write(*,'(A32,1X,F16.10,A3)') ' Electronic   energy: ',ERHF,' au'
  write(*,'(A32,1X,F16.10,A3)') ' Nuclear   repulsion: ',ENuc,' au'
  write(*,'(A32,1X,F16.10,A3)') ' RHF          energy: ',ERHF + ENuc,' au'
  write(*,'(A50)')           '-----------------------------------------'
  write(*,'(A32,1X,F16.6,A3)')  ' HF HOMO      energy: ',eHF(HOMO)*HaToeV,' eV'
  write(*,'(A32,1X,F16.6,A3)')  ' HF LUMO      energy: ',eHF(LUMO)*HaToeV,' eV'
  write(*,'(A32,1X,F16.6,A3)')  ' HF HOMO-LUMO gap   : ',Gap*HaToeV,' eV'
  write(*,'(A50)')           '-----------------------------------------'
  write(*,'(A35)')           ' Dipole moment (Debye)    '
  write(*,'(10X,4A10)')      'X','Y','Z','Tot.'
  write(*,'(10X,4F10.6)')    (dipole(ixyz)*auToD,ixyz=1,ncart),norm2(dipole)*auToD
  write(*,'(A50)')           '-----------------------------------------'
  write(*,*)

! Print results

  write(*,'(A50)')  '---------------------------------------'
  write(*,'(A32)') 'MO coefficients'
  write(*,'(A50)')  '---------------------------------------'
  call matout(nBas,nBas,cHF)
  write(*,*)
  write(*,'(A50)') '---------------------------------------'
  write(*,'(A32)') 'MO energies'
  write(*,'(A50)') '---------------------------------------'
  call matout(nBas,1,eHF)
  write(*,*)

end subroutine 
