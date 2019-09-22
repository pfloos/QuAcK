subroutine print_RHF(nBas,nO,eHF,cHF,ENuc,ET,EV,EJ,EK,ERHF)

! Print one-electron energies and other stuff for G0W0

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas,nO
  double precision,intent(in)        :: eHF(nBas),cHF(nBas,nBas),ENuc,ET,EV,EJ,EK,ERHF

  integer                            :: HOMO,LUMO
  double precision                   :: Gap

! HOMO and LUMO

  HOMO = nO
  LUMO = HOMO + 1
  Gap = eHF(LUMO)-eHF(HOMO)

! Dump results


  write(*,*)
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A32)')           ' Summary              '
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A32,1X,F16.10)') ' One-electron energy  ',ET + EV
  write(*,'(A32,1X,F16.10)') ' Kinetic      energy  ',ET
  write(*,'(A32,1X,F16.10)') ' Potential    energy  ',EV
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A32,1X,F16.10)') ' Two-electron energy  ',EJ + EK
  write(*,'(A32,1X,F16.10)') ' Coulomb      energy  ',EJ
  write(*,'(A32,1X,F16.10)') ' Exchange     energy  ',EK
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A32,1X,F16.10)') ' Electronic   energy  ',ERHF
  write(*,'(A32,1X,F16.10)') ' Nuclear   repulsion  ',ENuc
  write(*,'(A32,1X,F16.10)') ' Hartree-Fock energy  ',ERHF + ENuc
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A36,F13.6)')  ' HF HOMO      energy (eV):',eHF(HOMO)*HaToeV
  write(*,'(A36,F13.6)')  ' HF LUMO      energy (eV):',eHF(LUMO)*HaToeV
  write(*,'(A36,F13.6)')  ' HF HOMO-LUMO gap    (eV):',Gap*HaToeV
  write(*,'(A50)')           '---------------------------------------'
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

end subroutine print_RHF


