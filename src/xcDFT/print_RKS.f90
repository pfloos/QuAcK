subroutine print_RKS(nBas,nO,e,C,ENuc,ET,EV,EJ,Ex,Ec,EKS)

! Print one- and two-electron energies and other stuff for RKS calculation

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas,nO
  double precision,intent(in)        :: e(nBas),c(nBas,nBas),ENuc,ET,EV,EJ,Ex,Ec,EKS

  integer                            :: HOMO,LUMO
  double precision                   :: Gap

! HOMO and LUMO

  HOMO = nO
  LUMO = HOMO + 1
  Gap = e(LUMO) - e(HOMO)

! Dump results


  write(*,*)
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A32)')           ' Summary              '
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A32,1X,F16.10)') ' One-electron energy  ',ET + EV
  write(*,'(A32,1X,F16.10)') ' Kinetic      energy  ',ET
  write(*,'(A32,1X,F16.10)') ' Potential    energy  ',EV
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A32,1X,F16.10)') ' Two-electron energy  ',EJ + Ex + Ec
  write(*,'(A32,1X,F16.10)') ' Coulomb      energy  ',EJ
  write(*,'(A32,1X,F16.10)') ' Exchange     energy  ',Ex
  write(*,'(A32,1X,F16.10)') ' Correlation  energy  ',Ec
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A32,1X,F16.10)') ' Electronic   energy  ',EKS
  write(*,'(A32,1X,F16.10)') ' Nuclear   repulsion  ',ENuc
  write(*,'(A32,1X,F16.10)') ' Kohn-Sham    energy  ',EKS + ENuc
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A36,F13.6)')     ' KS HOMO      energy (eV):',e(HOMO)*HatoeV
  write(*,'(A36,F13.6)')     ' KS LUMO      energy (eV):',e(LUMO)*Hatoev
  write(*,'(A36,F13.6)')     ' KS HOMO-LUMO gap    (eV):',Gap*Hatoev
  write(*,'(A50)')           '---------------------------------------'
  write(*,*)

! Print results

  write(*,'(A50)') '---------------------------------------'
  write(*,'(A50)') 'Kohn-Sham orbital coefficients         '
  write(*,'(A50)') '---------------------------------------'
  call matout(nBas,nBas,C)
  write(*,*)
  write(*,'(A50)') '---------------------------------------'
  write(*,'(A50)') ' Kohn-Sham orbital energies            '
  write(*,'(A50)') '---------------------------------------'
  call matout(nBas,1,e)
  write(*,*)

end subroutine print_RKS


