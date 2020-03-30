subroutine print_RKS(nBas,nO,eps,c,ENuc,ET,EV,EJ,Ex,Ec,Ew)

! Print one- and two-electron energies and other stuff for KS calculation

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nO
  double precision,intent(in)        :: eps(nBas)
  double precision,intent(in)        :: c(nBas,nBas)
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ET
  double precision,intent(in)        :: EV
  double precision,intent(in)        :: EJ
  double precision,intent(in)        :: Ex
  double precision,intent(in)        :: Ec
  double precision,intent(in)        :: Ew

  integer                            :: HOMO
  integer                            :: LUMO
  double precision                   :: Gap

! HOMO, LUMO, and Gap

  HOMO = nO

  LUMO = HOMO + 1

  Gap = eps(LUMO) - eps(HOMO)

! Dump results


  write(*,*)
  write(*,'(A60)')              '-------------------------------------------------'
  write(*,'(A40)')              ' Summary              '
  write(*,'(A60)')              '-------------------------------------------------'
  write(*,'(A40,1X,F16.10,A3)') ' One-electron    energy: ',ET + EV,' au'
  write(*,'(A40,1X,F16.10,A3)') ' Kinetic         energy: ',ET,' au'
  write(*,'(A40,1X,F16.10,A3)') ' Potential       energy: ',EV,' au'
  write(*,'(A60)')              '-------------------------------------------------'
  write(*,'(A40,1X,F16.10,A3)') ' Two-electron    energy: ',EJ + Ex + Ec,' au'
  write(*,'(A40,1X,F16.10,A3)') ' Coulomb         energy: ',EJ,' au'
  write(*,'(A40,1X,F16.10,A3)') ' Exchange        energy: ',Ex,' au'
  write(*,'(A40,1X,F16.10,A3)') ' Correlation     energy: ',Ec,' au'
  write(*,'(A60)')              '-------------------------------------------------'
  write(*,'(A40,1X,F16.10,A3)') ' Electronic      energy: ',Ew       ,' au'
  write(*,'(A40,1X,F16.10,A3)') ' Nuclear      repulsion: ',     ENuc,' au'
  write(*,'(A40,1X,F16.10,A3)') ' Kohn-Sham       energy: ',Ew + ENuc,' au'
  write(*,'(A60)')              '-------------------------------------------------'
  write(*,'(A40,F13.6,A3)')     ' KS HOMO      energy:',eps(HOMO),' au'
  write(*,'(A40,F13.6,A3)')     ' KS LUMO      energy:',eps(LUMO),' au'
  write(*,'(A40,F13.6,A3)')     ' KS HOMO-LUMO    gap:',Gap      ,' au'
  write(*,'(A60)')              '-------------------------------------------------'
  write(*,'(A40,F13.6,A3)')     ' KS HOMO      energy:',eps(HOMO)*HatoeV,' eV'
  write(*,'(A40,F13.6,A3)')     ' KS LUMO      energy:',eps(LUMO)*HatoeV,' eV'
  write(*,'(A40,F13.6,A3)')     ' KS HOMO-LUMO    gap:',Gap      *HatoeV,' eV'
  write(*,'(A60)')              '-------------------------------------------------'
  write(*,*)

! Print results

  write(*,'(A50)') '-----------------------------------------'
  write(*,'(A50)') '     Kohn-Sham orbital coefficients      '
  write(*,'(A50)') '-----------------------------------------'
  call matout(nBas,nBas,c(:,:))
  write(*,*)
  write(*,'(A50)') '---------------------------------------'
  write(*,'(A50)') '     Kohn-Sham orbital energies        '
  write(*,'(A50)') '---------------------------------------'
  call matout(nBas,1,eps(:))
  write(*,*)

end subroutine print_RKS
