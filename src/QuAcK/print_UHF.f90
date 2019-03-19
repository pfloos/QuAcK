subroutine print_UHF(nBas,nO,eps,c,ENuc,ET,EV,EJ,Ex,Ec,Ew)

! Print one- and two-electron energies and other stuff for UHF calculation

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nO(nspin)
  double precision,intent(in)        :: eps(nBas,nspin)
  double precision,intent(in)        :: c(nBas,nBas,nspin)
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ET(nspin)
  double precision,intent(in)        :: EV(nspin)
  double precision,intent(in)        :: EJ(nsp)
  double precision,intent(in)        :: Ex(nspin)
  double precision,intent(in)        :: Ec(nsp)
  double precision,intent(in)        :: Ew

  integer                            :: HOMO(nspin)
  integer                            :: LUMO(nspin)
  double precision                   :: Gap(nspin)

! HOMO and LUMO

  HOMO(:) = nO(:)

  LUMO(:) = HOMO(:) + 1

  Gap(1) = eps(LUMO(1),1) - eps(HOMO(1),1)
  Gap(2) = eps(LUMO(2),2) - eps(HOMO(2),2)

! Dump results


  write(*,*)
  write(*,'(A60)')              '-------------------------------------------------'
  write(*,'(A40)')              ' Summary              '
  write(*,'(A60)')              '-------------------------------------------------'
  write(*,'(A40,1X,F16.10,A3)') ' One-electron    energy: ',sum(ET(:))  + sum(EV(:)),' au'
  write(*,'(A40,1X,F16.10,A3)') ' One-electron a  energy: ',ET(1) + EV(1),' au'
  write(*,'(A40,1X,F16.10,A3)') ' One-electron b  energy: ',ET(2) + EV(2),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Kinetic         energy: ',sum(ET(:)),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Kinetic      a  energy: ',ET(1),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Kinetic      b  energy: ',ET(2),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Potential       energy: ',sum(EV(:)),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Potential    a  energy: ',EV(1),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Potential    b  energy: ',EV(2),' au'
  write(*,'(A60)')              '-------------------------------------------------'
  write(*,'(A40,1X,F16.10,A3)') ' Two-electron a  energy: ',sum(EJ(:)) + sum(Ex(:))  + sum(Ec(:)),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Two-electron aa energy: ',EJ(1) + Ex(1) + Ec(1),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Two-electron ab energy: ',EJ(2) + Ec(2),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Two-electron bb energy: ',EJ(3) + Ex(2) + Ec(3),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Coulomb         energy: ',sum(EJ(:)),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Coulomb      aa energy: ',EJ(1),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Coulomb      ab energy: ',EJ(2),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Coulomb      bb energy: ',EJ(3),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Exchange        energy: ',sum(Ex(:)),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Exchange     a  energy: ',Ex(1),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Exchange     b  energy: ',Ex(2),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Correlation     energy: ',sum(Ec(:)),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Correlation  aa energy: ',Ec(1),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Correlation  ab energy: ',Ec(2),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Correlation  bb energy: ',Ec(3),' au'
  write(*,'(A60)')              '-------------------------------------------------'
  write(*,'(A40,1X,F16.10,A3)') ' Electronic      energy: ',Ew,' au'
  write(*,'(A40,1X,F16.10,A3)') ' Nuclear      repulsion: ',ENuc,' au'
  write(*,'(A40,1X,F16.10,A3)') ' UHF             energy: ',Ew + ENuc,' au'
  write(*,'(A60)')              '-------------------------------------------------'
  write(*,'(A40,F13.6,A3)')     ' UHF HOMO a    energy:',eps(HOMO(1),1)*HatoeV,' eV'
  write(*,'(A40,F13.6,A3)')     ' UHF LUMO a    energy:',eps(LUMO(1),1)*HatoeV,' eV'
  write(*,'(A40,F13.6,A3)')     ' UHF HOMOa-LUMOa  gap:',Gap(1)*HatoeV,' eV'
  write(*,'(A60)')              '-------------------------------------------------'
  write(*,'(A40,F13.6,A3)')     ' UHF HOMO b    energy:',eps(HOMO(2),2)*HatoeV,' eV'
  write(*,'(A40,F13.6,A3)')     ' UHF LUMO b    energy:',eps(LUMO(2),2)*HatoeV,' eV'
  write(*,'(A40,F13.6,A3)')     ' UHF HOMOb-LUMOb gap :',Gap(2)*HatoeV,' eV'
  write(*,'(A60)')              '-------------------------------------------------'
  write(*,*)

! Print results

  write(*,'(A50)') '-----------------------------------------'
  write(*,'(A50)') 'UHF spin-up   orbital coefficients '
  write(*,'(A50)') '-----------------------------------------'
  call matout(nBas,nBas,c(:,:,1))
  write(*,'(A50)') '-----------------------------------------'
  write(*,'(A50)') 'UHF spin-down orbital coefficients '
  write(*,'(A50)') '-----------------------------------------'
  call matout(nBas,nBas,c(:,:,2))
  write(*,*)
  write(*,'(A50)') '---------------------------------------'
  write(*,'(A50)') ' UHF spin-up   orbital energies  '
  write(*,'(A50)') '---------------------------------------'
  call matout(nBas,1,eps(:,1))
  write(*,*)
  write(*,'(A50)') '---------------------------------------'
  write(*,'(A50)') ' UHF spin-down orbital energies  '
  write(*,'(A50)') '---------------------------------------'
  call matout(nBas,1,eps(:,2))
  write(*,*)

end subroutine print_UHF
