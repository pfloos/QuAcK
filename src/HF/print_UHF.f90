subroutine print_UHF(nBas,nO,e,c,ENuc,ET,EV,EJ,Ex,EUHF)

! Print one- and two-electron energies and other stuff for UHF calculation

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nO(nspin)
  double precision,intent(in)        :: e(nBas,nspin)
  double precision,intent(in)        :: c(nBas,nBas,nspin)
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ET(nspin)
  double precision,intent(in)        :: EV(nspin)
  double precision,intent(in)        :: EJ(nsp)
  double precision,intent(in)        :: Ex(nspin)
  double precision,intent(in)        :: EUHF

  integer                            :: ispin
  double precision                   :: HOMO(nspin)
  double precision                   :: LUMO(nspin)
  double precision                   :: Gap(nspin)

! HOMO and LUMO

  do ispin=1,nspin
    if(nO(ispin) > 0) then 
      HOMO(ispin) = e(nO(ispin),ispin)
      LUMO(ispin) = e(nO(ispin)+1,ispin)
      Gap(ispin)  = LUMO(ispin) - HOMO(ispin)
    else
      HOMO(ispin) = 0d0
      LUMO(ispin) = e(1,ispin)
      Gap(ispin)  = 0d0
    end if
  end do

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
  write(*,'(A40,1X,F16.10,A3)') ' Two-electron a  energy: ',sum(EJ(:)) + sum(Ex(:)),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Two-electron aa energy: ',EJ(1) + Ex(1),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Two-electron ab energy: ',EJ(2),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Two-electron bb energy: ',EJ(3) + Ex(2),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Coulomb         energy: ',sum(EJ(:)),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Coulomb      aa energy: ',EJ(1),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Coulomb      ab energy: ',EJ(2),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Coulomb      bb energy: ',EJ(3),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Exchange        energy: ',sum(Ex(:)),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Exchange     a  energy: ',Ex(1),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Exchange     b  energy: ',Ex(2),' au'
  write(*,'(A60)')              '-------------------------------------------------'
  write(*,'(A40,1X,F16.10,A3)') ' Electronic      energy: ',EUHF,' au'
  write(*,'(A40,1X,F16.10,A3)') ' Nuclear      repulsion: ',ENuc,' au'
  write(*,'(A40,1X,F16.10,A3)') ' UHF             energy: ',EUHF + ENuc,' au'
  write(*,'(A60)')              '-------------------------------------------------'
  write(*,'(A40,F13.6,A3)')     ' UHF HOMO a    energy:',HOMO(1)*HatoeV,' eV'
  write(*,'(A40,F13.6,A3)')     ' UHF LUMO a    energy:',LUMO(1)*HatoeV,' eV'
  write(*,'(A40,F13.6,A3)')     ' UHF HOMOa-LUMOa  gap:',Gap(1)*HatoeV,' eV'
  write(*,'(A60)')              '-------------------------------------------------'
  write(*,'(A40,F13.6,A3)')     ' UHF HOMO b    energy:',HOMO(2)*HatoeV,' eV'
  write(*,'(A40,F13.6,A3)')     ' UHF LUMO b    energy:',LUMO(2)*HatoeV,' eV'
  write(*,'(A40,F13.6,A3)')     ' UHF HOMOb-LUMOb gap :',Gap(2)*HatoeV,' eV'
  write(*,'(A60)')              '-------------------------------------------------'
  write(*,*)

! Print results

  write(*,'(A50)') '-----------------------------------------'
  write(*,'(A50)') 'UHF spin-up   orbital coefficients '
  write(*,'(A50)') '-----------------------------------------'
  call matout(nBas,nBas,c(:,:,1))
  write(*,*)
  write(*,'(A50)') '-----------------------------------------'
  write(*,'(A50)') 'UHF spin-down orbital coefficients '
  write(*,'(A50)') '-----------------------------------------'
  call matout(nBas,nBas,c(:,:,2))
  write(*,*)
  write(*,'(A50)') '---------------------------------------'
  write(*,'(A50)') ' UHF spin-up   orbital energies  '
  write(*,'(A50)') '---------------------------------------'
  call matout(nBas,1,e(:,1))
  write(*,*)
  write(*,'(A50)') '---------------------------------------'
  write(*,'(A50)') ' UHF spin-down orbital energies  '
  write(*,'(A50)') '---------------------------------------'
  call matout(nBas,1,e(:,2))
  write(*,*)

end subroutine print_UHF