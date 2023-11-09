subroutine print_ROHF(nBas,nO,e,c,ENuc,ET,EV,EJ,Ex,EHF,dipole)

! Print one- and two-electron energies and other stuff for RoHF calculation

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nO(nspin)
  double precision,intent(in)        :: e(nBas)
  double precision,intent(in)        :: c(nBas,nBas)
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ET(nspin)
  double precision,intent(in)        :: EV(nspin)
  double precision,intent(in)        :: EJ(nsp)
  double precision,intent(in)        :: Ex(nspin)
  double precision,intent(in)        :: EHF
  double precision,intent(in)        :: dipole(ncart)

  integer                            :: ixyz
  integer                            :: ispin
  double precision                   :: HOMO(nspin)
  double precision                   :: LUMO(nspin)
  double precision                   :: Gap(nspin)
  double precision                   :: S,S2

! HOMO and LUMO

  do ispin=1,nspin
    if(nO(ispin) > 0) then 
      HOMO(ispin) = e(nO(ispin))
      if(nO(ispin) < nBas) then
        LUMO(ispin) = e(nO(ispin)+1)
      else
        LUMO(ispin) = 0d0
      end if
      Gap(ispin)  = LUMO(ispin) - HOMO(ispin)
    else
      HOMO(ispin) = 0d0
      LUMO(ispin) = e(1)
      Gap(ispin)  = 0d0
    end if
  end do

  S2 = dble(nO(1) - nO(2))/2d0*(dble(nO(1) - nO(2))/2d0 + 1d0) 
  S  = 0.5d0*dble(nO(1) - nO(2))

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
  write(*,'(A40,1X,F16.10,A3)') ' Two-electron    energy: ',sum(EJ(:)) + sum(Ex(:)),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Two-electron aa energy: ',EJ(1) + Ex(1),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Two-electron ab energy: ',EJ(2),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Two-electron bb energy: ',EJ(3) + Ex(2),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Hartree         energy: ',sum(EJ(:)),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Hartree      aa energy: ',EJ(1),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Hartree      ab energy: ',EJ(2),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Hartree      bb energy: ',EJ(3),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Exchange        energy: ',sum(Ex(:)),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Exchange     a  energy: ',Ex(1),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Exchange     b  energy: ',Ex(2),' au'
  write(*,'(A60)')              '-------------------------------------------------'
  write(*,'(A40,1X,F16.10,A3)') ' Electronic      energy: ',EHF,' au'
  write(*,'(A40,1X,F16.10,A3)') ' Nuclear      repulsion: ',ENuc,' au'
  write(*,'(A40,1X,F16.10,A3)') ' ROHF            energy: ',EHF + ENuc,' au'
  write(*,'(A60)')              '-------------------------------------------------'
  write(*,'(A40,1X,F16.6,A3)')  ' ROHF HOMO a   energy:',HOMO(1)*HatoeV,' eV'
  write(*,'(A40,1X,F16.6,A3)')  ' ROHF LUMO a   energy:',LUMO(1)*HatoeV,' eV'
  write(*,'(A40,1X,F16.6,A3)')  ' ROHF HOMOa-LUMOa gap:',Gap(1)*HatoeV,' eV'
  write(*,'(A60)')              '-------------------------------------------------'
  write(*,'(A40,1X,F16.6,A3)')  ' ROHF HOMO b   energy:',HOMO(2)*HatoeV,' eV'
  write(*,'(A40,1X,F16.6,A3)')  ' ROHF LUMO b   energy:',LUMO(2)*HatoeV,' eV'
  write(*,'(A40,1X,F16.6,A3)')  ' ROHF HOMOb-LUMOb gap:',Gap(2)*HatoeV,' eV'
  write(*,'(A60)')              '-------------------------------------------------'
  write(*,'(A40,1X,F16.6)')     '  S                  :',2d0*S       + 1d0
  write(*,'(A40,1X,F16.6)')     ' <S**2>              :',S2
  write(*,'(A60)')              '-------------------------------------------------'
  write(*,'(A45)')              ' Dipole moment (Debye)    '
  write(*,'(19X,4A10)')         'X','Y','Z','Tot.'
  write(*,'(19X,4F10.4)')       (dipole(ixyz)*auToD,ixyz=1,ncart),norm2(dipole)*auToD
  write(*,'(A60)')              '-------------------------------------------------'
  write(*,*)

! Print results

  write(*,'(A50)') '-----------------------------------------'
  write(*,'(A50)') 'ROHF orbital coefficients '
  write(*,'(A50)') '-----------------------------------------'
  call matout(nBas,nBas,c)
  write(*,*)
  write(*,'(A50)') '---------------------------------------'
  write(*,'(A50)') ' ROHF orbital energies (au) '
  write(*,'(A50)') '---------------------------------------'
  call matout(nBas,1,e)
  write(*,*)

end subroutine 
