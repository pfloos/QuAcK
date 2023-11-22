subroutine print_UHF(nBas,nO,Ov,eHF,c,ENuc,ET,EV,EJ,Ex,EUHF,dipole)

! Print one- and two-electron energies and other stuff for UHF calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nO(nspin)
  double precision,intent(in)        :: Ov(nBas,nBas)
  double precision,intent(in)        :: eHF(nBas,nspin)
  double precision,intent(in)        :: c(nBas,nBas,nspin)
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ET(nspin)
  double precision,intent(in)        :: EV(nspin)
  double precision,intent(in)        :: EJ(nsp)
  double precision,intent(in)        :: Ex(nspin)
  double precision,intent(in)        :: EUHF
  double precision,intent(in)        :: dipole(ncart)

! Local variables

  integer                            :: ixyz
  integer                            :: ispin
  double precision                   :: HOMO(nspin)
  double precision                   :: LUMO(nspin)
  double precision                   :: Gap(nspin)
  double precision                   :: Sz
  double precision                   :: Sx2,Sy2,Sz2

  logical                            :: dump_orb = .false.

! HOMO and LUMO

  do ispin=1,nspin
    if(nO(ispin) > 0) then 
      HOMO(ispin) = eHF(nO(ispin),ispin)
      if(nO(ispin) < nBas) then
        LUMO(ispin) = eHF(nO(ispin)+1,ispin)
      else
        LUMO(ispin) = 0d0
      end if
      Gap(ispin)  = LUMO(ispin) - HOMO(ispin)
    else
      HOMO(ispin) = 0d0
      LUMO(ispin) = eHF(1,ispin)
      Gap(ispin)  = 0d0
    end if
  end do

  Sz =  0.5d0*dble(nO(1) - nO(2))
  Sx2 = 0.25d0*dble(nO(1) - nO(2)) + 0.5d0*nO(2) - 0.5d0*sum(matmul(transpose(c(:,1:nO(1),1)),matmul(Ov,c(:,1:nO(2),2)))**2)
  Sy2 = 0.25d0*dble(nO(1) - nO(2)) + 0.5d0*nO(2) - 0.5d0*sum(matmul(transpose(c(:,1:nO(1),1)),matmul(Ov,c(:,1:nO(2),2)))**2)
  Sz2 = 0.25d0*dble(nO(1) - nO(2))**2

! Dump results

  write(*,*)
  write(*,'(A60)')              '---------------------------------------------'
  write(*,'(A40)')              ' Summary                  '
  write(*,'(A60)')              '---------------------------------------------'
  write(*,'(A40,1X,F16.10,A3)') ' One-electron    energy = ',sum(ET)  + sum(EV),' au'
  write(*,'(A40,1X,F16.10,A3)') ' One-electron a  energy = ',ET(1) + EV(1),' au'
  write(*,'(A40,1X,F16.10,A3)') ' One-electron b  energy = ',ET(2) + EV(2),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Kinetic         energy = ',sum(ET),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Kinetic      a  energy = ',ET(1),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Kinetic      b  energy = ',ET(2),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Potential       energy = ',sum(EV),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Potential    a  energy = ',EV(1),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Potential    b  energy = ',EV(2),' au'
  write(*,'(A60)')              '---------------------------------------------'
  write(*,'(A40,1X,F16.10,A3)') ' Two-electron    energy = ',sum(EJ) + sum(Ex),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Two-electron aa energy = ',EJ(1) + Ex(1),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Two-electron ab energy = ',EJ(2),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Two-electron bb energy = ',EJ(3) + Ex(2),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Hartree         energy = ',sum(EJ),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Hartree      aa energy = ',EJ(1),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Hartree      ab energy = ',EJ(2),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Hartree      bb energy = ',EJ(3),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Exchange        energy = ',sum(Ex),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Exchange     a  energy = ',Ex(1),' au'
  write(*,'(A40,1X,F16.10,A3)') ' Exchange     b  energy = ',Ex(2),' au'
  write(*,'(A60)')              '---------------------------------------------'
  write(*,'(A40,1X,F16.10,A3)') ' Electronic      energy = ',EUHF,' au'
  write(*,'(A40,1X,F16.10,A3)') ' Nuclear      repulsion = ',ENuc,' au'
  write(*,'(A40,1X,F16.10,A3)') ' UHF             energy = ',EUHF + ENuc,' au'
  write(*,'(A60)')              '---------------------------------------------'
  write(*,'(A40,1X,F16.6,A3)')  ' UHF HOMO a    energy   = ' ,HOMO(1)*HatoeV,' eV'
  write(*,'(A40,1X,F16.6,A3)')  ' UHF LUMO a    energy   = ' ,LUMO(1)*HatoeV,' eV'
  write(*,'(A40,1X,F16.6,A3)')  ' UHF HOMOa-LUMOa  gap   = ' ,Gap(1)*HatoeV,' eV'
  write(*,'(A60)')              '---------------------------------------------'
  write(*,'(A40,1X,F16.6,A3)')  ' UHF HOMO b    energy   = ',HOMO(2)*HatoeV,' eV'
  write(*,'(A40,1X,F16.6,A3)')  ' UHF LUMO b    energy   = ',LUMO(2)*HatoeV,' eV'
  write(*,'(A40,1X,F16.6,A3)')  ' UHF HOMOb-LUMOb  gap   = ',Gap(2)*HatoeV,' eV'
  write(*,'(A60)')              '---------------------------------------------'
  write(*,'(A40,1X,F10.6)')     ' <Sz>                   = ',Sz
  write(*,'(A40,1X,F10.6)')     ' <S^2>                  = ',Sx2+Sy2+Sz2
  write(*,'(A60)')              '---------------------------------------------'
  write(*,'(A45)')              ' Dipole moment (Debye)    '
  write(*,'(19X,4A10)')         'X','Y','Z','Tot.'
  write(*,'(19X,4F10.4)')       (dipole(ixyz)*auToD,ixyz=1,ncart),norm2(dipole)*auToD
  write(*,'(A60)')              '---------------------------------------------'
  write(*,*)

! Print results

  if(dump_orb) then
    write(*,'(A40)') '-----------------------------------------'
    write(*,'(A40)') 'UHF spin-up   orbital coefficients '
    write(*,'(A40)') '-----------------------------------------'
    call matout(nBas,nBas,c(:,:,1))
    write(*,*)
    write(*,'(A40)') '-----------------------------------------'
    write(*,'(A40)') 'UHF spin-down orbital coefficients '
    write(*,'(A40)') '-----------------------------------------'
    call matout(nBas,nBas,c(:,:,2))
    write(*,*)
  end if
  write(*,'(A40)') '---------------------------------------'
  write(*,'(A40)') ' UHF spin-up   orbital energies  '
  write(*,'(A40)') '---------------------------------------'
  call vecout(nBas,eHF(:,1))
  write(*,*)
  write(*,'(A40)') '---------------------------------------'
  write(*,'(A40)') ' UHF spin-down orbital energies  '
  write(*,'(A40)') '---------------------------------------'
  call vecout(nBas,eHF(:,2))
  write(*,*)

end subroutine 
