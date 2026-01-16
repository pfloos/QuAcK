subroutine print_MOM_UHF(nBas,nO,S,eHF,c,P,ENuc,ET,EV,EJ,Ex,EUHF,dipole,occupations)

! Print one- and two-electron energies and other stuff for UHF calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nO(nspin)
  double precision,intent(in)        :: S(nBas,nBas)
  double precision,intent(in)        :: eHF(nBas,nspin)
  double precision,intent(in)        :: c(nBas,nBas,nspin)
  double precision,intent(in)        :: P(nBas,nBas,nspin)
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ET(nspin)
  double precision,intent(in)        :: EV(nspin)
  double precision,intent(in)        :: EJ(nsp)
  double precision,intent(in)        :: Ex(nspin)
  double precision,intent(in)        :: EUHF
  double precision,intent(in)        :: dipole(ncart)
  integer,intent(in)                 :: occupations(maxval(nO),nspin)

! Local variables

  integer                            :: ixyz
  integer                            :: ispin
  double precision                   :: eHOMO(nspin)
  double precision                   :: eLUMO(nspin)
  double precision                   :: Gap
  double precision                   :: Sz
  double precision                   :: Sx2,Sy2,Sz2
  integer                            :: mu,nu
  integer,allocatable                :: unoccupied(:,:)

  logical                            :: dump_orb = .false.

allocate(unoccupied(nBas - minval(nO),nspin))

! HOMO and LUMO
  do ispin=1,nspin
    call non_occupied(nO(ispin),nBas,occupations(1:nO(ispin),ispin),unoccupied(1:nBas - nO(ispin),ispin))
    eHOMO(ispin) = maxval(eHF(occupations(1:nO(ispin),ispin),ispin))
    eLUMO(ispin) = minval(eHF(unoccupied(1:nBas-nO(ispin),ispin),ispin))
  end do
  Gap = minval(eLUMO)  -maxval(eHOMO)

  Sz =  0.5d0*dble(nO(1) - nO(2))
  Sx2 = 0.25d0*dble(nO(1) - nO(2)) + 0.5d0*nO(2) - &
        0.5d0*sum(matmul(transpose(c(:,occupations(1:nO(1),1),1)),&
                  matmul(S,c(:,occupations(1:nO(2),2),2)))**2)
  Sy2 = 0.25d0*dble(nO(1) - nO(2)) + 0.5d0*nO(2) - &
        0.5d0*sum(matmul(transpose(c(:,occupations(1:nO(1),1),1)),&
                  matmul(S,c(:,occupations(1:nO(2),2),2)))**2)
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
  write(*,'(A40,1X,F16.6,A3)')  ' UHF HOMO      energy   = ' ,maxval(eHOMO)*HatoeV,' eV'
  write(*,'(A40,1X,F16.6,A3)')  ' UHF LUMO      energy   = ' ,minval(eLUMO)*HatoeV,' eV'
  write(*,'(A40,1X,F16.6,A3)')  ' UHF HOMO-LUMO    gap   = ' ,Gap*HatoeV,' eV'
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
  write(*,'(A40)') ' UHF Occupied alpha orbital energies  '
  write(*,'(A40)') '---------------------------------------'
  call vecout(nO(1),eHF(occupations(1:nO(1),1),1))
  write(*,*)
  write(*,'(A40)') '---------------------------------------'
  write(*,'(A40)') ' UHF Unoccupied alpha orbital energies '
  write(*,'(A40)') '---------------------------------------'
  call vecout(nBas - nO(1),eHF(unoccupied(1:nBas - nO(1),1),1))
  write(*,*)
  write(*,'(A40)') '---------------------------------------'
  write(*,'(A40)') ' UHF Occupied beta orbital energies  '
  write(*,'(A40)') '---------------------------------------'
  call vecout(nO(2),eHF(occupations(1:nO(2),2),2))
  write(*,*)
  write(*,'(A40)') '---------------------------------------'
  write(*,'(A40)') ' UHF Unoccupied beta orbital energies '
  write(*,'(A40)') '---------------------------------------'
  call vecout(nBas - nO(2),eHF(unoccupied(1:nBas-nO(2),2),2))
  write(*,*)

  print *, "Orbital occupations for MOMUHF:"
  print *, "Alpha:"
  print *, occupations(1:nO(1),1)
  print *, "Beta:"
  print *, occupations(1:nO(2),2)
  print *, ""

  deallocate(unoccupied)
end subroutine 
