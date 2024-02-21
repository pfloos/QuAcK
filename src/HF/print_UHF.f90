subroutine print_UHF(nBas,nO,S,eHF,c,P,ENuc,ET,EV,EJ,Ex,EUHF,dipole)

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

! Local variables

  integer                            :: ixyz
  integer                            :: ispin
  double precision                   :: eHOMO(nspin)
  double precision                   :: eLUMO(nspin)
  double precision                   :: Gap
  double precision                   :: Sz
  double precision                   :: Sx2,Sy2,Sz2
  integer                            :: mu,nu
  double precision,allocatable       :: qa(:),qb(:)

  logical                            :: dump_orb = .false.

! HOMO and LUMO

      
  do ispin=1,nspin
    eHOMO(ispin) = maxval(eHF(1:nO(ispin),ispin))
    eLUMO(ispin) = minval(eHF(nO(ispin)+1:nBas,ispin))
  end do
  Gap = minval(eLUMO)  -maxval(eHOMO)

  Sz =  0.5d0*dble(nO(1) - nO(2))
  Sx2 = 0.25d0*dble(nO(1) - nO(2)) + 0.5d0*nO(2) - 0.5d0*sum(matmul(transpose(c(:,1:nO(1),1)),matmul(S,c(:,1:nO(2),2)))**2)
  Sy2 = 0.25d0*dble(nO(1) - nO(2)) + 0.5d0*nO(2) - 0.5d0*sum(matmul(transpose(c(:,1:nO(1),1)),matmul(S,c(:,1:nO(2),2)))**2)
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
  write(*,'(A40)') ' UHF spin-up   orbital energies  '
  write(*,'(A40)') '---------------------------------------'
  call vecout(nBas,eHF(:,1))
  write(*,*)
  write(*,'(A40)') '---------------------------------------'
  write(*,'(A40)') ' UHF spin-down orbital energies  '
  write(*,'(A40)') '---------------------------------------'
  call vecout(nBas,eHF(:,2))
  write(*,*)

  allocate(qa(nBas),qb(nBas))

  qa(:) = 0d0
  qb(:) = 0d0
  do mu=1,nBas
    do nu=1,nBas
      qa(mu) = qa(mu) + P(mu,nu,1)*S(nu,mu)
      qb(mu) = qb(mu) + P(mu,nu,2)*S(nu,mu)
    end do
  end do

  call vecout(nBas,qa)
  call vecout(nBas,qb)

end subroutine 
