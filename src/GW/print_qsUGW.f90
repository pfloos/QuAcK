subroutine print_qsUGW(nBas,nO,nSCF,Conv,thresh,eHF,eGW,c,Ov,ENuc,ET,EV,EJ,Ex,EcGM,EcRPA,EqsGW,SigC,Z,dipole)

! Print one-electron energies and other stuff for qsUGW

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nO(nspin)
  integer,intent(in)                 :: nSCF
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ET(nspin)
  double precision,intent(in)        :: EV(nspin)
  double precision,intent(in)        :: EJ(nsp)
  double precision,intent(in)        :: Ex(nspin)
  double precision,intent(in)        :: EcGM(nspin)
  double precision,intent(in)        :: EcRPA
  double precision,intent(in)        :: EqsGW
  double precision,intent(in)        :: Conv
  double precision,intent(in)        :: thresh
  double precision,intent(in)        :: eHF(nBas,nspin)
  double precision,intent(in)        :: eGW(nBas,nspin)
  double precision,intent(in)        :: c(nBas,nBas,nspin)
  double precision,intent(in)        :: Ov(nBas,nBas)
  double precision,intent(in)        :: SigC(nBas,nBas,nspin)
  double precision,intent(in)        :: Z(nBas,nspin)
  double precision,intent(in)        :: dipole(ncart)

! Local variables

  logical                            :: dump_orb = .false.
  integer                            :: p
  integer                            :: ispin,ixyz
  double precision                   :: HOMO(nspin)
  double precision                   :: LUMO(nspin)
  double precision                   :: Gap(nspin)
  double precision                   :: Sz
  double precision                   :: Sx2,Sy2,Sz2
  double precision,external          :: trace_matrix

! HOMO and LUMO

  do ispin=1,nspin
    if(nO(ispin) > 0) then
      HOMO(ispin) = eGW(nO(ispin),ispin)
      LUMO(ispin) = eGW(nO(ispin)+1,ispin)
      Gap(ispin)  = LUMO(ispin) - HOMO(ispin)
    else
      HOMO(ispin) = 0d0
      LUMO(ispin) = eGW(1,ispin)
      Gap(ispin)  = 0d0
    end if
  end do

  Sz =  0.5d0*dble(nO(1) - nO(2))
  Sx2 = 0.25d0*dble(nO(1) - nO(2)) + 0.5d0*nO(2) - 0.5d0*sum(matmul(transpose(c(:,1:nO(1),1)),matmul(Ov,c(:,1:nO(2),2)))**2)
  Sy2 = 0.25d0*dble(nO(1) - nO(2)) + 0.5d0*nO(2) - 0.5d0*sum(matmul(transpose(c(:,1:nO(1),1)),matmul(Ov,c(:,1:nO(2),2)))**2)
  Sz2 = 0.25d0*dble(nO(1) - nO(2))**2

! Dump results

  write(*,*)'----------------------------------------------------------------'// &
            '----------------------------------------------------------------'
  if(nSCF < 10) then
    write(*,'(1X,A20,I1,A1,I1,A16)')' Self-consistent qsG',nSCF,'W',nSCF,'@UHF calculation'
  elseif(nSCF < 100) then
    write(*,'(1X,A20,I2,A1,I2,A16)')' Self-consistent qsG',nSCF,'W',nSCF,'@UHF calculation'
  else
    write(*,'(1X,A20,I3,A1,I3,A16)')' Self-consistent qsG',nSCF,'W',nSCF,'@UHF calculation'
  end if
  write(*,*)'----------------------------------------------------------------'// &
            '----------------------------------------------------------------'
  write(*,'(A1,A3,A1,A30,A1,A30,A1,A30,A1,A30,A1)') &
            '|',' ','|','e_HF (eV)         ','|','Sig_GW (eV)        ','|','Z             ','|','e_GW (eV)         ','|'
  write(*,'(A1,A3,A1,2A15,A1,2A15,A1,2A15,A1,2A15,A1)') &
            '|','#','|','up     ','dw     ','|','up     ','dw     ','|','up     ','dw     ','|','up     ','dw     ','|'
  write(*,*)'----------------------------------------------------------------'// &
            '----------------------------------------------------------------'

  do p=1,nBas
    write(*,'(A1,I3,A1,2F15.6,A1,2F15.6,A1,2F15.6,A1,2F15.6,A1)') &
    '|',p,'|',eHF(p,1)*HaToeV,eHF(p,2)*HaToeV,'|',SigC(p,p,1)*HaToeV,SigC(p,p,2)*HaToeV,'|', &
              Z(p,1),Z(p,2),'|',eGW(p,1)*HaToeV,eGW(p,2)*HaToeV,'|'
  end do

  write(*,*)'----------------------------------------------------------------'// &
            '----------------------------------------------------------------'
  write(*,'(2X,A10,I3)')   'Iteration ',nSCF
  write(*,'(2X,A14,F15.5)')'Convergence = ',Conv
  write(*,*)'----------------------------------------------------------------'// &
            '----------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'qsGW@UHF HOMO      energy = ',maxval(HOMO)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'qsGW@UHF LUMO      energy = ',minval(LUMO)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'qsGW@UHF HOMO-LUMO gap    = ',(minval(LUMO)-maxval(HOMO))*HaToeV,' eV'
  write(*,*)'----------------------------------------------------------------'// &
            '----------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') '    qsGW@UHF total       energy = ',ENuc + EqsGW,' au'
  write(*,'(2X,A60,F15.6,A3)') '    qsGW@UHF exchange    energy = ',sum(Ex),' au'
  write(*,'(2X,A60,F15.6,A3)') ' GM@qsGW@UHF correlation energy = ',sum(EcGM),' au'
  write(*,'(2X,A60,F15.6,A3)') 'RPA@qsGW@UHF correlation energy = ',EcRPA,' au'
  write(*,*)'----------------------------------------------------------------'// &
            '----------------------------------------------------------------'
  write(*,*)

! Dump results for final iteration

  if(Conv < thresh) then

    write(*,*)
    write(*,'(A60)')              '-------------------------------------------------'
    write(*,'(A40)')              ' Summary              '
    write(*,'(A60)')              '-------------------------------------------------'
    write(*,'(A40,1X,F16.10,A3)') ' One-electron    energy = ',sum(ET)  + sum(EV),' au'
    write(*,'(A40,1X,F16.10,A3)') ' One-electron a  energy = ',ET(1) + EV(1),' au'
    write(*,'(A40,1X,F16.10,A3)') ' One-electron b  energy = ',ET(2) + EV(2),' au'
    write(*,*)
    write(*,'(A40,1X,F16.10,A3)') ' Kinetic         energy = ',sum(ET),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Kinetic      a  energy = ',ET(1),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Kinetic      b  energy = ',ET(2),' au'
    write(*,*)
    write(*,'(A40,1X,F16.10,A3)') ' Potential       energy = ',sum(EV),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Potential    a  energy = ',EV(1),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Potential    b  energy = ',EV(2),' au'
    write(*,'(A60)')              '-------------------------------------------------'
    write(*,'(A40,1X,F16.10,A3)') ' Two-electron    energy = ',sum(EJ) + sum(Ex),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Two-electron aa energy = ',EJ(1) + Ex(1),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Two-electron ab energy = ',EJ(2),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Two-electron bb energy = ',EJ(3) + Ex(2),' au'
    write(*,*)
    write(*,'(A40,1X,F16.10,A3)') ' Hartree         energy = ',sum(EJ),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Hartree      aa energy = ',EJ(1),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Hartree      ab energy = ',EJ(2),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Hartree      bb energy = ',EJ(3),' au'
    write(*,*)
    write(*,'(A40,1X,F16.10,A3)') ' Exchange        energy = ',sum(Ex),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Exchange     a  energy = ',Ex(1),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Exchange     b  energy = ',Ex(2),' au'
    write(*,*)
    write(*,'(A40,1X,F16.10,A3)') ' Correlation     energy = ',sum(EcGM),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Correlation  aa energy = ',EcGM(1),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Correlation  bb energy = ',EcGM(2),' au'
    write(*,'(A60)')              '-------------------------------------------------'
    write(*,'(A40,1X,F16.10,A3)') ' Electronic      energy = ',EqsGW,' au'
    write(*,'(A40,1X,F16.10,A3)') ' Nuclear      repulsion = ',ENuc,' au'
    write(*,'(A40,1X,F16.10,A3)') ' qsUGW           energy = ',ENuc + EqsGW,' au'
    write(*,'(A60)')              '-------------------------------------------------'
    write(*,'(A40,1X,F10.6)')     ' <Sz>                   = ',Sz
    write(*,'(A40,1X,F10.6)')     ' <S^2>                  = ',Sx2+Sy2+Sz2
    write(*,'(A60)')              '-------------------------------------------------'
    write(*,'(A45)')              ' Dipole moment (Debye)    '
    write(*,'(19X,4A10)')         'X','Y','Z','Tot.'
    write(*,'(19X,4F10.4)')       (dipole(ixyz)*auToD,ixyz=1,ncart),norm2(dipole)*auToD
    write(*,'(A60)')              '-------------------------------------------------'
    write(*,*)

    ! Print orbitals

    if(dump_orb) then
      write(*,'(A50)') '-----------------------------------------'
      write(*,'(A50)') 'qsUGW spin-up   orbital coefficients '
      write(*,'(A50)') '-----------------------------------------'
      call matout(nBas,nBas,c(:,:,1))
      write(*,*)
      write(*,'(A50)') '-----------------------------------------'
      write(*,'(A50)') 'qsUGW spin-down orbital coefficients '
      write(*,'(A50)') '-----------------------------------------'
      call matout(nBas,nBas,c(:,:,2))
      write(*,*)
    end if

  end if

end subroutine 
