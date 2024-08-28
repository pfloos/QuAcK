
! ---

subroutine print_qsRGTeh(nBas_AOs, nBas_MOs, nO, nSCF, Conv, thresh, eHF, eGT, c, SigC, &
                         Z, ENuc, ET, EV, EJ, Ex, EcGM, EcRPA, EqsGT, dipole)

! Print one-electron energies and other stuff for qsGTeh

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                 :: nBas_AOs, nBas_MOs
  integer,intent(in)                 :: nO
  integer,intent(in)                 :: nSCF
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ET
  double precision,intent(in)        :: EV
  double precision,intent(in)        :: EJ
  double precision,intent(in)        :: Ex
  double precision,intent(in)        :: EcGM
  double precision,intent(in)        :: EcRPA(nspin)
  double precision,intent(in)        :: Conv
  double precision,intent(in)        :: thresh
  double precision,intent(in)        :: eHF(nBas_MOs)
  double precision,intent(in)        :: eGT(nBas_MOs)
  double precision,intent(in)        :: c(nBas_AOs,nBas_MOs)
  double precision,intent(in)        :: SigC(nBas_MOs,nBas_MOs)
  double precision,intent(in)        :: Z(nBas_MOs)
  double precision,intent(in)        :: EqsGT
  double precision,intent(in)        :: dipole(ncart)

! Local variables

  logical                            :: dump_orb = .false.

  integer                            :: p,ixyz,HOMO,LUMO
  double precision                   :: Gap
  double precision,external          :: trace_matrix

! Output variables

! HOMO and LUMO

  HOMO = nO
  LUMO = HOMO + 1
  Gap = eGT(LUMO)-eGT(HOMO)

! Compute energies

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  if(nSCF < 10) then
    write(*,'(1X,A21,I1,A3,I1,A12)')'  Self-consistent qsG',nSCF,'Teh',nSCF,' calculation'
  else
    write(*,'(1X,A21,I2,A3,I2,A12)')'  Self-consistent qsG',nSCF,'Teh',nSCF,' calculation'
  end if
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','Sig_GTeh (eV)','|','Z','|','e_GTeh (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'

  do p=1,nBas_MOs
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',eHF(p)*HaToeV,'|',SigC(p,p)*HaToeV,'|',Z(p),'|',eGT(p)*HaToeV,'|'
  end do

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A10,I3)')   'Iteration ',nSCF
  write(*,'(2X,A14,F15.5)')'Convergence = ',Conv
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'qsRGTeh HOMO      energy =',eGT(HOMO)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'qsRGTeh LUMO      energy =',eGT(LUMO)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'qsRGTeh HOMO-LUMO gap    =',Gap*HaToeV,' eV'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') '      qsRGTeh total       energy =',ENuc + EqsGT,' au'
  write(*,'(2X,A60,F15.6,A3)') '      qsRGTeh exchange    energy =',Ex,' au'
  write(*,'(2X,A60,F15.6,A3)') '   GM@qsRGTeh correlation energy =',EcGM,' au'
  write(*,'(2X,A60,F15.6,A3)') 'ppRPA@qsRGTeh correlation energy =',sum(EcRPA),' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

! Dump results for final iteration

  if(Conv < thresh) then

    write(*,*)
    write(*,'(A50)')           '---------------------------------------'
    write(*,'(A32)')           ' Summary              '
    write(*,'(A50)')           '---------------------------------------'
    write(*,'(A32,1X,F16.10,A3)') ' One-electron energy = ',ET + EV,' au'
    write(*,'(A32,1X,F16.10,A3)') ' Kinetic      energy = ',ET,' au'
    write(*,'(A32,1X,F16.10,A3)') ' Potential    energy = ',EV,' au'
    write(*,'(A50)')           '---------------------------------------'
    write(*,'(A32,1X,F16.10,A3)') ' Two-electron energy = ',EJ + Ex,' au'
    write(*,'(A32,1X,F16.10,A3)') ' Hartree      energy = ',EJ,' au'
    write(*,'(A32,1X,F16.10,A3)') ' Exchange     energy = ',Ex,' au'
    write(*,'(A32,1X,F16.10,A3)') ' Correlation  energy = ',EcGM,' au'
    write(*,'(A50)')           '---------------------------------------'
    write(*,'(A32,1X,F16.10,A3)') ' Electronic   energy = ',EqsGT,' au'
    write(*,'(A32,1X,F16.10,A3)') ' Nuclear   repulsion = ',ENuc,' au'
    write(*,'(A32,1X,F16.10,A3)') ' qsGTeh       energy = ',ENuc + EqsGT,' au'
    write(*,'(A50)')           '---------------------------------------'
    write(*,'(A35)')           ' Dipole moment (Debye)    '
    write(*,'(10X,4A10)')      'X','Y','Z','Tot.'
    write(*,'(10X,4F10.6)')    (dipole(ixyz)*auToD,ixyz=1,ncart),norm2(dipole)*auToD
    write(*,'(A50)')           '-----------------------------------------'
    write(*,*)
 
    if(dump_orb) then
      write(*,'(A50)') '---------------------------------------'
      write(*,'(A32)') ' qsGTeh MO coefficients'
      write(*,'(A50)') '---------------------------------------'
      call matout(nBas_AOs, nBas_MOs, c)
      write(*,*)
    end if
    write(*,'(A50)') '---------------------------------------'
    write(*,'(A32)') ' qsGTeh MO energies'
    write(*,'(A50)') '---------------------------------------'
    call vecout(nBas_MOs, eGT)
    write(*,*)

  end if

end subroutine 
