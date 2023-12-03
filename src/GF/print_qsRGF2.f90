subroutine print_qsRGF2(nBas,nO,nSCF,Conv,thresh,eHF,eGF,c,SigC,Z,ENuc,ET,EV,EJ,Ex,Ec,EqsGF,dipole)

! Print one-electron energies and other stuff for qsGF2

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nO
  integer,intent(in)                 :: nSCF
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: Conv
  double precision,intent(in)        :: thresh
  double precision,intent(in)        :: eHF(nBas)
  double precision,intent(in)        :: eGF(nBas)
  double precision,intent(in)        :: c(nBas)
  double precision,intent(in)        :: SigC(nBas,nBas)
  double precision,intent(in)        :: Z(nBas)
  double precision,intent(in)        :: ET
  double precision,intent(in)        :: EV
  double precision,intent(in)        :: EJ
  double precision,intent(in)        :: Ex
  double precision,intent(in)        :: Ec
  double precision,intent(in)        :: EqsGF
  double precision,intent(in)        :: dipole(ncart)

! Local variables

  integer                            :: q,ixyz,HOMO,LUMO
  double precision                   :: Gap
  double precision,external          :: trace_matrix

! Output variables

! HOMO and LUMO

  HOMO = nO
  LUMO = HOMO + 1
  Gap = eGF(LUMO)-eGF(HOMO)

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  if(nSCF < 10) then
    write(*,'(1X,A21,I1,A2,A12)')'  Self-consistent qsG',nSCF,'F2',' calculation'
  else
    write(*,'(1X,A21,I2,A2,A12)')'  Self-consistent qsG',nSCF,'F2',' calculation'
  end if
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','Sig_c (eV)','|','Z','|','e_QP (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'

  do q=1,nBas
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',q,'|',eHF(q)*HaToeV,'|',SigC(q,q)*HaToeV,'|',Z(q),'|',eGF(q)*HaToeV,'|'
  end do

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A10,I3)')   'Iteration ',nSCF
  write(*,'(2X,A14,F15.5)')'Convergence = ',Conv
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'qsGF2 HOMO      energy =',eGF(HOMO)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'qsGF2 LUMO      energy =',eGF(LUMO)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'qsGF2 HOMO-LUMO gap    =',Gap*HaToeV,' eV'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') '    qsGF2 total       energy =',ENuc + EqsGF,' au'
  write(*,'(2X,A60,F15.6,A3)') '    qsGF2 exchange    energy =',Ex,' au'
  write(*,'(2X,A60,F15.6,A3)') '    qsGF2 correlation energy =',Ec,' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

! Dump results for final iteration

  if(Conv < thresh) then

    write(*,*)
    write(*,'(A50)')           '---------------------------------------'
    write(*,'(A32)')           ' Summary              '
    write(*,'(A50)')           '---------------------------------------'
    write(*,'(A32,1X,F16.10,A3)') ' One-electron energy: ',ET + EV,' au'
    write(*,'(A32,1X,F16.10,A3)') ' Kinetic      energy: ',ET,' au'
    write(*,'(A32,1X,F16.10,A3)') ' Potential    energy: ',EV,' au'
    write(*,'(A50)')           '---------------------------------------'
    write(*,'(A32,1X,F16.10,A3)') ' Two-electron energy: ',EJ + Ex + Ec,' au'
    write(*,'(A32,1X,F16.10,A3)') ' Hartree      energy: ',EJ,' au'
    write(*,'(A32,1X,F16.10,A3)') ' Exchange     energy: ',Ex,' au'
    write(*,'(A32,1X,F16.10,A3)') ' Correlation  energy: ',Ec,' au'
    write(*,'(A50)')           '---------------------------------------'
    write(*,'(A32,1X,F16.10,A3)') ' Electronic   energy: ',EqsGF,' au'
    write(*,'(A32,1X,F16.10,A3)') ' Nuclear   repulsion: ',ENuc,' au'
    write(*,'(A32,1X,F16.10,A3)') ' qsGF2        energy: ',ENuc + EqsGF,' au'
    write(*,'(A50)')           '---------------------------------------'
    write(*,'(A35)')           ' Dipole moment (Debye)    '
    write(*,'(10X,4A10)')      'X','Y','Z','Tot.'
    write(*,'(10X,4F10.6)')    (dipole(ixyz)*auToD,ixyz=1,ncart),norm2(dipole)*auToD
    write(*,'(A50)')           '-----------------------------------------'
    write(*,*)
 
    write(*,'(A50)') '---------------------------------------'
    write(*,'(A32)') ' qsGF2 MO coefficients'
    write(*,'(A50)') '---------------------------------------'
    call matout(nBas,nBas,c)
    write(*,*)
    write(*,'(A50)') '---------------------------------------'
    write(*,'(A32)') ' qsGF2 MO energies'
    write(*,'(A50)') '---------------------------------------'
    call matout(nBas,1,eGF)
    write(*,*)

  end if


end subroutine 
