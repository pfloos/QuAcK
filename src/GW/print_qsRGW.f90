subroutine print_qsRGW(nBas,nO,nSCF,Conv,thresh,eHF,eGW,c,SigC,Z,ENuc,ET,EV,EJ,EK,EcGM,EcRPA,EqsGW,dipole)

! Print useful information about qsRGW calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nO
  integer,intent(in)                 :: nSCF
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ET
  double precision,intent(in)        :: EV
  double precision,intent(in)        :: EJ
  double precision,intent(in)        :: EK
  double precision,intent(in)        :: EcGM
  double precision,intent(in)        :: EcRPA
  double precision,intent(in)        :: Conv
  double precision,intent(in)        :: thresh
  double precision,intent(in)        :: eHF(nBas)
  double precision,intent(in)        :: eGW(nBas)
  double precision,intent(in)        :: c(nBas)
  double precision,intent(in)        :: SigC(nBas,nBas)
  double precision,intent(in)        :: Z(nBas)
  double precision,intent(in)        :: EqsGW
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
  Gap = eGW(LUMO)-eGW(HOMO)

! Compute energies

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  if(nSCF < 10) then
    write(*,'(1X,A20,I1,A1,I1,A16)')' Self-consistent qsG',nSCF,'W',nSCF,'@RHF calculation'
  elseif(nSCF < 100) then
    write(*,'(1X,A20,I2,A1,I2,A16)')' Self-consistent qsG',nSCF,'W',nSCF,'@RHF calculation'
  else
    write(*,'(1X,A20,I3,A1,I3,A16)')' Self-consistent qsG',nSCF,'W',nSCF,'@RHF calculation'
  end if
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','Sig_GW (eV)','|','Z','|','e_GW (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'

  do p=1,nBas
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',eHF(p)*HaToeV,'|',SigC(p,p)*HaToeV,'|',Z(p),'|',eGW(p)*HaToeV,'|'
  enddo

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A10,I3)')   'Iteration ',nSCF
  write(*,'(2X,A14,F15.5)')'Convergence = ',Conv
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'qsGW@RHF HOMO      energy = ',eGW(HOMO)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'qsGW@RHF LUMO      energy = ',eGW(LUMO)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'qsGW@RHF HOMO-LUMO gap    = ',Gap*HaToeV,' eV'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') '      qsGW@RHF total       energy = ',ENuc + EqsGW,' au'
  write(*,'(2X,A60,F15.6,A3)') '      qsGW@RHF exchange    energy = ',EK,' au'
  write(*,'(2X,A60,F15.6,A3)') '   GM@qsGW@RHF correlation energy = ',EcGM,' au'
  write(*,'(2X,A60,F15.6,A3)') 'phRPA@qsGW@RHF correlation energy = ',EcRPA,' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

! Dump results for final iteration

  if(Conv < thresh) then

    write(*,*)
    write(*,'(A50)')           '---------------------------------------'
    write(*,'(A33)')           ' Summary               '
    write(*,'(A50)')           '---------------------------------------'
    write(*,'(A33,1X,F16.10,A3)') ' One-electron energy = ',ET + EV,' au'
    write(*,'(A33,1X,F16.10,A3)') ' Kinetic      energy = ',ET,' au'
    write(*,'(A33,1X,F16.10,A3)') ' Potential    energy = ',EV,' au'
    write(*,'(A50)')           '---------------------------------------'
    write(*,'(A33,1X,F16.10,A3)') ' Two-electron energy = ',EJ + EK,' au'
    write(*,'(A33,1X,F16.10,A3)') ' Hartree      energy = ',EJ,' au'
    write(*,'(A33,1X,F16.10,A3)') ' Exchange     energy = ',EK,' au'
    write(*,'(A33,1X,F16.10,A3)') ' Correlation  energy = ',EcGM,' au'
    write(*,'(A50)')           '---------------------------------------'
    write(*,'(A33,1X,F16.10,A3)') ' Electronic   energy = ',EqsGW,' au'
    write(*,'(A33,1X,F16.10,A3)') ' Nuclear   repulsion = ',ENuc,' au'
    write(*,'(A33,1X,F16.10,A3)') ' qsRGW        energy = ',ENuc + EqsGW,' au'
    write(*,'(A50)')           '---------------------------------------'
    write(*,'(A36)')           ' Dipole moment (Debye)    '
    write(*,'(10X,4A10)')      'X','Y','Z','Tot.'
    write(*,'(10X,4F10.4)')    (dipole(ixyz)*auToD,ixyz=1,ncart),norm2(dipole)*auToD
    write(*,'(A50)')           '---------------------------------------'
    write(*,*)
 
    if(dump_orb) then
      write(*,'(A50)') '---------------------------------------'
      write(*,'(A50)') ' Restricted qsGW orbital coefficients'
      write(*,'(A50)') '---------------------------------------'
      call matout(nBas,nBas,c)
      write(*,*)
    end if
    write(*,'(A50)') '---------------------------------------'
    write(*,'(A50)') ' Restricted qsGW orbital energies (au) '
    write(*,'(A50)') '---------------------------------------'
    call vecout(nBas,eGW)
    write(*,*)

  endif

end subroutine 
