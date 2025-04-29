
! ---

subroutine print_complex_qsRGF2(nBas, nOrb, nO, nSCF, Conv, thresh, eHF, eGW, c, SigC, &
                       Z, ENuc, ET, EV,EW, EJ, EK, EcGM, EcRPA, EqsGW, dipole)

! Print useful information about qsRGW calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                 :: nBas, nOrb
  integer,intent(in)                 :: nO
  integer,intent(in)                 :: nSCF
  double precision,intent(in)        :: ENuc
  complex*16,intent(in)              :: ET
  complex*16,intent(in)              :: EV
  complex*16,intent(in)              :: EW
  complex*16,intent(in)              :: EJ
  complex*16,intent(in)              :: EK
  complex*16,intent(in)              :: EcGM
  complex*16,intent(in)              :: EcRPA
  double precision,intent(in)        :: Conv
  double precision,intent(in)        :: thresh
  complex*16,intent(in)              :: eHF(nOrb)
  complex*16,intent(in)              :: eGW(nOrb)
  complex*16,intent(in)              :: c(nBas,nOrb)
  complex*16,intent(in)              :: SigC(nOrb,nOrb)
  complex*16,intent(in)              :: Z(nOrb)
  complex*16,intent(in)              :: EqsGW
  complex*16,intent(in)              :: dipole(ncart)

! Local variables

  logical                            :: dump_orb = .false.
  integer                            :: p,ixyz,HOMO,LUMO
  complex*16                         :: Gap
  double precision,external          :: complex_trace_matrix

! Output variables

! HOMO and LUMO

  HOMO = maxloc(real(eGW(1:nO)),1)
  LUMO = minloc(real(eGW(nO+1:nBas)),1) + nO
  Gap = eGW(LUMO)-eGW(HOMO)

! Compute energies

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)' Self-consistent qsGF2 calculation'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','Re(e_HF (eV))','|','Re(Sig_GW) (eV)','|','Re(Z)','|','Re(e_GW) (eV)','|'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','Im(e_HF (eV))','|','Im(Sig_GW) (eV)','|','Im(Z)','|','Im(e_GW) (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'

  do p=1,nOrb
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
        '|',p,'|',real(eHF(p))*HaToeV,'|',real(SigC(p,p))*HaToeV,'|',real(Z(p)),'|',real(eGW(p))*HaToeV,'|'
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
        '|',p,'|',aimag(eHF(p))*HaToeV,'|',aimag(SigC(p,p))*HaToeV,'|',aimag(Z(p)),'|',aimag(eGW(p))*HaToeV,'|'
    write(*,*)'-------------------------------------------------------------------------------'
    if(p==nO) then
      write(*,*)'-------------------------------------------------------------------------------'
    end if

  end do

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A10,I3)')   'Iteration ',nSCF
  write(*,'(2X,A14,F15.5)')'Convergence = ',Conv
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'qsGW@RHF HOMO real energy = ',real(eGW(HOMO))*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'qsGW@RHF HOMO imag energy = ',aimag(eGW(HOMO))*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'qsGW@RHF LUMO real energy = ',real(eGW(LUMO))*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'qsGW@RHF LUMO imag energy = ',aimag(eGW(LUMO))*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'qsGW@RHF HOMO-LUMO gap    = ',real(Gap)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'qsGW@RHF HOMO-LUMO gap    = ',aimag(Gap)*HaToeV,' eV'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') '      qsGW@RHF total  real energy = ',ENuc + real(EqsGW),' au'
  write(*,'(2X,A60,F15.6,A3)') '      qsGW@RHF total  imag energy = ',aimag(EqsGW),' au'
  write(*,'(2X,A60,F15.6,A3)') '      qsGW@RHF exchange    energy = ',real(EK),' au'
  write(*,'(2X,A60,F15.6,A3)') '      qsGW@RHF exchange    energy = ',aimag(EK),' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

! Dump results for final iteration

  if(Conv < thresh) then

    write(*,*)
    write(*,'(A50)')           '---------------------------------------'
    write(*,'(A33)')           ' Summary               '
    write(*,'(A50)')           '---------------------------------------'
    write(*,'(A33,1X,F16.10,A3)') ' One-electron energy = ',real(ET) + real(EV) + real(EW),' au'
    write(*,'(A33,1X,F16.10,A3)') ' One-electron energy = ',aimag(ET) + aimag(EV) + aimag(EW),' au'
    write(*,'(A33,1X,F16.10,A3)') ' Kinetic      energy = ',real(ET),' au'
    write(*,'(A33,1X,F16.10,A3)') ' Kinetic      energy = ',aimag(ET),' au'
    write(*,'(A33,1X,F16.10,A3)') ' Potential    energy = ',real(EV),' au'
    write(*,'(A33,1X,F16.10,A3)') ' Potential    energy = ',aimag(EV),' au'
    write(*,'(A33,1X,F16.10,A3)') ' CAP          energy = ',real(EW),' au'
    write(*,'(A33,1X,F16.10,A3)') ' CAP          energy = ',aimag(EW),' au'
    write(*,'(A50)')           '---------------------------------------'
    write(*,'(A33,1X,F16.10,A3)') ' Two-electron energy = ',real(EJ + EK),' au'
    write(*,'(A33,1X,F16.10,A3)') ' Two-electron energy = ',aimag(EJ + EK),' au'
    write(*,'(A33,1X,F16.10,A3)') ' Hartree      energy = ',real(EJ),' au'
    write(*,'(A33,1X,F16.10,A3)') ' Hartree      energy = ',aimag(EJ),' au'
    write(*,'(A33,1X,F16.10,A3)') ' Exchange     energy = ',real(EK),' au'
    write(*,'(A33,1X,F16.10,A3)') ' Exchange     energy = ',aimag(EK),' au'
    write(*,'(A33,1X,F16.10,A3)') ' Correlation  energy = ',real(EcGM),' au'
    write(*,'(A33,1X,F16.10,A3)') ' Correlation  energy = ',aimag(EcGM),' au'
    write(*,'(A50)')           '---------------------------------------'
    write(*,'(A33,1X,F16.10,A3)') ' Electronic   energy = ',real(EqsGW),' au'
    write(*,'(A33,1X,F16.10,A3)') ' Electronic   energy = ',aimag(EqsGW),' au'
    write(*,'(A33,1X,F16.10,A3)') ' Nuclear   repulsion = ',ENuc,' au'
    write(*,'(A33,1X,F16.10,A3)') ' qsRGW        energy = ',ENuc + real(EqsGW),' au'
    write(*,'(A33,1X,F16.10,A3)') ' qsRGW        energy = ',aimag(EqsGW),' au'
    write(*,'(A50)')           '---------------------------------------'
    write(*,*)
 
    if(dump_orb) then
      write(*,'(A50)') '---------------------------------------'
      write(*,'(A50)') ' Restricted qsGW orbital coefficients'
      write(*,'(A50)') '---------------------------------------'
      call complex_matout(nBas, nOrb, c)
      write(*,*)
    end if
    write(*,'(A50)') '---------------------------------------'
    write(*,'(A50)') ' Restricted qsGW orbital energies (au) '
    write(*,'(A50)') '---------------------------------------'
    call complex_vecout(nOrb, eGW)
    write(*,*)

  end if

end subroutine 
