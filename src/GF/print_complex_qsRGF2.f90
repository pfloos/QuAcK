
! ---

subroutine print_complex_qsRGF2(nBas, nOrb, nO, nSCF, Conv, thresh, eHF, eGF, c, SigC, &
                       Z, ENuc, ET, EV,EW, EJ, EK, EcGM, EcRPA, EqsGF, dipole)

! Print useful information about qsRGF calculation

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
  complex*16,intent(in)              :: eGF(nOrb)
  complex*16,intent(in)              :: c(nBas,nOrb)
  complex*16,intent(in)              :: SigC(nOrb,nOrb)
  complex*16,intent(in)              :: Z(nOrb)
  complex*16,intent(in)              :: EqsGF
  complex*16,intent(in)              :: dipole(ncart)

! Local variables

  logical                            :: dump_orb = .false.
  integer                            :: p,ixyz,HOMO,LUMO
  complex*16                         :: Gap
  double precision,external          :: complex_trace_matrix

! Output variables

! HOMO and LUMO

  HOMO = maxloc(real(eGF(1:nO)),1)
  LUMO = minloc(real(eGF(nO+1:nBas)),1) + nO
  Gap = eGF(LUMO)-eGF(HOMO)

! Compute energies

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)' Self-consistent qsGF2 calculation'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','Re(e_HF (eV))','|','Re(Sig_GF) (eV)','|','Re(Z)','|','Re(e_GF) (eV)','|'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','Im(e_HF (eV))','|','Im(Sig_GF) (eV)','|','Im(Z)','|','Im(e_GF) (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'

  do p=1,nOrb
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
        '|',p,'|',real(eHF(p))*HaToeV,'|',real(SigC(p,p))*HaToeV,'|',real(Z(p)),'|',real(eGF(p))*HaToeV,'|'
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
        '|',p,'|',aimag(eHF(p))*HaToeV,'|',aimag(SigC(p,p))*HaToeV,'|',aimag(Z(p)),'|',aimag(eGF(p))*HaToeV,'|'
    write(*,*)'-------------------------------------------------------------------------------'
    if(p==nO) then
      write(*,*)'-------------------------------------------------------------------------------'
    end if

  end do

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A10,I3)')   'Iteration ',nSCF
  write(*,'(2X,A14,F15.5)')'Convergence = ',Conv
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'qsGF@RHF HOMO real energy = ',real(eGF(HOMO))*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'qsGF@RHF HOMO imag energy = ',aimag(eGF(HOMO))*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'qsGF@RHF LUMO real energy = ',real(eGF(LUMO))*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'qsGF@RHF LUMO imag energy = ',aimag(eGF(LUMO))*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'qsGF@RHF HOMO-LUMO gap    = ',real(Gap)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'qsGF@RHF HOMO-LUMO gap    = ',aimag(Gap)*HaToeV,' eV'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') '      qsGF@RHF total  real energy = ',ENuc + real(EqsGF),' au'
  write(*,'(2X,A60,F15.6,A3)') '      qsGF@RHF total  imag energy = ',aimag(EqsGF),' au'
  write(*,'(2X,A60,F15.6,A3)') '      qsGF@RHF exchange    energy = ',real(EK),' au'
  write(*,'(2X,A60,F15.6,A3)') '      qsGF@RHF exchange    energy = ',aimag(EK),' au'
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
    write(*,'(A33,1X,F16.10,A3)') ' Electronic   energy = ',real(EqsGF),' au'
    write(*,'(A33,1X,F16.10,A3)') ' Electronic   energy = ',aimag(EqsGF),' au'
    write(*,'(A33,1X,F16.10,A3)') ' Nuclear   repulsion = ',ENuc,' au'
    write(*,'(A33,1X,F16.10,A3)') ' qsRGF        energy = ',ENuc + real(EqsGF),' au'
    write(*,'(A33,1X,F16.10,A3)') ' qsRGF        energy = ',aimag(EqsGF),' au'
    write(*,'(A50)')           '---------------------------------------'
    write(*,*)
 
    if(dump_orb) then
      write(*,'(A50)') '---------------------------------------'
      write(*,'(A50)') ' Restricted qsGF orbital coefficients'
      write(*,'(A50)') '---------------------------------------'
      call complex_matout(nBas, nOrb, c)
      write(*,*)
    end if
    write(*,'(A50)') '---------------------------------------'
    write(*,'(A50)') ' Restricted qsGF orbital energies (au) '
    write(*,'(A50)') '---------------------------------------'
    call complex_vecout(nOrb, eGF)
    write(*,*)

  end if

end subroutine 
