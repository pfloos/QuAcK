subroutine print_evGW(nBas,nO,nSCF,Conv,e,ENuc,EHF,SigmaC,Z,eGW,EcRPA,EcGM)

! Print one-electron energies and other stuff for evGW

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas,nO,nSCF
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: EHF
  double precision,intent(in)        :: EcRPA
  double precision,intent(in)        :: EcGM
  double precision,intent(in)        :: Conv,e(nBas),SigmaC(nBas),Z(nBas),eGW(nBas)

  integer                            :: x,HOMO,LUMO
  double precision                   :: Gap

! HOMO and LUMO

  HOMO = nO
  LUMO = HOMO + 1
  Gap = eGW(LUMO)-eGW(HOMO)

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  if(nSCF < 10) then
    write(*,'(1X,A21,I1,A1,I1,A12)')'  Self-consistent evG',nSCF,'W',nSCF,' calculation'
  else
    write(*,'(1X,A21,I2,A1,I2,A12)')'  Self-consistent evG',nSCF,'W',nSCF,' calculation'
  endif
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','Sigma_c (eV)','|','Z','|','e_QP (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'

  do x=1,nBas
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',x,'|',e(x)*HaToeV,'|',SigmaC(x)*HaToeV,'|',Z(x),'|',eGW(x)*HaToeV,'|'
  enddo

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A10,I3)')   'Iteration ',nSCF
  write(*,'(2X,A14,F15.5)')'Convergence = ',Conv
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A30,F15.6)') 'evGW HOMO      energy (eV):',eGW(HOMO)*HaToeV
  write(*,'(2X,A30,F15.6)') 'evGW LUMO      energy (eV):',eGW(LUMO)*HaToeV
  write(*,'(2X,A30,F15.6)') 'evGW HOMO-LUMO gap    (eV):',Gap*HaToeV
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A30,F15.6)') 'RPA@evGW total energy       =',ENuc + EHF + EcRPA
  write(*,'(2X,A30,F15.6)') 'RPA@evGW correlation energy =',EcRPA
  write(*,'(2X,A30,F15.6)') 'GM@evGW  total energy       =',ENuc + EHF + EcGM
  write(*,'(2X,A30,F15.6)') 'GM@evGW  correlation energy =',EcGM
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

end subroutine print_evGW
