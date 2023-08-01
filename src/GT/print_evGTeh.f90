subroutine print_evGTeh(nBas,nO,nSCF,Conv,eHF,ENuc,ERHF,SigC,Z,eGT,EcRPA,EcGM)

! Print one-electron energies and other stuff for evGTeh

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas,nO,nSCF
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ERHF
  double precision,intent(in)        :: Conv
  double precision,intent(in)        :: eHF(nBas)
  double precision,intent(in)        :: SigC(nBas)
  double precision,intent(in)        :: Z(nBas)
  double precision,intent(in)        :: eGT(nBas)
  double precision,intent(in)        :: EcRPA
  double precision,intent(in)        :: EcGM

  integer                            :: p,HOMO,LUMO
  double precision                   :: Gap

! HOMO and LUMO

  HOMO = nO
  LUMO = HOMO + 1
  Gap = eGT(LUMO)-eGT(HOMO)

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  if(nSCF < 10) then
    write(*,'(1X,A21,I1,A3,I1,A12)')'  Self-consistent evG',nSCF,'Teh',nSCF,' calculation'
  else
    write(*,'(1X,A21,I2,A3,I2,A12)')'  Self-consistent evG',nSCF,'Teh',nSCF,' calculation'
  endif
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','Sig_GTeh (eV)','|','Z','|','e_GTeh (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'

  do p=1,nBas
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',eHF(p)*HaToeV,'|',SigC(p)*HaToeV,'|',Z(p),'|',eGT(p)*HaToeV,'|'
  enddo

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A10,I3)')   'Iteration ',nSCF
  write(*,'(2X,A14,F15.5)')'Convergence = ',Conv
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'evGTeh HOMO      energy =',eGT(HOMO)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'evGTeh LUMO      energy =',eGT(LUMO)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'evGTeh HOMO-LUMO gap    =',Gap*HaToeV,' eV'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'phRPA@evGTeh total energy       =',ENuc + ERHF + EcRPA,' au'
  write(*,'(2X,A60,F15.6,A3)') 'phRPA@evGTeh correlation energy =',EcRPA,' au'
  write(*,'(2X,A60,F15.6,A3)') '   GM@evGTeh total energy       =',ENuc + ERHF + EcGM,' au'
  write(*,'(2X,A60,F15.6,A3)') '   GM@evGTeh correlation energy =',EcGM,' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

end subroutine
