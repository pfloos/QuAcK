subroutine print_evRGW(nBas,nO,nSCF,Conv,eHF,ENuc,ERHF,SigC,Z,eGW,EcRPA,EcGM)

! Print one-electron energies and other stuff for evGW

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas,nO,nSCF
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ERHF
  double precision,intent(in)        :: Conv
  double precision,intent(in)        :: eHF(nBas)
  double precision,intent(in)        :: SigC(nBas)
  double precision,intent(in)        :: Z(nBas)
  double precision,intent(in)        :: eGW(nBas)
  double precision,intent(in)        :: EcRPA
  double precision,intent(in)        :: EcGM

  integer                            :: p,HOMO,LUMO
  double precision                   :: Gap

! HOMO and LUMO

  HOMO = nO
  LUMO = HOMO + 1
  Gap = eGW(LUMO)-eGW(HOMO)

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  if(nSCF < 10) then
    write(*,'(1X,A20,I1,A1,I1,A16)')' Self-consistent evG',nSCF,'W',nSCF,'@RHF calculation'
  elseif(nSCF < 100) then
    write(*,'(1X,A20,I2,A1,I2,A16)')' Self-consistent evG',nSCF,'W',nSCF,'@RHF calculation'
  else
    write(*,'(1X,A20,I3,A1,I3,A16)')' Self-consistent evG',nSCF,'W',nSCF,'@RHF calculation'
  end if
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','Sig_GW (eV)','|','Z','|','e_GW (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'

  do p=1,nBas
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',eHF(p)*HaToeV,'|',SigC(p)*HaToeV,'|',Z(p),'|',eGW(p)*HaToeV,'|'
  end do

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A10,I3)')   'Iteration ',nSCF
  write(*,'(2X,A14,F15.5)')'Convergence = ',Conv
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'evGW@RHF HOMO      energy = ',eGW(HOMO)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'evGW@RHF LUMO      energy = ',eGW(LUMO)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'evGW@RHF HOMO-LUMO gap    = ',Gap*HaToeV,' eV'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'phRPA@evGW@RHF total       energy = ',ENuc + ERHF + EcRPA,' au'
  write(*,'(2X,A60,F15.6,A3)') 'phRPA@evGW@RHF correlation energy = ',EcRPA,' au'
  write(*,'(2X,A60,F15.6,A3)') '   GM@evGW@RHF total       energy = ',ENuc + ERHF + EcGM,' au'
  write(*,'(2X,A60,F15.6,A3)') '   GM@evGW@RHF correlation energy = ',EcGM,' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

end subroutine 
