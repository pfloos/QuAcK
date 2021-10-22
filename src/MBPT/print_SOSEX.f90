subroutine print_SOSEX(nBas,nO,eHF,ENuc,ERHF,SigC,Z,eSOSEX,EcRPA,EcGM)

! Print one-electron energies and other stuff for SOSEX

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas,nO
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ERHF
  double precision,intent(in)        :: EcRPA
  double precision,intent(in)        :: EcGM
  double precision,intent(in)        :: eHF(nBas)
  double precision,intent(in)        :: SigC(nBas)
  double precision,intent(in)        :: Z(nBas)
  double precision,intent(in)        :: eSOSEX(nBas)

  integer                            :: p,HOMO,LUMO
  double precision                   :: Gap

! HOMO and LUMO

  HOMO = nO
  LUMO = HOMO + 1
  Gap = eSOSEX(LUMO)-eSOSEX(HOMO)

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)'  One-shot SOSEX calculation  '
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','Sig_c (eV)','|','Z','|','e_QP (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'

  do p=1,nBas
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',eHF(p)*HaToeV,'|',SigC(p)*HaToeV,'|',Z(p),'|',eSOSEX(p)*HaToeV,'|'
  enddo

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A30,F15.6,A3)') 'SOSEX HOMO      energy:',eSOSEX(HOMO)*HaToeV,' eV'
  write(*,'(2X,A30,F15.6,A3)') 'SOSEX LUMO      energy:',eSOSEX(LUMO)*HaToeV,' eV'
  write(*,'(2X,A30,F15.6,A3)') 'SOSEX HOMO-LUMO gap   :',Gap*HaToeV,' eV'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A30,F15.6,A3)') 'RPA@SOSEX total energy      :',ENuc + ERHF + EcRPA,' au'
  write(*,'(2X,A30,F15.6,A3)') 'RPA@SOSEX correlation energy:',EcRPA,' au'
  write(*,'(2X,A30,F15.6,A3)') 'GM@SOSEX  total energy      :',ENuc + ERHF + EcGM,' au'
  write(*,'(2X,A30,F15.6,A3)') 'GM@SOSEX  correlation energy:',EcGM,' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

end subroutine print_SOSEX


