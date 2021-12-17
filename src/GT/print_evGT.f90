subroutine print_evGT(nBas,nO,nSCF,Conv,eHF,SigT,Z,eGT)

! Print one-electron energies and other stuff for evGT

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nO
  integer,intent(in)                 :: nSCF
  double precision,intent(in)        :: Conv
  double precision,intent(in)        :: eHF(nBas)
  double precision,intent(in)        :: SigT(nBas)
  double precision,intent(in)        :: Z(nBas)
  double precision,intent(in)        :: eGT(nBas)

  integer                            :: x,HOMO,LUMO
  double precision                   :: Gap

! HOMO and LUMO

  HOMO = nO
  LUMO = HOMO + 1
  Gap = eGT(LUMO)-eGT(HOMO)

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  if(nSCF < 10) then
    write(*,'(1X,A21,I1,A1,I1,A12)')'  Self-consistent evG',nSCF,'T',nSCF,' calculation'
  else
    write(*,'(1X,A21,I2,A1,I2,A12)')'  Self-consistent evG',nSCF,'T',nSCF,' calculation'
  endif
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','Sigma_T (eV)','|','Z','|','e_QP (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'

  do x=1,nBas
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',x,'|',eHF(x)*HaToeV,'|',SigT(x)*HaToeV,'|',Z(x),'|',eGT(x)*HaToeV,'|'
  enddo

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A10,I3)')   'Iteration ',nSCF
  write(*,'(2X,A14,F15.5)')'Convergence = ',Conv
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A30,F15.6)') 'evGT HOMO      energy (eV):',eGT(HOMO)*HaToeV
  write(*,'(2X,A30,F15.6)') 'evGT LUMO      energy (eV):',eGT(LUMO)*HaToeV
  write(*,'(2X,A30,F15.6)') 'evGT HOMO-LUMO gap    (eV):',Gap*HaToeV
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

end subroutine print_evGT
