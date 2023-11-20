subroutine print_evRGTpp(nBas,nO,nSCF,Conv,eHF,ENuc,ERHF,SigT,Z,eGT,EcGM,EcRPA)

! Print one-electron energies and other stuff for evGT

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nO
  integer,intent(in)                 :: nSCF
  double precision,intent(in)        :: Conv
  double precision,intent(in)        :: eHF(nBas)
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ERHF
  double precision,intent(in)        :: EcGM
  double precision,intent(in)        :: EcRPA(nspin)
  double precision,intent(in)        :: SigT(nBas)
  double precision,intent(in)        :: Z(nBas)
  double precision,intent(in)        :: eGT(nBas)

  integer                            :: p,HOMO,LUMO
  double precision                   :: Gap

! HOMO and LUMO

  HOMO = nO
  LUMO = HOMO + 1
  Gap = eGT(LUMO)-eGT(HOMO)

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  if(nSCF < 10) then
    write(*,'(1X,A21,I1,A3,I1,A12)')'  Self-consistent evG',nSCF,'Tpp',nSCF,' calculation'
  else
    write(*,'(1X,A21,I2,A3,I2,A12)')'  Self-consistent evG',nSCF,'Tpp',nSCF,' calculation'
  endif
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','Sig_GTpp (eV)','|','Z','|','e_GTpp (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'

  do p=1,nBas
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',eHF(p)*HaToeV,'|',SigT(p)*HaToeV,'|',Z(p),'|',eGT(p)*HaToeV,'|'
  enddo

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A10,I3)')   'Iteration ',nSCF
  write(*,'(2X,A14,F15.5)')'Convergence = ',Conv
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'evRGTpp HOMO      energy =',eGT(HOMO)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'evRGTpp LUMO      energy =',eGT(LUMO)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'evRGTpp HOMO-LUMO gap    =',Gap*HaToeV,' eV'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'ppRPA@evRGTpp correlation energy (singlet) =',EcRPA(1),' au'
  write(*,'(2X,A60,F15.6,A3)') 'ppRPA@evRGTpp correlation energy (triplet) =',EcRPA(2),' au'
  write(*,'(2X,A60,F15.6,A3)') 'ppRPA@evRGTpp correlation energy           =',sum(EcRPA),' au'
  write(*,'(2X,A60,F15.6,A3)') 'ppRPA@evRGTpp total energy                 =',ENuc + ERHF + sum(EcRPA),' au'
  write(*,'(2X,A60,F15.6,A3)') '   GM@evRGTpp correlation energy           =',EcGM,' au'
  write(*,'(2X,A60,F15.6,A3)') '   GM@evRGTpp total energy                 =',ENuc + ERHF + EcGM,' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

end subroutine 
