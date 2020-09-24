subroutine print_evUGW(nBas,nO,nSCF,Conv,e,ENuc,EHF,SigC,Z,eGW,EcRPA)

! Print one-electron energies and other stuff for evGW

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nO(nspin)
  integer,intent(in)                 :: nSCF
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: EHF
  double precision,intent(in)        :: EcRPA
  double precision,intent(in)        :: Conv
  double precision,intent(in)        :: e(nBas,nspin)
  double precision,intent(in)        :: SigC(nBas,nspin)
  double precision,intent(in)        :: Z(nBas,nspin)
  double precision,intent(in)        :: eGW(nBas,nspin)

  integer                            :: p
  double precision                   :: HOMO
  double precision                   :: LUMO
  double precision                   :: Gap

! HOMO and LUMO

  HOMO = max(eGW(nO(1),1),eGW(nO(2),2))
  LUMO = min(eGW(nO(1)+1,1),eGW(nO(2)+1,2))
  Gap = LUMO - HOMO

! Dump results

  write(*,*)'-------------------------------------------------------------------------------& 
              -------------------------------------------------'
  if(nSCF < 10) then
    write(*,'(1X,A21,I1,A1,I1,A12)')'  Self-consistent evG',nSCF,'W',nSCF,' calculation'
  else
    write(*,'(1X,A21,I2,A1,I2,A12)')'  Self-consistent evG',nSCF,'W',nSCF,' calculation'
  endif
  write(*,*)'-------------------------------------------------------------------------------& 
              -------------------------------------------------'
  write(*,'(A1,A3,A1,A30,A1,A30,A1,A30,A1,A30,A1)') &
            '|',' ','|','e_HF          ','|','Sig_c          ','|','Z          ','|','e_QP          ','|'
  write(*,'(A1,A3,A1,2A15,A1,2A15,A1,2A15,A1,2A15,A1)') &
            '|','#','|','up     ','dw     ','|','up     ','dw     ','|','up     ','dw     ','|','up     ','dw     ','|'
  write(*,*)'-------------------------------------------------------------------------------& 
              -------------------------------------------------'

  do p=1,nBas
    write(*,'(A1,I3,A1,2F15.6,A1,2F15.6,A1,2F15.6,A1,2F15.6,A1)') &
    '|',p,'|',e(p,1)*HaToeV,e(p,2)*HaToeV,'|',SigC(p,1)*HaToeV,SigC(p,2)*HaToeV,'|', &
              Z(p,1),Z(p,2),'|',eGW(p,1)*HaToeV,eGW(p,2)*HaToeV,'|'
  enddo

  write(*,*)'-------------------------------------------------------------------------------& 
              -------------------------------------------------'
  write(*,'(2X,A10,I3)')   'Iteration ',nSCF
  write(*,'(2X,A14,F15.5)')'Convergence = ',Conv
  write(*,*)'-------------------------------------------------------------------------------& 
              -------------------------------------------------'
  write(*,'(2X,A30,F15.6)') 'evGW HOMO      energy (eV):',HOMO*HaToeV
  write(*,'(2X,A30,F15.6)') 'evGW LUMO      energy (eV):',LUMO*HaToeV
  write(*,'(2X,A30,F15.6)') 'evGW HOMO-LUMO gap    (eV):',Gap*HaToeV
  write(*,*)'-------------------------------------------------------------------------------& 
              -------------------------------------------------'
  write(*,'(2X,A30,F15.6)') 'RPA@evGW total energy       =',ENuc + EHF + EcRPA
  write(*,'(2X,A30,F15.6)') 'RPA@evGW correlation energy =',EcRPA
  write(*,*)'-------------------------------------------------------------------------------& 
              -------------------------------------------------'
  write(*,*)

end subroutine print_evUGW
