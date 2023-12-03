subroutine print_evUGW(nBas,nO,nSCF,Conv,eHF,ENuc,EUHF,SigC,Z,eGW,EcRPA,EcGM)

! Print one-electron energies and other stuff for evGW

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nO(nspin)
  integer,intent(in)                 :: nSCF
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: EUHF
  double precision,intent(in)        :: EcRPA
  double precision,intent(in)        :: EcGM(nspin)
  double precision,intent(in)        :: Conv
  double precision,intent(in)        :: eHF(nBas,nspin)
  double precision,intent(in)        :: SigC(nBas,nspin)
  double precision,intent(in)        :: Z(nBas,nspin)
  double precision,intent(in)        :: eGW(nBas,nspin)

  integer                            :: p
  integer                            :: ispin
  double precision                   :: HOMO(nspin)
  double precision                   :: LUMO(nspin)
  double precision                   :: Gap(nspin)

! HOMO and LUMO

  do ispin=1,nspin
    if(nO(ispin) > 0) then
      HOMO(ispin) = eGW(nO(ispin),ispin)
      LUMO(ispin) = eGW(nO(ispin)+1,ispin)
      Gap(ispin)  = LUMO(ispin) - HOMO(ispin)
    else
      HOMO(ispin) = 0d0
      LUMO(ispin) = eGW(1,ispin)
      Gap(ispin)  = 0d0
    end if
  end do

! Dump results

  write(*,*)'----------------------------------------------------------------'// &
            '----------------------------------------------------------------'
  if(nSCF < 10) then
    write(*,'(1X,A20,I1,A1,I1,A16)')' Self-consistent evG',nSCF,'W',nSCF,'@UHF calculation'
  elseif(nSCF < 100) then
    write(*,'(1X,A20,I2,A1,I2,A16)')' Self-consistent evG',nSCF,'W',nSCF,'@UHF calculation'
  else
    write(*,'(1X,A20,I3,A1,I3,A16)')' Self-consistent evG',nSCF,'W',nSCF,'@UHF calculation'
  end if
  write(*,*)'----------------------------------------------------------------'// &
            '----------------------------------------------------------------'
  write(*,'(A1,A3,A1,A30,A1,A30,A1,A30,A1,A30,A1)') &
            '|',' ','|','e_HF (eV)         ','|','Sig_GW (eV)        ','|','Z             ','|','e_GW (eV)         ','|'
  write(*,'(A1,A3,A1,2A15,A1,2A15,A1,2A15,A1,2A15,A1)') &
            '|','#','|','up     ','dw     ','|','up     ','dw     ','|','up     ','dw     ','|','up     ','dw     ','|'
  write(*,*)'----------------------------------------------------------------'// &
            '----------------------------------------------------------------'

  do p=1,nBas
    write(*,'(A1,I3,A1,2F15.6,A1,2F15.6,A1,2F15.6,A1,2F15.6,A1)') &
    '|',p,'|',eHF(p,1)*HaToeV,eHF(p,2)*HaToeV,'|',SigC(p,1)*HaToeV,SigC(p,2)*HaToeV,'|', &
              Z(p,1),Z(p,2),'|',eGW(p,1)*HaToeV,eGW(p,2)*HaToeV,'|'
  end do

  write(*,*)'----------------------------------------------------------------'// &
            '----------------------------------------------------------------'
  write(*,'(2X,A10,I3)')   'Iteration ',nSCF
  write(*,'(2X,A14,F15.5)')'Convergence = ',Conv
  write(*,*)'----------------------------------------------------------------'// &
            '----------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'evGW@UHF HOMO      energy = ',maxval(HOMO)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'evGW@UHF LUMO      energy = ',minval(LUMO)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'evGW@UHF HOMO-LUMO gap    = ',(minval(LUMO)-maxval(HOMO))*HaToeV,' eV'
  write(*,*)'----------------------------------------------------------------'// &
            '----------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'RPA@evG@UHFW total       energy = ',ENuc + EUHF + EcRPA,' au'
  write(*,'(2X,A60,F15.6,A3)') 'RPA@evG@UHFW correlation energy = ',EcRPA,' au'
  write(*,'(2X,A60,F15.6,A3)') ' GM@evG@UHFW total       energy = ',ENuc + EUHF + sum(EcGM),' au'
  write(*,'(2X,A60,F15.6,A3)') ' GM@evG@UHFW correlation energy = ',sum(EcGM),' au'
  write(*,*)'----------------------------------------------------------------'// &
            '----------------------------------------------------------------'
  write(*,*)

end subroutine 
