subroutine print_G0T0(nBas,nO,e,ENuc,ERHF,SigT,Z,eGW,EcRPA)

! Print one-electron energies and other stuff for G0T0

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas,nO
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ERHF
  double precision,intent(in)        :: EcRPA(nspin)
  double precision,intent(in)        :: e(nBas)
  double precision,intent(in)        :: SigT(nBas)
  double precision,intent(in)        :: Z(nBas)
  double precision,intent(in)        :: eGW(nBas)

  integer                            :: p,HOMO,LUMO
  double precision                   :: Gap

! HOMO and LUMO

  HOMO = nO
  LUMO = HOMO + 1
  Gap = eGW(LUMO)-eGW(HOMO)

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)'  One-shot G0T0 calculation (T-matrix self-energy)  '
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','Sigma_T (eV)','|','Z','|','e_QP (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'

  do p=1,nBas
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',e(p)*HaToeV,'|',SigT(p)*HaToeV,'|',Z(p),'|',eGW(p)*HaToeV,'|'
  enddo

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A40,F15.6)') 'G0T0 HOMO      energy (eV)            :',eGW(HOMO)*HaToeV
  write(*,'(2X,A40,F15.6)') 'G0T0 LUMO      energy (eV)            :',eGW(LUMO)*HaToeV
  write(*,'(2X,A40,F15.6)') 'G0T0 HOMO-LUMO gap    (eV)            :',Gap*HaToeV
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A40,F15.6)') 'RPA@G0T0 correlation energy (singlet) =',EcRPA(1)
  write(*,'(2X,A40,F15.6)') 'RPA@G0T0 correlation energy (triplet) =',EcRPA(2)
  write(*,'(2X,A40,F15.6)') 'RPA@G0T0 correlation energy           =',EcRPA(1) + EcRPA(2)
  write(*,'(2X,A40,F15.6)') 'RPA@G0T0 total energy                 =',ENuc + ERHF + EcRPA(1) + EcRPA(2)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

end subroutine print_G0T0


