subroutine print_UG0W0(nBas,nO,e,ENuc,EHF,SigC,Z,eGW,EcRPA)

! Print one-electron energies and other stuff for G0W0

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nO(nspin)
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: EHF
  double precision,intent(in)        :: EcRPA
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
  write(*,*)'  Unrestricted one-shot G0W0 calculation (eV)'
  write(*,*)'-------------------------------------------------------------------------------& 
              -------------------------------------------------'
  write(*,'(A1,A3,A1,2A15,A1,2A15,A1,2A15,A1,2A15,A1)') &
            '|','#','|','e_HF up','e_HF dw','|','Sig_c up','Sig_c dw','|', & 
                        'Z up','Z dw','|','e_QP up','e_QP dw','|'
  write(*,*)'-------------------------------------------------------------------------------& 
              -------------------------------------------------'

  do p=1,nBas
    write(*,'(A1,I3,A1,2F15.6,A1,2F15.6,A1,2F15.6,A1,2F15.6,A1)') &
    '|',p,'|',e(p,1)*HaToeV,e(p,2)*HaToeV,'|',SigC(p,1)*HaToeV,SigC(p,2)*HaToeV,'|', & 
              Z(p,1),Z(p,2),'|',eGW(p,1)*HaToeV,eGW(p,2)*HaToeV,'|'
  enddo

  write(*,*)'-------------------------------------------------------------------------------& 
              -------------------------------------------------'
  write(*,'(2X,A30,F15.6)') 'G0W0 HOMO      energy (eV):',HOMO*HaToeV
  write(*,'(2X,A30,F15.6)') 'G0W0 LUMO      energy (eV):',LUMO*HaToeV
  write(*,'(2X,A30,F15.6)') 'G0W0 HOMO-LUMO gap    (eV):',Gap*HaToeV
  write(*,*)'-------------------------------------------------------------------------------& 
              -------------------------------------------------'
  write(*,'(2X,A30,F15.6)') 'RPA@HF   total energy       =',ENuc + EHF + EcRPA
  write(*,'(2X,A30,F15.6)') 'RPA@HF   correlation energy =',EcRPA
  write(*,*)'-------------------------------------------------------------------------------& 
              -------------------------------------------------'
  write(*,*)

end subroutine print_UG0W0


