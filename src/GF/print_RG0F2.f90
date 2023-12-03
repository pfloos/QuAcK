subroutine print_RG0F2(nBas,nO,eHF,Sig,eGF,Z,ENuc,ERHF,Ec)

! Print one-electron energies and other stuff for G0F2

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nO
  double precision,intent(in)        :: eHF(nBas)
  double precision,intent(in)        :: Sig(nBas)
  double precision,intent(in)        :: eGF(nBas)
  double precision,intent(in)        :: Z(nBas)
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ERHF
  double precision,intent(in)        :: Ec

  integer                            :: p
  integer                            :: HOMO
  integer                            :: LUMO
  double precision                   :: Gap

! HOMO and LUMO

  HOMO = nO
  LUMO = HOMO + 1
  Gap = eGF(LUMO) - eGF(HOMO)

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)' One-shot G0F2 calculation'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','Sig_GF2 (eV)','|','Z','|','e_GF2 (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'

  do p=1,nBas
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',eHF(p)*HaToeV,'|',Sig(p)*HaToeV,'|',Z(p),'|',eGF(p)*HaToeV,'|'
  end do

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'G0F2  HOMO      energy =',eGF(HOMO)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'G0F2  LUMO      energy =',eGF(LUMO)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'G0F2  HOMO-LUMO gap    =',Gap*HaToeV,' eV'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'G0F2 total energy       =',ENuc + ERHF + Ec,' au'
  write(*,'(2X,A60,F15.6,A3)') 'G0F2 correlation energy =',Ec,' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

end subroutine 
