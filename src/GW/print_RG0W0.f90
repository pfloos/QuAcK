subroutine print_RG0W0(nBas,nO,eHF,ENuc,ERHF,SigC,Z,eGW,EcRPA,EcGM)

! Print one-electron energies and other stuff for G0W0

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
  double precision,intent(in)        :: eGW(nBas)

  integer                            :: p
  double precision                   :: eHOMO,eLUMO,Gap

! HOMO and LUMO

  eHOMO = maxval(eGW(1:nO))
  eLUMO = minval(eGW(nO+1:nBas))
  Gap = eLUMO-eHOMO

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)' G0W0@RHF calculation '
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','Sig_GW (eV)','|','Z','|','e_GW (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'

  do p=1,nBas
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',eHF(p)*HaToeV,'|',SigC(p)*HaToeV,'|',Z(p),'|',eGW(p)*HaToeV,'|'
  end do

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'G0W0@RHF HOMO      energy = ',eHOMO*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'G0W0@RHF LUMO      energy = ',eLUMO*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'G0W0@RHF HOMO-LUMO gap    = ',Gap*HaToeV,' eV'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'phRPA@G0W0@RHF total       energy = ',ENuc + ERHF + EcRPA,' au'
  write(*,'(2X,A60,F15.6,A3)') 'phRPA@G0W0@RHF correlation energy = ',EcRPA,' au'
  write(*,'(2X,A60,F15.6,A3)') '   GM@G0W0@RHF total       energy = ',ENuc + ERHF + EcGM,' au'
  write(*,'(2X,A60,F15.6,A3)') '   GM@G0W0@RHF correlation energy = ',EcGM,' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

end subroutine 
