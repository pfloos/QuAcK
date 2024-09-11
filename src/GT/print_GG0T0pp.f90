subroutine print_GG0T0pp(nBas,nO,eHF,ENuc,ERHF,SigT,Z,eGT,EcGM,EcRPA)

! Print one-electron energies and other stuff for G0T0

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nO
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ERHF
  double precision,intent(in)        :: EcGM
  double precision,intent(in)        :: EcRPA
  double precision,intent(in)        :: eHF(nBas)
  double precision,intent(in)        :: SigT(nBas)
  double precision,intent(in)        :: Z(nBas)
  double precision,intent(in)        :: eGT(nBas)

  integer                            :: p,HOMO,LUMO
  double precision                   :: Gap

! HOMO and LUMO

  HOMO = nO
  LUMO = HOMO + 1
  if(nBas >= LUMO) then

    Gap = eGT(LUMO) - eGT(HOMO)

  else

    Gap = 0d0

  end if

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)' G0T0pp@GHF calculation '
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','Sig_GTpp (eV)','|','Z','|','e_GTpp (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'

  do p=1,nBas
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',eHF(p)*HaToeV,'|',SigT(p)*HaToeV,'|',Z(p),'|',eGT(p)*HaToeV,'|'
  end do

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'G0T0pp@RHF HOMO      energy                 = ',eGT(HOMO)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'G0T0pp@RHF LUMO      energy                 = ',eGT(LUMO)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'G0T0pp@RHF HOMO-LUMO gap                    = ',Gap*HaToeV,' eV'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'ppRPA@G0T0pp@RHF correlation energy           = ',EcRPA,' au'
  write(*,'(2X,A60,F15.6,A3)') 'ppRPA@G0T0pp@RHF total       energy           = ',ENuc + ERHF + EcRPA,' au'
  write(*,'(2X,A60,F15.6,A3)') '   GM@G0T0pp@RHF correlation energy           = ',EcGM,' au'
  write(*,'(2X,A60,F15.6,A3)') '   GM@G0T0pp@RHF total       energy           = ',ENuc + ERHF + EcGM,' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

end subroutine 
