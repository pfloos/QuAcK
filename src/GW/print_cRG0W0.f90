subroutine print_cRG0W0(nBas,nO,eHF,ENuc,ERHF,Re_SigC,Im_SigC,Re_Z,Im_Z,Re_eGW,Im_eGW,EcRPA,EcGM,CAP)

! Print one-electron energies and other stuff for G0W0

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas,nO
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ERHF
  double precision,intent(in)        :: EcRPA
  double precision,intent(in)        :: EcGM
  double precision,intent(in)        :: eHF(nBas)
  double precision,intent(in)        :: Re_SigC(nBas)
  double precision,intent(in)        :: Im_SigC(nBas)
  double precision,intent(in)        :: Re_Z(nBas)
  double precision,intent(in)        :: Im_Z(nBas)
  double precision,intent(in)        :: Re_eGW(nBas)
  double precision,intent(in)        :: Im_eGW(nBas)
  double precision,intent(in)        :: CAP(nBas,nBas)

  integer                            :: p
  double precision                   :: eHOMO,eLUMO,Gap

! HOMO and LUMO

  eHOMO = maxval(Re_eGW(1:nO))
  eLUMO = minval(Re_eGW(nO+1:nBas))
  Gap = eLUMO-eHOMO

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)' G0W0@RHF calculation '
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A11,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','Re(Sig_GW) (eV)','|','Re(Z)','|','Re(e_GW) (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'
  do p=1,nBas
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',eHF(p)*HaToeV,'|',Re_SigC(p)*HaToeV,'|',Re_Z(p),'|',Re_eGW(p)*HaToeV,'|'
  end do
  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A11,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','CAP (eV)','|','Im(Sig_GW) (eV)','|','Im(Z)','|','Im(e_GW) (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'
  do p=1,nBas
    write(*,'(1X,A11,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',CAP(p,p)*HaToeV,'|',Im_SigC(p)*HaToeV,'|',Im_Z(p),'|',Im_eGW(p)*HaToeV,'|'
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
