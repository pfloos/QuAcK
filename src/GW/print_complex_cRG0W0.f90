subroutine print_complex_cRG0W0(nBas,nO,Re_eHF,Im_eHF,ENuc,ERHF,Re_SigC,Im_SigC,Re_Z,Im_Z,Re_eGW,Im_eGW,EcRPA,EcGM)

! Print one-electron energies and other stuff for G0W0

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas,nO
  double precision,intent(in)        :: ENuc
  complex*16,intent(in)              :: ERHF
  complex*16,intent(in)              :: EcRPA
  complex*16,intent(in)              :: EcGM
  double precision,intent(in)        :: Re_eHF(nBas)
  double precision,intent(in)        :: Im_eHF(nBas)
  double precision,intent(in)        :: Re_SigC(nBas)
  double precision,intent(in)        :: Im_SigC(nBas)
  double precision,intent(in)        :: Re_Z(nBas)
  double precision,intent(in)        :: Im_Z(nBas)
  double precision,intent(in)        :: Re_eGW(nBas)
  double precision,intent(in)        :: Im_eGW(nBas)

  integer                            :: p,index_homo,index_lumo
  double precision                   :: Re_eHOMO,Re_eLUMO,Im_eHOMO,Im_eLUMO,Re_Gap,Im_Gap
 
! HOMO and LUMO

  index_homo = maxloc(Re_eGW(1:nO),1)
  Re_eHOMO = Re_eGW(index_homo)
  Im_eHOMO = Im_eGW(index_homo)
  index_lumo = minloc(Re_eGW(nO+1:nBas),1) + nO
  Re_eLUMO = Re_eGW(index_lumo)
  Im_eLUMO = Im_eGW(index_lumo)
  Re_Gap = Re_eLUMO-Re_eHOMO
  Im_Gap = Im_eLUMO-Im_eHOMO

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)' cG0W0@RHF calculation '
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','Re(e_HF) (eV)','|','Re(Sig_GW) (eV)','|','Re(Z)','|','Re(e_GW) (eV)','|'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','Im(e_HF) (eV)','|','Im(Sig_GW) (eV)','|','Im(Z)','|','Im(e_GW) (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'
  do p=1,nBas
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',Re_eHF(p)*HaToeV,'|',Re_SigC(p)*HaToeV,'|',Re_Z(p),'|',Re_eGW(p)*HaToeV,'|'
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',Im_eHF(p)*HaToeV,'|',Im_SigC(p)*HaToeV,'|',Im_Z(p),'|',Im_eGW(p)*HaToeV,'|'
  write(*,*)'-------------------------------------------------------------------------------'
  end do
  write(*,*)
  write(*,*)'---------------------------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A5,F15.6,A3)') 'cG0W0@RHF HOMO      energy = ',Re_eHOMO*HaToeV,' + i*',Im_eHOMO*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A5,F15.6,A3)') 'cG0W0@RHF LUMO      energy = ',Re_eLUMO*HaToeV,' + i*',Im_eLUMO*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A5,F15.6,A3)') 'cG0W0@RHF HOMO-LUMO gap    = ',Re_Gap*HaToeV,' + i*',Im_Gap*HaToeV,' eV'
  write(*,*)'---------------------------------------------------------------------------------------------------'
!  write(*,'(2X,A60,F15.6,A3)') 'phRPA@cG0W0@RHF total       energy = ',ENuc + ERHF + EcRPA,' au'
!  write(*,'(2X,A60,F15.6,A3)') 'phRPA@cG0W0@RHF correlation energy = ',EcRPA,' au'
!  write(*,'(2X,A60,F15.6,A3)') '   GM@cG0W0@RHF total       energy = ',ENuc + ERHF + EcGM,' au'
!  write(*,'(2X,A60,F15.6,A3)') '   GM@cG0W0@RHF correlation energy = ',EcGM,' au'
!  write(*,*)'-------------------------------------------------------------------------------'
!  write(*,*)
!
end subroutine 
