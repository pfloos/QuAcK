subroutine print_complex_evRGW(nBas,nO,nSCF,Conv,Re_eHF,Im_eHF,ENuc,ERHF,Re_SigC,Im_SigC,Re_Z,Im_Z,Re_eGW,Im_eGW,EcRPA,EcGM)

! Print one-electron energies and other stuff for G0W0

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas,nO,nSCF
  double precision,intent(in)        :: Conv
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
  if(nSCF < 10) then
    write(*,'(1X,A20,I1,A1,I1,A17)')' Self-consistent evG',nSCF,'W',nSCF,'@cRHF calculation'
  elseif(nSCF < 100) then
    write(*,'(1X,A20,I2,A1,I2,A17)')' Self-consistent evG',nSCF,'W',nSCF,'@cRHF calculation'
  else
    write(*,'(1X,A20,I3,A1,I3,A17)')' Self-consistent evG',nSCF,'W',nSCF,'@cRHF calculation'
  end if
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
    if(p==nO) then
      write(*,*)'-------------------------------------------------------------------------------'
    end if
  write(*,*)'-------------------------------------------------------------------------------'
  end do
  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A10,I3)')   'Iteration ',nSCF
  write(*,'(2X,A14,F15.5)')'Convergence = ',Conv
  write(*,*)
end subroutine 
