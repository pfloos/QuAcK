subroutine print_complex_evRGF2(nBas,nO,nSCF,Conv,Re_eHF,Im_eHF,ENuc,ERHF,Re_SigC,Im_SigC,Re_Z,Im_Z,Re_eGF,Im_eGF,EcRPA,EcGM)

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
  double precision,intent(in)        :: Re_eGF(nBas)
  double precision,intent(in)        :: Im_eGF(nBas)

  integer                            :: p,index_homo,index_lumo
  double precision                   :: Re_eHOMO,Re_eLUMO,Im_eHOMO,Im_eLUMO,Re_Gap,Im_Gap
 
! HOMO and LUMO

  index_homo = maxloc(Re_eGF(1:nO),1)
  Re_eHOMO = Re_eGF(index_homo)
  Im_eHOMO = Im_eGF(index_homo)
  index_lumo = minloc(Re_eGF(nO+1:nBas),1) + nO
  Re_eLUMO = Re_eGF(index_lumo)
  Im_eLUMO = Im_eGF(index_lumo)
  Re_Gap = Re_eLUMO-Re_eHOMO
  Im_Gap = Im_eLUMO-Im_eHOMO

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)' Self-consistent evGF2 calculation'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','Re(e_HF) (eV)','|','Re(Sig_GF) (eV)','|','Re(Z)','|','Re(e_GF) (eV)','|'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','Im(e_HF) (eV)','|','Im(Sig_GF) (eV)','|','Im(Z)','|','Im(e_GF) (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'
  do p=1,nBas
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',Re_eHF(p)*HaToeV,'|',Re_SigC(p)*HaToeV,'|',Re_Z(p),'|',Re_eGF(p)*HaToeV,'|'
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',Im_eHF(p)*HaToeV,'|',Im_SigC(p)*HaToeV,'|',Im_Z(p),'|',Im_eGF(p)*HaToeV,'|'
    if(p==nO) then
      write(*,*)'-------------------------------------------------------------------------------'
    end if
  write(*,*)'-------------------------------------------------------------------------------'
  end do
  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A10,I3)')   'Iteration ',nSCF
  write(*,'(2X,A14,F15.9)')'Convergence = ',Conv
  write(*,*)
end subroutine 
