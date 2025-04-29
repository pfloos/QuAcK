subroutine print_cRG0F2(nBas,nO,eHF,e_cap,Re_Sig,Im_Sig,Re_eGF,Im_eGF,Re_Z,Im_Z,ENuc,ERHF,Ec)

! Print one-electron energies and other stuff for G0F2

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nO
  double precision,intent(in)        :: eHF(nBas)
  double precision,intent(in)        :: e_cap(nBas)
  double precision,intent(in)        :: Re_Sig(nBas)
  double precision,intent(in)        :: Im_Sig(nBas)
  double precision,intent(in)        :: Re_eGF(nBas)
  double precision,intent(in)        :: Im_eGF(nBas)
  double precision,intent(in)        :: Re_Z(nBas)
  double precision,intent(in)        :: Im_Z(nBas)
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ERHF
  double precision,intent(in)        :: Ec

  integer                            :: p
  integer                            :: HOMO
  integer                            :: LUMO
  double precision                   :: Gap



! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)' One-shot G0F2 calculation'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','Re(Sig_GF2) (eV)','|','Re(Z)','|','Re(e_GF2) (eV)','|'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','CAP(p,p) (eV)','|','Im(Sig_GF2) (eV)','|','Im(Z)','|','Im(e_GF2) (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'

  do p=1,nBas
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',eHF(p)*HaToeV,'|',Re_Sig(p)*HaToeV,'|',Re_Z(p),'|',Re_eGF(p)*HaToeV,'|'
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',e_cap(p)*HaToeV,'|',Im_Sig(p)*HaToeV,'|',Im_Z(p),'|',Im_eGF(p)*HaToeV,'|'
  write(*,*)'-------------------------------------------------------------------------------'
  end do
end subroutine 
