subroutine print_parquet_1b(nOrb,nO,eHF,SigC,eQP,Z,n_it_1b,err_1b,ENuc,ERHF,Ec)

! Print one-electron energies and other stuff for G0F2

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nOrb
  integer,intent(in)                 :: nO
  double precision,intent(in)        :: eHF(nOrb)
  double precision,intent(in)        :: SigC(nOrb)
  double precision,intent(in)        :: eQP(nOrb)
  double precision,intent(in)        :: Z(nOrb)
  integer,intent(in)                 :: n_it_1b
  double precision,intent(in)        :: err_1b
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
  Gap = eQP(LUMO) - eQP(HOMO)

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)' Parquet self-energy '
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','Sig_c (eV)','|','Z','|','e_QP (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'

  do p=1,nOrb
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',eHF(p)*HaToeV,'|',SigC(p)*HaToeV,'|',Z(p),'|',eQP(p)*HaToeV,'|'
  end do

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,I15)')   'One-body iteration # ',n_it_1b
  write(*,'(2X,A60,F15.6)') 'One-body convergence ',err_1b
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'Parquet HOMO    energy = ',eQP(HOMO)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'Parquet LUMO    energy = ',eQP(LUMO)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'Parquet HOMO-LUMO gap  = ',Gap*HaToeV,' eV'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'Parquet total       energy = ',ENuc + ERHF + Ec,' au'
  write(*,'(2X,A60,F15.6,A3)') 'Parquet correlation energy = ',Ec,' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

end subroutine 
