subroutine R_print_parquet_1b(nOrb,nC,nO,nV,nR,eHF,SigC,eQP,Z,n_it_1b,err_1b,ENuc,EGHF,EcGM,Ec_eh,Ec_pp)

! Print one-electron energies and other stuff for G0F2

  implicit none
  include 'parameters.h'

  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: SigC(nOrb)
  double precision,intent(in)   :: eQP(nOrb)
  double precision,intent(in)   :: Z(nOrb)
  integer,intent(in)            :: n_it_1b
  double precision,intent(in)   :: err_1b
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EGHF
  double precision,intent(in)   :: EcGM
  double precision,intent(in)   :: Ec_eh(nspin)
  double precision,intent(in)   :: Ec_pp(nspin)

  integer                       :: i,a
  double precision              :: eHOMO
  double precision              :: eLUMO
  double precision              :: Gap

! HOMO and LUMO

  eHOMO = maxval(eQP(1:nO))
  eLUMO = minval(eQP(nO+1:nOrb))
  Gap = eLUMO - eHOMO
  
! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)' Parquet self-energy '
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','Sig_c (eV)','|','Z','|','e_QP (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'

  ! Occupied states

  do i=nC+1,nO
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',i,'|',eHF(i)*HaToeV,'|',SigC(i)*HaToeV,'|',Z(i),'|',eQP(i)*HaToeV,'|'
  end do

  ! Fermi level

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)'-------------------------------------------------------------------------------'

  ! Vacant states

  do a=nO+1,nOrb-nR
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',a,'|',eHF(a)*HaToeV,'|',SigC(a)*HaToeV,'|',Z(a),'|',eQP(a)*HaToeV,'|'
  end do

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,I15)')   'One-body iteration # ',n_it_1b
  write(*,'(2X,A60,F15.6)') 'One-body convergence ',err_1b
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'Parquet HOMO    energy = ',eHOMO*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'Parquet LUMO    energy = ',eLUMO*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'Parquet HOMO-LUMO gap  = ',Gap*HaToeV,' eV'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') '     Parquet total       energy = ',ENuc + EGHF + EcGM,' au'
  write(*,'(2X,A60,F15.6,A3)') '     Parquet correlation energy = ',EcGM,' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') '     eh-RPA  correlation energy = ',Ec_eh(1)+3d0*Ec_eh(2),' au'
  write(*,'(2X,A60,F15.6,A3)') '     pp-RPA  correlation energy = ',Ec_pp(1)+3d0*Ec_pp(2),' au'
  !write(*,'(2X,A60,F15.6,A3)') '(eh+pp)-RPA  correlation energy = ',Ec_pp,' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

end subroutine 
