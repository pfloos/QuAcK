subroutine print_R_2SOSEX(nOrb,nC,nO,nV,nR,eHF,ENuc,ERHF,SigC,Z,eQP,EcRPA,EcGM)

! Print one-electron energies and other stuff for 2SOSEX-psd

  implicit none
  include 'parameters.h'

  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: EcRPA
  double precision,intent(in)   :: EcGM
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: SigC(nOrb)
  double precision,intent(in)   :: Z(nOrb)
  double precision,intent(in)   :: eQP(nOrb)

  integer                       :: i,a
  double precision              :: eHOMO,eLUMO,Gap

! HOMO and LUMO

  eHOMO = maxval(eQP(1:nO))
  eLUMO = minval(eQP(nO+1:nOrb))
  Gap   = eLUMO - eHOMO

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)' 2SOSEX@RHF calculation '
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','Sig_QP (eV)','|','Z','|','e_QP (eV)','|'
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
  write(*,'(2X,A60,F15.6,A3)') '2SOSEX@RHF HOMO      energy = ',eHOMO*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') '2SOSEX@RHF LUMO      energy = ',eLUMO*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') '2SOSEX@RHF HOMO-LUMO gap    = ',Gap*HaToeV,' eV'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'phRPA@2SOSEX@RHF total       energy = ',ENuc + ERHF + EcRPA,' au'
  write(*,'(2X,A60,F15.6,A3)') 'phRPA@2SOSEX@RHF correlation energy = ',EcRPA,' au'
  write(*,'(2X,A60,F15.6,A3)') '   GM@2SOSEX@RHF total       energy = ',ENuc + ERHF + EcGM,' au'
  write(*,'(2X,A60,F15.6,A3)') '   GM@2SOSEX@RHF correlation energy = ',EcGM,' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

end subroutine 
