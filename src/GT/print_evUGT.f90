subroutine print_evUGT(nBas,nO,nSCF,Conv,eHF,ENuc,EUHF,SigT,Z,eGT,EcGM,EcRPA)

! Print one-electron energies and other stuff for UG0T0

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nO(nspin)
  integer,intent(in)                 :: nSCF
  double precision,intent(in)        :: Conv
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: EUHF
  double precision,intent(in)        :: EcGM
  double precision,intent(in)        :: EcRPA(nspin)
  double precision,intent(in)        :: eHF(nBas,nspin)
  double precision,intent(in)        :: SigT(nBas,nspin)
  double precision,intent(in)        :: Z(nBas,nspin)
  double precision,intent(in)        :: eGT(nBas,nspin)

  integer                            :: p
  integer                            :: ispin 
  double precision                   :: HOMO(nspin)
  double precision                   :: LUMO(nspin)
  double precision                   :: Gap(nspin)

! HOMO and LUMO
  do ispin=1,nspin
    if(nO(ispin) > 0) then
      HOMO(ispin) = eGT(nO(ispin),ispin)
      LUMO(ispin) = eGT(nO(ispin)+1,ispin)
      Gap(ispin)  = LUMO(ispin) - HOMO(ispin)
    else
      HOMO(ispin) = 0d0
      LUMO(ispin) = eGT(1,ispin)
      Gap(ispin)  = 0d0
    end if
  end do

! Dump results
 
  write(*,*)'-------------------------------------------------------------------------------'
  if(nSCF < 10) then
    write(*,'(1X,A21,I1,A1,I1,A12)')'  Self-consistent evG',nSCF,'T',nSCF,' calculation'
  else
    write(*,'(1X,A21,I2,A1,I2,A12)')'  Self-consistent evG',nSCF,'T',nSCF,' calculation'
  end if
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','Sigma_T (eV)','|','Z','|','e_QP (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'

  do p=1,nBas
    write(*,'(A1,I3,A1,2F15.6,A1,2F15.6,A1,2F15.6,A1,2F15.6,A1)') &
 '|',p,'|',eHF(p,1)*HaToeV,eHF(p,2)*HaToeV,'|',SigT(p,1)*HaToeV,SigT(p,2)*HaToeV,'|', &
 Z(p,1),Z(p,2),'|',eGT(p,1)*HaToeV,eGT(p,2)*HaToeV,'|'
  end do

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A10,I3)')   'Iteration ',nSCF
  write(*,'(2X,A14,F15.5)')'Convergence = ',Conv
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F15.6,A3)') 'evUGT HOMO      energy (eV)            =',maxval(HOMO(:))*HaToeV,' eV'
  write(*,'(2X,A50,F15.6,A3)') 'evUGT LUMO      energy (eV)            =',minval(LUMO(:))*HaToeV,' eV'
  write(*,'(2X,A50,F15.6,A3)') 'evUGT HOMO-LUMO gap    (eV)            =',(minval(LUMO(:))-maxval(HOMO(:)))*HaToeV,' eV'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F20.10,A3)') ' Tr@ppRPA@evUGT correlation energy (singlet) =',EcRPA(1),' au'
  write(*,'(2X,A50,F20.10,A3)') ' Tr@ppRPA@evUGT correlation energy (triplet) =',EcRPA(2),' au'
  write(*,'(2X,A50,F20.10,A3)') ' Tr@ppRPA@evUGT correlation energy           =',EcRPA(1) + EcRPA(2),' au'
  write(*,'(2X,A50,F20.10,A3)') ' Tr@ppRPA@evUGT total energy                 =',ENuc + EUHF + EcRPA(1) + EcRPA(2),' au'
  write(*,'(2X,A50,F20.10,A3)') '       GM@evUGT correlation energy           =',EcGM,' au'
  write(*,'(2X,A50,F20.10,A3)') '       GM@evUGT total energy                 =',ENuc + EUHF + EcGM,' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

end subroutine 
