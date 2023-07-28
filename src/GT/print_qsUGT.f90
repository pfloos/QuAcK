subroutine print_qsUGT(nBas,nO,nSCF,Conv,thresh,eHF,eGT,c,SigT,Z,ENuc,ET,EV,EJ,Ex,EcGM,EcRPA,EqsGT,dipole)

! Print one-electron energies and other stuff for UG0T0

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nO(nspin)
  integer,intent(in)                 :: nSCF
  double precision,intent(in)        :: Conv
  double precision,intent(in)        :: thresh
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ET(nspin)
  double precision,intent(in)        :: EV(nspin)
  double precision,intent(in)        :: EJ(nsp)
  double precision,intent(in)        :: Ex(nspin)
  double precision,intent(in)        :: EcGM(nspin)
  double precision,intent(in)        :: EcRPA(nspin)
  double precision,intent(in)        :: EqsGT
  double precision,intent(in)        :: eHF(nBas,nspin)
  double precision,intent(in)        :: SigT(nBas,nBas,nspin)
  double precision,intent(in)        :: Z(nBas,nspin)
  double precision,intent(in)        :: eGT(nBas,nspin)
  double precision,intent(in)        :: c(nBas,nBas,nspin)
  double precision,intent(in)        :: dipole(ncart)

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
    write(*,'(1X,A21,I1,A1,I1,A12)')'  Self-consistent qsG',nSCF,'T',nSCF,' calculation'
  else
    write(*,'(1X,A21,I2,A1,I2,A12)')'  Self-consistent qsG',nSCF,'T',nSCF,' calculation'
  endif
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','Sigma_T (eV)','|','Z','|','e_QP (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'

  do p=1,nBas
    write(*,'(A1,I3,A1,2F15.6,A1,2F15.6,A1,2F15.6,A1,2F15.6,A1)') &
 '|',p,'|',eHF(p,1)*HaToeV,eHF(p,2)*HaToeV,'|',SigT(p,p,1)*HaToeV,SigT(p,p,2)*HaToeV,'|', &
 Z(p,1),Z(p,2),'|',eGT(p,1)*HaToeV,eGT(p,2)*HaToeV,'|'
  enddo

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A10,I3)')   'Iteration ',nSCF
  write(*,'(2X,A14,F15.5)')'Convergence = ',Conv
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F15.6,A3)') 'qsUGT HOMO      energy (eV)            =',maxval(HOMO(:))*HaToeV,' eV'
  write(*,'(2X,A50,F15.6,A3)') 'qsUGT LUMO      energy (eV)            =',minval(LUMO(:))*HaToeV,' eV'
  write(*,'(2X,A50,F15.6,A3)') 'qsUGT HOMO-LUMO gap    (eV)            =',(minval(LUMO(:))-maxval(HOMO(:)))*HaToeV,' eV'
  write(*,*)'-------------------------------------------------------------------------------'
write(*,'(2X,A50,F15.6,A3)') '      qsGT total       energy:',ENuc + EqsGT,' au'
write(*,'(2X,A50,F15.6,A3)') '      qsGT exchange    energy:',sum(Ex(:)),' au'
write(*,'(2X,A50,F15.6,A3)') '   GM@qsGT correlation energy:',sum(EcGM(:)),' au'
write(*,'(2X,A50,F15.6,A3)') 'ppRPA@qsGT correlation energy:',sum(EcRPA(:)),' au'
write(*,*)'-------------------------------------------'
write(*,*) 

end subroutine 
