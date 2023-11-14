subroutine print_UG0F2(nBas,nO,eHF,ENuc,EUHF,SigC,Z,eGF2,Ec)

! Print one-electron energies and other stuff for G0W0

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nO(nspin)
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: EUHF
  double precision,intent(in)        :: Ec(nsp)
  double precision,intent(in)        :: eHF(nBas,nspin)
  double precision,intent(in)        :: SigC(nBas,nspin)
  double precision,intent(in)        :: Z(nBas,nspin)
  double precision,intent(in)        :: eGF2(nBas,nspin)

  integer                            :: p
  integer                            :: ispin
  double precision                   :: HOMO(nspin)
  double precision                   :: LUMO(nspin)
  double precision                   :: Gap(nspin)

! HOMO and LUMO

  do ispin=1,nspin
    if(nO(ispin) > 0) then
      HOMO(ispin) = eGF2(nO(ispin),ispin)
      LUMO(ispin) = eGF2(nO(ispin)+1,ispin)
      Gap(ispin)  = LUMO(ispin) - HOMO(ispin)
    else
      HOMO(ispin) = 0d0
      LUMO(ispin) = eGF2(1,ispin)
      Gap(ispin)  = 0d0
    end if
  end do

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'&
              '------------------------------------------------'
  write(*,*)'  Unrestricted one-shot G0F2 calculation (eV)'
  write(*,*)'-------------------------------------------------------------------------------'&
              '------------------------------------------------'
  write(*,'(A1,A3,A1,A30,A1,A30,A1,A30,A1,A30,A1)') &
            '|',' ','|','e_HF            ','|','Sig_c            ','|','Z            ','|','e_QP            ','|'
  write(*,'(A1,A3,A1,2A15,A1,2A15,A1,2A15,A1,2A15,A1)') &
            '|','#','|','up     ','dw     ','|','up     ','dw     ','|','up     ','dw     ','|','up     ','dw     ','|'
  write(*,*)'-------------------------------------------------------------------------------'&
              '------------------------------------------------'

  do p=1,nBas
    write(*,'(A1,I3,A1,2F15.6,A1,2F15.6,A1,2F15.6,A1,2F15.6,A1)') &
    '|',p,'|',eHF(p,1)*HaToeV,eHF(p,2)*HaToeV,'|',SigC(p,1)*HaToeV,SigC(p,2)*HaToeV,'|', & 
              Z(p,1),Z(p,2),'|',eGF2(p,1)*HaToeV,eGF2(p,2)*HaToeV,'|'
  enddo

  write(*,*)'-------------------------------------------------------------------------------'&
              '------------------------------------------------'
  write(*,'(2X,A30,F15.6,A3)') 'UG0F2 HOMO      energy:',maxval(HOMO(:))*HaToeV,' eV'
  write(*,'(2X,A30,F15.6,A3)') 'UG0F2 LUMO      energy:',minval(LUMO(:))*HaToeV,' eV'
  write(*,'(2X,A30,F15.6,A3)') 'UG0F2 HOMO-LUMO gap   :',(minval(LUMO(:))-maxval(HOMO(:)))*HaToeV,' eV'
  write(*,*)'-------------------------------------------------------------------------------'&
              '------------------------------------------------'
  write(*,'(2X,A30,F15.6,A3)') '    UG0F2 total energy      :',ENuc + EUHF + sum(Ec(:)),' au'
  write(*,'(2X,A30,F15.6,A3)') '    UG0F2 correlation energy:',sum(Ec(:)),' au'
  write(*,*)'-------------------------------------------------------------------------------'&
              '------------------------------------------------'
  write(*,*)

end subroutine 
