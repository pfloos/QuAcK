subroutine print_UG0W0(nBas,nO,e,ENuc,EHF,SigC,Z,eGW,EcRPA)

! Print one-electron energies and other stuff for G0W0

  implicit none
  include 'parameters.h'

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nO(nspin)
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: EHF
  double precision,intent(in)        :: EcRPA
  double precision,intent(in)        :: e(nBas,nspin)
  double precision,intent(in)        :: SigC(nBas,nspin)
  double precision,intent(in)        :: Z(nBas,nspin)
  double precision,intent(in)        :: eGW(nBas,nspin)

  integer                            :: p
  integer                            :: ispin
  double precision                   :: HOMO(nspin)
  double precision                   :: LUMO(nspin)
  double precision                   :: Gap(nspin)

! HOMO and LUMO

  do ispin=1,nspin
    if(nO(ispin) > 0) then
      HOMO(ispin) = eGW(nO(ispin),ispin)
      LUMO(ispin) = eGW(nO(ispin)+1,ispin)
      Gap(ispin)  = LUMO(ispin) - HOMO(ispin)
    else
      HOMO(ispin) = 0d0
      LUMO(ispin) = e(1,ispin)
      Gap(ispin)  = 0d0
    end if
  end do

! Dump results

  write(*,*)'-------------------------------------------------------------------------------&
              -------------------------------------------------'
  write(*,*)'  Unrestricted one-shot G0W0 calculation (eV)'
  write(*,*)'-------------------------------------------------------------------------------& 
              -------------------------------------------------'
  write(*,'(A1,A3,A1,A30,A1,A30,A1,A30,A1,A30,A1)') &
            '|',' ','|','e_HF            ','|','Sig_c            ','|','Z            ','|','e_QP            ','|'
  write(*,'(A1,A3,A1,2A15,A1,2A15,A1,2A15,A1,2A15,A1)') &
            '|','#','|','up     ','dw     ','|','up     ','dw     ','|','up     ','dw     ','|','up     ','dw     ','|'
  write(*,*)'-------------------------------------------------------------------------------& 
              -------------------------------------------------'

  do p=1,nBas
    write(*,'(A1,I3,A1,2F15.6,A1,2F15.6,A1,2F15.6,A1,2F15.6,A1)') &
    '|',p,'|',e(p,1)*HaToeV,e(p,2)*HaToeV,'|',SigC(p,1)*HaToeV,SigC(p,2)*HaToeV,'|', & 
              Z(p,1),Z(p,2),'|',eGW(p,1)*HaToeV,eGW(p,2)*HaToeV,'|'
  enddo

  write(*,*)'-------------------------------------------------------------------------------& 
              -------------------------------------------------'
  write(*,'(2X,A30,F15.6)') 'G0W0 HOMO      energy (eV):',maxval(HOMO(:))*HaToeV
  write(*,'(2X,A30,F15.6)') 'G0W0 LUMO      energy (eV):',minval(LUMO(:))*HaToeV
  write(*,'(2X,A30,F15.6)') 'G0W0 HOMO-LUMO gap    (eV):',(minval(LUMO(:))-maxval(HOMO(:)))*HaToeV
  write(*,*)'-------------------------------------------------------------------------------& 
              -------------------------------------------------'
  write(*,'(2X,A30,F15.6)') 'RPA@G0W0 total energy       =',ENuc + EHF + EcRPA
  write(*,'(2X,A30,F15.6)') 'RPA@G0W0 correlation energy =',EcRPA
  write(*,*)'-------------------------------------------------------------------------------& 
              -------------------------------------------------'
  write(*,*)

end subroutine print_UG0W0

