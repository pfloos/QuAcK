subroutine print_qsGW(nBas,nO,nSCF,Conv,thresh,eHF,eGW,c,ENuc,P,T,V,Hc,J,K,F,SigmaC,Z,EcRPA,EcGM)


! Print one-electron energies and other stuff for qsGW

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                 :: nBas,nO,nSCF
  double precision,intent(in)        :: ENuc,EcRPA,EcGM,Conv,thresh
  double precision,intent(in)        :: eHF(nBas),eGW(nBas),c(nBas),P(nBas,nBas) 
  double precision,intent(in)        :: T(nBas,nBas),V(nBas,nBas),Hc(nBas,nBas)
  double precision,intent(in)        :: J(nBas,nBas),K(nBas,nBas),F(nBas,nBas)
  double precision,intent(in)        :: Z(nBas),SigmaC(nBas,nBas)

! Local variables

  integer                            :: x,HOMO,LUMO
  double precision                   :: Gap,ET,EV,EJ,Ex,Ec,EqsGW
  double precision,external          :: trace_matrix


! HOMO and LUMO

  HOMO = nO
  LUMO = HOMO + 1
  Gap = eGW(LUMO)-eGW(HOMO)

  ET = trace_matrix(nBas,matmul(P,T))
  EV = trace_matrix(nBas,matmul(P,V))
  EJ = 0.5d0*trace_matrix(nBas,matmul(P,J))
  Ex = 0.5d0*trace_matrix(nBas,matmul(P,K))
  EqsGW = ET + EV + EJ + Ex
  Ec = 0d0

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  if(nSCF < 10) then
    write(*,'(1X,A21,I1,A1,I1,A12)')'  Self-consistent qsG',nSCF,'W',nSCF,' calculation'
  else
    write(*,'(1X,A21,I2,A1,I2,A12)')'  Self-consistent qsG',nSCF,'W',nSCF,' calculation'
  endif
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','e_QP-e_HF (eV)','|','Z','|','e_QP (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'

  do x=1,nBas
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',x,'|',eHF(x)*HaToeV,'|',(eGW(x)-eHF(x))*HaToeV,'|',Z(x),'|',eGW(x)*HaToeV,'|'
  enddo


  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A10,I3)')   'Iteration ',nSCF
  write(*,'(2X,A19,F15.5)')'max(|FPS - SPF|) = ',Conv
  write(*,*)'-------------------------------------------'
  write(*,'(2X,A27,F15.6)') 'qsGW HOMO      energy (eV):',eGW(HOMO)*HaToeV
  write(*,'(2X,A27,F15.6)') 'qsGW LUMO      energy (eV):',eGW(LUMO)*HaToeV
  write(*,'(2X,A27,F15.6)') 'qsGW HOMO-LUMO gap    (eV):',Gap*HaToeV
  write(*,*)'-------------------------------------------'
  write(*,'(2X,A27,F15.6)') 'qsGW total       energy   =',EqsGW + ENuc
  write(*,'(2X,A27,F15.6)') 'qsGW GM total    energy   =',EqsGW + ENuc + EcGM
  write(*,'(2X,A27,F15.6)') 'qsGW exchange    energy   =',Ex
  write(*,'(2X,A27,F15.6)') 'qsGW correlation energy   =',Ec
  write(*,'(2X,A27,F15.6)') 'RPA  correlation energy   =',EcRPA
  write(*,'(2X,A27,F15.6)') 'GM  correlation energy    =',EcGM
  write(*,*)'-------------------------------------------'
  write(*,*)

! Dump results for final iteration

  if(Conv < thresh) then

    write(*,*)
    write(*,'(A50)')           '---------------------------------------'
    write(*,'(A32)')           ' Summary              '
    write(*,'(A50)')           '---------------------------------------'
    write(*,'(A32,1X,F16.10)') ' One-electron energy  ',ET + EV
    write(*,'(A32,1X,F16.10)') ' Kinetic      energy  ',ET
    write(*,'(A32,1X,F16.10)') ' Potential    energy  ',EV
    write(*,'(A50)')           '---------------------------------------'
    write(*,'(A32,1X,F16.10)') ' Two-electron energy  ',EJ + Ex
    write(*,'(A32,1X,F16.10)') ' Coulomb      energy  ',EJ
    write(*,'(A32,1X,F16.10)') ' Exchange     energy  ',Ex
    write(*,'(A32,1X,F16.10)') ' Correlation  energy  ',Ec
    write(*,'(A50)')           '---------------------------------------'
    write(*,'(A32,1X,F16.10)') ' Electronic   energy  ',EqsGW
    write(*,'(A32,1X,F16.10)') ' Nuclear   repulsion  ',ENuc
    write(*,'(A32,1X,F16.10)') ' qsGW         energy  ',ENuc + EqsGW
    write(*,'(A32,1X,F16.10)') ' RPA corr.    energy  ',EcRPA
    write(*,'(A50)')           '---------------------------------------'
    write(*,*)
 
    write(*,'(A50)') '---------------------------------------'
    write(*,'(A32)') ' qsGW MO coefficients'
    write(*,'(A50)') '---------------------------------------'
    call matout(nBas,nBas,c)
    write(*,*)
    write(*,'(A50)') '---------------------------------------'
    write(*,'(A32)') ' qsGW MO energies'
    write(*,'(A50)') '---------------------------------------'
    call matout(nBas,1,eGW)
    write(*,*)

  endif


end subroutine print_qsGW
