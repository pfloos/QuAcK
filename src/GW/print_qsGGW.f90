subroutine print_qsGGW(nBas,nBas2,nO,nSCF,Conv,thresh,eHF,eGW,c,S,SigC,Z,ENuc,ET,EV,EJ,Ex,EcGM,EcRPA,EqsGW,dipole)

! Print information for the generalized version of qsGW

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nBas2
  integer,intent(in)                 :: nO
  integer,intent(in)                 :: nSCF
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ET
  double precision,intent(in)        :: EV
  double precision,intent(in)        :: EJ
  double precision,intent(in)        :: Ex
  double precision,intent(in)        :: EcGM
  double precision,intent(in)        :: EcRPA
  double precision,intent(in)        :: Conv
  double precision,intent(in)        :: thresh
  double precision,intent(in)        :: eHF(nBas2)
  double precision,intent(in)        :: eGW(nBas2)
  double precision,intent(in)        :: c(nBas2,nBas2)
  double precision,intent(in)        :: S(nBas,nBas)
  double precision,intent(in)        :: SigC(nBas2,nBas2)
  double precision,intent(in)        :: Z(nBas2)
  double precision,intent(in)        :: EqsGW
  double precision,intent(in)        :: dipole(ncart)

! Local variables

  logical                            :: dump_orb = .false.

  integer                            :: i,j
  integer                            :: p,ixyz,HOMO,LUMO
  double precision                   :: Gap
  double precision,external          :: trace_matrix

  double precision                   :: Sx,Sy,Sz
  double precision                   :: SmSp,SpSm,Sz2,S2

  double precision,allocatable       :: Ca(:,:)
  double precision,allocatable       :: Cb(:,:)
  double precision,allocatable       :: Paa(:,:)
  double precision,allocatable       :: Pab(:,:)
  double precision,allocatable       :: Pba(:,:)
  double precision,allocatable       :: Pbb(:,:)

! Output variables

! HOMO and LUMO

  HOMO = nO
  LUMO = HOMO + 1
  Gap = eGW(LUMO)-eGW(HOMO)

! Density matrices

  allocate(Paa(nO,nO),Pab(nO,nO),Pba(nO,nO),Pbb(nO,nO))

  allocate(Ca(nBas,nO),Cb(nBas,nO))

  Ca(:,:) = C(     1:nBas ,1:nO)
  Cb(:,:) = C(nBas+1:nBas2,1:nO)

  Paa = matmul(transpose(Ca),matmul(S,Ca))
  Pab = matmul(transpose(Ca),matmul(S,Cb))
  Pba = matmul(transpose(Cb),matmul(S,Ca))
  Pbb = matmul(transpose(Cb),matmul(S,Cb))

! Compute components of S = (Sx,Sy,Sz)

  Sx = 0.5d0*(trace_matrix(nO,Pab) + trace_matrix(nO,Pba))
  Sy = 0.5d0*(trace_matrix(nO,Pab) - trace_matrix(nO,Pba))
  Sz = 0.5d0*(trace_matrix(nO,Paa) - trace_matrix(nO,Pbb))

! Compute <S^2> = <Sx^2> + <Sy^2> + <Sz^2>

  SpSm = 0d0
  do i=1,nO
    do j=1,nO
      SpSm = SpSm + Pab(i,i)*Pba(j,j) - Pab(i,j)*Pba(j,i)
    end do
  end do
  SpSm = trace_matrix(nO,Paa) + SpSm

  SmSp = 0d0
  do i=1,nO
    do j=1,nO
      SmSp = SmSp + Pba(i,i)*Pab(j,j) - Pba(i,j)*Pab(j,i)
    end do
  end do
  SmSp = trace_matrix(nO,Pbb) + SmSp

  Sz2 = 0d0
  do i=1,nO
    do j=1,nO
      Sz2 = Sz2 + (Paa(i,i) - Pbb(i,i))*(Paa(j,j) - Pbb(j,j)) - (Paa(i,j) - Pbb(i,j))**2
    end do
  end do
  Sz2 = 0.25d0*(dble(nO) + Sz2)

! Compute <S^2> from Sz^2, S^+S^- and S^-S^+

  S2 = Sz2 + 0.5d0*(SpSm + SmSp)

  call print_GHF_spin(nBas,nBas2,nO,C,S)

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  if(nSCF < 10) then
    write(*,'(1X,A20,I1,A1,I1,A16)')' Self-consistent qsG',nSCF,'W',nSCF,'@GHF calculation'
  elseif(nSCF < 100) then
    write(*,'(1X,A20,I2,A1,I2,A16)')' Self-consistent qsG',nSCF,'W',nSCF,'@GHF calculation'
  else
    write(*,'(1X,A20,I3,A1,I3,A16)')' Self-consistent qsG',nSCF,'W',nSCF,'@GHF calculation'
  end if
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','Sig_GW (eV)','|','Z','|','e_GW (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'

  do p=1,nBas
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',eHF(p)*HaToeV,'|',SigC(p,p)*HaToeV,'|',Z(p),'|',eGW(p)*HaToeV,'|'
  end do

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A10,I3)')   'Iteration ',nSCF
  write(*,'(2X,A14,F15.5)')'Convergence = ',Conv
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') 'qsGW@GHF HOMO     energy = ',eGW(HOMO)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'qsGW@GHF LUMO     energy = ',eGW(LUMO)*HaToeV,' eV'
  write(*,'(2X,A60,F15.6,A3)') 'qsGW@GHF HOMO-LUMO gap   = ',Gap*HaToeV,' eV'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') '      qsGW@GHF total       energy = ',ENuc + EqsGW,' au'
  write(*,'(2X,A60,F15.6,A3)') '      qsGW@GHF exchange    energy = ',Ex,' au'
  write(*,'(2X,A60,F15.6,A3)') '   GM@qsGW@GHF correlation energy = ',EcGM,' au'
  write(*,'(2X,A60,F15.6,A3)') 'phRPA@qsGW@GHF correlation energy = ',EcRPA,' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

! Dump results for final iteration

  if(Conv < thresh) then

    write(*,*)
    write(*,'(A50)')           '---------------------------------------'
    write(*,'(A33)')           ' Summary              '
    write(*,'(A50)')           '---------------------------------------'
    write(*,'(A33,1X,F16.10,A3)') ' One-electron energy = ',ET + EV,' au'
    write(*,'(A33,1X,F16.10,A3)') ' Kinetic      energy = ',ET,' au'
    write(*,'(A33,1X,F16.10,A3)') ' Potential    energy = ',EV,' au'
    write(*,'(A50)')           '---------------------------------------'
    write(*,'(A33,1X,F16.10,A3)') ' Two-electron energy = ',EJ + Ex,' au'
    write(*,'(A33,1X,F16.10,A3)') ' Hartree      energy = ',EJ,' au'
    write(*,'(A33,1X,F16.10,A3)') ' Exchange     energy = ',Ex,' au'
    write(*,'(A33,1X,F16.10,A3)') ' Correlation  energy = ',EcGM,' au'
    write(*,'(A50)')           '---------------------------------------'
    write(*,'(A33,1X,F16.10,A3)') ' Electronic   energy = ',EqsGW,' au'
    write(*,'(A33,1X,F16.10,A3)') ' Nuclear   repulsion = ',ENuc,' au'
    write(*,'(A33,1X,F16.10,A3)') ' qsGGW        energy = ',ENuc + EqsGW,' au'
    write(*,'(A50)')           '---------------------------------------'
    write(*,'(A36)')           ' Dipole moment (Debye)    '
    write(*,'(10X,4A10)')      'X','Y','Z','Tot.'
    write(*,'(10X,4F10.4)')    (dipole(ixyz)*auToD,ixyz=1,ncart),norm2(dipole)*auToD
    write(*,'(A50)')           '-----------------------------------------'
    write(*,*)
 
    if(dump_orb) then
      write(*,'(A50)') '---------------------------------------'
      write(*,'(A32)') ' Generalized qsGW orbital coefficients '
      write(*,'(A50)') '---------------------------------------'
      call matout(nBas,nBas,c)
      write(*,*)
    end if
    write(*,'(A50)') '---------------------------------------'
    write(*,'(A32)') ' Generalized qsGW orbital energies (au)'
    write(*,'(A50)') '---------------------------------------'
    call vecout(nBas,eGW)
    write(*,*)

  end if


end subroutine 
