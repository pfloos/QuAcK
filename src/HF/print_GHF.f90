subroutine print_GHF(nBas,nBas2,nO,eHF,C,S,ENuc,ET,EV,EJ,EK,EGHF,dipole)


! Print one-electron energies and other stuff for GHF

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nBas2
  integer,intent(in)                 :: nO
  double precision,intent(in)        :: eHF(nBas2)

  double precision,intent(in)        :: C(nBas2,nBas2)
  double precision,intent(in)        :: S(nBas,nBas)
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ET
  double precision,intent(in)        :: EV
  double precision,intent(in)        :: EJ
  double precision,intent(in)        :: EK
  double precision,intent(in)        :: EGHF
  double precision,intent(in)        :: dipole(ncart)

! Local variables

  logical                            :: dump_orb = .false.

  integer                            :: i,j
  integer                            :: ixyz

  integer                            :: HOMO
  integer                            :: LUMO
  double precision                   :: Gap
  double precision                   :: Sx,Sy,Sz
  double precision                   :: SmSp,SpSm,Sz2,S2

  double precision,allocatable       :: Ca(:,:)
  double precision,allocatable       :: Cb(:,:)
  double precision,allocatable       :: Paa(:,:)
  double precision,allocatable       :: Pab(:,:)
  double precision,allocatable       :: Pba(:,:)
  double precision,allocatable       :: Pbb(:,:)

  double precision,external          :: trace_matrix

! HOMO and LUMO

  HOMO = nO
  LUMO = HOMO + 1
  Gap = eHF(LUMO)-eHF(HOMO)

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

  write(*,*)
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33)')           ' Summary               '
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33,1X,F16.10,A3)') ' One-electron energy = ',ET + EV,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Kinetic      energy = ',ET,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Potential    energy = ',EV,' au'
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33,1X,F16.10,A3)') ' Two-electron energy = ',EJ + EK,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Hartree      energy = ',EJ,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Exchange     energy = ',EK,' au'
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33,1X,F16.10,A3)') ' Electronic   energy = ',EGHF,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Nuclear   repulsion = ',ENuc,' au'
  write(*,'(A33,1X,F16.10,A3)') ' GHF          energy = ',EGHF + ENuc,' au'
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33,1X,F16.6,A3)')  ' GHF HOMO     energy = ',eHF(HOMO)*HaToeV,' eV'
  write(*,'(A33,1X,F16.6,A3)')  ' GHF LUMO     energy = ',eHF(LUMO)*HaToeV,' eV'
  write(*,'(A33,1X,F16.6,A3)')  ' GHF HOMO-LUMO gap   = ',Gap*HaToeV,' eV'
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33,1X,F16.6)')     ' <Sx>                = ',Sx
  write(*,'(A33,1X,F16.6)')     ' <Sy>                = ',Sy
  write(*,'(A33,1X,F16.6)')     ' <Sz>                = ',Sz
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33,1X,F16.6)')     ' <Sx^2+Sy^2>         = ',S2 - Sz2
  write(*,'(A33,1X,F16.6)')     ' <Sz^2>              = ',Sz2
  write(*,'(A33,1X,F16.6)')     ' <S^2>               = ',S2
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A36)')           ' Dipole moment (Debye)    '
  write(*,'(10X,4A10)')      'X','Y','Z','Tot.'
  write(*,'(10X,4F10.4)')    (dipole(ixyz)*auToD,ixyz=1,ncart),norm2(dipole)*auToD
  write(*,'(A50)')           '---------------------------------------'
  write(*,*)

! Print results

  if(dump_orb) then
    write(*,'(A50)') '---------------------------------------'
    write(*,'(A50)') ' GHF orbital coefficients '
    write(*,'(A50)') '---------------------------------------'
    call matout(nBas2,nBas2,C)
    write(*,*)
  end if
  write(*,'(A50)') '---------------------------------------'
  write(*,'(A50)') ' GHF orbital energies (au) '
  write(*,'(A50)') '---------------------------------------'
  call vecout(nBas2,eHF)
  write(*,*)

end subroutine 
