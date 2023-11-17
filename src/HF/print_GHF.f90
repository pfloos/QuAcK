subroutine print_GHF(nBas,nBas2,nO,eHF,C,P,S,ENuc,ET,EV,EJ,EK,EGHF,dipole)


! Print one-electron energies and other stuff for GHF

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nBas2
  integer,intent(in)                 :: nO
  double precision,intent(in)        :: eHF(nBas2)

  double precision,intent(in)        :: C(nBas2,nBas2)
  double precision,intent(in)        :: P(nBas2,nBas2)
  double precision,intent(in)        :: S(nBas,nBas)
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ET
  double precision,intent(in)        :: EV
  double precision,intent(in)        :: EJ
  double precision,intent(in)        :: EK
  double precision,intent(in)        :: EGHF
  double precision,intent(in)        :: dipole(ncart)

! Local variables

  integer                            :: i,j
  integer                            :: ixyz

  integer                            :: mu,nu
  integer                            :: HOMO
  integer                            :: LUMO
  double precision                   :: Gap
  double precision                   :: Sx,Sy,Sz
  double precision                   :: SmSp,SpSm,Sz2,S2
! double precision                   :: na, nb
! double precision                   :: nonco_z, contam_uhf, xy_perp, contam_ghf

  double precision,allocatable       :: Ca(:,:)
  double precision,allocatable       :: Cb(:,:)
  double precision,allocatable       :: Paa(:,:)
  double precision,allocatable       :: Pab(:,:)
  double precision,allocatable       :: Pba(:,:)
  double precision,allocatable       :: Pbb(:,:)
  double precision,allocatable       :: tmp(:,:)

  double precision,allocatable       :: Mx(:,:)
  double precision,allocatable       :: My(:,:)
  double precision,allocatable       :: Mz(:,:)
  double precision,allocatable       :: PP(:,:)
  double precision                   :: T(3,3)
  double precision                   :: vec(3,3)
  double precision                   :: val(3)
  double precision                   :: lambda

  double precision,external          :: trace_matrix

  logical                            :: dump_orb = .false.

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




! deallocate(Paa,Pab,Pba,Pbb)

! Check collinearity and coplanarity 

! allocate(PP(nO,nO),Mx(nO,nO),My(nO,nO),Mz(nO,nO))

! PP(:,:) = 0.5d0*(Paa(:,:) + Pbb(:,:))
! Mx(:,:) = 0.5d0*(Pba(:,:) + Pab(:,:))
! My(:,:) = 0.5d0*(Pba(:,:) - Pab(:,:))
! Mz(:,:) = 0.5d0*(Paa(:,:) - Pbb(:,:))

! T(1,1) = trace_matrix(nO,matmul(Mx,Mx))
! T(1,2) = trace_matrix(nO,matmul(Mx,My))
! T(1,3) = trace_matrix(nO,matmul(Mx,Mz))
! T(2,1) = trace_matrix(nO,matmul(My,Mx))
! T(2,2) = trace_matrix(nO,matmul(My,My))
! T(2,3) = trace_matrix(nO,matmul(My,Mz))
! T(3,1) = trace_matrix(nO,matmul(Mz,Mx))
! T(3,2) = trace_matrix(nO,matmul(Mz,My))
! T(3,3) = trace_matrix(nO,matmul(Mz,Mz))

! lambda = trace_matrix(nO,PP - matmul(PP,PP))
! write(*,'(A,F10.6)') 'Tr(P - P^2) = ',lambda

! vec(:,:) = T(:,:)
! call diagonalize_matrix(3,vec,val)
! write(*,'(A,3F10.6)') 'Eigenvalues of T = ',val

! T(1,1) = - T(1,1) + lambda
! T(2,2) = - T(2,2) + lambda 
! T(3,3) = - T(3,3) + lambda

! vec(:,:) = T(:,:)
! call diagonalize_matrix(3,vec,val)
! write(*,'(A,3F10.6)') 'Eigenvalues of A = ',val

! deallocate(PP,Mx,My,Mz)


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
  write(*,'(A33,1X,F16.6)')     ' <Sz**2>             = ',Sz2
  write(*,'(A33,1X,F16.6)')     ' <S**2>              = ',S2
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

  call print_GHFspin(nBas, nBas2, nO, C, S)

end subroutine 

! ---

subroutine print_GHFspin(nBas, nBas2, nO, C, S)

  implicit none

  integer,          intent(in)  :: nBas, nBas2, nO
  double precision, intent(in)  :: C(nBas2,nBas2), S(nBas,nBas)

  integer                       :: i, j
  double precision              :: Na, Nb
  double precision              :: nonco_z, contam_ghf
  double precision              :: S2, Sz, Sz2
  double precision, allocatable :: Ca(:,:), Cb(:,:)
  double precision, allocatable :: Paa(:,:), Pab(:,:), Pba(:,:), Pbb(:,:)
  double precision, allocatable :: Mc(:,:), Eigc(:)

  print *, ' Spin properties for GHF WF:'

  allocate(Ca(nBas,nO), Cb(nBas,nO))
  do i = 1, nO
    do j = 1, nBas
      Ca(j,i) = C(     j,i)
      Cb(j,i) = C(nBas+j,i)
    enddo
  enddo

  ! TODO DGEMM
  allocate(Paa(nO,nO), Pab(nO,nO), Pba(nO,nO), Pbb(nO,nO))
  Paa = matmul(transpose(Ca), matmul(S, Ca))
  Pab = matmul(transpose(Ca), matmul(S, Cb))
  Pba = matmul(transpose(Cb), matmul(S, Ca))
  Pbb = matmul(transpose(Cb), matmul(S, Cb))

  deallocate(Ca, Cb)

  Na = 0.d0
  Nb = 0.d0
  do i = 1, nO
    Na = Na + Paa(i,i)
    Nb = Nb + Pbb(i,i)
  enddo

  nonco_z = dble(nO)
  do j = 1, nO
    do i = 1, nO
      nonco_z = nonco_z - (Paa(i,j) - Pbb(i,j))**2
    enddo
  enddo
  nonco_z = 0.25d0 * nonco_z

  contam_ghf = 0.d0
  do j = 1, nO
    do i = 1, nO
      contam_ghf = contam_ghf - (Pab(i,i)*Pba(j,j) - Pab(i,j)*Pba(j,i))
    enddo
  enddo

  Sz  = 0.5d0 * (Na - Nb)
  Sz2 = Sz*Sz + nonco_z
  S2  = Sz * (Sz + 1.d0) + nonco_z + contam_ghf

  print *, 'Sz, Sz^2 = ', Sz, Sz2
  print *, 'S^2      = ', S2
  

  ! --- --- --- --- --- --- --- --- ---
  ! calculate the axis of Collinearity
  ! --- --- --- --- --- --- --- --- ---

  allocate(Mc(3,3), Eigc(3))

  Mc(:,:) = 0.d0
  Mc(1,1) = 0.25d0 * dble(nO)
  Mc(2,2) = 0.25d0 * dble(nO)
  Mc(3,3) = 0.25d0 * dble(nO)
  do j = 1, nO
    do i = 1, nO
      Mc(1,1) = Mc(1,1) - 0.25d0 * (Pba(i,j) + Pab(i,j))**2
      Mc(2,2) = Mc(2,2) - 0.25d0 * (Pba(i,j) - Pab(i,j))**2
      Mc(3,3) = Mc(3,3) - 0.25d0 * (Paa(i,j) - Pbb(i,j))**2
      Mc(1,3) = Mc(1,3) - 0.25d0 * (Pab(i,j) + Pba(i,j))*(Paa(j,i) - Pbb(j,j))
    enddo
  enddo
  Mc(3,1) = Mc(1,3)

  call diagonalize_matrix(3, Mc, Eigc)
  print *, ' eigenvalues of Collinearity matrix:', Eigc
  deallocate(Mc, Eigc)

  ! --- --- --- --- --- --- --- --- ---
  ! --- --- --- --- --- --- --- --- ---

  deallocate(Paa, Pab, Pba, Pbb)

end



