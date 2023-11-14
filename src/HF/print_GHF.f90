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

  integer                            :: ixyz
  integer                            :: mu,nu
  integer                            :: HOMO
  integer                            :: LUMO
  double precision                   :: Gap
  double precision                   :: Sx ,Sy ,Sz
  double precision                   :: Sx2,Sy2,Sz2

  double precision,allocatable       :: Ca(:,:)
  double precision,allocatable       :: Cb(:,:)
  double precision,allocatable       :: Paa(:,:)
  double precision,allocatable       :: Pab(:,:)
  double precision,allocatable       :: Pba(:,:)
  double precision,allocatable       :: Pbb(:,:)

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

! Paa(:,:) = P(     1:nBas ,     1:nBas )
! Pab(:,:) = P(     1:nBas ,nBas+1:nBas2)
! Pba(:,:) = P(nBas+1:nBas2,     1:nBas )
! Pbb(:,:) = P(nBas+1:nBas2,nBas+1:nBas2)

  allocate(Ca(nBas,nBas2),Cb(nBas,nBas2))

  Ca(:,:) = C(     1:nBas ,1:nO)
  Cb(:,:) = C(nBas+1:nBas2,1:nO)

  Paa = matmul(transpose(Ca),matmul(S,Ca))
  Pab = matmul(transpose(Ca),matmul(S,Cb))
  Pba = matmul(transpose(Cb),matmul(S,Ca))
  Pbb = matmul(transpose(Cb),matmul(S,Cb))

! Compute <S>

  Sx = 0.5d0*(trace_matrix(nO,Pab) + trace_matrix(nO,Pba))
  Sy = 0.5d0*(trace_matrix(nO,Pab) - trace_matrix(nO,Pba))
  Sz = 0.5d0*(trace_matrix(nO,Paa) - trace_matrix(nO,Pbb))

! Compute <S^2>

  Sx2 = 0.25d0*trace_matrix(nO,Paa+Pbb) + 0.25d0*trace_matrix(nO,Pab+Pba)**2 &
      - 0.5d0*trace_matrix(nO,matmul(Paa,Pbb) + matmul(Pab,Pab))

  Sy2 = 0.25d0*trace_matrix(nO,Paa+Pbb) - 0.25d0*trace_matrix(nO,Pab+Pba)**2 &
      - 0.5d0*trace_matrix(nO,matmul(Paa,Pbb) - matmul(Pab,Pab))

  Sz2 = 0.25d0*trace_matrix(nO,Paa+Pbb) + 0.25d0*trace_matrix(nO,Pab-Pba)**2 &
      - 0.25d0*trace_matrix(nO,matmul(Paa,Paa) - matmul(Pbb,Pbb)) &
      + 0.25d0*trace_matrix(nO,matmul(Pab,Pba) - matmul(Pba,Pab))

! deallocate(Paa,Pab,Pba,Pbb)

! Check collinearity and coplanarity 

  allocate(PP(nO,nO),Mx(nO,nO),My(nO,nO),Mz(nO,nO))

  PP(:,:) = 0.5d0*(Paa(:,:) + Pbb(:,:))
  Mx(:,:) = 0.5d0*(Pba(:,:) + Pab(:,:))
  My(:,:) = 0.5d0*(Pba(:,:) - Pab(:,:))
  Mz(:,:) = 0.5d0*(Paa(:,:) - Pbb(:,:))

! T(1,1) = trace_matrix(nBas,matmul(Mx,transpose(Mx)))
! T(1,2) = - trace_matrix(nBas,matmul(Mx,transpose(My)))
! T(1,3) = trace_matrix(nBas,matmul(Mx,transpose(Mz)))
! T(2,1) = - trace_matrix(nBas,matmul(My,transpose(Mx)))
! T(2,2) = + trace_matrix(nBas,matmul(My,transpose(My)))
! T(2,3) = - trace_matrix(nBas,matmul(My,transpose(Mz)))
! T(3,1) = trace_matrix(nBas,matmul(Mz,transpose(Mx)))
! T(3,2) = - trace_matrix(nBas,matmul(Mz,transpose(My)))
! T(3,3) = trace_matrix(nBas,matmul(Mz,transpose(Mz)))

  print*,2d0*trace_matrix(nO,PP)
! print*,'Value of Tr(P - P^2)'
! lambda = trace_matrix(nBas,PP - matmul(PP,transpose(PP)))
! print*,lambda

! print*,'Eigenvalues of T'
! vec(:,:) = T(:,:)
! call diagonalize_matrix(3,vec,val)
! print*,val

! T(1,1) = - T(1,1) + lambda
! T(2,2) = - T(2,2) + lambda 
! T(3,3) = - T(3,3) + lambda

! print*,'Eigenvalues of A'
! vec(:,:) = T(:,:)
! call diagonalize_matrix(3,vec,val)
! print*,val

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
  write(*,'(A33,1X,F16.6)')     ' <S>                 = ',Sx+Sy+Sz
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33,1X,F16.6)')     ' <Sx**2>             = ',Sx2
  write(*,'(A33,1X,F16.6)')     ' <Sy**2>             = ',Sy2
  write(*,'(A33,1X,F16.6)')     ' <Sz**2>             = ',Sz2
  write(*,'(A33,1X,F16.6)')     ' <S**2>              = ',Sx2+Sy2+Sz2
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
  call matout(nBas2,1,eHF)
  write(*,*)

end subroutine 
