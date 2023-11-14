subroutine print_GHF(nBas,nBas2,nO,e,C,P,ENuc,ET,EV,EJ,EK,EHF,dipole)

! Print one-electron energies and other stuff for GHF

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nBas2
  integer,intent(in)                 :: nO
  double precision,intent(in)        :: e(nBas2)
  double precision,intent(in)        :: C(nBas2,nBas2)
  double precision,intent(in)        :: P(nBas2,nBas2)
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ET
  double precision,intent(in)        :: EV
  double precision,intent(in)        :: EJ
  double precision,intent(in)        :: EK
  double precision,intent(in)        :: EHF
  double precision,intent(in)        :: dipole(ncart)

! Local variables

  integer                            :: ixyz
  integer                            :: mu,nu
  integer                            :: HOMO
  integer                            :: LUMO
  double precision                   :: Gap
  double precision                   :: Sx2,Sy2,Sz2,S2

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
  Gap = e(LUMO)-e(HOMO)

! Density matrices

  allocate(Paa(nBas,nBas),Pab(nBas,nBas),Pba(nBas,nBas),Pbb(nBas,nBas))

  Paa(:,:) = P(     1:nBas ,     1:nBas )
  Pab(:,:) = P(     1:nBas ,nBas+1:nBas2)
  Pba(:,:) = P(nBas+1:nBas2,     1:nBas )
  Pbb(:,:) = P(nBas+1:nBas2,nBas+1:nBas2)

! allocate(Ca(nBas,nBas2),Cb(nBas,nBas2))

! Ca(:,:) = C(     1:nBas ,1:nBas2)
! Cb(:,:) = C(nBas+1:nBas2,1:nBas2)

! Compute expectation values of S^2 (WRONG!)

! Sx2 = 0.25d0*trace_matrix(nBas,Paa+Pbb) + 0.25d0*trace_matrix(nBas,Pab+Pba)**2
! do mu=1,nBas
!   do nu=1,nBas
!       Sx2 = Sx2 - 0.5d0*(Paa(mu,nu)*Pbb(nu,mu) + Pab(mu,nu)*Pab(nu,mu))
!   end do
! end do

! Sy2 = 0.25d0*trace_matrix(nBas,Paa+Pbb) - 0.25d0*trace_matrix(nBas,Pab+Pba)**2
! do mu=1,nBas
!   do nu=1,nBas
!       Sy2 = Sy2 - 0.5d0*(Paa(mu,nu)*Pbb(nu,mu) - Pab(mu,nu)*Pab(nu,mu))
!   end do
! end do

! Sz2 = 0.25d0*trace_matrix(nBas,Paa+Pbb) + 0.25d0*trace_matrix(nBas,Pab-Pba)**2
! do mu=1,nBas
!   do nu=1,nBas
!       Sz2 = Sz2 - 0.25d0*(Paa(mu,nu)*Pbb(nu,mu) - Pab(mu,nu)*Pab(nu,mu))
!       Sz2 = Sz2 + 0.25d0*(Pab(mu,nu)*Pba(nu,mu) - Pba(mu,nu)*Pab(nu,mu))
!   end do
! end do
! 
! S2 = Sx2 + Sy2 + Sz2

! Checl collinearity and coplanarity 

  allocate(PP(nBas,nBas),Mx(nBas,nBas),My(nBas,nBas),Mz(nBas,nBas))

  PP(:,:) = 0.5d0*(Paa(:,:) + Pbb(:,:))
  Mx(:,:) = 0.5d0*(Pba(:,:) + Pab(:,:))
  My(:,:) = 0.5d0*(Pba(:,:) - Pab(:,:))
  Mz(:,:) = 0.5d0*(Paa(:,:) - Pbb(:,:))

  T(1,1) = trace_matrix(nBas,matmul(Mx,transpose(Mx)))
  T(1,2) = trace_matrix(nBas,matmul(Mx,transpose(My)))
  T(1,3) = trace_matrix(nBas,matmul(Mx,transpose(Mz)))
  T(2,1) = trace_matrix(nBas,matmul(My,transpose(Mx)))
  T(2,2) = trace_matrix(nBas,matmul(My,transpose(My)))
  T(2,3) = trace_matrix(nBas,matmul(My,transpose(Mz)))
  T(3,1) = trace_matrix(nBas,matmul(Mz,transpose(Mx)))
  T(3,2) = trace_matrix(nBas,matmul(Mz,transpose(My)))
  T(3,3) = trace_matrix(nBas,matmul(Mz,transpose(Mz)))

  print*,'Value of Tr(P - P^2)'
  lambda = trace_matrix(nBas,PP - matmul(PP,transpose(PP)))
  print*,lambda

  print*,'Eigenvalues of T'
  vec(:,:) = T(:,:)
  call diagonalize_matrix(3,vec,val)
  print*,val

  T(1,1) = - T(1,1) + lambda
  T(2,2) = - T(2,2) + lambda 
  T(3,3) = - T(3,3) + lambda

  print*,'Eigenvalues of A'
  vec(:,:) = T(:,:)
  call diagonalize_matrix(3,vec,val)
  print*,val

  deallocate(PP,Mx,My,Mz)

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
  write(*,'(A33,1X,F16.10,A3)') ' Electronic   energy = ',EHF,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Nuclear   repulsion = ',ENuc,' au'
  write(*,'(A33,1X,F16.10,A3)') ' GHF          energy = ',EHF + ENuc,' au'
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33,1X,F16.6,A3)')  ' GHF HOMO     energy = ',e(HOMO)*HaToeV,' eV'
  write(*,'(A33,1X,F16.6,A3)')  ' GHF LUMO     energy = ',e(LUMO)*HaToeV,' eV'
  write(*,'(A33,1X,F16.6,A3)')  ' GHF HOMO-LUMO gap   = ',Gap*HaToeV,' eV'
  write(*,'(A50)')           '---------------------------------------'
! write(*,'(A32,1X,F16.6)')     ' <S**2>             :',S2
! write(*,'(A50)')           '---------------------------------------'
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
    call matout(nBas2,nBas2,c)
    write(*,*)
  end if
  write(*,'(A50)') '---------------------------------------'
  write(*,'(A50)') ' GHF orbital energies (au) '
  write(*,'(A50)') '---------------------------------------'
  call matout(nBas2,1,e)
  write(*,*)

end subroutine 
