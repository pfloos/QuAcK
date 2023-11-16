subroutine print_GHF(nBas,nBas2,nO,e,Sao,C,P,ENuc,ET,EV,EJ,EK,EHF,dipole)

! Print one-electron energies and other stuff for GHF

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nBas2
  integer,intent(in)                 :: nO
  double precision,intent(in)        :: e(nBas2)
  ! TODO
  ! add AO overlap as input
  double precision,intent(in)        :: Sao(nBas,nBas)
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

  integer                            :: i, j, ixyz
  integer                            :: mu,nu
  integer                            :: HOMO
  integer                            :: LUMO
  double precision                   :: Gap
  double precision                   :: Sz,Sx2,Sy2,Sz2,S2
  double precision                   :: na, nb
  double precision                   :: nonco_z, contam_uhf, xy_perp, contam_ghf

  double precision,allocatable       :: Ca(:,:)
  double precision,allocatable       :: Cb(:,:)
  double precision,allocatable       :: Paa(:,:), Saa(:,:)
  double precision,allocatable       :: Pab(:,:), Sab(:,:)
  double precision,allocatable       :: Pba(:,:), Sba(:,:)
  double precision,allocatable       :: Pbb(:,:), Sbb(:,:)
  double precision,allocatable       :: tmp(:,:)

  double precision,external          :: trace_matrix

  logical                            :: dump_orb = .false.

! HOMO and LUMO

  HOMO = nO
  LUMO = HOMO + 1
  Gap = e(LUMO)-e(HOMO)

! Density matrices

  allocate(Paa(nBas2,nBas2),Pab(nBas2,nBas2),Pba(nBas2,nBas2),Pbb(nBas2,nBas2))

  Paa(:,:) = P(     1:nBas ,     1:nBas )
  Pab(:,:) = P(     1:nBas ,nBas+1:nBas2)
  Pba(:,:) = P(nBas+1:nBas2,     1:nBas )
  Pbb(:,:) = P(nBas+1:nBas2,nBas+1:nBas2)

  ! TODO
  ! check C size
  allocate(Ca(nBas,nBas), Cb(nBas,nBas))
  do i = 1, nBas
    do j = 1, nBas
      Ca(j,i) = C(j,     i)
      Cb(j,i) = C(j,nBas+i)
    enddo
  enddo

  allocate(Saa(nBas,nBas),Sab(nBas,nBas),Sba(nBas,nBas),Sbb(nBas,nBas))
  allocate(tmp(nBas,nBas))

  ! Saa = Ca x Sao x Ca.T
  call dgemm("N", "N", nBas, nBas, nBas, 1.d0,  Ca,  size(Ca, 1), Sao, size(Sao, 1), 0.d0, tmp, size(tmp, 1))
  call dgemm("N", "T", nBas, nBas, nBas, 1.d0, tmp, size(tmp, 1),  Ca,  size(Ca, 1), 0.d0, Saa, size(Saa, 1))

  ! Sab = Ca x Sao x Cb.T
  call dgemm("N", "N", nBas, nBas, nBas, 1.d0,  Ca,  size(Ca, 1), Sao, size(Sao, 1), 0.d0, tmp, size(tmp, 1))
  call dgemm("N", "T", nBas, nBas, nBas, 1.d0, tmp, size(tmp, 1),  Cb,  size(Cb, 1), 0.d0, Sab, size(Sab, 1))

  ! Sba = Cb x Sao x Ca.T
  !     = Sab.T
  Sba = transpose(Sab)

  ! Sbb = Cb x Sao x Cb.T
  call dgemm("N", "N", nBas, nBas, nBas, 1.d0,  Cb,  size(Cb, 1), Sao, size(Sao, 1), 0.d0, tmp, size(tmp, 1))
  call dgemm("N", "T", nBas, nBas, nBas, 1.d0, tmp, size(tmp, 1),  Cb,  size(Cb, 1), 0.d0, Sbb, size(Sbb, 1))

  deallocate(tmp)

  ! TODO
  ! nO = nb of electrons ?
  na = 0.d0
  nb = 0.d0
  do i = 1, nO
    na = na + Saa(i,i)
    nb = nb + Sbb(i,i)
  enddo
  
  nonco_z = dble(nO)
  do j = 1, nO
    do i = 1, nO
      nonco_z = nonco_z - (Saa(i,j) - Sbb(i,j))**2
    enddo
  enddo
  nonco_z = 0.25d0 * nonco_z

  Sz  = 0.5d0 * (na - nb)
  Sz2 = Sz*Sz + nonco_z

  ! If Na > Nb
  !contam_uhf = nb
  !do j = 1, nO
  !  do i = 1, nO
  !    contam_uhf = contam_uhf - (Sab(i,j) - Sba(j,i))
  !  enddo
  !enddo
  !xy_perp = 0.d0
  !do i = 1, nO
  !  xy_perp = xy_perp + (Sba(i,i))**2
  !enddo
  !S2 = Sz * (Sz + 1.d0) + nonco_z + contam_uhf + xy_perp

  contam_ghf = 0.d0
  do j = 1, nO
    do i = 1, nO
      contam_ghf = contam_ghf - (Sab(i,i)*Sba(j,j) - Sab(i,j)*Sba(j,i))
    enddo
  enddo
  S2 = Sz * (Sz + 1.d0) + nonco_z + contam_ghf
  

! Compute expectation values of S^2 (WRONG!)

!  Sx2 = 0.25d0*trace_matrix(nBas,Paa+Pbb) + 0.25d0*trace_matrix(nBas,Pab+Pba)**2
!  do mu=1,nBas
!    do nu=1,nBas
!        Sx2 = Sx2 - 0.5d0*(Paa(mu,nu)*Pbb(nu,mu) + Pab(mu,nu)*Pab(nu,mu))
!    end do
!  end do
!
!  Sy2 = 0.25d0*trace_matrix(nBas,Paa+Pbb) - 0.25d0*trace_matrix(nBas,Pab+Pba)**2
!  do mu=1,nBas
!    do nu=1,nBas
!        Sy2 = Sy2 - 0.5d0*(Paa(mu,nu)*Pbb(nu,mu) - Pab(mu,nu)*Pab(nu,mu))
!    end do
!  end do
!
!  Sz2 = 0.25d0*trace_matrix(nBas,Paa+Pbb) + 0.25d0*trace_matrix(nBas,Pab-Pba)**2
!  do mu=1,nBas
!    do nu=1,nBas
!        Sz2 = Sz2 - 0.25d0*(Paa(mu,nu)*Pbb(nu,mu) - Pab(mu,nu)*Pab(nu,mu))
!        Sz2 = Sz2 + 0.25d0*(Pab(mu,nu)*Pba(nu,mu) - Pba(mu,nu)*Pab(nu,mu))
!    end do
!  end do
!  
!  S2 = Sx2 + Sy2 + Sz2

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
