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
  double precision                   :: Sx ,Sy ,Sz
  double precision                   :: Sx2,Sy2,Sz2
  double precision                   :: SmSp,SpSm,S2
  double precision                   :: na, nb
  double precision                   :: nonco_z, contam_uhf, xy_perp, contam_ghf

  double precision,allocatable       :: Ca(:,:)
  double precision,allocatable       :: Cb(:,:)
  double precision,allocatable       :: Paa(:,:), Saa(:,:)
  double precision,allocatable       :: Pab(:,:), Sab(:,:)
  double precision,allocatable       :: Pba(:,:), Sba(:,:)
  double precision,allocatable       :: Pbb(:,:), Sbb(:,:)
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
  SpSm = 0.5d0*(trace_matrix(nO,Paa) + SpSm)

! Sx2 = 0.25d0*trace_matrix(nO,Paa+Pbb) + 0.25d0*trace_matrix(nO,Pab+Pba)**2 &
!     - 0.5d0*trace_matrix(nO,matmul(Paa,Pbb) + matmul(Pab,Pab))

  SmSp = 0d0
  do i=1,nO
    do j=1,nO
      SmSp = SmSp + Pba(i,i)*Pab(j,j) - Pba(i,j)*Pab(j,i)
    end do
  end do
  SmSp = 0.5d0*(trace_matrix(nO,Pbb) + SmSp)

! Sy2 = 0.25d0*trace_matrix(nO,Paa+Pbb) - 0.25d0*trace_matrix(nO,Pab-Pba)**2 &
!     - 0.5d0*trace_matrix(nO,matmul(Paa,Pbb) - matmul(Pab,Pab))

  Sz2 = 0d0
  do i=1,nO
    do j=1,nO
      Sz2 = Sz2 + (Paa(i,i) - Pbb(i,i))*(Paa(j,j) - Pbb(j,j)) - (Paa(i,j) - Pbb(i,j))**2
    end do
  end do
  Sz2 = 0.25d0*(dble(nO) + Sz2)
  print*,'<Sz^2> = ',Sz2

  Sz2 = 0.25d0*trace_matrix(nO,Paa+Pbb)    &
      - 0.25d0*trace_matrix(nO,Paa-Pbb)**2 &
      - 0.25d0*trace_matrix(nO,matmul(Paa-Pbb,Paa-Pbb))

  S2 = Sz*(Sz+1d0) + trace_matrix(nO,Pbb) + 0.25d0*trace_matrix(nO,Paa+Pbb)
  do i=1,nO
    do j=1,nO
      S2 = S2 - 0.25d0*(Paa(i,j) - Pbb(i,j))**2 &
              + (Pba(i,i)*Pab(j,j) - Pba(i,j)*Pab(j,i))
    end do   
  end do   
  print*,'<S^2> = ',S2

  Sx2 = 0.5d0*(S2 - Sz2 + 0.5d0*(SmSp + SpSm))
  print*,'<Sx^2> = ',Sx2
  Sy2 = 0.5d0*(S2 - Sz2 - 0.5d0*(SmSp + SpSm))
  print*,'<Sy^2> = ',Sy2
!  Sx2 = 0.25d0*trace_matrix(nO,Paa+Pbb) + 0.25d0*trace_matrix(nO,Pab+Pba)**2 &
!      - 0.5d0*trace_matrix(nO,matmul(Paa,Pbb) + matmul(Pab,Pab))

!  Sx2 = trace_matrix(

!  Sy2 = 0.25d0*trace_matrix(nO,Paa+Pbb) - 0.25d0*trace_matrix(nO,Pab-Pba)**2 &
!      - 0.5d0*trace_matrix(nO,matmul(Paa,Pbb) - matmul(Pab,Pab))


!  Sz2 = 0.25d0*trace_matrix(nO,Paa+Pbb) + 0.25d0*trace_matrix(nO,Paa-Pbb)**2 &
!      - 0.25d0*trace_matrix(nO,matmul(Paa,Paa) + matmul(Pbb,Pbb)) &
!      + 0.25d0*trace_matrix(nO,matmul(Pab,Pba) + matmul(Pba,Pab))

!  S2 = Sz*(Sz+1d0) + trace_matrix(nO,Pbb) + 0.25d0*trace_matrix(nO,Paa+Pbb)

!  do i=1,nO
!    do j=1,nO
!      S2 = S2 - 0.25d0*(Paa(i,j) - Pbb(i,j))**2 &
!              + (Pba(i,i)*Pab(j,j) - Pba(i,j)*Pab(j,i))
!    end do   
!  end do   
! print*,'<S^2> = ',S2

  ! TODO
  ! check C size
  allocate(Ca(nBas,nBas), Cb(nBas,nBas))
  do i = 1, nBas
    do j = 1, nBas
      Ca(j,i) = C(j,     i)
      Cb(j,i) = C(j,nBas+i)
    enddo
  enddo

! allocate(Saa(nBas,nBas),Sab(nBas,nBas),Sba(nBas,nBas),Sbb(nBas,nBas))
! allocate(tmp(nBas,nBas))

  ! Saa = Ca x Sao x Ca.T
! call dgemm("N", "N", nBas, nBas, nBas, 1.d0,  Ca,  size(Ca, 1), Sao, size(Sao, 1), 0.d0, tmp, size(tmp, 1))
! call dgemm("N", "T", nBas, nBas, nBas, 1.d0, tmp, size(tmp, 1),  Ca,  size(Ca, 1), 0.d0, Saa, size(Saa, 1))

  ! Sab = Ca x Sao x Cb.T
! call dgemm("N", "N", nBas, nBas, nBas, 1.d0,  Ca,  size(Ca, 1), Sao, size(Sao, 1), 0.d0, tmp, size(tmp, 1))
! call dgemm("N", "T", nBas, nBas, nBas, 1.d0, tmp, size(tmp, 1),  Cb,  size(Cb, 1), 0.d0, Sab, size(Sab, 1))

  ! Sba = Cb x Sao x Ca.T
  !     = Sab.T
! Sba = transpose(Sab)

  ! Sbb = Cb x Sao x Cb.T
! call dgemm("N", "N", nBas, nBas, nBas, 1.d0,  Cb,  size(Cb, 1), Sao, size(Sao, 1), 0.d0, tmp, size(tmp, 1))
! call dgemm("N", "T", nBas, nBas, nBas, 1.d0, tmp, size(tmp, 1),  Cb,  size(Cb, 1), 0.d0, Sbb, size(Sbb, 1))

! deallocate(tmp)

  ! TODO
  ! nO = nb of electrons ?
! na = 0.d0
! nb = 0.d0
! do i = 1, nO
!   na = na + Saa(i,i)
!   nb = nb + Sbb(i,i)
! enddo
  
! nonco_z = dble(nO)
! do j = 1, nO
!   do i = 1, nO
!     nonco_z = nonco_z - (Saa(i,j) - Sbb(i,j))**2
!   enddo
! enddo
! nonco_z = 0.25d0 * nonco_z

! Sz  = 0.5d0 * (na - nb)
! Sz2 = Sz*Sz + nonco_z

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

! contam_ghf = 0.d0
! do j = 1, nO
!   do i = 1, nO
!     contam_ghf = contam_ghf - (Sab(i,i)*Sba(j,j) - Sab(i,j)*Sba(j,i))
!   enddo
! enddo
! S2 = Sz * (Sz + 1.d0) + nonco_z + contam_ghf
  



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
  write(*,'(A33,1X,F16.6)')     ' <Sx**2>             = ',Sx2
  write(*,'(A33,1X,F16.6)')     ' <Sy**2>             = ',Sy2
  write(*,'(A33,1X,F16.6)')     ' <Sz**2>             = ',Sz2
  write(*,'(A33,1X,F16.6)')     ' <S**2>              = ',Sx2+Sy2+Sz2
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

end subroutine 
