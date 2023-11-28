subroutine print_GHF_spin(nBas,nBas2,nO,C,S)

  implicit none

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nBas2
  integer,intent(in)            :: nO
  double precision, intent(in)  :: C(nBas2,nBas2)
  double precision, intent(in)  :: S(nBas,nBas)

  integer                       :: i, j
  double precision              :: Na, Nb
  double precision              :: nonco_z, contam_ghf
  double precision              :: S2, Sz, Sz2
  complex*16                    :: Sc_x, Sc_y, Sc_z
  complex*16                    :: Sc_xx, Sc_xy, Sc_xz
  complex*16                    :: Sc_yx, Sc_yy, Sc_yz
  complex*16                    :: Sc_zx, Sc_zy, Sc_zz
  double precision, allocatable :: Ca(:,:), Cb(:,:)
  double precision, allocatable :: Paa(:,:), Pab(:,:), Pba(:,:), Pbb(:,:)
  double precision, allocatable :: Mc(:,:), Eigc(:)

  write(*,*)  
  write(*,*) '****************************************'
  write(*,*) '* Spin properties of GHF wave function *'
  write(*,*) '****************************************'
  write(*,*)  

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
      contam_ghf = contam_ghf + (Pab(i,i)*Pba(j,j) - Pab(i,j)*Pba(j,i))
    enddo
  enddo

  Sz  = 0.5d0 * (Na - Nb)
  Sz2 = Sz*Sz + nonco_z
  S2  = Sz2 + 0.5d0 * (Na + Nb) + contam_ghf

  write(*,'(A15,2F10.6)') ' < Sz > = ', Sz
  write(*,'(A15,2F10.6)') ' < Sz^2 > = ', Sz2
  write(*,'(A15,2F10.6)') ' < S^2 > = ', S2
  write(*,*)

  ! --- --- --- --- --- --- --- --- ---
  ! Compute <Si> & <SiSj> for all i, j
  ! --- --- --- --- --- --- --- --- ---

  Sc_x = (0.d0,0.d0)
  Sc_y = (0.d0,0.d0)
  Sc_z = (0.d0,0.d0)
  do i = 1, nO
    Sc_x = Sc_x + (+0.5d0,0.d0) * (Pab(i,i) + Pba(i,i))
    Sc_y = Sc_y + (0.d0,-0.5d0) * (Pab(i,i) - Pba(i,i))
    Sc_z = Sc_z + (+0.5d0,0.d0) * (Paa(i,i) - Pbb(i,i))
  enddo
  write(*,'(A15,2F10.6)') ' < Sx > = ',Sc_x
  write(*,'(A15,2F10.6)') ' < Sy > = ',Sc_y
  write(*,'(A15,2F10.6)') ' < Sz > = ',Sc_z
  write(*,*)

  Sc_xx = Sc_x * Sc_x + 0.25d0*dble(nO)*(1.d0,0.d0)
  Sc_yy = Sc_y * Sc_y + 0.25d0*dble(nO)*(1.d0,0.d0)
  Sc_zz = Sc_z * Sc_z + 0.25d0*dble(nO)*(1.d0,0.d0)
  do i = 1, nO
    do j = 1, nO
      Sc_xx = Sc_xx - zabs((+0.5d0,0.d0) * (Pab(i,j) + Pba(i,j)))**2
      Sc_yy = Sc_yy - zabs((0.d0,-0.5d0) * (Pab(i,j) - Pba(i,j)))**2
      Sc_zz = Sc_zz - zabs((+0.5d0,0.d0) * (Paa(i,j) - Pbb(i,j)))**2
    enddo
  enddo
  write(*,'(A15,2F10.6)') ' < Sx^2 > = ',Sc_xx
  write(*,'(A15,2F10.6)') ' < Sy^2 > = ',Sc_yy
  write(*,'(A15,2F10.6)') ' < Sz^2 > = ',Sc_zz
  write(*,*)

  Sc_xy = Sc_x * Sc_y
  Sc_yx = Sc_x * Sc_y
  do i = 1, nO
    Sc_xy = Sc_xy + (0.d0,0.5d0) * (+0.5d0,0.0d0) * (Paa(i,i) - Pbb(i,i))
    Sc_yx = Sc_yx - (0.d0,0.5d0) * (+0.5d0,0.0d0) * (Paa(i,i) - Pbb(i,i))
    do j = 1, nO
      Sc_xy = Sc_xy - (+0.5d0,0.d0) * (Pab(i,j) + Pba(i,j)) * (0.d0,-0.5d0) * (Pab(j,i) - Pba(j,i))
      Sc_yx = Sc_yx - (+0.5d0,0.d0) * (Pab(j,i) + Pba(j,i)) * (0.d0,-0.5d0) * (Pab(i,j) - Pba(i,j))
    enddo
  enddo
  write(*,'(A15,2F10.6)') ' < Sx.Sy > = ',Sc_xy
  write(*,'(A15,2F10.6)') ' < Sy.Sx > = ',Sc_yx

  Sc_xz = Sc_x * Sc_z
  Sc_zx = Sc_x * Sc_z
  do i = 1, nO
    Sc_xz = Sc_xz - (0.d0,0.5d0) * (0.d0,-0.5d0) * (Pab(i,i) - Pba(i,i))
    Sc_zx = Sc_zx + (0.d0,0.5d0) * (0.d0,-0.5d0) * (Pab(i,i) - Pba(i,i))
    do j = 1, nO
      Sc_xz = Sc_xz - (+0.5d0,0.d0) * (Pab(i,j) + Pba(i,j)) * (+0.5d0,0.d0) * (Paa(j,i) - Pbb(j,i))
      Sc_zx = Sc_zx - (+0.5d0,0.d0) * (Pab(j,i) + Pba(j,i)) * (+0.5d0,0.d0) * (Paa(i,j) - Pbb(i,j))
    enddo
  enddo
  write(*,'(A15,2F10.6)') ' < Sx.Sz > = ',Sc_xz
  write(*,'(A15,2F10.6)') ' < Sz.Sx > = ',Sc_zx

  Sc_yz = Sc_y * Sc_z
  Sc_zy = Sc_y * Sc_z
  do i = 1, nO
    Sc_yz = Sc_yz + (0.d0,0.5d0) * (+0.5d0,0.d0) * (Pab(i,i) + Pba(i,i))
    Sc_zy = Sc_zy - (0.d0,0.5d0) * (+0.5d0,0.d0) * (Pab(i,i) + Pba(i,i))
    do j = 1, nO
      Sc_yz = Sc_yz - (0.d0,-0.5d0) * (Pab(i,j) - Pba(i,j)) * (+0.5d0,0.d0) * (Paa(j,i) - Pbb(j,i))
      Sc_zy = Sc_zy - (0.d0,-0.5d0) * (Pab(j,i) - Pba(j,i)) * (+0.5d0,0.d0) * (Paa(i,j) - Pbb(i,j))
    enddo
  enddo
  write(*,'(A15,2F10.6)') ' < Sy.Sz > = ',Sc_yz
  write(*,'(A15,2F10.6)') ' < Sz.Sy > = ', Sc_zy
  write(*,*)  


  ! --- --- --- --- --- !
  !  Collinearity Test  !
  ! --- --- --- --- --- !

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
      Mc(1,3) = Mc(1,3) - 0.25d0 * (Pab(i,j) + Pba(i,j))*(Paa(j,i) - Pbb(j,i))
    enddo
  enddo
  Mc(3,1) = Mc(1,3)

  write(*,*) 'The collinearity matrix is'
  call matout(3,3,Mc)

  call diagonalize_matrix(3,Mc,Eigc)
  write(*,*)
  write(*,'(A40,3F10.6)') 'Eigenvalues of collinearity matrix:', Eigc
  write(*,'(A40,1F10.6)') 'Smallest eigenvalue:',Eigc(1)
  write(*,'(A40)')        '(0 iff wave function collinear)'
  deallocate(Mc,Eigc)

  ! --- --- --- --- --- --- --- --- ---
  ! --- --- --- --- --- --- --- --- ---

  deallocate(Paa, Pab, Pba, Pbb)

end
