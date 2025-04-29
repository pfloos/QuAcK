subroutine diagonalize_nonsym_matrix(N, A, L, e_re, thr_d, thr_nd, thr_deg, imp_bio, verbose)

  ! Diagonalize a non-symmetric matrix A
  ! 
  ! Output
  !   right-eigenvectors are saved in A
  !    left-eigenvectors are saved in L
  !         eigenvalues  are saved in e = e_re + i e_im

  implicit none

  integer,          intent(in)    :: N
  logical,          intent(in)    :: imp_bio, verbose
  double precision, intent(in)    :: thr_d, thr_nd, thr_deg
  double precision, intent(inout) :: A(N,N)
  double precision, intent(out)   :: e_re(N), L(N,N)

  integer                         :: i, j, ii
  integer                         :: lwork, info
  double precision                :: accu_d, accu_nd
  integer,          allocatable   :: iorder(:), deg_num(:)
  double precision, allocatable   :: Atmp(:,:), Ltmp(:,:), work(:), e_im(:)
  double precision, allocatable   :: S(:,:)

  if(verbose) then
    print*, ' Starting a non-Hermitian diagonalization ...'
    print*, ' Good Luck ;)'
    print*, ' imp_bio = ', imp_bio
  endif

  ! ---
  ! diagonalize

  allocate(Atmp(N,N), e_im(N))
  Atmp(1:N,1:N) = A(1:N,1:N)

  allocate(work(1))
  lwork = -1
  call dgeev('V', 'V', N, Atmp, N, e_re, e_im, L, N, A, N, work, lwork, info)
  if(info .gt. 0) then
    print*,'dgeev failed !!', info
    stop
  endif

  lwork = max(int(work(1)), 1)
  deallocate(work)
  allocate(work(lwork))

  call dgeev('V', 'V', N, Atmp, N, e_re, e_im, L, N, A, N, work, lwork, info)
  if(info .ne. 0) then
    print*,'dgeev failed !!', info
    stop
  endif

  deallocate(Atmp, WORK)


  ! ---
  ! check if eigenvalues are real

  i  = 1
  ii = 0
  do while(i .le. N)
    if(dabs(e_im(i)) .gt. 1.d-12) then
      ii = ii + 1
      if(verbose) then
        print*, ' Warning: complex eigenvalue !'
        print*, i, e_re(i), e_im(i)
        if(dabs(e_im(i)/e_re(i)) .lt. 1.d-6) then
          print*, ' small enouph to be igored'
        else
          print*, ' IMAGINARY PART IS SIGNIFANT !!!'
        endif
      endif
    endif
    i = i + 1
  enddo

  if(verbose) then
    if(ii .eq. 0) print*, ' congratulations :) eigenvalues are real-valued !!'
  endif


  ! ---
  ! track & sort the real eigenvalues 

  allocate(Atmp(N,N), Ltmp(N,N), iorder(N))

  do i = 1, N
    iorder(i) = i
  enddo
  call quick_sort(e_re, iorder, N)

  Atmp(:,:) = A(:,:)
  Ltmp(:,:) = L(:,:)
  do i = 1, N
    do j = 1, N
      A(j,i) = Atmp(j,iorder(i))
      L(j,i) = Ltmp(j,iorder(i))
    enddo
  enddo

  deallocate(Atmp, Ltmp, iorder)




  ! ---
  ! check bi-orthog

  allocate(S(N,N))
  call check_biorthog(N, N, L, A, accu_d, accu_nd, S, thr_d, thr_nd, .false., verbose)

  if((accu_nd .lt. thr_nd) .and. (dabs(accu_d-dble(N))/dble(N) .lt. thr_d)) then

    if(verbose) then
      print *, ' lapack vectors are normalized and bi-orthogonalized'
    endif

  elseif((accu_nd .lt. thr_nd) .and. (dabs(accu_d - dble(N)) .gt. thr_d)) then

    if(verbose) then
      print *, ' lapack vectors are not normalized but bi-orthogonalized'
    endif

    call check_biorthog_binormalize(N, N, L, A, thr_d, thr_nd, .true.)
    call check_biorthog(N, N, L, A, accu_d, accu_nd, S, thr_d, thr_nd, .true., verbose)

  else

    if(verbose) then
      print *, ' lapack vectors are not normalized neither bi-orthogonalized'
    endif

    allocate(deg_num(N))
    call reorder_degen_eigvec(N, thr_deg, deg_num, e_re, L, A)
    call impose_biorthog_degen_eigvec(N, deg_num, e_re, L, A)
    deallocate(deg_num)

    call check_biorthog(N, N, L, A, accu_d, accu_nd, S, thr_d, thr_nd, .false., verbose)
    if((accu_nd .lt. thr_nd) .and. (dabs(accu_d-dble(N))/dble(N) .lt. thr_d)) then
      if(verbose) then
        print *, ' lapack vectors are now normalized and bi-orthogonalized'
      endif
    elseif((accu_nd .lt. thr_nd) .and. (dabs(accu_d - dble(N)) .gt. thr_d)) then
      if(verbose) then
        print *, ' lapack vectors are now not normalized but bi-orthogonalized'
      endif
      call check_biorthog_binormalize(N, N, L, A, thr_d, thr_nd, .true.)
      call check_biorthog(N, N, L, A, accu_d, accu_nd, S, thr_d, thr_nd, .true., verbose)
    else
      if(verbose) then
        print*, ' bi-orthogonalization failed !'
      endif
      if(imp_bio) then
        print*, ' bi-orthogonalization failed !'
        deallocate(S)
        stop
      endif
    endif

  endif

  deallocate(S)
  return

end

! ---

subroutine check_biorthog(n, m, Vl, Vr, accu_d, accu_nd, S, thr_d, thr_nd, stop_ifnot, verbose)

  implicit none

  integer,          intent(in)  :: n, m
  logical,          intent(in)  :: stop_ifnot, verbose
  double precision, intent(in)  :: Vl(n,m), Vr(n,m)
  double precision, intent(in)  :: thr_d, thr_nd
  double precision, intent(out) :: accu_d, accu_nd, S(m,m)

  integer                       :: i, j
  double precision, allocatable :: SS(:,:)

  if(verbose) then
    print *, ' check bi-orthogonality'
  endif

  ! ---

  call dgemm( 'T', 'N', m, m, n, 1.d0          &
            , Vl, size(Vl, 1), Vr, size(Vr, 1) &
            , 0.d0, S, size(S, 1) )

  accu_d  = 0.d0
  accu_nd = 0.d0
  do i = 1, m
    do j = 1, m
      if(i==j) then
        accu_d = accu_d + dabs(S(i,i))
      else
        accu_nd = accu_nd + S(j,i) * S(j,i)
      endif
    enddo
  enddo
  accu_nd = dsqrt(accu_nd) / dble(m)

  if(verbose) then
    if((accu_nd .gt. thr_nd) .or. dabs(accu_d-dble(m))/dble(m) .gt. thr_d) then
      print *, ' non bi-orthogonal vectors !'
      print *, ' accu_nd = ', accu_nd
      print *, ' accu_d  = ', dabs(accu_d-dble(m))/dble(m)
    else
      print *, ' vectors are bi-orthogonals'
    endif
  endif

  ! ---

  if(stop_ifnot .and. ((accu_nd .gt. thr_nd) .or. dabs(accu_d-dble(m))/dble(m) .gt. thr_d)) then
    print *, ' non bi-orthogonal vectors !'
    print *, ' accu_nd = ', accu_nd
    print *, ' accu_d  = ', dabs(accu_d-dble(m))/dble(m)
    stop
  endif

end

! ---

subroutine check_biorthog_binormalize(n, m, Vl, Vr, thr_d, thr_nd, stop_ifnot)

  implicit none

  integer,          intent(in)    :: n, m
  logical,          intent(in)    :: stop_ifnot
  double precision, intent(in)    :: thr_d, thr_nd
  double precision, intent(inout) :: Vl(n,m), Vr(n,m)

  integer                         :: i, j
  double precision                :: accu_d, accu_nd, s_tmp
  double precision, allocatable   :: S(:,:)

  ! ---

  allocate(S(m,m))
  call dgemm( 'T', 'N', m, m, n, 1.d0          &
            , Vl, size(Vl, 1), Vr, size(Vr, 1) &
            , 0.d0, S, size(S, 1) )

  do i = 1, m
    if(S(i,i) .lt. 0.d0) then
      do j = 1, n
        Vl(j,i) = -1.d0 * Vl(j,i)
      enddo
      S(i,i) = -S(i,i)
    endif
  enddo

  accu_d  = 0.d0
  accu_nd = 0.d0
  do i = 1, m
    do j = 1, m
      if(i==j) then
        accu_d = accu_d + S(i,i)
      else
        accu_nd = accu_nd + S(j,i) * S(j,i)
      endif
    enddo
  enddo
  accu_nd = dsqrt(accu_nd) / dble(m)

  ! ---

  if( (accu_nd .lt. thr_nd) .and. (dabs(accu_d-dble(m))/dble(m) .gt. thr_d) ) then

    do i = 1, m
      if(S(i,i) <= 0.d0) then
        print *, ' negative overlap !'
        print *, i, S(i,i)
        exit
      endif
      if(dabs(S(i,i) - 1.d0) .gt. thr_d) then
        s_tmp = 1.d0 / dsqrt(S(i,i))
        do j = 1, n
          Vl(j,i) = Vl(j,i) * s_tmp
          Vr(j,i) = Vr(j,i) * s_tmp
        enddo
      endif

    enddo

  endif

  ! ---

  call dgemm( 'T', 'N', m, m, n, 1.d0          &
            , Vl, size(Vl, 1), Vr, size(Vr, 1) &
            , 0.d0, S, size(S, 1) )

  accu_d  = 0.d0
  accu_nd = 0.d0
  do i = 1, m
    do j = 1, m
      if(i==j) then
        accu_d = accu_d + S(i,i)
      else
        accu_nd = accu_nd + S(j,i) * S(j,i)
      endif
    enddo
  enddo
  accu_nd = dsqrt(accu_nd) / dble(m)

  deallocate(S)

  ! ---

  if( stop_ifnot .and. ((accu_nd .gt. thr_nd) .or. (dabs(accu_d-dble(m))/dble(m) .gt. thr_d)) ) then
    print *, accu_nd, thr_nd
    print *, dabs(accu_d-dble(m))/dble(m), thr_d
    print *, ' biorthog_binormalize failed !'
    stop
  endif

end

! ---

subroutine reorder_degen_eigvec(n, thr_deg, deg_num, e0, L0, R0)

  implicit none

  integer,          intent(in)    :: n
  double precision, intent(in)    :: thr_deg
  double precision, intent(inout) :: e0(n), L0(n,n), R0(n,n)
  integer,          intent(out)   :: deg_num(n)

  logical                         :: complex_root
  integer                         :: i, j, k, m, ii, j_tmp
  double precision                :: ei, ej, de
  double precision                :: accu_d, accu_nd
  double precision                :: e0_tmp, L0_tmp(n), R0_tmp(n)
  double precision, allocatable   :: L(:,:), R(:,:), S(:,:), S_inv_half(:,:)

  do i = 1, n
    deg_num(i) = 1
  enddo

  do i = 1, n-1
    ei = e0(i)

    ! already considered in degen vectors
    if(deg_num(i) .eq. 0) cycle

    ii = 0
    do j = i+1, n
      ej = e0(j)
      de = dabs(ei - ej)

      if(de .lt. thr_deg) then
        ii = ii + 1

        j_tmp = i + ii

        deg_num(j_tmp) = 0

        e0_tmp    = e0(j_tmp)
        e0(j_tmp) = e0(j)
        e0(j)     = e0_tmp

        L0_tmp(1:n)   = L0(1:n,j_tmp)
        L0(1:n,j_tmp) = L0(1:n,j)
        L0(1:n,j)     = L0_tmp(1:n)

        R0_tmp(1:n)   = R0(1:n,j_tmp)
        R0(1:n,j_tmp) = R0(1:n,j)
        R0(1:n,j)     = R0_tmp(1:n)
      endif
    enddo

    deg_num(i) = ii + 1
  enddo

  ii = 0
  do i = 1, n
    if(deg_num(i) .gt. 1) then
      ii = ii + 1
    endif
  enddo

  if(ii .eq. 0) then
    print*, ' WARNING: bi-orthogonality is lost but there is no degeneracies'
    print*, ' rotations may change energy'
    stop
  endif

end

! ---

subroutine impose_biorthog_degen_eigvec(n, deg_num, e0, L0, R0)

  implicit none

  integer,          intent(in)    :: n, deg_num(n)
  double precision, intent(in)    :: e0(n)
  double precision, intent(inout) :: L0(n,n), R0(n,n)

  logical                         :: complex_root
  integer                         :: i, j, k, m
  double precision                :: ei, ej, de
  double precision                :: accu_d, accu_nd
  double precision, allocatable   :: L(:,:), R(:,:), S(:,:), S_inv_half(:,:)

  !do i = 1, n
  !  if(deg_num(i) .gt. 1) then
  !    print *, ' degen on', i, deg_num(i), e0(i)
  !  endif
  !enddo

  ! ---

  do i = 1, n
    m = deg_num(i)

    if(m .gt. 1) then

      allocate(L(n,m), R(n,m), S(m,m))

      do j = 1, m
        L(1:n,j) = L0(1:n,i+j-1)
        R(1:n,j) = R0(1:n,i+j-1)
      enddo

      ! ---

      call dgemm( 'T', 'N', m, m, n, 1.d0      &
                , L, size(L, 1), R, size(R, 1) &
                , 0.d0, S, size(S, 1) )

      accu_nd = 0.d0
      do j = 1, m
        do k = 1, m
          if(j==k) cycle
          accu_nd = accu_nd + dabs(S(j,k))
        enddo
      enddo

      if(accu_nd .lt. 1d-12) then
        deallocate(S, L, R)
        cycle
      endif

      call impose_biorthog_svd(n, m, L, R)

      call dgemm( 'T', 'N', m, m, n, 1.d0      &
                , L, size(L, 1), R, size(R, 1) &
                , 0.d0, S, size(S, 1) )
      accu_nd = 0.d0
      do j = 1, m
        do k = 1, m
          if(j==k) cycle
          accu_nd = accu_nd + dabs(S(j,k))
        enddo
      enddo
      if(accu_nd .gt. 1d-7) then
        print*, ' accu_nd =', accu_nd
        print*, ' your strategy for degenerates orbitals failed !'
        print*, m, 'deg on', i
        stop
      endif

      deallocate(S)

      ! ---

      do j = 1, m
        L0(1:n,i+j-1) = L(1:n,j)
        R0(1:n,i+j-1) = R(1:n,j)
      enddo

      deallocate(L, R)

    endif
  enddo

end

! ---

subroutine impose_biorthog_svd(n, m, L, R)

  implicit none

  integer,          intent(in)    :: n, m
  double precision, intent(inout) :: L(n,m), R(n,m)

  integer                         :: i, j, num_linear_dependencies
  double precision                :: threshold
  double precision, allocatable   :: S(:,:), tmp(:,:)
  double precision, allocatable   :: U(:,:), V(:,:), Vt(:,:), D(:)

  allocate(S(m,m))

  call dgemm( 'T', 'N', m, m, n, 1.d0      &
            , L, size(L, 1), R, size(R, 1) &
            , 0.d0, S, size(S, 1) )

  ! ---

  allocate(U(m,m), Vt(m,m), D(m))

  call svd_local(S, m, U, m, D, Vt, m, m, m)

  deallocate(S)

  threshold               = 1.d-6
  num_linear_dependencies = 0
  do i = 1, m
    if(dabs(D(i)) <= threshold) then
      D(i) = 0.d0
      num_linear_dependencies = num_linear_dependencies + 1
    else
      D(i) = 1.d0 / dsqrt(D(i))
    endif
  enddo
  if(num_linear_dependencies > 0) then
    write(*,*) ' linear dependencies = ', num_linear_dependencies
    write(*,*) ' m                   = ', m
    stop
  endif

  allocate(V(m,m))
  do i = 1, m
    do j = 1, m
      V(j,i) = Vt(i,j)
    enddo
  enddo
  deallocate(Vt)

  ! ---

  ! R <-- R x V x D^{-0.5}
  ! L <-- L x U x D^{-0.5}

  do i = 1, m
    do j = 1, m
      V(j,i) = V(j,i) * D(i)
      U(j,i) = U(j,i) * D(i)
    enddo
  enddo

  allocate(tmp(n,m))
  tmp(:,:) = R(:,:)
  call dgemm( 'N', 'N', n, m, m, 1.d0          &
            , tmp, size(tmp, 1), V, size(V, 1) &
            , 0.d0, R, size(R, 1))

  tmp(:,:) = L(:,:)
  call dgemm( 'N', 'N', n, m, m, 1.d0          &
            , tmp, size(tmp, 1), U, size(U, 1) &
            , 0.d0, L, size(L, 1))

  deallocate(tmp, U, V, D)

end

! ---

subroutine svd_local(A, LDA, U, LDU, D, Vt, LDVt, m, n)

  implicit none

  ! Compute A = U D Vt

  integer,          intent(in)  :: LDA, LDU, LDVt, m, n
  double precision, intent(in)  :: A(LDA,n)
  double precision, intent(out) :: U(LDU, min(m, n))
  double precision, intent(out) :: Vt(LDVt,n)
  double precision, intent(out) :: D(min(m, n))

  integer                       :: info, lwork, i, j, k
  double precision, allocatable :: work(:)

  double precision,allocatable  :: A_tmp(:,:)

  allocate(A_tmp(LDA,n))
  do k = 1, n
    do i = 1, m
      A_tmp(i,k) = A(i,k)
    enddo
  enddo

  ! Find optimal size for temp arrays
  allocate(work(1))
  lwork = -1
  call dgesvd('S','S', m, n, A_tmp, LDA, &
              D, U, LDU, Vt, LDVt, work, lwork, info)
  ! /!\ int(WORK(1)) becomes negative when WORK(1) > 2147483648
  lwork = max(int(work(1)), 10*MIN(M,N))
  deallocate(work)

  allocate(work(lwork))
  call dgesvd('S','S', m, n, A_tmp, LDA, &
              D, U, LDU, Vt, LDVt, work, lwork, info)

  deallocate(A_tmp, work)

  if (info /= 0) then
    print *,  info, ': SVD failed'
    stop
  endif

  do j = 1, min(m, n)
    do i = 1, m
      if (dabs(U(i,j)) < 1.d-14)  U(i,j) = 0.d0
    enddo
  enddo

  do j = 1, n
    do i = 1, n
      if(dabs(Vt(i,j)) < 1.d-14) Vt(i,j) = 0.d0
    enddo
  enddo

end
