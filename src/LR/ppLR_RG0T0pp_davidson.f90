
! ---

subroutine ppLR_RG0T0pp_davidson(ispin, TDA, nC, nO, nR, nOrb, nOO, nVV, lambda, e, eF, ERI, Om, R, n_states, n_states_diag)

  ! 
  ! Extract the low n_states 
  ! Om(i) (eigenvalues) and 
  ! R(:,i) (right-eigenvectors) 
  ! of the pp-RPA matrix
  !
  !   (+C    +B) 
  !   (        )  
  !   (-B.T  -D)
  !

  implicit none

  logical,          intent(in)  :: TDA
  integer,          intent(in)  :: ispin
  integer,          intent(in)  :: nC, nO, nR, nOrb, nOO, nVV
  integer,          intent(in)  :: n_states      ! nb of physical states
  integer,          intent(in)  :: n_states_diag ! nb of states used to get n_states
  double precision, intent(in)  :: lambda, eF
  double precision, intent(in)  :: e(nOrb)
  double precision, intent(in)  :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision, intent(out) :: Om(n_states)
  double precision, intent(out) :: R(nOO+nVV,n_states_diag)
  
  integer                       :: N, M
  integer                       :: iter, itermax, itertot
  integer                       :: shift1, shift2
  integer                       :: i, j, ij, k, l, kl, a, b, ab, c, d, cd
  integer                       :: i_omax(n_states)
  logical                       :: converged
  character*(16384)             :: write_buffer
  double precision              :: r1, r2, dtwo_pi
  double precision              :: lambda_tmp
  double precision              :: to_print(2,n_states)
  double precision              :: mem
  double precision, allocatable :: H_diag(:)
  double precision, allocatable :: W(:,:)
  double precision, allocatable :: U(:,:)
  double precision, allocatable :: h(:,:), h_vec(:,:), h_val(:)
  double precision, allocatable :: residual_norm(:)
  double precision, allocatable :: overlap(:)
  double precision, allocatable :: S_check(:,:)

  double precision, external    :: u_dot_u

  dtwo_pi = 6.283185307179586d0

  N = nOO + nVV
  itermax = 8
  M = n_states_diag * itermax

  if(M .ge. N) then
    print*, 'N = ', N
    print*, 'M = ', M
    print*, ' use Lapack or decrease n_states and/or itermax '
    stop
  endif

  write(6,'(A)') ''
  write(6,'(A)') 'Davidson Diagonalization'
  write(6,'(A)') '------------------------'
  write(6,'(A)') ''

  write(*,'(A40, I12)') 'Number of states = ', n_states
  write(*,'(A40, I12)') 'Number of states in diagonalization = ', n_states_diag
  write(*,'(A40, I12)') 'Number of basis functions = ', N




  allocate(H_diag(N))
  allocate(U(N,M))
  allocate(W(N,M))
  allocate(h(M,M), h_vec(M,M), h_val(M))
  allocate(overlap(n_states_diag))
  allocate(residual_norm(n_states_diag))

  mem = 8.d0 * dble(nOrb + nOrb**4 + N*n_states) / 1d6
  write(*,'(A40, F12.4)') 'I/O mem (MB) = ', mem

  mem = 8.d0 * dble(N + N*M + N*M + M*M + M*M + M + n_states_diag + n_states_diag) / 1d6
  write(*,'(A40, F12.4)') 'tmp mem (MB) = ', mem


  call ppLR_H_diag(ispin, nOrb, nC, nO, nR, nOO, nVV, lambda, e, eF, ERI, H_diag)

  !print*, "H_diag:"
  !do ab = 1, N
  !  print*, ab, H_diag(ab)
  !enddo

  ! initialize guess
  R = 0.d0
  do k = 1, n_states
    R(k,k) = 1.d0
  enddo
  do k = n_states+1, n_states_diag
    do i = 1, N
      call random_number(r1)
      call random_number(r2)
      r1 = dsqrt(-2.d0*dlog(r1))
      r2 = dtwo_pi*r2
      R(i,k) = r1*dcos(r2)
    enddo
    R(k,k) = R(k,k) + 10.d0
  enddo

  do k = 1, n_states_diag
    call normalize(R(1,k), N)
  enddo

  !print*, 'guess'
  !do k = 1, N
  !  write(*,'(100(F15.7,2X))') (R(k,i), i = 1, n_states_diag)
  !enddo

  ! working vectors
  do k = 1, n_states_diag
    do i = 1, N
      U(i,k) = R(i,k)
    enddo
  enddo

  !print*, 'working vectors'
  !do k = 1, N
  !  write(*,'(100(F15.7,2X))') (U(k,i), i = 1, n_states_diag)
  !enddo

  write(6,'(A)') ''
  write_buffer = '====='
  do i = 1, n_states
    write_buffer = trim(write_buffer)//' ================  ==========='
  enddo
  write(6,'(A)') write_buffer(1:6+41*n_states)
  write_buffer = 'Iter'
  do i = 1, n_states
    write_buffer = trim(write_buffer)//'       Energy         Residual '
  enddo
  write(6,'(A)') write_buffer(1:6+41*n_states)
  write_buffer = '====='
  do i = 1, n_states
    write_buffer = trim(write_buffer)//' ================  ==========='
  enddo
  write(6,'(A)') write_buffer(1:6+41*n_states)


  converged = .False.
  itertot = 0

  do while (.not.converged)

    itertot = itertot + 1
    if(itertot == itermax) then
      print*, 'exit before convergence !'
      print*, 'itertot == itermax', itertot
      exit
    endif

    do iter = 1, itermax-1

      shift1 = n_states_diag * (iter - 1)
      shift2 = shift1 + n_states_diag
      !print*, iter, shift1, shift2

      if((iter > 1) .or. (itertot == 1)) then

        call ortho_qr(U(1,1), size(U, 1), N, shift2)
        !call ortho_qr(U(1,1), size(U, 1), N, shift2)

        !print*, 'working vectors after qr'
        !do k = 1, N
        !  write(*,'(100(F15.7,2X))') (U(k,i), i = 1, n_states_diag)
        !enddo
        !allocate(S_check(shift2,shift2))
        !call dgemm("T", "N", shift2, shift2, N, 1.d0,      &
        !           U(1,1), size(U, 1), U(1,1), size(U, 1), &
        !           0.d0, S_check(1,1), size(S_check, 1))
        !do k = 1, shift2
        !  write(*,'(100(F15.7,2X))') (S_check(k,i), i = 1, shift2)
        !enddo
        !deallocate(S_check)

        call ppLR_HR_calc(ispin, nOrb, nC, nO, nR, nOO, nVV, lambda, e, eF, n_states_diag, &
                          ERI(1,1,1,1), U(1,shift1+1), W(1,shift1+1))

      else

        ! computed below
        continue
      endif

      ! h = U.T H U
      call dgemm('T', 'N', shift2, shift2, N, 1.d0,      &
                 U(1,1), size(U, 1), W(1,1), size(W, 1), &
                 0.d0, h(1,1), size(h, 1))

      ! h h_vec = h_val h_vec
      call diag_nonsym_right(shift2, h(1,1), size(h, 1), h_vec(1,1), size(h_vec, 1), h_val(1), size(h_val, 1))
      !print*, 'h_val', h_val(1:shift2)

      ! U1 = U0 h_vec
      call dgemm('N', 'N', N, n_states_diag, shift2, 1.d0,       &
                 U(1,1), size(U, 1), h_vec(1,1), size(h_vec, 1), &
                 0.d0, U(1,shift2+1), size(U, 1))

      do k = 1, n_states_diag
        call normalize(U(1,shift2+k), N)
      enddo

      do l = 1, n_states
        !do k = 1, n_states_diag
        !  overlap(k) = 0.d0
        !  do i = 1, N
        !    overlap(k) = overlap(k) + U(i,shift2+k) * R(i,l)
        !  enddo
        !  overlap(k) = dabs(overlap(k))
        !  !print *, ' overlap =', k, overlap(k)
        !enddo
        !lambda_tmp = 0.d0
        !do k = 1, n_states_diag
        !  if(overlap(k) .gt. lambda_tmp) then
        !    i_omax(l) = k
        !    lambda_tmp = overlap(k)
        !  endif
        !enddo
        !if(lambda_tmp .lt. 0.7d0) then
        !  print *, ' small overlap ...', l, i_omax(l)
        !  print *, ' max overlap =', lambda_tmp
        !  !stop
        !endif
        !if(i_omax(l) .ne. l) then
        !  print *, ' !!! WARNING !!!'
        !  print *, ' index of state', l, i_omax(l)
        !endif
      enddo

      ! W1 = W0 h_vec
      call dgemm('N', 'N', N, n_states_diag, shift2, 1.d0,       &
                 W(1,1), size(W, 1), h_vec(1,1), size(h_vec, 1), &
                 0.d0, W(1,shift2+1), size(W, 1))

      ! check if W1 = U1 h_val
      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP          PRIVATE(i, k) &
      !$OMP          SHARED(n_states, n_states_diag, N, shift2, U, h_val, W, H_diag, residual_norm, to_print)
      !$OMP DO 
      do k = 1, n_states_diag
        do i = 1, N
          U(i,shift2+k) = (h_val(k) * U(i,shift2+k) - W(i,shift2+k)) / max(H_diag(i) - h_val(k), 1.d-2)
        enddo
        if(k <= n_states) then
          residual_norm(k) = u_dot_u(U(1,shift2+k), N)
          to_print(1,k) = h_val(k)
          to_print(2,k) = residual_norm(k)
        endif
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !print*, " to_print", to_print

      if((itertot > 1) .and. (iter == 1)) then
        continue
      else
        write(*,'(1X, I3, 1X, 100(1X, F16.10, 1X, F12.6))') iter-1, to_print(1:2,1:n_states)
        !write(*, '(1X, I3, 1X, 100(1X, F16.10, 1X, F16.10, 1X, F16.10))') iter-1, to_print(1:2,1:n_states)
      endif

      !print*, 'iter = ', iter
      if(iter > 1) then
        converged = dabs(maxval(residual_norm(1:n_states))) < 1d-15
      endif

      do k = 1, n_states
        if(residual_norm(k) > 1.d8) then
          print *, 'Davidson failed'
          stop -1
        endif
      enddo

      if(converged) exit

    enddo ! loop over iter


    ! Re-contract U and update W
    ! --------------------------------

    call dgemm('N', 'N', N, n_states_diag, shift2, 1.d0,       &
               W(1,1), size(W, 1), h_vec(1,1), size(h_vec, 1), &
               0.d0, R, size(R, 1))
    do k = 1, n_states_diag
      do i = 1, N
        W(i,k) = R(i,k)
      enddo
    enddo

    call dgemm('N', 'N', N, n_states_diag, shift2, 1.d0,       &
               U(1,1), size(U, 1), h_vec(1,1), size(h_vec, 1), &
               0.d0, R(1,1), size(R, 1))

    do k = 1, n_states_diag
      do i = 1, N
        U(i,k) = R(i,k)
      enddo
    enddo

    call ortho_qr(U(1,1), size(U, 1), N, n_states_diag)
    !call ortho_qr(U(1,1), size(U, 1), N, n_states_diag)

    do j = 1, n_states_diag
      k = 1
      do while((k < N) .and. (U(k,j) == 0.d0))
        k = k+1
      enddo
      if(U(k,j) * R(k,j) < 0.d0) then
        do i = 1, N
          W(i,j) = -W(i,j)
        enddo
      endif
    enddo

  enddo ! loop over while

  ! ---

  write_buffer = '====='
  do i = 1, n_states
    write_buffer = trim(write_buffer)//' ================  ==========='
  enddo
  write(6,'(A)') trim(write_buffer)
  write(6,'(A)') ''


  print*, " Davidson eigenvalues"
  do k = 1, n_states
    Om(k) = h_val(k)
    print*, k, Om(k)
  enddo

  deallocate(H_diag)
  deallocate(U)
  deallocate(W)
  deallocate(h)
  deallocate(h_vec)
  deallocate(h_val)
  deallocate(overlap)
  deallocate(residual_norm)

  return
end

! ---

subroutine ppLR_HR_calc(ispin, nOrb, nC, nO, nR, nOO, nVV, lambda, e, eF, n_states_diag, ERI, U, W)

  implicit none

  integer,          intent(in)  :: ispin
  integer,          intent(in)  :: n_states_diag
  integer,          intent(in)  :: nOO, nVV, nOrb, nC, nO, nR
  double precision, intent(in)  :: lambda, eF
  double precision, intent(in)  :: e(nOrb)
  double precision, intent(in)  :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision, intent(in)  :: U(nOO+nVV,n_states_diag)
  double precision, intent(out) :: W(nOO+nVV,n_states_diag)

  integer                       :: i, j, ij, k, l, kl
  integer                       :: a, b, c, d, ab, cd
  integer                       :: state
  double precision              :: mat_tmp
  double precision              :: diff_loc, diff_tot
  double precision, allocatable :: M_ref(:,:), W_ref(:,:)
  double precision, allocatable :: Cpp_ref(:,:), Dpp_ref(:,:), Bpp_ref(:,:)

  double precision, external    :: Kronecker_delta

  if(ispin .eq. 1) then

    ab = 0
    do a = nO+1, nOrb-nR
      do b = a, nOrb-nR
        ab = ab + 1
  
        do state = 1, n_states_diag
  
          W(ab,state) = 0.d0
  
          cd = 0
          do c = nO+1, nOrb-nR
            do d = c, nOrb-nR
              cd = cd + 1
    
              mat_tmp = (e(a) + e(b) - eF) * Kronecker_delta(a, c) * Kronecker_delta(b, d) &
                      + lambda * (ERI(a,b,c,d) + ERI(a,b,d,c)) / dsqrt( (1.d0 + Kronecker_delta(a, b)) &
                                                                      * (1.d0 + Kronecker_delta(c, d)))

              W(ab,state) = W(ab,state) + mat_tmp * U(cd,state)
            enddo
          enddo
  
          ij = nVV
          do i = nC+1, nO
            do j = i, nO
              ij = ij + 1
  
              mat_tmp = lambda * (ERI(a,b,i,j) + ERI(a,b,j,i)) / dsqrt( (1.d0 + Kronecker_delta(a, b)) &
                                                                      * (1.d0 + Kronecker_delta(i, j)))
  
              W(ab,state) = W(ab,state) - mat_tmp * U(ij,state)
            enddo
          enddo
  
        enddo ! state
      enddo ! b
    enddo ! a

    ! ---
  
    ij = nVV
    do i = nC+1, nO
      do j = i, nO
        ij = ij + 1
  
        do state = 1, n_states_diag
  
          W(ij,state) = 0.d0
  
          cd = 0
          do c = nO+1, nOrb-nR
            do d = c, nOrb-nR
              cd = cd + 1
    
              mat_tmp = lambda * (ERI(c,d,i,j) + ERI(c,d,j,i)) / dsqrt( (1.d0 + Kronecker_delta(c, d)) &
                                                                      * (1.d0 + Kronecker_delta(i, j)))

              W(ij,state) = W(ij,state) + mat_tmp * U(cd,state)
            enddo
          enddo
  
          kl = nVV
          do k = nC+1, nO
            do l = k, nO
              kl = kl + 1
    
              mat_tmp = - (e(i) + e(j) - eF) * Kronecker_delta(i, k) * Kronecker_delta(j, l) &
                      + lambda * (ERI(i,j,k,l) + ERI(i,j,l,k)) / dsqrt( (1.d0 + Kronecker_delta(i, j)) &
                                                                      * (1.d0 + Kronecker_delta(k, l)))

              W(ij,state) = W(ij,state) - mat_tmp * U(kl,state)
            enddo
          enddo
  
        enddo ! state
      enddo ! j
    enddo ! i

  elseif((ispin .eq. 2) .or. (ispin .eq. 4)) then

    ab = 0
    do a = nO+1, nOrb-nR
      do b = a+1, nOrb-nR
        ab = ab + 1
  
        do state = 1, n_states_diag
  
          W(ab,state) = 0.d0

          cd = 0
          do c = nO+1, nOrb-nR
            do d = c+1, nOrb-nR
              cd = cd + 1
    
              mat_tmp = (e(a) + e(b) - eF) * Kronecker_delta(a, c) * Kronecker_delta(b, d) &
                      + lambda * (ERI(a,b,c,d) - ERI(a,b,d,c))

              W(ab,state) = W(ab,state) + mat_tmp * U(cd,state)
            enddo
          enddo

          ij = nVV
          do i = nC+1, nO
            do j = i+1, nO
              ij = ij + 1
  
              mat_tmp = lambda * (ERI(a,b,i,j) - ERI(a,b,j,i))
  
              W(ab,state) = W(ab,state) - mat_tmp * U(ij,state)
            enddo
          enddo
  
        enddo ! state
      enddo ! b
    enddo ! a

    ! ---
  
    ij = nVV
    do i = nC+1, nO
      do j = i+1, nO
        ij = ij + 1
  
        do state = 1, n_states_diag
  
          W(ij,state) = 0.d0
  
          cd = 0
          do c = nO+1, nOrb-nR
            do d = c+1, nOrb-nR
              cd = cd + 1
    
              mat_tmp = lambda * (ERI(c,d,i,j) - ERI(c,d,j,i))

              W(ij,state) = W(ij,state) + mat_tmp * U(cd,state)
            enddo
          enddo

          kl = nVV
          do k = nC+1, nO
            do l = k+1, nO
              kl = kl + 1
    
              mat_tmp = - (e(i) + e(j) - eF) * Kronecker_delta(i, k) * Kronecker_delta(j, l) &
                      + lambda * (ERI(i,j,k,l) - ERI(i,j,l,k))

              W(ij,state) = W(ij,state) - mat_tmp * U(kl,state)
            enddo
          enddo
  
        enddo ! state
      enddo ! j
    enddo ! i

  elseif(ispin .eq. 3) then

    ab = 0
    do a = nO+1, nOrb-nR
      do b = nO+1, nOrb-nR
        ab = ab + 1
  
        do state = 1, n_states_diag
  
          W(ab,state) = 0.d0
  
          cd = 0
          do c = nO+1, nOrb-nR
            do d = nO+1, nOrb-nR
              cd = cd + 1

              mat_tmp = (e(a) + e(b) - eF) * Kronecker_delta(a, c) * Kronecker_delta(b, d) &
                      + lambda * ERI(a,b,c,d)

              W(ab,state) = W(ab,state) + mat_tmp * U(cd,state)
            enddo
          enddo

          ij = nVV
          do i = nC+1, nO
            do j = nC+1, nO
              ij = ij + 1
  
              mat_tmp = lambda * ERI(a,b,i,j)
  
              W(ab,state) = W(ab,state) - mat_tmp * U(ij,state)
            enddo
          enddo
  
        enddo ! state
      enddo ! b
    enddo ! a
  
    ! ---
  
    ij = nVV
    do i = nC+1, nO
      do j = nC+1, nO
        ij = ij + 1
  
        do state = 1, n_states_diag
  
          W(ij,state) = 0.d0
  
          cd = 0
          do c = nO+1, nOrb-nR
            do d = nO+1, nOrb-nR
              cd = cd + 1
    
              mat_tmp = lambda * ERI(c,d,i,j)

              W(ij,state) = W(ij,state) + mat_tmp * U(cd,state)
            enddo
          enddo
  
          kl = nVV
          do k = nC+1, nO
            do l = nC+1, nO
              kl = kl + 1

              mat_tmp = - (e(i) + e(j) - eF) * Kronecker_delta(i, k) * Kronecker_delta(j, l) &
                      + lambda * ERI(i,j,k,l)

              W(ij,state) = W(ij,state) - mat_tmp * U(kl,state)
            enddo
          enddo
  
        enddo ! state
      enddo ! j
    enddo ! i

  else

    print*, ' ispin is not supported'
    print*, ' ispin = ', ispin
    stop

  endif


!  print*, ' debug ppLR_HR_calc:'
!  print*, ispin, nOO, nVV
!  allocate(M_ref(nOO+nVV,nOO+nVV))
!  allocate(Bpp_ref(nVV,nOO), Cpp_ref(nVV,nVV), Dpp_ref(nOO,nOO))
!  allocate(W_ref(nOO+nVV,n_states_diag))
!
!  call ppLR_C(ispin, nOrb, nC, nO, nOrb-nO, nR, nVV, 1d0, e, ERI, Cpp_ref)
!  call ppLR_D(ispin, nOrb, nC, nO, nOrb-nO, nR, nOO, 1d0, e, ERI, Dpp_ref)
!  call ppLR_B(ispin, nOrb, nC, nO, nOrb-nO, nR, nOO, nVV, 1d0, ERI, Bpp_ref)
!  M_ref = 0.d0
!  M_ref(    1:nVV    ,    1:nVV)     = + Cpp_ref(1:nVV,1:nVV)
!  M_ref(nVV+1:nVV+nOO,nVV+1:nVV+nOO) = - Dpp_ref(1:nOO,1:nOO)
!  M_ref(    1:nVV    ,nVV+1:nOO+nVV) = - Bpp_ref(1:nVV,1:nOO)
!  M_ref(nVV+1:nOO+nVV,    1:nVV)     = + transpose(Bpp_ref(1:nVV,1:nOO))
!
!  call dgemm('N', 'N', nOO+nVV, n_states_diag, nOO+nVV, 1.d0, &
!             M_ref(1,1), size(M_ref, 1), U(1,1), size(U, 1),  &
!             0.d0, W_ref(1,1), size(W_ref, 1))
!
!  diff_tot = 0.d0
!  do state = 1, n_states_diag
!    do ab = 1, nOO
!      diff_loc = dabs(W(ab,state) - W_ref(ab,state))
!      if(diff_loc .gt. 1d-12) then
!        print*, ' important diff on:', ab, state
!        print*, W(ab,state), W_ref(ab,state)
!        stop
!      endif
!      diff_tot = diff_tot + diff_loc
!    enddo
!    do ij = nVV+1, nVV+nOO
!      diff_loc = dabs(W(ij,state) - W_ref(ij,state))
!      if(diff_loc .gt. 1d-12) then
!        print*, ' important diff on:', ij, state
!        print*, W(ij,state), W_ref(ij,state)
!        stop
!      endif
!      diff_tot = diff_tot + diff_loc
!    enddo
!  enddo
!  print*, 'diff_tot = ', diff_tot
!
!  deallocate(M_ref)
!  deallocate(Bpp_ref, Cpp_ref, Dpp_ref)
!  deallocate(W_ref)

  return
end

! ---

subroutine ppLR_H_diag(ispin, nOrb, nC, nO, nR, nOO, nVV, lambda, e, eF, ERI, H_diag)

  implicit none

  integer,          intent(in)  :: ispin
  integer,          intent(in)  :: nOO, nVV, nOrb, nC, nO, nR
  double precision, intent(in)  :: lambda, eF
  double precision, intent(in)  :: e(nOrb)
  double precision, intent(in)  :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision, intent(out) :: H_diag(nOO+nVV)

  integer                       :: i, j, ij, k, l, kl
  integer                       :: a, b, c, d, ab, cd
  double precision              :: diff_loc, diff_tot
  double precision, allocatable :: M_ref(:,:)
  double precision, allocatable :: Cpp_ref(:,:), Dpp_ref(:,:), Bpp_ref(:,:)

  double precision, external    :: Kronecker_delta


  if(ispin .eq. 1) then

    ab = 0
    do a = nO+1, nOrb-nR
      do b = a, nOrb-nR
        ab = ab + 1
        cd = 0
        do c = nO+1, nOrb-nR
          do d = c, nOrb-nR
            cd = cd + 1
            if(a .ne. c) cycle
            if(b .ne. d) cycle
            H_diag(ab) = e(a) + e(b) - eF &
                       + lambda * (ERI(a,b,c,d) + ERI(a,b,d,c)) / dsqrt( (1.d0 + Kronecker_delta(a, b)) &
                                                                       * (1.d0 + Kronecker_delta(c, d)))
          enddo
        enddo
      enddo ! b
    enddo ! a
  
    ij = nVV
    do i = nC+1, nO
      do j = i, nO
        ij = ij + 1
        kl = 0
        do k = nC+1, nO
          do l = k, nO
            kl = kl + 1
            if(i .ne. k) cycle
            if(j .ne. l) cycle
            H_diag(ij) = e(i) + e(j) - eF &
                       - lambda * (ERI(i,j,k,l) + ERI(i,j,l,k)) / dsqrt( (1.d0 + Kronecker_delta(i, j)) &
                                                                       * (1.d0 + Kronecker_delta(k, l)))
          enddo
        enddo
      enddo ! j
    enddo ! i

  elseif((ispin .eq. 2) .or. (ispin .eq. 4)) then

    ab = 0
    do a = nO+1, nOrb-nR
      do b = a+1, nOrb-nR
        ab = ab + 1
        cd = 0
        do c = nO+1, nOrb-nR
          do d = c+1, nOrb-nR
            cd = cd + 1
            if(a .ne. c) cycle
            if(b .ne. d) cycle
            H_diag(ab) = e(a) + e(b) - eF + lambda * (ERI(a,b,c,d) - ERI(a,b,d,c))
          enddo
        enddo
      enddo ! b
    enddo ! a
  
    ij = nVV
    do i = nC+1, nO
      do j = i+1, nO
        ij = ij + 1
        kl = 0
        do k = nC+1, nO
          do l = k+1, nO
            kl = kl + 1
            if(i .ne. k) cycle
            if(j .ne. l) cycle
            H_diag(ij) = e(i) + e(j) - eF - lambda * (ERI(i,j,k,l) - ERI(i,j,l,k))
          enddo
        enddo
      enddo ! j
    enddo ! i

  elseif(ispin .eq. 3) then

    ab = 0
    do a = nO+1, nOrb-nR
      do b = nO+1, nOrb-nR
        ab = ab + 1
        cd = 0
        do c = nO+1, nOrb-nR
          do d = nO+1, nOrb-nR
            cd = cd + 1
            if(a .ne. c) cycle
            if(b .ne. d) cycle
            H_diag(ab) = (e(a) + e(b) - eF) + lambda * ERI(a,b,c,d)
          enddo
        enddo
      enddo ! b
    enddo ! a

    ij = nVV
    do i = nC+1, nO
      do j = nC+1, nO
        ij = ij + 1 
        kl = 0
        do k = nC+1, nO
          do l = nC+1, nO
            kl = kl + 1
            if(i .ne. k) cycle
            if(j .ne. l) cycle
            H_diag(ij) = (e(i) + e(j) - eF) - lambda * ERI(i,j,k,l)
          enddo
        enddo
      enddo ! j
    enddo ! i

  else

    print*, ' ispin is not supported'
    print*, ' ispin = ', ispin
    stop

  endif


!  print*, ' debug ppLR_H_diag:'
!  print*, ispin, nOO, nVV
!  allocate(M_ref(nOO+nVV,nOO+nVV))
!  allocate(Bpp_ref(nVV,nOO), Cpp_ref(nVV,nVV), Dpp_ref(nOO,nOO))
!
!  call ppLR_C(ispin, nOrb, nC, nO, nOrb-nO, nR, nVV, 1d0, e, ERI, Cpp_ref)
!  call ppLR_D(ispin, nOrb, nC, nO, nOrb-nO, nR, nOO, 1d0, e, ERI, Dpp_ref)
!  call ppLR_B(ispin, nOrb, nC, nO, nOrb-nO, nR, nOO, nVV, 1d0, ERI, Bpp_ref)
!  M_ref = 0.d0
!  M_ref(    1:nVV    ,    1:nVV)     = + Cpp_ref(1:nVV,1:nVV)
!  M_ref(nVV+1:nVV+nOO,nVV+1:nVV+nOO) = - Dpp_ref(1:nOO,1:nOO)
!  M_ref(    1:nVV    ,nVV+1:nOO+nVV) = - Bpp_ref(1:nVV,1:nOO)
!  M_ref(nVV+1:nOO+nVV,    1:nVV)     = + transpose(Bpp_ref(1:nVV,1:nOO))
!
!  diff_tot = 0.d0
!  do ab = 1, nOO
!    diff_loc = dabs(H_diag(ab) - M_ref(ab,ab))
!    if(diff_loc .gt. 1d-12) then
!      print*, ' important diff on:', ab
!      print*, H_diag(ab), M_ref(ab,ab)
!      stop
!    endif
!    diff_tot = diff_tot + diff_loc
!  enddo
!  do ij = nVV+1, nVV+nOO
!    diff_loc = dabs(H_diag(ij) - M_ref(ij,ij))
!    if(diff_loc .gt. 1d-12) then
!      print*, ' important diff on:', ij
!      print*, H_diag(ij), M_ref(ij,ij)
!      stop
!    endif
!    diff_tot = diff_tot + diff_loc
!  enddo
!  print*, 'diff_tot = ', diff_tot
!
!  deallocate(M_ref)
!  deallocate(Bpp_ref, Cpp_ref, Dpp_ref)

  return
end

! ---

subroutine diag_nonsym_right(n, A, A_ldim, V, V_ldim, energy, E_ldim)

  implicit none

  integer,          intent(in)  :: n, A_ldim, V_ldim, E_ldim
  double precision, intent(in)  :: A(A_ldim,n)
  double precision, intent(out) :: energy(E_ldim), V(V_ldim,n)

  character*1                   :: JOBVL, JOBVR, BALANC, SENSE
  integer                       :: i, j
  integer                       :: ILO, IHI, lda, ldvl, ldvr, LWORK, INFO
  double precision              :: ABNRM
  integer,          allocatable :: iorder(:), IWORK(:)
  double precision, allocatable :: WORK(:), SCALE_array(:), RCONDE(:), RCONDV(:)
  double precision, allocatable :: Atmp(:,:), WR(:), WI(:), VL(:,:), VR(:,:), Vtmp(:)
  double precision, allocatable :: energy_loc(:), V_loc(:,:)

  !print*, " n = ", n
  allocate(Atmp(n,n))
  allocate(WI(n))
  allocate(WR(n))
  allocate(VR(n,n))
  allocate(VL(1,1))

  do i = 1, n
    do j = 1, n
      Atmp(j,i) = A(j,i)
    enddo
  enddo

  JOBVL  = "N" ! computes the left  eigenvectors
  JOBVR  = "V" ! computes the right eigenvectors
  BALANC = "B" ! Diagonal scaling and Permutation for optimization
  SENSE  = "V" ! Determines which reciprocal condition numbers are computed
  lda  = n
  ldvr = n
  ldvl = 1

  allocate(WORK(1), SCALE_array(n), RCONDE(n), RCONDV(n), IWORK(2*n-2))

  LWORK = -1 ! to ask for the optimal size of WORK
  call dgeevx(BALANC, JOBVL, JOBVR, SENSE,                  & ! CHARACTERS
              n, Atmp, lda,                                 & ! MATRIX TO DIAGONALIZE
              WR, WI,                                       & ! REAL AND IMAGINARY PART OF EIGENVALUES
              VL, ldvl, VR, ldvr,                           & ! LEFT AND RIGHT EIGENVECTORS
              ILO, IHI, SCALE_array, ABNRM, RCONDE, RCONDV, & ! OUTPUTS OF OPTIMIZATION
              WORK, LWORK, IWORK, INFO)

  if(INFO .ne. 0) then
    print*, 'dgeevx failed !!', INFO
    stop
  endif

  LWORK = max(int(work(1)), 1) ! this is the optimal size of WORK
  deallocate(WORK)
  allocate(WORK(LWORK))
  call dgeevx(BALANC, JOBVL, JOBVR, SENSE,                  &
              n, Atmp, lda,                                 &
              WR, WI,                                       &
              VL, ldvl, VR, ldvr,                           &
              ILO, IHI, SCALE_array, ABNRM, RCONDE, RCONDV, &
              WORK, LWORK, IWORK, INFO)

  if(INFO .ne. 0) then
    print*, 'dgeevx failed !!', INFO
    stop
  endif

  deallocate(WORK, SCALE_array, RCONDE, RCONDV, IWORK)
  deallocate(VL, Atmp)


  allocate(energy_loc(n), V_loc(n,n))
  energy_loc = 0.d0
  V_loc = 0.d0

  i = 1
  do while(i .le. n)

    !print*, i, WR(i), WI(i)

    if(dabs(WI(i)) .gt. 1d-7) then

      print*, ' Found an imaginary component to eigenvalue'
      print*, ' Re(i) + Im(i)', i, WR(i), WI(i)

      energy_loc(i) = WR(i)
      do j = 1, n
        V_loc(j,i) = WR(i) * VR(j,i) - WI(i) * VR(j,i+1)
      enddo
      energy_loc(i+1) = WI(i)
      do j = 1, n
        V_loc(j,i+1) = WR(i) * VR(j,i+1) + WI(i) * VR(j,i)
      enddo
      i = i + 2

    else

      energy_loc(i) = WR(i)
      do j = 1, n
        V_loc(j,i) = VR(j,i)
      enddo
      i = i + 1

    endif

  enddo

  deallocate(WR, WI, VR)


  ! ordering
!  do j = 1, n
!    write(444, '(100(1X, F16.10))') (V_loc(j,i), i=1,5)
!  enddo
  allocate(iorder(n))
  do i = 1, n
    iorder(i) = i
  enddo
  call quick_sort(energy_loc, iorder, n)
  do i = 1, n
    energy(i) = energy_loc(i)
    do j = 1, n
      V(j,i) = V_loc(j,iorder(i))
    enddo
  enddo
  deallocate(iorder)
!  do j = 1, n
!    write(445, '(100(1X, F16.10))') (V_loc(j,i), i=1,5)
!  enddo
  deallocate(V_loc, energy_loc)

end

! ---

