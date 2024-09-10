

! ---

subroutine ppLR_GW_HR_calc(ispin, nOrb, nC, nO, nR, nOO, nVV, nS, lambda, e, eF, n_states_diag, &
                           ERI, eta, rho, Om, U, W)


  implicit none

  integer,          intent(in)  :: ispin
  integer,          intent(in)  :: n_states_diag
  integer,          intent(in)  :: nOO, nVV, nOrb, nC, nO, nR, nS
  double precision, intent(in)  :: lambda, eF, eta
  double precision, intent(in)  :: e(nOrb)
  double precision, intent(in)  :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision, intent(in)  :: rho(nS,nOrb,nOrb), Om(nS)
  double precision, intent(in)  :: U(nOO+nVV,n_states_diag)
  double precision, intent(out) :: W(nOO+nVV,n_states_diag)

  integer                       :: i, j, ij
  integer                       :: ab
  integer                       :: m
  integer                       :: state
  double precision              :: eta2
  double precision              :: t1, t2
  double precision              :: diff_tot, diff_loc
  double precision, allocatable :: M_ref(:,:)
  double precision, allocatable :: Bpp_ref(:,:), Cpp_ref(:,:), Dpp_ref(:,:)
  double precision, allocatable :: KB_sta(:,:), KC_sta(:,:), KD_sta(:,:)
  double precision, allocatable :: W_ref(:,:)
  double precision, allocatable :: rho_t(:,:,:)

!  call wall_time(t1)

  if((nOO+nVV) .le. 20000) then

    call ppLR_GW_HR_calc_oneshot(ispin, nOrb, nC, nO, nR, nOO, nVV, nS, lambda, e(1), eF, n_states_diag, &
                                 ERI(1,1,1,1), eta, rho(1,1,1), Om(1), U(1,1), W(1,1))

  else

    call ppLR_GW_HR_calc_batches(ispin, nOrb, nC, nO, nR, nOO, nVV, nS, lambda, e(1), eF, n_states_diag, &
                                 ERI(1,1,1,1), eta, rho(1,1,1), Om(1), U(1,1), W(1,1))

  endif


!  print*, ' debug ppLR_GW_H_diag:'
!  allocate(M_ref(nOO+nVV,nOO+nVV))
!  allocate(Bpp_ref(nVV,nOO), Cpp_ref(nVV,nVV), Dpp_ref(nOO,nOO))
!  allocate(KB_sta(nVV,nOO), KC_sta(nVV,nVV), KD_sta(nOO,nOO))
!  allocate(rho_t(nOrb,nOrb,nS))
!  allocate(W_ref(nOO+nVV,n_states_diag))
!
!  call ppLR_C(ispin, nOrb, nC, nO, nOrb-nO, nR, nVV,      1d0, e, ERI, Cpp_ref)
!  call ppLR_D(ispin, nOrb, nC, nO, nOrb-nO, nR, nOO,      1d0, e, ERI, Dpp_ref)
!  call ppLR_B(ispin, nOrb, nC, nO, nOrb-nO, nR, nOO, nVV, 1d0,    ERI, Bpp_ref)
!
!  do j = 1, nOrb
!    do i = 1, nOrb
!      do m = 1, nS
!        rho_t(i,j,m) = rho(m,i,j) 
!      enddo
!    enddo
!  enddo
!
!  call RGW_ppBSE_static_kernel_C(ispin, eta, nOrb, nC, nO, nOrb-nO, nR, nS, nVV,      1.d0, ERI, Om, rho_t, KC_sta)
!  call RGW_ppBSE_static_kernel_D(ispin, eta, nOrb, nC, nO, nOrb-nO, nR, nS, nOO,      1.d0, ERI, Om, rho_t, KD_sta)
!  call RGW_ppBSE_static_kernel_B(ispin, eta, nOrb, nC, nO, nOrb-nO, nR, nS, nOO, nVV, 1.d0, ERI, Om, rho_t, KB_sta)
!
!  Cpp_ref = Cpp_ref + KC_sta
!  Dpp_ref = Dpp_ref + KD_sta
!  Bpp_ref = Bpp_ref + KB_sta
!
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
!      if(diff_loc .gt. 1d-10) then
!        print*, ' important diff on:', ab, state
!        print*, W(ab,state), W_ref(ab,state)
!        stop
!      endif
!      diff_tot = diff_tot + diff_loc
!    enddo
!    do ij = nVV+1, nVV+nOO
!      diff_loc = dabs(W(ij,state) - W_ref(ij,state))
!      if(diff_loc .gt. 1d-10) then
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
!  deallocate(KB_sta, KC_sta, KD_sta)
!  deallocate(W_ref)
!  deallocate(rho_t)


!  call wall_time(t2)
!  write(*,'(A50, F12.4)') 'total wall time for ppLR_GW_HR_calc (sec): ', t2-t1

  return
end

! ---

subroutine ppLR_GW_H_diag(ispin, nOrb, nC, nO, nR, nOO, nVV, nS, lambda, e, eF, ERI, eta, rho, Om, H_diag)

  implicit none

  integer,          intent(in)  :: ispin
  integer,          intent(in)  :: nOO, nVV, nOrb, nC, nO, nR, nS
  double precision, intent(in)  :: lambda, eF, eta
  double precision, intent(in)  :: e(nOrb)
  double precision, intent(in)  :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision, intent(in)  :: rho(nS,nOrb,nOrb), Om(nS)
  double precision, intent(out) :: H_diag(nOO+nVV)

  integer                       :: i, j, ij, k, l, kl
  integer                       :: a, b, c, d, ab, cd
  integer                       :: m
  double precision              :: chi, eps
  double precision              :: t1, t2

  double precision, external    :: Kronecker_delta

  call wall_time(t1)

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
            chi = 0.d0
            do m = 1, nS
              eps = Om(m)**2 + eta**2
              chi = chi - rho(m,a,c) * rho(m,b,d) * Om(m) / eps &
                        - rho(m,a,d) * rho(m,b,c) * Om(m) / eps
            end do
            H_diag(ab) = H_diag(ab) + 4.d0 * lambda * chi / dsqrt( (1.d0 + Kronecker_delta(a, b)) &
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
            chi = 0.d0
            do m = 1, nS
              eps = Om(m)**2 + eta**2
              chi = chi - rho(m,i,k) * rho(m,j,l) * Om(m) / eps &
                        - rho(m,i,l) * rho(m,j,k) * Om(m) / eps
            enddo
            H_diag(ij) = H_diag(ij) - 4.d0 * lambda * chi / dsqrt( (1.d0 + Kronecker_delta(i, j)) &
                                                                 * (1.d0 + Kronecker_delta(k, l)))
          enddo
        enddo
      enddo ! j
    enddo ! i

  elseif(ispin .eq. 2) then

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
            chi = 0.d0
            do m = 1, nS
              eps = Om(m)**2 + eta**2
              chi = chi - rho(m,a,c) * rho(m,b,d) * Om(m) / eps &
                        + rho(m,a,d) * rho(m,b,c) * Om(m) / eps
            enddo
            H_diag(ab) = H_diag(ab) + 4.d0 * lambda * chi
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
            chi = 0.d0
            do m = 1, nS
              eps = Om(m)**2 + eta**2
              chi = chi - rho(m,i,k) * rho(m,j,l) * Om(m) / eps &
                        + rho(m,i,l) * rho(m,j,k) * Om(m) / eps
            end do
            H_diag(ij) = H_diag(ij) - 4.d0 * lambda * chi
          enddo
        enddo
      enddo ! j
    enddo ! i

  else

    print*, ' Error in ppLR_GW_H_diag'
    print*, ' ispin is not supported'
    print*, ' ispin = ', ispin
    stop

  endif

  call wall_time(t2)
  write(*,'(A50, F12.4)') 'total wall time for ppLR_GW_H_diag (sec): ', t2-t1

  return
end

! ---

subroutine ppLR_GW_HR_calc_oneshot(ispin, nOrb, nC, nO, nR, nOO, nVV, nS, lambda, e, eF, n_states_diag, &
                                   ERI, eta, rho, Om, U, W)

  implicit none

  integer,          intent(in)  :: ispin
  integer,          intent(in)  :: n_states_diag
  integer,          intent(in)  :: nOO, nVV, nOrb, nC, nO, nR, nS
  double precision, intent(in)  :: lambda, eF, eta
  double precision, intent(in)  :: e(nOrb)
  double precision, intent(in)  :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision, intent(in)  :: rho(nS,nOrb,nOrb), Om(nS)
  double precision, intent(in)  :: U(nOO+nVV,n_states_diag)
  double precision, intent(out) :: W(nOO+nVV,n_states_diag)

  integer                       :: i, j, ij, k, l, kl
  integer                       :: a, b, c, d, ab, cd
  integer                       :: a0, aa, i0, ii
  integer                       :: m
  integer                       :: state
  double precision              :: mat_tmp, chi, eps
  double precision              :: eta2
  double precision              :: tmp_e, tmp_ab, tmp_ij
  double precision, allocatable :: Om_tmp(:), H_mat(:,:)


  if(ispin .eq. 1) then

    allocate(Om_tmp(nS))

    a0 = nOrb - nR - nO

    eta2 = eta * eta
    do m = 1, nS
      Om_tmp(m) = Om(m) / (Om(m) * Om(m) + eta2)
    enddo

    allocate(H_mat(nVV,nOO+nVV))

    !$OMP PARALLEL                                            &
    !$OMP DEFAULT(NONE)                                       &
    !$OMP PRIVATE(a, b, aa, ab, c, d, cd, i, j, ij, m, state, &
    !$OMP         tmp_e, tmp_ab, chi, mat_tmp)                &
    !$OMP SHARED(nC, nO, nOrb, nR, nS, n_states_diag, nVV,    &
    !$OMP        nOO, a0, eF, lambda, e, Om_tmp, rho, ERI, H_mat)
    !$OMP DO SCHEDULE(GUIDED)
    do a = nO+1, nOrb-nR
      aa = a0 * (a - nO - 1) - (a - nO - 1) * (a - nO) / 2 - nO
      do b = a, nOrb-nR
        ab = aa + b

        tmp_e = e(a) + e(b) - eF

        tmp_ab = lambda
        if(a .eq. b) then
          tmp_ab = 0.7071067811865475d0 * tmp_ab
        endif

        cd = 0
        do c = nO+1, nOrb-nR
          do d = c, nOrb-nR
            cd = cd + 1

            chi = 0.d0
            do m = 1, nS
              chi = chi - Om_tmp(m) * (rho(m,a,c) * rho(m,b,d) + rho(m,a,d) * rho(m,b,c))
            enddo

            mat_tmp = tmp_ab
            if(c .eq. d) then
              mat_tmp = 0.7071067811865475d0 * mat_tmp
            endif

            mat_tmp = mat_tmp * (4.d0 * chi + ERI(d,c,b,a) + ERI(d,c,a,b))
            if((a .eq. c) .and. (b .eq. d)) then
              mat_tmp = mat_tmp + tmp_e
            endif

            H_mat(ab,cd) = mat_tmp
          enddo ! d
        enddo ! c

        ij = nVV
        do i = nC+1, nO
          do j = i, nO
            ij = ij + 1

            chi = 0.d0
            do m = 1, nS
              chi = chi - Om_tmp(m) * (rho(m,i,a) * rho(m,j,b) + rho(m,i,b) * rho(m,a,j))
            enddo

            mat_tmp = tmp_ab
            if(i .eq. j) then
              mat_tmp = 0.7071067811865475d0 * mat_tmp
            endif

            mat_tmp = mat_tmp * (4.d0 * chi + ERI(j,i,b,a) + ERI(j,i,a,b))

            H_mat(ab,ij) = -mat_tmp
          enddo ! j
        enddo ! i

      enddo ! b
    enddo ! a
    !$OMP END DO
    !$OMP END PARALLEL

    call dgemm("N", "N", nVV, n_states_diag, nOO+nVV, 1.d0,    &
               H_mat(1,1), size(H_mat, 1), U(1,1), size(U, 1), &
               0.d0, W(1,1), size(W, 1))

    deallocate(H_mat)

    ! ---

    allocate(H_mat(nOO,nOO+nVV))

    i0 = nO - nC

    !$OMP PARALLEL                                            &
    !$OMP DEFAULT(NONE)                                       &
    !$OMP PRIVATE(i, j, ii, ij, a, b, ab, k, l, kl, m, state, &
    !$OMP         tmp_e, tmp_ij, chi, mat_tmp)                &
    !$OMP SHARED(nC, nO, nOrb, nR, nS, nVV, i0,               &
    !$OMP        eF, lambda, e, Om_tmp, rho, ERI, U, H_mat)
    !$OMP DO SCHEDULE(GUIDED)
    do i = nC+1, nO
      ii = i0 * (i - nC - 1) - (i - nC - 1) * (i - nC) / 2 - nC
      do j = i, nO
        ij = ii + j

        tmp_e = e(i) + e(j) - eF

        tmp_ij = lambda
        if(i .eq. j) then
          tmp_ij = 0.7071067811865475d0 * tmp_ij
        endif
  
        ab = 0
        do a = nO+1, nOrb-nR
          do b = a, nOrb-nR
            ab = ab + 1

            chi = 0.d0
            do m = 1, nS
              chi = chi - Om_tmp(m) * (rho(m,i,a) * rho(m,j,b) + rho(m,i,b) * rho(m,a,j))
            enddo

            mat_tmp = tmp_ij
            if(a .eq. b) then
              mat_tmp = 0.7071067811865475d0 * mat_tmp
            endif
    
            mat_tmp = mat_tmp * (4.d0 * chi + ERI(b,a,j,i) + ERI(b,a,i,j))

            H_mat(ij,ab) = mat_tmp
          enddo ! b
        enddo ! a

        kl = nVV
        do k = nC+1, nO
          do l = k, nO
            kl = kl + 1
    
            chi = 0.d0
            do m = 1, nS
              chi = chi - Om_tmp(m) * (rho(m,i,k) * rho(m,j,l) + rho(m,i,l) * rho(m,j,k))
            enddo

            mat_tmp = tmp_ij
            if(k .eq. l) then
              mat_tmp = 0.7071067811865475d0 * mat_tmp
            endif

            mat_tmp = mat_tmp * (4.d0 * chi + ERI(l,k,j,i) + ERI(l,k,i,j))

            if((i .eq. k) .and. (j .eq. l)) then
              mat_tmp = mat_tmp - tmp_e
            endif

            H_mat(ij,kl) = -mat_tmp
          enddo ! l
        enddo ! k
      enddo ! j
    enddo ! i
    !$OMP END DO
    !$OMP END PARALLEL

    call dgemm("N", "N", nOO, n_states_diag, nOO+nVV, 1.d0,    &
               H_mat(1,1), size(H_mat, 1), U(1,1), size(U, 1), &
               0.d0, W(nVV+1,1), size(W, 1))
    deallocate(H_mat)

    deallocate(Om_tmp)

    ! ---

  elseif(ispin .eq. 2) then

    ! ---

    allocate(Om_tmp(nS))

    a0 = nOrb - nR - nO - 1

    eta2 = eta * eta
    do m = 1, nS
      Om_tmp(m) = Om(m) / (Om(m) * Om(m) + eta2)
    enddo

    allocate(H_mat(nVV,nOO+nVV))

    !$OMP PARALLEL                                            &
    !$OMP DEFAULT(NONE)                                       &
    !$OMP PRIVATE(a, b, aa, ab, c, d, cd, i, j, ij, m, state, &
    !$OMP         tmp_e, tmp_ab, chi, mat_tmp)                &
    !$OMP SHARED(nC, nO, nOrb, nR, nS, n_states_diag, nVV,    &
    !$OMP        nOO, a0, eF, lambda, e, Om_tmp, rho, ERI, H_mat)
    !$OMP DO SCHEDULE(GUIDED)
    do a = nO+1, nOrb-nR
      aa = a0 * (a - nO - 1) - (a - nO - 1) * (a - nO) / 2 - nO - 1
      do b = a+1, nOrb-nR
        ab = aa + b

        tmp_e = e(a) + e(b) - eF

        cd = 0
        do c = nO+1, nOrb-nR
          do d = c+1, nOrb-nR
            cd = cd + 1
    
            mat_tmp = lambda * (ERI(d,c,b,a) - ERI(d,c,a,b))

            chi = 0.d0
            do m = 1, nS
              chi = chi - Om_tmp(m) * (rho(m,a,c) * rho(m,b,d) - rho(m,a,d) * rho(m,b,c))
            enddo
            mat_tmp = mat_tmp + 4.d0 * lambda * chi

            if((a .eq. c) .and. (b .eq. d)) then
              mat_tmp = mat_tmp + tmp_e
            endif

            H_mat(ab,cd) = mat_tmp
          enddo ! d
        enddo ! c

        ij = nVV
        do i = nC+1, nO
          do j = i+1, nO
            ij = ij + 1

            mat_tmp = lambda * (ERI(j,i,b,a) - ERI(j,i,a,b))

            chi = 0.d0
            do m = 1, nS
              chi = chi - Om_tmp(m) * (rho(m,i,a) * rho(m,j,b) - rho(m,i,b) * rho(m,a,j))
            end do
            mat_tmp = mat_tmp + 4.d0 * lambda * chi

            H_mat(ab,ij) = -mat_tmp
          enddo ! j
        enddo ! i
      enddo ! b
    enddo ! a
    !$OMP END DO
    !$OMP END PARALLEL

    call dgemm("N", "N", nVV, n_states_diag, nOO+nVV, 1.d0,    &
               H_mat(1,1), size(H_mat, 1), U(1,1), size(U, 1), &
               0.d0, W(1,1), size(W, 1))

    deallocate(H_mat)

    ! ---

    i0 = nO - nC - 1
    allocate(H_mat(nOO,nOO+nVV))

    !$OMP PARALLEL                                            &
    !$OMP DEFAULT(NONE)                                       &
    !$OMP PRIVATE(i, j, ii, ij, a, b, ab, k, l, kl, m, state, &
    !$OMP         tmp_e, tmp_ij, chi, mat_tmp)                &
    !$OMP SHARED(nC, nO, nOrb, nR, nS, nVV, i0,               &
    !$OMP        eF, lambda, e, Om_tmp, rho, ERI, U, H_mat)
    !$OMP DO SCHEDULE(GUIDED)
    do i = nC+1, nO
      ii = i0 * (i - nC - 1) - (i - nC - 1) * (i - nC) / 2 - nC - 1
      do j = i+1, nO
        ij = ii + j

        tmp_e = e(i) + e(j) - eF
  
        ab = 0
        do a = nO+1, nOrb-nR
          do b = a+1, nOrb-nR
            ab = ab + 1
    
            mat_tmp = lambda * (ERI(b,a,j,i) - ERI(b,a,i,j))

            chi = 0.d0
            do m = 1, nS
              chi = chi - Om_tmp(m) * (rho(m,i,a) * rho(m,j,b) - rho(m,i,b) * rho(m,a,j))
            enddo
            mat_tmp = mat_tmp + 4.d0 * lambda * chi

            H_mat(ij,ab) = mat_tmp
          enddo ! b
        enddo ! a

        kl = nVV
        do k = nC+1, nO
          do l = k+1, nO
            kl = kl + 1
   
            mat_tmp = lambda * (ERI(l,k,j,i) - ERI(l,k,i,j))

            chi = 0.d0
            do m = 1, nS
              chi = chi - Om_tmp(m) * (rho(m,i,k) * rho(m,j,l) - rho(m,i,l) * rho(m,j,k))
            enddo
            mat_tmp = mat_tmp + 4.d0 * lambda * chi

            if((i .eq. k) .and. (j .eq. l)) then
              mat_tmp = mat_tmp - tmp_e
            endif

            H_mat(ij,kl) = -mat_tmp
          enddo ! l
        enddo ! k
      enddo ! j
    enddo ! i
    !$OMP END DO
    !$OMP END PARALLEL

    call dgemm("N", "N", nOO, n_states_diag, nOO+nVV, 1.d0,    &
               H_mat(1,1), size(H_mat, 1), U(1,1), size(U, 1), &
               0.d0, W(nVV+1,1), size(W, 1))
    deallocate(H_mat)

    deallocate(Om_tmp)


  else

    print*, ' Error in ppLR_GW_HR_calc_oneshot'
    print*, ' ispin is not supported'
    print*, ' ispin = ', ispin
    stop

  endif

  return
end

! ---

subroutine ppLR_GW_HR_calc_batches(ispin, nOrb, nC, nO, nR, nOO, nVV, nS, lambda, e, eF, n_states_diag, &
                                   ERI, eta, rho, Om, U, W)


  use omp_lib

  implicit none

  integer,          intent(in)  :: ispin
  integer,          intent(in)  :: n_states_diag
  integer,          intent(in)  :: nOO, nVV, nOrb, nC, nO, nR, nS
  double precision, intent(in)  :: lambda, eF, eta
  double precision, intent(in)  :: e(nOrb)
  double precision, intent(in)  :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision, intent(in)  :: rho(nS,nOrb,nOrb), Om(nS)
  double precision, intent(in)  :: U(nOO+nVV,n_states_diag)
  double precision, intent(out) :: W(nOO+nVV,n_states_diag)

  integer                       :: i, j, ij, k, l, kl
  integer                       :: a, b, c, d, ab, cd
  integer                       :: a0, aa, i0, ii, bb
  integer                       :: m
  integer                       :: state
  double precision              :: mat_tmp, chi, eps
  double precision              :: eta2
  double precision              :: tmp_e, tmp_ab, tmp_ij
  double precision, allocatable :: Om_tmp(:), H_mat(:,:)

  double precision, external    :: Kronecker_delta


  if(ispin .eq. 1) then

    call omp_set_max_active_levels(1)

    allocate(Om_tmp(nS))

    a0 = nOrb - nR - nO

    eta2 = eta * eta
    do m = 1, nS
      Om_tmp(m) = Om(m) / (Om(m) * Om(m) + eta2)
    enddo

    !$OMP PARALLEL                                            &
    !$OMP DEFAULT(NONE)                                       &
    !$OMP PRIVATE(a, b, aa, ab, c, d, cd, i, j, ij, m, state, &
    !$OMP         bb, tmp_e, tmp_ab, chi, mat_tmp, H_mat)     &
    !$OMP SHARED(nC, nO, nOrb, nR, nS, n_states_diag, nVV,    &
    !$OMP        nOO, a0, eF, lambda, e, Om_tmp, rho, ERI, U, W)

    allocate(H_mat(nOO+nVV,a0))

    !$OMP DO SCHEDULE(GUIDED)
    do a = nO+1, nOrb-nR
      aa = a0 * (a - nO - 1) - (a - nO - 1) * (a - nO) / 2 - nO
      do b = a, nOrb-nR
        ab = aa + b

        bb = b - a + 1

        tmp_e = e(a) + e(b) - eF

        tmp_ab = lambda
        if(a .eq. b) then
          tmp_ab = 0.7071067811865475d0 * tmp_ab
        endif

        cd = 0
        do c = nO+1, nOrb-nR
          do d = c, nOrb-nR
            cd = cd + 1

            chi = 0.d0
            do m = 1, nS
              chi = chi - Om_tmp(m) * (rho(m,a,c) * rho(m,b,d) + rho(m,a,d) * rho(m,b,c))
            enddo

            mat_tmp = tmp_ab
            if(c .eq. d) then
              mat_tmp = 0.7071067811865475d0 * mat_tmp
            endif

            mat_tmp = mat_tmp * (4.d0 * chi + ERI(d,c,b,a) + ERI(d,c,a,b))
            if((a .eq. c) .and. (b .eq. d)) then
              mat_tmp = mat_tmp + tmp_e
            endif

            H_mat(cd,bb) = mat_tmp
          enddo ! d
        enddo ! c

        ij = nVV
        do i = nC+1, nO
          do j = i, nO
            ij = ij + 1

            chi = 0.d0
            do m = 1, nS
              chi = chi - Om_tmp(m) * (rho(m,i,a) * rho(m,j,b) + rho(m,i,b) * rho(m,a,j))
            enddo

            mat_tmp = tmp_ab
            if(i .eq. j) then
              mat_tmp = 0.7071067811865475d0 * mat_tmp
            endif

            mat_tmp = mat_tmp * (4.d0 * chi + ERI(j,i,b,a) + ERI(j,i,a,b))

            H_mat(ij,bb) = -mat_tmp
          enddo ! j
        enddo ! i
      enddo ! b

      call dgemm("T", "N", nOrb-nR-a+1, n_states_diag, nOO+nVV, 1.d0, &
                 H_mat(1,1), size(H_mat, 1), U(1,1), size(U, 1),      &
                 0.d0, W(aa+a,1), size(W, 1))

    enddo ! a
    !$OMP END DO

    deallocate(H_mat)

    !$OMP END PARALLEL

    ! ---

    i0 = nO - nC
    allocate(H_mat(nOO,nOO+nVV))

    !$OMP PARALLEL                                            &
    !$OMP DEFAULT(NONE)                                       &
    !$OMP PRIVATE(i, j, ii, ij, a, b, ab, k, l, kl, m, state, &
    !$OMP         tmp_e, tmp_ij, chi, mat_tmp)                &
    !$OMP SHARED(nC, nO, nOrb, nR, nS, nVV, i0,               &
    !$OMP        eF, lambda, e, Om_tmp, rho, ERI, U, H_mat)
    !$OMP DO SCHEDULE(GUIDED)
    do i = nC+1, nO
      ii = i0 * (i - nC - 1) - (i - nC - 1) * (i - nC) / 2 - nC
      do j = i, nO
        ij = ii + j

        tmp_e = e(i) + e(j) - eF

        tmp_ij = lambda
        if(i .eq. j) then
          tmp_ij = 0.7071067811865475d0 * tmp_ij
        endif
  
        ab = 0
        do a = nO+1, nOrb-nR
          do b = a, nOrb-nR
            ab = ab + 1

            chi = 0.d0
            do m = 1, nS
              chi = chi - Om_tmp(m) * (rho(m,i,a) * rho(m,j,b) + rho(m,i,b) * rho(m,a,j))
            enddo

            mat_tmp = tmp_ij
            if(a .eq. b) then
              mat_tmp = 0.7071067811865475d0 * mat_tmp
            endif
    
            mat_tmp = mat_tmp * (4.d0 * chi + ERI(b,a,j,i) + ERI(b,a,i,j))

            H_mat(ij,ab) = mat_tmp
          enddo ! b
        enddo ! a

        kl = nVV
        do k = nC+1, nO
          do l = k, nO
            kl = kl + 1
    
            chi = 0.d0
            do m = 1, nS
              chi = chi - Om_tmp(m) * (rho(m,i,k) * rho(m,j,l) + rho(m,i,l) * rho(m,j,k))
            enddo

            mat_tmp = tmp_ij
            if(k .eq. l) then
              mat_tmp = 0.7071067811865475d0 * mat_tmp
            endif

            mat_tmp = mat_tmp * (4.d0 * chi + ERI(l,k,j,i) + ERI(l,k,i,j))

            if((i .eq. k) .and. (j .eq. l)) then
              mat_tmp = mat_tmp - tmp_e
            endif

            H_mat(ij,kl) = -mat_tmp
          enddo ! l
        enddo ! k
      enddo ! j
    enddo ! i
    !$OMP END DO
    !$OMP END PARALLEL

    call dgemm("N", "N", nOO, n_states_diag, nOO+nVV, 1.d0,    &
               H_mat(1,1), size(H_mat, 1), U(1,1), size(U, 1), &
               0.d0, W(nVV+1,1), size(W, 1))
    deallocate(H_mat)

    deallocate(Om_tmp)

    ! ---

  elseif(ispin .eq. 2) then

    call omp_set_max_active_levels(1)

    allocate(Om_tmp(nS))

    a0 = nOrb - nR - nO - 1

    eta2 = eta * eta
    do m = 1, nS
      Om_tmp(m) = Om(m) / (Om(m) * Om(m) + eta2)
    enddo

    ! ---

    !$OMP PARALLEL                                            &
    !$OMP DEFAULT(NONE)                                       &
    !$OMP PRIVATE(a, b, aa, ab, c, d, cd, i, j, ij, m, state, &
    !$OMP         bb, tmp_e, tmp_ab, chi, mat_tmp, H_mat)     &
    !$OMP SHARED(nC, nO, nOrb, nR, nS, n_states_diag, nVV,    &
    !$OMP        nOO, a0, eF, lambda, e, Om_tmp, rho, ERI, U, W)

    allocate(H_mat(nOO+nVV,a0))

    !$OMP DO SCHEDULE(GUIDED)
    do a = nO+1, nOrb-nR
      aa = a0 * (a - nO - 1) - (a - nO - 1) * (a - nO) / 2 - nO - 1
      do b = a+1, nOrb-nR
        ab = aa + b

        bb = b - a

        tmp_e = e(a) + e(b) - eF

        cd = 0
        do c = nO+1, nOrb-nR
          do d = c+1, nOrb-nR
            cd = cd + 1
    
            mat_tmp = lambda * (ERI(d,c,b,a) - ERI(d,c,a,b))

            chi = 0.d0
            do m = 1, nS
              chi = chi - Om_tmp(m) * (rho(m,a,c) * rho(m,b,d) - rho(m,a,d) * rho(m,b,c))
            enddo
            mat_tmp = mat_tmp + 4.d0 * lambda * chi

            if((a .eq. c) .and. (b .eq. d)) then
              mat_tmp = mat_tmp + tmp_e
            endif

            H_mat(cd,bb) = mat_tmp
          enddo ! d
        enddo ! c

        ij = nVV
        do i = nC+1, nO
          do j = i+1, nO
            ij = ij + 1

            mat_tmp = lambda * (ERI(j,i,b,a) - ERI(j,i,a,b))

            chi = 0.d0
            do m = 1, nS
              chi = chi - Om_tmp(m) * (rho(m,i,a) * rho(m,j,b) - rho(m,i,b) * rho(m,a,j))
            end do
            mat_tmp = mat_tmp + 4.d0 * lambda * chi

            H_mat(ij,bb) = -mat_tmp
          enddo ! j
        enddo ! i
      enddo ! b

      call dgemm("T", "N", nOrb-nR-a, n_states_diag, nOO+nVV, 1.d0, &
                 H_mat(1,1), size(H_mat, 1), U(1,1), size(U, 1),    &
                 0.d0, W(aa+a+1,1), size(W, 1))

    enddo ! a
    !$OMP END DO

    deallocate(H_mat)

    !$OMP END PARALLEL

    ! ---

    i0 = nO - nC - 1
    allocate(H_mat(nOO,nOO+nVV))

    !$OMP PARALLEL                                            &
    !$OMP DEFAULT(NONE)                                       &
    !$OMP PRIVATE(i, j, ii, ij, a, b, ab, k, l, kl, m, state, &
    !$OMP         tmp_e, tmp_ij, chi, mat_tmp)                &
    !$OMP SHARED(nC, nO, nOrb, nR, nS, nVV, i0,               &
    !$OMP        eF, lambda, e, Om_tmp, rho, ERI, U, H_mat)
    !$OMP DO SCHEDULE(GUIDED)
    do i = nC+1, nO
      ii = i0 * (i - nC - 1) - (i - nC - 1) * (i - nC) / 2 - nC - 1
      do j = i+1, nO
        ij = ii + j

        tmp_e = e(i) + e(j) - eF
  
        ab = 0
        do a = nO+1, nOrb-nR
          do b = a+1, nOrb-nR
            ab = ab + 1
    
            mat_tmp = lambda * (ERI(b,a,j,i) - ERI(b,a,i,j))

            chi = 0.d0
            do m = 1, nS
              chi = chi - Om_tmp(m) * (rho(m,i,a) * rho(m,j,b) - rho(m,i,b) * rho(m,a,j))
            enddo
            mat_tmp = mat_tmp + 4.d0 * lambda * chi

            H_mat(ij,ab) = mat_tmp
          enddo ! b
        enddo ! a

        kl = nVV
        do k = nC+1, nO
          do l = k+1, nO
            kl = kl + 1
   
            mat_tmp = lambda * (ERI(l,k,j,i) - ERI(l,k,i,j))

            chi = 0.d0
            do m = 1, nS
              chi = chi - Om_tmp(m) * (rho(m,i,k) * rho(m,j,l) - rho(m,i,l) * rho(m,j,k))
            enddo
            mat_tmp = mat_tmp + 4.d0 * lambda * chi

            if((i .eq. k) .and. (j .eq. l)) then
              mat_tmp = mat_tmp - tmp_e
            endif

            H_mat(ij,kl) = -mat_tmp
          enddo ! l
        enddo ! k
      enddo ! j
    enddo ! i
    !$OMP END DO
    !$OMP END PARALLEL

    call dgemm("N", "N", nOO, n_states_diag, nOO+nVV, 1.d0,    &
               H_mat(1,1), size(H_mat, 1), U(1,1), size(U, 1), &
               0.d0, W(nVV+1,1), size(W, 1))
    deallocate(H_mat)

    deallocate(Om_tmp)

    ! ---

  else

    print*, ' Error in ppLR_GW_HR_calc_batches'
    print*, ' ispin is not supported'
    print*, ' ispin = ', ispin
    stop

  endif

  return
end

! ---

