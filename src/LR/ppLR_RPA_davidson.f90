

! ---

subroutine ppLR_RPA_HR_calc(ispin, nOrb, nC, nO, nR, nOO, nVV, lambda, e, eF, n_states_diag, ERI, U, W)

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

subroutine ppLR_RPA_H_diag(ispin, nOrb, nC, nO, nR, nOO, nVV, lambda, e, eF, ERI, H_diag)

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

