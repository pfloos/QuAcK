

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
  double precision, intent(in)  :: rho(nOrb,nOrb,nS), Om(nS)
  double precision, intent(in)  :: U(nOO+nVV,n_states_diag)
  double precision, intent(out) :: W(nOO+nVV,n_states_diag)

  integer                       :: i, j, ij, k, l, kl
  integer                       :: a, b, c, d, ab, cd
  integer                       :: m
  integer                       :: state
  double precision              :: mat_tmp, chi, eps

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
              chi = 0.d0
              do m = 1, nS
                eps = Om(m)**2 + eta**2
                chi = chi - rho(a,c,m) * rho(b,d,m) * Om(m) / eps &
                          - rho(a,d,m) * rho(b,c,m) * Om(m) / eps
              enddo
              mat_tmp = mat_tmp + 4.d0 * lambda * chi / dsqrt( (1.d0 + Kronecker_delta(a, b)) &
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

              chi = 0.d0
              do m = 1, nS
                eps = Om(m)**2 + eta**2
                chi = chi - rho(i,a,m) * rho(j,b,m) * Om(m) / eps &
                          - rho(i,b,m) * rho(a,j,m) * Om(m) / eps
              enddo
              mat_tmp = mat_tmp + 4.d0 * lambda * chi / dsqrt( (1.d0 + Kronecker_delta(a, b)) &
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
  
          ab = 0
          do a = nO+1, nOrb-nR
            do b = a, nOrb-nR
              ab = ab + 1
    
              mat_tmp = lambda * (ERI(a,b,i,j) + ERI(a,b,j,i)) / dsqrt( (1.d0 + Kronecker_delta(a, b)) &
                                                                      * (1.d0 + Kronecker_delta(i, j)))

              chi = 0.d0
              do m = 1, nS
                eps = Om(m)**2 + eta**2
                chi = chi - rho(i,a,m) * rho(j,b,m) * Om(m) / eps &
                          - rho(i,b,m) * rho(a,j,m) * Om(m) / eps
              enddo
              mat_tmp = mat_tmp + 4.d0 * lambda * chi / dsqrt( (1.d0 + Kronecker_delta(a, b)) &
                                                             * (1.d0 + Kronecker_delta(i, j)))

              W(ij,state) = W(ij,state) + mat_tmp * U(ab,state)
            enddo
          enddo

          kl = nVV
          do k = nC+1, nO
            do l = k, nO
              kl = kl + 1
    
              mat_tmp = - (e(i) + e(j) - eF) * Kronecker_delta(i, k) * Kronecker_delta(j, l) &
                      + lambda * (ERI(i,j,k,l) + ERI(i,j,l,k)) / dsqrt( (1.d0 + Kronecker_delta(i, j)) &
                                                                      * (1.d0 + Kronecker_delta(k, l)))

              chi = 0.d0
              do m = 1, nS
                eps = Om(m)**2 + eta**2
                chi = chi - rho(i,k,m) * rho(j,l,m) * Om(m) / eps &
                          - rho(i,l,m) * rho(j,k,m) * Om(m) / eps
              enddo
              mat_tmp = mat_tmp + 4.d0 * lambda * chi / dsqrt( (1.d0 + Kronecker_delta(i, j)) &
                                                             * (1.d0 + Kronecker_delta(k, l)))

              W(ij,state) = W(ij,state) - mat_tmp * U(kl,state)
            enddo
          enddo
  
        enddo ! state
      enddo ! j
    enddo ! i


  elseif(ispin .eq. 2) then

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

              chi = 0.d0
              do m = 1, nS
                eps = Om(m)**2 + eta**2
                chi = chi - rho(a,c,m) * rho(b,d,m) * Om(m) / eps &
                          + rho(a,d,m) * rho(b,c,m) * Om(m) / eps
              enddo
              mat_tmp = mat_tmp + 4.d0 * lambda * chi

              W(ab,state) = W(ab,state) + mat_tmp * U(cd,state)
            enddo
          enddo

          ij = nVV
          do i = nC+1, nO
            do j = i+1, nO
              ij = ij + 1
  
              mat_tmp = lambda * (ERI(a,b,i,j) - ERI(a,b,j,i))

              chi = 0.d0
              do m = 1, nS
                eps = Om(m)**2 + eta**2
                chi = chi - rho(i,a,m) * rho(j,b,m) * Om(m) / eps &
                          + rho(i,b,m) * rho(a,j,m) * Om(m) / eps
              end do
              mat_tmp = mat_tmp + 4.d0 * lambda * chi
  
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
  
          ab = 0
          do a = nO+1, nOrb-nR
            do b = a+1, nOrb-nR
              ab = ab + 1
    
              mat_tmp = lambda * (ERI(a,b,i,j) - ERI(a,b,j,i))

              chi = 0.d0
              do m = 1, nS
                eps = Om(m)**2 + eta**2
                chi = chi - rho(i,a,m) * rho(j,b,m) * Om(m) / eps &
                          + rho(i,b,m) * rho(a,j,m) * Om(m) / eps
              end do
              mat_tmp = mat_tmp + 4.d0 * lambda * chi

              W(ij,state) = W(ij,state) + mat_tmp * U(ab,state)
            enddo
          enddo

          kl = nVV
          do k = nC+1, nO
            do l = k+1, nO
              kl = kl + 1
    
              mat_tmp = - (e(i) + e(j) - eF) * Kronecker_delta(i, k) * Kronecker_delta(j, l) &
                      + lambda * (ERI(i,j,k,l) - ERI(i,j,l,k))

              chi = 0.d0
              do m = 1, nS
                eps = Om(m)**2 + eta**2
                chi = chi - rho(i,k,m) * rho(j,l,m) * Om(m) / eps &
                          + rho(i,l,m) * rho(j,k,m) * Om(m) / eps
              enddo
              mat_tmp = mat_tmp + 4.d0 * lambda * chi

              W(ij,state) = W(ij,state) - mat_tmp * U(kl,state)
            enddo
          enddo

        enddo ! state
      enddo ! j
    enddo ! i

  else

    print*, ' Error in ppLR_GW_HR_calc'
    print*, ' ispin is not supported'
    print*, ' ispin = ', ispin
    stop

  endif

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
  double precision, intent(in)  :: rho(nOrb,nOrb,nS), Om(nS)
  double precision, intent(out) :: H_diag(nOO+nVV)

  integer                       :: i, j, ij, k, l, kl
  integer                       :: a, b, c, d, ab, cd
  integer                       :: m
  double precision              :: chi, eps

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
            chi = 0.d0
            do m = 1, nS
              eps = Om(m)**2 + eta**2
              chi = chi - rho(a,c,m) * rho(b,d,m) * Om(m) / eps &
                        - rho(a,d,m) * rho(b,c,m) * Om(m) / eps
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
              chi = chi - rho(i,k,m) * rho(j,l,m) * Om(m) / eps &
                        - rho(i,l,m) * rho(j,k,m) * Om(m) / eps
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
              chi = chi - rho(a,c,m) * rho(b,d,m) * Om(m) / eps &
                        + rho(a,d,m) * rho(b,c,m) * Om(m) / eps
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
              chi = chi - rho(i,k,m) * rho(j,l,m) * Om(m) / eps &
                        + rho(i,l,m) * rho(j,k,m) * Om(m) / eps
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

  return
end

! ---

