subroutine complex_sort_eigenvalues_RPA(n, evals, evecs)
  implicit none
  integer, intent(in)           :: n
  complex*16, intent(inout)     :: evals(n)
  complex*16, intent(inout)     :: evecs(n, n)
  integer                       :: i, j, k, nhalf, counter
  complex*16                    :: temp_val
  complex*16                    :: temp_vec(n)
  double precision              :: etanorm
  double precision              :: threshold = 1d-6

  ! -------------------------------
  ! Step 1: initial bubble sort
  ! Sort in descending order by real part, if two real parts are the same than first the one with lower imag part
  ! -------------------------------
  do i = 1, n-1
   do j = i+1, n

      if (abs(real(evals(i)) - real(evals(j))) > threshold) then
         ! Primary sort: real part ascending
         if (real(evals(i)) < real(evals(j))) then
            ! swap
            temp_val = evals(i)
            evals(i) = evals(j)
            evals(j) = temp_val

            temp_vec = evecs(:, i)
            evecs(:, i) = evecs(:, j)
            evecs(:, j) = temp_vec
         end if
      else
         ! Secondary sort: imaginary part descending
         if (aimag(evals(i)) > aimag(evals(j))) then
            ! swap
            temp_val = evals(i)
            evals(i) = evals(j)
            evals(j) = temp_val

            temp_vec = evecs(:, i)
            evecs(:, i) = evecs(:, j)
            evecs(:, j) = temp_vec
         end if
      end if

   end do
end do

  ! -------------------------------
  ! Step 2: classify first half
  !  - eta-norm > 0  -> excitation
  !  - eta-norm < 0  -> swap
  !  - if Re(omega) = 0 -> use imaginary part (this criteria has to be checked)
  ! -------------------------------
  nhalf = n/2
  counter = 0
  print *,"OmOmminus"
  call complex_vecout(n,evals)

  do i = 1, nhalf

     etanorm = sum(abs(evecs(1:nhalf,i))**2) - &
               sum(abs(evecs(nhalf+1:n,i))**2)
     print *, i, etanorm
     ! tolerance for zero real part
     if (abs(real(evals(i))) < threshold) then

        ! purely imaginary pair
        if (aimag(evals(i)) > 0d0) then
           ! positive imaginary part -> move to second half
           j = nhalf + i

           temp_val = evals(i)
           evals(i) = evals(j)
           evals(j) = temp_val

           temp_vec = evecs(:, i)
           evecs(:, i) = evecs(:, j)
           evecs(:, j) = temp_vec

        end if

     else

        ! normal case: classify by eta-norm
        if (abs(etanorm) < threshold) then
           print *, 'Mode with vanishing eta-norm has been found !'
           print *, 'Omega:',i, evals(i)

        else if (etanorm < threshold) then

           j = nhalf + i

           temp_val = evals(i)
           evals(i) = evals(j)
           evals(j) = temp_val

           temp_vec = evecs(:, i)
           evecs(:, i) = evecs(:, j)
           evecs(:, j) = temp_vec

           counter = counter + 1

        end if

     end if

  end do

  if(counter /= 0) then
    print *, "Attention:", counter, "excitations with negative excitation energies have been found!" 
  end if

  ! -------------------------------
  ! Step 3: separate positives and negatives
  ! -------------------------------
  ! Assuming after step 2, first nhalf are excitations (positive eta norm) and last nhalf are de-excitations (negative eta norm)
  
  ! Sort positive excitations ascending
  do i = 1, nhalf-1
     do j = i+1, nhalf
        if (real(evals(i)) > real(evals(j))) then
          temp_val = evals(i)
          evals(i) = evals(j)
          evals(j) = temp_val
          temp_vec = evecs(:, i)
          evecs(:, i) = evecs(:, j)
          evecs(:, j) = temp_vec
        end if
     end do
  end do

  ! Sort negative de-excitations descending
  do i = nhalf+1, n-1
     do j = i+1, n
        if (real(evals(i)) < real(evals(j))) then
          temp_val = evals(i)
          evals(i) = evals(j)
          evals(j) = temp_val
          temp_vec = evecs(:, i)
          evecs(:, i) = evecs(:, j)
          evecs(:, j) = temp_vec
        end if
     end do
  end do
  print *,"OmOmminus"
  call complex_vecout(n,evals)
!  print *,'eta norm of omomminus'
!  do i=1,n
!    print *,i,sum(abs(evecs(1:nhalf,i))**2) - sum(abs(evecs(nhalf+1:n,i))**2)
!    call complex_vecout(n,evecs(:,i))
!  end do

end subroutine
