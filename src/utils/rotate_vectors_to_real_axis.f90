subroutine rotate_vectors_to_real_axis(n, m, V)
    implicit none
    integer, intent(in) :: n        ! length of eigenvectors
    integer, intent(in) :: m        ! number of eigenvectors
    complex*16, intent(inout) :: V(n,m)

    integer :: i, k
    double precision :: phi, threshold
    complex*16 :: phase

    threshold = 1.0d-10

    do i = 1, m

        ! Find index of largest component (for stability)
        k = maxloc(abs(V(:,i)), dim=1)

        if (abs(V(k,i)) < threshold) cycle

        ! Remove global phase
        phi = atan2(aimag(V(k,i)), real(V(k,i)))
        phase = cmplx(cos(-phi), sin(-phi), kind=8)
        V(:,i) = V(:,i) * phase
    
    end do

end subroutine
