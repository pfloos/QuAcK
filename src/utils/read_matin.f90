subroutine read_matin(m, n, A, fpath)
    implicit none
    integer, intent(in)            :: m, n
    double precision, intent(out) :: A(m, n)
    character(len=*), intent(in)  :: fpath
    integer :: i, j, unit, ios

    open(newunit=unit, file=fpath, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, 'Error opening file: ', fpath
        return
    end if

    do i = 1, m
        read(unit, *, iostat=ios) (A(i, j), j=1, n)
        if (ios /= 0) then
            print *, 'Error reading row ', i
            close(unit)
            return
        end if
    end do

    close(unit)
end subroutine

