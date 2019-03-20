subroutine GenerateShell(atot,nShellFunction,ShellFunction)

  implicit none

! Input variables

  integer,intent(in)   :: atot,nShellFunction

! Local variables

  integer              :: ax,ay,az,ia

! Output variables

  integer,intent(out)  :: ShellFunction(nShellFunction,3)

  ia = 0
  do ax=atot,0,-1
    do az=0,atot
      ay = atot - ax - az
      if(ay >= 0) then
        ia = ia + 1
        ShellFunction(ia,1) = ax
        ShellFunction(ia,2) = ay
        ShellFunction(ia,3) = az
      end if
    end do
  end do

end subroutine GenerateShell
