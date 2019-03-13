subroutine select_rung(rung,DFA)

! Select rung of Jacob's ladder

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: rung
  character(len=12),intent(in)  :: DFA

  select case (rung)

! Hartree calculation
    case(0)
      write(*,*) "* 0th rung of Jacob's ladder: Hartree calculation *"

! LDA functionals
    case(1)
      write(*,*) "* 1st rung of Jacob's ladder: local-density approximation (LDA)   *"

! GGA functionals
    case(2)
      write(*,*) "* 2nd rung of Jacob's ladder: generalized gradient approximation (GGA) *"

! meta-GGA functionals
    case(3)
      write(*,*) "* 3rd rung of Jacob's ladder: meta-GGA functional (MGGA) *"

! Hybrid functionals
    case(4)
      write(*,*) "* 4th rung of Jacob's ladder: hybrid functional *"

! Hartree-Fock calculation
    case(666)
      write(*,*) "* rung 666: Hartree-Fock calculation *"

! Default    
    case default
      write(*,*) "!!! rung not available !!!"
      stop

  end select

  ! Print selected functional

  write(*,*) '* You have selected the following functional: ',DFA,'        *'
  write(*,*) '*******************************************************************'
  write(*,*)


end subroutine select_rung
