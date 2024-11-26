module cu_quack_module

  use, intrinsic :: iso_c_binding

  implicit none

  ! ---

  interface

    subroutine ph_drpa(nO, nBas, nS, eps, ERI, &
                       Omega, XpY, XmY) bind(C, name = "ph_drpa")

      import c_int, c_double
      integer(c_int), intent(in), value :: nO, nBas, nS
      real(c_double), intent(in)        :: eps(nBas)
      real(c_double), intent(in)        :: ERI(nBas,nBas,nBas,nBas)
      real(c_double), intent(out)       :: Omega(nS)
      real(c_double), intent(out)       :: XpY(nS,nS)
      real(c_double), intent(out)       :: XmY(nS,nS)

    end subroutine ph_drpa

  end interface

  ! ---

  contains

    subroutine cu_quack_module_test()
        implicit none
        print*, ' hello from mod_test'
    end subroutine cu_quack_module_test

  ! ---

end module cu_quack_module


