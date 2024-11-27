module cu_quack_module

  use, intrinsic :: iso_c_binding

  implicit none

  ! ---

  interface

    subroutine ph_drpa_tda(nO, nBas, nS, eps, ERI, &
                           Omega, X) bind(C, name = "ph_drpa_tda")

      import c_int, c_double
      integer(c_int), intent(in), value :: nO, nBas, nS
      real(c_double), intent(in)        :: eps(nBas)
      real(c_double), intent(in)        :: ERI(nBas,nBas,nBas,nBas)
      real(c_double), intent(out)       :: Omega(nS)
      real(c_double), intent(out)       :: X(nS,nS)

    end subroutine ph_drpa_tda

  end interface

  ! ---

  contains

    subroutine cu_quack_module_test()
        implicit none
        print*, ' hello from cu_quack_module'
    end subroutine cu_quack_module_test

  ! ---

end module cu_quack_module


