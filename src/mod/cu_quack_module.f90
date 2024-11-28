module cu_quack_module

  use, intrinsic :: iso_c_binding

  implicit none

  ! ---

  interface

    subroutine ph_drpa_tda_sing(nO, nBas, nS, eps, ERI, &
        Omega, X) bind(C, name = "ph_drpa_tda_sing")

      import c_int, c_double
      integer(c_int), intent(in), value :: nO, nBas, nS
      real(c_double), intent(in)        :: eps(nBas)
      real(c_double), intent(in)        :: ERI(nBas,nBas,nBas,nBas)
      real(c_double), intent(out)       :: Omega(nS)
      real(c_double), intent(out)       :: X(nS,nS)

    end subroutine ph_drpa_tda_sing

    ! ---

    subroutine ph_drpa_tda_trip(nO, nBas, nS, eps, ERI, &
        Omega, X) bind(C, name = "ph_drpa_tda_trip")

      import c_int, c_double
      integer(c_int), intent(in), value :: nO, nBas, nS
      real(c_double), intent(in)        :: eps(nBas)
      real(c_double), intent(in)        :: ERI(nBas,nBas,nBas,nBas)
      real(c_double), intent(out)       :: Omega(nS)
      real(c_double), intent(out)       :: X(nS,nS)

    end subroutine ph_drpa_tda_trip

    ! ---

    subroutine ph_drpa_sing(nO, nBas, nS, eps, ERI, &
        Omega, X) bind(C, name = "ph_drpa_sing")

      import c_int, c_double
      integer(c_int), intent(in), value :: nO, nBas, nS
      real(c_double), intent(in)        :: eps(nBas)
      real(c_double), intent(in)        :: ERI(nBas,nBas,nBas,nBas)
      real(c_double), intent(out)       :: Omega(nS)
      real(c_double), intent(out)       :: X(nS,nS)

    end subroutine ph_drpa_sing

    ! ---

    subroutine ph_drpa_trip(nO, nBas, nS, eps, ERI, &
        Omega, X) bind(C, name = "ph_drpa_trip")

      import c_int, c_double
      integer(c_int), intent(in), value :: nO, nBas, nS
      real(c_double), intent(in)        :: eps(nBas)
      real(c_double), intent(in)        :: ERI(nBas,nBas,nBas,nBas)
      real(c_double), intent(out)       :: Omega(nS)
      real(c_double), intent(out)       :: X(nS,nS)

    end subroutine ph_drpa_trip

    ! ---

  end interface

  ! ---

  contains

    subroutine cu_quack_module_test()
        implicit none
        print*, ' hello from cu_quack_module'
    end subroutine cu_quack_module_test

  ! ---

end module cu_quack_module


