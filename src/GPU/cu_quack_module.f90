module cu_quack_module

  use, intrinsic :: iso_c_binding

  implicit none

!#ifdef USE_GPU
!  interface
!    subroutine ph_drpa_tda_sing(nO, nBas, nS, eps, ERI, &
!        Omega, X) bind(C, name = "ph_drpa_tda_sing")
!
!      import c_int, c_double
!      integer(c_int), intent(in), value :: nO, nBas, nS
!      real(c_double), intent(in)        :: eps(nBas)
!      real(c_double), intent(in)        :: ERI(nBas,nBas,nBas,nBas)
!      real(c_double), intent(out)       :: Omega(nS)
!      real(c_double), intent(out)       :: X(nS,nS)
!
!    end subroutine ph_drpa_tda_sing
!  end interface
!#else
!  interface
!    subroutine ph_drpa_tda_sing(nO, nBas, nS, eps, ERI, Omega, X)
!      integer, intent(in)           :: nO, nBas, nS
!      double precision, intent(in)  :: eps(nBas)
!      double precision, intent(in)  :: ERI(nBas,nBas,nBas,nBas)
!      double precision, intent(out) :: Omega(nS)
!      double precision, intent(out) :: X(nS,nS)
!    end subroutine ph_drpa_tda_sing
!  end interface
!#endif

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

  end interface

  ! ---

end module cu_quack_module


