module cu_quack_module

  use, intrinsic :: iso_c_binding

  implicit none

  interface

    ! ---

    subroutine ph_drpa(nO, nBas, eps, ERI) bind(C, name = "cutc_int")

      import c_int, c_double
      integer(c_int), intent(in), value :: nO, nBas
      real(c_double), intent(in)        :: eps(nBas)
      real(c_double), intent(in)        :: ERI(nBas,nBas,nBas,nBas)

    end subroutine ph_drpa

    ! ---

  end interface

end module cu_quack_module


