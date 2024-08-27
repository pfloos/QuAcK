module pyscf_module

    use, intrinsic :: iso_c_binding

    implicit none

    interface

        subroutine call_mol_prop(xyz, input_basis, charge, multiplicity, unit_type, cartesian, &
                                 natoms, nalpha, nbeta, Enuc) bind(C, name="call_mol_prop")

            import c_char, c_int, c_double, c_ptr

            character(kind=c_char), intent(in) :: xyz(*)
            character(kind=c_char), intent(in) :: input_basis(*)
            integer(c_int), intent(in), value  :: charge
            integer(c_int), intent(in), value  :: multiplicity
            character(kind=c_char), intent(in) :: unit_type(*)
            integer(c_int), intent(in), value  :: cartesian
            integer(c_int), intent(out)        :: natoms
            integer(c_int), intent(out)        :: nalpha
            integer(c_int), intent(out)        :: nbeta
            real(c_double), intent(out)        :: Enuc
        end subroutine call_mol_prop

    end interface

end module pyscf_module

