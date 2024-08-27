program test_python_call

    use pyscf_module
    implicit none

    character(len = 100) :: xyz, input_basis, unit_type
    integer :: charge, multiplicity, cartesian
    integer :: natoms, nalpha, nbeta
    double precision :: Enuc

    xyz = '/users/p18005/ammar/tmpdir/QuAcK/mol/H2O.xyz' // C_NULL_CHAR
    input_basis = 'sto-3g' // C_NULL_CHAR
    charge = 0
    multiplicity = 1
    unit_type = 'Angstrom' // C_NULL_CHAR
    cartesian = 0

    call call_mol_prop(xyz, input_basis, charge, multiplicity, unit_type, cartesian, &
                       natoms, nalpha, nbeta, Enuc)

    !call call_mol_prop(trim(adjustl(xyz)), trim(adjustl(input_basis)), charge, multiplicity, trim(adjustl(unit_type)), cartesian, &
    !                   natoms, nalpha, nbeta, Enuc)

    print *, "Number of atoms:", natoms
    print *, "Number of alpha electrons:", nalpha
    print *, "Number of beta electrons:", nbeta
    print *, "Nuclear energy:", Enuc

end program test_python_call


