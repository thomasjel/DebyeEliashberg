program tc

    use global
    use constants
    use io
    use tools
    use self_energy

    implicit none

    type(parameters) :: x
    type(matsubara) :: im

    !============================================================
    ! Load parameters
    !============================================================

    call load_parameters(x)

    !============================================================
    ! Allocate Matsubara arrays
    !============================================================

    call allocate_matsubara(x, im)

    !============================================================
    ! Print parameters
    !============================================================

    call print_parameters(x)

    !============================================================
    ! Solve Eliashberg equations
    !============================================================

    call solve_eliashberg(x, im)

    !============================================================
    ! Final output
    !============================================================

    print *
    print *, '========================================='
    print *, ' Results'
    print *, '========================================='
    print *

    print *, 'Temperature [K] : ', x%T
    print *, 'Delta(0) [meV]  : ', &
        im%Delta(0) * eV_to_meV

    print *

end program tc