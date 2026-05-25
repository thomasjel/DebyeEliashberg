program gap

    use global
    use constants
    use io
    use tools
    use self_energy

    implicit none

    type(parameters) :: x
    type(matsubara) :: im

    integer :: n
    integer :: unit

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
    ! Save Matsubara gap
    !============================================================

    open(newunit=unit, file='gap.dat', status='replace')

    write(unit, '(A)') '# n omega_n[eV] Delta_n[meV] Z_n'

    do n = 0, x%nMatsubara

        write(unit,'(I6,3E20.10)') &
            n, &
            im%omega(n), &
            im%Delta(n) * eV_to_meV, &
            im%Z(n)

    end do

    close(unit)

    !============================================================
    ! Final output
    !============================================================

    print *
    print *, '========================================='
    print *, ' Gap calculation completed'
    print *, '========================================='
    print *

    print *, 'Output file : gap.dat'
    print *, 'Delta(0) [meV] : ', &
        im%Delta(0) * eV_to_meV

    print *

end program gap