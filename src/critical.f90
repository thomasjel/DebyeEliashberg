program critical

    use global
    use constants
    use io
    use tools
    use eigenvalue

    implicit none

    type(parameters) :: x
    type(matsubara) :: im

    real(dp) :: lambdaMax

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
    ! Solve linearized Eliashberg equations
    !============================================================

    call solve_eigenvalue(x, im, lambdaMax)

    !============================================================
    ! Output
    !============================================================

    print *
    print *, '========================================='
    print *, ' Critical temperature analysis'
    print *, '========================================='
    print *

    print *, 'Temperature [K] : ', x%T
    print *, 'Largest eigenvalue : ', lambdaMax
    print *

    if (abs(lambdaMax - 1.0_dp) < 1.0e-2_dp) then

        print *, 'System is close to Tc'
        print *

    else if (lambdaMax > 1.0_dp) then

        print *, 'Superconducting state expected'
        print *

    else

        print *, 'Normal state expected'
        print *

    end if

end program critical