module tools

    use global
    use constants

    implicit none

contains

    !============================================================
    ! Build Matsubara frequencies
    !
    ! omega_n = (2n+1) pi kB T
    !============================================================

    subroutine build_matsubara(x, im)

        type(parameters), intent(in)   :: x
        type(matsubara), intent(inout) :: im

        integer :: n

        do n = 0, x%nMatsubara

            im%omega(n) = &
                (2*n + 1) * pi * kB * x%T

        end do

    end subroutine build_matsubara

    !============================================================
    ! Debye kernel
    !
    ! lambda_nm =
    !
    ! lambda / [1 + ((n-m)/nD)^2]
    !
    !============================================================

    pure function debye_kernel(x, dn) result(kernel)

        type(parameters), intent(in) :: x

        integer, intent(in) :: dn

        real(dp) :: kernel

        real(dp) :: nD

        nD = x%wd / (2.0_dp * pi * kB * x%T)

        kernel = x%lambda / &
            (1.0_dp + (real(dn,dp)/nD)**2)

    end function debye_kernel

    !============================================================
    ! Max norm difference
    !
    ! Used for convergence
    !============================================================

    pure function norm_diff(a, b) result(diff)

        real(dp), intent(in) :: a(:)
        real(dp), intent(in) :: b(:)

        real(dp) :: diff

        diff = maxval(abs(a - b))

    end function norm_diff

    !============================================================
    ! Trapezoidal integration
    !============================================================

    pure function trapz(x, y, n) result(integral)

        integer, intent(in) :: n

        real(dp), intent(in) :: x(n)
        real(dp), intent(in) :: y(n)

        real(dp) :: integral

        integer :: i

        integral = 0.0_dp

        do i = 1, n-1

            integral = integral + &
                0.5_dp * (x(i+1)-x(i)) * &
                (y(i+1)+y(i))

        end do

    end function trapz

    !============================================================
    ! Fermi-Dirac distribution
    !============================================================

    pure function fermi(E, T) result(f)

        real(dp), intent(in) :: E
        real(dp), intent(in) :: T

        real(dp) :: f

        f = 1.0_dp / &
            (exp(E/(kB*T)) + 1.0_dp)

    end function fermi

    !============================================================
    ! Bose-Einstein distribution
    !============================================================

    pure function bose(E, T) result(b)

        real(dp), intent(in) :: E
        real(dp), intent(in) :: T

        real(dp) :: b

        b = 1.0_dp / &
            (exp(E/(kB*T)) - 1.0_dp)

    end function bose

    !============================================================
    ! Linear interpolation
    !============================================================

    pure function interpolate(x, y, n, x0) result(y0)

        integer, intent(in) :: n

        real(dp), intent(in) :: x(n)
        real(dp), intent(in) :: y(n)

        real(dp), intent(in) :: x0

        real(dp) :: y0

        integer :: i

        y0 = 0.0_dp

        !--------------------------------------------------------
        ! Boundary conditions
        !--------------------------------------------------------

        if (x0 <= x(1)) then

            y0 = y(1)
            return

        end if

        if (x0 >= x(n)) then

            y0 = y(n)
            return

        end if

        !--------------------------------------------------------
        ! Linear interpolation
        !--------------------------------------------------------

        do i = 1, n-1

            if (x0 >= x(i) .and. x0 <= x(i+1)) then

                y0 = y(i) + &
                    (y(i+1)-y(i)) * &
                    (x0-x(i)) / &
                    (x(i+1)-x(i))

                return

            end if

        end do

    end function interpolate

    !============================================================
    ! Save vector to file
    !============================================================

    subroutine save_vector(filename, x, y, n)

        character(len=*), intent(in) :: filename

        integer, intent(in) :: n

        real(dp), intent(in) :: x(n)
        real(dp), intent(in) :: y(n)

        integer :: i
        integer :: unit

        open(newunit=unit, file=filename, status='replace')

        do i = 1, n

            write(unit,'(2E20.10)') x(i), y(i)

        end do

        close(unit)

    end subroutine save_vector

end module tools