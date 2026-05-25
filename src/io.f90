module io

    use global
    use constants

    implicit none

contains

    !============================================================
    ! Read command-line parameters
    !============================================================

    subroutine load_parameters(x)

        type(parameters), intent(inout) :: x

        integer :: i
        integer :: nargs
        integer :: pos

        character(len=256) :: arg
        character(len=128) :: key
        character(len=128) :: value

        nargs = command_argument_count()

        do i = 1, nargs

            call get_command_argument(i, arg)

            pos = index(arg, '=')

            if (pos == 0) cycle

            key   = adjustl(trim(arg(:pos-1)))
            value = adjustl(trim(arg(pos+1:)))

            select case (trim(key))

            !----------------------------------------------------
            ! Physical parameters
            !----------------------------------------------------

            case ('T')
                read(value, *) x%T

            case ('lambda')
                read(value, *) x%lambda

            case ('mu')
                read(value, *) x%muStar

            case ('wd')
                read(value, *) x%wd

            !----------------------------------------------------
            ! Numerical parameters
            !----------------------------------------------------

            case ('nMatsubara')
                read(value, *) x%nMatsubara

            case ('maxIter')
                read(value, *) x%maxIter

            case ('tolerance')
                read(value, *) x%tolerance

            case ('mixing')
                read(value, *) x%mixing

            !----------------------------------------------------
            ! Kernel selection
            !----------------------------------------------------

            case ('debye')
                read(value, *) x%useDebye

            case ('a2f')
                call load_a2F(trim(value), x)

            end select

        end do

    end subroutine load_parameters

    !============================================================
    ! Load alpha2F spectral function
    !
    ! File format:
    !
    ! omega[eV]   alpha2F
    !============================================================

    subroutine load_a2F(filename, x)

        character(len=*), intent(in) :: filename

        type(parameters), intent(inout) :: x

        integer :: i
        integer :: n
        integer :: unit
        integer :: ios

        real(dp) :: w
        real(dp) :: a

        !--------------------------------------------------------
        ! Count lines
        !--------------------------------------------------------

        n = 0

        open(newunit=unit, file=filename, &
             status='old', action='read')

        do

            read(unit, *, iostat=ios) w, a

            if (ios /= 0) exit

            n = n + 1

        end do

        rewind(unit)

        !--------------------------------------------------------
        ! Allocate arrays
        !--------------------------------------------------------

        call allocate_a2F(x, n)

        !--------------------------------------------------------
        ! Read data
        !--------------------------------------------------------

        do i = 1, n

            read(unit, *) &
                x%a2F_omega(i), &
                x%a2F(i)

        end do

        close(unit)

        !--------------------------------------------------------
        ! Activate spectral kernel
        !--------------------------------------------------------

        x%useA2F   = .true.
        x%useDebye = .false.

    end subroutine load_a2F

    !============================================================
    ! Print parameters
    !============================================================

    subroutine print_parameters(x)

        type(parameters), intent(in) :: x

        print *
        print *, '========================================='
        print *, ' DebyeEliashberg'
        print *, '========================================='
        print *

        print *, 'Temperature [K]      : ', x%T
        print *, 'Lambda               : ', x%lambda
        print *, 'mu*                  : ', x%muStar
        print *, 'Debye frequency [eV] : ', x%wd
        print *, 'Matsubara frequencies: ', x%nMatsubara

        print *

        if (x%useDebye) then

            print *, 'Kernel : Debye'

        end if

        if (x%useA2F) then

            print *, 'Kernel : alpha2F'
            print *, 'a2F points : ', x%nA2F

        end if

        print *

    end subroutine print_parameters

end module io