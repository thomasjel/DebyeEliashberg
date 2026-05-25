module self_energy

    use global
    use constants
    use tools
    use spectral_function

    implicit none

contains

    !============================================================
    ! Solve isotropic Eliashberg equations
    !
    ! Imaginary-axis self-consistent solver
    !
    ! Solves:
    !
    !   Z_n
    !   Delta_n
    !
    !============================================================

    subroutine solve_eliashberg(x, im)

        type(parameters), intent(in)   :: x
        type(matsubara), intent(inout) :: im

        integer :: n
        integer :: m
        integer :: nMats
        integer :: iter

        real(dp) :: error

        real(dp) :: kernel
        real(dp) :: denom

        real(dp) :: Znew
        real(dp) :: Dnew

        real(dp), allocatable :: lambda_nm(:,:)

        real(dp), allocatable :: Zold(:)
        real(dp), allocatable :: Dold(:)

        !--------------------------------------------------------
        ! Matsubara size
        !--------------------------------------------------------

        nMats = x%nMatsubara

        !--------------------------------------------------------
        ! Allocate temporary arrays
        !--------------------------------------------------------

        allocate(Zold(0:nMats))
        allocate(Dold(0:nMats))

        !--------------------------------------------------------
        ! Build Matsubara frequencies
        !--------------------------------------------------------

        call build_matsubara(x, im)

        !--------------------------------------------------------
        ! Build spectral kernel matrix
        !--------------------------------------------------------

        if (x%useA2F) then

            call build_lambda_matrix(x, im, lambda_nm)

        end if

        !--------------------------------------------------------
        ! Initial guess
        !--------------------------------------------------------

        im%Z     = 1.0_dp
        im%Delta = 1.0e-4_dp

        !========================================================
        ! Self-consistency loop
        !========================================================

        do iter = 1, x%maxIter

            Zold = im%Z
            Dold = im%Delta

            !----------------------------------------------------
            ! Matsubara equations
            !----------------------------------------------------

            do n = 0, nMats

                Znew = 1.0_dp
                Dnew = 0.0_dp

                do m = 0, nMats

                    !--------------------------------------------
                    ! Kernel
                    !--------------------------------------------

                    if (x%useDebye) then

                        kernel = debye_kernel(x, n-m)

                    else if (x%useA2F) then

                        kernel = lambda_nm(n,m)

                    else

                        kernel = 0.0_dp

                    end if

                    !--------------------------------------------
                    ! Denominator
                    !--------------------------------------------

                    denom = sqrt( &
                        im%omega(m)**2 + &
                        im%Delta(m)**2 )

                    !--------------------------------------------
                    ! Z equation
                    !--------------------------------------------

                    Znew = Znew + &
                        (pi * kB * x%T / im%omega(n)) * &
                        kernel * &
                        im%omega(m) / denom

                    !--------------------------------------------
                    ! Gap equation
                    !--------------------------------------------

                    Dnew = Dnew + &
                        pi * kB * x%T * &
                        (kernel - x%muStar) * &
                        im%Delta(m) / denom

                end do

                !--------------------------------------------
                ! Normalize gap equation
                !--------------------------------------------

                Dnew = Dnew / Znew

                !--------------------------------------------
                ! Mixing
                !--------------------------------------------

                im%Z(n) = &
                    x%mixing * Znew + &
                    (1.0_dp - x%mixing) * Zold(n)

                im%Delta(n) = &
                    x%mixing * Dnew + &
                    (1.0_dp - x%mixing) * Dold(n)

            end do

            !----------------------------------------------------
            ! Convergence
            !----------------------------------------------------

            error = max( &
                norm_diff(im%Z, Zold), &
                norm_diff(im%Delta, Dold) )

            if (error < x%tolerance) exit

        end do

        !========================================================
        ! Output
        !========================================================

        print *
        print *, '========================================='
        print *, ' Eliashberg converged'
        print *, '========================================='
        print *

        print *, 'Iterations : ', iter
        print *, 'Error      : ', error
        print *

        print *, 'Delta(0) [meV] : ', &
            im%Delta(0) * eV_to_meV

        print *

        !--------------------------------------------------------
        ! Clean memory
        !--------------------------------------------------------

        deallocate(Zold)
        deallocate(Dold)

        if (allocated(lambda_nm)) then
            deallocate(lambda_nm)
        end if

    end subroutine solve_eliashberg

end module self_energy