module eigenvalue

    use global
    use constants
    use tools
    use spectral_function

    implicit none

contains

    !============================================================
    ! Linearized Eliashberg equations
    !
    ! Power iteration for largest eigenvalue
    !
    ! Tc occurs when:
    !
    ! lambdaMax = 1
    !============================================================

    subroutine solve_eigenvalue(x, im, lambdaMax)

        type(parameters), intent(in)   :: x
        type(matsubara), intent(inout) :: im

        real(dp), intent(out) :: lambdaMax

        integer :: n
        integer :: m
        integer :: nMats
        integer :: iter

        real(dp) :: kernel
        real(dp) :: denom

        real(dp) :: error

        real(dp), allocatable :: lambda_nm(:,:)

        real(dp), allocatable :: DeltaOld(:)
        real(dp), allocatable :: DeltaNew(:)

        !--------------------------------------------------------
        ! Matsubara size
        !--------------------------------------------------------

        nMats = x%nMatsubara

        !--------------------------------------------------------
        ! Allocate arrays
        !--------------------------------------------------------

        allocate(DeltaOld(0:nMats))
        allocate(DeltaNew(0:nMats))

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

        DeltaOld = 1.0_dp

        !========================================================
        ! Power iteration
        !========================================================

        do iter = 1, x%maxIter

            DeltaNew = 0.0_dp

            do n = 0, nMats

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
                    ! Linearized denominator
                    !--------------------------------------------

                    denom = abs(im%omega(m))

                    !--------------------------------------------
                    ! Linearized gap equation
                    !--------------------------------------------

                    DeltaNew(n) = DeltaNew(n) + &
                        pi * kB * x%T * &
                        (kernel - x%muStar) * &
                        DeltaOld(m) / denom

                end do

            end do

            !----------------------------------------------------
            ! Largest eigenvalue estimate
            !----------------------------------------------------

            lambdaMax = maxval(abs(DeltaNew))

            !----------------------------------------------------
            ! Normalize eigenvector
            !----------------------------------------------------

            if (lambdaMax > tiny) then

                DeltaNew = DeltaNew / lambdaMax

            end if

            !----------------------------------------------------
            ! Convergence
            !----------------------------------------------------

            error = norm_diff(DeltaNew, DeltaOld)

            if (error < x%tolerance) exit

            DeltaOld = DeltaNew

        end do

        !--------------------------------------------------------
        ! Store eigenvector
        !--------------------------------------------------------

        im%Delta = DeltaNew

        !========================================================
        ! Output
        !========================================================

        print *
        print *, '========================================='
        print *, ' Linearized Eliashberg solver'
        print *, '========================================='
        print *

        print *, 'Iterations       : ', iter
        print *, 'Largest eigenvalue : ', lambdaMax
        print *, 'Error            : ', error

        print *

        !--------------------------------------------------------
        ! Clean memory
        !--------------------------------------------------------

        deallocate(DeltaOld)
        deallocate(DeltaNew)

        if (allocated(lambda_nm)) then
            deallocate(lambda_nm)
        end if

    end subroutine solve_eigenvalue

end module eigenvalue