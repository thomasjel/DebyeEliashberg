module spectral_function

    use global
    use constants
    use tools

    implicit none

contains

    !============================================================
    ! Compute spectral kernel
    !
    ! lambda(n,m) =
    !
    ! 2 ∫ dΩ α²F(Ω) Ω /
    !   [ Ω² + (ω_n-ω_m)² ]
    !
    !============================================================

    pure function spectral_kernel(x, wn) result(kernel)

        type(parameters), intent(in) :: x

        real(dp), intent(in) :: wn

        real(dp) :: kernel

        integer :: i
        integer :: nPoints

        real(dp), allocatable :: integrand(:)

        !--------------------------------------------------------
        ! Number of alpha2F points
        !--------------------------------------------------------

        nPoints = x%nA2F

        allocate(integrand(nPoints))

        !--------------------------------------------------------
        ! Build integrand
        !--------------------------------------------------------

        do i = 1, nPoints

            integrand(i) = &
                2.0_dp * &
                x%a2F(i) * &
                x%a2F_omega(i) / &
                (x%a2F_omega(i)**2 + wn**2)

        end do

        !--------------------------------------------------------
        ! Integrate
        !--------------------------------------------------------

        kernel = trapz( &
            x%a2F_omega, &
            integrand, &
            nPoints)

        deallocate(integrand)

    end function spectral_kernel

    !============================================================
    ! Build lambda(n,m) matrix
    !============================================================

    subroutine build_lambda_matrix(x, im, lambda_nm)

        type(parameters), intent(in) :: x
        type(matsubara), intent(in)  :: im

        real(dp), allocatable, intent(out) :: lambda_nm(:,:)

        integer :: n
        integer :: m
        integer :: nMats

        real(dp) :: wn

        !--------------------------------------------------------
        ! Matsubara size
        !--------------------------------------------------------

        nMats = x%nMatsubara

        allocate(lambda_nm(0:nMats,0:nMats))

        !--------------------------------------------------------
        ! Build kernel matrix
        !--------------------------------------------------------

        do n = 0, nMats

            do m = 0, nMats

                wn = abs(im%omega(n) - im%omega(m))

                lambda_nm(n,m) = &
                    spectral_kernel(x, wn)

            end do

        end do

    end subroutine build_lambda_matrix

    !============================================================
    ! Compute total electron-phonon coupling
    !
    ! lambda =
    !
    ! 2 ∫ dΩ α²F(Ω) / Ω
    !
    !============================================================

    pure function compute_lambda(x) result(lambda)

        type(parameters), intent(in) :: x

        real(dp) :: lambda

        integer :: i
        integer :: nPoints

        real(dp), allocatable :: integrand(:)

        !--------------------------------------------------------
        ! Number of alpha2F points
        !--------------------------------------------------------

        nPoints = x%nA2F

        allocate(integrand(nPoints))

        !--------------------------------------------------------
        ! Build integrand
        !--------------------------------------------------------

        do i = 1, nPoints

            if (x%a2F_omega(i) > tiny) then

                integrand(i) = &
                    2.0_dp * &
                    x%a2F(i) / &
                    x%a2F_omega(i)

            else

                integrand(i) = 0.0_dp

            end if

        end do

        !--------------------------------------------------------
        ! Integrate
        !--------------------------------------------------------

        lambda = trapz( &
            x%a2F_omega, &
            integrand, &
            nPoints)

        deallocate(integrand)

    end function compute_lambda

end module spectral_function