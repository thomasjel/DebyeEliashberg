module global

    implicit none

    !============================================================
    ! Precision
    !============================================================

    integer, parameter :: dp = selected_real_kind(15, 307)

    !============================================================
    ! Parameters structure
    !============================================================

    type parameters

        !--------------------------------------------------------
        ! Physical parameters
        !--------------------------------------------------------

        real(dp) :: T        = 10.0_dp

        real(dp) :: lambda   = 1.0_dp
        real(dp) :: muStar   = 0.10_dp

        real(dp) :: wd       = 0.020_dp

        !--------------------------------------------------------
        ! Matsubara parameters
        !--------------------------------------------------------

        integer :: nMatsubara = 512

        !--------------------------------------------------------
        ! Self-consistency parameters
        !--------------------------------------------------------

        integer  :: maxIter   = 10000

        real(dp) :: tolerance = 1.0e-10_dp
        real(dp) :: mixing    = 0.5_dp

        !--------------------------------------------------------
        ! Kernel selection
        !--------------------------------------------------------

        logical :: useDebye = .true.
        logical :: useA2F   = .false.

        !--------------------------------------------------------
        ! alpha2F spectral function
        !--------------------------------------------------------

        integer :: nA2F = 0

        real(dp), allocatable :: a2F_omega(:)
        real(dp), allocatable :: a2F(:)

    end type parameters

    !============================================================
    ! Matsubara quantities
    !============================================================

    type matsubara

        real(dp), allocatable :: omega(:)

        real(dp), allocatable :: Z(:)
        real(dp), allocatable :: Delta(:)

    end type matsubara

contains

    !============================================================
    ! Allocate Matsubara arrays
    !============================================================

    subroutine allocate_matsubara(x, im)

        type(parameters), intent(in)   :: x
        type(matsubara), intent(inout) :: im

        integer :: N

        N = x%nMatsubara

        allocate(im%omega(0:N))
        allocate(im%Z(0:N))
        allocate(im%Delta(0:N))

        im%omega = 0.0_dp
        im%Z     = 1.0_dp
        im%Delta = 1.0e-4_dp

    end subroutine allocate_matsubara

    !============================================================
    ! Allocate alpha2F arrays
    !============================================================

    subroutine allocate_a2F(x, n)

        type(parameters), intent(inout) :: x

        integer, intent(in) :: n

        x%nA2F = n

        allocate(x%a2F_omega(n))
        allocate(x%a2F(n))

        x%a2F_omega = 0.0_dp
        x%a2F       = 0.0_dp

    end subroutine allocate_a2F

end module global