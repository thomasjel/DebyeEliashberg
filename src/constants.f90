module constants

    use global, only: dp

    implicit none

    !============================================================
    ! Mathematical constants
    !============================================================

    real(dp), parameter :: pi = 3.14159265358979323846_dp
    real(dp), parameter :: twopi = 2.0_dp * pi

    !============================================================
    ! Physical constants
    !============================================================

    ! Boltzmann constant [eV/K]
    real(dp), parameter :: kB = 8.617333262145e-5_dp

    ! Reduced Planck constant [eV*s]
    real(dp), parameter :: hbar = 6.582119569e-16_dp

    ! Electron volt to meV
    real(dp), parameter :: eV_to_meV = 1000.0_dp

    ! Kelvin to meV conversion
    real(dp), parameter :: K_to_meV = kB * eV_to_meV

    !============================================================
    ! Numerical constants
    !============================================================

    real(dp), parameter :: zero = 0.0_dp
    real(dp), parameter :: one  = 1.0_dp
    real(dp), parameter :: two  = 2.0_dp

    real(dp), parameter :: tiny = 1.0e-14_dp

end module constants