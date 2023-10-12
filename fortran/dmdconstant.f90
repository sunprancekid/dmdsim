! AUTHOR : Matthew A Dorsey
! DATE : 2022-12-20
! FILENAME : dmdconstant.f90
! PURPOSE : contains constants and essential methods for dmd simulations

module dmdconstant
implicit none

! ** constants ***********************************************
integer, parameter :: dbl = selected_real_kind (32) ! integer which determines precision of real numbers
real(kind=dbl), parameter :: pi = 3.141592653589793238
real(kind=dbl), parameter :: twopi = 2. * pi
real(kind=dbl), parameter :: halfpi = pi / 2.
real(kind=dbl), parameter :: bigtime = 1e10 ! unreasonably large time
! ** colors (rbg format) *************************************
character(len=12), parameter :: red = ' 1 0.15 0.15'
character(len=12), parameter :: blue = ' 0.1 1 0.3'
character(len=12), parameter :: green = ' 0 0 1'
character(len=12), parameter :: orange = ' 1 0 0.64'
character(len=12), parameter :: purple = ' 1 1 0'
character(len=15), parameter :: white = ' 0.75 0.75 0.75'

contains


! This function generates a random integer based on a 
! Poisson distribution. The poisson variate is generated 
! using an algorithm that simulates a poisson distrubution
! directly.
!
! SOURCE: https://hpaulkeeler.com/simulating-poisson-random-variables-direct-method/
integer function poissondist (lambda)
    implicit none
    real(kind=dbl), intent(in) :: lambda ! parameter used in poisson PDF
    ! average number of times an event occurs per unit time
    real(kind=dbl) :: exp_lambda
    real(kind=dbl) :: randUni ! random uniform variable
    real(kind=dbl) :: prodUni ! product of random uniform variables
    real(kind=dbl) :: randPoisson ! random poisson variate returned by function

    ! initialize variables
    exp_lambda = exp(-lambda)
    randPoisson = -1
    prodUni = 1
    do 
        call random_number(randUni)
        prodUni = prodUni * randUni
        randPoisson = randPoisson + 1
        if (exp_lambda > prodUni) exit
    end do 
    poissondist = randPoisson
end function poissondist

end module dmdconstant