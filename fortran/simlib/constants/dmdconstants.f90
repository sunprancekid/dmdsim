! AUTHOR : Matthew Dorsey
! DATE : 2022-07-22
! FILENAME : dmdconstants.f90
! PURPOSE : This module contains consants that are used in other DMD modules

module dmdconstants
implicit none

! simulation constants
integer, parameter :: dim2 = 2 ! two dimensions
integer, parameter :: dim3 = 3 ! three dimensions
integer, parameter :: dbl = selected_real_kind(16) ! double precision, used for real numbers
integer, parameter :: quad = selected_real_kind(32) ! quad precision, used for real number
real(kind=quad), parameter :: tol = 0.001 ! tolerance used to as acceptance criteria for numerical integration
real, parameter :: pi = 4.D0*DATAN(1.D0) ! defines PI 
real, parameter :: halfpi = pi / 2. ! half of pi
real, parameter :: twopi = pi * 2. ! two times pi

! ** colors (rbg format) *************************************
character(len=12), parameter :: red = ' 1 0.15 0.15'
character(len=12), parameter :: blue = ' 0.1 1 0.3'
character(len=12), parameter :: green = ' 0 0 1'
character(len=12), parameter :: orange = ' 1 0 0.64'
character(len=12), parameter :: purple = ' 1 1 0'
character(len=15), parameter :: white = ' 0.75 0.75 0.75'

! ** I / O ****************************************************
integer, parameter :: saveiounit = 11 
integer, parameter :: simiounit = 12
integer, parameter :: coorsphiounit = 13
integer, parameter :: coorsquiounit = 14
integer, parameter :: reportiounit = 15
integer, parameter :: annealiounit = 16
integer, parameter :: opiounit = 17
integer, parameter :: equiliounit = 18

end module dmdconstants