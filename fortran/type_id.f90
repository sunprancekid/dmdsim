! AUTHOR : Matthew A Dorsey
! DATE : 2022-12-20
! FILENAME : type_id.f90
! PURPOSE : contains type id class and methods

module type_id
use dmdconstant
implicit none 

! define id class
type id 
    integer :: one ! group number id
    integer :: two ! particle number id
end type

contains 

! function that initializes id object with null pointers
! stored with id type module
type(id) function nullset()
    implicit none
    nullset%one = 0
    nullset%two = 0
end function nullset

! function that converts two number representing an id type
! to a molecular identifier
! stored with id type module
integer function id2mol(i, m, mer)
    implicit none
    integer, intent(in) :: i ! group id number 
    integer, intent(in) :: m ! particle id number
    integer, intent(in) :: mer ! number of particles in one group
    id2mol = (i-1)*mer + m
end function id2mol

! functrion that convers a molecular id to a type id
! stored with id type module 
type(id) function mol2id(i, mer)
    implicit none
    integer, intent(in) :: i ! molecular id
    integer, intent(in) :: mer ! number of particles in one group

    mol2id%one = ((i - 1) / mer) + 1
    mol2id%two = mod(i, mer)
    if (mol2id%two == 0) mol2id%two = mer
end function mol2id

! function that determines if two ids are the same
! stored with id type module
logical function idequiv(id1, id2)
    implicit none
    ! ** calling variables ***********************************
    type(id), intent(in) :: id1, id2
    ! ** local variables *************************************

    idequiv = (id1%one == id2%one) .and. (id1%two == id2%two)
end function idequiv

! function that generates a random group and sub-particle
! according to the number of groups in the simulation and 
! subparticles in the group
! stored with type id module
function random_group (g, p) result (rand)
    implicit none 
    integer, intent(in) :: g ! max number of groups
    integer, intent(in) :: p ! max number of sub particles per group 
    type(id) :: rand 
    real(kind=dbl) :: randnum 

    call random_number(randnum)
    rand%one = ceiling(randnum * real(g))

    call random_number(randnum)
    rand%two = ceiling(randnum * real(p))
end function random_group

end module type_id