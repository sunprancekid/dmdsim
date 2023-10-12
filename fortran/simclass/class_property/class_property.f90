! AUTHOR : Matthew Dorsey
! DATE : 2022-07-22
! FILENAME : class_property.f90
! PURPOSE : used to track the values of properties over the course of molecular simulations

module class_property
use dmdconstants
implicit none 

! property class description
type :: property
    private 
    real(kind=quad) :: val, sum, sum2
    integer :: count
end type property


contains 


! CONSTRUCTOR

! construct property object
function construct_property () result (p)
    type(property) :: p 
    p%val = 0.
    p%sum = 0.
    p%sum2 = 0.
    p%count = 0
end function construct_property


! GETTERS AND SETTERS

! set property to value passed as parameter 
subroutine set_property_value (p, val)
    type(property), intent(out) :: p
    real(kind=quad), intent(in) :: val
    p%val = val
end subroutine set_property_value

! get property average
function get_property_average (p) result (avg)
    type(property), intent(in) :: p 
    real(kind=quad) :: avg
    avg = p%sum / p%count
end function get_property_average

! get stadard deviation of property
function get_property_stddev (p) result (stddev)
    type(property), intent(in) :: p
    real(kind=quad) :: avg, stddev
    avg = get_property_average (p)
    stddev = sqrt((p%sum2 / p%count) - (avg ** 2))
end function get_property_stddev


! HELPER METHODS

! initilizes property values, count, sums to zero
subroutine reset_property (p)
    type(property), intent(inout) :: p 
    p%val = 0.
    p%sum = 0.
    p%sum2 = 0.
    p%count = 0
end subroutine reset_property

! adds value to property sums and increments count
subroutine accumulate_property (p, val)
    type(property), intent(inout) :: p ! property to accumulate
    real(kind=quad), intent(in) :: val ! defines opertation to preform on proepry
    p%val = val
    p%sum = p%sum + p%val
    p%sum2 = p%sum2 + (p%val ** 2)
    p%count = p%count + 1
end subroutine accumulate_property

end module class_property