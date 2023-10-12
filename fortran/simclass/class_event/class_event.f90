! AUTHOR : Matthew Dorsey
! DATE : 2022-07-21
! FILENAME : class_event.f90
! PURPOSE : the event class is used to describe events that occur in DMD simulations.
!           each event has a tiem when it will occur, a type, and a partner.    

module class_event
use dmdconstants
implicit none
real(kind=quad), parameter :: bigtime = 1e10 ! unreasonably large time for events that will not occur

! define event class
type :: event
     private
     real(kind=quad) :: time ! time until event occurs (in seconds)
     integer :: type ! number representing the type of event that will occur
     integer :: up, dn ! uplist and downlist partneres associated with event
end type event

contains

! CONSTRUCTORS

! constructs event object
function construct_event (time, type, up, dn) result (e)
    real(kind=quad), intent(in) :: time 
    integer, intent(in) :: type, up, dn
    type(event) :: e 

    ! assign parameters to new event
    e%time = time 
    e%type = type
    e%up = up
    e%dn = dn
end function construct_event

! returns an event with big time, and no partner or type 
function null_event () result (e)
    type(event) :: e 

    ! create null event 
    e%time = bigtime 
    e%type = 0 
    e%up = 0
    e%dn = 0
end function null_event

! GETTERS AND SETTERS 

! returns time of event passed to function 
function get_event_time (e) result (t)
    type(event), intent(in) :: e 
    real(kind=quad) :: t
    t = e%time
end function get_event_time

! return type of event passed to partner
function get_event_type (e) result (t) 
    type(event), intent(in) :: e 
    integer :: t 
    t = e%type
end function get_event_type

! return uplist partner of event passed to function
function get_event_uppartner (e) result (up)
    type(event), intent(in) :: e
    integer :: up
    up = e%up
end function get_event_uppartner

! return dnlist partner of event passed to function
function get_event_dnpartner (e) result (dn)
    type(event), intent(in) :: e 
    integer :: dn
    dn = e%dn 
end function get_event_dnpartner

! I / O

! reports status of event to CLT 
subroutine report_event_status (e)
    type(event), intent(in) :: e

    if (e%up == 0) then
        write (*,'("No event scheduled.")')
    else
        write(*,'("Particle ", I4," will collid with particle ", I4," in ", F5.4," seconds. (Type ",I2,").")') &
            e%up, e%dn, e%time, e%type
    endif
end subroutine report_event_status

end module class_event