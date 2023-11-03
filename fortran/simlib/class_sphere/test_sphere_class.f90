! AUTHOR : Matthew Dorsey
! DATE : 2022-07-21
! FILENAME : test_sphere_class.f90
! PURPOSE : Program that tests the sphere class and functions

program test_sphere_class
use dmdconstants
use class_event
use class_sphere
implicit none

! simulation constants 
real(kind=quad), parameter :: region = 10.
real(kind=quad), parameter :: temp = 1.
integer, parameter :: mols = 100
real(kind=quad), parameter :: sigma1 = 1.

! local variables
type(sphere), dimension(mols) :: s, s2
type(event), dimension(mols) :: schedule
type(event) :: nxtevent ! used for predicting events
integer :: i, j ! indexing parameters 

! initialization
if (.not. generate_random_positions(s, region)) then 
    write (*,*)"Unable to generate random positions with region and number of spheres provided."
    call exit()
endif
call generate_random_velocities(s, temp)

! predict hard sphere collisions
call schedule_spheres (s, schedule, region)

! test that sphere array saving and loading works
s2 = s 
if(.not. save_spheres(s, 'sphereSAVE.dat')) call exit()
if(.not. load_spheres(s, 'sphereSAVE.dat')) call exit()

! report to user
write(*,*)
do i = 1, mols
    call report_sphere_status(s(i))
    call report_event_status(schedule(i))
    if (.not. compare_spheres(s(i), s2(i))) then
        write (*,'("Data transfer for sphere number ", I4," in array unsuccessful.")') i
        call exit() 
    endif 
    write(*,*)
enddo

end program test_sphere_class