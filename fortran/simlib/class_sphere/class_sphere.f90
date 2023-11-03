! AUTHOR : Matthew Dorsey
! DATE : 2022-07-21
! FILENAME : class_sphere.f90
! PURPOSE : sphere class represents and isotropic object in DMD simulations.
!           every sphere is defined by position and velocity vectors.
!           the dimension of vectors, typically 2 or 3, is set in the dmdconstants module.


module class_sphere
use dmdconstants
use class_event
implicit none

! class constants
integer, parameter :: ndim = dim3 ! number of dimensions that define the sphere class
!character(len=20) :: '_posSAVE.dat' ! tag used when saving and loading spheres
real(kind=quad), parameter :: hard_sphere_diameter = 1. ! default diameter of all hard spheres, unleses otherwise specified
real(kind=quad), parameter :: hard_sphere_excluded_volume = (4. / 3.) * pi * (hard_sphere_diameter / 2.)

! define sphere class
type :: sphere
     integer :: id
     real(kind=quad) :: diameter = hard_sphere_diameter ! diameter of sphere
     real(kind=quad), dimension(ndim) :: position, velocity ! position and velocity vectors for modeling dynamics
end type sphere

contains


! CONSTRUCTORS

! generates a random distribution of positions such that no spheres are overlapping
function generate_random_positions (s, region) result (success)
    implicit none
    type(sphere), dimension(:), intent(inout) :: s
    integer :: l ! length of sphere array
    real(kind=quad) :: region ! length of simulation box
    integer, parameter :: attemplimit = 1e6 ! maximum number of attemps at generating random positions
    logical :: success ! boolean that defines whether positions were randomly generated
    integer :: i, j, count ! used for indexing

    ! initailize incrementing variables
    l = size(s)
    count = 0
    i = 01
    do
        count = count + 1

        ! assign the sphere a random position
        call random_position (s(i)%position, region)

        ! check that the position does not overlap with any s that have already been assigned
        success = .true.
        if (i /= 1) then  
            do j = 1, i - 1
                if (is_overlapping(s(i), s(j), region)) success = .false.
            enddo
        endif

        ! if the sphere is not overlapping with any others
        if (success) then
            s(i)%id = i
            i = i + 1 ! increment the loop
            count = 0 ! reset the count
            if (i > l) return ! exit the subroutine
        else ! otherwise
            ! check the count
            if (count == attemplimit) return ! exit the subroutine
        endif
    enddo
end function generate_random_positions

! assigns sphere random position according to region length
subroutine random_position (pos, region)
    implicit none
    real(kind=quad), dimension(ndim), intent(inout) :: pos ! position of sphere
    real(kind=quad), intent(in) :: region ! maximum length of simulation box
    integer :: m ! used for indexing

    do m = 1, ndim
        call random_number (pos(m))
        pos(m) = pos(m) * region
    enddo
end subroutine random_position

! generates position distribution of sphere along a fcc lattice
!function generate_fcc_lattice (s, region) result (success)
!    implicit none 

!end function generate_fcc_lattice

! generates a MB distribution of velcoties, whose linear momentum along each axis is zero
subroutine generate_random_velocities (s, temp)
    implicit none
    type(sphere), dimension(:), intent(inout) :: s ! array of spheres
    integer :: l ! length of sphere array
    real(kind=quad), intent(in) :: temp ! temperature to assign velocity distribution along
    real(kind=quad), dimension(ndim) :: lm ! linear momentum
    integer :: i, m ! used for indexing 

    ! asssign velocities
    l = size(s)
    do i = 1, l 
        call random_velocity (s(i)%velocity, temp)
    enddo

    ! calculate the linear momentum contributions and reduce
    do m = 1, ndim 
        lm(m) = sum(s(:)%velocity(m)) / real(l)
        ! adjust the average linear momentum contrubtions to zero
        s(:)%velocity(m) = s(:)%velocity(m) - lm(m)
    enddo
end subroutine generate_random_velocities

! assigns sphere random velocity accorinding to a MB distribution
subroutine random_velocity (vel, temp)
    implicit none
    real(kind=quad), dimension(ndim), intent(out) :: vel ! velocity vector of sphere
    real(kind=quad), intent(in) :: temp ! system temperature
    real(kind=quad) :: mu, sigma, u1, u2 ! random numbers used for distribution generation
    integer :: m ! used for indexing

    mu = 0. 
    sigma = sqrt(temp)

    do m = 1, ndim
        call random_number (u1)
        call random_number (u2)
        vel(m) = mu + sigma * sqrt( -2. * log(u1)) * sin (2 * pi * u2)
    enddo
end subroutine random_velocity


! I / O

! prints the current status of a sphere to the CLT
subroutine report_sphere_status (i)
    type(sphere), intent(in) :: i
    write(*,'("SPHERE: ", I4," (DIAMETER: ", F4.2,")")') i%id, i%diameter
    write(*, '("POSITION: ", 3(F6.2))') i%position(:)
    write(*, '("VELOCITY: ", 3(F6.2))') i%velocity(:)
end subroutine report_sphere_status

! prints the status of the sphere in xyz file format 

! save an array of spheres to a file 
function save_spheres (s, filename) result (iobool)
    implicit none 
    type(sphere), dimension(:), intent(in) :: s ! array of spheres 
    character(len=*), intent(in) :: filename ! location and file name to store information
    logical :: iobool ! determines if data transfer was successful
    integer :: ioerr ! used to determine success for information transfer
    integer :: i, m ! used for indexing

    iobool = .false.
    open (unit = saveiounit, file = trim(filename), status = 'REPLACE', action = 'WRITE', iostat = ioerr)
    if (ioerr == 0) then ! write information to file
        do i = 1, size(s)
            write(saveiounit, *) s%id 
            write(saveiounit, *) s%diameter  
            do m = 1, ndim 
                write(saveiounit,*) s%position(m)
            enddo
            do m = 1, ndim
                write(saveiounit,*) s%velocity(m)
            enddo
        enddo
        close (unit = saveiounit, status = 'KEEP')
        iobool = .true.
    endif
    ! if the operation was unsuccessful inform the used 
    if (.not. iobool) write(*,*) 'Unable to save sphere array to ', trim(filename),'.'
end function save_spheres

! load list of spheres from save file
function load_spheres (s, filename) result (iobool)
    implicit none 
    type(sphere), dimension(:), intent(out) :: s ! array of spheres
    character(len=*), intent(in) :: filename
    logical :: iobool ! used to denote the success of the load operation
    integer :: ioerr ! used to track success of data transfer 
    integer :: i, m ! used for data transfer

    iobool = .false.
    open (unit = saveiounit, file = trim(filename), status = 'OLD', action = 'READ', iostat = ioerr)
    if (ioerr == 0) then ! read information from the file
        do i = 1, size(s)
            read(saveiounit, *) s%id 
            read(saveiounit, *) s%diameter 
            do m = 1, ndim 
                read(saveiounit, *) s%position(m)
            enddo
            do m = 1, ndim 
                read(saveiounit, *) s%velocity(m)
            enddo
        enddo
        close (unit = saveiounit, status = 'KEEP')
        iobool = .true.
    endif
    ! if the operation was unsuccessful, inform the user
    if (.not. iobool) write(*,*) 'Unable to load sphere array from ', trim(filename),'.'
end function load_spheres

! HELPER METHODS

! logical function that determines if two sphere objects contain the same information
function compare_spheres (i, j) result (samebool)
    type(sphere), intent(in) :: i, j 
    logical :: samebool ! boolean that determines whether two spheres are the same
    integer :: m ! used for indexing
    samebool = (i%id == j%id) .and. (i%diameter == j%diameter)
    do m = 1, ndim 
        samebool = samebool .and. (abs(i%position(m) - j%position(m)) < tol)
        samebool = samebool .and. (abs(i%velocity(m) - j%velocity(m)) < tol)
    enddo
end function compare_spheres

! logical function that determines if two spheres are overlapping based on their diameters
function is_overlapping(i, j, region) result (bool)
    implicit none 
    type(sphere), intent(in) :: i, j ! uplist and downlist spheres
    real(kind=quad), intent(in) :: region ! length of simulation box
    real(kind=quad) :: dist ! minimum distance between two sphere i and j according to their diameter
    logical :: bool
    dist = (i%diameter + j%diameter) / 2.
    bool = (calculate_distance(i, j, region) < dist)
end function is_overlapping

! calculates the distance between two spheres
function calculate_distance (i, j, region) result (d)
    type(sphere), intent(in) :: i, j ! two spheres to calculate distance between
    real(kind=quad), intent(in) :: region ! length of simulation box, used to apply MIC
    real(kind=quad) :: d ! distance between two spheres 
    real(kind=quad), dimension(ndim) :: rij ! position vector pointing from i to j
    integer :: m ! indexing paramter
    rij = calculate_distance_vector(i, j, region)
    d = sqrt(dot_product(rij, rij))
end function calculate_distance

! determines the vector pointing from sphere i to sphere j
function calculate_distance_vector (i, j, region) result (rij) 
    type(sphere), intent(in) :: i, j ! two spheres to calculate relative position vector between
    real(kind=quad), intent(in) :: region ! length of simulation box, used to apply MIC
    real(kind=quad), dimension(ndim) :: rij ! position vector pointing from i to j 
    rij = i%position - j%position
    call apply_minimum_image_convention (rij, region)
end function calculate_distance_vector

! applies minimum image convention to position vector pointing between two spheres
subroutine apply_minimum_image_convention (rij, region)
    real(kind=quad), dimension(ndim), intent(inout) :: rij ! position vector pointing from i to j 
    real(kind=quad), intent(in) :: region ! length of simulation box 
    integer :: m ! used for indexing 
    do m = 1, ndim
        if (rij(m) >= (0.5 * region)) rij(m) = rij(m) - region 
        if (rij(m) < (-0.5 * region)) rij(m) = rij(m) + region
    enddo
end subroutine apply_minimum_image_convention

! TODO create interaction object (?) and move functions like this to that object type (hard, well, bond)

! create an event calander for a list of spheres
subroutine schedule_spheres (s, e, region)
    implicit none 
    type(sphere), dimension(:), intent(inout) :: s ! list of spheres
    type(event), dimension(:), intent(out) :: e ! event calander
    real(kind=quad), intent(in) :: region ! length of simulation box
    integer :: mols ! length of sphere and event arrays
    integer :: i, j ! used for indexing
    type(event) :: nxtevnt

    ! predict hard sphere collisions for all spheres in list
    mols = size(s)
    e = null_event()
    do i = 1, mols - 1
        do j = i + 1, mols 
            nxtevnt = predict_hard_collision (s(i), s(j), region)
            if (get_event_time(nxtevnt) < get_event_time(e(i))) then 
                e(i) = nxtevnt
            endif
        enddo
    enddo
end subroutine schedule_spheres

! calculates the time until two spheres will occur
! if collision will not occur, returns event with large time
function predict_hard_collision (i, j, region) result (e) 
    type(sphere), intent(in) :: i, j ! two spheres to compare, i is downlist of j and vice versa
    real(kind=quad), intent(in) :: region ! length of simulation box 
    type(event) :: e ! event predicted by algorithm
    real(kind=quad) :: colldist ! collision distance between two spheres
    real(kind=quad), dimension(ndim) :: rij, vij ! position and velocity vectors between i and j
    real(kind=quad) :: aij, bij, cij, discr

    ! set the events uplist partner
    e = null_event()
    colldist = (i%diameter + j%diameter) / 2.

    ! calculate position and velocity vectors between i and j
    rij = calculate_distance_vector (i, j, region)
    vij = i%velocity - j%velocity

    ! calculate quadratic parameters 
    aij = dot_product (vij, vij)
    bij = dot_product (vij, rij)
    cij = dot_product (rij, rij) - colldist ** 2 
    discr = (bij ** 2) - (aij * cij)

    ! predict if the spheres will collide
    if ((discr > 0.) .and. (bij < 0.)) then 
        ! an event will occur at the distance passed to the method 
        e = construct_event((-bij - sqrt(discr)) / aij, 1, i%id, j%id)
    endif ! otherwise no event will occur
end function predict_hard_collision


! GETTERS AND SETTERS

! returns the position vector of a sphere
function get_sphere_position (i) result (pos)
    type(sphere), intent(in) :: i
    real(quad), dimension(ndim) :: pos 
    pos = i%position
end function get_sphere_position

! returns the velocity of a sphere
function get_sphere_velocity (i) result (vel)
    type(sphere), intent(in) :: i ! sphere 
    real(quad), dimension(ndim) :: vel 
    vel = i%velocity
end function get_sphere_velocity
    
end module class_sphere