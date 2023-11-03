!*************************************************************
!** Programmer: Matthew A Dorsey 
!**				Dr. Carol Hall Research Laboratory
!**				Chemical and Biomolecular Engineering 
!**				North Carolina State University 
!**
!** 
!** Date: 11/18/2022
!** 
!**	Purpose: Uses DMD algorithm to simulate hard, active discs.
!**
!*************************************************************

module active_disk_module
use dmdconstant
use binarytree
use type_id

!*************************************************************
!** GLOBAL CONSTANTS
!*************************************************************

! ** ensemble ************************************************
integer, parameter :: ndim = 2 ! number of dimensions 
integer, parameter :: cell = 15 ! number of lattice cells
integer, parameter :: mer = 1 ! number of hardspheres which make one disc grouping
real(kind=dbl), parameter :: xa = 1.0 ! number fraction of active particles
real(kind=dbl), parameter :: eta = 0.20 ! packing fraction: ration of total sphere area to total area [RW MAX: ~ 0.4..]
real(kind=dbl), parameter :: tempstart = 3.0 ! starting temperature of the simulation 
real(kind=dbl), parameter :: tempfinal = 0.01 ! final system temperature set point 
real(kind=dbl), parameter :: sigma1 = 1.0 ! distance from sphere center to first discontinuity (diameter of hardsphere)
real(kind=dbl), parameter :: sigma2 = 0.6 * sigma1 ! distance from sphere center to second discontinuity
real(kind=dbl), parameter :: sigma3 = 0.85 * sigma1 ! distance from sphere center to third discontinuity
real(kind=dbl), parameter :: sigma4 = 1.4 * sigma1 ! distance from sphere center to fourth discontinuity
real(kind=dbl), parameter :: epsilon2 = 2.0 ! depth of first well, should be greater than epsilon 3
real(kind=dbl), parameter :: epsilon3 = 0.5 ! depth of second well 
real(kind=dbl), parameter :: epsilon4 = 0.25 ! depth of third well
real(kind=dbl), parameter :: delta = 0.03 ! half the bond length
!real(kind=dbl), parameter :: exittemp = 0.01 ! final system temperature 
!*************************************************************
integer, parameter :: disk = cell ** ndim ! number of discs
integer, parameter :: mols = disk * mer ! number of hard spheres
integer, parameter :: na = disk * xa ! number of active discs
real(kind=dbl), parameter :: excluded_area = pi * ((sigma1 / 2.) ** 2) ! area occupied by one hard disk
real(kind=dbl), parameter :: area = (excluded_area * real(disk)) / eta ! area of simulation box
real(kind=dbl), parameter :: region = sqrt (area) ! length of simulation box wall
real(kind=dbl), parameter :: density = real(disk) / area ! number density of dumbbells
real(kind=dbl), parameter :: sg1sq = sigma1 ** 2
real(kind=dbl), parameter :: sg2sq = sigma2 ** 2
real(kind=dbl), parameter :: sg3sq = sigma3 ** 2
real(kind=dbl), parameter :: sg4sq = sigma4 ** 2
real(kind=dbl), parameter :: inbond = (sigma1 / 4.) * (1. - delta) ! inner bond distance
real(kind=dbl), parameter :: onbond = (sigma1 / 4.) * (1. + delta) ! outer bond distance
real(kind=dbl), parameter :: icbond = (sigma1 / 2.) * (1. - delta) ! inner across bond distance
real(kind=dbl), parameter :: ocbond = (sigma1 / 2.) * (1. + delta) ! outer across bond distance


! ** ANDERSON ************************************************
! ** anderson thermostat *************************************
logical, parameter :: thermostat = .true. ! anderson thermostat status: .false. == off, .true. == on
real(kind=dbl), parameter :: thermal_conductivity = 50.0 ! the thermal conductivity of the hardsphere system [THIS VALUE HAS NOT BEEN VERIFIED!!]
real(kind=dbl), parameter :: thermo_period = (density ** (1./2.)) / thermal_conductivity ! average number of collisions occuring per second
! ** anderson field ******************************************
logical, parameter :: field = .false. ! status of magnetic field: .false. == off, .true. == on
real(kind=dbl), parameter :: magforce = sqrt(2. * tempstart) ! magnitude of the force that is exerted by the magnetic field
real(kind=dbl), parameter :: mag_period = (1. / 600.) ! amount of time between collision events, in seconds
! ** field rotation **
logical, parameter :: field_rotation = .false. ! determines whether magnetic field rotates, field must be on
real(kind=dbl), parameter :: field_angvel = 50.0 ! rate at which the angle of the field changes in rads / sec
! ** field actuation **
logical, parameter :: field_actuation = .false. ! determines whether the magnetic field switches on and off, field must be on
real(kind=dbl), parameter :: field_actfreq = 1. / 50.0 ! rate at which switchs from on to off (and vice versa) in herz
! ** anderson activity ****************************************
logical, parameter :: activity = .false.
real(kind=dbl), parameter :: active_force = sqrt(2. * tempstart)
real(kind=dbl) :: active_period = (1. / 600.)

! ** SETTINGS ************************************************
! ** events **************************************************
real(kind=dbl), parameter :: frac = 0.95 ! fractional amount by which to reduce the temperarture of the annealing simulation
integer, parameter :: total_events = 200000000 ! total number of events in the simulation at the temperature set point
integer, parameter :: event_equilibrium = 0.1 * total_events ! number of steps during which the system is considered at equilibrium, equilibrium property values are taken
integer, parameter :: event_equilibriate = 0.9 * total_events ! number of steps during which the system is allowed to equilibriate
integer, parameter :: event_reschedule = 300 ! number of steps by which to reschedule entire event calander, adjust temperatures
integer, parameter :: event_average = 1000000 ! number of steps between property calculations and complete rescheduling
integer, parameter :: propfreq = 1000000 ! frequency of property calculations, when the system is not at equilibrium
real(kind=dbl), parameter :: tol = 0.01 ! amount by which to allow mistakes from numberical integration (i.e. overlap and bcs)
integer, parameter :: debug = 0 ! debugging status: 0 == off, 1 == on, 2 == extra on 
! ** order parameters ****************************************
integer, parameter :: orderlength = 15
real(kind=dbl), parameter :: orderwidth = sigma2
! ** file management *****************************************
character(len=15), parameter :: simtitle = 'actdisk'
character(len=10), parameter :: simid = 'TEST'
character(len=40), parameter :: group_save_file = trim(simtitle) // trim(simid) // '_groupSAVE.dat' ! save file containing all false position vectors 
character(len=40), parameter :: state_save_file = trim(simtitle) // trim(simid) // '_dmdSAVE.dat' ! save file containing all simulation state 
integer, parameter :: saveiounit = 11 
integer, parameter :: simiounit = 12
integer, parameter :: sphcooriounit = 13
integer, parameter :: objcooriounit = 14
integer, parameter :: reportiounit = 15
integer, parameter :: annealiounit = 16
integer, parameter :: opiounit = 17
! ** animation settings **************************************
logical, parameter :: movie = .true. ! movie making status: 0 == off, 1 == on
real(kind=dbl), parameter :: movfreq = 0.1 ! frequency to take snapshots of movies [reduced seconds]

! ** EFFICIENCY TECHNIQUES ***********************************
! ** cell + neighbor list ************************************
integer, parameter :: nbrListSize = 40 ! average number of neighbors per particle 
integer, parameter :: nbrListSizeMax = 100 ! maximum number of particles accessible to list 
real(kind=dbl), parameter :: nbrRadiusMinInt = 5.0 ! integer used for determing the max displacement required for a neighborlist update
real(kind=dbl), parameter :: nbrRadiusMin = (nbrRadiusMinInt / (nbrRadiusMinInt - 1)) * (sigma4 - (sigma1 / nbrRadiusMinInt)) ! minimum required radius
real(kind=dbl), parameter :: nbrRadius = max(sqrt((real(nbrListSize) / (real(mer) * density)) * (1.0 / pi)), nbrRadiusMin) ! radius of neighborlist
real(kind=dbl), parameter :: nbrDispMax = (nbrRadius - sigma1) / nbrRadiusMinInt ! max particle displacement before a neighborlist update is required
integer, parameter :: nCells = floor (region / nbrRadius) ! number of cells in one dimension, cell length cannot be shorter than the nerighbor radius
real(kind=dbl), parameter :: lengthCell = region / real (nCells) ! legnth of each cell in one dimension, must be greater than the square well length (sig2)



!*************************************************************
!** GLOBAL TYPES: data constructions for OOP
!*************************************************************

type :: position 
    real(kind=dbl), dimension(ndim) :: r ! vector describing the position of one sphere
end type position 

type :: velocity
    real(kind=dbl), dimension(ndim) :: v ! vector describing the velocity of one sphere 
end type 

type :: event 
    real(kind=dbl) :: time 
    integer :: type
    type(id) :: dnpart, uppart
end type event 

type :: particle
    type(position) :: fpos ! false location of sphere 
    type(velocity) :: vel ! velocity of sphere
    integer :: pol ! polarization of sphere (0 = neutral, -1 = negative, 1 = positive)
    type(event) :: collision ! next collision event for sphere 
    type(id), dimension(nbrListSizeMax) :: upnab, dnnab ! uplist and downlist neightbors with nbrRadius distance used for event scheduling
    type(id), dimension(orderlength) :: orderlist ! list of oppositely charged pairs within orderwidth distance used for calculating order parameters 
end type particle

type :: group
    type(particle), dimension(mer) :: circle ! each group is made of mer circles 
    integer :: chai ! integer describing the chiraliry of each dumbbell group
    ! TODO add activity boolean
    ! TODO: add string describing each circle?
end type group

type :: dmd  
    integer :: n_cell ! specifices the number of lattice cells in each dimension of the simulation box
    real(kind=dbl) :: temp ! system temperature
    real(kind=dbl) :: den ! system density as an area fraction
    integer :: anneal ! number of times the simulation has reduced the temperature
    real(kind=dbl) :: timenow ! current length of simulation 
    real(kind=dbl) :: timeperiod ! current length of period 
    integer :: n_events, n_col, n_ghost, n_bond, n_hard, n_well, n_mag ! event counting
end type dmd

type :: property
    real(kind=dbl) :: value, sum, sum2, equilavg, equilstd
    integer :: count, equilibcount
end type property


!*************************************************************
!** GLOBAL VARIABLES
!*************************************************************
type(group), dimension(disk) :: d ! list of discs
type(dmd) :: sim ! state of dmd simulation
real(kind=dbl) :: tsl ! time since the system was last updated, encapsulated in list of groups
! ** event scheduling ****************************************
!type(event), dimension(mols + 4) :: schedule  ! TODO can these be pointers?
integer, parameter :: ghost = 5 ! number of ghost events
type(event), dimension(ghost) :: ghost_event
type(node), dimension(mols + ghost) :: eventTree ! binary tree list used for scheduling collision events
! type(event) :: ghost_event, mag_event, switch_event, movie_event ! periodic simulation events
! ** simulation properties ***********************************
type(property) :: te ! total energy 
type(property) :: ke ! kinetic energy 
type(property) :: pot ! potential energy 
type(property) :: temp ! temperature
type(property) :: nem ! nematic order parameter 
!type(property) :: perc ! percolation order parameter
! ** efficiency methods ** 
logical :: nbrnow
real(kind=dbl) :: dispTotal
integer :: rootnode ! pointer to first node of binary tree using for scheduling collision events


contains 

! ** SIMULATION **********************************************

! ** initialization ******************************************

! method called by main program, initialize the system 
! belongs to ... ? (simulation, module ... or object?)
subroutine initialize_system ()

    ! load, save state
    if (.not. load_state(sim, state_save_file, saveiounit)) then
        call set_dmd(sim) ! belongs to dmd object
        if (.not. save_state(sim, state_save_file, saveiounit)) then
            call exit(1)
        endif
    endif

    ! load, save groups
    if (.not. load_group(d, group_save_file, saveiounit)) then 
        call set_group(d, region, sim%temp) ! belogs to group object 
    	if (.not. save_group(d, group_save_file, saveiounit)) then
            call exit(1)
        endif
    endif

    ! initialize properties
    call open_files()
    call initialize_binarytree(eventTree, rootnode)
    call build_neighborlist(d, region)
    call update_properties(0)
end subroutine initialize_system

! method saves the state of the simulation
! called by simulation
! stored with simulation 
subroutine save() 
    implicit none 

    ! simulation
    if (.not. save_state(sim, state_save_file, saveiounit)) &
        write (*, '(" save: unable to save the simulation state.")')

    ! save groups 
    if (.not. save_group(d, group_save_file, saveiounit)) &
        write (*, '(" save: unable to save group state.")')
end subroutine save

! ** dynamics ************************************************

! method that processes a single event in DMD simulation
function single_step() result (stop)
    implicit none 
    logical :: stop
    logical :: one_step
    type(id) :: next

    stop = .false.
    ! determine if the neighbor list should be updates
    if (nbrnow) call build_neighborlist(d, region)

    ! find the next event
    next = find_next()

    ! move the system "forward in time" 
    call forward(next)

    ! process the next event
    if (next%one <= disk) then 

        ! if the event is a collision 
        call collide (next)

    elseif (next%one > disk) then 
        if ((next%one - disk) == 1) then 
            ! ghost thermostat collision
            call thermal_ghost_event (d, sim%temp, region, sim%n_ghost)

            ! reschedule ghost thermostat event
            ghost_event(1) = predict_ghost(thermo_period)
            call addbranch(eventTree, mols+1)

        elseif ((next%one - disk) == 2) then 
            ! ghost magnetic collision 
            call field_ghost_event (d, sim%temp, region, sim%n_mag)

            ! reschedule the ghost magnetic event
            ghost_event(2) = predict_ghost(mag_period)
            call addbranch(eventTree, mols+2)

        elseif ((next%one - disk) == 4) then 
            ! record position
            call record_position_circles(d, sphcooriounit, region)
            call record_position_objects(d, objcooriounit, region)

            ! reschedule ghost event
            ghost_event(4) = predict_ghost(movfreq)
            call addbranch(eventTree, mols + 4)
        elseif ((next%one - disk) == 5) then 
            ! ghost activity collision
            call activity_ghost_event (d, sim%temp, region, sim%n_ghost)

            ! reschedule acitivty ghost event 
            ghost_event(5) = predict_ghost(active_period)
            call addbranch(eventTree, mols+5)

        else
            ! if the event is a ghost 
            write (*,'(" single_step: ghost event ", I2, " called.")') (next%one - disk)
            call exit()
        endif
    endif
    ! event based decisions 
    ! stop simulation

    ! 
    if (mod(sim%n_events, event_reschedule) == 1) then 
        call complete_reschedule()
        write (*,'(" single_step: complete reschedule called after ", I10, " events at time ", F4.1)') sim%n_events, sim%timenow
    endif

    if (sim%n_events >= 100000) stop = .true.
end function single_step

! function that finds the next event
! using logic gates can switch between using an event tree
! and serial searching
! called by simulation
! stored with simulation 
function find_next() result (next)
    implicit none 
    type(id) :: next 
    integer :: i 

    ! TODO: if binary tree method is on
    i = findnextevent (eventTree, rootnode)
    if (i <= mols) then 
        next = mol2id(i, mer)
    elseif (i > mols) then 
        next%one = disk + (i - mols)
        next%two = 0
    endif

    ! else find the next event using serial search method
end function find_next

! method moves the system forward in time to the next event
! called by simulation
! stored with simulation
subroutine forward (next) 
    implicit none 
    type(id), intent(in) :: next 
    type(id) :: a 
    integer :: i 
    real(kind=dbl) :: tnext 
    type(event) :: enext

    sim%n_events = sim%n_events + 1

    ! check time
    if (next%one <= disk) then 
        enext = d(next%one)%circle(next%two)%collision
    elseif(next%one > disk) then 
        enext = ghost_event(next%one - disk)
    end if 
    tnext = enext%time

    if (tnext < 0) then
        write (*,'(" forward: negative time calculated at ", I9)') sim%n_events
        write (*,*) ''
        write (*,*) '** N E G A T I V E   T I M E   C A L C U L A T I O N **'
        write (*,*) '' 
        call report_schedule()
        !call report_neighborlist(d)
        call exit (1)
    endif

    ! report the event to the user
    if (debug >= 1) then 
        call report_event (enext)
    endif

    ! TODO: if the false position method is off
    ! integrate the position of the particles forward in time
    ! according to their velocity

    ! remove the time til the next event from event calander
    sim%timenow = sim%timenow + tnext 
    sim%timeperiod = sim%timeperiod + tnext
    tsl = tsl + tnext
    do i = 1, mols
        a = mol2id(i, mer)
        d(a%one)%circle(a%two)%collision%time = &
            d(a%one)%circle(a%two)%collision%time - tnext
    enddo
    if (thermostat) ghost_event(1)%time = ghost_event(1)%time - tnext
    if (field) ghost_event(2)%time = ghost_event(2)%time - tnext 
    if (field_actuation) ghost_event(3)%time = ghost_event(3)%time - tnext 
    if (movie) ghost_event(4)%time = ghost_event(4)%time - tnext
    if (activity) ghost_event(5)%time = ghost_event(5)%time - tnext

    ! check boundaries 
    if ((debug >= 1) .and. check_boundaries()) then 
        write (*, '(" forward: boundaries overlapping at event", I9)') sim%n_events
        write (*,*) ''
        write (*,*) '** B O U N D A R I E S   O V E R L A P **'
        write (*,*) '' 
        call report_schedule()
        call exit(1)
    endif
end subroutine forward

! method that determines how the particles veloccities should
! change when a collision event occurs 
! called by simulation
! stored with group_list
subroutine collide (a)
    implicit none
    type(id), intent(in) :: a ! downlist particle
    type(id) :: b ! uplist particle 
    type(event) :: e ! event for participating particles
    real(kind=dbl), parameter :: smdistance = 1e-11
    logical, parameter :: debug_collide = (debug >= 1) .and. (.true.)
    type(position) :: rij
    type(velocity) :: vij
    real(kind=dbl), dimension(ndim) :: impulse_a, impulse_b ! the transfer of energy between the two atoms 
    type(velocity) :: vel_a, vel_b
    real(kind=dbl) :: bij, distance, bump, discr2, discr3, discr4, delep ! dot prroduct of moment and velocity vectors, distance between atoms
    real(kind=dbl) :: dispa, dispb ! displacement of a and b particles
    integer :: q ! indexing parameter 
    real(kind=dbl) :: scale_a, scale_b ! values to scale impulse by to keep the velocity constant

    ! incriment the number of collisions by 1
    sim%n_col = sim%n_col + 1

    !initialize parameters 
    b = d(a%one)%circle(a%two)%collision%uppart
    e = d(a%one)%circle(a%two)%collision
    distance = 0.
    bij = 0.
    dispa = 0.
    dispb = 0.

    ! determine the real position of the particle, calulate distance
    do q = 1, ndim 
        ! calculate the false position of the particle pair
        rij%r(q) = d(a%one)%circle(a%two)%fpos%r(q) - d(b%one)%circle(b%two)%fpos%r(q)
        vij%v(q) = d(a%one)%circle(a%two)%vel%v(q) - d(b%one)%circle(b%two)%vel%v(q)
        ! calculate the real position of the particle pair
        rij%r(q) = rij%r(q) + vij%v(q) * tsl
        ! apply minimum image convention
        if (rij%r(q) >= 0.5*region) rij%r(q) = rij%r(q) - region 
        if (rij%r(q) < -0.5*region) rij%r(q) = rij%r(q) + region 
        distance = distance + (rij%r(q) ** 2)
        bij = bij + (rij%r(q) * vij%v(q))
        ! determine the displacement of either particle based on their pre-collision velocities
        dispa = dispa + (d(a%one)%circle(a%two)%vel%v(q) * tsl) ** 2
        dispb = dispb + (d(b%one)%circle(b%two)%vel%v(q) * tsl) ** 2
    end do 
    distance = sqrt(distance)

    ! TODO: do not need to do this if the neighbor list method is not being used
    ! determine if the pre-collision displacement of either particle is greater than
    ! the maximum displacement required for a neighborlist update
    if (dispa >= dispb) then
        dispa = sqrt(dispa)
        if (dispa >= nbrDispMax) nbrnow = .true.
    else
        dispb = sqrt(dispb)
        if (dispb >= nbrDispMax) nbrnow = .true.
    endif

    ! notify user about location of the particles at the time of the event 
    if (debug_collide) then 
        write(*, 120) a%two, a%one, b%two, b%one, &
            distance, d(a%one)%circle(a%two)%collision%type
        120 format (' collision: The distance between ', I3, ' of ',I5, &
            ' and ', I3, ' of ', I5,' is ', F10.8,' (Type ',I2,').')
    end if 

    ! calculate collision dynamics based on event type 
    if (e%type == 1) then ! the event is at the repulsive core
        do q = 1, ndim 
            impulse_a(q) = -rij%r(q) * (bij / sg1sq)
        end do 
        bump = 0.0
        sim%n_hard = sim%n_hard + 1
    else 
        ! there are no other types of collisions
        write (*,*) 'collide: non-hard sphere interaction called. Abort.'
        call exit()
    end if 

    ! debugging
    if (debug >= 1) then
        write (*,*) ''
        write (*,*) 'pre scaling'
        !write (*,*) ''
        call report_circles (d(a%one)%circle)
        call report_circles (d(b%one)%circle)
    endif

    ! store values for calculation
    impulse_b = -impulse_a
    vel_a = d(a%one)%circle(a%two)%vel
    vel_b = d(b%one)%circle(b%two)%vel

    ! calculate scale values
    scale_a = scale_velocity (d(a%one)%circle(a%two)%vel, impulse_a)
    scale_b = scale_velocity (d(b%one)%circle(b%two)%vel, impulse_b)
    !scale_a = 1.
    scale_b = 1.

    ! adjust the velocity vector of each atom accordingly
    do q = 1, ndim 
        ! adjust the velocity of both colliding atoms 
        d(a%one)%circle(a%two)%vel%v(q) = (vel_a%v(q) + impulse_a(q)) * scale_a
        d(b%one)%circle(b%two)%vel%v(q) = (vel_b%v(q) + impulse_b(q)) * scale_b
        ! move the particles a small distance away from the collision point
        d(a%one)%circle(a%two)%fpos%r(q) = d(a%one)%circle(a%two)%fpos%r(q) + &
            bump * smdistance * rij%r(q) * sigma2
        d(b%one)%circle(b%two)%fpos%r(q) = d(b%one)%circle(b%two)%fpos%r(q) - &
            bump * smdistance * rij%r(q) * sigma2
        ! calculate the new false position of each particle based on their new velocities
        d(a%one)%circle(a%two)%fpos%r(q) = d(a%one)%circle(a%two)%fpos%r(q) + tsl * &
            (vel_a%v(q) - d(a%one)%circle(a%two)%vel%v(q))
        d(b%one)%circle(b%two)%fpos%r(q) = d(b%one)%circle(b%two)%fpos%r(q) + tsl * &
            (vel_b%v(q) - d(b%one)%circle(b%two)%vel%v(q))
    end do  

    ! apply periodic boundary conditions
    call periodic_boundary_conditions(d(a%one)%circle(a%two), region)
    call periodic_boundary_conditions(d(b%one)%circle(b%two), region)

    ! calculate the dispalcement of either particle based on their post-collision, scaled velocities
    dispa = 0.
    dispb = 0.
    do q = 1, ndim 
        dispa = dispa + (d(a%one)%circle(a%two)%vel%v(q) * tsl) ** 2
        dispb = dispb + (d(b%one)%circle(b%two)%vel%v(q) * tsl) ** 2
    enddo

    ! TODO: only need to do this when false position method is .on.
    ! determine if the post-collision displacement of either particle is greater than
    ! the maximum displacement required for a neighborlist update
    if (dispa >= dispb) then
        dispa = sqrt(dispa)
        if (dispa >= nbrDispMax) nbrnow = .true.
    else
        dispb = sqrt(dispb)
        if (dispb >= nbrDispMax) nbrnow = .true.
    endif

    ! reschedule events for the particles that participated in the event
    call collision_reschedule(a, b) 

    ! debugging
    if (debug >= 1) then 
        write (*,*) ''
        write (*,*) 'post scaling, reschedule'
        !write (*,*) ''
        call report_circles (d(a%one)%circle)
        call report_circles (d(b%one)%circle)
    endif
end subroutine collide

! method that scales the velocity a particle to a set value
! called by group list
! stored with DMD methods
function scale_velocity (vel, i) result (scale)
    implicit none 
    type(velocity), intent(in) :: vel ! velocity vectory
    real(kind=dbl), dimension(ndim), intent(in) :: i ! impulse vector
    real(kind=dbl) :: scale ! amount to scale impulse vector by
    logical, parameter :: debug_scale = (debug >= 1)
    real(kind=dbl) :: mag, new_mag ! magnitude to scale particle's velocity to
    type(velocity) :: new_velocity
    integer :: q 

    ! target magnitude of the new vector
    mag = sqrt(2. * tempstart)

    ! calculate the new velocity
    do q = 1, ndim 
        new_velocity%v(q) = vel%v(q) + i(q)
    enddo

    ! determine how much to scale the new velocity by
    new_mag = 0.
    do q = 1, ndim 
        new_mag = new_mag + new_velocity%v(q) ** 2
    enddo
    new_mag = sqrt(new_mag)

    scale = mag / new_mag
end function scale_velocity

! completes uplist and down list rescheduling for
! both particles that participate in a collision 
! event
! called by group_list object / module
! stored with group_list object / module 
subroutine collision_reschedule(a, b)
    implicit none
    ! ** calling variables ***********************************
    type(id), intent(in) :: a, b
    ! ** local variables *************************************
    type(id) :: j1, j2 ! event partners 
    integer :: i, m, n  ! indexing parameters

    ! loop through all particles 
    do i = 1, disk
        do m = 1, mer 
            ! if any particles just participated in an event or were scheduled to participate
            ! in an event with the event particles, reschedule their uplist partner 
            if (((i == a%one) .and. (m == a%two)) .or. ((d(i)%circle(m)%collision%uppart%one == a%one) .and.&
                (d(i)%circle(m)%collision%uppart%two == a%two)) .or. ((i == b%one) .and. (m == b%two)) .or.&
                ((d(i)%circle(m)%collision%uppart%one == b%one) .and. &
                (d(i)%circle(m)%collision%uppart%two == b%two))) then 
                j1%one = i
                j1%two = m
                ! reset its event and find next event 
                d(i)%circle(m)%collision = reset_event()
                ! loop through their uplist neighbors to find next event partner
                n = 0
                uplist: do 
                    n = n + 1
                    j2 = d(i)%circle(m)%upnab(n)
                    if (j2%one == 0) exit 
                    if (idequiv(j1,a) .and. idequiv(j2, b)) then 
                        ! if j2 just participated in an even with j1
                        ! skip
                        cycle uplist
                    endif
                    call predict_collision (d(j1%one)%circle(j1%two), &
                        d(j2%one)%circle(j2%two), j1, j2, tsl, region)
                enddo uplist
                call addbranch (eventTree, id2mol(j1%one, j1%two, mer)) ! add the event to the event tree
            endif 
            ! if the particles are event partners 
            if (((i == a%one) .and. (m == a%two)) .or. ((i == b%one) .and. (m == b%two))) then 
                ! loop through downlist partners 
                j2%one = i 
                j2%two = m 
                n = 0
                downlist: do 
                    n = n + 1 
                    j1 = d(i)%circle(m)%dnnab(n)
                    if (j1%one == 0) exit 
                    if (idequiv(j1, a) .and. idequiv(j2,b)) then 
                        ! if j2 just participated in an event with j1
                        ! skip
                        cycle downlist
                    endif
                    call predict_collision (d(j1%one)%circle(j1%two), &
                        d(j2%one)%circle(j2%two), j1, j2, tsl, region)
                    if (idequiv(d(j1%one)%circle(j1%two)%collision%uppart,j2)) then ! if j1 is scheduled for an event with j2
                        call addbranch (eventTree, id2mol(j1%one, j1%two, mer))
                    endif
                enddo downlist
            endif 
        enddo
    enddo
end subroutine collision_reschedule

! completes uplist and downlist rescheduling for
! the one particle that participates in a ghost
! collision 
! called by group_list object / module
! stored with the group_list object / module
subroutine ghost_reschedule(ghost)
    implicit none
    integer, intent(in) :: ghost
    type(id) :: j1, j2 ! prediction event partners 
    integer :: i, m, n ! indexing parameters 

    ! loop through all particles 
    do i = 1, disk
        do m = 1, mer 
            ! if the partcle either participated in the ghost event, or was scheduled 
            ! to collide with a particle that participated in the ghost event 
            if ((i == ghost) .or. (d(i)%circle(m)%collision%uppart%one == ghost)) then 
                j1%one = i
                j1%two = m 
                ! reset and find its next event with an uplist particle
                d(j1%one)%circle(j1%two)%collision = reset_event()
                n = 0 
                uplist: do 
                    n = n + 1
                    j2 = d(j1%one)%circle(j1%two)%upnab(n)
                    if (j2%one == 0) exit 
                    call predict_collision (d(j1%one)%circle(j1%two), &
                        d(j2%one)%circle(j2%two), j1, j2, tsl, region)
                enddo uplist 
                call addbranch (eventTree, id2mol(j1%one, j1%two, mer)) ! add the event to the event tree
            endif
            ! if the particle participated in the ghost event 
            if (i == ghost) then 
                j2%one = i 
                j2%two = m 
                ! find any downlist particles in its collision path 
                n = 0
                downlist: do 
                    n = n + 1
                    j1 = d(j2%one)%circle(j2%two)%dnnab(n)
                    if (j1%one == 0) exit
                    call predict_collision (d(j1%one)%circle(j1%two), &
                        d(j2%one)%circle(j2%two), j1, j2, tsl, region)
                    if (idequiv(d(j1%one)%circle(j1%two)%collision%uppart,j2)) then ! if j1 is scheduled for an event with j2
                        call addbranch (eventTree, id2mol(j1%one, j1%two, mer))
                    endif
                enddo downlist 
            endif 
        enddo
    enddo
end subroutine ghost_reschedule

! ** GROUP INITIALIZATION ************************************ 

! method called by simulation module
! method that initializes the disks
! belongs to group object
subroutine set_group(d, region, temp)
    implicit none
    type(group), dimension(:), intent(out) :: d ! array of groups, that are discs 
    real(kind=dbl), intent(in) :: region !	 length of simulation box in each dimension
    real(kind=dbl), intent(in) :: temp ! value to initialize velocities

    write (*,*) 'set_group: initializing disks to initial state'
    call set_chirality(d)
    call set_polarity(d)
    call set_velocity(d, temp, region)
    call set_position(d, region)

    ! method owned by group (simulation?) that dictates how to apply the dmd algorithm 
    ! particles have interactions, which are stored in the list of disks
    ! the simulation should take those interactions, and incorperate ghost interactions
    ! as well as other periodic events (i.e have another list)
end subroutine set_group

! method that initializes the initial position of the disks
! initialization is either random, or ordered
! called by group object
! owned by group object (sphere object?)
subroutine set_position(d, region)
    implicit none
    type(group), dimension(:), intent(out) :: d ! array of groups, that are disks
    real(kind=dbl), intent(in) :: region !	 length of simulation box in each dimension
    integer, parameter :: max_attempts = 2
    integer :: attempts

    ! attept to randomly distribute the discs
    attempts = 0
    write (*,*) 'set_position: attempting to generate random configuration.'
    do 
	    if (random_walk(d, region)) then 
            write (*,*) 'set_position: random configuration generated'
            return 
	    else
            attempts = attempts + 1
	        if (attempts >= max_attempts) then
                write (*,*) 'set_position: unable to generate random configuration, ',&
                    'generating lattice configuration.'
	            if (lattice(d, region)) then 
                    write (*,*) 'set_position: lattice configuration generated.'
                    return
	            else
                    write (*,*) 'set_position: unable to initialize positions.'
                    call exit(1)
	            endif
	        endif
	    endif
    enddo
end subroutine set_position

! method that initializes the discs randomly
! if the method fails, false is returned as a value, else true
! called by group object
! owned by group object (sphere object)
function random_walk (d, region) result (success)
    implicit none
    type(group), dimension(:), intent(out) :: d ! array of groups, that are disks
    real(kind=dbl), intent(in) :: region !	 length of simulation box in each dimension
    logical, parameter :: debug_random_walk = (debug >= 1) .and. (.false.)
    integer, parameter :: max_attempts = 1000
    real(kind=dbl), parameter :: zero_tol = 0. 
    real(kind=dbl) :: dist ! distance between two spheres
    real(kind=dbl) :: theta ! orientation
    type(position) :: center ! central position of particle position (random)
    logical :: success ! whether the random walk algorithm generates positions for all discs
    logical :: overlap ! boolean used for determining particle overlap 
    integer :: i, m, j, n, q, count

    success = .false.
    tsl = 0.
    count = 0
    i = 1 ! for every disk
    do 
        if (debug_random_walk) write (*,"(' set_position: generating position for disk ', I4)") i
        ! give the disc a random position and random orientation
        call random_number(theta)
        theta = theta * twopi
        do q = 1, ndim
            call random_number(center%r(q))
            center%r(q) = center%r(q) * region
        enddo
        call generate_position(d(i), center, theta, region)

        ! determine if new position overlaps with other previous ones
        overlap = check_position(d(i), region)
        if (i > 1) then 
            do j = 1, i - 1
                overlap = check_overlap(d(i), d(j), region, zero_tol)
                if (overlap) exit
            enddo
        endif

        if (overlap) then 
            ! the group is overlapping with a different group
            ! try again 
            count = count + 1
            if (count > max_attempts) return 
        else
            ! the algorithm was successful
            if (i == size(d)) then 
                ! the algorithm was successful
                success = .true.
                return 
            else
                ! increment and do again
                i = i + 1
                count = 0
            endif
        endif
    enddo
end function random_walk

! function that generates the disc at a random position in the box
! called by group_list object / module
! contained by group_object / module 
subroutine generate_position(d, center, theta, region)
    implicit none
    type(group), intent(out) :: d
    real(kind=dbl), intent(in) :: region ! length of simulation box
    type(position), intent(in) :: center ! center of disc
    real(kind=dbl), intent(in) :: theta ! orientation of disk
    integer :: m, q 

    ! assign the disc the central position
    do m = 1, size(d%circle)
        do q = 1, ndim 
            d%circle(m)%fpos%r(q) = center%r(q)
        enddo
    enddo
end subroutine generate_position

! method that checks if the position of the disc is correct
! called by group_list
! store in group_object / list
function check_position (d, region) result (overlap)
    implicit none 
    type(group), intent(in) :: d 
    real(kind=dbl), intent(in) :: region ! length of simulation box
    logical :: overlap

    overlap = .false.
end function check_position

! function that calculates the distance between all 
! circles in an arrangement and determines if they are 
! overlapping. (NOTE: polarized spheres do not have hard 
! sphere interactions atm)
function check_overlap (di, dj, region, tol) result (overlap)
    implicit none
    type(group), intent(in) :: di, dj 
    real(kind=dbl), intent(in) :: region
    real(kind=dbl), intent(in) :: tol ! acceptance criteria
    logical :: overlap
    logical, parameter :: verbose = (debug >= 1)
    real(kind=dbl), parameter :: SIGMAI = 1.0 ! diameter of ith sphere
    real(kind=dbl), parameter :: SIGMAJ = 1.0 ! diameter of jth sphere
    real(kind=dbl) :: dist
    integer :: m, n 

    overlap = .true.
    do m = 1, mer
        do n = 1, mer 
            if ((di%circle(m)%pol == 0) .and. (dj%circle(n)%pol == 0)) then 
                dist = distance(di%circle(m), dj%circle(n), region)
                if ((dist + tol) < ((SIGMAI + SIGMAJ) / 2.)) then
                    if (verbose) then
                        write(*,*) 'check_overlap: distance between pair is', dist
                    endif
                    return
                endif
            endif
        enddo
    enddo
    overlap = .false.
end function check_overlap

! method that initializes the discs in an ordered arrangement
! method cannot fail, unless the area fraction is above the 
! max packing fraction.
! called by group object
! owned by group object (sphere object?)
function lattice (d, region) result (success)
    implicit none
    type(group), dimension(:), intent(out) :: d ! array of groups, that are disks
    real(kind=dbl), intent(in) :: region !	 length of simulation box in each dimension
    type(position) :: center
    real(kind=dbl) :: ori 
    logical :: success ! whether the random walk algorithm generates positions for all discs
    real(kind=dbl) :: n, l ! used for lattice configuration
    integer :: x, y, i, j

    success = .false.
    tsl = 0.

    n = sqrt(real(size(d))) ! number of cells in 1d
    l = (region / ceiling(n)) - tol ! length of each cell
    do x = 1, ceiling(n)
        center%r(1) = l * (real(x) - 0.75)
        do y = 1, ceiling(n)
            center%r(2) = l * (real(y) - 0.75)
            call random_number (ori)
            ori = ori * twopi
            i = (x - 1) * ceiling(n) + y
            call generate_position (d(i), center, ori, region)
        enddo
    enddo

    ! check that none of the discs overlap
    if (check_boundaries()) return
    success = .true.
end function lattice

! method that sets the velocity of all discs along a MB dist
! method utilizes a disc level utility, that sets the MB for the entire part
! called by group object
! owned by sphere object (disc module?)
subroutine set_velocity(d, temp, region)
    implicit none
    type(group), dimension(:), intent(out) :: d ! array of groups, that are disks
    real(kind=dbl), intent(in) :: temp ! used to initialize velocities
    real(kind=dbl), intent(in) :: region
    integer :: i ! used for indexing
    
    tsl = 0.
    do i = 1, size(d)
        call active_velocity(d(i), temp, region)
    enddo

    ! update the circles to their current position
    call update_positions(d, tsl, region)
    ! remove any exess momentum from the system
    call center_momentum(d)
end subroutine set_velocity

! method that sets the velocity of a disc to a velocity with
! constant magnitude in a random direction
! method is called by group_list
! method is stored with group
subroutine active_velocity(d, temp, region)
    implicit none
    type(group), intent(inout) :: d
    real(kind=dbl), intent(in) :: temp
    real(kind=dbl), intent(in) :: region
    real(kind=dbl) :: mag ! magnitude of the vector
    real(kind=dbl) :: theta ! orientation of vector
    integer :: i, q

    ! the magnitude of the velocity is proportional to the temperature
    mag = sqrt(2. * temp)

    ! give the velocity vector a random orientation
    call random_number(theta)
    theta = theta * twopi

    do i = 1, mer 

        ! bring the position of the circles forward in time
        call fastforward(d%circle(i), tsl, region)

        ! assign the velocity vector according to the random orientation
        d%circle(i)%vel%v(1) = cos(theta) * mag 
        d%circle(i)%vel%v(2) = sin(theta) * mag 

        ! return circle with new velocity to tsl position 
        call fastbackward(d%circle(i), tsl, region)
    enddo
end subroutine active_velocity

! method that sets the velocities of all circles in disc along
! a Maxwell-Boltzmann distribution
! method is called by group
! method is stored in disk module
subroutine MB(d, temp, region) 
    implicit none
    type(group), intent(inout) :: d
    real(kind=dbl), intent(in) :: temp
    real(kind=dbl), intent(in) :: region
    real(kind=dbl) :: mu = 0. ! average of gaussian 
    real(kind=dbl) :: sigma ! distribution of Gaussian
    real(kind=dbl) :: u_1, u_2 ! Box-Mueller algorithm parameters
    integer :: i, q

    sigma = sqrt(temp)
    do i = 1, mer 

        ! bring the position of the circles forward in time
        call fastforward(d%circle(i), tsl, region)

        ! use box-mueller algorithm to generate random numbers
        ! along a Gaussian in each dimension
        do q = 1, ndim
            d%circle(i)%vel%v(q) = boxmueller(mu, sigma)
        enddo

        ! return circle with new velocity to tsl position 
        call fastbackward(d%circle(i), tsl, region)
    enddo
end subroutine MB

! function that returns a random number along a box mueller 
! called by module level disk
! stored in math (util??) (dmd constants)
function boxmueller (mu, sigma) result (rand)
    implicit none
    real(kind=dbl) :: rand ! random number returned by algorithm
    real(kind=dbl), intent(in) :: mu ! average of Gaussian
    real(kind=dbl), intent(in) :: sigma ! dist of Gaussian
    real(kind=dbl) :: u_1, u_2 ! random numbers
    
    ! Generate a pseudo random vector components based on a 
    ! Gaussian distribution according the Box-Mueller algorithm

    ! generate
    call random_number (u_1)
    call random_number (u_2)

    ! distribute, assign
    rand = mu + sigma * sqrt( -2. * log(u_1)) * sin (twopi * u_2)
end function boxmueller

! method that determines if a circle resides in the simulation
! box. if it does not, the periodic boundaries are applied
! called by disk module
! stored in circle module
subroutine periodic_boundary_conditions(c, region)
    implicit none
    type(particle), intent(inout) :: c 
    real(kind=dbl), intent(in) :: region
    integer :: q ! used for indexing

    do q = 1, ndim 
        if (c%fpos%r(q) >= region) c%fpos%r(q) = c%fpos%r(q) - region
        if (c%fpos%r(q) < 0.) c%fpos%r(q) = c%fpos%r(q) + region 
    enddo
end subroutine periodic_boundary_conditions

! adjusts the velocities of all groups by setting the 
! linear momentum in each dimension to zero
subroutine center_momentum (d)
    implicit none
    type(group), dimension(:), intent(inout) :: d 
    real(kind=dbl), dimension(ndim) :: vsum
    integer :: i, m, q 

    ! WARNING: do not call outside of setting_velocity method,
    ! when tsl is 0. Otherwise, the particles should be
    ! move forward in time and the periodic boundaries 
    ! considered

    vsum = 0.

    ! calculate the drift of the system
    do i = 1, size(d)
        do m = 1, mer
            do q = 1, ndim 
                vsum(q) = vsum(q) + (d(i)%circle(m)%vel%v(q))
            enddo
        enddo
    enddo

    ! average the linear momentum
    vsum = vsum / real(size(d))

    ! remove the drift from the systen 
    do i = 1, size(d)
        do m = 1, mer 
            do q = 1, ndim 
                d(i)%circle(m)%vel%v(q) = d(i)%circle(m)%vel%v(q) - vsum(q)
            enddo
        enddo
    enddo

    ! WARNING: subroutine assumes that tsl is zero, 
    ! othwise the false position of all of the partcles
    ! should be updated
end subroutine center_momentum

! sets the chirality all all discs according to their
! number fraction
subroutine set_chirality(d)
    implicit none 
    type(group), dimension(:), intent(out) :: d 
    real(kind=dbl), parameter :: xa = 1.
    real(kind=dbl) :: rand 
    integer :: i

    do i = 1, size(d) 
        ! randomly generate a number, check against acceptance
        ! criteria
        call random_number (rand)
        if (rand > xa) then 
            d(i)%chai = 2
        else 
            d(i)%chai = 1
        endif
    enddo
end subroutine set_chirality

! sets the polarity of all the discs according to their order
! called by group (list?)
! stored in disc module
subroutine set_polarity(d)
    implicit none
    type(group), dimension(:), intent(out) :: d
    integer :: i, m 

    do i = 1, size(d) 
        do m = 1, mer 
            d(i)%circle(m)%pol = 0
        enddo
    enddo
end subroutine set_polarity

! method that scheduled collisions for all groups and subpartiles
! according to their interactions 
! called by group_list module / object 
! stored in group_list module / object
subroutine schedule_collisions (d, region)
    implicit none 
    type(group), dimension(:), intent(inout) :: d 
    real(kind=dbl), intent(in) :: region 
    type(id) :: a, b 
    integer :: i, m, n, j1, j2

    ! loop through up list pairs 
    do i = 1, size(d)
        do m = 1, mer
            n = 0
            j1 = id2mol(i, m, mer)
            a = mol2id(j1, mer)
            d(a%one)%circle(a%two)%collision = reset_event()
            upnab: do
                n = n + 1 
                if ((n == (nbrListSizeMax + 1)) .or. &
                    (d(a%one)%circle(a%two)%upnab(n)%one == 0)) exit
                b = d(a%one)%circle(a%two)%upnab(n)
                call predict_collision(d(a%one)%circle(a%two), &
                    d(b%one)%circle(b%two), a, b, tsl, region)
            enddo upnab
        enddo
    enddo
end subroutine schedule_collisions

! method that predicts collisions between two circles that
! are within a disk object
! called by group_list module / object 
! sotre in group module / object
subroutine predict_collision (ci, cj, a, b, tsl, l) 
    implicit none 
    type(particle), intent(inout) :: ci, cj ! uplist and downlist particles 
    type(id), intent(in) :: a, b ! uplist and downlist partners
    real(kind=dbl), intent(in) :: tsl ! time since last update 
    real(kind=dbl), intent(in) :: l ! length of simulation box
    type(position) :: rij
    type(velocity) :: vij 
    type(event) :: prediction
    integer :: j, n, q

    prediction = reset_event()
    if ((a%one > b%one) .or. ((a%one == b%one) .and. &
        (a%two >= b%two))) return

    ! calculate the real position of the pair 
    rij = calculate_pair_position(ci, cj, tsl, l)
    do q = 1, ndim 
        vij%v(q) = ci%vel%v(q) - cj%vel%v(q)
    enddo

    ! determine the type of interaction that should occur
    ! determine when that interaction will occur 
    if (a%one == b%one) then ! if the pair are in the same group
        ! this should not happen
        write (*,*) 'pedict_collision: called bonding methods but no discs are bonded.'
        call exit()
    elseif (b%one > a%one) then ! if the pair are not in the same group
        ! the pair have a hard sphere interaction
        prediction = hardsphere_event (rij, vij, a, b)
    endif

    ! determine if that interaction is sooner than the previous
    if (sooner(prediction, ci%collision)) then 
        ci%collision = prediction
    endif
end subroutine predict_collision

! method that checks if any of the boundaries between groups
! are overlapping
! called by simulation
! stored with group_list 
function check_boundaries () result(overlap)
    implicit none 
    logical :: overlap 
    integer :: i, j
    logical :: position, group 

    position = .false.
    group = .false.

    single: do i = 1, (disk - 1) 

        ! check whether the internal position of the disc
        ! is correct
        if (check_position(d(i), region)) then
            position = .true.
            write (*,010) i 
            010 format(" check_boundaries: overlap in internal position of group",&
                 I4, ".")
        endif
        ! if (position) exit

        pair: do j = (i + 1), disk 

            ! check whether two disks are overlapping 
            if (check_overlap(d(i), d(j), region, tol)) then 
                group = .true.
                write (*, 011) i, j 
                011 format (" check_boundaries: groups ", I4,&
                     " and ", I4, " are overlapping.")
            endif
            ! if (group) exit single

        enddo pair
    enddo single

    ! report to user 
    !if (position) then 
        !write (*,010) i 
        !010 format(" check_overlap: overlap in internal position of group", I4, ".")
    !endif 

    !if (group) then 
        !write (*, 011) i, j 
        !011 format (" check_overlap: groups ", I4, " and ", I4, " are overlapping.")
    !endif

    overlap = position .or. group 
end function check_boundaries


! STATE INITIALIZATION 

! method that initializes the state according to defaults
! called by simulation object / module
! every thing is intialized to zero
subroutine set_dmd(sim)
    implicit none 
    type(dmd), intent(out) :: sim
    
    write (*,*) 'set_dmd: initializing dmd simulation to default, initial state.'
    sim%n_cell = cell
    sim%temp = tempstart
    sim%den = eta 
    sim%anneal = 0
    sim%timenow = 0.
    sim%timeperiod = 0.
    sim%n_events = 0
    sim%n_col = 0
    sim%n_hard = 0
    sim%n_well = 0
    sim%n_bond = 0
    sim%n_ghost = 0
    sim%n_mag = 0
end subroutine set_dmd


! ** PROPERTIES **********************************************

! method used for calculating property values and equilibrium
! owned by property module / type
! called by group object
subroutine accumulate_properties (prop, number)
    implicit none
    ! ** calling variables ***********************************
    type(property), intent(inout) :: prop
    integer, intent(in) :: number
    ! ** local variables *************************************

	if (number == 0) then
        prop%sum = 0.
        prop%sum2 = 0.
        prop%equilavg = 0.
        prop%count = 0
        prop%equilibcount = 0
	else if (number == 1) then 
        prop%sum = 0.
        prop%sum2 = 0.
        prop%count = 0
	else if (number == 2) then
        prop%count = prop%count + 1
        prop%sum = prop%sum + prop%value
        prop%sum2 = prop%sum2 + (prop%value ** 2)
	else if (number == 3) then
        prop%sum = prop%sum / prop%count
        prop%sum2 = sqrt(prop%sum2 / prop%count - (prop%sum ** 2))
	else if (number == 4) then
        prop%equilavg = prop%equilavg + prop%sum 
        prop%equilstd = prop%equilstd + (prop%sum ** 2)
        prop%equilibcount = prop%equilibcount + 1
	else if (number == 5) then 
        prop%equilavg = prop%equilavg / prop%equilibcount
        prop%equilstd = sqrt(prop%equilstd / prop%equilibcount - (prop%equilavg ** 2))
	else
        write(*,*) 'Error in execution of accumulate_properties subroutine'
	end if
end subroutine accumulate_properties

! method that performs the same operation on all properties
! owned by group object
! called by group object
subroutine update_properties(number)
    implicit none
    integer, intent(in) :: number ! operation to perform on
    ! all properties

    call accumulate_properties(te, number)
    call accumulate_properties(pot, number)
    call accumulate_properties(ke, number)
    call accumulate_properties(temp, number)
    call accumulate_properties(nem,number)
    ! call accumulate_properties(perc, number)
end subroutine update_properties

! method that calculates systems properties
! called by simulation object 
! stored in object module
subroutine calculate_properties (d, region)
    implicit none 
    type(group), dimension(:), intent(in) :: d 
    real(kind=dbl), intent(in) :: region ! length of simulation box
    real(kind=dbl) :: vsum, vvsum
    integer :: i, m, q ! indexing parameters

    ! calculate all properties 
    vsum = 0.
    vvsum = 0.

    ! sum the velocities
    do i = 1, size(d)
        do m = 1, mer
            do q = 1, ndim
                !vsum = vsum + d(i)%circle(m)%vel%v(q)
                vvsum = vvsum + (d(i)%circle(m)%vel%v(q) ** 2)
            enddo
        enddo
    enddo

    ! kinetic energy
    ke%value = (0.5 * vvsum) / (real(size(d)) * real(mer))

    ! potential energy
    pot%value = accumulate_potential(d, region)

    ! total energy
    te%value = ke%value + pot%value 

    ! temperature
    temp%value = 2. * ke%value / real(ndim)

    ! nematic op
    nem%value = calculate_nematic(d, region)

    ! percolation op
    ! perc%value = calculate_percolation(d)

    ! update
    call update_properties(2)
end subroutine calculate_properties

! method that calculates the potential of a system of discs
! called by group_list object 
! stored in group objkect (disc module)
function accumulate_potential(d, region) result(p)
    implicit none 
    type(group), dimension(:), intent(in) :: d
    real(kind=dbl), intent(in) :: region ! length of simulation box
    real(kind=dbl) :: p ! system potential
    real(kind=dbl) :: rij ! distance between two interacting circles
    integer :: i, j, m, n, o

    ! initialize
    p = 0.

    ! loop through and calculate the potential betwee each disc
    do i = 1, size(d) - 1
        do m = 1, mer 

            ! loop through all uplist atoms 
            o = 0 
			uplist: do 
                o = o + 1
                j = d(i)%circle(m)%upnab(o)%one 
                if (j == 0) exit 
                if (i /= j) then 
                    n = d(i)%circle(m)%upnab(o)%two 
                    ! TODO group_object method that calcualtes potential between two objects 
                    ! calculate the real distance between the pair
                    rij = distance(d(i)%circle(m), d(j)%circle(n), region)
					if (rij < sigma2) then ! if the pair is within the first, innermost well
						if ((d(i)%circle(m)%pol * d(j)%circle(n)%pol) == -1) then ! if the spheres are attracted to
                            p = p - epsilon2
						else if ((d(i)%circle(m)%pol * d(j)%circle(n)%pol) == 1) then ! if the spheres are repulsed from 
                            p = p + epsilon2
						end if 
					else if (rij < sigma3) then ! if the pair is within the second, middle well 
						if ((d(i)%circle(m)%pol * d(j)%circle(n)%pol) == -1) then ! if the spheres are attracted to
                            p = p - epsilon3
						else if ((d(i)%circle(m)%pol * d(j)%circle(n)%pol) == 1) then ! if the spheres are repulsed from
                            p = p + epsilon3
						end if
					else if (rij < sigma4) then ! the pair is within the third, outermost well
						if ((d(i)%circle(m)%pol * d(j)%circle(n)%pol) == -1) then ! if the spheres are attracted to
                            p = p - epsilon4
						else if ((d(i)%circle(m)%pol * d(j)%circle(n)%pol) == 1) then ! if the spheres are repulsed from
                            p = p + epsilon4
						end if
					end if
				endif
			enddo uplist
        enddo
    enddo

    ! average potential 
    p = p / (real(size(d)) * real(mer))
end function accumulate_potential

! method that calculates the nematic order of the system
! called by group_list module
! store in group_object module
function calculate_nematic (d, region) result(n)
    implicit none 
    type(group), dimension(:), intent(in) :: d 
    real(kind=dbl), intent(in) :: region ! length of simulation box
    real(kind=dbl), dimension(size(d)) :: phi 
    real(kind=dbl) :: n ! nematic order of system 
    integer :: i, j, count

    ! initialize
    n = 0.
    count = 0

    ! calculate angle for each disc
    do i = 1, size(d)
        phi(i) = calculate_orientation(d(i), region, 0.)
    enddo


    ! pair wise comparison of all groups 
    do i = 1, size(d) - 1
        do j = i, size(d)
            n = n + cos(phi(i) - phi(j)) ** 2
            count = count + 1
        enddo
    enddo

    n = n / real(count)
    n = (3. * n - 1.) / 2.
end function calculate_nematic

! method that calcualtes the orientation of a disc relative to
! an angle provided by the user
! called by group_object
! called in group_object module
function calculate_orientation (d, region, theta) result (phi)
    implicit none
    type(group), intent(in) :: d 
    real(kind=dbl), intent(in) :: region ! length of simulation box
    real, intent(in) :: theta ! reduced rotation (0 = 0, 1 = 2pi, etc.)
    integer, parameter :: CIRCLE_1 = mer - 1
    integer, parameter :: CIRCLE_2 = mer
    type(particle), dimension(2) :: c_ori ! two positions used to calculate orientation
    type(position) :: dr 
    real(kind=dbl) :: phi 
    integer :: m, q

    ! the orientation of the particle is the direction of its velocity
    phi = 0.
    do q = 1, ndim 
        phi = phi + (d%circle(1)%vel%v(q) ** 2)
    enddo
    phi = sqrt(phi)
    phi = (d%circle(1)%vel%v(1) / phi)
    phi = acos(phi) ! bounds [-1, 1], range [0, pi] 
    if (d%circle(1)%vel%v(2) < 0.) phi = -phi
    phi = phi + theta * twopi
end function calculate_orientation

! method that calulcates the distance between two circles
! method called by group
! stored in partcile type / module
function distance(ci, cj, region) result(dist)
    implicit none 
    type(particle), intent(in) :: ci, cj ! two circles 
    real(kind=dbl), intent(in) :: region ! length of simulation box
    real(kind=dbl) :: dist ! distance between two circles
    type(position) :: rij
    type(velocity) :: vij
    integer :: q ! indexing parameter 

    ! calculate the real distance between the particle pair
    dist = 0. 
	do q = 1, ndim 
        rij%r(q) = (ci%fpos%r(q) - cj%fpos%r(q))
        vij%v(q) = (ci%vel%v(q) - cj%vel%v(q))
        rij%r(q) = rij%r(q) + vij%v(q) * tsl
        call minimum_image_convention (rij, region)
        dist = dist + (rij%r(q) ** 2)
	end do 
    dist = sqrt (dist)
end function distance

! method the applies the minimum image convention to comparisons
! called by particle module
! store in particle module
subroutine minimum_image_convention (rij, region)
    implicit none 
    type(position), intent(inout) :: rij
    real(kind=dbl), intent(in) :: region ! length of simulation box
    integer :: q

    do q = 1, size(rij%r)
        if (rij%r(q) >= (0.5 * region)) rij%r(q) = rij%r(q) - region
        if (rij%r(q) < (-0.5 * region)) rij%r(q) = rij%r(q) + region
    enddo
end subroutine minimum_image_convention


! ** EVENT SHEDULING *****************************************

! method that reschedules the entire event calander, including
! particle collisions and ghost events 
! called by simulation object
! stored with simulation object
subroutine complete_reschedule ()
    implicit none 
    integer :: i
    type(id) :: a, b 

    ! schedule collisions between particles
    call schedule_collisions(d, region)

    ! schedule when the ghost collisions will occur

    ! add collisions to the event tree
    call initialize_binarytree(eventTree, rootnode)
    do i = 1, mols 
        a = mol2id(i, mer)
        call addbranch(eventTree, i)
    enddo

    ! add ghost events
    if (thermostat) call addbranch (eventTree, mols+1)
    if (field) call addbranch (eventTree, mols+2)
    if (field_actuation) call addbranch (eventTree, mols+3)
    if (movie) call addbranch (eventTree, mols+4)
    if (activity) call addbranch (eventTree, mols+5)
end subroutine complete_reschedule

! method that generates the ghost events
! called by simulation
! stored with simulation
subroutine schedule_ghostcollisions(ghost_event)
    implicit none
    type(event), dimension(ghost), intent(out) :: ghost_event
    integer :: i 

    ghost_event = reset_event()
    do i = 1, ghost
        if (i == 1) then
            ! thermostat
            ghost_event(i) = predict_ghost(thermo_period)
        elseif (i == 2) then 
            ! field 
            ghost_event(i) = predict_ghost(mag_period)
        elseif (i == 3) then 
            ! switch
            ghost_event(i) = predict_ghost(1. / field_actfreq)
        elseif (i == 4) then 
            ! movie
            ghost_event(i) = predict_ghost(movfreq)
        elseif (i == 5) then 
            ! activity 
            ghost_event(i) = predict_ghost(active_period)
        endif
    enddo
end subroutine schedule_ghostcollisions

! function that predicts when a ghost collision will occur
! stored with simulation
! called by simulation
type(event) function predict_ghost(period)
    implicit none
    ! ** calling variables ***********************************
    real(kind=dbl), intent(in) :: period ! average time between collision events
    ! ** local variables *************************************

    predict_ghost%time = period
    predict_ghost%uppart = nullset()
    predict_ghost%dnpart = random_group(disk, mer)
    predict_ghost%type = -1
end function predict_ghost

! method that lists the order of events in the calander
! called by group_list
! stored with event module / object
subroutine report_schedule()
    implicit none 
    type(node), dimension(mols + ghost) :: tree 
    integer :: i, e, incint, root
    type(id) :: a 

    root = rootnode
    tree = eventTree
    incint = mols
    if (thermostat) incint = incint + 1 
    if (field) incint = incint + 1
    if (field_actuation) incint = incint + 1
    if (movie) incint = incint + 1
    if (activity) incint = incint + 1
    do i = 1, incint
        e = findnextevent(tree, root)
        if (e <= mols) then 
            a = mol2id(e, mer)
            call report_event(d(a%one)%circle(a%two)%collision)
            call delbranch(tree, root, e)
        elseif (e > mols) then 
            if ((thermostat .and. ((e - mols) == 1)) .or. (field .and. &
                ((e - mols) == 2)) .or. (field_actuation .and. &
                ((e - mols) == 3)) .or. (movie .and. ((e - mols) == 4)) &
                .or. (activity .and. ((e - mols) == 5))) then
                    call report_ghostevent ((e-mols), ghost_event(e - mols))
                    call delbranch(tree, root, e)
            endif
        endif
    enddo
end subroutine report_schedule

! method that processes a ghost thermostat collision
! called by simulation
! stored with group_list
! TODO: transition part of the routine to group
subroutine thermal_ghost_event(d, temp, region, n_ghost)
    implicit none
    type(group), dimension(:), intent(inout) :: d ! array of group of particles
    real(kind=dbl), intent(in) :: temp ! system temp 
    real(kind=dbl), intent(in) :: region ! length of simultion box
    integer, intent(inout) :: n_ghost ! number of times ghost events have happened
    real(kind=dbl) :: disp ! used for measuring particle displacement 
    type(id) :: a
    integer :: i, j, n, m, q ! indexing parameter 
    type(position), dimension(mer) :: rg ! stores the real position of each particle in ghost coliision

    ! determine how many ghost collisions could take place
    ! the number of particles that experience a ghost collisions
    ! is determined accoring to a random Poisson variate
    n = poissondist(real(1., dbl))
    n_ghost = n_ghost + n
    if (debug >= 1) write (*, 140) n
    140 format (' ', I2, ' thermal ghost collisions occur.')

    if (n > 0) then ! the number of particles that should experience ghost collisions is 1 or greater
        ! TODO: create sub method with group that processes the random collision
        do j = 1, n

            ! select a random particle
            a = random_group(disk, mer)
            i = a%one
            if (debug >= 1) write (*, 141) i
            141 format (' particle ', I5, ' experiences a thermal ghost collision ')

            ! reassign the velocity according to a maxwell boltzmann distribution
            call active_velocity (d(i), temp, region) 

            ! measure the dispalcement of the partcle
            do m = 1, mer 
                disp = 0 
                do q = 1, ndim 
                    disp = disp + (d(i)%circle(m)%vel%v(q) * tsl) ** 2
                enddo 
                if (sqrt(disp) > nbrDispMax) nbrnow = .true.
            end do 

            ! reschedule the event calander for that particle
            call ghost_reschedule(i)

        end do
    end if
end subroutine thermal_ghost_event

! method that processes a ghost field event
! called by simultion
! stored with group_list
subroutine field_ghost_event (d, temp, region, n_mag)
    implicit none 
    type(group), dimension(:), intent(inout) :: d ! array of group of particles
    real(kind=dbl), intent(in) :: temp ! system temp 
    real(kind=dbl), intent(in) :: region ! length of simultion box
    integer, intent(inout) :: n_mag ! number of magnetic events that have occured
    real(kind=dbl) :: disp ! used for measuring particle displacement
    type(id) :: a 
    integer :: n, m, q, j, i ! indexing parameters

    ! determine how many ghost collisions should take place
    ! the number of particles that experience a ghost collisions
    ! is determined accoring to a random Poisson variate
    n = poissondist(real(1., dbl))
    n_mag = n_mag + n
    if (debug >= 1) write (*, 160) n
    160 format (' ', I2, ' field ghost collisions occur.')

    if (n > 0) then 

        do j = 1, n 

            ! select a random group
            a = random_group(disk, mer)
            i = a%one
            if (debug >= 1) write (*, 161) i
            161 format (' particle ', I5, ' experiences a field ghost collision ')

            ! perform the field collision
            call field_group_collision (d(i), temp, region)

            ! determine the displacement
            do m = 1, mer 
                disp = 0 
                do q = 1, ndim 
                    disp = disp + (d(i)%circle(m)%vel%v(q) * tsl) ** 2
                enddo 
                if (sqrt(disp) > nbrDispMax) nbrnow = .true.
            end do 

            ! reschedule the event calander for the particle 
            ! that participated in the event
            call ghost_reschedule(i)
        enddo
    endif
end subroutine field_ghost_event

! method that processes ghost field collision for a group of particles
! called by group_list
! stored with group
subroutine field_group_collision (d, temp, region)
    implicit none
    type(group), intent(inout) :: d ! group
    real(kind=dbl), intent(in) :: temp ! temperature setpoint of the system
    real(kind=dbl), intent(in) :: region ! length of the simulation box 
    real(kind=dbl), dimension(ndim) :: fieldvec ! vector describing the direction of the external field
    integer :: m, q

    ! determine the direction that the field is pointing
    if (field_rotation) then ! if the field is rotating
        ! the direction that the field points depends on the rotational velocity of the field
        fieldvec(1) = cos((sim%timenow / field_angvel) * twopi) * magforce
        fieldvec(2) = sin((sim%timenow / field_angvel) * twopi) * magforce
    else ! otherwise the direction of the field is constant
        ! the field points in the direction of the x axis
        fieldvec(1) = magforce
        fieldvec(2) = 0.
    endif

    ! for each particle in the group
    do m = 1, mer 

        ! only polarized particles participate in the event
        if (d%circle(m)%pol == 0) cycle

        ! update the real position of the particle
        call fastforward (d%circle(m), tsl, region)

        ! assign the field vector to the particle
        ! based on its polarity
        do q = 1, ndim 
            if (d%circle(m)%pol == 1) then 
                d%circle(m)%vel%v(q) = fieldvec(q)
            else
                d%circle(m)%vel%v(q) = -fieldvec(q)
            endif
        enddo

        ! update the false position of the particle
        call fastbackward (d%circle(m), tsl, region)
    enddo
end subroutine field_group_collision

! method that process a ghost activity event
! called by simulation
! stored with group_list
! TODO: transition part of the routine to group
subroutine activity_ghost_event (d, temp, region, n_ghost)
    implicit none 
    type(group), dimension(:), intent(inout) :: d ! array of group of particles
    real(kind=dbl), intent(in) :: temp ! system temp 
    real(kind=dbl), intent(in) :: region ! length of simultion box
    integer, intent(inout) :: n_ghost ! number of times ghost events have happened
    real(kind=dbl) :: disp ! used for measuring particle displacement 
    type(id) :: a
    integer :: n, m, j, i, q ! indexing 

    ! determine how many ghost collisions should take place
    ! the number of particles that experience a ghost collisions
    ! is determined accoring to a random Poisson variate
    n = poissondist(real(1., dbl))
    n_ghost = n_ghost + n
    if (debug >= 1) write (*, 150) n
    150 format (' ', I2, ' acitivity ghost collisions occur.')

    if (n > 0) then 

        do j = 1, n 

            ! select a random group 
            a = random_group(disk, mer)
            i = a%one
            if (debug >= 1) write (*, 151) i
            151 format (' particle ', I5, ' experiences an active ghost collision ')

            ! perform the active collision
            call active_group_collision (d(i), temp, region)

            ! measure the dispalcement of the partcle
            do m = 1, mer 
                disp = 0 
                do q = 1, ndim 
                    disp = disp + (d(i)%circle(m)%vel%v(q) * tsl) ** 2
                enddo 
                if (sqrt(disp) > nbrDispMax) nbrnow = .true.
            end do 

            ! reschedule the event calander for that particle
            call ghost_reschedule(i)

        enddo

    endif
end subroutine activity_ghost_event

! method that processes an active ghost collision for a group of spheres
! stored with group
! called by group_list
subroutine active_group_collision (d, temp, region)
    implicit none 
    type(group), intent(inout) :: d ! group
    real(kind=dbl), intent(in) :: temp ! temperature setpoint of the system
    real(kind=dbl), intent(in) :: region ! length of the simulation box 
    real(kind=dbl) :: magnitude
    type(velocity) :: velvec ! velocity vector
    real(kind=dbl) :: phi ! orientation of the group
    integer :: m, q ! indexing


    ! calculate the orientation of the dipole
    ! relative to the x-axis
    phi = calculate_orientation(d, region, 0.)

    ! apply the activity motion in an orientation that is 
    ! orthogonal to the orientation of the dipole
    magnitude = sqrt(2. * temp)
    velvec%v(1) = cos(phi - halfpi) * magnitude
    velvec%v(2) = sin(phi - halfpi) * magnitude

    if (debug >= 1) then 
        write (*,*) 'The orientation of the partcile is ', phi
        write (*,*) 'The applied vector is: ', velvec%v(1),', ', velvec%v(2)
    endif

    ! apply the vector each the particles in the group
    ! that are apart of the dipole
    do m = 1, mer 
        ! the first particle is not a part of the dipole
        if (m == 1) cycle

        ! update the real position of the particle
        call fastforward (d%circle(m), tsl, region)

        ! assign the velocity vector to the particle
        do q = 1, ndim 
            d%circle(m)%vel%v(q) = velvec%v(q)
        enddo

        ! update the false position of the particle
        call fastbackward (d%circle(m), tsl, region)
    enddo
end subroutine active_group_collision



! ** I / O ***************************************************

! method that opens the files used during the molecular simulation
! owned by simulation module / object
! called by the simulaition module / object
subroutine open_files ()
    implicit none 
    character(len=10), parameter :: format = "(I0.3)"
    character(len=10) :: rp
    character(len=40), parameter :: sphmovfile_end = '_sphmov.xyz' 
    character(len=40), parameter :: objmovfile_end = '_objmov.xyz'
    character(len=40) :: sphmovfile, objmovfile, simfile, reportfile, annealfile, opfile

    write (rp, format) sim%anneal

    ! open the sphere movie file
    sphmovfile = trim(simtitle) // trim(simid) // trim(rp) // trim(sphmovfile_end)
    open (unit = sphcooriounit, file = trim(sphmovfile), status = 'REPLACE')

    ! open the object movie file
    objmovfile = trim(simtitle) // trim(simid) // trim(rp) // trim(objmovfile_end)
    open (unit = objcooriounit, file = trim(objmovfile), status = 'REPLACE')
end subroutine open_files

! method that closes that files used during the molecular simulation
! owned by the simulation module / object
! called by the simulaion module / object
subroutine close_files ()
    implicit none 

    close (unit = sphcooriounit, status = 'KEEP')
    close (unit = objcooriounit, status = 'KEEP')
end subroutine close_files

! ** group - save and load *********************************

! method that saves position vectors in the system to the save file
! belongs to group type
function save_group(d, filename, iounit) result (iobool)
    implicit none
    type(group), dimension(:), intent(in) :: d ! array of discs (groupings of spheres)
    character(len=*), intent(in) :: filename ! path to location to store information
    integer, intent(in) :: iounit ! port to use while saving data
    logical :: iobool ! indicates successful data transfer
    integer :: ierror ! used to determine read/write error status 
    integer :: i, m, q ! used for indexing 

    iobool = .false.
    open (unit = iounit, file = trim(filename), status = 'REPLACE', action = 'WRITE', iostat = ierror)
	if (ierror == 0) then 

        write (*,*) 'save_group: spheres saved to ', trim(filename)

        ! time since last false position update
        write (iounit, *) tsl

		do i = 1, size(d) ! for each object 

            ! write the assigned chirality of the disk
            write(iounit,*) d(i)%chai

			do m = 1, mer ! for each particle 

                ! write the assigned polarity of the disk
                write(iounit, *) d(i)%circle(m)%pol

                ! write the assigned false position
				do q = 1 , ndim ! for each dimension 
                    write (iounit, *) d(i)%circle(m)%fpos%r(q)
                enddo

                ! write the assigned velocity of the particle
                do q = 1, ndim ! for each dimension
                    write (iounit, *) d(i)%circle(m)%vel%v(q)
                enddo

			enddo 
		enddo 
        close (unit = saveiounit, status = 'KEEP')
        iobool = .true.
	else 
        write (*,*) 'save_group: unable to open file. failed to record sphere data.'
	end if 
end function save_group

! method that loads the position vectors for the system from the save file
! belongs to group type
function load_group(d, filename, iounit) result (iobool)
    implicit none
    type(group), dimension(:), intent(out) :: d ! array of discs (groupings of spheres)
    character(len=*), intent(in) :: filename ! path to location to store information
    integer, intent(in) :: iounit ! port to use while saving data
    logical :: iobool ! indicates successful data transfer 
    integer :: ierror ! used to read / write error status
    integer :: i, m, q ! used for indexing

    iobool = .false.
    open (unit = iounit, file = trim(filename), status = 'OLD', action = 'READ', iostat = ierror)
	if (ierror == 0) then 

        write(*,*) 'load_group: spheres were read from', trim(filename)

        ! read the time since the false position of the particles has been updated
        read (iounit, *) tsl

		do i = 1, size(d) ! for each disk

            ! read the assigned chirality of the disk
            read(iounit, *) d(i)%chai

			do m = 1, mer ! for each particle 

                ! read the assigned polarity of the disk
                read (iounit, *) d(i)%circle(m)%pol

                ! read the false position of the particles
    	        do q = 1 , ndim ! for each dimension
                    read (iounit, *) d(i)%circle(m)%fpos%r(q)
    			end do 

                ! read the velocity of the particles
                do q = 1, ndim ! for each dimension
                    read(iounit, *) d(i)%circle(m)%vel%v(q)
                enddo

    		end do 
    	end do
        close (unit = iounit, status = 'KEEP')
        iobool = .true.
    else
        write (*,*) 'load_group: unable to access file. sphere data not loaded.'
    endif
end function load_group

! ** state - save and load ******************************

! method that saves the state of the simulation to the save file
! belongs to dmd object
function save_state (sim, filename, iounit) result (iobool)
    implicit none
    type(dmd), intent(in) :: sim ! simulation state
    character(len=*), intent(in) :: filename ! path to location to store information
    integer, intent(in) :: iounit ! port to use while saving data
    logical :: iobool ! indicates successful data transfer 
    integer :: ioerror ! used to read / write error status

    iobool = .false.
    open (unit = iounit, file = trim(filename), status = 'REPLACE', action = 'WRITE', iostat = ioerror)
    write (*,*) 'save_state: saving state of dmd simulation to ', trim(filename)
    if (ioerror == 0) then 
        write(iounit, *) sim%n_cell
        write(iounit, *) sim%temp
        write(iounit, *) sim%den 
        write(iounit, *) sim%anneal 
        write(iounit, *) sim%timenow
        write(iounit, *) sim%timeperiod
        write(iounit, *) sim%n_events
        write(iounit, *) sim%n_col
        write(iounit, *) sim%n_hard
        write(iounit, *) sim%n_well 
        write(iounit, *) sim%n_bond
        write(iounit, *) sim%n_ghost
        write(iounit, *) sim%n_mag

        close(unit = iounit, status = 'KEEP')
        iobool = .true.
    else
        write(*,*) 'save_state: unable to open file. failed to record state of dmd simulation.'
    endif
end function save_state

! method that loads the state of the simulation from the save file
! belongs to dmd object
function load_state (sim, filename, iounit) result (iobool)
    implicit none
    type(dmd), intent(out) :: sim ! simulation state
    character(len=*), intent(in) :: filename ! path to location to store information
    integer, intent(in) :: iounit ! port to use while saving data
    logical :: iobool ! indicates successful data transfer 
    integer :: ioerror ! used to read / write error status

    iobool = .false.
    open (unit = iounit, file = trim(filename), status = 'OLD', action = 'READ', iostat = ioerror)
    write (*,*) 'load_state: loading state of dmd simulation from ', trim(filename)
    if (ioerror == 0) then 
        read(iounit, *) sim%n_cell
        read(iounit, *) sim%temp
        read(iounit, *) sim%den 
        read(iounit, *) sim%anneal 
        read(iounit, *) sim%timenow
        read(iounit, *) sim%timeperiod
        read(iounit, *) sim%n_events
        read(iounit, *) sim%n_col
        read(iounit, *) sim%n_hard
        read(iounit, *) sim%n_well 
        read(iounit, *) sim%n_bond
        read(iounit, *) sim%n_ghost
        read(iounit, *) sim%n_mag

        close(unit = iounit, status = 'KEEP')
        iobool = .true.
    else
        write(*,*) 'load_state: unable to open file. failed to load state of dmd simulation.'
    endif
end function load_state

! ** snapshot generation *************************************

! method that generates snapshot of all dmd circles
! called by group_list (?)
! owned by circles (?)
subroutine record_position_circles(d, iounit, region)
    implicit none 
    type(group), dimension(:), intent(in) :: d 
    integer, intent(in) :: iounit 
    real(kind=dbl), intent(in) :: region ! length of simulation box
    type(particle) :: ri
    character(len=15), parameter :: positive = blue
    character(len=15), parameter :: negative = red 
    character(len=15), parameter :: neutral = white 
    integer :: i, m, q ! indexing parameters
    character(len=20) :: format, string, charge, rad

    format = "(2(' ', F7.3))"

    write(iounit,*) mols 
    write(iounit,*) ' x, y, polcolor (3)' ! comment line frame number starting from 1

    do i = 1, size(d)
        do m = 1, mer 
            ! calculate the real position of the particle
            do q = 1, ndim
                ri%fpos%r(q) = d(i)%circle(m)%fpos%r(q) + d(i)%circle(m)%vel%v(q) * tsl
            enddo
            call periodic_boundary_conditions (ri, region)

            ! format the position
            write(string, format) ri%fpos%r(1), ri%fpos%r(2)
            ! determine the chairality of the sphere  
            ! determine the polarization of the sphere 
            select case (d(i)%circle(m)%pol)
			case (-1)
                charge = negative
                write (rad, '(" ", F5.3)') (sigma1 / 2.)
			case (1) 
                charge = positive
                write (rad, '(" ", F5.3)') (sigma1 / 2.)
			case default 
                charge = neutral
                write (rad, '(" ", F5.3)') (sigma1 / 1.)
			end select 
            ! write the description 
            write (iounit, *) trim(string), trim(charge), trim(rad)
        enddo
    enddo
end subroutine record_position_circles

! method that writes xyz style snapshot of groups
! method called by group_list object / module
! method stored in group_object / module
subroutine record_position_objects(d, iounit, region)
    implicit none 
    type(group), dimension(:), intent(in) :: d 
    integer, intent(in) :: iounit
    real(kind=dbl), intent(in) :: region ! length of simulation box
    type(particle) :: center
    real(kind=dbl) :: phi
    real(kind=dbl) :: x, y, z, w ! components of a quaternion
    integer :: i, m, q
    character(len=20) :: num, format, string, typenum, typecol, qxy, qzw, rad, latx, laty

    format = "(2(' ', F7.3))"

    write(iounit,*) size(d)

    ! write simulation box lengths to xyz file
    write (latx, format) region, 0.
    write (laty, format) 0., region
    write(iounit,*) 'Lattice="', trim(latx), trim(laty),'"' ! comment line frame number starting from 1
    do i = 1, size(d)
        ! calculate the real, central postion of the particle
        do q = 1, ndim 
            center%fpos%r(q) = d(i)%circle(1)%fpos%r(q) + d(i)%circle(1)%vel%v(q) * tsl
        enddo
        call periodic_boundary_conditions (center, region)

        ! calculate the orientation of the particle
        phi = calculate_orientation(d(i), region, 0.)

        ! calculate the quaternion of the particle
        call calculate_quaternion2D(x, y, z, w, phi)

        ! write to file
        write (num, '(" ", I4)') i 
        write (string, format) center%fpos%r(1), center%fpos%r(2)
        write (qxy, format) x, y 
        write (qzw, format) z, w 
        write (rad, '(" ", F5.3)') sigma1
        write (typenum, '(" ", I4)') d(i)%chai
        write (iounit, *) trim(num), trim(string), trim(qxy), trim(qzw), trim(typenum), trim(rad)
    enddo
end subroutine record_position_objects

! function that calculates the two dimensional quaternion
! of an objet provided its angle 
! called by the simulation
! stored with the quaternion object (?)
subroutine calculate_quaternion2D (x, y, z, w, phi)
    implicit none 
    real(kind=dbl), intent(inout) :: w, x, y, z ! quaternion of object 
    real(kind=dbl), intent(in) :: phi ! angle of object
    real(kind=dbl) :: t ! used for quaternion calculation
    real(kind=dbl), dimension(3,3) :: rotate ! rotation matrix

    ! calculate the rotation matrix of the square using the single angle
    rotate(1,1) = cos(phi)
    rotate(1,2) = sin(phi)
    rotate(1,3) = 0.
    rotate(2,1) = -sin(phi)
    rotate(2,2) = cos(phi)
    rotate(2,3) = 0.
    rotate(3,1) = 0. 
    rotate(3,2) = 0.
    rotate(3,3) = 1.

    ! calculate the quaternion using the rotation matrix
    x = 0.
    y = 0.
    z = 0.
    w = 0. 
    if (rotate(3,3) < 0.) then 
        ! won't occur
    else
        if (rotate(1,1) < -rotate(2,2)) then 
            t = 1 - rotate(1,1) - rotate(2,2) + rotate(3,3)
            x = rotate(3,1) + rotate(1,3)
            y = rotate(2,3) + rotate(3,2)
            z = t
            w = rotate(1,2) - rotate(2,1)
        else
            t = 1 + rotate(1,1) + rotate(2,2) + rotate(3,3)
            x = rotate(2,3) - rotate(3,2)
            y = rotate(3,1) - rotate(1,3)
            z = rotate(1,2) - rotate(2,1)
            w = t
        endif
    endif
    t = 0.5 / sqrt(t)
    x = t * x 
    y = t * y 
    z = t * z 
    w = t * w 
end subroutine calculate_quaternion2D


! ** report functions **

! method that reports the state of a grouping of circles
! used primarily for debugging
! called by group_list object / module
! stored with group object / module
subroutine report_groups (d)
    implicit none 
    type(group), dimension(:), intent(in) :: d 
    integer :: i 

    ! write out the state of all groups 
    do i = 1, size(d)
        write (*,'(" DISK: ", I4)') i 
        write (*, '(" chirality: ", I2)') d(i)%chai
        write (*,*)
        call report_circles(d(i)%circle)
        write (*,*) 
    enddo
end subroutine report_groups

! method that reports the state of a list of circles
! used primarily for debugging
! called by group object / module
! stored with circle object / module
! TODO create circle and false pos modules of the same
! subroutine, one that requires tsl and one that doesnt
subroutine report_circles (c)
    implicit none
    type(particle), dimension(:), intent(in) :: c 
    integer :: i

    ! write out the state of all circles in the list 
    do i = 1, size(c)
        write (*, '( " CIRCLE: ", I4)') i
        call report_realposition(c(i)%fpos, c(i)%vel, tsl)
        call report_velocity(c(i)%vel)
        call report_polarization(c(i)%pol)
        call report_event(c(i)%collision)
        write (*,*)
    enddo
end subroutine report_circles

! method that reports the false position of a circle
! according to its false position 
! called by circle method
! stored within position object / module 
subroutine report_falseposition (pos)
    implicit none 
    type(position), intent(in) :: pos 
    character(len=20), parameter :: header = "position:"
    character(len=20), parameter :: format_string = "(2(' ', F7.3))"
    character(len=20) :: pos_string
    
    write (pos_string,format_string) pos%r(1), pos%r(2)
    write (*,*) trim(header) // trim(pos_string)
end subroutine report_falseposition

subroutine report_realposition (fpos, vel, tsl)
    implicit none 
    type(position), intent(in) :: fpos 
    type(velocity), intent(in) :: vel
    real(kind=dbl), intent(in) :: tsl 
    character(len=20), parameter :: header = "position:"
    character(len=20), parameter :: format_string = "(2(' ', F7.3))"
    type(position) :: realpos
    character(len=20) :: pos_string
    integer :: q

    ! TODO turn fastforward / backward subroutines into functions
    do q = 1, ndim 
        realpos%r(q) = fpos%r(q) + vel%v(q) * tsl
    enddo
    write (pos_string,format_string) realpos%r(1), realpos%r(2)
    write (*,*) trim(header) // trim(pos_string)
end subroutine report_realposition

! method that reports the velocity of a circle
! called by circle / method
! stored within velocity module / object
subroutine report_velocity (vel)
    implicit none 
    type(velocity), intent(in) :: vel 
    character(len=20), parameter :: header = "velocity:"
    character(len=20), parameter :: format_string = "(2(' ', F7.3))"
    character(len=20) :: vel_string 

    write (vel_string, format_string) vel%v(1), vel%v(2)
    write (*,*) trim(header) // trim(vel_string)
end subroutine report_velocity

! method that reports the polarization of a square 
! method called by circle module / object
! method store in the circle module / object
subroutine report_polarization (pol)
    implicit none 
    integer, intent(in) :: pol 
    character(len=20), parameter :: header = "polarization:"
    character(len=20), parameter :: format_string = "(' ', I2)"
    character(len=20) :: pol_string 

    write (pol_string, format_string) pol
    write (*,*) trim(header) // trim(pol_string)
end subroutine report_polarization

! method that reports the scheduled event of one particle
! owned by circle object / method
! stored within circle object / method
subroutine report_event (e)
    implicit none 
    type(id) :: id ! id corresponding to downlist event partner
    type(event), intent(in) :: e ! event
    character(len=20), parameter :: header = "collision:"
    character(len=120) :: event_string 

    if (e%type > 0) then 
        write (event_string, "(' event (dnlist: ', I3,' of group ', &
            I4,'; uplist: ', I3,' of group ', I4') in ',F8.5, &
            ' seconds (collision type ', I2,').')") e%dnpart%two, &
            e%dnpart%one, e%uppart%two, e%uppart%one, e%time, e%type 
        write (*,*) trim(header) // trim(event_string)
    elseif (e%type < 0) then 
        write (event_string, "(' ghost event in ',F8.5,' seconds (ghost collision type ', &
            I2,').')") e%time, e%type 
        write (*,*) trim(header) // trim(event_string)
    endif
end subroutine report_event

! method that reports any ghost event passed to the method
! owned by the simulation module?
! stored with simulation / ghost event module
subroutine report_ghostevent (id, ghost_event) 
    implicit none 
    integer, intent(in) :: id
    type(event), intent(in) :: ghost_event 
    character(len=20), parameter :: header = "ghost_event:"
    character(len=100) :: event_string 

    if (id == 1) then 
        write (event_string, "(' thermal ghost collision in ', F8.5 &
            , ' seconds (ghost collision type ', I2,').')") ghost_event%time, id
        write (*,*) trim(header) // trim(event_string)
    else if (id == 2) then 
        write (event_string, "(' magnetic ghost collision in ', F8.5 &
            , ' seconds (ghost collision type ', I2,').')") ghost_event%time, id
        write (*,*) trim(header) // trim(event_string)
    else if (id == 3) then 
        write (event_string, "(' field actuation ghost collision in ', F8.5 &
            , ' seconds (ghost collision type ', I2,').')") ghost_event%time, id
        write (*,*) trim(header) // trim(event_string)
    else if (id == 4) then 
        write (event_string, "(' movie snapshot event in ', F8.5 &
            , ' seconds (ghost collision type ', I2,').')") ghost_event%time, id
        write (*,*) trim(header) // trim(event_string)
    else if (id == 5) then 
        write (event_string, "(' active particle ghost event in ', F8.5 &
            , ' seconds (ghost collision type ', I2,').')") ghost_event%time, id
        write (*,*) trim(header) // trim(event_string)
    else
        write (*,*) 'TODO: report ghost event ', id 
    endif
end subroutine report_ghostevent

! write simulations state to csv

! write properties to csv (should this be done individually for each property?)


! ** EFFICIENCY METHODS **************************************

! ** binary tree **

subroutine addbranch(tree, newnode)
    implicit none
    ! ** calling variables ***********************************
    type(node), dimension(mols+ghost), intent(inout) :: tree ! binary tree
    integer, intent(in) :: newnode ! node to be added to tree
    ! ** local variables *************************************
    type(id) :: a ! id of particle corresponding to node
    real(kind=dbl) :: tnew, tcomp ! collision time of new and compartison nodes
    integer :: ncomp ! current insertion node, which is being compared to the new node
    logical :: nfound ! true if insert position of new branch has been found, else false

    if (rootnode == 0) then ! if the tree is empty
        ! establish the root node 
        rootnode = newnode 
    else ! the tree has already begun 
        if ((tree(newnode)%pnode /= 0) .or. (newnode == rootnode)) then ! the node is already in the tree
            ! remove the node from the tree
            call delbranch (tree, rootnode, newnode)
        endif
        ! save the time of the next event
        if (newnode <= mols) then ! for collision events
            a = mol2id(newnode, mer)
            tnew = d(a%one)%circle(a%two)%collision%time
        elseif (newnode > mols) then ! ghost events
            tnew = ghost_event(newnode - mols)%time
        endif

        ! search through all branches to find the proper insert position
        ncomp = rootnode
        nfound = .false.
        do 
            if (nfound) exit
            ! calculate the time of comparison node's event
            if (ncomp <= mols) then 
                a = mol2id(ncomp, mer) ! for collision events
                tcomp = d(a%one)%circle(a%two)%collision%time
            elseif (ncomp > mols) then ! for ghost events
                tcomp = ghost_event(ncomp - mols)%time
            endif

            if (tnew <= tcomp) then ! if the new event is sooner than the current node event
                ! go to the left 
                if (tree(ncomp)%lnode /= 0) then ! if the current node has a left node 
                    ncomp = tree(ncomp)%lnode ! compare the new node to the left branch of the current node
                else ! connect the new node to the left of current node
                    tree(ncomp)%lnode = newnode
                    nfound = .true.
                    exit
                endif
            else ! the new event is later than the current node event
                ! go to the right 
                if (tree(ncomp)%rnode /= 0) then ! if the current node has a right node 
                    ncomp = tree(ncomp)%rnode ! compare the new node to the right branch of the current node
                else ! connect the new node to the right of current node
                    tree(ncomp)%rnode = newnode
                    nfound = .true.
                    exit
                end if
            endif
        enddo 
        ! link the new node to the previous node
        tree(newnode)%pnode = ncomp
    endif
end subroutine addbranch

! ** neighbor + cell list method **

! builds linked list of particles
subroutine build_linkedlist (d, list)
    implicit none 
    type(group), dimension(:), intent(inout) :: d
    integer, dimension(mols + (nCells ** ndim)), intent(out) :: list 
    integer :: i, m, q ! indexing parameters 
    integer, dimension (ndim) :: c ! used to record the cell number of a particle in each dimenion 
    integer :: cellindex ! used to record the absolute cell number of each particle 

    ! reset the pointer for each cell to -1 
    do i = (mols + 1), (mols + (nCells ** ndim))
        list(i) = -1
    enddo

    ! generate linked list according to the position 
    ! of each particle
    do i = 1, size(d)
        do m = 1, mer 
            do q = 1, ndim 
                c(q) = floor (d(i)%circle(m)%fpos%r(q) / lengthCell)
            enddo
            cellindex = (c(2) * nCells) + c(1) + mols + 1
            list(id2mol(i, m,mer)) = list(cellindex)
            list(cellindex) = id2mol(i, m, mer)
        enddo
    enddo
end subroutine build_linkedlist

! builds list of neighboring particles
subroutine build_neighborlist(d, region)
    implicit none 
    type(group), dimension(:), intent(inout) :: d 
    real(kind=dbl), intent(in) :: region ! length of simulation box
    logical, parameter :: debug_neighborlist = (debug >= 1) .and. (.false.)
    integer, dimension(mols + (nCells ** ndim)) :: cellList ! linked list
    integer :: m1x, m1y, m1cell ! integers used to determine the reference cell 
    integer :: delx, dely ! used to shift 
    integer :: m2x, m2y, m2cell ! integers used to determine the search cell 
    real(kind=dbl) :: rij
    integer :: j1, j2 ! particle pair 
    type(id) :: a, b
    integer, dimension (mols) :: nup, ndn ! down and uplist indexes 
    integer :: i, m ! used for indexing

    ! update the position of all particles
    call update_positions (d, tsl, region)

    ! build the linked list
    call build_linkedlist(d, cellList)

    ! reset neighbor list for each sphere
    do i = 1, size(d)
        do m = 1, mer
            d(i)%circle(m)%upnab = nullset()
            d(i)%circle(m)%dnnab = nullset()
        enddo
    enddo 

    ! intialize list index to the first position for all particles 
    nup = 1
    ndn = 1
    do m1y = 1, nCells 
        do m1x = 1, nCells
            ! determine the reference cell 
            m1cell = ((m1y - 1) * nCells) + m1x + mols
            ! establish reference atoms 
            j1 = cellList (m1cell)
            reference_cell: do ! loop through each atom in reference cell 
                if (j1 <= 0) exit ! until the end of the linked list (-1) is reached 
                ! calculate groupid from linked list index 
                a = mol2id(j1, mer)
                ! loop through all neighboring cells 
                do dely = -1, 1 
                    ! determine m2y including wrap around effects 
                     m2y = m1y + dely 
                    if (m2y > nCells) then 
                        m2y = 1
                    else if (m2y < 1) then 
                        m2y = nCells 
                    end if 
                    do delx = -1, 1 
                        ! determine m2x including wrap around effects 
                        m2x = m1x + delx 
                        if (m2x > nCells) then 
                            m2x = 1
                        else if (m2x < 1) then 
                            m2x = nCells
                        end if 
                        ! calculate the neighvoring cell to search through 
                        m2cell = ((m2y - 1) * nCells) + m2x + mols
                        ! establish atom in neighboring cell to compare to j1
                        j2 = cellList (m2cell)
                        compare: do 
                            if (j2 <= 0) exit 
                            if (j2 > j1) then ! for all uplist atoms 
                                b = mol2id(j2, mer)
                                rij = distance(d(a%one)%circle(a%two), d(b%one)%circle(b%two), region) ! calculate the distance between pair i and j 
                                ! determine if the pair are neighbors
                                if (rij < nbrRadius) then ! if the distance between pair i and j is less than the max search radius
                                    if ((nup(j1) > nbrListSizeMax) .or. (ndn(j2) > nbrListSizeMax)) then 
                                        write (*,*) 'build_neighborlist: too many neighbors in neighborlist.',&
                                            ' increase size of list for this desnity.'
                                        call report_neighborlist(d)
                                        call exit ()
                                    endif
                                    ! j is uplist of i 
                                    d(a%one)%circle(a%two)%upnab(nup(j1)) = b
                                    nup (j1) = nup(j1) + 1
                                    ! i is downlist of j 
                                    d(b%one)%circle(b%two)%dnnab(ndn(j2)) = a
                                    ndn (j2) = ndn(j2) + 1
                                endif
                            endif
                            j2 = cellList(j2)
                        enddo compare
                    enddo
                enddo
                j1 = cellList(j1)
            enddo reference_cell
        enddo
    enddo
    nbrnow = .false.
    dispTotal = 0.

    if (debug_neighborlist) call list_neighborlist(d)
end subroutine build_neighborlist

! method that lists the up and downlist neighbors of every particle
! called by build_neighborlist, meant for debugging
! owner by neighborlist efficiency method 
subroutine report_neighborlist (d) 
    implicit none
    type(group), dimension(:), intent(in) :: d 
    logical, parameter :: verbose = (debug >= 2)
    type(id) :: a, b 
    integer :: i, m, up, dn
    real(kind=dbl) :: avg

    ! calculate the average number of particles in each spheres neighborlist
    avg = 0.

    do i = 1, size(d)
        do m = 1, size(d(i)%circle)
            up = 0
            do 
                avg = avg + 1.
                up = up + 1 
                if ((up == nbrListSizeMax) .or. (d(i)%circle(m)%upnab(up)%one == 0)) exit
                b = d(i)%circle(m)%upnab(up)
                if (verbose) write (*, 001) b%two, b%one, m, i
            enddo

            dn = 0
            do
                avg = avg + 1.
                dn = dn + 1 
                if ((dn == nbrListSizeMax) .or. (d(i)%circle(m)%dnnab(dn)%one == 0)) exit 
                b = d(i)%circle(m)%dnnab(dn)
                if (verbose) write (*, 002) b%two, b%one, m, i
            enddo
            write (*,003) m, i, (up - 1), (dn - 1) 
        enddo
    enddo

    ! calculate the average number of neighbors per particle
    avg = avg / real(size(d))
    write (*,004) avg

    001 format (' list_neighborlist: circle ', I2,' of disk ', I4,&
        ' uplist neighbor of circle ', I2, ' of disk ',  I4)
    002 format (' list_neighborlist: circle ', I2,' of disk ', I4,&
        ' downlist neighbor of circle ', I2, ' of disk ',  I4)
    003 format (' list_neighborlist: circle ', I2, ' of disk ', I4, ' has ',&
        I4, ' uplist neighbors and ', I4, ' downlist neighbors.')
    004 format ('report_neighborlist: ', F5.2, ' neighbors per particle.')
end subroutine report_neighborlist

! outputs the neightbor list parameters to the use
! called by build neighborlist, debugging
! owned by neighborlist efficiency method
subroutine state_neighborlist_parameters ()
    implicit none
end subroutine state_neighborlist_parameters

! ** false position method **

! TODO: add method that checks the displacement of any
! particle and then updates the nbrnow parameter if
! it surpasses the max dist.

! returns the position off all the circles to their current
! position, updates time since last (tsl) to 0
! called by group object
! stored in group object (disc module? efficiency module)
subroutine update_positions(d, tsl, region)
    implicit none 
    type(group), dimension(:), intent(inout) :: d 
    real(kind=dbl), intent(inout) :: tsl ! time since last update
    real(kind=dbl), intent(in) :: region
    integer :: i, m ! used for indexing

    do i = 1, size(d)
        do m = 1, mer
            call fastforward(d(i)%circle(m), tsl, region)
            call periodic_boundary_conditions(d(i)%circle(m), region)
        enddo
    enddo
    tsl = 0.
end subroutine update_positions

! subroutine that calculates the pair position vectors
! between two circles based on their false positions
! method called by group 
! method stored in the false position module
function calculate_pair_position (ci, cj, tsl, l) result (rij)
    implicit none 
    type(particle), intent(in) :: ci, cj ! pair of circles to compare
    real(kind=dbl), intent(in) :: tsl ! time since last update
    real(kind=dbl), intent(in) :: l ! length of the simulation box 
    type(position) :: rij ! real position vector between two circles
    type(velocity) :: vij 
    integer :: q ! used for indexing

    do q = 1, ndim 
        rij%r(q) = ci%fpos%r(q) - cj%fpos%r(q)
        vij%v(q) = ci%vel%v(q) - cj%vel%v(q)
        rij%r(q) = rij%r(q) + vij%v(q) * tsl 
        call minimum_image_convention(rij, l)
    enddo
end function calculate_pair_position

! method that returns the current position of a circle to its
! tsl position 
! called by disk module, 
! stored in circle module
subroutine fastbackward(c, tsl, region)
    implicit none
    type(particle), intent(inout) :: c 
    real(kind=dbl), intent(in) :: tsl 
    real(kind=dbl), intent(in) :: region
    integer :: i ! used for indexing

    do i = 1, ndim 
        c%fpos%r(i) = c%fpos%r(i) - c%vel%v(i) * tsl
    enddo
    call periodic_boundary_conditions(c, region)
end subroutine fastbackward

! method that brings the position of a circle from its tsl
! position to the current position
! called by disk module, 
! stored in circle module
subroutine fastforward(c, tsl, region)
    implicit none
    type(particle), intent(inout) :: c
    real(kind=dbl), intent(in) :: tsl 
    real(kind=dbl), intent(in) :: region
    integer :: i ! used for indexing

    do i = 1, ndim 
        c%fpos%r(i) = c%fpos%r(i) + c%vel%v(i) * tsl
    enddo
    call periodic_boundary_conditions(c, region)
end subroutine fastforward

! ** EVENT FUNCTIONS *****************************************

type(event) function reset_event () 
    implicit none 
    ! ** calling variables ***********************************
    ! ** local variables *************************************

    reset_event%time = bigtime
    reset_event%type = 0
    reset_event%dnpart = nullset()
    reset_event%uppart = nullset()
end function reset_event

logical function sooner (newevent, oldevent) 
    implicit none 
    ! ** calling variables ***********************************
    type(event), intent(in) :: newevent
    type(event), intent(in) :: oldevent
    ! ** local variables *************************************

    sooner = newevent%time < oldevent%time 
end function sooner

! ** disk interactions ***************************************

! function that determines whether two circles are neighbor bonded
! called by group object / module
! owned by group??
logical function neighborbonded (a, b)
    implicit none 
    type(id), intent(in) :: a, b 

    ! assume false 
    neighborbonded = .false.
    if ((a%one == b%one) .and. (a%two == 1)) then 
        ! if the two spheres are in the same group
        ! and next to one another
        neighborbonded = .true.
    end if 
end function neighborbonded

! function that determines whether two circles are acrossed bonded
! called by group object / module
! owned by group ??
logical function acrossbonded (a, b)
    implicit none 
    type(id), intent(in) :: a, b 

    ! assume false 
    acrossbonded = .false.
    if ((a%one == b%one) .and. (a%two /= 1)) then 
        ! if the two spheres are in the same group
        ! and across from one another
        acrossbonded = .true.
    end if 
end function acrossbonded

logical function hardsphereoverlap (distance)
    implicit none 
    ! ** calling variables ***********************************
    real(kind=dbl), intent(in) :: distance
    ! ** local variables *************************************

    hardsphereoverlap = distance < (sigma1 - tol)
end function hardsphereoverlap

logical function acrossbondoverlap (distance) 
    implicit none 
    ! ** calling variables ***********************************
    real(kind=dbl), intent(in) :: distance 
    ! ** local variables *************************************

    acrossbondoverlap = (distance > (ocbond + tol)) .or. (distance < (icbond - tol))
end function acrossbondoverlap

logical function neighborbondoverlap (distance) 
    implicit none 
    ! ** calling variables ***********************************
    real(kind=dbl), intent(in) :: distance 
    ! ** local variables *************************************

    neighborbondoverlap = (distance > (onbond + tol)) .or. (distance < (inbond - tol))
end function neighborbondoverlap

! method that predicts when a neighborbonded event will occur
! assumes that the pair meet the conditions for a neighbor
! bonded event to occur
! method called by group object / module
! method owned by group object / module 
type(event) function neighborbond_event (rij, vij, a, b)
    implicit none  
    type(position), intent(in) :: rij 
    type(velocity), intent(in) :: vij 
    type(id), intent(in) :: a, b ! downlist, uplist partners
    real(kind=dbl) :: aij, bij, icij, ocij ! quadratic equation constants wrt inner and outer neighbor bond 
    real(kind=dbl) :: idiscr, odiscr ! discrimenant wrt icij and ocij respectively

    ! assign event partners
    neighborbond_event%dnpart = a
    neighborbond_event%uppart = b

    ! calculate quadratic parameters 
    aij = (vij%v(1) ** 2) + (vij%v(2) ** 2)
    bij = (rij%r(1) * vij%v(1)) + (rij%r(2) * vij%v(2))
    icij = (rij%r(1) ** 2) + (rij%r(2) ** 2) - (inbond ** 2)
    ocij = (rij%r(1) ** 2) + (rij%r(2) ** 2) - (onbond ** 2)
    idiscr = (bij ** 2) - (aij * icij)
    odiscr = (bij ** 2) - (aij * ocij)

    ! determine the event type and time 
    if (bij < 0.) then ! the centers are approaching
        if (idiscr > 0.) then ! the centers will collide 
            neighborbond_event%type = 8 ! an event will occur at the inner bond length
            neighborbond_event%time = (-bij - sqrt(idiscr)) / aij
        else ! the repulsive centers will miss each other 
            neighborbond_event%type = 9 ! an event will take place at the outer bond length
            neighborbond_event%time = (-bij + sqrt(odiscr)) / aij
        end if 
    else ! the centers are receding
        neighborbond_event%type = 9 ! an event will occur at the outer bond length 
        neighborbond_event%time = (-bij + sqrt(odiscr)) / aij
    end if 
end function neighborbond_event

! method that predicts when an across bonded event will occur
! assumes that the pair meet the conditions for a cross
! bonded event to occur
! method called by group object / module
! method owned by group objet / module 
type(event) function acrossbond_event (rij, vij, a, b)
    implicit none  
    type(position), intent(in) :: rij 
    type(velocity), intent(in) :: vij 
    type(id), intent(in) :: a, b ! dnlist, uplist partners
    real(kind=dbl) :: aij, bij, icij, ocij ! quadratic equation constants wrt inner and outer neighbor bond 
    real(kind=dbl) :: idiscr, odiscr ! discrimenant wrt icij and ocij respectively

    ! determine the event partner
    acrossbond_event%dnpart = a 
    acrossbond_event%uppart = b

    ! calculate quadratic parameters 
    aij = (vij%v(1) ** 2) + (vij%v(2) ** 2)
    bij = (rij%r(1) * vij%v(1)) + (rij%r(2) * vij%v(2))
    icij = (rij%r(1) ** 2) + (rij%r(2) ** 2) - (icbond ** 2)
    ocij = (rij%r(1) ** 2) + (rij%r(2) ** 2) - (ocbond ** 2)
    idiscr = (bij ** 2) - (aij * icij)
    odiscr = (bij ** 2) - (aij * ocij)

    ! determine the event type and time 
    if (bij < 0.) then ! the centers are approaching
        if (idiscr > 0.) then ! the centers will collide 
            acrossbond_event%type = 10 ! an event will occur at the inner bond length
            acrossbond_event%time = (-bij - sqrt(idiscr)) / aij
        else ! the repulsive centers will miss each other 
            acrossbond_event%type = 11 ! an event will take place at the outer bond length
            acrossbond_event%time = (-bij + sqrt(odiscr)) / aij
        end if 
    else ! the centers are receding
        acrossbond_event%type = 11 ! an event will occur at the outer bond length 
        acrossbond_event%time = (-bij + sqrt(odiscr)) / aij
    end if 
end function acrossbond_event

! method that predicts when a hard sphere event will
! occur between two circles
! method called by group object / module
! method owned by circle object / module
type(event) function hardsphere_event (rij, vij, a, b)
    implicit none 
    type(position), intent(in) :: rij 
    type(velocity), intent(in) :: vij 
    type(id), intent(in) :: a, b ! dnlist, uplist partners
    real(kind=dbl) :: aij, bij, cij, discr ! quadratic equation constants, discrimenant

    ! determine the event partner
    hardsphere_event%dnpart = a
    hardsphere_event%uppart = b

    ! calculate quadratic parameters 
    aij = (vij%v(1) ** 2) + (vij%v(2) ** 2)
    bij = (rij%r(1) * vij%v(1)) + (rij%r(2) * vij%v(2))
    cij = (rij%r(1) ** 2) + (rij%r(2) ** 2) - sg1sq
    discr = (bij ** 2) - (aij * cij)

    ! predict if the hard spheres will collide
    if ((discr > 0.) .and. (bij < 0.)) then
        ! calculate the time until collision
        hardsphere_event%type = 1 ! an event will occur at sigma1
        hardsphere_event%time = (-bij - sqrt(discr)) / aij
    else 
        hardsphere_event%type = 0 ! no event will occur
        hardsphere_event%time = bigtime
    end if 
end function hardsphere_event

! methd that predits when a polarized attractice or repulsive
! collision will occur according to a square potential.
! method is called by group object / module
! method stored with circle object / module 
type(event) function polsphere_event (rij, vij, a, b)
    implicit none 
    type(position), intent(in) :: rij 
    type(velocity), intent(in) :: vij 
    type(id), intent(in) :: a, b ! dnlist, uplist partners
    real(kind=dbl) :: aij, bij, cij_1, cij_2, cij_3, cij_4 ! quadratic equation constants wrt to discontinuities 1, 2, 3, and 4
    real(kind=dbl) :: discr1, discr2, discr3, discr4 ! discrimenant wrt cij_1, _2, _3, and _4

    ! determine the event partner
    polsphere_event%dnpart = a
    polsphere_event%uppart = b

    ! calculate quadratic parameters 
    aij = (vij%v(1) ** 2) + (vij%v(2) ** 2)
    bij = (rij%r(1) * vij%v(1)) + (rij%r(2) * vij%v(2))
    cij_1 = (rij%r(1) ** 2) + (rij%r(2) ** 2) - sg1sq ! first discontinuity
    cij_2 = (rij%r(1) ** 2) + (rij%r(2) ** 2) - sg2sq ! second discontinuity
    cij_3 = (rij%r(1) ** 2) + (rij%r(2) ** 2) - sg3sq ! third discontinuity
    cij_4 = (rij%r(1) ** 2) + (rij%r(2) ** 2) - sg4sq ! fourth discontinuity
    discr1 = (bij ** 2) - (aij * cij_1)
    discr2 = (bij ** 2) - (aij * cij_2)
    discr3 = (bij ** 2) - (aij * cij_3)
    discr4 = (bij ** 2) - (aij * cij_4)

    ! predict the next event
    if (bij < 0.0) then ! the centers are approaching
        if (cij_2 < 0.0) then ! if rij is within the first well 
            if (discr1 > 0.0) then ! if the cores will collide
                ! no hard sphere events can occur between polarized charges, as they are located inside the hard sphere
                polsphere_event = reset_event ()
            else ! the cores will miss
                polsphere_event%type = 2 ! an event will take place at sigma2-
                polsphere_event%time = (-bij + sqrt(discr2)) / aij
            end if 
        else if (cij_3 < 0.0) then ! if rij is within the second well
            if (discr2 > 0.0) then ! if the cores will collide
                polsphere_event%type = 3 ! an event will take place at sigma2+
                polsphere_event%time = (-bij - sqrt(discr2)) / aij
            else ! the cores will miss 
                polsphere_event%type = 4 ! an event will take place at sigma3-
                polsphere_event%time = (-bij + sqrt(discr3)) / aij
            end if 
        else if (cij_4 < 0.0) then ! if rij is within the third well
            if (discr3 > 0.0) then ! if the cores will collide
                polsphere_event%type = 5 ! an event will take place at sigma3+
                polsphere_event%time = (-bij - sqrt(discr3)) / aij
            else ! the cores will miss 
                polsphere_event%type = 6 ! an event will take place at sigma4-
                polsphere_event%time = (-bij + sqrt(discr4)) / aij
            end if 
        else ! rij is outside the potential
            if (discr4 > 0.0) then ! if the cores will collide
                polsphere_event%type = 7 ! an event will take place at sigma4+
                polsphere_event%time = (-bij - sqrt(discr4)) / aij
            else ! the outermost cores will miss
                ! no event will take place 
                polsphere_event = reset_event ()
            end if 
        end if 
    else ! the centers are receding
        if (cij_2 < 0.0) then ! if rij is within the first well
            polsphere_event%type = 2 ! an event will take place at sigma2-
            polsphere_event%time = (-bij + sqrt(discr2)) / aij
        else if (cij_3 < 0.0) then ! if rij is within the second well
            polsphere_event%type = 4 ! en event will take place at sigma3-
            polsphere_event%time = (-bij + sqrt(discr3)) / aij
        else if (cij_4 < 0.0) then !if rij is within the third well
            polsphere_event%type = 6 ! an event will take place at sigma4-
            polsphere_event%time = (-bij + sqrt(discr4)) / aij
        else ! rij is outside the potential
            ! no event will take place
            polsphere_event = reset_event ()
        end if 
    end if 
end function polsphere_event


end module active_disk_module

program pol_disk_simulation
    use active_disk_module
    implicit none

    call initialize_system()
    call schedule_ghostcollisions(ghost_event)
    call complete_reschedule() ! move to single_step or initialization

    call record_position_circles(d, sphcooriounit, region)
    call record_position_objects(d, objcooriounit, region)

    do 
        if (single_step()) exit
    enddo

    !call record_position_circles(d, sphcooriounit, region)
    !call record_position_objects(d, objcooriounit, region)


    ! deallocation sequence
    call save()
    call close_files()
end program pol_disk_simulation
