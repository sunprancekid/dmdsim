! AUTHOR : Matthew Dorsey
! DATE : 2022-07-22
! FILENAME : dmdmodule.f90
! PURPOSE : [class description]

module dmdmodule
use dmdconstants
implicit none
integer, parameter :: default_collision_limit = 1e9 ! unless otherwise specified, collision limit is set to this value

type :: dmdnvt
    character(len=40) :: simid
    integer :: mols ! number of particles in dmd simulation 
    real(kind=quad) :: eta, temp ! volume fraction and temperature of dmd simulation
    integer :: num_coll, equil_coll, lim_coll ! number of collisions, number of equilibrium collisions, simulation collision limit
    real(kind=quad) :: time ! elapsed time since start of simultion, in seconds
endtype dmdnvt

contains

! CONSTRUCTORS

! constructs dmd object using command line terminal propmts
subroutine construct_dmdnvt_clt (sim, particle_string)
    implicit none
    type(dmdnvt), intent(out) :: sim
    character(len=*), intent(in) :: particle_string
    integer :: numparticles
    real(kind=quad) :: density
    real(kind=quad) :: temperature

    write(*,'("Enter integer number of ", A, "s to simulate:")') trim(particle_string)
    read(*,*) numparticles
    write(*,'("")')
    write(*,'("Enter simulation desnity as area / volume fraction: ")')
    read(*,*) density 
    write(*,'("")')
    write(*,'("Enter simulation temperature (reduced): ")')
    read(*,*) temperature
    write(*,'("")')

    call initialize_dmdnvt (sim, numparticles, density, temperature)
end subroutine construct_dmdnvt_clt

! constructs dmd object using information stored in file
subroutine construct_dmdnvt (sim, filename)
    implicit none 
    type(dmdnvt), intent(in) :: sim 
    character(len=*), intent(in) :: filename
end subroutine construct_dmdnvt

subroutine initialize_dmdnvt (sim, mols, eta, temp)
    type(dmdnvt), intent(out) :: sim 
    integer, intent(in) :: mols 
    real(kind=quad), intent(in) :: eta, temp
    sim%mols = mols 
    sim%eta = eta
    sim%temp = temp
    sim%num_coll = 0
    sim%lim_coll = default_collision_limit
    sim%time = 0.
end subroutine initialize_dmdnvt


! I / O

subroutine report_status_dmdnvt (sim)
    implicit none 
    type(dmdnvt), intent(in) :: sim

    write(*,'("NVT SIMULATION: ID = ", A)') sim%simid
    write(*,'("NVT SIMULATION: N = ", I8)') sim%mols
    write(*,'("NVT SIMULATION: V = ", F4.3,"% of total")') sim%eta
    write(*,'("NVT SIMULATION: T = ", F4.2)') sim%temp 

! save state of dmd simulation to filename
function save_dmdnvt (sim, filename) result (iobool)
    implicit none 
    type(dmdnvt), intent(in) :: sim
    character(len=*), intent(in) :: filename ! location to save file to
    logical :: iobool ! indicated whether the save operation was successful
    integer :: ioerr ! used to check the success of file opening operation

    iobool = .false.
    open (unit = saveiounit, file = trim(filename), status = 'REPLACE', action = 'WRITE', iostat = ioerr)
    if (ioerr == 0) then ! then write to file 
        write(saveiounit, *) sim%simid
        write(saveiounit, *) sim%mols
        write(saveiounit, *) sim%eta 
        write(saveiounit, *) sim%temp 
        write(saveiounit, *) sim%num_coll
        write(saveiounit, *) sim%equil_coll
        write(saveiounit, *) sim%lim_coll 
        write(saveiounit, *) sim%time 
        close (unit = saveiounit, status = 'KEEP')
        iobool = .true.
    endif
    ! if the operation failed, inform the user 
    if (.not. iobool) write(*,*) 'Unable to save dmd state to ', trim(filename), '.'
end function save_dmdnvt

! load state of dmd simulation to filename 
function load_dmdnvt (sim, filename) result (iobool)
    implicit none 
    type(dmdnvt), intent(out) :: sim 
    character(len=*), intent(in) :: filename 
    logical :: iobool ! indicates whether the load operation was successful
    integer :: ioerr ! used to determine whether file opening operation was successful

    iobool = .false.
    open (unit = saveiounit, file = trim(filename), status = 'OLD', action = 'READ', iostat = ioerr)
    if (ioerr == 0) then ! read from save file
        read(saveiounit, *) sim%simid
        read(saveiounit, *) sim%mols
        read(saveiounit, *) sim%eta 
        read(saveiounit, *) sim%temp
        read(saveiounit, *) sim%num_coll
        read(saveiounit, *) sim%equil_coll
        read(saveiounit, *) sim%lim_coll
        read(saveiounit, *) sim%time
        close(unit = saveiounit, status = 'KEEP')
        iobool = .true.
    endif
    ! if the operation failed, inform the user
    if (.not. iobool) write(*,*) 'Unable to load dmd state from ', trim(filename), '.'
end function load_dmdnvt

! GETTERS AND SETTERS

! returns the number of particles in the simulation
function get_dmd_particles (sim) result (n)
    implicit none
    type(dmdnvt), intent(in) :: sim 
    integer :: n 
    n = sim%mols
end function get_dmd_particles

! returns simulation desnity as volume / area fraction
function get_dmd_density (sim) result (density)
    implicit none 
    type(dmdnvt), intent(in) :: sim 
    real(kind=quad) :: density 
    density = sim%eta 
end function get_dmd_density

! sets maximum number of collisions for dmd simulation to value passed as parameter
! values passed to subroutine that are equal to zero or less are ignored
subroutine set_dmd_collision_limit (sim, limit)
    implicit none
    type(dmdnvt), intent(inout) :: sim
    integer, intent(in) :: limit
    if (limit > 0) sim%lim_coll = limit
end subroutine set_dmd_collision_limit

! returns the number of collisions that have occured in the simulation so far
function get_dmd_collisions (sim) result (num_coll)
    implicit none
    type(dmdnvt), intent(in) :: sim
    integer :: num_coll
    num_coll = sim%num_coll
end function get_dmd_collisions



! LOAD AND SAVE FUNCTIONS



! HELPER METHODS

! logical function that determines if the maximum number of events in the simulation have occurred
function simulation_limit_reached (sim) result (bool)
    implicit none
    type(dmdnvt), intent(in) :: sim
    logical :: bool 
    bool = sim%num_coll >= sim%lim_coll
end function simulation_limit_reached

! increments the number of collisions in the simulation by one
subroutine increment_collisions (sim)
    implicit none
    type(dmdnvt), intent(inout) :: sim
    sim%num_coll = sim%num_coll + 1
end subroutine increment_collisions

end module dmdmodule