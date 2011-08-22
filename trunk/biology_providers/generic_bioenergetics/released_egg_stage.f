cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Generic egg stage template. 
c
c     Provides module released_egg_stage and class released_egg
c     for embedding in stage manager module (do not provide
c     independent particle_state interface)
c     ---------------------------------------------------
c     $Rev: $
c     $LastChangedDate: $
c     $LastChangedBy:  $ 
c
c     TODO:
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module released_egg_stage  

      use egg_properties        ! import egg_devel_rate
      use time_tools            ! import clock type
      use particle_tracking     ! space types/methods 
      use physical_fields
      use run_context, only: ctrlfile => simulation_file
      use input_parser
      use output                ! access polytype for get_prop_state/get_metadata_state


      implicit none            
      private     

c     ----------------------------------------------------------------    
c     Class released_egg is for the embedded egg stage interface 
c     and contains attributes and state variables relating to the 
c     egg stage alone
c     ----------------------------------------------------------------

      type released_egg
      private
         real :: completion   ! 0 < completion < 1. Hatch at 1
      end type
      public :: released_egg ! make type state_attributes visible outside


c.....Particle state interface: (stand-alone)
      public :: init_released_egg_stage  ! module operator
      public :: close_released_egg_stage ! module operator
      public :: init_state_attributes_egg
      public :: get_active_velocity_egg
      public :: update_particle_state_egg
      public :: delete_state_attributes_egg
c      public :: write_state_attributes
c      public :: get_particle_version ! hide this, because it is only for embedded usage
c      public :: get_property         ! currently void
c      public :: get_metadata_state   ! currently void


      

c.....All other subroutine than these below are private to the module

      interface get_property
      module procedure get_prop_state
      end interface
    
       
c     ==============================================================
c     ========             module data section              ========
c     ==============================================================
      integer, parameter       :: verbose = 0 ! control output volume
      character*(*), parameter :: title   = "egg stage model"
c     ===============================================================
                               contains
c     ===============================================================
      
      
      subroutine init_released_egg_stage () ! module operator
c     ---------------------------------------------------------------
c     This subroutine should initialize this module 
c     ---------------------------------------------------------------
      write (*,*) "initializing ", trim(get_particle_version())

      end subroutine init_released_egg_stage 



      character*100 function get_particle_version()  
c     ---------------------------------------------------------------
c     ---------------------------------------------------------------
      get_particle_version = title // ": $Rev: $"
      end function



      subroutine close_released_egg_stage () ! module operator
c     ---------------------------------------------------------------
c     Module clean up stuff here
c     (e.g. memory deallocation, MPI finalisation, final dumps etc.)
c     Only deallocate own data, i.e. data allocated in this module
c     ---------------------------------------------------------------
      write(*,*) "close_released_egg_stage ():"
      end subroutine close_released_egg_stage  
      


      subroutine init_state_attributes_egg(state, space, time_dir,             
     +                                 initdata)
c     ---------------------------------------------------------------
c     This subroutine should initialize state attributes (and 
c     possibly also space attributes like boundary conditions and mobility)
c
c     Their space part shas already been initialized so that e.g. 
c     local physical conditions
c     can be assessed using their position. 
c    
c     Currently initdata is unused
c     --------------------------------------------------------------- 
      type(released_egg),intent(out)        :: state
      type(spatial_attributes),intent(inout) :: space
      real,intent(in)                        :: time_dir            
      character*(*),intent(in)               :: initdata
c     ---------------------------------------------------------------
      if (time_dir<0) then
         write(*,*) "init_state_attributes: time_dir<0 not implemented" 
         stop
      endif
      
c     parse initdata string here, if needed

      call set_egg_spawned(state, space) ! currently only support this init state
  
c     ----------------------------------------------------
      end subroutine init_state_attributes_egg

      


      subroutine get_active_velocity_egg(state, space, v_active)
c     --------------------------------------------------------------
c     Currently no migration pattern implemented (set v_active = 0 0 0)
c     If available, this subroutine should determine
c     v_active from (state, space) of the particle
c     --------------------------------------------------------------
      type(released_egg), intent(in)       :: state
      type(spatial_attributes), intent(in) :: space      
      real, intent(out)                    :: v_active(:) ! meter/second
c 
      v_active = 0.
c     --------------------------------------------------------------
      end subroutine get_active_velocity_egg


      subroutine update_particle_state_egg(state, space, dt, 
     +                       mortality_rate, die, next)
c     --------------------------------------------------------------
c     This subroutine should integrate state forward for
c     time period dt.
c     Also set auxillary parameters:
c        mortality_rate: current net mortality rate [1/sec]
c        die           : if .true. kill this entity (avoid having mortality_rate=infinity)
c        next          : if .true. advance (or regress) to ontogenetic stage (depending on sign of dt)
c     The subroutines below call a potentially deep hierachy
c     of other subroutines; to avoid duplicated interpolations of
c     the local environment, a flat container structure local_environment
c     is created which contains potential look-ups, and this 
c     container structure is passed around to increase flexibility
c     and avoid messy argument relaying 
c     --------------------------------------------------------------
      type(released_egg), intent(inout)       :: state
      type(spatial_attributes), intent(inout) :: space
      real,intent(in)                         :: dt ! in seconds
      real,intent(out)                        :: mortality_rate
      logical,intent(out)                     :: die, next   
c
      real                                    :: dcomp_dt, xyz(3),temp
      integer                                 :: istat
      real,parameter                          :: epsil = 1.e-4 ! avoid flagging if at stage boundary
c     --------------------------------------------------------------
      if (dt<0) stop "negative dt values currently blocked"

      call get_tracer_position(space, xyz)
      call interpolate_temp(xyz,temp,istat)
      call egg_devel_rate(temp, dcomp_dt)
      
      state%completion = state%completion + dcomp_dt*dt  ! OK for dt<0
c
c     --- set mortality_rate + die ---
c              
      call evaluate_survival_chance(state, space, dt,
     +                                    mortality_rate, die)
c
c     --- set stage promotion/regression request (next) ---
c   
      if ((state%completion > 1.0+epsil).or.
     +     state%completion < -epsil) then
         next = .true.
      else
         next = .false.
      endif
c     --------------------------------------------------------------
      end subroutine update_particle_state_egg


      subroutine evaluate_survival_chance(state, space, dt,
     +                                    mortality_rate, die)
c     ---------------------------------------------------
c     Evaluate mortality_rate and death request (die) for
c     this egg stage corresponding to time interval dt
c     ---------------------------------------------------
      type(released_egg),intent(inout)       :: state
      type(spatial_attributes),intent(inout) :: space
      real,intent(in)                        :: dt             ! [sec]
      real,intent(out)                       :: mortality_rate ! [1/sec]
      logical,intent(out)                    :: die
c     ---------------------------------------------------
      mortality_rate = 0.0 ! elaborate later
      die            = .false.
      end subroutine evaluate_survival_chance


      subroutine delete_state_attributes_egg(state) 
      type(released_egg), intent(inout) :: state 
      end subroutine delete_state_attributes_egg


      subroutine write_state_attributes(state)
c     --------------------------------------------------------------
c     log output of state
c     --------------------------------------------------------------
      type(released_egg), intent(in) :: state 
      end subroutine 

c     ==============================================================
c     ========      internal (non public) subroutines       ========
c     ==============================================================

      subroutine set_egg_spawned(state, space) 
c     ---------------------------------------------------
c     Initialize egg in newly spawned state
c     make no assumptions about previous history
c     ---------------------------------------------------
      type(released_egg), intent(inout)      :: state
      type(spatial_attributes),intent(inout) :: space
c     ---------------------------------------------------
      call set_tracer_mobility_free(space)     
      state%completion = 0.0         
      end subroutine set_egg_spawned

c     ==============================================================
c     output interface - to be filled in later
c     ==============================================================

      subroutine get_prop_state(state,var,bucket,status)
c------------------------------------------------------------  
      type(released_egg),intent(in) :: state
      type(variable),intent(in)         :: var
      type(polytype), intent(out)       :: bucket
      integer, intent(out)              :: status
c------------------------------------------------------------  
      end subroutine


      subroutine get_metadata_state(var_name,var,status)
c------------------------------------------------------------  
      character(*), intent(in)   :: var_name
      type(variable),intent(out) :: var
      integer, intent(out)       :: status
c------------------------------------------------------------  
      end subroutine get_metadata_state

      end module
