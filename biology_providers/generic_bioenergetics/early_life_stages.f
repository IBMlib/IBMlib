cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Multi stage manager for simulating early life stages
c     with several ontogenetic stages 
c     ---------------------------------------------------
c     $Rev: $
c     $LastChangedDate: $
c     $LastChangedBy:  $ 
c
c     This module assembles several stages expressed as classes
c     imported from other modules into a multi stage container class and
c     implements the particle_state interface for this assembly class
c     and manages transitions between different ontogenetic stages
c     This module allows both forward/backward time simulation
c
c     TODO: 
c          consider meaning and update of state%survival for dt<0 
c          more flexible init for internal stage shift
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module particle_state

      use time_tools            ! import clock type
      use particle_tracking     ! space types/methods 
      use physical_fields
      use run_context, only: ctrlfile => simulation_file
      use input_parser
      use output                ! access polytype for get_prop_state/get_metadata_state
      
c
c     import with ONLY modifier, because modules released_egg_stage,
c     yolksac_larval_stage, feeding_larval_stage may provide full 
c     particle_state interface themself
c
      use released_egg_stage,    ONLY: released_egg,
     +                                 init_released_egg_stage,
     +                                 close_released_egg_stage,
     +                                 init_state_attributes_egg,
     +                                 get_active_velocity_egg,
     +                                 update_particle_state_egg,
     +                                 delete_state_attributes_egg
      
      use yolksac_larval_stage,  ONLY: yolksac_larvae,
     +                                 init_yolksac_larval_stage,
     +                                 close_yolksac_larval_stage,
     +                                 init_state_attributes_YL,
     +                                 get_active_velocity_YL,
     +                                 update_particle_state_YL,
     +                                 delete_state_attributes_YL    

      use feeding_larval_stage,  ONLY: feeding_larvae,
     +                                 init_feeding_larval_stage,
     +                                 close_feeding_larval_stage,
     +                                 init_state_attributes_FL,
     +                                 get_active_velocity_FL,
     +                                 update_particle_state_FL,
     +                                 delete_state_attributes_FL
      
      implicit none            
      private     

c     -----------------------------------------------    
c     Assembly class of different ontogenetic stages
c
c     Class state_attributes contains all attributes 
c     related to the state of a particle beyond 
c     spatial aspects
c     -----------------------------------------------
      type state_attributes
      private

      type(released_egg)   :: egg      ! active_stage == 1
      type(yolksac_larvae) :: ysaclarv ! active_stage == 2
      type(feeding_larvae) :: feedlarv ! active_stage == 3 
      integer              :: active_stage

      real        :: survival   ! 0 < survival < 1 of this particle 
      logical     :: alive                
      integer     :: death_day  ! Julian day
      
      integer     :: orig_boxID ! spatial release box
      integer     :: particleID ! ensemble ID - negative is unset
      
      end type
      public :: state_attributes ! make type state_attributes visible outside


c.....All other subroutine than these below are private to the module
      public :: init_particle_state ! module operator
      public :: close_particle_state ! module operator
      public :: init_state_attributes
      public :: get_active_velocity
      public :: update_particle_state
      public :: delete_state_attributes 
      public :: write_state_attributes
      public :: get_particle_version
      public :: get_property    ! currently void
      public :: get_metadata_state ! currently void


      interface get_property
      module procedure get_prop_state
      end interface
      
c     ------ internal interfaces ------

      interface init_state_attributes_
      module procedure init_state_attributes_egg   
      module procedure init_state_attributes_YL   
      module procedure init_state_attributes_FL   
      end interface 

      interface get_active_velocity_
      module procedure get_active_velocity_egg  
      module procedure get_active_velocity_YL   
      module procedure get_active_velocity_FL  
      end interface 

      interface update_particle_state_
      module procedure update_particle_state_egg  
      module procedure update_particle_state_YL  
      module procedure update_particle_state_FL  
      end interface 

      interface delete_state_attributes_
      module procedure delete_state_attributes_egg  
      module procedure delete_state_attributes_YL     
      module procedure delete_state_attributes_FL 
      end interface 

c     ==============================================================
c     ========             module data section              ========
c     ==============================================================

      integer, parameter    :: verbose = 0      ! control output volume
      integer               :: particle_counter ! module counter for generating IDs for particles
      integer, parameter    :: nstages = 3      ! in class state_attributes
      
c     ===============================================================
                               contains
c     ===============================================================
      
      
      subroutine init_particle_state() ! module operator
c     ---------------------------------------------------------------
c     This subroutine should initialize this module 
c     ---------------------------------------------------------------
      write (*,*) trim(get_particle_version())
c     
c     --- pass init to sub modules
c     
      call init_released_egg_stage()    ! module init        
      call init_yolksac_larval_stage()  ! module init   
      call init_feeding_larval_stage () ! module init  

      particle_counter = 1      ! ID for next larvae   
    
cc      ----- test section -----

      end subroutine 


      character*100 function get_particle_version()  
c     ---------------------------------------------------------------
c     ---------------------------------------------------------------
      get_particle_version = "multi stage manager  : $Rev: $"
      end function


      subroutine close_particle_state() ! module operator
c     ---------------------------------------------------------------
c     Module clean up stuff here
c     (e.g. memory deallocation, MPI finalisation, final dumps etc.)
c     Only deallocate own data, i.e. data allocated in this module
c     ---------------------------------------------------------------
      write(*,*) "close_particle_state():"
      particle_counter               = 0
c
c     --- pass init to sub modules reversely
c
      call close_released_egg_stage()    ! module close down    
      call close_yolksac_larval_stage()  ! module close down   
      call close_feeding_larval_stage () ! module close down
c      
      end subroutine 
      


      subroutine init_state_attributes(state, space, time_dir,             
     +                                 initdata, emitboxID)
c     ---------------------------------------------------------------
c     This subroutine should initialize state attributes (and 
c     possibly also space attributes like boundary conditions and mobility)
c
c     Their space part shas already been initialized so that e.g. 
c     local physical conditions can be assessed using their position. 
c    
c     The first 3 letters of initdata indicate in which stage
c     the particle is initialized; these three letters are chopped
c     off initdata, before initdata it is passed on to stage
c     initialization. The following tags apply for picking 
c     initial stages
c        egg_stage
c        yolksac_larvae
c        feeding_larvae) :: feedlarv ! active_stage == 3 
c     --------------------------------------------------------------- 
      type(state_attributes),intent(out)     :: state
      type(spatial_attributes),intent(inout) :: space
      real,intent(in)                        :: time_dir            
      character*(*),intent(in)               :: initdata
      integer,intent(in)                     :: emitboxID
c
      character*3                            :: stage_tag
      character(len=max(1,len(initdata)))    :: rest
      integer                                :: leni
c     ---------------------------------------------------------------
      
c      
c     --- parse initdata(1:3) here to resolve which stage to initialize ---
c
      leni = len(initdata)
      if (leni < 3) then
         stop "init_state_attributes: no stage specified"
      else
         stage_tag = initdata(1:3)
         if (leni > 3) then
            rest = initdata(4:)
         else
            rest = ""
         endif
      endif
c      
c     --- activate requested sub stage ---
c     
      if     (stage_tag == "egg") then
         call init_state_attributes_(state%egg, space, time_dir, rest)
         state%active_stage = 1
      elseif (stage_tag == "ylv") then
         call init_state_attributes_(state%ysaclarv,space,time_dir,rest)
         state%active_stage = 2
      elseif (stage_tag == "flv") then  
         call init_state_attributes_(state%feedlarv,space,time_dir,rest)
         state%active_stage = 3
      else
         write(*,*) "init_state_attributes: unsupported stage tag ", 
     +               stage_tag
         stop 
      endif
c
c     --- set particle administrative sttributes ---
c
      state%orig_boxID = emitboxID   ! store releasing box
      state%particleID = particle_counter    
      particle_counter = particle_counter + 1 ! module data
      state%survival   = 1.0   
      state%alive      = .true.
c     ----------------------------------------------------
      end subroutine 



      subroutine get_active_velocity(state, space, v_active)
c     --------------------------------------------------------------
c     --------------------------------------------------------------
      type(state_attributes), intent(in)   :: state
      type(spatial_attributes), intent(in) :: space      
      real, intent(out)                    :: v_active(:) ! meter/second
c     --------------------------------------------------------------
      if      (state%active_stage == 1) then
         call get_active_velocity_(state%egg,      space, v_active)
      elseif  (state%active_stage == 2) then
         call get_active_velocity_(state%ysaclarv, space, v_active)  
      elseif  (state%active_stage == 3) then   
         call get_active_velocity_(state%feedlarv, space, v_active)
      else
         stop "get_active_velocity: unexpected"
      endif
c     --------------------------------------------------------------
      end subroutine get_active_velocity



      subroutine update_particle_state(state, space, dt)
c     --------------------------------------------------------------
c     This subroutine should integrate state forward for
c     time period dt 
c
c     The subroutines below call a potentially deep hierachy
c     of other subroutines; to avoid duplicated interpolations of
c     the local environment, a flat container structure local_environment
c     is created which contains potential look-ups, and this 
c     container structure is passed around to increase flexibility
c     and avoid messy argument relaying 
c
c     When changing stage, keep old data, becauuse it may contain 
c     needed life history information
c
c     TODO:
c        proper inits for dt<0 when changing stages
c     --------------------------------------------------------------
      type(state_attributes), intent(inout)   :: state
      type(spatial_attributes), intent(inout) :: space
      real,intent(in)                         :: dt ! in seconds
c
      real                                    :: mortality_rate
      logical                                 :: die, next
c     --------------------------------------------------------------
      if (dt<0) stop "dt<0 not released"
      if (.not.state%alive) return
c
c     --- fanout update to active sub stage ---
c
      if      (state%active_stage == 1) then
         call update_particle_state_(state%egg, space, dt, 
     +                       mortality_rate, die, next)
      elseif  (state%active_stage == 2) then
         call update_particle_state_(state%ysaclarv, space, dt, 
     +                       mortality_rate, die, next)
      elseif  (state%active_stage == 3) then   
         call update_particle_state_(state%feedlarv, space, dt, 
     +                       mortality_rate, die, next)
      else
         stop "get_active_velocity: unexpected"
      endif
c
c     --- update survival/death/stage attributes ---
c
      state%survival = state%survival * exp(-mortality_rate*dt)
      if (die) call set_state_dead(state, space)  

      if (.not.state%alive) return ! do not change stage if particle is dead
      
      if (next .eqv. .true.) then
         if     (dt > 0) then
            state%active_stage = state%active_stage + 1
         elseif (dt < 0) then
            state%active_stage = state%active_stage - 1
         endif

c
c    Here init is with initdata = "". We need to differentiate this and to
c    cover both forward/backward time cases (which should be initialized differently)
c    apply time_dir=dt 
         if      (state%active_stage == 1) then
            call init_state_attributes_(state%egg,space,      dt, "")            
         elseif  (state%active_stage == 2) then
            call init_state_attributes_(state%ysaclarv,space, dt, "")  
         elseif  (state%active_stage == 3) then   
            call init_state_attributes_(state%feedlarv,SPACE, DT, "")  
         else
            stop "update_particle_state: unknown stage"
         endif

      endif    ! if (next == .true.)
  
c     --------------------------------------------------------------
      end subroutine 



      subroutine delete_state_attributes(state) 
      type(state_attributes), intent(inout) :: state 
      
      call delete_state_attributes_(state%egg)      ! has potentially allocated data
      call delete_state_attributes_(state%ysaclarv) ! has potentially allocated data
      call delete_state_attributes_(state%feedlarv) ! has potentially allocated data

      end subroutine 


      subroutine write_state_attributes(state)
c     --------------------------------------------------------------
c     log output of state
c     --------------------------------------------------------------
      type(state_attributes), intent(in) :: state 
      end subroutine 

c     ==============================================================
c     ========      internal (non public) subroutine        ========
c     ==============================================================


      subroutine set_state_dead(state, space)
c     ---------------------------------------------------
c     Set particle state to dead
c     leave variables at value just before larvae died
c     make no assumptions about previous history
c     ---------------------------------------------------
      type(state_attributes), intent(inout)  :: state
      type(spatial_attributes),intent(inout) :: space
      type(clock), pointer                   :: current_time
c     ---------------------------------------------------
      call set_tracer_mobility_stop(space)           
      state%survival   = 0.0    
      state%alive      = .false.
      current_time     => get_master_clock()
      call get_julian_day(current_time, state%death_day)
      end subroutine set_state_dead
  
c     ==============================================================
c     output interface - to be filled in later
c     ==============================================================

      subroutine get_prop_state(state,var,bucket,status)
c------------------------------------------------------------  
      type(state_attributes),intent(in) :: state
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
