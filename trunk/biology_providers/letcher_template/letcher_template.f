cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Letcher template 
c     ---------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
c     Preliminary skeleton for Letcher type larval bioenergetics
c     with daily growth increments
c
c     Phoenix version of IBMlib
c
c     asc:July2010 - made rough skeleton for Letcher model
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module particle_state

      use time_tools           ! import clock type
      use particle_tracking    ! space types/methods 
      use physical_fields
      use run_context, only: ctrlfile => simulation_file
      use input_parser

      implicit none            
      private     

c     -----------------------------------------------    
c     Class state_attributes contains all attributes 
c     related to the state of a particle beyond 
c     spatial aspects
c     -----------------------------------------------
      type state_attributes
      private
      real        :: length               ! length in mm
      real        :: weight               
      real        :: survival             ! 0 < survival < 1 of this larvae
      real        :: max_weight           ! previously max weight of this larvae

      logical     :: alive
      integer     :: hatch_day            ! Julian day
      integer     :: death_day            ! Julian day

      integer     :: orig_boxID ! spatial release box
      integer     :: larvID     ! ensemble ID - negative is unset
      
      end type
      public :: state_attributes     ! make type state_attributes visible outside


c.....All other subroutine than these below are private to the module
      public :: init_particle_state  ! module operator
      public :: close_particle_state ! module operator
      public :: init_state_attributes
      public :: get_active_velocity
      public :: update_particle_state
      public :: delete_state_attributes 
      public :: write_state_attributes
      public :: get_particle_version


c     ==============================================================
c     ========             module data section              ========
c     ==============================================================
      integer, parameter    :: verbose = 0    ! control output volume

      integer               :: larv_counter   ! for generating IDs for larvae

      
c     here one should define data and parameters relating to the
c     species, not the individual, e.g hatch length

c     ===============================================================
                                  contains
c     ===============================================================
     
 
      subroutine init_particle_state()  ! module operator
c     ---------------------------------------------------------------
c     This subroutine should initialize this module 
c     ---------------------------------------------------------------
      write (*,*) trim(get_particle_version())
      larv_counter = 1 ! ID for next larvae

c     you may read parameters like this:
c
c     call read_control_data(ctrlfile,"my_parameter_name", value)
      end subroutine 

      character*100 function get_particle_version()  
      get_particle_version = "Letcher biology provider : $Rev$"
      end function


      subroutine close_particle_state() ! module operator
c     ---------------------------------------------------------------
c     Module clean up stuff here
c     (e.g. memory deallocation, MPI finalisation, final dumps etc.)
c     Only deallocate own data, i.e. data allocated in this module
c     ---------------------------------------------------------------
      write(*,*) "close_particle_state():"
      larv_counter               = 0
      end subroutine 
      


      subroutine init_state_attributes(state, space, time_dir,             
     +                                 initdata, emitboxID)
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
      type(state_attributes),intent(out)     :: state
      type(spatial_attributes),intent(inout) :: space
      real,intent(in)                        :: time_dir            
      character*(*),intent(in)               :: initdata
      integer,intent(in)                     :: emitboxID
c     ---------------------------------------------------------------
      if (time_dir<0) then
         write(*,*) "init_state_attributes: time_dir<0 not implemented" 
         stop
      endif
      
c     parse initdata string here, if needed

      state%orig_boxID = emitboxID   ! store releasing box
      state%larvID     = larv_counter    
      larv_counter     = larv_counter + 1 ! module data
      call set_state_hatched(state, space)
c     ----------------------------------------------------
      end subroutine 


      subroutine get_active_velocity(state, space, v_active)
c     --------------------------------------------------------------
c     Currently no migration pattern implemented (set v_active = 0 0 0)
c     If available, this subroutine should determine
c     v_active from (state, space) of the particle

c     --------------------------------------------------------------
      type(state_attributes), intent(in)   :: state
      type(spatial_attributes), intent(in) :: space      
      real, intent(out)                    :: v_active(:) ! meter/second
c 
      v_active = 0.
c     --------------------------------------------------------------
      end subroutine 



      subroutine update_particle_state(state, space, dt)
c     --------------------------------------------------------------
c     This subroutine should integrate state forward for
c     time period dt (need not be small)
c     Currently only dt == 1 day is accepted
c     --------------------------------------------------------------
      type(state_attributes), intent(inout)   :: state
      type(spatial_attributes), intent(inout) :: space
      real,intent(in)                         :: dt ! in seconds
      real                                    :: xyz(3),zbiomass(1)
      integer                                 :: status
c     --------------------------------------------------------------
      
c     --------------------------------------------------------------
      if (dt<0) stop "negative dt values currently blocked"
      if (abs(dt-86400.) > 1.e-3) then
         write(*,*) "update_particle_state: currently only daily steps" 
         stop 
      endif
      
      call get_tracer_position(space, xyz)
      call interpolate_zooplankton (xyz, zbiomass, status) 
c      assess food availability (zbiomass -> predens)

      call calculate_ingestion(predens, ingestion)
      call grow_larvae(state,ingestion)
      call update_survival_chance(state, space)
c     --------------------------------------------------------------
      end subroutine 



      subroutine delete_state_attributes(state) 
      type(state_attributes), intent(inout) :: state 
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


      subroutine length_to_weight(length, weight)
c     ---------------------------------------------------
c     Convert larval length to larval weight
c     ---------------------------------------------------
      real, intent(in)  :: length
      real, intent(out) :: weight
      weight = ...
      end subroutine length_to_weight


      subroutine weight_to_length(weight, length)
c     ---------------------------------------------------
c     Convert weight to length 
c     ---------------------------------------------------
      real, intent(in)  :: weight
      real, intent(out) :: length
      length = ...
      end subroutine length_to_weight


      subroutine set_state_hatched(state, space) 
c     ---------------------------------------------------
c     Initialize larvae in newly hatched state
c     make no assumptions about previous history
c     ---------------------------------------------------
      call set_tracer_mobility_free(space) 
      state%length     =           
      state%weight     =                       
      state%survival   = 1.0 ! fresh start         
      state%max_weight = state%weight
      state%alive      = .true.
      call get_julian_day(current_time, state%hatch_day)
      end subroutine set_state_hatched



      subroutine set_state_dead(state, space)
c     ---------------------------------------------------
c     Set larval state to dead
c     leave most variable at value just before larvae died
c     make no assumptions about previous history
c     ---------------------------------------------------
      call set_tracer_mobility_stop(space)           
      state%survival   = 0.0    
      state%alive      = .false.
      call get_julian_day(current_time, state%death_day)
      end subroutine set_state_dead



      subroutine calculate_ingestion(predens, ingestion)
c     ---------------------------------------------------
c     Letcher Eqs 3-13
c     ---------------------------------------------------
      compute encounter rates for each prey class
      compute capture success for each prey class
      compute handling time   for each prey class
      rank prey -> diet composition -> ingestion
      ingestion = max(ingestion, Cmax)
      end subroutine calculate_ingestion


      subroutine grow_larvae(state,ingestion)
c     ---------------------------------------------------
c     Letcher Eqs 14-17
c     ---------------------------------------------------
      call calculate_TC(...)
      compute AE
      growth = ingestion*AE - TC
      update state%weight 
      update state%length

      state%max_weight = max(state%max_weight, state%weight)

      end subroutine grow_larvae



      subroutine calculate_TC(...)
c     ---------------------------------------------------
c     Letcher Eqs 16+17
c     ---------------------------------------------------
      end subroutine calculate_TC




      subroutine update_survival_chance(state, space)
c     ---------------------------------------------------
c     Update the survivalk counter for this larvae 
c     corresponding to life this day 
c     Letcher Eqs 18-22
c     ---------------------------------------------------

      calculate starvation survival   
      call calculate_predation_survival(state,space,pred_surv)  
      
      tot_surv = starv_surv*pred_surv            ! for this day
      state%survival = state%survival * tot_surv ! for life time
      
      if (died) call set_state_dead(state, space) 

      end subroutine update_survival_chance

      

      subroutine calculate_predation_survival(state,space,pred_surv)
c     --------------------------------------------------- 
c     Letcher Eqs 19-22
c     At some point we may also overlay predator maps (space available)
c     ---------------------------------------------------
      end subroutine calculate_predation_survival
      



      end module
