cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Based on generic_bioenergetics/yolksac_larval_stage.f
c     MEECE specific ecotox parameterization for sprat yolksac larvae
c
c     Provides module yolksac_larval_stage and class yolksac_larvae
c     for embedding in stage manager module (do not provide
c     independent particle_state interface)
c
c     hatch length is set from ambient conditions (assumed constant during egg development)
c     ---------------------------------------------------
c     $Rev: $
c     $LastChangedDate: $
c     $LastChangedBy:  $ 
c
c     TODO:
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module yolksac_larval_stage  

      use time_tools            ! import clock type
      use particle_tracking     ! space types/methods 
      use physical_fields
      use run_context, only: ctrlfile => simulation_file
      use input_parser
      use output                ! access polytype for get_prop_state/get_metadata_state


      implicit none            
      private     

c     ----------------------------------------------------------------    
c     Class yolksac_larvae is for the embedded stage interface 
c     and contains attributes and state variables relating to the 
c     yolksac stage alone
c     ----------------------------------------------------------------

      type yolksac_larvae
      private
         real :: completion   ! 0 < completion < 1. Hatch at 1
         real :: hatch_length ! before yolk absorbtion
      end type
      public :: yolksac_larvae ! make type state_attributes visible outside


c.....Particle sub state interface: (for multi stage manager)
      public :: init_yolksac_larval_stage  ! module operator
      public :: close_yolksac_larval_stage ! module operator
      public :: init_state_attributes_YL
      public :: get_active_velocity_YL
      public :: update_particle_state_YL
      public :: delete_state_attributes_YL
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
      character*(*), parameter :: title   = "yolksac_larval_stage model"

      real,parameter :: devel_prefac = 0.27/100./86400. ! daewel etal 2008
      real,parameter :: devel_expo   = 1.495            ! daewel etal 2008

c     ---- ecotox parameters for hatch length ----
      real,parameter :: L0_hatch   = 7.0 ! mm
      real,parameter :: Linf_hatch = 5.0 ! mm
      real,parameter :: K1_prefac  = 0.00049 ! yg/l
      real,parameter :: K1_expo    = 5.95  !

c     ===============================================================
                               contains
c     ===============================================================
      
      
      subroutine init_yolksac_larval_stage () ! module operator
c     ---------------------------------------------------------------
c     This subroutine should initialize this module 
c     ---------------------------------------------------------------
      write (*,*) "initializing ", trim(get_particle_version())

      end subroutine init_yolksac_larval_stage



      character*100 function get_particle_version()  
c     ---------------------------------------------------------------
c     ---------------------------------------------------------------
      get_particle_version = title // ": $Rev: $"
      end function



      subroutine close_yolksac_larval_stage () ! module operator
c     ---------------------------------------------------------------
c     Module clean up stuff here
c     (e.g. memory deallocation, MPI finalisation, final dumps etc.)
c     Only deallocate own data, i.e. data allocated in this module
c     ---------------------------------------------------------------
      write(*,*) "close_yolksac_larval_stage ():"
      end subroutine close_yolksac_larval_stage 
      


      subroutine init_state_attributes_YL(state, space, time_dir,             
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
      type(yolksac_larvae),intent(out)       :: state
      type(spatial_attributes),intent(inout) :: space
      real,intent(in)                        :: time_dir            
      character*(*),intent(in)               :: initdata
c     ---------------------------------------------------------------
      if (time_dir<0) then
         write(*,*) "init_state_attributes: time_dir<0 not implemented" 
         stop
      endif
      
c     
      if (len(trim(adjustl(initdata))) == 0) then
         call set_larvae_hatched(state, space) ! currently only support this init state
      else
         write(*,*) "init_state_attributes_YL: unknown init request="//
     +              trim(adjustl(initdata))
      endif
c     ----------------------------------------------------
      end subroutine init_state_attributes_YL

      


      subroutine get_active_velocity_YL(state, space, v_active)
c     --------------------------------------------------------------
c     Currently no migration pattern implemented (set v_active = 0 0 0)
c     If available, this subroutine should determine
c     v_active from (state, space) of the particle
c     --------------------------------------------------------------
      type(yolksac_larvae), intent(in)       :: state
      type(spatial_attributes), intent(in) :: space      
      real, intent(out)                    :: v_active(:) ! meter/second
c 
      v_active = 0.
c     --------------------------------------------------------------
      end subroutine get_active_velocity_YL


      subroutine update_particle_state_YL(state, space, dt, 
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
      type(yolksac_larvae), intent(inout)       :: state
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

c     inlined biology_providers/sprat/sprat_yolksac_larv.f
      
      dcomp_dt = devel_prefac * temp**devel_expo ! [1/sec]
      state%completion = state%completion + dcomp_dt*dt  ! OK for dt<0
c
c     --- set mortality_rate + die ---
c              
      mortality_rate = 0.0 ! no information
      die            = .false.
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
      end subroutine update_particle_state_YL


      subroutine evaluate_survival_chance(state, space, dt,
     +                                    mortality_rate, die)
c     ---------------------------------------------------
c     Evaluate mortality_rate and death request (die) for
c     this yolk-sac larval corresponding to time interval dt
c     ---------------------------------------------------
      type(yolksac_larvae),intent(inout)       :: state
      type(spatial_attributes),intent(inout) :: space
      real,intent(in)                        :: dt             ! [sec]
      real,intent(out)                       :: mortality_rate ! [1/sec]
      logical,intent(out)                    :: die
c     ---------------------------------------------------
      mortality_rate = 0.0 ! elaborate later
      die            = .false.
      end subroutine evaluate_survival_chance


      subroutine delete_state_attributes_YL(state) 
      type(yolksac_larvae), intent(inout) :: state 
      end subroutine delete_state_attributes_YL


      subroutine write_state_attributes(state)
c     --------------------------------------------------------------
c     log output of state
c     --------------------------------------------------------------
      type(yolksac_larvae), intent(in) :: state 
      end subroutine 

c     ==============================================================
c     ========      internal (non public) subroutines       ========
c     ==============================================================

      subroutine set_larvae_hatched(state, space) 
c     ---------------------------------------------------
c     Initializelarvae in newly hatched state
c     make no assumptions about previous history
c     ---------------------------------------------------
      type(yolksac_larvae), intent(inout)    :: state
      type(spatial_attributes),intent(inout) :: space
      real                                    :: cu,xyz(3),temp, K1
      integer                                 :: istat
c     ---------------------------------------------------
      call set_tracer_mobility_free(space)     
      state%completion   = 0.0    
c     --- set hatch length from ambient conditions (assumed constant during egg development)
      call get_tracer_position(space, xyz)
      call interpolate_temp(xyz,temp,istat)
      call interpolate_copper(xyz,cu,istat)  ! unit = kg/m3
      cu = cu*1.0e-6                         ! convert to unit = yg/l
      K1                 = K1_prefac * (temp**K1_expo) ! yg/l
      state%hatch_length = L0_hatch * exp(- cu / K1) + Linf_hatch
      end subroutine set_larvae_hatched

c     ==============================================================
c     output interface - to be filled in later
c     ==============================================================

      subroutine get_prop_state(state,var,bucket,status)
c------------------------------------------------------------  
      type(yolksac_larvae),intent(in) :: state
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
