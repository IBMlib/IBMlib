cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Growth/survival of larvae based one bioenergetics
c     ---------------------------------------------------
c     $Rev: $
c     $LastChangedDate: $
c     $LastChangedBy:  $ 
c
c     Preliminary skeleton for Letcher type larval bioenergetics
c     with daily growth increments based on optimal fouraging.
c     This module is generic and does not contain species-specific parameters,
c     but deals with generic characteristics as length and weight.
c     This module provides both the particle state interface
c     and the embedded larvae interface (for inclusion in a stage resolving upper module)
c
c     Species-specific proterties are obtained through 
c     subroutines in module larval_properties. The prey availability 
c     is assessed through module prey_community. This module composes the primary 
c     class state_attributes of species and state-specific subclasses
c     which are inherited from sub modules. Sub classes are tight, i.e
c     they have public scope. This means they should only be inherited 
c     locally.
c     Module association structure (classes in paranthesis):
c
c                     optimal_forager (-> state_attributes + feeding_larvae)
c       
c                          |                                |
c
c      larval_properties (-> larval_physiology)   prey_community (-> no classes)
c                                       
c                          |              |                 |             
c
c                     particle_state_base (-> local_environment)
c
c     Growth representation: growth is in mass. Length and weight
c     is stored separately. If weight > nominal weight(current length)
c     then current length is updated so weight = nominel weight(length)
c     Length can never decrease. Weight can increase/decrease.
c    
c     Currently, handling time (handling_time) and capture success (capture_sucess) does only
c     depend on predator/prey size ratio. This is the basis for 
c     precalculating foraging windows, depending only on larval size
c     and feeding profitability - otherwise the optimal diet determination
c     must be generalized (slower).
c     The specific functions capture_sucess, handling_time, search_volume  
c     must be parameterized as continuous functions with a consistent derivative
c     for prey sizes lp in the range:
c       lp_min_search_value <= lp <= lpmaxfactor*llarv
c     
c
c     Units: weight  = myg DW
c            length  = mm
c            time    = seconds
c
c     TODO
c       
c       stabilize calculate_ingestion_rate to situation where max is close to p=0
c       egg sub model
c       FF larv sub model
c       settlement sub model
c       more general configuration protocol for ingestion_ref_DB (shape,content)
c       poisson process for encounter over average encounter
c       signal interpolation exceeding lplow_2Dspline/lphigh_2Dspline
c       reset intent for dummies in subroutine calculate_ingestion_rate when testing is done
c
c     Questions to Ute Daewel etal:
c         MDMIN in beta ==  weight0??(currently set so)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module feeding_larval_stage  ! use name change module to switch to particle_state

      use time_tools            ! import clock type
      use particle_tracking     ! space types/methods 
      use physical_fields
      use run_context, only: ctrlfile => simulation_file
      use input_parser
      use output                ! access polytype for get_prop_state/get_metadata_state

      use particle_state_base  
      use larval_properties     ! import: type(larval_physiology) + capture_sucess+ handling_time + search_volume
      use prey_community        ! import: prey_mass + prey_spectrum_density
      use numerical_1d_integrals ! import trapez_integral

      implicit none            
      private     

c     ----------------------------------------------------------------    
c     Class feeding_larvae is for the embedded larvae interface 
c     and contains attributes relating to the larval stage alone
c     along with attributes relating to optimal foraging
c
c     Currently time bundling of growth processes are disabled
c     ----------------------------------------------------------------

      
      type feeding_larvae
      private
         type(larval_physiology) :: larvstate         ! defined in module larval_properties 
      end type
      public :: feeding_larvae ! make type state_attributes visible outside

c     -----------------------------------------------    
c     Class state_attributes is for the 
c     particle state interface and contains all attributes 
c     related to the state of a particle beyond 
c     spatial aspects.
c     If used in a stage manager class that provides state_attributes
c     itself, state_attributes in this module should be excluded at
c     use association like:
c        use feeding_larval_stage,  ONLY: ...
c     or diverted like:
c        use feeding_larval_stage,  state_attributes_FL => state_attributes
c     to avoid name space conflict
c     -----------------------------------------------

      type state_attributes
      private

         type(feeding_larvae) :: larvae
c
c        ---- administrative attributes relating to existence context ----   
c 
         real        :: survival   ! 0 < survival < 1 of this particle 
         logical     :: alive                
         integer     :: death_day  ! Julian day     
         integer     :: sourceBox ! spatial release box
         integer     :: tracerID ! ensemble ID - negative is unset

      end type
      public :: state_attributes ! make type state_attributes visible outside


c.....Particle state interface: (stand-alone)
   
      public :: init_particle_state  ! module operator
      public :: close_particle_state ! module operator
      public :: init_state_attributes
      public :: get_active_velocity
      public :: update_particle_state
      public :: delete_state_attributes 
      public :: write_state_attributes
      public :: get_particle_version
      public :: get_property         ! currently void
      public :: get_metadata_state   ! currently void

c.....Particle sub state interface: (for multi stage manager)
      public :: init_feeding_larval_stage
      public :: close_feeding_larval_stage
      public :: init_state_attributes_FL
      public :: get_active_velocity_FL
      public :: update_particle_state_FL
      public :: delete_state_attributes_FL
c.....All other subroutine than these below are private to the module

      interface get_property
      module procedure get_prop_SA
      module procedure get_prop_FL
      end interface
      
      interface init_state_attributes
      module procedure init_state_attributes_SA
      end interface

      interface get_active_velocity
      module procedure get_active_velocity_SA
      end interface
      
      
      interface update_particle_state
      module procedure update_particle_state_SA
      end interface

      interface set_state_hatched   ! internal interface
      module procedure set_state_hatched_FL
      module procedure set_state_hatched_SA
      end interface

      interface delete_state_attributes
      module procedure delete_state_attributes_SA
      end interface   
  

c     ==============================================================
c     ========             module data section              ========
c     ==============================================================
      integer, parameter    :: verbose = 0 ! control output volume

      integer               :: particle_counter ! module counter for generating IDs for particles

      integer               :: nforpt ! sampling point for ingration integrals
      

      real,parameter        :: lp_min_search_value = 0.0 ! [mm] lower prey size, for optimality searches
      REAL,parameter        :: lpmaxfactor         = 1.0 ! max prey/pred size ratio, for optimality searches 
c
      logical               :: precalc_ingestion_DB = .false. ! later move to input file
      real,allocatable      :: ingestion_ref_DB(:,:,:)   ! (llarv,logZ,julian_day)
      real,allocatable      :: llarv_grid(:)
      real,allocatable      :: logZ_grid(:)              ! log base = e
      real,allocatable      :: jday_grid(:)              ! currently no periodic BC applies
c     ===============================================================
                               contains
c     ===============================================================
     
 
      subroutine init_particle_state()
      call init_feeding_larval_stage()
      end subroutine init_particle_state

      subroutine init_feeding_larval_stage() ! module operator
c     ---------------------------------------------------------------
c     This subroutine should initialize this module 
c     ---------------------------------------------------------------
c     
      type(local_environment) :: local_env ! for testing
      type(larval_physiology) :: ref_larv  ! for testing
      real                    :: lp,llarv,dlp ,drp,rprof
      real                    :: irate,irateDB,irate_max
      real                    :: lp_rmax, rmax, svol, dum
      real                    :: dllarv, dlogZ, djday,z
      integer                 :: m,mpt,i0,j0,k0,i1,j1,k1,i
c     ---------------------------------------------------------------
      write (*,*) trim(get_particle_version())
c     
c     --- pass init to sub modules
c     
      call init_larval_properties() ! module init
      particle_counter = 1      ! ID for next larvae (inactive when embedded)


c     
      call read_control_data(ctrlfile,"ingestion_integral_sampling",
     +     nforpt)
      write(*,*) "init_particle_state: ingestion integral sampling =",
     +            nforpt
      
      call read_control_data(ctrlfile,"precalc_ingestion_DB",
     +                       precalc_ingestion_DB)
      if (precalc_ingestion_DB) then
         call setup_ingestion_db()
      else
         write(*,*) "init_particle_state: optimal ingestion rate "//
     +              "calculated dynamically"
      endif
c     
cc     ==================  test section ==================
c
c      call probe_local_environment((/0.0, 0.0, 1.0/),local_env)
c      lp    = 0.5
c      llarv = 10.0
c      dlp   = 0.01
c      call test_lp_derivatives(lp, llarv, local_env, dlp)
c      rprof = 0.56
c      drp   = 0.01
c      call test_rprof_derivatives(llarv, rprof, local_env, drp)
c      stop 555
c
c      ------------ foraging window ------------
c
c      local_env%zbiomass = 100./1e6/0.32 ! kgDW/m3 converted from mgC/m3
c      local_env%temp     = 18.0
c      local_env%julday   = 100  
c      call test_foraging_window_llarv(ref_larv,local_env,7.0,20.0,0.1)
c
c      ------------ convergence test on nforpt ------------
c
c      local_env%zbiomass = 200./1e6/0.32 ! kgDW/m3 converted from mgC/m3
c      local_env%temp     = 10.0
c      local_env%julday   = 100  
c      ref_larv%length    = 15.0
c      do nforpt=2,30
c         call calculate_ingestion_rate(ref_larv,local_env,irate,rprof) 
c         write(*,*) nforpt,irate,rprof
c      enddo
c      ------------ test DB interpolation ------------
c      local_env%zbiomass = exp(-10.0) ! kgDW/m3 converted from mgC/m3
c      local_env%temp     = 10.0
c      local_env%julday   = 101
c      ref_larv%length    = 15.0
c      call calculate_ingestion_rate(ref_larv,local_env,irate,rprof) 
c      call interpolate_ingestion_rate(ref_larv,local_env,irateDB) 
c      write(*,*) irate, irateDB
c
c      ------------ food utilization degree ------------
c      local_env%zbiomass = 100./1e6/0.32  ! kgDW/m3 converted from mgC/m3
c      local_env%temp     = 10.0
c      local_env%julday   = 101
c      ref_larv%length    = 15.0
c      call find_rank_maximum(ref_larv%length, local_env, lp_rmax, rmax)
c      call search_volume(lp_rmax, ref_larv%length, local_env, svol, dum) !  (svol unit = mm3/sec) 
c      irate_max = local_env%zbiomass * svol    ! kgDW/m3 * mm3/sec = mygDW/sec
c      call calculate_ingestion_rate(ref_larv,local_env,irate,rprof) 
c      write(*,*) irate/irate_max
      

      end subroutine init_feeding_larval_stage



      character*100 function get_particle_version()  
c     ---------------------------------------------------------------
c     ---------------------------------------------------------------
      get_particle_version = "sprat letcher model : $Rev: 147 $"
      end function



      subroutine close_particle_state()
      call close_feeding_larval_stage()
      end subroutine close_particle_state


      subroutine close_feeding_larval_stage() ! module operator
c     ---------------------------------------------------------------
c     Module clean up stuff here
c     (e.g. memory deallocation, MPI finalisation, final dumps etc.)
c     Only deallocate own data, i.e. data allocated in this module
c     ---------------------------------------------------------------
      write(*,*) "close_particle_state():"

      particle_counter               = 0

      if (allocated(ingestion_ref_DB)) deallocate(ingestion_ref_DB)
      if (allocated(llarv_grid))       deallocate(llarv_grid)
      if (allocated(logZ_grid))        deallocate(logZ_grid)  
c
c     --- pass init to sub modules reversely
c
      call close_larval_properties()  ! module close down

      end subroutine close_feeding_larval_stage
      


      subroutine init_state_attributes_SA(state, space, time_dir,             
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
      state%survival   = 1.0
      state%alive      = .true.
      state%sourceBox = emitboxID   ! store releasing box
      state%tracerID = particle_counter    
      particle_counter = particle_counter + 1 ! module data

      call init_state_attributes_FL(state%larvae, space, time_dir,             
     +                                 initdata) 
c     ----------------------------------------------------
      end subroutine init_state_attributes_SA




      subroutine init_state_attributes_FL(state, space, time_dir,             
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
      type(feeding_larvae),intent(out)       :: state
      type(spatial_attributes),intent(inout) :: space
      real,intent(in)                        :: time_dir            
      character*(*),intent(in)               :: initdata
c     ---------------------------------------------------------------
      if (time_dir<0) then
         write(*,*) "init_state_attributes: time_dir<0 not implemented" 
         stop
      endif
      
c     parse initdata string here, if needed

      call set_larvae_hatched(state%larvstate, space) ! currently only support this init state
  
c     ----------------------------------------------------
      end subroutine init_state_attributes_FL

      


      subroutine get_active_velocity_SA(state, space, v_active)
c     --------------------------------------------------------------
c     Currently no migration pattern implemented (set v_active = 0 0 0)
c     If available, this subroutine should determine
c     v_active from (state, space) of the particle
c
c     This should really be in species properties
c     --------------------------------------------------------------
      type(state_attributes), intent(in)   :: state
      type(spatial_attributes), intent(in) :: space      
      real, intent(out)                    :: v_active(:) ! meter/second
c 
      v_active = 0.
c     --------------------------------------------------------------
      end subroutine get_active_velocity_SA


      subroutine get_active_velocity_FL(state, space, v_active)
c     --------------------------------------------------------------
c     Currently no migration pattern implemented (set v_active = 0 0 0)
c     If available, this subroutine should determine
c     v_active from (state, space) of the particle
c
c     This should really be in species properties
c     --------------------------------------------------------------
      type(feeding_larvae), intent(in)     :: state
      type(spatial_attributes), intent(in) :: space      
      real, intent(out)                    :: v_active(:) ! meter/second
c 
      v_active = 0.
c     --------------------------------------------------------------
      end subroutine get_active_velocity_FL



      subroutine update_particle_state_SA(state, space, dt)
c     --------------------------------------------------------------
c     --------------------------------------------------------------
      type(state_attributes), intent(inout)   :: state
      type(spatial_attributes), intent(inout) :: space
      real,intent(in)                         :: dt   ! in seconds
c
      real                                    :: mortality_rate
      logical                                 :: die, next ! stage shift/death request
c     --------------------------------------------------------------
      if (.not.state%alive) return
      call update_particle_state_FL(state%larvae, space, dt,  
     +                       mortality_rate, die, next)
c
c     Update administrative attributes - currently ignore next, if set
c
      state%survival = state%survival * exp(-mortality_rate*dt)
      if (die) call set_state_dead(state, space)  
c     --------------------------------------------------------------
      end subroutine update_particle_state_SA




      subroutine update_particle_state_FL(state, space, dt, 
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
      type(feeding_larvae), intent(inout)     :: state
      type(spatial_attributes), intent(inout) :: space
      real,intent(in)                         :: dt ! in seconds
      real,intent(out)                        :: mortality_rate
      logical,intent(out)                     :: die, next
c
      type(local_environment)                 :: local_env
      real                                    :: xyz(3)
      integer                                 :: status
      real                                    :: irate, rprof
c     --------------------------------------------------------------
      if (dt<0) stop "negative dt values currently blocked"

      call get_tracer_position(space, xyz)
      call probe_local_environment(xyz, local_env) ! cache to avoid repeated interpolations
     
c     --- determine ingestion rate - assume it is a visual predator     
      if (local_env%light) then 
         if (precalc_ingestion_DB) then
            call interpolate_ingestion_rate(state%larvstate,local_env,
     +                                      irate)         ! currently just average, no Poisson  
         else
            call calculate_ingestion_rate(state%larvstate, local_env, 
     +                                    irate, rprof)   
         endif
      else
         irate = 0.0   ! no light -> no feeding
      endif
      
      call grow_larvae(state%larvstate, local_env, irate, dt) ! grow in mass, sync in length
      call inquire_stage_change(state%larvstate, next) 
      call evaluate_mortality_rate(state%larvstate, local_env, 
     +                             mortality_rate, die)

c      write(55,*) state%larvstate%length,  state%larvstate%weight  ! hack

      call clear_local_environment(local_env) ! in case memory were allocated
c     --------------------------------------------------------------
      end subroutine update_particle_state_FL




      subroutine delete_state_attributes_SA(state) 
      type(state_attributes), intent(inout) :: state 
      call delete_state_attributes_FL(state%larvae)
      end subroutine 

      subroutine delete_state_attributes_FL(state) 
      type(feeding_larvae), intent(inout) :: state 
      call close_larval_physiology(state%larvstate)
      end subroutine 

      subroutine write_state_attributes(state)
      type(state_attributes), intent(in) :: state 
      end subroutine 

      

c     ==============================================================
c     ========      internal (non public) subroutines       ========
c     ==============================================================



      subroutine get_condition_number(state, phi)
c     ---------------------------------------------------
      type(feeding_larvae), intent(in)     :: state 
      real, intent(out)                    :: phi
      real                                 :: nomweight
c     ---------------------------------------------------
      call length_to_nominal_weight(state%larvstate%length, nomweight)
      phi = state%larvstate%weight / nomweight
      end subroutine get_condition_number

      


      subroutine set_state_hatched_FL(state, space) 
c     ---------------------------------------------------
c     Initialize larvae in newly hatched state
c     make no assumptions about previous history
c     ---------------------------------------------------
      type(feeding_larvae), intent(inout)    :: state
      type(spatial_attributes),intent(inout) :: space
c     ---------------------------------------------------
      call set_tracer_mobility_free(space)   
      call set_larvae_hatched(state%larvstate,space)
 
      end subroutine set_state_hatched_FL



      subroutine set_state_hatched_SA(state, space) 
c     ---------------------------------------------------
c     Initialize larvae in newly hatched state
c     make no assumptions about previous history
c     ---------------------------------------------------
      type(state_attributes), intent(inout)  :: state
      type(spatial_attributes),intent(inout) :: space
c     ---------------------------------------------------       
      call set_state_hatched(state%larvae, space)  
      state%survival          = 1.0   
      state%alive             = .true.
      end subroutine set_state_hatched_SA




      subroutine set_state_dead(state, space)
c     ---------------------------------------------------
c     Set larval state to dead (there is no overloaded variant for 
c     type feeding_larvae, because it does not know if it is dead or not)
c     leave most variable at value just before larvae died
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





      subroutine eval_gain(llarv,rprof,local_env, g, dgdp)
c     ------------------------------------------------------------------
c     Evaluate the energetic gain g [myg/sec] that a larvae of length llarv
c     have when feeding at profitability level rprof
c     Also evaluate profitability derivatives dgdp
c
c     Notice that the analytic rule for derivatives of an integral 
c     can not be used when a finite numerical sampling scheme is applied
c     (there will be a small discrepancy between analytic rule and a 
c      numerical differentiation, which will prevent convergence of 
c      optimization)
c
c     nforpt (sampling points for foraging integrals) is module data setup 
c     at module initialization
c     ------------------------------------------------------------------
      real,intent(in)                    :: llarv
      real,intent(in)                    :: rprof
      type(local_environment),intent(in) :: local_env ! ambient conditions
      real,intent(out)                   :: g
      real,intent(out)                   :: dgdp
c
      integer :: i
      real    :: lp_intv,lp
      real    :: y(nforpt)   ! automatic array to represent ingestion integrands
      real    :: dy(nforpt)  ! automatic array to represent integrand derivatives
      real    :: w(nforpt)   ! automatic array to represent handling time integrands
      real    :: dw(nforpt)  ! automatic array to represent integrand derivatives

      real    :: lp0, lp1, h
      real    :: ing, ding_dlp0, ding_dlp1
      real    :: xht, dxht_dlp0, dxht_dlp1
      real    :: ding_dp, dxht_dp,  b, db_dp 
      real    :: dlp0_dp, dlp1_dp  
      real    :: enc, pval, htime, denc_dlp, dpval_dlp, dhtime_dlp
c     ------------------------------------------------------------------ 
c
c     --- retreieve foraging window [lp0,lp1] for this (llarv,rprof)
c         and the derivatives with respect to profitability
       
      call find_foraging_window(llarv,local_env,rprof, 
     +                          lp0, dlp0_dp, lp1, dlp1_dp)
c
c     --- setup ingestion / handling time integrands and evaluate integrals
c         along with derivatives
c      
      lp_intv = lp1-lp0
      h       = lp_intv/(nforpt-1)
      do i=1,nforpt
         lp = lp0 + lp_intv*(i-1)/(nforpt-1)
         call encounter_rate(lp, llarv, local_env, enc, denc_dlp)
         call prey_value(lp, llarv, local_env, pval, dpval_dlp)
         call handling_time(lp, llarv, local_env, htime, dhtime_dlp)
c      
c        y  =  encounter_rate * prey_value        
c        w  =  encounter_rate * handling_time
c
         y(i)   = enc * pval
         dy(i)  = denc_dlp * pval  +  enc * dpval_dlp  
c
         w(i)   = enc * htime 
         dw(i)  = denc_dlp * htime  +  enc * dhtime_dlp       
c
      enddo
c
      call trapez_integral(y, dy, h, ing, ding_dlp0, ding_dlp1)
      call trapez_integral(w, dw, h, xht, dxht_dlp0, dxht_dlp1) 
      xht = xht + 1.0  ! per unit time (does not affect derivatives)
c
c     --- compute end point derivatives consistently
c         with the numerical integration scheme using chain rule
c        
c
      ding_dp   = ding_dlp0*dlp0_dp + ding_dlp1*dlp1_dp
      dxht_dp   = dxht_dlp0*dlp0_dp + dxht_dlp1*dlp1_dp
      b         = ding_dp*xht - ing*dxht_dp
c
c     --- finally collect gain g and derivates (dgdp,d2gdp2) from
c         computed integrals using chain rule
c     
      g      = ing/xht
      dgdp   = b/xht**2

      end subroutine eval_gain

 

      subroutine calculate_ingestion_rate(ref_larv, local_env, irate, 
     +                                    rprof) 
c     ------------------------------------------------------------------
c     Calculate the intake rate irate [myg/sec] based on optimal foraging theory
c     for a reference larvae in state ref_larv when subject to 
c     local environment local_env.
c     0 <  rprof < 1 is the relative profitability corresponding to the optimal feeding rate irate
c     This algorithm localizes the relative feeding profitability level rprof corresponding
c     to optimal gain using bisection of gain function value/slope
c     ------------------------------------------------------------------
      type(larval_physiology),intent(in) :: ref_larv
      type(local_environment),intent(in) :: local_env 
      REAL,intent(out)                   :: irate ! optimal feeding rate
      REAL,intent(out)                   :: rprof ! relative profitability for optimal feeding rate
c
      REAL                               :: llarv
      real                               :: p, g,dgdp
      real                               :: pleft,pright,gleft,gright
      real                               :: wdist
      integer                            :: ist, nsteps,i,ib
      
      real, parameter                    :: wdist0 = 0.1   ! initial left/right search bracket 
      
      real, parameter                    :: multifac      = 0.303
      integer, parameter                 :: max_loc_steps = 20
      real, parameter                    :: p_resol      = 1.e-7  ! resolution limit on gain maximum
      logical                            :: printinfo
c     ------------------------------------------------------------------
      llarv     = ref_larv%length 
c      
c     --- establish profitability search bracket [pleft,pright] 
c         for gain maximum by expanding toward 
c         absolute limits p=0 and p=1
c      
      wdist = wdist0
      do ib = 1, max_loc_steps
         pleft = wdist
         call eval_gain(llarv,pleft,local_env, gleft, dgdp)  
         if (dgdp > 0) exit ! pleft is a valid left delimiter
         wdist = wdist*multifac
      enddo
      if (dgdp <= 0) then
         write(*,211) "left"
         stop
      endif

      wdist = wdist0
      do ib = 1, max_loc_steps
         pright = 1.0 - wdist
         call eval_gain(llarv,pright,local_env, gright, dgdp)    
         if (dgdp < 0) exit ! pright is a valid right delimiter
         wdist = wdist*multifac
      enddo
      if (dgdp >= 0) then
         write(*,211) "right"
         stop
      endif
 211  format("calculate_ingestion_rate: unable to establish ",a,
     +       " bracket for gain maximum ")
   
c
c     Bisect maximum gain interval 
c
      nsteps = int(log((pright-pleft)/p_resol) / log(2.))
      
      do ist = 1, nsteps
         p = 0.5*(pright + pleft)
         call eval_gain(llarv, p, local_env, g, dgdp)   
         if (dgdp > 0) then
            pleft  = p
            gleft  = g
         else
            pright = p
            gright = g
         endif
c         write(*,*) "[",pleft,">",g,"<",pright, "]"
      enddo     
c
c     Finally set irate = max gain 
c
      rprof   = 0.5*(pleft + pright) ! final estimator 
      call eval_gain(llarv, rprof, local_env, irate, dgdp) 
      

      end subroutine calculate_ingestion_rate
     



      subroutine interpolate_ingestion_rate(larvstate, local_env, irate) 
c     ------------------------------------------------------------------
c     Interpolate (upper limit) optimal ingestion rate from ingestion_ref_DB
c     which was initialized in setup_ingestion_db
c     This interpolator assumes regular llarv,logZ,jday grids
c     Extract grid spacing/size dynamically rather than storing it as module data
c
c     subroutine is validated
c     ------------------------------------------------------------------
      type(larval_physiology),intent(in) :: larvstate
      type(local_environment),intent(in) :: local_env 
      REAL,intent(out)                   :: irate ! optimal feeding rate
c
      REAL                               :: llarv,logZ,jday,vc(8) ! 8=2*2*2
      REAL                               :: sllarv,slogZ,sjday
      REAL                               :: dllarv,dlogZ,djday
      integer                            :: i0,i1,j0,j1,k0,k1
      integer                            :: nllarv,nlogZ,njday,deriv
      real, parameter                    :: zbiomass_min = 1.0e-20 ! avoid log() error
c     ------------------------------------------------------------------
c     
c     locate grid coordinates in ingestion_ref_DB(nllarv,nlogZ,njday)
c
      if (.not.precalc_ingestion_DB) then
         write(*,*) "interpolate_ingestion_rate: precalc_ingestion_DB=F"
         stop 
      endif

      llarv = larvstate%length
      nllarv = size(llarv_grid)
      if ((llarv < llarv_grid(1)).or.(llarv > llarv_grid(nllarv))) then
         write(*,*) "interpolate_ingestion_rate: " //
     +              "extrapolating for dimension llarv, llarv=",llarv
      endif
      dllarv = llarv_grid(2)-llarv_grid(1) ! regular grid
      i0     = 1 + int((llarv-llarv_grid(1)) / dllarv)
      i1     = i0 + 1
      i0     = max(1,min(i0,nllarv))
      i1     = max(1,min(i1,nllarv))
      sllarv = (llarv - llarv_grid(i0))/dllarv
c      
      logZ   = log(max(local_env%zbiomass, zbiomass_min)) ! avoid log() error
      nlogZ  = size(logZ_grid)
      if ((logZ < logZ_grid(1)).or.(logZ > logZ_grid(nlogZ))) then
         write(*,*) "interpolate_ingestion_rate: " //
     +              "extrapolating for dimension logZ, logZ=",logZ
      endif
      dlogZ  = logZ_grid(2)-logZ_grid(1)   ! regular grid
      j0     = 1 + int((logZ-logZ_grid(1)) / dlogZ)
      j1     = j0 + 1
      j0     = max(1,min(j0,nlogZ))
      j1     = max(1,min(j1,nlogZ))
      slogZ = (logZ - logZ_grid(j0))/dlogZ
c      
      jday  = float(local_env%julday)  
      njday = size(jday_grid)              ! integer
      if ((jday < jday_grid(1)).or.(jday > jday_grid(njday))) then
         write(*,*) "interpolate_ingestion_rate: " //
     +              "extrapolating for dimension jday, jday=",jday
      endif
      djday = jday_grid(2)-jday_grid(1) ! regular grid
      k0     = 1 + int((jday-jday_grid(1)) / djday)
      k1     = k0 + 1
      k0     = max(1,min(k0,njday))
      k1     = max(1,min(k1,njday))
      sjday = (jday - jday_grid(k0))/djday
c      
c     vc(1:8) = v000,v001,v010,v011, v100,v101,v110,111
c
      vc(1) = ingestion_ref_DB(i0,j0,k0)
      vc(2) = ingestion_ref_DB(i0,j0,k1)
      vc(3) = ingestion_ref_DB(i0,j1,k0)
      vc(4) = ingestion_ref_DB(i0,j1,k1)
      vc(5) = ingestion_ref_DB(i1,j0,k0)
      vc(6) = ingestion_ref_DB(i1,j0,k1)
      vc(7) = ingestion_ref_DB(i1,j1,k0)
      vc(8) = ingestion_ref_DB(i1,j1,k1)
      deriv = 0        ! request value interpoaltion
c
      call interp_3Dbox_data(sllarv,slogZ,sjday,vc,deriv,irate)
c     
      end subroutine interpolate_ingestion_rate


      subroutine encounter_rate(lp, llarv, local_env, enc, denc_dlp)
c     --------------------------------------
c     Evaluate the encounter rate density enc [unit == individuals/mm/sec] and derivatives
c     of prey at length lp [mm] of a larval predator of length llarv [mm]
c     The encounter rate of prey in a length interval [lp - 0.5*dL; lp + 0.5*dL]
c     corresponds to dN = enc_rate*dL  [unit == individuals/sec]
c     Encounter rate density is evaluated as clearence volume [unit == mm3/sec]
c     times density of prey [unit == individuals/mm3/mm] 
c
c     Applies prey_spectrum_density (from module prey_community)
c     and search_volume (from module larval_properties)          
c     --------------------------------------
      REAL,intent(in)                    :: lp        ! prey length [mm]
      REAL,intent(in)                    :: llarv     ! larval length [mm]
      type(local_environment),intent(in) :: local_env ! ambient conditions
      REAL,intent(out)                   :: enc       ! encounter rate density
      REAL,intent(out),optional          :: denc_dlp  ! (d/dlp) encounter rate density   

      REAL                               :: dens, ddens_dlp
      REAL                               :: svol, dsvol_dlp
c     --------------------------------------
      if ( present(denc_dlp) ) then
         call prey_spectrum_density(lp, local_env, dens, ddens_dlp)
         call search_volume(lp, llarv, local_env, svol, dsvol_dlp)
         enc      = dens * svol
         denc_dlp = ddens_dlp * svol  +  dens * dsvol_dlp
      else
         call prey_spectrum_density(lp, local_env, dens)
         call search_volume(lp, llarv, local_env, svol)
         enc = dens * svol
      endif
     
      end subroutine encounter_rate

 

      subroutine prey_value(lp,llarv,local_env,pval,dpval_dlp)
c     --------------------------------------
c     Evaluate prey value = prey mass * capture success  [myg] 
c     of prey length lp      
c     --------------------------------------
      REAL,intent(in)                    :: lp        ! prey length [mm]
      REAL,intent(in)                    :: llarv     ! larval length [mm]
      type(local_environment),intent(in) :: local_env 
      REAL,intent(out)                   :: pval      ! prey value
      REAL,intent(out),optional          :: dpval_dlp ! (d/dlp) prey value

      REAL                               :: pmass, dpmass_dlp
      REAL                               :: csucc, dcsucc_dlp
c     --------------------------------------
      if ( present(dpval_dlp) ) then
         call prey_mass(lp, pmass, dpmass_dlp)
         call capture_sucess(lp, llarv, local_env, csucc, dcsucc_dlp)
         pval      = pmass * csucc
         dpval_dlp = dpmass_dlp * csucc  +  pmass * dcsucc_dlp
      else   
         call prey_mass(lp, pmass)
         call capture_sucess(lp, llarv, local_env, csucc)
         pval = pmass * csucc
      endif

      end subroutine prey_value


      
      subroutine rank(lp, llarv, local_env, rnk, drnk_dlp)
c     ----------------------------------------------------------------------------
c     Evaluate            
c                           prey_value    [myg] 
c         prey rank =   -----------------------------
c                        handling time    [seconds]
c   
c     of prey length lp, where prey_value = prey_mass * capture success
c     and optionally (d/dlp) prey rank 
c     ----------------------------------------------------------------------------
      REAL,intent(in)                    :: lp        ! prey length [mm]
      REAL,intent(in)                    :: llarv     ! larval length [mm]
      type(local_environment),intent(in) :: local_env ! ambient conditions
      REAL,intent(out)                   :: rnk       ! prey rank 
      REAL,intent(out),optional          :: drnk_dlp  ! (d/dlp) prey rank 
 
      REAL                               :: pval, dpval_dlp 
      REAL                               :: htime, dhtime_dlp 
c     ----------------------------------------------------------------------------
      if ( present(drnk_dlp) ) then
         call prey_value(lp, llarv, local_env, pval, dpval_dlp)
         call handling_time(lp, llarv, local_env, htime, dhtime_dlp)
         rnk      = pval/htime
         drnk_dlp = (dpval_dlp*htime - pval*dhtime_dlp)/htime**2
      else
         call prey_value(lp, llarv, local_env, pval)
         call handling_time(lp, llarv, local_env, htime)
         rnk = pval/htime
      endif
   
      end subroutine rank


      
      subroutine find_rank_maximum(llarv, local_env, lp, rmax)
c     --------------------------------------
c     Locate prey length lp corresponding to the 
c     prey rank maximum rmax for a larvae 
c     with larval length llarv using rank slope bisection
c     to within accuracy lp_resolution
c
c     Assert lp_min_search_value <= lp <= llarv to establish a search bracket
c     --------------------------------------
      REAL,intent(in)                    :: llarv ! larval length [mm]
      type(local_environment),intent(in) :: local_env ! ambient conditions
      REAL,intent(out)                   :: lp    ! prey length at rank maximum [mm]
      REAL,intent(out)                   :: rmax  ! rank maximum [myg/second]
c
       
      REAL,parameter     :: lp_resolution = 1.e-6 ! [mm]
   
      REAL               :: lpleft,lpright,rnk,drnk_dlp
      integer            :: nsteps, ist
c     -------------------------------------
      lpleft  = lp_min_search_value  ! assert (d/dlp) rank (lleft)  > 0
      lpright = lpmaxfactor*llarv    ! assert (d/dlp) rank (lright) < 0
c     --- check assertions ---
      call rank(lpleft, llarv, local_env, rnk, drnk_dlp)   
      if (drnk_dlp < 0) then   ! 
         stop "find_rank_maximum: invalid left bracket"
      endif
      call rank(lpright, llarv, local_env, rnk, drnk_dlp) 
      if (drnk_dlp > 0) then
         stop "find_rank_maximum: invalid right bracket"
      endif
      nsteps = int(log((lpright-lpleft)/lp_resolution) / log(2.))
      do ist = 1, nsteps
         lp = 0.5*(lpleft + lpright)
         call rank(lp, llarv, local_env, rnk, drnk_dlp) 
         if (drnk_dlp > 0) then
            lpleft  = lp
         else
            lpright = lp
         endif
      enddo

      lp   = 0.5*(lpleft + lpright) ! final estimator 
      call rank(lp, llarv, local_env, rmax)
      
      end subroutine find_rank_maximum



      subroutine find_foraging_window(llarv,local_env,rprof,
     +                  lplow, dlplow_drp,
     +                  lphigh, dlphigh_drp)
c     -------------------------------------------------------
c     Locate the larval foraging prey size interval 
c     window [lplow,lphigh] where rprof < rank(lp)/max_rank < 1
c     along with derivatives (dlplow_drp, dlphigh_drp) of (lplow,lphigh) with respect to rprof
c
c     These derivatives are determined exactly by perturbation theory from 
c     the derivatives of rank(lp). In the singular limit rprof->1
c     derivatives are assigned +/- infty, with sign from left side limit of rprof->1
c
c     Assume optimal foraging window is a contiguous interval 
c     Assume optimal foraging window is within the
c     search interval [lp_min_search_value, lpmaxfactor*llarv]
c 
c     Intercept singuar limit rprof -> 1
c     Currently do not intercept limit rprof -> 0
c     -------------------------------------------------------- 
      REAL,intent(in)                    :: llarv                   ! larval length [mm]
      type(local_environment),intent(in) :: local_env ! ambient conditions
      REAL,intent(in)                    :: rprof                   ! minimal relative profitability [1/second]
      REAL,intent(out)                   :: lplow, lphigh           ! foraging window [mm]
      REAL,intent(out)                   :: dlplow_drp, dlphigh_drp ! derivatives of foraging window [mm]
c
      REAL               :: lp_rmax,rmax,minrank,rnk,drnk_dlp
      REAL               :: lplow_limit, lphigh_limit 

      REAL,parameter     :: rprof_singular = 1.e-4 ! interception limit when rprof too close to 1
      REAL,parameter     :: infty = 1.0e30         ! deriv value returned in singular limit
c     -------------------------------------------------------- 
c
c     locate rank maximum as delimiter for locating lplow,lphigh
c
      call find_rank_maximum(llarv, local_env, lp_rmax, rmax) 

      if ((rprof>1.0).or.(rprof<0.0)) then
         write(*,*) "find_foraging_window: rprof=", rprof
         stop "find_foraging_window: rprof range error"
      endif
c
c     Capture singular limit rprof=1, where lplow=lphigh=lp_rmax
c
      if (abs(rprof-1.0) < rprof_singular) then
         lplow         = lp_rmax
         lphigh        = lp_rmax
c        --- set derivatives to infinity with sign corresponding to left side limit of rprof->1
         dlplow_drp    = +infty
         dlphigh_drp   = -infty         
         return
      else
c
c     else locate lplow/lphigh by newton_search
c  
         minrank      = rmax*rprof 
         lplow_limit  = lp_min_search_value
         lphigh_limit = lpmaxfactor*llarv
         call newton_search(lplow_limit, lp_rmax,      minrank, lplow)
         call newton_search(lp_rmax,     lphigh_limit, minrank, lphigh)
c        --- set derivatives (rescale to relative profitability = profitability/rmax)
         call rank(lplow, llarv, local_env, rnk, drnk_dlp)
         dlplow_drp    = rmax/drnk_dlp ! exact limit from perturbation theory 
         call rank(lphigh, llarv, local_env, rnk, drnk_dlp)
         dlphigh_drp   = rmax/drnk_dlp ! exact limit from perturbation theory
      endif

      contains

      subroutine newton_search(xleft, xright, r0, x)
c     ------------------------------------------------
c     Solve rank(x,llarv,value)=r0 in [xleft,xright] in local scope
c     using a bracketed Newton search algorithm.
c     The algorithm uses a bisection fall back on a safe bracket,
c     in case the Newton search does not converge satisfactory
c     
c     use llarv and local_env from embedding scope
c     ------------------------------------------------
      real,intent(in)    :: xleft, xright ! input solution bracket
      real,intent(in)    :: r0            ! rank value to solve for
      real,intent(out)   :: x             ! solution 
      real               :: s0,s1,s,f,df,dx,rnk0,rnk1,drnk_dlp
      real               :: x0,x1,resi,relresid      ! dynamical bracket
      logical            :: do_bisec, converged
      character*2        :: steptyp
      REAL,parameter     :: dx_limit       = 1.0e-6 ! resolution limit on x
      REAL,parameter     :: relresid_limit = 1.0e-4 ! resolution limit on rank==r0
      integer            :: iter,i
      integer,parameter  :: max_iter = 20 ! suspend newton search, if exceeded
      logical,parameter  :: log = .false. ! switch on to generate convergence trace
c     ------------------------------------------------
      iter = 1
      x0 = xleft   ! copy initial bracket
      x1 = xright  ! copy initial bracket
      if (log) write(*,317)

c     ---- check assertions ----
      if (x1<x0) stop "newton_search: invalid bracket"
      call rank(x0, llarv, local_env, rnk0)
      s0 = sign(1.0, rnk0-r0)
      call rank(x1, llarv, local_env, rnk1)
      s1 = sign(1.0, rnk1-r0)
      if ((s0*s1)>0) stop "newton_search: root not bracketed"
c     -------------------------------------------------------- 
c     locate the solution x: xleft <= x0 < x < x1 < xright
c       x0,x1 : dynamic bracket on root
c       x     : current value on root
c       dx    : current estimate on proximity to root
c     --------------------------------------------------------
      
c
c        ---- 0) loop initialization ----
c
      x       = 0.5*(x0 + x1)   
      dx      = x1-x0           ! set value in case already abs(dx) < dxlimit
c
      do                ! hybrid Newton/bisection loop
c
c        ---- 1) evaluate current point x ----
c
         call rank(x, llarv, local_env, rnk, drnk_dlp)
c
c        ---- 2) check convergence / print log  ----
c         
         resi      = rnk-r0
         relresid  = abs(resi/max(r0, 1.e-4))
         if (log) write(*,318) iter,x0,x,x1,abs(dx),relresid,steptyp

         converged = ((abs(dx)<dx_limit).and.(relresid<relresid_limit))
         if (converged) then
c            write(32,*) iter
            exit
         endif
c
c        ---- 3) update bracket from last evaluation  ----
c
         s     = sign(1.0, rnk-r0)
         if (s0*s > 0) then     ! update left bracket + sign
            x0 = x
            s0 = s
         else                   ! update right bracket + sign
            x1 = x
            s1 = s
         endif
c
c        ---- 4) generate next point x and step monitor dx ----
c
         if ((drnk_dlp /= 0.0).and.(iter <= max_iter)) then
            dx       = -(rnk-r0)/drnk_dlp   ! Newton step  
            x        = x + dx
            steptyp  = "nw"
            do_bisec = ((x<x0).or.(x>x1))   ! discard Newton step beyond bracket
         else
            do_bisec = .true.   ! Newton step undefined, apply bisection fallback
         endif
c     
         if (do_bisec .or. (iter>max_iter)) then ! bisection fallback 
            x        = 0.5*(x0 + x1)      
            dx       = x1-x0
            steptyp  = "bi"
         endif
         iter            = iter+1

      enddo ! while     

 317  format(20("-"),"newton_search",20("-"))
 318  format(i3," [",3(f12.7),"]  adx=",f12.7," rres=",f12.7,2x,a2)

      end subroutine newton_search ! local function to find_foraging_window
      end subroutine find_foraging_window





      subroutine setup_ingestion_db()
c     ----------------------------------------------------------------
c     Setup upper ingestion limit data base for reference ranges
c     Implicit range for jday_grid is 1..365 (ignore leap years in this context)
c     ----------------------------------------------------------------
      real             :: llarv_len_max, llarv_len_min  ! llarv_grid
      real             :: logZ_min, logZ_max            ! logZ_grid
      real,parameter   :: jday_min = 1.0          ! currently hard code this
      real,parameter   :: jday_max = 365.0        ! currently hard code this
      integer          :: nllarv, nlogZ, njday, i,j,k       
      real             :: dllarv, dlogZ, djday, irate, rprof
      type(local_environment) :: local_env 
      type(larval_physiology) :: ref_larv
c     ----------------------------------------------------------------
      write(*,*)   "init_particle_state: setup_ingestion_db"
c
c     Read reference grid specifications
c
      call read_control_data(ctrlfile,"larval_length_grid_min", 
     +     llarv_len_min) 
      call read_control_data(ctrlfile,"larval_length_grid_max", 
     +     llarv_len_max)     
      call read_control_data(ctrlfile,"larval_length_grid_points", 
     +     nllarv)  
      write(*,732) "larval length",  llarv_len_min,llarv_len_max, nllarv
c     
      call read_control_data(ctrlfile,"logZ_grid_min", logZ_min) 
      call read_control_data(ctrlfile,"logZ_grid_max", logZ_max)
      call read_control_data(ctrlfile,"logZ_grid_points", nlogZ)
      write(*,732) "logZ",  logZ_min, logZ_max, nlogZ
c
      call read_control_data(ctrlfile,"julian_day_grid_points", njday)
      write(*,732) "Julian_day", jday_min, jday_max, njday
      
 732  format("setup_ingestion_db: ",a," grid = [",2f12.7,"] with ",
     +                           i6," points")

c
c     Setup reference grids
c
      allocate( ingestion_ref_DB(nllarv,nlogZ,njday)  )
      allocate( llarv_grid(nllarv) )
      allocate( logZ_grid(nlogZ)   )
      allocate( jday_grid(njday)   )     

      dllarv = llarv_len_max - llarv_len_min
      do i=1,nllarv
         llarv_grid(i) = llarv_len_min + dllarv*(i-1)/(nllarv-1)
      enddo

      dlogZ = logZ_max - logZ_min
      do j=1,nlogZ
         logZ_grid(j) = logZ_min + dlogZ*(j-1)/(nlogZ-1)
      enddo

      djday = jday_max - jday_min
      do k=1,njday
         jday_grid(k) = jday_min + djday*(k-1)/(njday-1)
      enddo
c
c     ------ main loop over reference points ------
c
      local_env%light = .true.
      local_env%temp  = -1.0 ! should not be used
      do i = 1, nllarv
         ref_larv%length = llarv_grid(i)
         do j = 1, nlogZ
            local_env%zbiomass = exp(logZ_grid(j))  ! log base = e
            do k = 1, njday
               local_env%julday = jday_grid(k)
               call calculate_ingestion_rate(ref_larv, local_env, 
     +                                       irate, rprof) 
               if (verbose>0) write(*,532) llarv_grid(i), logZ_grid(j), 
     +                      int(jday_grid(k)), irate, rprof
               ingestion_ref_DB(i,j,k)  = irate
            enddo
         enddo
      enddo
      write(*,*)   "setup_ingestion_db: ingestion_ref_DB initialized"
 532  format("setup_ingestion_db: llarv=",e11.5," logZ=",e11.5,
     +       " jday=",i4," -> ingestion rate =",e11.5," @ p=",f12.8)   
c     ----------------------------------------------------------------
      end subroutine setup_ingestion_db
      


c     subroutine calculate_predation_survival(state,space,pred_surv)
c     --------------------------------------------------- 
c     Letcher Eqs 19-22
c     At some point we may also overlay predator maps (space available)
c     ---------------------------------------------------
c     end subroutine calculate_predation_survival
      

c     ==============================================================
c     testing/validation subroutines
c     ==============================================================


      subroutine test_foraging_window_llarv(ref_larv, local_env, 
     +                                  llarv_min, llarv_max, dllarv)
c     ---------------------------------------------------------------
c     For a larvae ref_larv subject to local environment local_env  
c     scan over larval length to assess optimal feeding conditions
c     larval length scan is on a regular grid [llarv_min,llarv_max] 
c     with sampling spacing dllarv
c     ---------------------------------------------------------------
      type(larval_physiology),intent(inout) :: ref_larv
      type(local_environment),intent(in)    :: local_env 
      real,intent(in)                       :: llarv_min,llarv_max  ! [mm]
      real,intent(in)                       :: dllarv               ! [mm]
c
      real                    :: llarv0,irate,rprof
      real                    :: lplow, dlplow_drp,lphigh, dlphigh_drp
c     ---------------------------------------------------------------
      llarv0          = ref_larv%length ! save orig value
      ref_larv%length = llarv_min
      do while (ref_larv%length < llarv_max)
         call calculate_ingestion_rate(ref_larv,local_env,irate,rprof)
         call find_foraging_window(ref_larv%length,local_env,rprof,
     +                  lplow, dlplow_drp,
     +                  lphigh, dlphigh_drp)
         write(*,*) ref_larv%length,1e3*lplow,1e3*lphigh,rprof,irate
         ref_larv%length = ref_larv%length + dllarv
      end do
      ref_larv%length = llarv0 !  restore orig value

      end subroutine test_foraging_window_llarv



      subroutine test_lp_derivatives(lp, llarv, local_env, dlp)
c     ---------------------------------------------------
c     Debugging/validation utility when testing a new setup
c     Test the consistency of analytical derivatives and
c     numerical derivatives of:
c
c         capture_sucess(lp, llarv, local_env, vector)   from larval_propertie
c         handling_time(lp,  llarv, local_env, vector)   from larval_propertie
c         search_volume(lp,  llarv, local_env, vector)   from larval_properties
c
c         prey_mass(lp, vector)                          from prey_community
c         prey_spectrum_density(lp, local_env, vector)   from prey_community
c
c         encounter_rate(lp, llarv, local_env, vector)   defined locally from above subroutines
c         prey_value(lp, llarv, local_env, vector)       defined locally from above subroutines 
c         rank local(lp, llarv, local_env, vector)       defined locally from above subroutines
c
c     where vector has length 1-3
c     Result is printed to stdout.
c     dlp is the step length for numerical differentiation
c     
c     It is recommended to apply auto doubling of real variables
c     to have exact validation of derivatives (set this in compiler_defaults.mk,
c     for the ifort compiler this is -r8)
c     ---------------------------------------------------
      REAL,intent(in)                    :: lp        ! prey length [mm]
      REAL,intent(in)                    :: llarv     ! larval length [mm]
      type(local_environment),intent(in) :: local_env ! ambient conditions
      REAL,intent(in)                    :: dlp       ! step length for numerical differentiation [mm]
      REAL                               :: val,slope,vm,vp
c     ---------------------------------------------------
      write(*,*) "--- test_lp_derivatives:begin --- "
      write(*,*)
      write(*,*) "lp    =", lp
      write(*,*) "llarv =", llarv
      write(*,*) "dlp   =", dlp
      write(*,*)
c     
      call capture_sucess(lp-dlp, llarv, local_env, vm) 
      call capture_sucess(lp+dlp, llarv, local_env, vp) 
      call capture_sucess(lp,     llarv, local_env, val, slope) ! extract analytic derivatives
      write(*,732) "capture_sucess", val
      write(*,733) "capture_sucess", slope,(vp-vm)/2/dlp
      
      write(*,*)
c
      call handling_time(lp-dlp, llarv, local_env, vm) 
      call handling_time(lp+dlp, llarv, local_env, vp)  
      call handling_time(lp,     llarv, local_env, val, slope) ! extract analytic derivatives
      write(*,732) "handling_time", val
      write(*,733) "handling_time", slope, (vp-vm)/2/dlp
      write(*,*)
c
      call search_volume(lp-dlp, llarv, local_env, vm) 
      call search_volume(lp+dlp, llarv, local_env, vp)   
      call search_volume(lp,     llarv, local_env, val, slope) ! extract analytic derivatives
      write(*,732) "search_volume", val
      write(*,733) "search_volume", slope, (vp-vm)/2/dlp
      write(*,*)
c
      call prey_mass(lp-dlp, vm) 
      call prey_mass(lp+dlp, vp) 
      call prey_mass(lp,     val, slope) ! extract analytic derivatives
      write(*,732) "prey_mass", val
      write(*,733) "prey_mass", slope, (vp-vm)/2/dlp
      write(*,*)
c
      call prey_spectrum_density(lp-dlp, local_env, vm) 
      call prey_spectrum_density(lp+dlp, local_env, vp) 
      call prey_spectrum_density(lp,     local_env, val, slope) ! extract analytic derivatives
      write(*,732) "prey_spectrum_density", val
      write(*,733) "prey_spectrum_density", slope, (vp-vm)/2/dlp
      write(*,*)
c
      call encounter_rate(lp-dlp, llarv, local_env, vm)     
      call encounter_rate(lp+dlp, llarv, local_env, vp)  
      call encounter_rate(lp,     llarv, local_env, val, slope) ! extract analytic derivatives     
      write(*,732) "encounter_rate", val
      write(*,733) "encounter_rate", slope, (vp-vm)/2/dlp 
      write(*,*)
c
      call prey_value(lp-dlp, llarv, local_env, vm)  
      call prey_value(lp+dlp, llarv, local_env, vp) 
      call prey_value(lp,     llarv, local_env, val, slope) ! extract analytic derivatives  
      write(*,732) "prey_value", val
      write(*,733) "prey_value", slope, (vp-vm)/2/dlp 
      write(*,*)
c
      call rank(lp-dlp, llarv, local_env, vm) 
      call rank(lp+dlp, llarv, local_env, vp) 
      call rank(lp,     llarv, local_env, val, slope) ! extract analytic derivatives  
      write(*,732) "rank", val
      write(*,733) "rank", slope, (vp-vm)/2/dlp 
      write(*,*)
      write(*,*) "--- test_lp_derivatives:end --- "
c     
 732  format("testing value   ",a,": value = ",e14.7)
 733  format("testing d/dlp   ",a,": analyt= ",e14.7," num= ",e14.7)    
      
    
      end subroutine test_lp_derivatives

      

      subroutine test_rprof_derivatives(llarv,rprof,local_env,drp)
c------------------------------------------------------------
c     Test the consistency of analytical derivatives and
c     numerical derivatives wrt. rprof
c
c     It is recommended to apply auto doubling of real variables
c     to have exact validation of derivatives (set this in compiler_defaults.mk,
c     for the ifort compiler this is -r8)
c------------------------------------------------------------
      REAL,intent(in)                    :: llarv     ! larval length [mm]
      REAL,intent(in)                    :: rprof     ! relative profitability (0 < rprof < 1)
      type(local_environment),intent(in) :: local_env ! ambient conditions
      REAL,intent(in)                    :: drp       ! step in rprof for numerical differentiation
      REAL                        :: g,dgdp
      REAL                        :: gm,gp,dum1,dum2
      REAL                        :: lplowm,lplowp,lplow,dlplow_drp
      REAL                        :: lphighm,lphighp,lphigh,dlphigh_drp
c------------------------------------------------------------
      write(*,*) "test_rprof_derivatives:"

      call eval_gain(llarv,rprof-drp,local_env, gm, dum1)
      call eval_gain(llarv,rprof+drp,local_env, gp, dum1)
      call eval_gain(llarv,rprof,    local_env, g,  dgdp)
      write(*,287) "eval_gain","g", g
      write(*,288) "eval_gain","dg_dp",dgdp,"dg_dp",(gp-gm)/2/drp 

      call find_foraging_window(llarv, local_env, rprof-drp,
     +                          lplowm, dum1, lphighm, dum2)
      call find_foraging_window(llarv, local_env, rprof+drp,
     +                          lplowp, dum1, lphighp, dum2)
      call find_foraging_window(llarv, local_env, rprof,
     +                          lplow, dlplow_drp, lphigh, dlphigh_drp)
      write(*,287) "find_foraging_window","lplow", lplow
      write(*,288) "find_foraging_window","dlplow_drp", dlplow_drp,
     +             "dlplow_drp", (lplowp-lplowm)/2/drp

      write(*,287) "find_foraging_window","lphigh", lphigh
      write(*,288) "find_foraging_window","dlphigh_drp", dlphigh_drp,
     +             "dlphigh_drp", (lphighp-lphighm)/2/drp
       

 287  format(a20,":",a12,"         ="e14.7)
 288  format(a20,":",a12,"(analyt) =",e14.7,a12,"(num) =",e14.7)

      end subroutine test_rprof_derivatives


c     ==============================================================
c     output interface - to be filled in later
c     ==============================================================

      
      subroutine get_prop_SA(state,var,bucket,status)
c--------------------------------------------------------------------------
c     Dig out current/virtual data from object using get_property interface
c     into bucket given the request var
c     Return status = 0, if data was delivered
c     Return status = 1, if data could not be resolved 
c
c     Pass unhandled data requests onto state%larvae
c--------------------------------------------------------------------------
      type(state_attributes),intent(in) :: state
      type(variable),intent(in)         :: var
      type(polytype), intent(out)       :: bucket
      integer, intent(out)              :: status
c
      integer,parameter                 :: namlen = 999
      character(len=namlen)             :: name
      character*(*),parameter           :: survival   = "survival"   ! query name
      character*(*),parameter           :: sourceBox  = "sourceBox"  ! query name
      character*(*),parameter           :: tracerID   = "tracerID" ! query name
c------------------------------------------------------------
      name = adjustl(get_name(var))
      if (name(1:len(survival)) == survival) then
         call construct(bucket, survival, state%survival)
         status = 0 ! data was extracted successfully
      elseif (name(1:len(sourceBox)) == sourceBox) then
         call construct(bucket, sourceBox, state%sourceBox)
         status = 0 ! data was extracted successfully
      elseif (name(1:len(tracerID)) == tracerID) then
         call construct(bucket, tracerID, state%tracerID)
         status = 0 ! data was extracted successfully
      else  ! pass other requests of larval sub component
         call get_prop_FL(state%larvae,var,bucket,status)
      endif
      end subroutine get_prop_SA


      subroutine get_prop_FL(state,var,bucket,status)
c------------------------------------------------------------  
c     Dig out current/virtual data from object using get_property interface
c     into bucket given the request var
c     Return status = 0, if data was delivered
c     Return status = 1, if data could not be resolved 
c
c     Pass unhandled data requests onto state%larvae
c--------------------------------------------------------------------------
      type(feeding_larvae),intent(in)   :: state
      type(variable),intent(in)         :: var
      type(polytype), intent(out)       :: bucket
      integer, intent(out)              :: status
c
      integer,parameter                 :: namlen = 999
      character(len=namlen)             :: name  
      character*(*),parameter           :: std_length = "std_length" ! query name
      character*(*),parameter           :: dry_weight = "dry_weight" ! query name
c------------------------------------------------------------
      name = adjustl(get_name(var))
      if      (name(1:len(std_length)) == std_length) then
         call construct(bucket, std_length, state%larvstate%length)
         status = 0 ! data was extracted successfully
      elseif  (name(1:len(dry_weight)) == dry_weight) then
         call construct(bucket, dry_weight, state%larvstate%weight)
         status = 0 ! data was extracted successfully
      else  ! unhandled, bail out
         status = 1 
      endif
      end subroutine get_prop_FL



      subroutine get_metadata_state(var_name,var,status)
c     ------------------------------------------------------------  
      character(*), intent(in)   :: var_name
      type(variable),intent(out) :: var
      integer, intent(out)       :: status
c     ------------------------------------------------------------  
      status=0 !Defaults to variable found
      select case (var_name)
      case ("tracerID")
         call construct(var,"tracerID","tracer ID number",
     +     units="-",fmt="(i6)",type="int")
      case ("sourceBox")
         call construct(var,"sourceBox","source box ID",
     +     units="-",fmt="(i6)",type="int")
      case ("std_length")
         call construct(var,"std_length","larval standard length",
     +     units="mm",fmt="(f8.3)",type="real")
      case ("dry_weight")
         call construct(var,"dry_weight","larval dry weight",
     +     units="micro gram",fmt="(f12.3)",type="real")
      case ("survival")
         call construct(var,"survival","survival rate",
     +     units="[0..1]",fmt="(f10.7)",type="real")   
      case default
        status=1  !Cannot find variable
      end select
      end subroutine get_metadata_state

      end module
