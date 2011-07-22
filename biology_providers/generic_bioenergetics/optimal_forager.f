cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Growth/survival based on bioenergetics
c     ---------------------------------------------------
c     $Rev: $
c     $LastChangedDate: $
c     $LastChangedBy:  $ 
c
c     Preliminary skeleton for Letcher type egg/larval bioenergetics
c     with daily growth increments based on optimal fouraging.
c     This module is generic and does not contain species-specific parameters,
c     but deals with generic characteristics as length and weight.
c     Species-specific proterties are obtained through 
c     subroutines in module larval_properties. The prey availability 
c     is assessed through module prey_community. This module composes the primary 
c     class state_attributes of species and state-specific subclasses
c     which are inherited from sub modules. Sub classes are tight, i.e
c     they have public scope. This means they should only be inherited 
c     locally.
c     Module association structure (classes in paranthesis):
c
c                     optimal_forager (-> state_attributes)
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
c     Currently, handling time (HT) and capture success (CS) does only
c     depend on predator/prey size ratio. This is the basis for 
c     precalculating foraging windows, depending only on larval size
c     and feeding profitability - otherwise the optimal diet determination
c     must be generalized (slower).
c
c     Units: weight  = myg DW
c            length  = mm
c            time    = seconds
c
c     TODO
c       
c       egg sub model
c       FF larv sub model
c       settlement sub model
c       poisson process for encounter over average encounter
c       
c       
c
c     Questions to Ute Daewel etal:
c         MDMIN in beta ==  weight0??(currently set so)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module particle_state

      use time_tools           ! import clock type
      use particle_tracking    ! space types/methods 
      use physical_fields
      use run_context, only: ctrlfile => simulation_file
      use input_parser
      use output               ! access polytype for get_prop_state/get_metadata_state

      use particle_state_base  
      use larval_properties   ! import type(physiology) 
      use prey_community
      use spline

      implicit none            
      private     

c     -----------------------------------------------    
c     Class state_attributes contains all attributes 
c     related to the state of a particle beyond 
c     spatial aspects
c     -----------------------------------------------
      type state_attributes
      private

        type(larval_physiology) :: larvstate   ! defined in module larval_properties 

        real        :: survival           ! 0 < survival < 1 of this particle 
        logical     :: alive                
        integer     :: death_day          ! Julian day
        real        :: time_since_growth  ! = sum(dt) - for growth aggregation
        real        :: tempxdt            ! = sum(temp*dt) - for averaging ambient temperature

        integer     :: orig_boxID   ! spatial release box
        integer     :: particleID   ! ensemble ID - negative is unset
      
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
      public :: get_property        ! currently void
      public :: get_metadata_state  ! currently void


      interface get_property
        module procedure get_prop_state
      end interface
      

c     ==============================================================
c     ========             module data section              ========
c     ==============================================================
      integer, parameter    :: verbose = 0    ! control output volume

      integer               :: particle_counter   ! module counter for generating IDs for particles

      type(spline_2D)       :: lplow_2Dspline   ! predation window, lower lp/llarv
      type(spline_2D)       :: lphigh_2Dspline  ! predation window, upper lp/llarv

      integer               :: nforpt           ! sampling point for ingeation integrals
      real, parameter       :: growth_increment_interval = 86400 ! later move to input

      logical, parameter    :: slope = .true.   ! for d/dlp functions
      logical, parameter    :: value = .false.  ! for d/dlp functions

      real,parameter        :: lp_min_search_value = 1.e-3 ! [mm] lower prey size, for optimality search
c     ===============================================================
                                  contains
c     ===============================================================
     
 
      subroutine init_particle_state()  ! module operator
c     ---------------------------------------------------------------
c     This subroutine should initialize this module 
c     ---------------------------------------------------------------
      real, allocatable :: rprof_grid(:), llarv_grid(:)
      integer           :: nrprof, nllarv, i
      real              :: rprofmin, dllarv,llarv_len_max,llarv_len_min
c     ---------------------------------------------------------------
      write (*,*) trim(get_particle_version())
      particle_counter = 1 ! ID for next larvae
c
c     --- setup optimal foraging auxillary grids
c
      call read_control_data(ctrlfile,"larval_rprof_grid_min", rprofmin)
      call read_control_data(ctrlfile,"larval_llarv_grid_min", 
     +                       llarv_len_min) 
      call read_control_data(ctrlfile,"larval_llarv_grid_max", 
     +                       llarv_len_max)
   
      write(*,*) "init_particle_state: larval interpolation grid = [",
     +           llarv_len_min,llarv_len_max,"]"
      write(*,*) "init_particle_state: larval profitability grid = [",
     +           rprofmin, " 1.0]"

      call read_control_data(ctrlfile,"larval_nprof_grid_points",nrprof)
      call read_control_data(ctrlfile,"larval_llarv_grid_points",nllarv)
      call read_control_data(ctrlfile,"ingestion_integral_sampling",
     +                       nforpt)

      write(*,*) "init_particle_state: larval size   sampling  =",nllarv
      write(*,*) "init_particle_state: profitability sampling  =",nrprof
      write(*,*) "init_particle_state: intake integral sampling=",nforpt

      allocate( rprof_grid(nrprof) )
      allocate( llarv_grid(nllarv) )
      do i=1,nrprof
         rprof_grid(i) = rprofmin + (1.0-rprofmin)*(i-1)/(nrprof-1)
      enddo
      dllarv = llarv_len_max - llarv_len_min 
      do i=1,nllarv
         llarv_grid(i) = llarv_len_min + dllarv*(i-1)/(nllarv-1)
      enddo
      write(*,155) "rprof_grid", rprofmin,      nrprof, 1.0
      write(*,155) "llarv_grid", llarv_len_min, nllarv, llarv_len_max
      write(*,*)   "init_particle_state: setup_foraging_window_db"
 155  format("init_particle_state:",a," =[",f9.5,",<",i4,">,",f9.5,"]")
      
      call setup_foraging_window_db(rprof_grid, llarv_grid)

c     grids llarv_grid and rprof_grid are copied by spline class
c     so they need not to be kept
      deallocate (rprof_grid)
      deallocate (llarv_grid)

      end subroutine 


      character*100 function get_particle_version()  
c     ---------------------------------------------------------------
c     ---------------------------------------------------------------
      get_particle_version = "sprat letcher model : $Rev: 147 $"
      end function



      subroutine close_particle_state() ! module operator
c     ---------------------------------------------------------------
c     Module clean up stuff here
c     (e.g. memory deallocation, MPI finalisation, final dumps etc.)
c     Only deallocate own data, i.e. data allocated in this module
c     ---------------------------------------------------------------
      write(*,*) "close_particle_state():"
      particle_counter               = 0

      call delete_spline(lplow_2Dspline)
      call delete_spline(lphigh_2Dspline)

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
      state%particleID = particle_counter    
      particle_counter = particle_counter + 1 ! module data

      call init_larval_physiology(state%larvstate)

      call set_state_hatched(state, space) 

      state%time_since_growth = 0 
      state%tempxdt           = 0    
      state%survival          = 1.0   

c     ----------------------------------------------------
      end subroutine 


      subroutine get_active_velocity(state, space, v_active)
c     --------------------------------------------------------------
c     Currently no migration pattern implemented (set v_active = 0 0 0)
c     If available, this subroutine should determine
c     v_active from (state, space) of the particle
c     This should really be in species properties
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
c     time period dt 
c
c     The subroutines below call a potentially deep hierachy
c     of other subroutines; to avoid duplicated interpolations of
c     the local environment, a flat container structure local_environment
c     is created which contains potential look-ups, and this 
c     container structure is passed around to increase flexibility
c     and avoid messy argument relaying 
c     --------------------------------------------------------------
      type(state_attributes), intent(inout)   :: state
      type(spatial_attributes), intent(inout) :: space
      real,intent(in)                         :: dt ! in seconds
c
      type(local_environment)                 :: local_env
      real                                    :: xyz(3)
      integer                                 :: status
      real                                    :: irate,tacc,tempavg
c     --------------------------------------------------------------
      if (dt<0) stop "negative dt values currently blocked"
      
      call get_tracer_position(space, xyz)
      call probe_local_environment(xyz,local_env)
      
      call add_routine_metabolic_costs(state%larvstate,local_env,dt) ! day/night
      if (local_env%light) then 
         call calculate_ingestion_rate(state%larvstate,local_env,irate)   ! visual predator
         call add_ingestion(state%larvstate, irate*dt, dt) ! currently just average, no Poisson       
      endif
      state%time_since_growth = state%time_since_growth + dt
      state%tempxdt           = state%tempxdt + dt*local_env%temp            

      if (state%time_since_growth > growth_increment_interval) then
         tacc    = state%time_since_growth
         tempavg = state%tempxdt/tacc
         call grow_larvae(state%larvstate, tempavg, tacc) ! growth in mass
         call update_larval_length(state) 
         call update_survival_chance(state, space, tacc) 
         state%time_since_growth = 0.0 ! reset time accumulator
         state%tempxdt           = 0.0 ! reset temp accumulator
      endif
      call clear_local_environment(local_env) ! in case memory were allocaetd
c     --------------------------------------------------------------
      end subroutine 



      subroutine delete_state_attributes(state) 
      type(state_attributes), intent(inout) :: state 
      call close_larval_physiology(state%larvstate)
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



      subroutine get_condition_number(state, phi)
c     ---------------------------------------------------
      type(state_attributes), intent(in)   :: state 
      real, intent(out)                    :: phi
      real                                 :: nomweight
c     ---------------------------------------------------
      call length_to_nominal_weight(state%larvstate%length, nomweight)
      phi = state%larvstate%weight / nomweight
      end subroutine get_condition_number
      



      subroutine update_larval_length(state)   
c     ---------------------------------------------------
c     Enforce that larval length can not decrease
c     ---------------------------------------------------
      type(state_attributes), intent(inout) :: state 
      real :: lnom
      call weight_to_nominal_length(state%larvstate%weight, lnom)
      if (lnom > state%larvstate%length) state%larvstate%length = lnom
      end subroutine update_larval_length


      


      subroutine set_state_hatched(state, space) 
c     ---------------------------------------------------
c     Initialize larvae in newly hatched state
c     make no assumptions about previous history
c     ---------------------------------------------------
      type(state_attributes), intent(inout)  :: state
      type(spatial_attributes),intent(inout) :: space
c     ---------------------------------------------------
      call set_tracer_mobility_free(space)       
      call set_larvae_hatched(state%larvstate,space)                         
      state%alive                = .true.
      
      end subroutine set_state_hatched



      subroutine set_state_dead(state, space)
c     ---------------------------------------------------
c     Set larval state to dead
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


c     ------ numerical integrals of continuous feeding intervals here ------


      subroutine trapez_integral(y, dy, h, intg, dintg_dx0, dintg_dx1) 
c     ------------------------------------------------------------------
c     Evaluate the trapez integral intg and the end points derivatives
c     (dintg_dx0, dintg_dx1) of the integrand y(x) on a regular sampling grid
c     x = [x0, x0+h, ..., x1](not in argument list - not needed explicitly) 
c     with regular grid spacing h. dy is the derivative of y, dydx(x), 
c     on same regular sampling grid.
c     
c     The end point derivatives (dintg_dx0, dintg_dx1) are consistent with
c     values obtained by numerical differentiation when using a 
c     finite sampling grid x. When a finite sampling grid is used, this
c     is slightly different from the analytic end point rule. When 
c     sampling grid x becomes infinitely dense, (dintg_dx0, dintg_dx1) 
c     should converge to the result of the analytic end point rule
c
c     Validated by numerical differentiation    
c     ------------------------------------------------------------------
      real,intent(in)  :: y(:)
      real,intent(in)  :: dy(:)
      real,intent(in)  :: h
      real,intent(out) :: intg, dintg_dx0, dintg_dx1
      integer          :: n,i
      real             :: w(size(y)),dhdx0,dhdx1,dsdx0,dsdx1
c     ------------------------------------------------------------------
      n     = size(y)
      w     = 1.0
      w(1)  = 0.5
      w(n)  = 0.5
      dhdx0 = -1.0/(n-1)/h
      dhdx1 =  1.0/(n-1)/h
c
      intg      = 0.0 
      dintg_dx0 = 0.0
      dintg_dx1 = 0.0
      do i=1,n
         dsdx0     = 1.0*(n-i)/(n-1)
         dsdx1     = 1.0*(i-1)/(n-1)
         intg      = intg + w(i)*y(i)
         dintg_dx0 = dintg_dx0 + w(i)*(dy(i)*dsdx0 + y(i)*dhdx0)
         dintg_dx1 = dintg_dx1 + w(i)*(dy(i)*dsdx1 + y(i)*dhdx1)
      enddo
      intg      = intg*h
      dintg_dx0 = dintg_dx0*h
      dintg_dx1 = dintg_dx1*h
      end subroutine trapez_integral



      subroutine eval_gain(llarv,rprof,local_env, g,dgdp,d2gdp2)
c     ------------------------------------------------------------------
c     Evaluate the energetic gain g [myg/sec] that a larvae of length llarv
c     have when feeding at profitability level rprof
c     Also evaluate profitability derivatives dgdp,d2gdp2
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
      real,intent(out)                   :: d2gdp2
c
      integer :: i
      real    :: lp_intv
      real    :: lp(nforpt)  ! automatic array to represent integration grids
      real    :: y(nforpt)   ! automatic array to represent ingestion integrands
      real    :: dy(nforpt)  ! automatic array to represent integrand derivatives

c     --- 2D spline derivative operator handlers ---
      integer,parameter :: deriv00(2) = (/0,0/)  ! just value 
      integer,parameter :: deriv01(2) = (/0,1/)  ! d/dp
      integer,parameter :: deriv02(2) = (/0,2/)  ! d2/dp2

      real    :: lp0, lp1, h
      real    :: ing, ding_dlp0, ding_dlp1
      real    :: xht, dxht_dlp0, dxht_dlp1
      real    :: ding_dp, dxht_dp, d2ing_dp2, d2xht_dp2, b, db_dp 
      real    :: dlp0_dp, d2lp0_dp2, dlp1_dp, d2lp1_dp2 

c     ------------------------------------------------------------------ 

c
c     --- retreieve foraging window [lp0,lp1] for this (llarv,rprof)
c
      call evaluate_spline(lplow_2Dspline, llarv,rprof,deriv00,lp0) 
      call evaluate_spline(lphigh_2Dspline,llarv,rprof,deriv00,lp1) 
      lp0     = lp0*llarv ! unscale 
      lp1     = lp1*llarv ! unscale 
c
c     --- setup integration sampling grid
c
      lp_intv = lp1-lp0
      h       = lp_intv/(nforpt-1)
      do i=1,nforpt
         lp(i) = lp0 + lp_intv*(i-1)/(nforpt-1)
      enddo
c
c     --- setup ingestion integrand and evaluate ingestion integral
c         along with derivatives
c      
      do i=1,nforpt
         y(i)  = enc_rate(lp(i), llarv, local_env, value)
     +           * epsil(lp(i), llarv, value)
         dy(i) = enc_rate(lp(i), llarv, local_env, slope)
     +          * epsil(lp(i), llarv, value) 
     +          + enc_rate(lp(i), llarv, local_env, value)
     +          * epsil(lp(i), llarv, slope)             
      enddo
      call trapez_integral(y, dy, h, ing, ding_dlp0, ding_dlp1) 
c
c     --- setup handling time (per unit time) integrand and evaluate 
c         handling time integral along with derivatives
c      
      do i=1,nforpt
          y(i) = enc_rate(lp(i), llarv, local_env, value)
     +          * HT(lp(i), llarv, value)
          dy(i) = enc_rate(lp(i), llarv, local_env, slope)
     +          * HT(lp(i), llarv, value) 
     +          + enc_rate(lp(i), llarv, local_env, value)
     +          * HT(lp(i), llarv, slope)        
      enddo
      call trapez_integral(y, dy, h, xht, dxht_dlp0, dxht_dlp1)
      xht = xht + 1.0  ! per unit time (does not affect derivatives)

c
c     --- compute end point derivatives consistently
c         with the numerical integration scheme using chain rule
c       
     
      call evaluate_spline(lplow_2Dspline, llarv,rprof,deriv01,dlp0_dp)    
      call evaluate_spline(lplow_2Dspline, llarv,rprof,deriv02,
     +                     d2lp0_dp2)  
      call evaluate_spline(lphigh_2Dspline,llarv,rprof,deriv01,dlp1_dp)   
      call evaluate_spline(lphigh_2Dspline,llarv,rprof,deriv02,
     +                     d2lp1_dp2) 

c
      ding_dp   = ding_dlp0*dlp0_dp + ding_dlp1*dlp1_dp
      dxht_dp   = dxht_dlp0*dlp0_dp + dxht_dlp1*dlp1_dp
      d2ing_dp2 =   ding_dlp0*(dlp0_dp**2 + d2lp0_dp2) 
     +            + ding_dlp1*(dlp1_dp**2 + d2lp1_dp2) 
      d2xht_dp2 =   dxht_dlp0*(dlp0_dp**2 + d2lp0_dp2) 
     +            + dxht_dlp1*(dlp1_dp**2 + d2lp1_dp2) 

      b     = ding_dp*xht - ing*dxht_dp
      db_dp = d2ing_dp2*xht - ing*d2xht_dp2

c
c     --- finally collect gain g and derivates (dgdp,d2gdp2) from
c         computed integrals using chain rule
c     
      g      = ing/xht
      dgdp   = b/xht**2
      d2gdp2 = (db_dp*xht**2 - 2.0*b*xht*dxht_dp)/xht**4


      end subroutine eval_gain



      subroutine calculate_ingestion_rate(self, local_env, irate) 
c     ------------------------------------------------------------------
c     Given the food biomass density zbiomass [g/l] calculate 
c     the intake rate irate [myg/sec] based on optimal foraging theory
c     for the larvae represented by self
c     This algorithm localizes the feeding profitability level corresponding
c     to optimal gain, using splines of continuous food selection intervals
c     Using a bracketed Newton process
c
c     Tests: 
c         1) make a plot of gain vs 0.1 < profitability < 1
c         2) check consistency of g,dgdp,d2gdp2 in eval_gain
c     ------------------------------------------------------------------
      type(larval_physiology),intent(in) :: self
      type(local_environment),intent(in) :: local_env ! ambient conditions
      REAL,intent(out)                   :: irate
c     ------------------------------------------------------------------
      irate = 99. ! silent compiler warnings
      end subroutine



      real function enc_rate(lp,llarv, local_env,deriv)
c     --------------------------------------
c     Evaluate the encounter rate density [unit == individuals/mm/sec] 
c     of prey at length lp [mm] of a larval predator of length llarv [mm]
c     The encounter rate of prey in a length interval [lp - 0.5*dL; lp + 0.5*dL]
c     corresponds to dN = enc_rate*dL  [unit == individuals/sec]
c     Encounter rate density is evaluated as clearence volume [unit == mm3/sec]
c     times density of prey [unit == individuals/mm3/mm] 
c
c     Applies spectrum_prey_density (from module prey_community)
c     and SV (from module larval_properties)      
c
c     deriv = false: return enc_rate   
c     deriv = true : return (d/dlp) enc_rate   
c     --------------------------------------
      REAL,intent(in)                    :: lp        ! prey length [mm]
      REAL,intent(in)                    :: llarv     ! larval length [mm]
      type(local_environment),intent(in) :: local_env ! ambient conditions
      logical,intent(in)                 :: deriv
c     --------------------------------------
      if (deriv) then
         enc_rate = 
     +         spectrum_prey_density(lp,local_env,slope)
     +         * SV(lp,llarv,local_env,value)
     +       + spectrum_prey_density(lp,local_env,value)
     +         * SV(lp,llarv,local_env,slope) 
      else
         enc_rate = 
     +         spectrum_prey_density(lp,local_env,value)
     +         * SV(lp,llarv,local_env,value)

      endif      
      end function enc_rate



      real function epsil(lp,llarv,deriv)
c     --------------------------------------
c     Evaluate prey value * capture success  [myg] 
c     of prey length lp
c     deriv = false: return epsil   
c     deriv = true : return (d/dlp) epsil   
c     --------------------------------------
      REAL,intent(in)                    :: lp        ! prey length [mm]
      REAL,intent(in)                    :: llarv     ! larval length [mm]
      logical,intent(in)                 :: deriv
c     --------------------------------------
      if (deriv) then
         epsil = prey_mass(lp,slope)*CS(lp,llarv,value)
     +         + prey_mass(lp,value)*CS(lp,llarv,slope)
      else
         epsil = prey_mass(lp,value)*CS(lp,llarv,value)
      endif      
      end function epsil


      
      real function rank(lp,llarv,deriv)
c     --------------------------------------
c     Evaluate            
c                       value * capture success   [myg] 
c         prey rank =  ---------------------------------
c                           handling time      [seconds]
c   
c     of prey length lp
c     deriv = false: return rank
c     deriv = true : return (d/dlp) rank
c     --------------------------------------
      REAL,intent(in)                    :: lp        ! prey length [mm]
      REAL,intent(in)                    :: llarv     ! larval length [mm]
c     type(local_environment),intent(in) :: local_env ! ambient conditions
      logical,intent(in)                 :: deriv
c     --------------------------------------
      if (deriv) then
         rank = (epsil(lp,llarv,slope)* HT(lp,llarv,value)
     +         - epsil(lp,llarv,value)* HT(lp,llarv,slope))
     +         / HT(lp,llarv,value)**2
      else
         rank =  epsil(lp,llarv,value)
     +         / HT(lp,llarv,value)
      endif      
      end function rank


      
      subroutine find_rank_maximum(llarv,lp,rmax)
c     --------------------------------------
c     Locate prey length lp corresponding to the 
c     prey rank maximum rmax for a larvae 
c     with larval length llarv using rank slope bisection
c
c     Assert lp_min_search_value < lp < llarv to establish a search bracket
c     --------------------------------------
      REAL,intent(in)    :: llarv ! larval length [mm]
      REAL,intent(out)   :: lp    ! prey length [mm]
      REAL,intent(out)   :: rmax  ! rank maximum [myg/second]
c
      REAL,parameter     :: lp_resolution = 1.e-6 ! [mm]
   
      REAL               :: lpleft,lpright,r
      integer            :: nsteps, ist
c     -------------------------------------
      lpleft  = lp_min_search_value  ! assert (d/dlp) rank (lleft)  > 0
      lpright = llarv                ! assert (d/dlp) rank (lright) < 0
c     --- check assertions ---
      if (rank(lpleft,llarv,slope) < 0) then
         stop "find_rank_maximum: invalid left bracket"
      endif
      if (rank(lpright,llarv,slope) > 0) then
         stop "find_rank_maximum: invalid right bracket"
      endif
      nsteps = int(log((lpright-lpleft)/lp_resolution) / log(2.))
      do ist=1,nsteps
         lp = 0.5*(lpleft + lpright)
         if (rank(lp,llarv,slope) > 0) then
            lpleft  = lp
         else
            lpright = lp
         endif
      enddo
      lp   = 0.5*(lpleft + lpright) ! final estimator
      rmax = rank(lp,llarv,value)
      end subroutine find_rank_maximum



      subroutine find_foraging_window(llarv,rprof,lplow,lphigh)
c     -------------------------------------------------------
c     Locate the larval foraging prey size interval 
c     window [lplow,lphigh] where rprof < rank(lp)/max_rank < 1
c     Can be accelerated using a Newton search process
c     Assume optimal foraging window is a contiguous interval 
c     Assume optimal foraging window is within the
c     search interval [lp_min_search_value, lpmaxfactor*llarv]
c     -------------------------------------------------------- 
      REAL,intent(in)    :: llarv ! larval length [mm]
      REAL,intent(in)    :: rprof ! minimal relative profitability [myg/second]
      REAL,intent(out)   :: lplow, lphigh ! foraging window [mm]
      
      REAL               :: lp_rmax,rmax,dlp,minrank
      REAL               :: lpleft,lpright
      REAL,parameter     :: lpmaxfactor = 1.0 ! max prey/pred size ratio 
      REAL               :: lplow_limit, lphigh_limit 
c     -------------------------------------------------------- 
c
c     1) locate rank maximum as delimiter for locating lplow,lphigh
c        Capture singular limit rprof=1, where lplow=lphigh=lp_rmax
      call find_rank_maximum(llarv, lp_rmax, rmax)
      if (abs(rprof-1.0)<1.e-4) then
         lplow  = lp_rmax
         lphigh = lp_rmax
         return
      endif
c
c     2) else locate lplow/lphigh
c  
      minrank = rmax*rprof
      lplow_limit  = lp_min_search_value
      lphigh_limit = lpmaxfactor*llarv
      call newton_search(lplow_limit, lp_rmax,      minrank, lplow)
      call newton_search(lp_rmax,     lphigh_limit, minrank, lphigh)

      contains

      subroutine newton_search(x0,x1,r0,x)
c     ------------------------------------------------
c     solve rank(x,llarv,value)=r0 in [x0,x1] in local scope
c     ------------------------------------------------
      real,intent(inout) :: x0,x1 ! input solution bracket
      real,intent(in)    :: r0    ! function offset
      real,intent(out)   :: x     ! solution 
      real               :: s0,s1,sm,xm,f,df,dx
      REAL,parameter     :: dxlimit = 1.0e-7
      integer            :: iter
      integer,parameter  :: max_iter = 20 ! suspend newton search, if exceeded
c     --------
      iter = 1
      x  = 0.5*(x0 + x1)
      s0 = sign(1.0, rank(x0,llarv,value)-r0)
      s1 = sign(1.0, rank(x1,llarv,value)-r0)
      dx = x1-x0 ! an upper bound 
      do while (abs(dx) > dxlimit) ! Newton loop     
         f  = rank(x,llarv,value)-r0
         df = rank(x,llarv,slope)
         dx = -f/df
         x  = x+dx
         if ((x<x0).or.(x>x1).or.(iter>max_iter)) then ! bisection fall back
            xm  = 0.5*(x0 + x1)
            sm  = sign(1.0, rank(xm,llarv,value)-r0)
            dx  = x1-x0 ! an upper bound
            if (s0*sm > 0) then
               x0 = xm
               s0 = sm
            else
               x1 = xm
               s1 = sm
            endif
            x  = 0.5*(x0 + x1) ! new mid point
         endif    ! bisection 
c         write(*,*) x0,x1,x,dx
         iter = iter+1
      enddo ! while     
      end subroutine newton_search ! local function to find_foraging_window
      end subroutine find_foraging_window



      subroutine setup_foraging_window_db(rprof_grid, llarv_grid)
c     ----------------------------------------------------------------
c     Setup 2D spline for lookup of foraging window, corresponding
c     to a given preying larval size and relative profitability.
c     The predation window is scaled to the preying larval size to 
c     condition the spline. Result is stored in module data 
c     lplow_2Dspline and lphigh_2Dspline
c 
c     If the relative profitability (rprof) is 1 (within numerical precision)
c     apply the exact limit to avoid numerical problem in solving 
c     near-singular problem. In the future, we may apply the
c     analytic expansion around rprof = 1
c     ----------------------------------------------------------------
      real, intent(in) :: rprof_grid(:)  ! grid of relative profitabilities
      real, intent(in) :: llarv_grid(:)  ! grid of preying larvae lengths
      integer          :: ip,il,nrprof,nllarv
      real             :: rprof, llarv, lplow,lphigh
      real             :: lp_rmax, rmax
      real,parameter   :: rprof_tol = 1.e-4
      real             :: lower(size(llarv_grid),size(rprof_grid))
      real             :: upper(size(llarv_grid),size(rprof_grid))
c     -----------------------------------------------------
      nllarv = size(llarv_grid)
      nrprof = size(rprof_grid)
      
      do il = 1, nllarv
         llarv = llarv_grid(il)
         do ip = 1, nrprof
            rprof = rprof_grid(ip) ! 
            if ((rprof < -rprof_tol).or.(rprof > 1.0+rprof_tol)) then
               write(*,*) "setup_foraging_window_db: invalid rprof", 
     +                     rprof
               stop
            endif
c           ---- define foraging window: [lplow,lphigh] ----       
            if (abs(rprof-1.0)<rprof_tol) then
               call find_rank_maximum(llarv, lp_rmax, rmax)
               lplow  = lp_rmax ! window collapses
               lphigh = lp_rmax ! window collapses
            else
               call find_foraging_window(llarv,rprof,lplow,lphigh)
            endif
            write(*,*) llarv,rprof, ":", lplow/llarv,lphigh/llarv
            lower(il,ip) = lplow/llarv
            upper(il,ip) = lphigh/llarv
         enddo
      enddo
c
c     ..... lplow_2Dspline/lphigh_2Dspline are module data .....
c           grids llarv_grid and rprof_grid are copied by spline class
c           so they need not to be kept
c
      call setup_spline(lplow_2Dspline, llarv_grid,rprof_grid,lower)
      call setup_spline(lphigh_2Dspline,llarv_grid,rprof_grid,upper)

      end subroutine setup_foraging_window_db
      



      subroutine update_survival_chance(state, space, dt)
c     ---------------------------------------------------
c     Update the survival counter state%survival for this organism
c     corresponding to time interval dt
c     Currently impose condition number hard limit
c     ---------------------------------------------------
      type(state_attributes),intent(inout)   :: state
      type(spatial_attributes),intent(inout) :: space
      real,intent(in)                        :: dt
      real                                   :: phi
c     ---------------------------------------------------
c     state%survival = ...

      call get_condition_number(state, phi)
      if (phi < min_condition_number) then
         call set_state_dead(state, space)  ! set state%survival to zero
      endif

      end subroutine update_survival_chance

      

c     subroutine calculate_predation_survival(state,space,pred_surv)
c     --------------------------------------------------- 
c     Letcher Eqs 19-22
c     At some point we may also overlay predator maps (space available)
c     ---------------------------------------------------
c     end subroutine calculate_predation_survival
      


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
