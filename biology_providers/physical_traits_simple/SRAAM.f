ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     -------------------------------------------------------------------------
c     Minimalistic biological module adapted for the SRAAM model from traits.f
c
c     This module addresses 3 major traits : 
c           1) pelagic period length (dynamic/fixed)
c           2) bouyancy 
c           3) start time of pelagic period (handled by calling modules)
c
c     Dynamics:   
c             horizontal passive particle with fixed sink_speed (buoyancy = 0) OR 
c             fixed object density (buoyancy = 1) in the latter case, the Stokes velocity is applied 
c             completely passive particle obtained by sink_speed = 0 and buoyancy = 0
c
c     Settlement: 
c             1a) settle_period_start < age < settle_period_end OR
c             1b) dynamic PLD + fixed length settlement window
c             2) reflective BC along shore
c     -------------------------------------------------------------------------    
c   
c     $Rev:  $
c     $LastChangedDate:  $
c     $LastChangedBy: asch $ 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module particle_state_bio ! NB, different module name

      use time_tools           ! import clock type
      use particle_tracking    ! space types/methods
      use physical_fields
      use run_context, only: ctrlfile => simulation_file
      use input_parser
     
      implicit none
      private    
      
c     ---------------------------------------------
c     This type contains only specific biological attributes
c     Generic aspects, like accumulated survival, origin etc
c     are maintained by embedding type
c     ---------------------------------------------
      type state_attributes
c      private   ! SRAAM public
        real    :: age          ! since instantiation [days]
        real    :: time_looking ! since first devel > 1 [days]
        real    :: devel        ! PLD completion (0 < devel < 1), for settlement_dynamics == 1
c       ------- logging data -------
        real    :: degree_days  ! Integrate(temperature, {t,0,now} [degreeC*days]
        real    :: min_temp     ! Min temperature encountered [degreeC]
        real    :: max_temp     ! Max temperature encountered [degreeC]
        real    :: min_salt     ! Min salinity encountered [PSU]
        real    :: max_salt     ! Max salinity encountered [PSU]
        
      end type
      public :: state_attributes

c     ---- export generic particle_state interface ----
      public :: init_particle_state  ! module operator
      public :: close_particle_state ! module operator
      public :: init_state_attributes
      public :: get_active_velocity
      public :: write_state_attributes
      public :: update_particle_state_ ! extended argument list
      

c     =================================================
c     ========       module data section       ========
c     =================================================
      integer, parameter    :: verbose = 0    ! control output volume      

      integer    :: settlement_dynamics
c     ---- settlement_dynamics == 0 : fixed settlement window ----      
      real       :: settle_period_start  ! since instantiation [days]
      real       :: settle_period_end    ! since instantiation [days]
c     ---- settlement_dynamics == 1 : dynamic PLD (OConnors2007)----
c          log(pld) = beta0 + beta1*log(temp/temp_c) + beta2*log(temp/temp_c)**2   
      real       :: beta0  ! correcponding to PLD being in unit days and T in degC
      real       :: beta1  ! correcponding to PLD being in unit days and T in degC
      real       :: beta2  ! correcponding to PLD being in unit days and T in degC
      real       :: temp_c       ! [deg C]
      real       :: settlement_window !  fixed [days]
      
c     buoyancy = 0: additive fixed sinking speed, in addition to turbulent velocity (default = 0.0)
c     buoyancy = 1: additive terminal Stokes velocity, corresponding to constant object density (no default) 

      integer    :: buoyancy             ! buoyancy model, default = 0
      real       :: sink_speed           ! [m/s],    applies for buoyancy=0, positive down
      real       :: rho_particle         ! [kg/m**3] fixed particle density, applies for buoyancy=1
      real       :: radius_particle ! [m]       for Stokes law, applies for buoyancy=1

c     active = 0: no active swimming
c     active = 1: keep within depth_min and depth_max by vertical swimming
c     active = 2: diurnal motion (respecting depth_min and depth_max)
c     active>1 leads to buoyancy = 0 with sinking speed = 0
      integer    :: active           ! active model, default = 0
      real       :: swim_depth_min   ! [m] upper confinement level > 0 (active>0)
      real       :: swim_depth_max   ! [m] lower confinement level (or sea bed, if shallower) (active>0
      real       :: swim_speed_light ! [m/s] swim speed when light, positive down (active=2)
      real       :: swim_speed_dark  ! [m/s] swim speed when dark, positive down  (active=2)
      real, parameter :: levelling_speed = 0.001 ! [m/s] swim speed > 0 to enforce swim_depth_min < z < swim_depth_min
c     ===============================================================
                                  contains
c     ===============================================================
     
      subroutine init_particle_state()  ! module operator
c     ------------------------------------------------------------------------
c
c     Resolve settlement dynamics
c
      if (count_tags(ctrlfile, "settle_period_start")>0) then ! precedence, if provided
         settlement_dynamics = 0
         call read_control_data(ctrlfile,"settle_period_start", 
     +                       settle_period_start)
         call read_control_data(ctrlfile,"settle_period_end",   
     +        settle_period_end)

         write(*,*)   "Fixed settlement window"
         write(*,388) "settle_period_start",settle_period_start, "days"
         write(*,388) "settle_period_end"  ,settle_period_end,   "days"
      else                      ! no other options
         settlement_dynamics = 1
         call read_control_data(ctrlfile,"beta0",  beta0)
         call read_control_data(ctrlfile,"beta1",  beta1)
         call read_control_data(ctrlfile,"beta2",  beta2)
         call read_control_data(ctrlfile,"temp_c", temp_c)
         call read_control_data(ctrlfile,"settlement_window",
     +                          settlement_window)
         write(*,*)   "Dynamics settlement window"
         write(*,388) "beta0",  beta0, ""
         write(*,388) "beta1",  beta1, ""
         write(*,388) "beta2",  beta2, ""
         write(*,388) "temp_c", beta2, "degC"
         write(*,388) "settlement_window", settlement_window, "days"
      endif
c
c     Resolve active mode
c         active = 0 : no active behavior (default)
c         active > 0 : enforce a confinement
c         active = 2 : light-driven vertical swimming
c
      active           = 0      ! default
      swim_depth_min   = 0      ! assign default (surface)
      swim_depth_max   = 1.0e8  ! assign default (bottom)
      swim_speed_light = 0.0    ! assign default (no swimming)
      swim_speed_dark  = 0.0    ! assign default (no swimming)
c      
      if (count_tags(ctrlfile, "swim_depth_min")>0) then
         call read_control_data(ctrlfile,"swim_depth_min",
     +                                    swim_depth_min) 
         active = 1             ! check confinement
         write(*,370) "upper", swim_depth_min
      endif
      if (count_tags(ctrlfile, "swim_depth_max")>0) then
         call read_control_data(ctrlfile,"swim_depth_max",
     +                                    swim_depth_max) 
         active = 1             ! check confinement
         write(*,370) "lower", swim_depth_max
      endif
      if (count_tags(ctrlfile, "swim_speed_light")>0) then
         call read_control_data(ctrlfile,"swim_speed_light",
     +                                    swim_speed_light)
         active = 2             
         write(*,371) "light", swim_speed_light 
      endif   
      if (count_tags(ctrlfile, "swim_speed_dark")>0) then
         call read_control_data(ctrlfile,"swim_speed_dark",
     +                                    swim_speed_dark)
         active = 2             
         write(*,371) "dark", swim_speed_light 
      endif
      
      if (swim_depth_min > swim_depth_max) then
         write(*,*) "init_particle_state: inconsistent parameters"
         write(*,277) swim_depth_min, swim_depth_max
         stop 45
      endif 
c
c     Resolve buoyancy mode, if active<2
c
c     fixed sinking speed takes precedence, if present (ignore rho_particle/radius_particle if present)
c
      buoyancy   = 0    ! default
      sink_speed = 0.0  ! default
      if (active<2) then
         if (count_tags(ctrlfile, "sink_speed")>0) then ! precedence, if provided

            call read_control_data(ctrlfile,"sink_speed",sink_speed) 
            buoyancy = 0
            write(*,390)  sink_speed

         elseif (count_tags(ctrlfile, "rho_particle")>0) then 

            call read_control_data(ctrlfile,"rho_particle",rho_particle) 
            call read_control_data(ctrlfile,"radius_particle",
     +                                       radius_particle) ! must be present, if rho_particle specified, no default
            buoyancy = 1
            write(*,391) rho_particle, radius_particle

         else                   ! default passive particle

            write(*,392) "particles are neutrally buoyant"

         endif
      endif ! active<2
 277  format("swim_depth_min =", f8.2, "> swim_depth_max =", f8.2)
      
 370  format(a, " confinement level = ", f8.2, " m")
 371  format("Swim speed at ", a, " = ", f8.3, " m/s")     
 388  format(a20, " = ", f8.3, 1x, a)
 390  format("particles have constant additive sinking speed = ",
     +       f8.3, " m/s")
 391  format("particles has constant density = ", 
     +       f8.3, " kg/m**3 and radius = ", e12.3, " m")
 392  format("particles are neutrally buoyant")

      end subroutine init_particle_state


      subroutine close_particle_state()  ! module operator
      end subroutine close_particle_state



      subroutine init_state_attributes(state, space, time_dir,             
     +                                 initdata, emitboxID)
c     ---------------------------------------------------------
      type(state_attributes),intent(out)     :: state
      type(spatial_attributes),intent(inout) :: space
      real,intent(in)                        :: time_dir            
      character*(*),intent(in)               :: initdata
      integer,intent(in)                     :: emitboxID
c     ---------------------------------------------------------
c     maintain reflective shore BC             ! FTH adaptation
c     ---------------------------------------------------------
      state%age           = 0.0
      state%time_looking  = 0.0
      state%devel         = 0.0 ! only updated for settlement_dynamics == 1
      state%degree_days   = 0.0
      state%min_temp      =  1.0e20 
      state%max_temp      = -1.0e20
      state%min_salt      =  1.0e20 
      state%max_salt      = -1.0e20
c      call set_shore_BC(space, BC_sticky)     ! FTH adaptation

      end subroutine init_state_attributes




      subroutine get_active_velocity(state, space, v_active)
c     --------------------------------------------------------------
c     Implement simple vertical velocities
c     buoyancy = 0: additive fixed sinking speed, in addition to turbulent velocity 
c     buoyancy = 1: additive terminal Stokes velocity, corresponding to constant object density 
c     --------------------------------------------------------------
      type(state_attributes), intent(in)   :: state
      type(spatial_attributes), intent(in) :: space      
      real, intent(out)                    :: v_active(:)

      real, parameter                      :: patm       = 1.01325    ! [bar] applied constant atmospheric pressure
      real, parameter                      :: rhow0      = 1027.0     ! [kg/m**3], applied constant water density for compressibility effects
      real, parameter                      :: g_x_Pa2bar = 9.81/1.e5 ! g * conversion from Pa to bar
      real, parameter                      :: g          = 9.81       ! [m/s2] gravitational constant
      real, parameter                      :: mydyn      = 1.002e-3   ! [kg/m/s] dynamic viscosity of water (at 20 oC)
      real                                 :: pressure, rho, salt, temp
      real                                 :: xyz(3)
      integer                              :: istat
      real, external                       :: rho_UNESCO ! defiend in ooceanography_providers/generic_elements/water_density.f
      type(clock),pointer                  :: now
c     --------------------------------------------------------------
      v_active(1:2) = 0             ! no horizontal component

      if      (buoyancy == 0) then  ! fixed sinking speed,

         v_active(3) = sink_speed   ! positive down

      elseif  (buoyancy == 1) then  ! fixed particle density

         call get_tracer_position(space, xyz)
         call interpolate_temp(xyz,   temp,  istat)
         call interpolate_salty(xyz,  salt,  istat)
         pressure    = patm + g_x_Pa2bar*xyz(3)*rhow0   ! zeroth order approx., unit == bar
         rho         = rho_UNESCO(salt,temp,pressure)   ! unit == kg/m**3
         v_active(3) = 2*(rho_particle-rho)*g*radius_particle**2/9/mydyn  ! [m/s], Stokes terminal velocity law  (positive down)

      else

         write(*,*) "get_active_velocity: unsupported buoyancy=",
     +              buoyancy
         stop 653
      endif
c
c     ---- case of light-triggered swimming
c
      if (active == 2) then
         now => get_master_clock()
         call get_tracer_position(space, xyz)
         if ( is_light(now,xyz) ) then ! apply default zenith
            v_active(3) = swim_speed_light
         else
            v_active(3) = swim_speed_dark
         endif
      endif
      
c     ---- confinement enforcement has precedence:
c      
c     soft enforcement approach, so that a swimming speed
c     will be added until vertical particle position z obeys
c           swim_depth_min < z < swim_depth_min
c     Notice the advection/diffusion components may lead to (small) violation of confinement
c 
      if (active > 0) then      ! check and enforce confinement
         call get_tracer_position(space, xyz)
         if     (xyz(3) < swim_depth_min) then  ! too high in water column
            v_active(3) = levelling_speed       ! go down
         elseif (xyz(3) > swim_depth_max) then  ! too deep in water column
            v_active(3) = -levelling_speed      ! go up
         endif
      endif
 
      end subroutine get_active_velocity



      subroutine write_state_attributes(state)
c     --------------------------------------------------------------
      type(state_attributes), intent(in)   :: state
c     --------------------------------------------------------------
      write(*,*) "age since init = ", state%age  
c
      end subroutine write_state_attributes



      subroutine update_particle_state_(state, space, dt, 
     +                       mortality_rate, die, next)
c     --------------------------------------------------------------
c     This subroutine should integrate state forward for
c     time period dt (extended argument list)
c     Also set auxillary parameters:
c        mortality_rate: current net mortality rate [1/sec]
c        die           : if .true. kill this entity (avoid having mortality_rate=infinity)
c        next          : if .true. advance (or regress) to ontogenetic stage (depending on sign of dt)
c
c     Settlement conditions: 
c          1) settle_period_start < age < settle_period_end
c     --------------------------------------------------------------
      type(state_attributes), intent(inout)   :: state
      type(spatial_attributes), intent(inout) :: space
      real,intent(in)                         :: dt ! in seconds
      real,intent(out)                        :: mortality_rate
      logical,intent(out)                     :: die, next 
      real                                    :: xyz(3), pld, y
      real                                    :: salt, temp
      integer                                 :: istat                               
c     --------------------------------------------------------------
      state%age      = state%age + dt/86400.
      mortality_rate = 0.0
      call get_tracer_position(space,xyz)
c     call interpolate_wdepth(xyz,cwd,istat)
c      
c     ====== update log vars ======
c  
      call interpolate_temp(xyz,   temp,  istat)
      call interpolate_salty(xyz,  salt,  istat)
      state%degree_days = state%degree_days + temp*dt/86400.
      state%min_temp    = min(state%min_temp, temp)
      state%max_temp    = max(state%max_temp, temp)
      state%min_salt    = min(state%min_salt, salt)
      state%max_salt    = max(state%max_salt, salt)
c      
c     ====== update settlement ======
c     
      if (settlement_dynamics == 0) then
c     ------ Fixed settlement window ------        
         if (state%age < settle_period_start) then
            die            = .false.
            next           = .false.
         elseif ((state%age > settle_period_start).and.
     +           (state%age < settle_period_end)) then
            die            = .false.
            next           = .true. ! try settling 
         else                   ! missed settlement
            die            = .true.   
            next           = .false.   
         endif
      elseif (settlement_dynamics == 1) then 
c     ------ dynamic settlement window --
         y   = log(temp/temp_c) ! temp assessed under log var update
         pld = exp(beta0 + beta1*y + beta2*y**2)  ! unit days, pld>0
         state%devel = state%devel + dt/86400/pld
         if (state%devel > 1) then
            state%time_looking = state%time_looking + dt/86400.
            if (state%time_looking<settlement_window) then
               die            = .false.! try settling
               next           = .true. ! try settling
            else
               die            = .true.  ! missed settlement
               next           = .false. ! missed settlement
            endif
         else
            die            = .false.   ! not ready to settle
            next           = .false.   ! not ready to settle
         endif
      else ! we should not end here
         write(*,*) "unexpected settlement_dynamics:",
     +               settlement_dynamics
         stop 744
      endif
c     ----------------------------------------
      end subroutine update_particle_state_


      end module
