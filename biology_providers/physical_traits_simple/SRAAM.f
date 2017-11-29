ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     -------------------------------------------------------------------------
c     Minimalistic biological module adapted for the SRAAM model from traits.f
c
c     This module addresses 3 major traits : 
c           1) pelagic period length 
c           2) bouyancy 
c           3) start time of pelagic period (handled by calling modules)
c
c     Dynamics:   
c             horizontal passive particle with fixed sink_speed (buoyancy = 0) OR 
c             fixed object density (buoyancy = 1) in the latter case, the Stokes velocity is applied 
c             completely passive particle obtained by sink_speed = 0 and buoyancy = 0
c
c     Settlement: 
c             1) settle_period_start < age < settle_period_end
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
      private
        real    :: age         ! since instantiation [days]
        real    :: devel       ! PLD completion (0 < devel < 1), for settlement_dynamics == 1
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
      real       :: beta0 ! correcponding to PLD being in unit days
      real       :: beta1
      real       :: beta2
      real       :: temp_c ! [deg C]
      
c     buoyancy = 0: additive fixed sinking speed, in addition to turbulent velocity (default = 0.0)
c     buoyancy = 1: additive terminal Stokes velocity, corresponding to constant object density (no default) 

      integer    :: buoyancy             ! buoyancy model, default = 0
      real       :: sink_speed           ! [m/s],    applies for buoyancy=0, positive down
      real       :: rho_particle         ! [kg/m**3] fixed particle density, applies for buoyancy=1
      real       :: radius_particle      ! [m]       for Stokes law, applies for buoyancy=1
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
         write(*,*)   "Dynamics settlement window"
         write(*,388) "beta0",  beta0, ""
         write(*,388) "beta1",  beta1, ""
         write(*,388) "beta2",  beta2, ""
         write(*,388) "temp_c", beta2, "degC"
      endif
c
c     Resolve buoyancy mode
c
c     fixed sinking speed takes precedence, if present (ignore rho_particle/radius_particle if present)
c
      if (count_tags(ctrlfile, "sink_speed")>0) then ! precedence, if provided

         call read_control_data(ctrlfile,"sink_speed",sink_speed) 
         buoyancy = 0
         write(*,390)  sink_speed

      elseif (count_tags(ctrlfile, "rho_particle")>0) then 

         call read_control_data(ctrlfile,"rho_particle",rho_particle) 
         call read_control_data(ctrlfile,"radius_particle",
     +                                   radius_particle)  ! must be present, if rho_particle specified, no default
         buoyancy = 1
         write(*,391) rho_particle, radius_particle

      else ! default passive particle

         sink_speed = 0.0
         buoyancy   = 0 
         write(*,392) "particles are neutrally buoyant"

      endif

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
      state%age    = 0.0
      state%deve   = 0.0  ! only updated for settlement_dynamics == 1
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
      real                                    :: xyz(3), pld, temp, y
      integer                                 :: status                               
c     --------------------------------------------------------------
      state%age      = state%age + dt/86400.
      mortality_rate = 0.0
      call get_tracer_position(space,xyz)
c     call interpolate_wdepth(xyz,cwd,status)
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
      elseif (settlement_dynamics == 1) then Fixed settlement window
c     ------ dynamic settlement window --
         call interpolate_temperature(xyz,temp,status)
         y   = log(temp/temp_c)
         pld = exp(beta0 + beta1*y + beta2*y**2)  ! unit days
         state%devel += dt/86400/pld
         if (state%devel > 1) then
            die            = .false.
            next           = .true. ! try settling
            TOHERE
      else
         write(*,*) "unexpected settlement_dynamics:",
     +               settlement_dynamics
         stop 744
      endif
c     ----------------------------------------
      end subroutine update_particle_state_


      end module
