ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     -------------------------------------------------------------------------
c     Simplified behavioral module for Plaice
c
c     Partially based on Hufnagl etal 2013: JSR 84, p26ff
c     stage_mortality rates are estimated from figure in Wennhage 1999 (M in figures wrong?)
c     -------------------------------------------------------------------------    
c   
c     $Rev: 181 $
c     $LastChangedDate: 2011-01-12 00:16:47 +0100 (Wed, 12 Jan 2011) $
c     $LastChangedBy: asch $ 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module particle_state_bio ! NB, different module name

      use particle_tracking     ! space types/methods
      use physical_fields
       

      implicit none
      private    
      
c     ---------------------------------------------
c     This type contains only specific biological attributes
c     Generic aspects, like accumulated survival, origin etc
c     are maintained by embedding type
c     ---------------------------------------------
      type state_attributes
      private
          real    :: age      ! days since instantiation
          integer :: otype    ! ontogenetic type: egg:0, pelagic:1, demersal:2
          real    :: completeness ! of this stage, 0 < completeness < 1
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

c     --------------- ontogeny (offset 1 to allow array mapping ) -------------------
      integer, parameter :: pelagic_egg   = 1
      integer, parameter :: pelagic_larv  = 2
      integer, parameter :: demersal_larv = 3
      character(len=15), parameter :: type_meaning(4) = 
     +                                          (/"pelagic egg    ",
     +                                            "pelagic larvae ",
     +                                            "demersal larvae",
     +                                            "beyond repres  "/)

c
c     stage_time(T[oC]) [days] = exp(x_0 - T[oC]*x_Tcoeff)
c
      real, parameter :: stage_time_0(3)      = (/3.95, 5.00, 3.40/)
      real, parameter :: stage_time_Tcoeff(3) = (/0.15, 0.15, 0.15/)
      real, parameter :: stage_mortality(3)   = (/33.4, 16.5, 5.86/) ! mortality per year, adapted from Wennhahe 1999 
      real, parameter :: seconds_per_year     = 365.0*86400.0
c     ===============================================================
                                  contains
c     ===============================================================
     
      subroutine init_particle_state()  ! module operator
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
      state%age          = 0.0  ! days since instantiation
      state%otype        = pelagic_egg 
      state%completeness = 0.0
      end subroutine init_state_attributes



      subroutine get_active_velocity(state, space, v_active)
c     --------------------------------------------------------------
c     Currently no horizontal migration pattern implemented
c
c     eggs:            stokes law, with drho= (petitgas 2006), all water column
c
c     pelagic  larvae: passive (v=0)
c
c     demersal larvae: passive (v=0) vdown_avg ~ 4 mm/s (Hufhagl;Ryland, 1963)
c                      If vertical diffusivity is constant D, then the larval aggregation layer
c                      has a thickness of approx L = D/v
c
c
c     Stokes terminal velocity = vt = (2/9) drho g R^2 / mu
c     mu is the dynamic viscosity (kg /m*s) water @ 10 oC mu = 1.3/1000 Pa s
c     R ~ 1mm = 1e-3 m
c     rho_egg ~ 1.02 g/cm3 = 1020 kg/m3 (assume constant, independent of pressure)
c     vt = (1020-rho)*1.7 mm/s (positive down)
c
c     function rho_UNESCO(ss,tt,pp) is defiend in oceanography_providers/generic_elements/water_density.f
c     only link to water_density.o, no module import
c     In this context we apply an avearage, static water density when assessing water compressibility
c     because this is a shallow water case, the small bias from this being much less than uncertainties 
c     from biological processes 
c     --------------------------------------------------------------
      type(state_attributes), intent(in)   :: state
      type(spatial_attributes), intent(in) :: space      
      real, intent(out)                    :: v_active(:)  
      real, parameter                      :: patm       = 1.01325    ! [bar] applied constant atmospheric pressure
      real, parameter                      :: rhow0      = 1027.0     ! [kg/m**3], applied constant water density for compressibility effects
      real, parameter                      :: g_x_Pa2bar =  9.81/1.e5 ! g * conversion from Pa to bar
      real                                 :: pressure, rho, salt, temp
      real                                 :: xyz(3)
      integer                              :: istat
      real, external                       :: rho_UNESCO ! defiend in ooceanography_providers/generic_elements/water_density.f
c     --------------------------------------------------------------
      v_active = 0 ! default, other cases than below
c
c     --- set vertical component in specific cases ---
c
      if     (state%otype == pelagic_egg) then

         call get_tracer_position(space, xyz)
         call interpolate_temp(xyz,   temp,  istat)
         call interpolate_salty(xyz,  salt,  istat)
         pressure = patm + g_x_Pa2bar*xyz(3)*rhow0   ! zeroth order approx., unit == bar
         rho      = rho_UNESCO(salt,temp,pressure)  ! unit == kg/m**3
         v_active(3) = (1020-rho)*1.7/1000             ! [m/s], stokes terminal velocity  (positive down)

      elseif (state%otype == demersal_larv) then

         v_active(3) = 4.0/1000. ! avg down velocity (i.e. > 0) Hufhagl;Ryland, 1963)

      endif

      end subroutine get_active_velocity



      subroutine write_state_attributes(state)
c     --------------------------------------------------------------
      type(state_attributes), intent(in)   :: state
      integer :: ist
c     --------------------------------------------------------------
      write(*,*) "age since instantiation = ", state%age, "days" 
      ist = min(1 + demersal_larv, state%otype) ! update_state should not exceed this, but anyway ...
      write(*,*) "stage type              = ", type_meaning(ist)
      write(*,*) "stage completeness      = ", state%completeness

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
c     This update does not update beyond end of demersal larvae stage,
c     so that age = age_at_settlement (if it settles)
c     --------------------------------------------------------------
      type(state_attributes), intent(inout)   :: state
      type(spatial_attributes), intent(inout) :: space
      real,intent(in)                         :: dt             ! in seconds
      real,intent(out)                        :: mortality_rate ! per second
      logical,intent(out)                     :: die, next  
c
      real    :: stage_time, ds
      integer :: ist, istat
      real    :: xyz(3), temp
      real,parameter :: do_die = 1000.0 ! just quite big
c     --------------------------------------------------------------
c     freeze state, if beyond demersal larvae (missed settlement)
      if (state%otype > demersal_larv) then
         die            = .true.
         next           = .false.   ! missed settlement
         mortality_rate = do_die
         return
      endif

      call get_tracer_position(space, xyz)
      call interpolate_temp(xyz, temp,  istat)

      state%age  = state%age + dt/86400.           ! unit = days
      ist        = min(state%otype, demersal_larv) ! index in stage maps
      stage_time = exp(stage_time_0(ist) - 
     +                 stage_time_Tcoeff(ist) * temp) ! unit = days
      state%completeness = state%completeness + dt/86400./stage_time
      
      mortality_rate = 0.5*stage_mortality(state%otype)
c     ---- advance to next stage represented here, completeness > 1 ----
      if (state%completeness > 1) then
         state%completeness = state%completeness - 1   ! better for medium time steps
         state%otype = state%otype + 1                 ! assumes stages are mapped consequtively
      endif
      ist = min(state%otype, demersal_larv)            ! index in stage maps
      mortality_rate = mortality_rate + 0.5*stage_mortality(ist) ! average stage_mortality at transition point  
      mortality_rate = mortality_rate/seconds_per_year ! convert to mortality per second
c
c     --- set transition flags depending on stage
c
      if (state%otype > demersal_larv) then ! ends here first time settlement missed
         die            = .true.
         next           = .false.           ! missed settlement
         mortality_rate = do_die
      elseif (state%otype == demersal_larv) then
         die            = .false.
         next           = .true.
      else                                  !  pelagic egg/larvae: not ready for settlement yet
         die            = .false.
         next           = .false.   
      endif
c     ----------------------------------------
      end subroutine update_particle_state_


      end module
