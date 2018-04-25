ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     -------------------------------------------------------------------------
c     Minimalistic Plastic module addressing 3 major traits :
c           1) Plastic fragmentation
c           2) Bouyancy
c           3) Biofouling
c
c     Dynamics:   
c             horizontal passive particle with fixed sink_speed (buoyancy = 0) OR 
c             fixed object density (buoyancy = 1) in the latter case, the Stokes velocity is applied 
c             completely passive particle obtained by sink_speed = 0 and buoyancy = 0
c
c     Settlement: 
c             1) settle_period_start < age < settle_period_end
c             2) set sticky BC along shore, and poll if ashore
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
c      private
          real    :: age            ! since instantiation [days]
          real    :: density        ! Density will  increase due to biofouling
          real    :: biofouling_r   ! Biofouling starts from 0 and grows exponentially until reaching carrying capacity?
          real    :: biofouling_v   !Biofouling volume
          real    :: P_volume
          real    :: nb_algae       ! number of algae on the surface
          real    :: ESD            ! Equivalent spherical diameter

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

      real       :: settle_period_start  ! since instantiation [days]
      real       :: settle_period_end    ! since instantiation [days]



c     buoyancy = 0: additive fixed sinking speed, in addition to turbulent velocity (default = 0.0)
c     buoyancy = 1: additive terminal Stokes velocity, corresponding to constant object density (no default) 

      integer    :: buoyancy             ! buoyancy model, default = 0
      real       :: sink_speed           ! [m/s],    applies for buoyancy=0, positive down
      real       :: rho_particle         ! [kg/m**3] fixed particle density, applies for buoyancy=1
      real       :: radius_particle      ! [m]       consider the radius of the DR for Dietrich Law
      real       :: length_particle      ! [m]       for Dietrich law
      real       :: biofouling_density(3)! will contain 3 types of density according to how long the plastic has been in the water
      real       :: PIV                  ! particle initial volume
      real       :: PIM                  ! Particle initial mass

c     ===============================================================
                                  contains
c     ===============================================================
     
      subroutine init_particle_state()  ! module operator
c     ------------------------------------------------------------------------
      call read_control_data(ctrlfile,"settle_period_start",
     +                       settle_period_start)
      call read_control_data(ctrlfile,"settle_period_end",
     +                       settle_period_end)
    
      write(*,388) "settle_period_start",settle_period_start, "days"
      write(*,388) "settle_period_end"  ,settle_period_end,   "days"
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
         call read_control_data(ctrlfile,"length_particle",
     +                                   length_particle)
         buoyancy = 1
         PIV      = 3.14*length_particle*radius_particle**2    !Particle initial volume
         PIM      =  PIV*rho_particle  !Particle initial mass
         write(*,391) rho_particle, radius_particle,
     +                            length_particle

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
      state%age           = 0.0
      state%biofouling_r  = 0.0
      state%biofouling_v  = 0.0
      state%density       = rho_particle !kg/m**3
      state%P_volume      =PIV
      state%nb_algae      = 0
      state%ESD           = 2*(PIV*3/(4*3.1416))**(1/3)! Equivalent spherical diameter


      call set_shore_BC(space, BC_sticky)
      call set_bottom_BC(space, BC_reflect)
      end subroutine init_state_attributes




      subroutine get_active_velocity(state, space, v_active)
c     --------------------------------------------------------------
c     Implement simple vertical velocities
c     buoyancy = 0: additive fixed sinking speed, in addition to turbulent velocity 
c     buoyancy = 1: additive terminal Dietrich velocity, corresponding to constant object density
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
      real                                 :: V_kin                ! kinematic viscosity of seawater
      real                                 ::W_star                ! dimentionmess settling velocity
      real                                 ::D_star                ! dimentionmess particle diameter
c      real                                 :: A ,B ,Ur ,Uw, U                  ! All variables to get the kinematik viscosity

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
c         v_active(3) = 2*(state%density-rho)*g*
c     +                 state%R_equivalent**2/9/mydyn  ! [m/s], Stokes terminal velocity law  (positive down) //// we dont want stokes
        V_kin       =  mydyn/rho
        D_star      = (state%density-rho)*9.81*
     +                state%ESD**3/(rho*V_kin**2)

             if (D_star  <  0.05) then
              W_star=1.74E-4 * D_star**2

             elseif (D_star<= 5E9) then
              W_star=exp(-3.76715 + 192944 * log(D_star)
     +                -0.09815*log(D_star)**2 - 0.00575*log(D_star)**3+
     +                 0.00056*log(D_star)**4)
             else
             write(*,*) "D_star out of range", D_star
             end if


        v_active(3)= ((state%density-rho)*9.81*W_star*V_kin/rho)**(1/3)   ! Dietrich

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
      write(78,*)  state%age,  state%biofouling_r,
     +             state%density,
     +             state%P_volume,  state%nb_algae
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
      real                                    :: xyz(3)
      integer                                 :: status
      real                                    :: Mp   ! Particule mass
      real                                    :: Vp !Particule volume
      real                                    :: biofilm_p
      real                                    :: biofouling_m
      real                                    :: rnd(5)
      integer, parameter                      :: nmax = 4
      real                                    :: size(2**nmax)
      integer                                 :: icut, n, i
      real                                    :: l1,l2
      real                                    :: Alg_v=0.000000001  !!! what units??? what range???

c     --------------------------------------------------------------
      state%age      = state%age + dt/86400.
      biofouling_density(1)=1040
      biofouling_density(2)=1100
      biofouling_density(3)=1300
      mortality_rate = 0.0
c      if (state%age < settle_period_start) then
c     die            = .false.
c      next           = .false.
c      elseif ((state%age > settle_period_start).and.
c     +        (state%age < settle_period_end)) then
c      die            = .false.
c      next           = .true.   ! try settling
c      else ! missed settlement
c      die            = .true.
c      next           = .false.
c      endif

      if (state%age < settle_period_start) then
      die            = .false.
      next           = .false.
      elseif ((state%age > settle_period_start).and.
     +        (state%age < settle_period_end)) then
      die            = .false.
      next           = .true.   ! try settling
      else ! missed settlement
      die            = .true.
      next           = .false.
      endif
      if (state%age <=1209600) then
        biofilm_p=biofouling_density(1)
      else if (state%age <= 3628800) then
        biofilm_p =biofouling_density(2)
      else
        biofilm_p=biofouling_density(3)
      endif
      state%nb_algae=state%nb_algae + dt/86400
      state%biofouling_v=state%nb_algae*Alg_v
      biofouling_m=state%biofouling_v*biofilm_p
      state%P_volume =PIV+state%biofouling_v
      Mp =PIM+biofouling_m
      state%density=Mp/state%P_volume

      write(*,*) state%density,Mp,state%P_volume
c      call get_tracer_position(space,xyz)
c      call interpolate_wdepth(xyz,cwd,status)
      

c     ----------------------------------------
      end subroutine update_particle_state_



      end module