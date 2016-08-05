ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     -------------------------------------------------------------------------
c     Simplified behavioral module for Nephrops 
c
c     Mainly based on:
c        Phelps, Jack J.C.; Polton, Jeff A.; Souza, Alejandro J.; Robinson, Leonie A. 2015 
c        The influence of behaviour on larval dispersal in shelf sea gyres: Nephrops norvegicus in the Irish Sea. 
c        Marine Ecology Progress Series, 518. 177-191.
c
c     Notes:
c        * Temperature dependence on stage development ignored (fixed to 9.23 oC), since Phelps et al
c          does not include this consistently (total PLD always fixed to 60 days)
c     -------------------------------------------------------------------------    
c     $Rev: 181 $
c     $LastChangedDate: 2011-01-12 00:16:47 +0100 (Wed, 12 Jan 2011) $
c     $LastChangedBy: asch $ 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module particle_state_bio ! NB, different module name

      use particle_tracking     ! space types/methods
      use physical_fields
      use time_tools
 

      implicit none
      private    
      
c     ---------------------------------------------
c     This type contains only specific biological attributes
c     Generic aspects, like accumulated survival, origin etc
c     are maintained by embedding type
c     ---------------------------------------------
      type state_attributes
      private
          real  :: age      ! days since instantiation
          real  :: stage    ! ontogenetic stage: pelagic in 1-4, 4-5 searching habitat, 5+ too late to settle
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
      integer, parameter :: pelagic_stage1 = 1
      integer, parameter :: pelagic_stage2 = 2
      integer, parameter :: pelagic_stage3 = 3
      integer, parameter :: try_settling   = 4
      integer, parameter :: beyond_repres  = 5
      character(len=15), parameter :: type_meaning(5) = 
     +                                          (/"pelagic stage 1",
     +                                            "pelagic stage 2",
     +                                            "pelagic stage 3",
     +                                            "try settling   ",
     +                                            "beyond repres  "/)

c
      real, parameter :: stage_duration(5)    = (/16.10243, 
     +                                            20.71276, 
     +                                            23.21990,
     +                                             1.00000,
     +                                             1.0e16/) ! days

      real, parameter :: swim_velocity(5) = (/1.0, 2.0, 3.0, 3.0, 0.0/) ! max vertical, mm/sec
     
      real, parameter :: homing_depth = 40.0   ! meters

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
      state%stage        = 1.0*pelagic_stage1
      end subroutine init_state_attributes



      subroutine get_active_velocity(state, space, v_active)
c     --------------------------------------------------------------
c     Currently no horizontal migration pattern implemented
c     --------------------------------------------------------------
      type(state_attributes), intent(in)   :: state
      type(spatial_attributes), intent(in) :: space      
      real, intent(out)                    :: v_active(:)  
      real                                 :: xyz(3)
      integer                              :: istat, istage, isec
      real                                 :: vervelo, harvest
      real                                 :: tau, lambda, phase, alpha
      type(clock),pointer                  :: current_time
c     --------------------------------------------------------------
      v_active = 0 ! default, other cases than below
c
c     --- set vertical component in specific cases ---
c
      istage = min(floor(state%stage), beyond_repres)  ! resolve current stage 1..5
      tau    = state%stage - istage                    ! intra stage completion
      call random_number(harvest)                      ! 0 < harvest < 1
      v_active = harvest*swim_velocity(istage)/1000.   ! unit = m/s, velocity not signed yet (positive = downward)

c     ... resolve time-in-day phase ...
      current_time => get_master_clock()
      call get_second_in_day(current_time, isec)
      phase        = 2.0*3.14159*isec/86400.0

c     ... above/below homing depth 
      call get_tracer_position(space, xyz)
      if (xyz(3) < homing_depth) then
         alpha = 0.0
      else
         alpha = 0.2  ! absorb factor 0.2 from Phelps here
      endif 
c
c     ..... resolve lambda, based on current stage .....
c
c           lambda < 0.5 -> up tendency of vertical motion
c           lambda > 0.5 -> down tendency of vertical motion
c
      if      (istage == pelagic_stage1) then
         if (tau < 0.5) then
            lambda = tau*(1.0 + sin(phase) - alpha)
         else
            lambda = 0.5*(1.0 + sin(phase) - alpha)
         endif
      elseif (istage == pelagic_stage2) then
         lambda = 0.5*(1.0 + sin(phase) - alpha)
      elseif (istage == pelagic_stage3) then
         lambda = 0.5*(1.0 + sin(phase))
      elseif (istage == try_settling) then
         lambda = (1.0-tau)*(1.0 + sin(phase))
      else
         lambda = 0.5
      endif
      call random_number(harvest)             ! 0 < harvest < 1
      if (harvest > lambda) then
         v_active = -v_active    ! direct step upward
      endif
      

      end subroutine get_active_velocity



      subroutine write_state_attributes(state)
c     --------------------------------------------------------------
      type(state_attributes), intent(in)   :: state
      integer :: istage
      real    :: tau
c     --------------------------------------------------------------
      write(*,*) "age since instantiation = ", state%age, "days" 
      istage = min(floor(state%stage), beyond_repres) 
      tau    = state%stage - istage
      write(*,*) "stage type              = ", type_meaning(istage)
      write(*,*) "stage completeness      = ", tau

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
      integer                                 :: istat, istage
      real,parameter                          :: do_die = 1000.0 ! just quite big
c     --------------------------------------------------------------
c     freeze state, if beyond demersal larvae (missed settlement)
      if (state%stage > beyond_repres) then
         die            = .true.
         next           = .false.   ! missed settlement
         mortality_rate = do_die
         return
      endif

c      call get_tracer_position(space, xyz)
c      call interpolate_temp(xyz, temp,  istat)

      mortality_rate = 0.0                             ! no mortality in this implementation
      state%age  = state%age + dt/86400.               ! unit = days
      istage = min(floor(state%stage), beyond_repres)  ! resolve current stage 1..5
      
      state%stage = state%stage + dt/86400.0/stage_duration(istage)     
c
c     --- set transition flags depending on stage
c
      if (istage == beyond_repres) then      ! ends here first time settlement missed
         die            = .true.
         next           = .false.           ! missed settlement
         mortality_rate = do_die
      elseif (istage == try_settling) then   ! settle, if habitat OK
         die            = .false.
         next           = .true.
      else                                  !  continue pelagic life
         die            = .false.
         next           = .false.   
      endif
c     ----------------------------------------
      end subroutine update_particle_state_


      end module
