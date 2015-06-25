ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     -------------------------------------------------------------------------
c     Sea grass dispersal model: 
c         dynamics:   passive particle with fixed sink_speed
c         settlement: 
c             1) settle_period_start < age < settle_period_end
c             2) set sticky BC along shore, and poll if ashore
c     -------------------------------------------------------------------------    
c   
c     $Rev: 181 $
c     $LastChangedDate: 2011-01-12 00:16:47 +0100 (Wed, 12 Jan 2011) $
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
          real    :: age            ! since instantiation [days]
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
      real       :: settle_height        ! meters above bottom

c     buoyancy = 0: additive fixed sinking speed (default = 0, in addition to turbulent velocity)
c     buoyancy = 1: additive fixed density driven (corresponding to constant object density - precedence with fixed sinking unresolved)

      integer    :: buoyancy             ! buoyancy model
      real       :: sink_speed           ! [m/s] applies for buoyancy=0, positive down

c     ===============================================================
                                  contains
c     ===============================================================
     
      subroutine init_particle_state()  ! module operator
      call read_control_data(ctrlfile,"settle_period_start", 
     +                       settle_period_start)
      call read_control_data(ctrlfile,"settle_period_end",   
     +                       settle_period_end)
      call read_control_data(ctrlfile,"settle_height",   
     +                       settle_height )
      write(*,388) "settle_period_start",settle_period_start, "days"
      write(*,388) "settle_period_end"  ,settle_period_end,   "days"
      write(*,388) "settle_height",      settle_height,"m over bott"

      if (count_tags(ctrlfile, "sink_speed")>0) then ! precedence, if provided
         call read_control_data(ctrlfile,"sink_speed",sink_speed) 
         buoyancy = 0
         write(*,390) "particles has constant sinking speed", sink_speed, "m/s"
      else ! default passive particle
         sink_speed = 0.0
         buoyancy   = 0 
         write(*,390) "particles are neutrally buoyant"
      endif

 388  format(a20, " = ", f8.3, 1x, a)
 390  format(a40, " = ", f8.3, 1x, a)  
 392  format(a40)     
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
      state%age      = 0.0
      call set_shore_BC(space, BC_sticky) 

      end subroutine init_state_attributes




      subroutine get_active_velocity(state, space, v_active)
c     --------------------------------------------------------------
c     Currently only constant sinking speed implemented
c     --------------------------------------------------------------
      type(state_attributes), intent(in)   :: state
      type(spatial_attributes), intent(in) :: space      
      real, intent(out)                    :: v_active(:)  
c     --------------------------------------------------------------
      if (buoyancy /= 0) then
         write(*,*) "get_active_velocity: unsupported buoyancy=",
     +              buoyancy
         stop 653
      endif
      v_active(1:2) = 0.0
      v_active(3)   = sink_speed ! positive down
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
c          2a) distance to sea bed < settle_height OR
c          2b) particle is ashore
c     --------------------------------------------------------------
      type(state_attributes), intent(inout)   :: state
      type(spatial_attributes), intent(inout) :: space
      real,intent(in)                         :: dt ! in seconds
      real,intent(out)                        :: mortality_rate
      logical,intent(out)                     :: die, next 
      real                                    :: cwd, xyz(3), h
      integer                                 :: status
c     --------------------------------------------------------------
      state%age      = state%age + dt/86400.
      mortality_rate = 0.0
      call get_tracer_position(space,xyz)
      call interpolate_wdepth(xyz,cwd,status)
      h = cwd-xyz(3) ! h = height above bottom > 0
c
      if (state%age < settle_period_start) then
         die            = .false.
         next           = .false.
      elseif ((state%age > settle_period_start).and.
     +        (state%age < settle_period_end)) then
         die            = .false.
         if ((h < settle_height).or.is_ashore(space)) then
            next           = .true.   ! try settling
         else
            next           = .false.  ! too high, keep trying
         endif  
      else ! missed settlement
         die            = .true.   
         next           = .false.   
      endif
c     ----------------------------------------
      end subroutine update_particle_state_


      end module
