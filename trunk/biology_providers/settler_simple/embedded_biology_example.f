ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     -------------------------------------------------------------------------
c     Example of an embedded biology provider with a fixed limited settlement 
c     time window
c     Only represent the presettlement period (this module assumes organism
c     is unsettled)
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
          real    :: age            ! since instantiation
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

      real, parameter       :: settle_period_start = 0.2*86400.0
      real, parameter       :: settle_period_end   = 0.5*86400.0

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
      state%age      = 0.0
      end subroutine init_state_attributes



      subroutine get_active_velocity(state, space, v_active)
c     --------------------------------------------------------------
c     Currently no migration pattern implemented
c     --------------------------------------------------------------
      type(state_attributes), intent(in)   :: state
      type(spatial_attributes), intent(in) :: space      
      real, intent(out)                    :: v_active(:)  
c     --------------------------------------------------------------
      v_active = 0
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
c     --------------------------------------------------------------
      type(state_attributes), intent(inout)   :: state
      type(spatial_attributes), intent(inout) :: space
      real,intent(in)                         :: dt ! in seconds
      real,intent(out)                        :: mortality_rate
      logical,intent(out)                     :: die, next   
c     --------------------------------------------------------------
      state%age      = state%age + dt
      mortality_rate = 0.0
      if (state%age < settle_period_start) then
         die            = .false.
         next           = .false.
      elseif ((state%age > settle_period_start).and.
     +        (state%age < settle_period_end)) then
         die            = .false.
         next           = .true.   
      else ! missed settlement
         die            = .true.   
         next           = .false.   
      endif
c     ----------------------------------------
      end subroutine update_particle_state_


      end module
