ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Passive particle biology provider 
c     ---------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
c     Particle_state interface of elementary passive tracer
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module particle_state

      use time_tools           ! import clock type
      use particle_tracking    ! space types/methods
      use output

      implicit none
      private     

c     -----------------------------------------------    
c     Class state_attributes contains all attributes 
c     related to the state of a particle beyond 
c     spatial aspects
c     Note: F90 does not allow an empty class, therefore
c     add a dummy variable
c     -----------------------------------------------
      type state_attributes
      private
        integer :: tracerID
      end type
      public :: state_attributes


      public :: init_particle_state  ! module operator
      public :: close_particle_state ! module operator
      public :: init_state_attributes
      public :: get_active_velocity
      public :: update_particle_state
      public :: delete_state_attributes 
      public :: write_state_attributes
      
      interface get_property
        module procedure get_prop_state
      end interface
      public :: get_property
      public :: get_metadata_state

c     -----------------
c     Module parameters
c     -----------------
      integer :: tracerID
      
c     --------
      contains
c     --------
     
      subroutine init_particle_state()  ! module operator
c     Display version
      write (*,*) "Passive particle biology provider : $Rev$"
      tracerID=1
      end subroutine 


      subroutine close_particle_state() ! module operator
      end subroutine 
      


      subroutine init_state_attributes(state, space, time_dir,             
     +                                 initdata, emitboxID)
c     ----------------------------------------------------
      type(state_attributes),intent(out)     :: state
      type(spatial_attributes),intent(inout) :: space
      real,intent(in)                        :: time_dir            
      character*(*),intent(in)               :: initdata
      integer,intent(in)                     :: emitboxID
c     ----------------------------------------------------
      state%tracerID = tracerID
      tracerID = tracerID +1
      call set_tracer_mobility_free(space)  
c     
      end subroutine 


      subroutine get_active_velocity(state, space, v_active)
      type(state_attributes), intent(in)   :: state
      type(spatial_attributes), intent(in) :: space      
      real, intent(out)                    :: v_active(:)  
      v_active = 0.
      end subroutine 

      subroutine update_particle_state(state, space, time_step)
      type(state_attributes), intent(inout) :: state
      type(spatial_attributes), intent(in)  :: space
      real,intent(in)                       :: time_step
c     -----------------
     
      end subroutine 


      subroutine delete_state_attributes(state) 
      type(state_attributes), intent(inout) :: state 
      end subroutine 


      subroutine write_state_attributes(state)
      type(state_attributes), intent(in) :: state 
      end subroutine 


      subroutine get_prop_state(state,var,bucket,status)
c------------------------------------------------------------  
      type(state_attributes),intent(in) :: state
      type(variable),intent(in) :: var
      type(polytype), intent(out) :: bucket
      integer, intent(out) :: status
c------------------------------------------------------------  
      status=0  !Variable is present
      select case (get_name(var))
      case ("tracerID")
        call construct(bucket,"tracerID",state%tracerID)
      case default
        status=1   !Cannont find variable name
      end select
      end subroutine


      subroutine get_metadata_state(var_name,var,status)
c------------------------------------------------------------  
      character(*), intent(in) :: var_name
      type(variable),intent(out) :: var
      integer, intent(out) :: status
c------------------------------------------------------------  
      status=0 !Defaults to variable found
      select case (var_name)
       case ("tracerID")
         call construct(var,"tracerID","tracer ID number",
     +     units="-",fmt="(i6)",type="int")
      case default
        status=1  !Cannot find variable
      end select
      end subroutine get_metadata_state



      end module
