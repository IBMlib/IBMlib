ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Constant Growth 
c     ---------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
c     Empty public interface of particle_state
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module particle_state

      use time_tools           ! import clock type
      use particle_tracking    ! space types/methods
      use output               ! access polytype for get_prop_state/get_metadata_state

      implicit none
      private     

c     -----------------------------------------------    
c     Class state_attributes contains all attributes 
c     related to the state of a particle beyond 
c     spatial aspects
c     -----------------------------------------------
      type state_attributes
      private
      real :: lenght ! unit = mm
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
      public :: get_particle_version


c     -------------------- module data --------------------  
      real, parameter :: growth_speed = 1./86400. ! mm/sec
      real, parameter :: init_lenght  = 5.0       ! mm
      

c     --------
      contains
c     --------
     
         
 
      subroutine init_particle_state()  ! module operator
      write (*,*) trim(get_particle_version())
      end subroutine 

    
      subroutine close_particle_state() ! module operator
      end subroutine 
      

      character*100 function get_particle_version()  
      get_particle_version = "Constant growth particle biology" //
     +     " provider : $Rev$"
      end function
      

      subroutine init_state_attributes(state, space, time_dir,             
     +                                 initdata, emitboxID)
c     ----------------------------------------------------
      type(state_attributes),intent(out)     :: state
      type(spatial_attributes),intent(inout) :: space
      real,intent(in)                        :: time_dir            
      character*(*),intent(in)               :: initdata
      integer,intent(in)                     :: emitboxID
c     ----------------------------------------------------
      state%lenght = init_lenght
      call set_tracer_mobility_free(space)
      
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
      state%lenght = state%lenght + growth_speed*time_step ! time_step is signed

      end subroutine 


      subroutine delete_state_attributes(state) 
      type(state_attributes), intent(inout) :: state 
      end subroutine 


      subroutine write_state_attributes(state)
      type(state_attributes), intent(in) :: state 
      write(*,*) "length = ", state%lenght, " mm"
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
      case ("lenght")
        call construct(bucket,"lenght",state%lenght)
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
         call construct(var,"lenght","lenght of organism",
     +     units="mm",fmt="(f12.4)",type="real")
      case default
        status=1  !Cannot find variable
      end select
      end subroutine get_metadata_state


      end module
