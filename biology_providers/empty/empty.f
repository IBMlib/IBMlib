ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Empty public interface of particle_state
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module particle_state
      use time_tools           ! import clock type
      use particle_tracking    ! space types/methods
      implicit none

c     -----------------------------------------------    
c     Class state_attributes contains all attributes 
c     related to the state of a particle beyond 
c     spatial aspects
c     -----------------------------------------------
      type state_attributes
      private
      logical :: dummy ! derived types can not be completely empty ...
      end type
      public :: state_attributes


      contains
         
 
      subroutine init_particle_state()  ! module operator
c     ----------------------------------------------------
      end subroutine 


      subroutine close_particle_state() ! module operator
c     ----------------------------------------------------
      end subroutine 
      

      subroutine init_state_attributes(state, space, time_dir,             
     +                                 initdata, emitboxID)
c     ----------------------------------------------------
      type(state_attributes),intent(out)     :: state
      type(spatial_attributes),intent(inout) :: space
      real,intent(in)                        :: time_dir            
      character*(*),intent(in)               :: initdata
      integer,intent(in)                     :: emitboxID
      end subroutine 


      subroutine get_active_velocity(state, space, v_active)
c     ----------------------------------------------------
      type(state_attributes), intent(in)   :: state
      type(spatial_attributes), intent(in) :: space      
      real, intent(out)                    :: v_active(:)  
      end subroutine 


      subroutine update_particle_state(state, space, time_step)
c     ----------------------------------------------------
      type(state_attributes), intent(inout) :: state
      type(spatial_attributes), intent(in)  :: space
      real,intent(in)                       :: time_step
      end subroutine 


      subroutine delete_state_attributes(state)
c     ----------------------------------------------------
      type(state_attributes), intent(inout) :: state 
      end subroutine 


      subroutine write_state_attributes(state)
c     ----------------------------------------------------
      type(state_attributes), intent(in) :: state 
      end subroutine 


      end module
