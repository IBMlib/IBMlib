ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Empty biology provider 
c     ---------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
c     Empty public interface of particle_state
c
c     This module may trigger compiler warnings because 
c     variables with intent(out) are not set
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module particle_state
      use time_tools           ! import clock type
      use particle_tracking    ! space types/methods
      use output               ! access polytype for get_prop_state/get_metadata_state

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

      interface get_property
        module procedure get_prop_state
      end interface

      contains
         
 
      subroutine init_particle_state()  ! module operator
c     ----------------------------------------------------
      end subroutine 

      subroutine close_particle_state() ! module operator
c     ----------------------------------------------------
      end subroutine 

      character*100 function get_particle_version()  
      end function


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

      

      subroutine get_prop_state(state,var,bucket,status)
c------------------------------------------------------------  
      type(state_attributes),intent(in) :: state
      type(variable),intent(in)         :: var
      type(polytype), intent(out)       :: bucket
      integer, intent(out)              :: status
c------------------------------------------------------------  
      end subroutine


      subroutine get_metadata_state(var_name,var,status)
c------------------------------------------------------------  
      character(*), intent(in)   :: var_name
      type(variable),intent(out) :: var
      integer, intent(out)       :: status
c------------------------------------------------------------  
      end subroutine get_metadata_state


      end module
