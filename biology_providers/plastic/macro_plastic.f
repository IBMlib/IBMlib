ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     -----------------------------------------------
c     Macroplastic template
c     -----------------------------------------------
c     
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
c
c     References:
c     Yoon2010 : Yoon JH, Kawano S, Igawa S
c                Modeling of marine litter drift and beaching in the Japan Sea.
c                Mar Pollut Bull. 2010 Mar; 60(3):448-63
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module particle_state

      use time_tools           ! import clock type
      use particle_tracking    ! space types/methods
      use output               ! access polytype for get_prop_state/get_metadata_state
      use input_parser       
      use run_context, only: ctrlfile => simulation_file
      use physical_fields
      
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
      real    :: age            ! [days]
      real    :: sink_rate      ! [1/sec]
      real    :: k_transfer     ! [dimless] 10m wind velocity transfer; k*sqrt(A/W) in Yoon2010
      end type
      public :: state_attributes


      public :: init_particle_state  ! module operator
      public :: close_particle_state ! module operator
      public :: init_state_attributes
      public :: get_active_velocity
      public :: update_particle_state
      public :: delete_state_attributes 
      public :: write_state_attributes

      public :: get_sink_rate        ! extended particle state API
      
c     --- I/O + version hand shake ---
      interface get_property
        module procedure get_prop_state
      end interface
      public :: get_particle_version
      public :: get_property
      public :: get_metadata_state


c     -----------------
c     Module parameters
c     -----------------
     
      real :: sink_rate   ! unit = per second
      real :: k_transfer  ! unitless  
      
c     --------
      contains
c     --------
     
      subroutine init_particle_state()  ! module operator
c     -------------------------------------------------------------    
      call read_control_data(ctrlfile,"sink_rate",  sink_rate)
      write(*,399) sink_rate
      sink_rate = sink_rate/86400. ! per day -> per sec     
 399  format("init_particle_state: sink_rate = ",f12.7," per day")

      call read_control_data(ctrlfile,"k_transfer",  k_transfer)
      write(*,400) k_transfer
 400  format("init_particle_state: k_transfer = ",f12.7)
      
      end subroutine 


      character*100 function get_particle_version()  
      get_particle_version = "hacked biology" //
     +     " provider : $Rev$"
      end function


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
      state%age = 0
      state%sink_rate  = sink_rate  ! currently same for all
      state%k_transfer = k_transfer ! currently same for all
      call set_tracer_mobility(space,1,1,0)   ! only horizontal motion
c     
      end subroutine 


      subroutine get_active_velocity(state, space, v_active)
c     ----------------------------------------------------
c     Add wind drag veolcity
c     ----------------------------------------------------      
      type(state_attributes), intent(in)   :: state
      type(spatial_attributes), intent(in) :: space      
      real, intent(out)                    :: v_active(:)  
      real                                 :: uv10(2),xyz(3)
      integer                              :: status
c     ----------------------------------------------------
      call get_tracer_position(space, xyz)
      call interpolate_wind_10m(xyz, uv10, status)
      v_active(1:2) = uv10 * state%k_transfer ! Yoon2010 eq.3
      v_active(3)   = 0.0
c     ----------------------------------------------------      
      end subroutine
      

      subroutine update_particle_state(state, space, time_step)
      type(state_attributes), intent(inout) :: state
      type(spatial_attributes), intent(in)  :: space
      real,intent(in)                       :: time_step
c     -----------------
      state%age = state%age + time_step/86400.
      end subroutine 


      real function get_sink_rate(state)
c     -----------------------------------------------    
c     Extended API needed by DRRS
c     Maintain object encapsulation      
c     -----------------------------------------------      
      type(state_attributes), intent(in) :: state
      get_sink_rate = state%sink_rate
      end function get_sink_rate

      
      
      subroutine delete_state_attributes(state) 
      type(state_attributes), intent(inout) :: state 
      end subroutine 

      

      subroutine write_state_attributes(state)
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
