ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Provide base classes/methods efficient representation
c     of bioenergetics 
c
c     local_environment is class to cache ambient
c       conditions at a particular time step to
c       avoid duplicated interpolations or extensive 
c       error prone argument lists (containing ambient parameters)
c       Avoid storing it for each particle - we only
c       need the content when updating the state of
c       a particle
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module particle_state_base
      use physical_fields
      use time_tools

      implicit none            
      private     
   
c     --- simple container class (public scope)
      type local_environment   
        real    :: temp     ! degrees Celcius
        real    :: zbiomass ! bulk zoo plankton density [myg DW/liter]
        logical :: light
        integer :: julday   ! Julian day, for seasonal patterns
      end type


c     ============  define public scope for this module =============

      public :: local_environment
      public :: probe_local_environment
      public :: clear_local_environment

c     ===============================================================
                                  contains
c     ===============================================================

      subroutine probe_local_environment(pos,local_env)
c     ---------------------------------------------------
c     ---------------------------------------------------
      real, intent(in)                     :: pos(:)
      type(local_environment), intent(out) :: local_env
      integer :: status
      type(clock), pointer :: current_time
      real                                 :: zdum(1)
c     ---------------------------------------------------
      call interpolate_temp(pos, local_env%temp, status) 
      local_env%light = is_light(current_time, pos) ! current_time from physical_fields
      call interpolate_zooplankton (pos, zdum, status) ! zdum: expects a vector
      local_env%zbiomass = zdum(1)
      current_time => get_master_clock()
      call get_julian_day(current_time, local_env%julday)

      end subroutine probe_local_environment
      



      subroutine clear_local_environment(local_env)
c     ---------------------------------------------------
c     Destructor - in case local_environment contains allocataed memory
c     Currently void
c     ---------------------------------------------------
      type(local_environment), intent(inout) :: local_env
      end subroutine clear_local_environment

      end 
