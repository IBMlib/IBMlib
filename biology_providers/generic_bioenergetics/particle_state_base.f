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
        real    :: zbiomass ! bulk zoo plankton density [kg DW / m3]
        logical :: light
        integer :: julday   ! Julian day, for seasonal patterns
      end type


c     ============  define public scope for this module =============

      public :: local_environment
      public :: probe_local_environment
      public :: write_local_environment
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
      call check_stat("interpolate_temp",status)
      current_time => get_master_clock()
      local_env%light = is_light(current_time, pos)    ! current_time from physical_fields
      call interpolate_zooplankton (pos, zdum, status) ! zdum: expects a vector
      call check_stat("interpolate_zooplankton",status)
      local_env%zbiomass = zdum(1)                     ! unit = kg DW / m3
      call get_julian_day(current_time, local_env%julday)
c     -----------
      contains   ! accesses parent scope
c     -----------      
      subroutine check_stat(who,intg)
      character*(*),intent(in) :: who
      integer, intent(in)      :: intg
      if (intg>0) then
         write(*,211) who,intg     
         write(*,212) pos ! access parent scope
      endif
 211  format("warning:probe_local_environment:",a,
     +       "failed with status = ",i2)
 212  format("interpolation position = ",3f12.7)
      end subroutine check_stat
      end subroutine probe_local_environment
      

      subroutine write_local_environment(local_env)
c     ---------------------------------------------------
c     For debugging
c     ---------------------------------------------------
      type(local_environment), intent(inout) :: local_env
      write(*,692) "temperature", local_env%temp
      write(*,692) "zoo biomass", local_env%zbiomass
      write(*,694) "light",       local_env%light
      write(*,693) "julian day",  local_env%julday
        
 692  format("local_environment:", a10, "=", f12.7)  
 693  format("local_environment:", a10, "=", i12) 
 694  format("local_environment:", a10, "=", l12)      

      end subroutine write_local_environment



      subroutine clear_local_environment(local_env)
c     ---------------------------------------------------
c     Destructor - in case local_environment contains allocataed memory
c     Currently void
c     ---------------------------------------------------
      type(local_environment), intent(inout) :: local_env
      end subroutine clear_local_environment

      end 
