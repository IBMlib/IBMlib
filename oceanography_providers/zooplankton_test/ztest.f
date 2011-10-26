ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Zooplankton access test module
c     ---------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
c     depth          = constant wdepth 
c     no coast line (all water, no land)
c     water temp     = seasonal trigonometric
c     turbulence     = 0
c     currents       = 0
c     all other interpolations = 0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module physical_fields
      use time_tools           ! import clock type
      use constants
      use geometry

      implicit none
      private     


      public :: init_physical_fields     
      public :: close_physical_fields
      public :: get_master_clock
      public :: set_master_clock
      public :: update_physical_fields
      public :: interpolate_turbulence
      public :: interpolate_turbulence_deriv
      public :: interpolate_currents
      public :: interpolate_zooplankton
      public :: interpolate_temp
      public :: interpolate_wdepth
      public :: is_wet    
      public :: is_land
      public :: horizontal_range_check
      public :: coast_line_intersection
      public :: get_pbi_version

c     -------------------- module data --------------------  
     
      

c     --- parameters for trigonometric seasonal temperature cycle ---
      real, parameter    :: avg_wtemp          = 10.5 ! yearly average - deg Celcius
      real, parameter    :: wtemp_amp          = 5.5  ! temperature oscillation, deg Celcius
      integer, parameter :: jday_with_mintemp  = 45   ! jday with lowest temperature in year


c     --- monthly averages of zooplankton abundance ---
c     Below are monthly averages over 1958-2007 of monthlty mean CPR count from 
c     Calanus_finmarchicus_C1.csv
c
      real, parameter    :: zoop_mon_avg(12) = 
     +    (/1.7842,  1.9752,  4.187, 
     +     24.3774, 26.7818, 18.122, 
     +     18.6544,  5.6966,  5.697, 
     +      8.2614,  4.4058,  2.2864/)
      
      >>>> please verify that units are kg DW/m3 <<<<


      real, parameter :: wdepth   = 40.0 ! water depth in meters
       
      

      type(clock), target  :: master_clock ! simulation master clock


c     --------
      contains
c     --------


c     ------------------------------------------    
      subroutine init_physical_fields(time)
      type(clock), intent(in),optional :: time 
      if (present(time)) master_clock = time
      write(*,*) trim(get_pbi_version()) 
      end subroutine 
c     ------------------------------------------ 
      character*100 function get_pbi_version()  
      get_pbi_version =  "Zooplankton pbi version: $Rev$"
      end function
c     ------------------------------------------ 
      subroutine close_physical_fields()
      end subroutine
c     ------------------------------------------ 
      function   get_master_clock()
      type(clock), pointer :: get_master_clock
      get_master_clock => master_clock
      end function 
c     ------------------------------------------ 
      subroutine set_master_clock(time)
      type(clock), intent(in) :: time
      master_clock = time ! local copy
      end subroutine 
c     ------------------------------------------ 
      subroutine update_physical_fields(time, dt)
      type(clock), intent(in),optional :: time
      integer, intent(in),optional     :: dt
      if (present(time)) then
         master_clock = time
      elseif (present(dt)) then
         call add_seconds_to_clock(master_clock, dt)
      endif
      end subroutine 
c     ------------------------------------------ 





      subroutine null_r3_interpolation(xyz, r3, status)
c     ------------------------------------------ 
c     flag sub sea bed quiries
c     ------------------------------------------ 
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
c     ------------------------------------------ 
      r3     = 0.
      if (is_wet(xyz)) then
         status = 0
      else
         status = 1
      endif
      end subroutine 
c     ------------------------------------------       
      subroutine interpolate_turbulence(xyz, r3, status)
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
      call null_r3_interpolation(xyz, r3, status)
      end subroutine 
c     ------------------------------------------ 
      subroutine interpolate_turbulence_deriv(xyz, r3, status)
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
      call null_r3_interpolation(xyz, r3, status)  
      end subroutine 
c     ------------------------------------------ 
      subroutine interpolate_currents(xyz, r3, status)
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
      call null_r3_interpolation(xyz, r3, status)  
      end subroutine 
c     ------------------------------------------ 



      subroutine interpolate_zooplankton (xyz, rvec, status) 
c     -------------------------------------------------- 
c     For testing, make lookup in monthly average table
c     Currently just define rvec(1), other components == 0
c     Ignore leap year (may cause minute discontinuity at new year)
c     -------------------------------------------------- 
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: rvec(:)
      integer, intent(out) :: status
      integer :: jday,ilow,iup
      real    :: s,scon,dz
c     -------------------------------------------------- 
      if (is_wet(xyz)) then
         call get_julian_day(master_clock, jday)
         scon = 0.5 + jday*12/365. ! maps mid point of months to 1,2,...12
         ilow = int(scon)
         s    = scon-ilow       ! relative distance between grid points
         iup  = ilow + 1
c        --- impose periodic BC ---     
         if (ilow<1) ilow = 12
         if (iup>12) iup  = 1
c        --- make linear time interpolation in zoop_mon_avg lookup table
         rvec    = 0.! default value
         dz      = zoop_mon_avg(iup)  - zoop_mon_avg(ilow) 
         rvec(1) = zoop_mon_avg(ilow) + dz*s
         status = 0  ! all OK
      else
         status = 1  ! non wet access 
         rvec   = 0. ! default value
         return
      endif
c
      end subroutine  


      subroutine interpolate_temp (xyz, r, status) 
c     ------------------------------------------ 
c     Trigonometric seasonal temperature cycle
c     No day-night difference
c     Ignore leap year (may cause minute discontinuity at new year)
c     ------------------------------------------ 
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
      integer :: jday
      call get_julian_day(master_clock, jday)
      r = avg_wtemp - wtemp_amp*cos(2*pi*(jday-jday_with_mintemp)/365.)
      status = 0    
      end subroutine  


c     ------------------------------------------ 
      subroutine interpolate_wdepth(xy, r, status) 
      real, intent(in)     :: xy(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
      r      = wdepth
      status = 0    
      end subroutine 
c     ------------------------------------------ 
      LOGICAL function is_wet(xyz)
      real, intent(in)     :: xyz(:) 
      if ((xyz(3)>0).and.(xyz(3)<wdepth)) then
         is_wet = .true.
      else
         is_wet = .false.
      endif
      end function
c     ------------------------------------------ 
      LOGICAL function is_land(xy)
      real, intent(in)     :: xy(:)
      is_land = .false.
      end function
c     ------------------------------------------ 
      LOGICAL function horizontal_range_check(xy)
      real, intent(in)     :: xy(:)
      horizontal_range_check = .true.
      end function
c     ------------------------------------------     
      subroutine coast_line_intersection(xyz1, xyz2, anycross, xyzref, 
     +                                   xyzhit) 
      real, intent(in)     :: xyz1(:),xyz2(:)
      logical, intent(out) :: anycross
      real, intent(out)    :: xyzref(:), xyzhit(:)
      anycross = .false. ! all water
      xyzhit   = xyz1    ! just assign a value to avoid compiler warnings
      xyzref   = xyz1    ! just assign a value to avoid compiler warnings
      end subroutine 
c     ------------------------------------------ 
   
      end module
