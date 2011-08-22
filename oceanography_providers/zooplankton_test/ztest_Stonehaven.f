ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Zooplankton access test module
c
c     Implements seasonal variability in the Stonehaven time series
c     partly following Marc Hufnagls parameterization
c     ---------------------------------------------------
c     $Rev: 361 $
c     $LastChangedDate: 2011-07-20 15:09:15 +0200 (Wed, 20 Jul 2011) $
c     $LastChangedBy: asch $ 
c
c     depth          = constant wdepth 
c     no coast line (all water, no land)
c     water temp     = seasonal trigonometric (no stratification)
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

      real, parameter :: avg_wtemp         = 10.0  ! yearly average - deg Celcius
      real, parameter :: wtemp_amp         = 3.9   ! oscillation amp - deg Celcius
      real, parameter :: jday_with_maxtemp = 215.0 ! jday with maximum 
      
c     --- parameters for seasonal cycle of total zooplanton abundance 
c
c      real, parameter :: avg_logN          = 3.23  ! yearly average - log10(N[1/m3])
c      real, parameter :: logN_amp          = 0.57  ! oscillation amp - log10(N[1/m3])
c      real, parameter :: jday_with_maxlogN = 195.0 ! jday with maximum 
c
c      phase  = day2rad*(jday - jday_with_maxlogN)
c      log10N = avg_logN + logN_amp*cos(phase)   ! log base == 10
c      N      = 10**log10N   ! unit = 1/m3
c
      


c     --- parameters for seasonal cycle of slope b of zooplanton spectrum
c         N[1/m3/mm] = a*lprey**b  (b<0)
      real, parameter :: SS_slope_avg = -1.85 
      real, parameter :: SS_slope_amp =  0.7
      real, parameter :: SS_dayoffset = -175+365/2. ! cosine offset     
      real, parameter :: lp_SSmin = 0.02 ! mm, lower Microplankton size limit
      real, parameter :: lp_SSmax = 2.0  ! mm, upper Mesoplankton size limit    

c     --- zooplanton length-weight key (assumed seasonally constant)
c         m = prey_mass0 *(lprey/prey_len0 )**prey_expo

      real, parameter :: prey_len0     = 1.0    ! [mm]
      real, parameter :: prey_expo     = 2.772 
      real, parameter :: prey_mass0    = 21.6197/10**9 ! NB: [kg] at prey_len0     


      real, parameter :: wdepth   = 40.0 ! water depth in meters
       
      real, parameter :: day2rad  = 2*pi/365.0
      
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
      get_pbi_version =  "Stonehaven zooplankton PBI vers: $Rev: 361 $"
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
c
c     zooplankton unit = kg dry weight / m3
c     -------------------------------------------------- 
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: rvec(:)
      integer, intent(out) :: status
      integer :: jday
      real    :: log10N,N,a,slope,phase,xp,s1,s0,xppr,zbio
c     -------------------------------------------------- 
      if (is_wet(xyz)) then
         call get_julian_day(master_clock, jday)
c        --- lookup seasonal slope in zooplankton spectrum (slope)---
         phase  = day2rad*(jday - SS_dayoffset)
         slope  = SS_slope_avg + SS_slope_amp*cos(phase) ! slope < 0
c        --- lookup seasonal zooplankton abundance (N [1/m3]) ---
         N = 0.005 + 160*( (1./(100.*sqrt(2*pi))**2)*   ! from Marc H paper 
     +            exp(-(jday-190.0)**2/2./1e4) )        ! [indv/ml] 
         N = N*1.e6                                      ! [indv/m3] 

c        --- convert N to biomass zbio, applying slope
         xp     = 1.0+slope
         xppr   = 1.0+slope+prey_expo
         s1     = lp_SSmax/prey_len0 
         s0     = lp_SSmin/prey_len0 
         a      = N*xp / prey_len0 / (s1**xp - s0**xp) ! unit = 1/m3/mm
c        
         zbio    = a*prey_mass0*prey_len0*(s1**xppr - s0**xppr)/xppr ! unit = kg/m3
         rvec(1) = zbio
         status  = 0  ! all OK

      else

         status  = 1  ! non wet access 
         rvec    = 0. ! default value
         
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
      r = avg_wtemp + wtemp_amp*cos(day2rad*(jday-jday_with_maxtemp))
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
      real,parameter       :: tol = 1.e-5
      if ((xyz(3)> -tol).and.(xyz(3)<wdepth+tol)) then
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
