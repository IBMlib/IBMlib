ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Zooplankton access test module 
c     Constant spatial/temporal environment
c
c     Read these tag from input file:
c         water_depth         =  < water depth in meters >
c         water_temperature   =  < water temperature in deg Celcius >
c         zooplankton_biomass =  < average zooplankton biomass in kg DW/m3 >
c
c     ---------------------------------------------------
c     $Rev: 361 $
c     $LastChangedDate: 2011-07-20 15:09:15 +0200 (Wed, 20 Jul 2011) $
c     $LastChangedBy: asch $ 
c    
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module physical_fields
      use run_context, only: ctrlfile => simulation_file
      use input_parser
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
        

      real :: wdepth   ! water depth in meters
      real :: wtemp    ! water temperature in deg Celcius
      real :: avgzoo   ! average zoo plankton biomass in kg DW/m3
      

      type(clock), target  :: master_clock ! simulation master clock


c     --------
      contains
c     --------


c     ------------------------------------------    
      subroutine init_physical_fields(time)
      type(clock), intent(in),optional :: time 
      if (present(time)) master_clock = time
      write(*,*) trim(get_pbi_version()) 
      call read_control_data(ctrlfile,"water_depth",         wdepth)
      call read_control_data(ctrlfile,"water_temperature",   wtemp)
      call read_control_data(ctrlfile,"zooplankton_biomass", avgzoo)
      write(*,*) "init_physical_fields: water depth =",wdepth,"m"
      write(*,*) "init_physical_fields: water temp  =",wtemp, "degC"
      write(*,*) "init_physical_fields: zooplankton =",avgzoo,"kgDW/m3"

      end subroutine 
c     ------------------------------------------ 
      character*100 function get_pbi_version()  
      get_pbi_version =  "Zooplankton pbi version: $Rev: 361 $"
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
         rvec(1) = avgzoo
         status  = 0
      else
         status = 1
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
c     ---------------------
      if (is_wet(xyz)) then
         r = wtemp
         status  = 0
      else
         status = 1
      endif
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
      if ((xyz(3) >= 0).and.(xyz(3) <= wdepth)) then
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
