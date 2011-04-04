ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Vortex1 pbi 
c     ---------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
c     Stationary gyre at (lon,lat) = (lonvc,latvc) with 
c     tangential current = cspeed m/s
c
c     topography: constant wdepth, 
c                 square (in lat-lon coordinates) water basin (land outside), 
c                 square (in lat-lon coordinates) domain box
c     turbulence: constant vertical/horizontal (no gradients)
c     water temp = constant wtemp
c     all other interpolations = 0
c     ------------------
c     Internal box coordinates: 0 < x,y < 1, boundaries gives by coast lines
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module physical_fields
      use run_context, only: simulation_file
      use time_tools           ! import clock type
      use constants
      use geometry
      use input_parser

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
      public :: interpolate_temp
      public :: interpolate_salty
      public :: interpolate_wind
      public :: interpolate_wdepth
      public :: is_wet    
      public :: is_land
      public :: horizontal_range_check
      public :: coast_line_intersection
      public :: get_pbi_version

c     -------------------- module data --------------------  
c
c     ...... lon-lat wet box parameters (land outside box, wet inside)
c  
      real :: wdepth       ! constant water depth, meters

      real :: coast_east   ! fixed eastern  coast, degrees East
      real :: coast_west   ! fixed western  coast, degrees East
      real :: coast_north  ! fixed northern coast, degrees North
      real :: coast_south  ! fixed southern coast, degrees North 
c
c     ...... lon-lat simulation box (boundary)   
c
      real :: box_east    ! fixed eastern  boundary, degrees East
      real :: box_west    ! fixed western  boundary, degrees East
      real :: box_north   ! fixed northern boundary, degrees North
      real :: box_south   ! fixed southern boundary, degrees North
c
c     ...... hydrodynamic parameters
c
      real :: latvc     ! vortex center, degrees North
      real :: lonvc     ! vortex center, degrees East
      real :: cspeed    ! clockwise current speed, m/s
      real :: wtemp     ! constant water temp, deg Celcius
      real :: vdiff     ! constant vertical   diffusivity, m2/s
      real :: hdiff     ! constant horizontal diffusivity, m2/s
c
c     ...... derived parameters
c
      real  :: lamvc    ! vortex center latitude  in radians
      real  :: phivc    ! vortex center longitude in radians
      real  :: c00(2)   ! SW corner point in degrees
      real  :: c01(2)   ! NW corner point in degrees
      real  :: c10(2)   ! SE corner point in degrees
      real  :: c11(2)   ! NE corner point in degrees

      type(clock), target  :: master_clock ! simulation master clock

      
c     --------
      contains
c     --------

   
      subroutine init_physical_fields(time)
c     ------------------------------------------ 
      type(clock), intent(in),optional :: time 
      if (present(time)) master_clock = time
c     ------------------------------------------ 
      write(*,*) trim(get_pbi_version()) 
c     ------- read input parameters -------
      call read_control_data(simulation_file,"wdepth",wdepth)
      write(*,*) "constant water depth = ",wdepth, "meters"
      call read_control_data(simulation_file,"coast_east",coast_east)
      call read_control_data(simulation_file,"coast_west",coast_west)
      write(*,*) "w/e coast = ",coast_west,coast_east,"degrees East"
      call read_control_data(simulation_file,"coast_north",coast_north)
      call read_control_data(simulation_file,"coast_south",coast_south)
      write(*,*) "s/n coast = ",coast_south,coast_north,"degrees North"
      call read_control_data(simulation_file,"box_east",box_east)
      call read_control_data(simulation_file,"box_west",box_west)
      write(*,*) "w/e boundary = ",box_west,box_east,"degrees East"
      call read_control_data(simulation_file,"box_north",box_north)
      call read_control_data(simulation_file,"box_south",box_south)
      write(*,*) "s/n  boundary",box_south,box_north,"degrees North"
      call read_control_data(simulation_file,"latvc",latvc)
      call read_control_data(simulation_file,"lonvc",lonvc)
      write(*,*) "vortex center = (",lonvc,latvc,") degrees E/N"
      call read_control_data(simulation_file,"cspeed",cspeed)
      write(*,*) "clockwise current speed",cspeed,"m/s"
      call read_control_data(simulation_file,"wtemp",wtemp)
      write(*,*) "constant water temp,",wtemp,"deg Celcius"
      call read_control_data(simulation_file,"vdiff",vdiff)
      write(*,*) "constant vertical   diffusivity",vdiff,"m2/s"
      call read_control_data(simulation_file,"hdiff",hdiff)
      write(*,*) "constant horizontal diffusivity",hdiff,"m2/s"
c     ------- post process input --------
      lamvc    = latvc*deg2rad   ! vortex center in radians
      phivc    = lonvc*deg2rad   ! vortex center in radians
c     --- set coast lines (c00,c01,c10,c11) ---      
      c00 = (/coast_west, coast_south/)   ! SW corner point in degrees
      c01 = (/coast_west, coast_north/)   ! NW corner point in degrees
      c10 = (/coast_east, coast_south/)   ! SE corner point in degrees
      c11 = (/coast_east, coast_north/)   ! NE corner point in degrees
      
      end subroutine 
c     ------------------------------------------ 
      character*100 function get_pbi_version()  
      get_pbi_version =  "Vortex 1 pbi version: $Rev$"
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
      type(clock), intent(in),optional :: time
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
      subroutine interpolate_turbulence(geo, r3, status)
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
      r3(1:2) = hdiff
      r3(3)   = vdiff
      status  = 0
      end subroutine 
c     ------------------------------------------ 
      subroutine interpolate_turbulence_deriv(geo, r3, status)
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
      r3     = 0.
      status = 0      
      end subroutine  
c     ------------------------------------------ 
      subroutine interpolate_currents(geo, r3, status)
c     ..............................................
c     Based on exact steamlines around vortex center
c     ..............................................
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
      real                 :: t(3),phi,lam
      real                 :: phidir(3),lamdir(3)

      if (.not.is_wet(geo)) then
         status = 1
         r3     = 0.
         return
      endif
      phi = geo(1)*deg2rad
      lam = geo(2)*deg2rad

c.....1) evaluate vortex steam line tangent t(3) in (phi,lam)
      t(1) = cos(lam)  * sin(phi)  * sin(lamvc) -
     +       cos(lamvc)* sin(phivc)* sin(lam)
      t(2) = cos(lamvc)* cos(phivc)* sin(lam) -
     +       cos(lam)  * cos(phi)  * sin(lamvc) 
      t(3) = cos(lam)  * cos(lamvc)*
     +       (cos(phi)*sin(phivc) - cos(phivc)*sin(phi))

      t    = t/sqrt(sum(t*t)) ! normalize

c.....2) resolve local tangent space (phidir,lamdir) in global frame
      phidir(1) = -sin(phi)  ! local WE vector in global frame
      phidir(2) =  cos(phi)  ! local WE vector in global frame
      phidir(3) =  0         ! local WE vector in global frame
      lamdir(1) = -sin(lam)*cos(phi) ! local SN vector in global frame
      lamdir(2) = -sin(lam)*sin(phi) ! local SN vector in global frame
      lamdir(3) =  cos(lam)          ! local SN vector in global frame

c.....3) project t onto local tangent space and 
c        scale with clock wise current speed (cspeed)
 
      r3(1)  = cspeed*sum(t*phidir)
      r3(2)  = cspeed*sum(t*lamdir)
      r3(3)  = 0.  ! no vertical currents
      status = 0

      end subroutine 
c     ------------------------------------------ 
      subroutine interpolate_temp (geo, r, status) 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
      r      = wtemp
      status = 0    
      end subroutine 
c     ------------------------------------------ 
      subroutine interpolate_salty(geo, r, status)
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
      r      = 0.
      status = 0    
      end subroutine  
c     ------------------------------------------ 
      subroutine interpolate_wind(geo, r2, status)
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r2(:)
      integer, intent(out) :: status
      r2     = 0.
      status = 0   
      end subroutine 
c     ------------------------------------------ 
      subroutine interpolate_wdepth(geo, r, status) 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
      r      = wdepth
      status = 0    
      end subroutine 
c     ------------------------------------------ 
      LOGICAL function is_wet(geo)
      real, intent(in)     :: geo(:) 
      if (is_land(geo).or.
     +    (geo(3)<0).or.
     +    (geo(3)>wdepth)) then
         is_wet = .false.
      else
         is_wet = .true.
      endif
      
      end function
c     ------------------------------------------ 
      LOGICAL function is_land(geo)
c     ..........................................
c     Impose square lon-lat water basin 
c     ..........................................
      real, intent(in)     :: geo(:)
      if  ((geo(1) < coast_west)  .or.
     +     (geo(1) > coast_east)  .or.
     +     (geo(2) > coast_north) .or.
     +     (geo(2) < coast_south)) then
         is_land = .true.
      else
         is_land = .false.
      endif
      end function
c     ------------------------------------------ 
      LOGICAL function horizontal_range_check(geo)
c     ..........................................
c     Impose horizontal lon-lat box domain  
c     ..........................................
      real, intent(in)     :: geo(:)
      if  ((geo(1) < box_west)  .or.
     +     (geo(1) > box_east)  .or.
     +     (geo(2) > box_north) .or.
     +     (geo(2) < box_south)) then
         horizontal_range_check = .false.
      else
         horizontal_range_check = .true.
      endif
      end function


c     ------------------------------------------ 
      subroutine coast_line_intersection(geo1, geo2, anycross, georef, 
     +                                   geohit) 
c     ..........................................
c     Represent square water box with land around
c     start position geo1 asserted valid 
c     ..........................................
      real, intent(in)     :: geo1(:),geo2(:)
      logical, intent(out) :: anycross
      real, intent(out)    :: georef(:), geohit(:)
      real                 :: direct(size(geo1)), reflect(size(geo1))
      real                 :: assured_wet(size(geo1))
      character*1          :: xface
      real                 :: s
      logical              :: ldum
c     ..........................................
      if (is_land(geo1)) stop "assertion failed: geo1 is on land"
      if (is_land(geo2)) then 
         anycross = .true. ! geo2 dry, we crossed coast line
         call assess_cell_crossing (geo1,geo2,s,xface)
         direct  = geo2-geo1
         reflect = direct
         if     (xface=="N" .or. xface =="S") then     
            reflect(2) = -reflect(2) ! horizontal box line reflection
         elseif (xface=="E" .or. xface =="W") then        
            reflect(1) = -reflect(1) ! vertical   box line reflection
         else
            stop "coast_line_intersection: unhandled reflection type"
         endif
c        |direct| == |reflect|        
         geohit = geo1   +     s*direct ! position where coast line is hit
         georef = geohit + (1-s)*reflect ! position for reflection in coast line
         ldum = point_on_land(geohit,geo1,assured_wet) ! waste logical result
         if (ldum) geohit = assured_wet
      else  
         anycross = .false.    ! georef, geohit not defined   
      endif

c     -----------------------------
      contains  ! local subroutines
c     -----------------------------


      subroutine assess_cell_crossing (g1,g2,s,exitface) ! internal to coast_line_intersection
c     --------------------------------------------------------------
c     Determine where the direct line from g1 through g2 
c     leaves the cell with center (ix,iy)
c     Return 0 < s  which is the coordinate along the vector (g2-g1) 
c     where the vector exits the cell. Notice s may exceed 1 
c     (which is the case, if cell (ix,iy) contains g2)
c     Also identify the exit face as one of "N","S","E","W" 
c     so the proper reflection can be calculated
c
c     Modified from assess_cell_leaving to avoid decision conflicts ASC/18Jan2011
c     --------------------------------------------------------------
      real, intent(in)         :: g1(:),g2(:)
      real, intent(out)        :: s
      character*1, intent(out) :: exitface
c
      logical                  :: yes
      real                     :: stest,tdum,dg(2)
c     --------------------------------------------------------------      
c
c     Determine how the step leaves the box defined by (c00,c01,c10,c11))
c     We can make efficiency gains by taking advantage of the direction
c     that the particle is travelling in - a particle travelling NE
c     can only leave by the North and East faces etc
c
      exitface = "@"    ! we should not pick this one up now, but ...
      s        = 1.0e12 ! must be set before first comparison
      dg     = g2(1:2)-g1(1:2)
c
c     First the east-west direction
c
      if (dg(1) > 0) then    !we're moving to the east

         call cross_2Dlines(g1,g2,c10,c11,stest,tdum,yes)
         if (yes.and.(stest<s)) then ! rely on right-to-left evaluation
            s         = stest
            exitface  = "E"   
         endif
      else   !we're moving to the west
         call cross_2Dlines(g1,g2,c00,c01,stest,tdum,yes)
         if (yes.and.(stest<s)) then ! rely on right-to-left evaluation
            s         = stest
            exitface  = "W"   
         endif
      endif
c
c     Then the north-south direction
c
      if (dg(2) > 0) then !we're moving to the North
         call cross_2Dlines(g1,g2,c01,c11,stest,tdum,yes)
         if (yes.and.(stest<s)) then ! rely on right-to-left evaluation
            s         = stest
            exitface  = "N"   
         endif
      else  !we're moving to the south
         call cross_2Dlines(g1,g2,c00,c10,stest,tdum,yes)
         if (yes.and.(stest<s)) then ! rely on right-to-left evaluation
            s         = stest
            exitface  = "S"   
         endif
      endif

      end subroutine assess_cell_crossing     ! local subroutine



      logical function point_on_land(maybewet, wetpt, assured_wetpt)
c     --------------------------------------------------------------
c     Ensure that maybewet is wet (i.e. is_land(maybewet) == .false.) 
c     by nudging it toward wetpt in very small steps to weed
c     out uncertainty in numerical solution of coast line intersection
c     equations. Perform algortihm in geo coordinates. Do not affect vertical
c     component of maybewet (if provided)
c     
c     Input:  maybewet      : may or may not be dry
c             wetpt         : assumed wet point
c
c     Output: assured_wetpt : a guarantied wet point if maybewet is dry
c                             undefined if maybewet tests wet          
c             point_on_land:  return .true. when maybewet was dry 
c                             and nudged into assured_wetpt    
c                             return .false. when maybewet was wet 
c                             and return assured_wetpt == maybewet 
c         
c     Assumptions: there exists wet points arbitrary close to wetpt
c                  is_land is the authoritative function for horizontal 
c                  wet/dry condition
c
c     ASC/09Feb2011
c     --------------------------------------------------------------
      real, intent(in)    :: maybewet(:), wetpt(:)
      real, intent(out)   :: assured_wetpt(:)      ! assumed same length as maybewet
      real                :: dg(2),scale,dgmax
      integer             :: ticks,istep
      logical             :: keep_going
      real,parameter      :: resol = 1.0e-6   ! approx last sign digit
      real,parameter      :: scale_default = 0.99    ! scale_min < scale_default < scale_max 
      real,parameter      :: scale_max     = 0.999 ! 0 < scale_min < scale_max < 1
      real,parameter      :: scale_min     = 0.5     ! 0 < scale_min < scale_max < 1
      integer,parameter   :: max_steps     = 1000
c     ------------------------------------------------------------      
      point_on_land = is_land(maybewet) ! .true. : needs nudging

      if (.not.point_on_land) return    ! no further processing
c     
c     resolve scaling factor (0 < scale_min < scale < scale_max < 1)
c     
      dg    = maybewet(1:2) - wetpt(1:2)
      dgmax = maxval(abs(dg))
      ticks = nint(dgmax/resol)
      if (ticks>0) then
         scale = abs(1.0*(ticks-1.0)/ticks)
         scale = min(max(scale, scale_min), scale_max)
      else
         scale = scale_default
      endif
c     
c     nudging loop: set assured_wetpt
c            
      istep         = 1
      assured_wetpt = maybewet ! load vertical component, if present
      keep_going    = .true. ! we know we need at least one step
      do while (keep_going)
         dg = dg*scale  ! scale down difference vector iteratively
         assured_wetpt(1:2) = wetpt(1:2) + dg ! 
         keep_going = is_land(assured_wetpt)
         istep = istep + 1
         if (istep > max_steps) then
            write(*,*) "point_on_land: max_steps exceeded"
            write(*,*) "maybewet=", maybewet
            write(*,*) "wetpt   =", wetpt
            write(*,*) "scale   =", scale
            stop ! fatal, no plan B available
         endif
      enddo
      
      end function point_on_land


      end subroutine


    
      end module
