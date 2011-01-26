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
c                 square water basin (land outside), 
c                 square domain box
c     turbulence: constant vertical/horizontal (no gradients)
c     water temp = constant wtemp
c     all other interpolations = 0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
      real, parameter :: wdepth      = 20.0 ! constant water depth, meters

      real, parameter :: coast_east  = 7.0  ! fixed eastern  coast, degrees East
      real, parameter :: coast_west  = 1.0  ! fixed western  coast, degrees East
      real, parameter :: coast_north = 60.0 ! fixed northern coast, degrees North
      real, parameter :: coast_south = 50.0 ! fixed southern coast, degrees North 
c
c     ...... lon-lat simulation box (boundary)   
c
      real, parameter :: box_east  =  8.0  ! fixed eastern  boundary, degrees East
      real, parameter :: box_west  =  0.0  ! fixed western  boundary, degrees East
      real, parameter :: box_north = 61.0  ! fixed northern boundary, degrees North
      real, parameter :: box_south = 49.0  ! fixed southern boundary, degrees North
c
c     ...... hydrodynamic parameters
c
      real, parameter :: latvc    = 55.0 ! vortex center, degrees North
      real, parameter :: lonvc    = 3.0  ! vortex center, degrees East
      real, parameter :: cspeed   = 1.0  ! clockwise current speed, m/s
      real, parameter :: wtemp    = 8.0  ! constant water temp, deg Celcius
      real, parameter :: vdiff    = 5.**2/86400. ! constant vertical   diffusivity, m2/s
      real, parameter :: hdiff    = 0.   ! constant horizontal diffusivity, m2/s
c
c     ...... derived parameters
c
      real, parameter :: lamvc    = latvc*deg2rad 
      real, parameter :: phivc    = lonvc*deg2rad  

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
      subroutine interpolate_turbulence(xyz, r3, status)
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
      r3(1:2) = hdiff
      r3(3)   = vdiff
      status  = 0
      end subroutine 
c     ------------------------------------------ 
      subroutine interpolate_turbulence_deriv(xyz, r3, status)
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
      r3     = 0.
      status = 0      
      end subroutine  
c     ------------------------------------------ 
      subroutine interpolate_currents(xyz, r3, status)
c     ..............................................
c     Based on exact steamlines around vortex center
c     ..............................................
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
      real                 :: t(3),phi,lam
      real                 :: phidir(3),lamdir(3)

      if (.not.is_wet(xyz)) then
         status = 1
         r3     = 0.
         return
      endif
      phi = xyz(1)*deg2rad
      lam = xyz(2)*deg2rad

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
      subroutine interpolate_temp (xyz, r, status) 
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
      r      = wtemp
      status = 0    
      end subroutine 
c     ------------------------------------------ 
      subroutine interpolate_salty(xyz, r, status)
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
      r      = 0.
      status = 0    
      end subroutine  
c     ------------------------------------------ 
      subroutine interpolate_wind(xy, r2, status)
      real, intent(in)     :: xy(:)
      real, intent(out)    :: r2(:)
      integer, intent(out) :: status
      r2     = 0.
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
      if (is_land(xyz).or.
     +    (xyz(3)<0).or.
     +    (xyz(3)>wdepth)) then
         is_wet = .false.
      else
         is_wet = .true.
      endif
      
      end function
c     ------------------------------------------ 
      LOGICAL function is_land(xy)
c     ..........................................
c     Impose square lon-lat water basin 
c     ..........................................
      real, intent(in)     :: xy(:)
      if  ((xy(1) < coast_west)  .or.
     +     (xy(1) > coast_east)  .or.
     +     (xy(2) > coast_north) .or.
     +     (xy(2) < coast_south)) then
         is_land = .true.
      else
         is_land = .false.
      endif
      end function
c     ------------------------------------------ 
      LOGICAL function horizontal_range_check(xy)
c     ..........................................
c     Impose horizontal lon-lat box domain  
c     ..........................................
      real, intent(in)     :: xy(:)
      if  ((xy(1) < box_west)  .or.
     +     (xy(1) > box_east)  .or.
     +     (xy(2) > box_north) .or.
     +     (xy(2) < box_south)) then
         horizontal_range_check = .false.
      else
         horizontal_range_check = .true.
      endif
      end function


c     ------------------------------------------ 
      subroutine coast_line_intersection(xyz1, xyz2, anycross, xyzref, 
     +                                   xyzhit) 
c     ..........................................
c     Represent square water box with land around
c     start position xyz1 asserted valid 
c     ..........................................
      real, intent(in)     :: xyz1(:),xyz2(:)
      logical, intent(out) :: anycross
      real, intent(out)    :: xyzref(:), xyzhit(:)
      real                 :: tx,ty,t,xr,yr
      real                 :: x1,y1,x2,y2
c     ..........................................
      if (is_land(xyz1)) stop "assertion failed: xyz1 is on land"
      if (is_land(xyz2)) then 
         anycross = .true. ! xyz2 dry, we crossed coast line
c        transform xys to box coordinates and compute reflections/crossings 
c        (x,y) reflections are handled independently and parallel in box coordinates         
         call lonlat2boxcoor(xyz1, x1,y1)
         call lonlat2boxcoor(xyz2, x2,y2)
         call umklapp(x1,x2,xr,tx)  ! in box coor
         call umklapp(y1,y2,yr,ty)  ! in box coor
         t         = min(tx,ty)
         xyzhit    = xyz1 + t*(xyz2-xyz1) ! in lonlat
         xyzref    = xyz2 ! copy z, if present
         call boxcoor2lonlat(xr,yr,xyzref) ! to lonlat
      else  
         anycross = .false.    ! xyzref, xyzhit not defined   
      endif
      end subroutine
c     ------------------------------------------ 
      subroutine lonlat2boxcoor(xy,x,y)
      real, intent(in)     :: xy(:)
      real, intent(out)    :: x,y
      x = (xy(1) - coast_west) /(coast_east  - coast_west)
      y = (xy(2) - coast_south)/(coast_north - coast_south)
      end subroutine
c     ------------------------------------------ 
      subroutine boxcoor2lonlat(x,y,xy)
      real, intent(in)    :: x,y
      real, intent(out)   :: xy(:)
      xy(1) = coast_west  + x*(coast_east  - coast_west)
      xy(2) = coast_south + y*(coast_north - coast_south)
      end subroutine boxcoor2lonlat

c     ----------------------------------------------------------
      subroutine umklapp(s1,s2,sr,t)
c     ----------------------------------------------------------
c     Compute (possibly multiply) reflection of s2 when moving
c     from s1 to s2 but confined to the unit interval 0 < s < 1
c
c     Assert 0 < s1 < 1
c     Output: sr = multiple reflection of s2 (if s2 inside  [0;1] sr undefined)
c             t  = relative position along s1->s2 where s=0/1
c                  is crossed. always t>0
c     note: s1 + (s2-s1)*t == 0 .or. 1
c           if s2 inside  [0;1], 1<t   and sr undefined
c           if s2 outside [0;1], 0<t<1 
c
c           mod(-q,1.) = -mod(q,1.)
c     ----------------------------------------------------------
      real, intent(in)  :: s1,s2
      real, intent(out) :: sr,t
c     
      integer           :: ibox
      real              :: ds
      logical           :: even, inside
      real,parameter    :: dslimit = 1.0e-12 ! lower step limit 
c     ----------------------------------------------------------'
      if ((s1>1.).or. (s1<0.)) stop "umklapp: s1 not in [0;1]"
c
c     ============== categorize interval where s2 is ==============
c 
      ibox = floor(s2)
      if (mod(ibox,2) == 0) then
         even = .true.
      else
         even = .false.
      endif

      if (ibox == 0) then
         inside = .true.
      else
         inside = .false.
      endif   
c
c     ============== generate reflection sr  ==============
c
      if (inside) then  
         sr = s2        ! s2 inside principal range -> no reflection
      else
c     ...... s2 is outside principal range
c     ...... generate reflected position sr ......    
c            mod(-q,1.) = -mod(q,1.)
c
         if (s2>0) then  
            if (even) then
               sr = mod(s2, 1.) ! even box order
            else
               sr = 1. - mod(s2, 1.) ! odd box order
            endif
         else                   ! s2<0
            if (even) then
               sr = 1. + mod(s2, 1.) ! even box order
            else
               sr = - mod(s2, 1.) ! odd box order
            endif
         endif
c
      endif

      if ((sr<0).or.(sr>1)) stop "umklapp: sr not in [0;1]"
c
c     ============== generate wall hit point t along ds  ==============
c
c     algorithm works for s2 both inside/outside principal range
c
      ds = s2-s1
c.....sign(a, b) returns the absolute value of a times the sign of b
      if (abs(ds)<dslimit) ds = sign(dslimit, ds)

      if (s2>s1) then
         t = (1.-s1)/ds
      else ! s2<s1
         t = -s1/ds
      endif
 
      end subroutine umklapp

    
      end module
