ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                 Time tool module                      c
c                                                       c
c                                                       c
c     recompiled asbjorn(Aug2005)                       c
c     default scope public                              c
c                                                       c
c     test compilation mode: (libtime must be compiled) c
c         ifort -e90 time_tools.f  -Llibtime -ltime77   c                       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module time_tools
      implicit none

c-----------------------------------------------------------
c     clock/period are interface classes to Thomas Forbriger's 
c     generic libtime (1997) that facilitate hydrodynamical/biological
c     flavor time handling
c
c     Make clock content private to allow changing implementation
c-----------------------------------------------------------
c        date(1)=year
c        date(2)=day of year
c        date(3)=hour
c        date(4)=minute
c        date(5)=second
c        date(6)=millisecond
c        date(7)=microsecond
c-----------------------------------------------------------
      type clock
      private
        integer :: date(7)
      end type

      interface set_clock
        module procedure set_clock_4i
        module procedure set_clock_from_clock
      end interface

      type time_period
      private
        integer :: date_start(7)
        integer :: date_end(7)
      end type

      interface set_period
        module procedure set_period_i8
        module procedure set_period_from_clocks
      end interface

      interface get_period_length_sec
        module procedure get_period_length_sec_period
        module procedure get_period_length_sec_clocks
      end interface

      interface get_period_length_hour
        module procedure get_period_length_hour_cc_i
        ! module procedure get_period_length_hour_cc_r 
        ! module procedure get_period_length_hour_p_i
        ! module procedure get_period_length_hour_p_r
      end interface

      interface get_relative_position
        module procedure get_relative_position_period
        module procedure get_relative_position_clocks
      end interface

      integer,parameter,private :: days_in_mon_noleap(12) = 
     *                     (/31,28,31,30,31,30,31,31,30,31,30,31/)


c-----------------------------------------------------------
c     Basic performance profiling tools
c     
c     Provides a timer class and associated functions, based around
c     the Fortran 90 function system_clock(). Timing precision is
c     system and compiler dependent, but is sufficient for most
c     simple tasks
c     Timing of a routine is achieved through two simple function calls
c        tic(timer)                        : Starts a timer
c        toc(timer,desc string,out_timer)  : Displays time elapsed
c     Inspiration for this style of timer comes from MATLAB
c-----------------------------------------------------------
      type timer
        integer counts
        integer count_rate
      end type
      
      public :: tic
      public :: toc


      contains


      subroutine tic(timer_out)
c     ------------------------------------------------------------
c     Starts a timer
c     The name is easy to remember - remember that a clock goes "tic, toc"
c     ------------------------------------------------------------
      type(timer),intent(out) :: timer_out
c     ------------------------------------------------------------
      call system_clock(timer_out%counts,timer_out%count_rate)
      end subroutine tic

      subroutine toc(timer_in,message,timer_out)
c     ------------------------------------------------------------
c     Calculates the time elapsed since the timer_in was initialised
c     Displays a message with that time elapsed. Optionally returns
c     the time at which the function was call.
c     The name is easy to remember - remember that a clock goes "tic, toc"
c     ------------------------------------------------------------
      type(timer),intent(in) :: timer_in
      type(timer),intent(out),optional :: timer_out
      type(timer) timer_now
      character (len=*), intent(in) :: message
      integer:: elapsed_counts
      real :: elapsed_time
      character*20 outstr
c     ------------------------------------------------------------
      !Calculate elapsed time
      call system_clock(timer_now%counts,timer_now%count_rate)
      elapsed_counts =timer_now%counts-timer_in%counts
      elapsed_time =real(elapsed_counts)/real(timer_now%count_rate)
      !Write output
      write(outstr,'(F12.4)') elapsed_time
      write(*,42) trim(message),trim(adjustl(outstr))
42    format ("toc : ",a," = ",a," s")
      !Return current time if requested
      if(present(timer_out)) timer_out=timer_now
      end subroutine toc
      
      subroutine set_clock_4i(aclock, year, month, day, second_of_day)
c     ------------------------------------------------------------
c     ------------------------------------------------------------
      type(clock)  :: aclock
      integer      :: year, month, day, second_of_day
c     ------------------------------------------------------------
      aclock%date     = 0  ! total reset
      aclock%date(1)  = year
      call time_setdoy(day, month, aclock%date)
      aclock%date(5)  = second_of_day
      call time_norm(aclock%date)
      end subroutine set_clock_4i
      

      subroutine set_clock_from_clock(aclock, src_clock)
c     ------------------------------------------------------------
c     variant to copy clocks
c     ------------------------------------------------------------
      type(clock)  :: aclock, src_clock
      aclock%date = src_clock%date
      end subroutine set_clock_from_clock




      subroutine set_period_i8(aperiod, i8)
c------------------------------------------------------------------------
c     i8 = (year, month, day, sec_of_day)| (year, month, day, sec_of_day)  
c            1      2     3       4      |   5      6     7       8
c                  period start          |         period end    
c------------------------------------------------------------------------
      type(time_period)  :: aperiod
      integer            :: i8(8)
c------------------------------------------------------------------------
c.....set date_start
      aperiod%date_start    = 0  ! total reset
      aperiod%date_start(1) = i8(1)
      call time_setdoy(i8(3), i8(2), aperiod%date_start)
      aperiod%date_start(5) = i8(4)
      call time_norm(aperiod%date_start)
c.....set date_end
      aperiod%date_end    = 0  ! total reset
      aperiod%date_end(1) = i8(5)
      call time_setdoy(i8(7), i8(6), aperiod%date_end)
      aperiod%date_end(5) = i8(8)
      call time_norm(aperiod%date_end)

      end subroutine set_period_i8


      subroutine set_period_from_clocks(aperiod, clock_start, clock_end)
c------------------------------------------------------------------------   
c------------------------------------------------------------------------
      type(time_period)  :: aperiod
      type(clock)        :: clock_start, clock_end
c------------------------------------------------------------------------  
      aperiod%date_start = clock_start%date
      call time_norm(aperiod%date_start)
      aperiod%date_end   = clock_end%date
      call time_norm(aperiod%date_end)

      end subroutine set_period_from_clocks






      subroutine add_seconds_to_clock(aclock, secs)
c     ------------------------------------------------
c     Works for both positive and negative increments secs
c     time_add includes a time normalization
c     ------------------------------------------------
      type(clock) :: aclock
      integer     :: old_date(7), secs
      integer     :: dt(7)
c     ------------------------------------------------
      old_date = aclock%date
      dt    = 0
      dt(5) = secs
      call time_add(old_date, dt, aclock%date)
      end subroutine add_seconds_to_clock


      integer function get_POSIXtime(aclock)
c     ------------------------------------------------
c     Converts a clock to POSIX (UNIX) time by doing 
c     a simple subtraction. The UNIX epoch is defined
c     as midnigth on 1st Jan 1970, so we simply set that
c     as a concstant and then subtract from the supplied
c     time.
c     ------------------------------------------------
      type(clock),intent(in) :: aclock
      type(clock) UNIX_epoch
      integer           :: secs
c     ------------------------------------------------
      call set_clock(UNIX_epoch,1970,1,1,0)
      call get_period_length_sec(UNIX_epoch, aclock, secs)
      get_POSIXtime=secs
      end function


c     TODO: removed julday - update interface

      subroutine get_date_from_clock(aclock, year, month, day)
c     -----------------------------------------------------
c     Extract traditional date components
c     -----------------------------------------------------
      type(clock), intent(in) :: aclock
      integer, intent(out)    :: year, month, day
c     -----------------------------------------------------
      call time_getdate(day, month, aclock%date)
      year   = aclock%date(1)
      end subroutine get_date_from_clock



c     TODO: renamed function get_second_of_day -> get_second_in_day - update interface

      subroutine get_second_in_day(aclock, isec)
c     -----------------------------------------------------
c     Compute accumulated number of seconds this day,
c     corresponding to time aclock, 
c     Return value satisfies 0 <= isec < 86400
c
c     ignores millisecond/microsecond, if set
c     -----------------------------------------------------
      type(clock), intent(in) :: aclock
      integer, intent(out)     :: isec
      isec = 3600*aclock%date(3) + 60*aclock%date(4) + aclock%date(5)
      end subroutine get_second_in_day


      subroutine get_julian_day(aclock, julday)
c     -----------------------------------------------------
      type(clock), intent(in)  :: aclock
      integer, intent(out)     :: julday
c     -----------------------------------------------------
      julday = aclock%date(2)
      end subroutine get_julian_day



      subroutine get_time_components(aclock, 
     +                               year, month_in_year, day_in_month, 
     +                               julian_day, hour_in_day,
     +                               minute_in_hour, second_in_minute)  
c     -----------------------------------------------------
c     Most flexible time access variant with complete optional 
c     output set
c     -----------------------------------------------------
      type(clock), intent(in)        :: aclock
      integer, intent(out), optional :: year,month_in_year,day_in_month 
      integer, intent(out), optional :: julian_day, hour_in_day
      integer, intent(out), optional :: minute_in_hour, second_in_minute

      integer                         :: month,day
c     -----------------------------------------------------
      call time_getdate(day, month, aclock%date)
      if (present(year))             year             = aclock%date(1)
      if (present(month_in_year))    month_in_year    = month
      if (present(day_in_month))     day_in_month     = day
      if (present(julian_day))       julian_day       = aclock%date(2)
      if (present(hour_in_day))      hour_in_day      = aclock%date(3)
      if (present(minute_in_hour))   minute_in_hour   = aclock%date(4)  
      if (present(second_in_minute)) second_in_minute = aclock%date(5)
      
      end subroutine get_time_components


      integer function compare_clocks(clock1, clock2)
c     -----------------------------------------------------
c     res = 0   for   clock1 = clock2
c     res = -1  for   clock1 < clock2
c     res = 1   for   clock1 > clock2
c     -----------------------------------------------------
      type(clock) :: clock1, clock2
      integer, external :: time_compare
c     -----------------------------------------------------
      compare_clocks = time_compare(clock1%date, clock2%date)    
      end function compare_clocks


      logical function time_inside_period(aclock, aperiod)
c     -----------------------------------------------------
      type(clock)       :: aclock
      type(time_period) :: aperiod
      integer, external :: time_compare
      integer           :: is, ie
c     -----------------------------------------------------
      is = time_compare(aclock%date,      aperiod%date_start) 
      ie = time_compare(aperiod%date_end, aclock%date) 
      time_inside_period = (is==1).and.(ie==1)

      end function time_inside_period


      subroutine get_period_length_sec_period(aperiod, numsec)
c     -----------------------------------------------------      
      type(time_period) :: aperiod 
      integer           :: numsec
      integer           :: reldiff(7) 
c     -----------------------------------------------------
      call time_sub(aperiod%date_start, aperiod%date_end, reldiff)
      numsec = 86400*reldiff(2) + 3600*reldiff(3) + 
     +         60*reldiff(4) + reldiff(5)
      end subroutine get_period_length_sec_period


      subroutine get_period_length_sec_clocks(clock_start, clock_end, 
     +                                        numsec)
c     -----------------------------------------------------      
      type(clock)       :: clock_start, clock_end   
      integer           :: numsec
      integer           :: reldiff(7)   
c     -----------------------------------------------------
      call time_sub(clock_start%date, clock_end%date, reldiff)
      numsec = 86400*reldiff(2) + 3600*reldiff(3) + 
     +         60*reldiff(4) + reldiff(5)
      end subroutine get_period_length_sec_clocks


      subroutine get_period_length_hour_cc_i(clock_start, clock_end, 
     +                                        numhours)
c     -----------------------------------------------------
c     For non integral hours, round to nearest hour
c     -----------------------------------------------------      
      type(clock)       :: clock_start, clock_end   
      integer           :: numhours
      integer           :: reldiff(7)   
c     -----------------------------------------------------
      call time_sub(clock_start%date, clock_end%date, reldiff)
      numhours = 24*reldiff(2) + reldiff(3) + 
     +         nint(reldiff(4)/60.0 + reldiff(5)/3600.0)
      end subroutine get_period_length_hour_cc_i



      subroutine get_relative_position_period(aperiod, aclock, rpos)
c     -----------------------------------------------------  
c     -----------------------------------------------------    
      type(time_period) :: aperiod 
      type(clock)       :: aclock
      real              :: rpos
      real              :: intvlen, t_since_start
      integer           :: reldiff(7), imulfac   
      integer, external :: time_compare   
c     -----------------------------------------------------
      call time_sub(aperiod%date_start, aperiod%date_end, reldiff)
      intvlen = 86400.*reldiff(2) + 3600.*reldiff(3) + 
     +         60.*reldiff(4) + reldiff(5)
c.....time_sub only evaluate the unsigned time diff
      call time_sub(aperiod%date_start, aclock%date, reldiff)
      t_since_start = 86400.*reldiff(2) + 3600.*reldiff(3) + 
     +                60.*reldiff(4) + reldiff(5)
      imulfac = time_compare(aclock%date,aperiod%date_start)
      t_since_start = t_since_start*imulfac ! now t_since_start is signed
      rpos = t_since_start/intvlen
      end subroutine get_relative_position_period


      subroutine get_relative_position_clocks(clock_start, clock_end,
     +                                        aclock, rpos)
c     -----------------------------------------------------  
c     Compute the relative position (rpos) of aclock in relation 
c     to a period (clock_start, clock_end)
c     aclock before aperiod: rpos < 0
c     aclock inside aperiod: 0 < rpos < 1
c     aclock after  aperiod: rpos > 1
c     -----------------------------------------------------     
      type(clock)       :: clock_start, clock_end, aclock
      real              :: rpos
      real              :: intvlen, t_since_start
      integer           :: reldiff(7), imulfac   
      integer, external :: time_compare
c     -----------------------------------------------------
      call time_sub(clock_start%date, clock_end%date, reldiff)
      intvlen = 86400.*reldiff(2) + 3600.*reldiff(3) + 
     +         60.*reldiff(4) + reldiff(5)
c.....time_sub only evaluate the unsigned time diff
      call time_sub(clock_start%date, aclock%date, reldiff)
      t_since_start = 86400.*reldiff(2) + 3600.*reldiff(3) + 
     +                60.*reldiff(4) + reldiff(5)
      imulfac = time_compare(aclock%date,clock_start%date)
      t_since_start = t_since_start*imulfac  ! now t_since_start is signed
      rpos = t_since_start/intvlen
      end subroutine get_relative_position_clocks


c result = 0   for   date1 = date2
c result = -1  for   date1 < date2
c result = 1   for   date1 > date2




      subroutine write_clock(aclock)
      type(clock) :: aclock
      write(*,354) aclock%date(1),aclock%date(2),
     +             aclock%date(3),aclock%date(4)
 354  format("year",i5," Julian day ",i4,i3,"h", i3,"m")
      end subroutine write_clock

      character*20 function  get_datetime(aclock)
      type(clock) :: aclock
      integer :: cyear, cmonth, cday

      call get_date_from_clock(aclock, cyear, cmonth, cday)
      write(get_datetime,355) aclock%date(1),cmonth,cday,
     +             aclock%date(3),aclock%date(4),aclock%date(5)
 355  format(i4,".",i2.2,".",i2.2,"__",i2.2,".", i2.2,".",i2.2)
      end function get_datetime


      subroutine write_time_period(aperiod)
      type(time_period) :: aperiod
      write(*,355) aperiod%date_start(1),aperiod%date_start(2),
     +             aperiod%date_start(3),aperiod%date_start(4)
      write(*,356) aperiod%date_end(1),aperiod%date_end(2),
     +             aperiod%date_end(3),aperiod%date_end(4)
 355  format("period start: year",i5," Julian day ",i4,i3,"h", i3,"m")
 356  format("period end  : year",i5," Julian day ",i4,i3,"h", i3,"m")
      end subroutine write_time_period


      function is_light(aclock,pos,zenith)
c     -------------------------------------------------------------------------
c     A wrapper to day_bounds subroutine that takes a date/time in the form of a
c     clock type together with a lat/long and returns a logical 
c     indicating if its light (true is light/day, false is night).
c     This function defaults to civil sunrise and sunset, which is a zenith 
c     angle of 96 degrees. However, other values can be specified if desired
c     All clocks are assuming to be running on GMT
c     -------------------------------------------------------------------------
      type(clock), intent(in) :: aclock
      real, intent(in) :: pos(:)
      logical :: is_light
      real, optional :: zenith
c     ----- locals -----
      integer :: day, month, year, secs, localOffset
      real    :: lat, lon
      real    :: sunrise, sunset 
      logical :: neverrise, neverset
c     -------------------------------------------------------------------------
      !Zenith angle defaults to civil sunrise,sunset, ie 96.0 degrees
      if(.NOT. present(zenith)) zenith= 96.0
      !Get the day bounds for the given day and position
      call get_date_from_clock(aclock,year,month,day)
      call get_second_in_day(aclock,secs)
      lon=pos(1)
      lat=pos(2)
      call day_bounds(day,month,year,0,lat,lon,zenith,
     &            sunrise,sunset,neverrise,neverset)
      if(neverset) then
        !its light all the time
        is_light = .TRUE.
      elseif(neverrise) then
        !its dark all the time
        is_light = .FALSE.
      else if(sunrise*3600<secs .AND. sunset*3600>secs) then  
        !its light
        is_light = .TRUE.
      else  !its dark
        is_light = .FALSE.
      endif

      end function
      
      

      subroutine day_bounds(day, month, year, localOffset,
     +                      latitude, longitude, zenith, 
     +                      sunrise, sunset, neverrise, neverset)
c -------------------------------------------------------------------------
c                        Sunrise/Sunset Algorithm
c  
c Direct source = http://williams.best.vwh.net/sunrise_sunset_algorithm.htm (Ed Williams)
c Transcribed to fortran 90 by Asbjorn Christensen, DIFRES, Oct. 5 2005.
c 
c
c Open questions: Summer time ??
c
c Source:
c        Almanac for Computers, 1990
c        published by Nautical Almanac Office
c        United States Naval Observatory
c        Washington, DC 20392
c
c Inputs:
c        day, month, year:      date of sunrise/sunset
c        localOffset:           local hour offset wrt. UTC
c        latitude, longitude:   location for sunrise/sunset
c        zenith:                Sun's zenith for sunrise/sunset
c          offical      = 90 degrees 50'
c          civil        = 96 degrees
c          nautical     = 102 degrees
c          astronomical = 108 degrees
c        
c        NOTE: longitude is positive for East and negative for West
c
c Outputs:
c        sunrise, sunset     : local time for event (after decimal as TIME FRACTION)
c        neverrise, neverset : true, if sun never rises/sets
c                              if true, returns: sunrise, sunset = -9.0d10
c Tests:
c        1) reproduces test example
c        
c        2) Central North Sea, Oct 5 2005: 
c            echo $day $month $year $localOffset $latitude $longitude $zenith     |sunrise
c            echo  5    10    2005      1          55.0         4.0    90.8333    |sunrise
c
c           http://aa.usno.navy.mil/data/docs/RS_OneDay.html: (assume $zenith = 90.8333)
c            Sunrise                  06:54                                
c            Sunset                   18:10   
c -------------------------------------------------------------------------
      integer, intent(in)  :: day, month, year, localOffset
      real,    intent(in)  :: latitude, longitude, zenith
      real,    intent(out) :: sunrise, sunset 
      logical, intent(out) :: neverrise, neverset
c     --- locals ---
c         index = 1 for sunrise
c         index = 2 for sunset
      real, parameter :: deg2rad = 1.7453292519943295d-2
      real, parameter :: rad2deg = 5.7295779513082323d1
      integer         :: N1, N2, N3, N, Lquadrant(2), RAquadrant(2)
      real            :: lngHour, tim(2), M(2), L(2), RA(2)
      real            :: sinDec(2),cosDec(2),cosH(2),H(2),T(2),UT(2)
c
c     1. first calculate the day of the year
c
      N1 = 275 * month / 9                              ! integer cast (floor)
      N2 = (month + 9) / 12                             ! integer cast (floor)
      N3 = 1 + (year - 4 * floor(year / 4.0d0) + 2) / 3 ! integer cast (floor)
      N  = N1 - (N2 * N3) + day - 30
c
c     2. convert the longitude to hour value and calculate an approximate time
c        
      lngHour = longitude / 15
      tim(1)    = N + (6  - lngHour) / 24.  ! index = 1 for rising time after here
      tim(2)    = N + (18 - lngHour) / 24.  ! index = 2 for setting time after here
c
c     3. calculate the Sun's mean anomaly
c        
      M = (0.9856 * tim) - 3.289
c
c     4. calculate the Sun's true longitude, confined to [0,360) 
c
      L = M + (1.916 * sin(deg2rad * M)) 
     +      + (0.020 * sin(2 * deg2rad * M)) + 282.634
      L = modulo(L, 360.0) 
c
c     5a. calculate the Sun's right ascension, confined to [0,360) 
c        
      RA = rad2deg * atan(0.91764 * tan(L*deg2rad))
      RA = modulo(RA, 360.0)  
c
c     5b. right ascension value needs to be in the same quadrant as L
c
      Lquadrant  = floor(L/90) * 90
      RAquadrant = floor(RA/90) * 90
      RA = RA + (Lquadrant - RAquadrant)
c
c     5c. right ascension value needs to be converted into hours
c
      RA = RA / 15
c
c     6. calculate the Sun's declination
c
      sinDec = 0.39782 * sin(L*deg2rad)
      cosDec = cos(asin(sinDec))    ! no deg/rad conversion
c
c     7a. calculate the Sun's local hour angle
c        
      cosH = (cos(zenith*deg2rad) - (sinDec * sin(latitude*deg2rad))) / 
     +                              (cosDec * cos(latitude*deg2rad))

c     never(1): the sun never rises on this location (on the specified date)
c     never(2): the sun never sets on this location  (on the specified date)
      
      !MPA 20090917: Fails on day when there is a sunrise, but no sunset, or
      !viceversa. Edited so that this is a neverset/neverrise 
      neverrise = ((cosH(1) >  1) .or. (cosH(2) >  1)) 
      neverset  = ((cosH(1) < -1) .or. (cosH(2) < -1)) 
c      if ((cosH(1) >  1) .and. (cosH(2) <  1)) stop  ! unexpected
c      if ((cosH(1) < -1) .and. (cosH(2) > -1)) stop  ! unexpected
c      if ((cosH(2) >  1) .and. (cosH(1) <  1)) stop  ! unexpected
c      if ((cosH(2) < -1) .and. (cosH(1) > -1)) stop  ! unexpected

      if (neverrise .or. neverset) then
         sunrise = -9.0d10
         sunset  = -9.0d10
         return
      endif
c
c     7b. finish calculating H and convert into hours 
c         (at this point, we know -1 < cosH < 1)
      H(1) = 360.0d0 - rad2deg*acos(cosH(1))   ! sunrise
      H(2) = rad2deg*acos(cosH(2))             ! sunset
      H    = H / 15
c
c     8. calculate local mean time of rising/setting
c        
      T = H + RA - (0.06571 * tim) - 6.622
c
c     9. adjust back to UTC adjusted into the range [0,24) 
c        
      UT = T - lngHour
      UT = modulo(UT, 24.0) 
c
c     10. convert UT value to local time zone of latitude/longitude
c         modulo 24h
c
      sunrise = modulo(UT(1) + localOffset, 24.0)
      sunset  = modulo(UT(2) + localOffset, 24.0)
      
      end subroutine day_bounds



      end module




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                        T E S T   C O D E                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c$$$      program test
c$$$      use time_tools
c$$$      implicit none
c$$$
c$$$      type(clock) :: clck
c$$$      integer     :: i1,i2
c$$$      call set_clock(clck, 1973, 5, 17, 456)
c$$$      call write_clock(clck)
c$$$      call get_time_components(clck, julian_day = i1)
c$$$      write(*,*) "jul day == ", i1
c$$$      call get_time_components(clck, minute_in_hour=i1, 
c$$$     +     second_in_minute=i2)
c$$$      write(*,*) "min == ", i1, "sec == ", i2
c$$$      end program
      

c$$$      program test
c$$$      use time_tools
c$$$      implicit none
c$$$
c$$$      type(clock)       :: clck1, clck2, now
c$$$      type(time_period) :: aperiod
c$$$      real              :: rcc,rp
c$$$      call set_clock(clck1, 1973, 5, 17, 33)
c$$$      call set_clock(clck2, 1973, 5, 18, 33)
c$$$      call set_period(aperiod, clck1, clck2)
c$$$c
c$$$      call set_clock(now, 1973, 5, 16, 33)
c$$$      call get_relative_position(clck1, clck2, now, rcc)
c$$$      call get_relative_position(aperiod, now, rp)
c$$$      write(*,*) rcc,rp
c$$$c      
c$$$      call set_clock(now, 1973, 5, 17, 33+3*86400/4)
c$$$      call get_relative_position(clck1, clck2, now, rcc)
c$$$      call get_relative_position(aperiod, now, rp)
c$$$      write(*,*) rcc,rp
c$$$c      
c$$$      call set_clock(now, 1973, 5, 20, 33)
c$$$      call get_relative_position(clck1, clck2, now, rcc)
c$$$      call get_relative_position(aperiod, now, rp)
c$$$      write(*,*) rcc,rp
c$$$      
c$$$      end program


c      integer function leap_day_this_year(iyea)
c********************************************************
c     Return 1 if this is a leap year, else 0
c********************************************************
c      integer :: iyea, iss
c      iss = 0
c      if((mod(iyea,4).eq.0).and.mod(iyea,100).ne.0) iss=1
c      if(mod(iyea,1000).eq.0) iss=1
c      leap_day_this_year = iss
c      end function leap_day_this_year
      

c      integer function days_this_year(iyea)
c********************************************************
c     Return 1 if this is a leap year, else 0
c********************************************************
c      integer :: iyea
c      days_this_year = 365 + leap_day_this_year(iyea)
c      end function days_this_year


c      subroutine julian_day(julday, iyear, imon, iday)
c********************************************************
c     Return Julian day this year, including leap year check
c********************************************************
c      integer :: julday, iyear, imon, iday
c      julday = iday
c      if (imon>1) then
c         julday = julday + sum(days_in_mon_noleap(1:imon-1))
c         if (imon>2) then
c            julday = julday + leap_day_this_year(iyear)
c         endif
c      endif
c      end subroutine julian_day


c      subroutine julian_day_hashtable(iyear, jdays, imode)
c********************************************************
c     Return Julian day of iyear (includes leap year check)
c     of first/mid day in each month, depending on imode. 
c     Return result into jdays
c     imode == 1: return Julian day of month start
c     imode == 2: return Julian day of middle of month 
c                 (as integer, truncated downward)
c********************************************************
c      integer :: iyear, jdays(*), imode
c      integer :: imon,halfway(12)     
c     --------------------------------------------------
c.....jdays is accumulated Julian day at month start
c     jdays == (1, 32, ...)
c      jdays(1) = 1
c      do imon = 2, 12
c         jdays(imon) = jdays(imon-1) + days_in_mon_noleap(imon-1)
c         if (imon==3) jdays(3) = jdays(3) + leap_day_this_year(iyear)   
c      enddo    
c.....vary result, depending on imode
c      if     (imode == 1) then
c         return
c      elseif (imode == 2) then  ! integer truncation of half day, OK for leap year
c         jdays(1:12) =  jdays(1:12) + days_in_mon_noleap/2  
c         return
c      else
c         stop "unknown imode"
c      endif
c
c      end subroutine julian_day_hashtable
c

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccc             quick tester           ccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      program test
c      use time_tools
c      implicit none
c      integer :: j, jdum(12), iyear
c      write(*,*) 2000, leap_day_this_year(2000)
c      write(*,*) 2001, leap_day_this_year(2001)
c      write(*,*) 2004, leap_day_this_year(2004)
c      write(*,*) 1900, leap_day_this_year(1900)
c      call julian_day(j, 1998, 1, 12)
c      write(*,*) 1998, 1, 12, j
c      call julian_day(j, 1998, 2, 12)
c      write(*,*) 1998, 2, 12, j
c      call julian_day(j, 1998, 3, 12)
c      write(*,*) 1998, 3, 12, j
c      call julian_day(j, 2004, 1, 12)
c      write(*,*) 2004, 1, 12, j
c      call julian_day(j, 2004, 2, 12)
c      write(*,*) 2004, 2, 12, j
c      call julian_day(j, 2004, 3, 12)
c      write(*,*) 2004, 3, 12, j
ccc----------------------------------------------
c      do iyear = 2000, 2004
c         call julian_day_hashtable(iyear, jdum, 1)
c         write(*,333) iyear, "b", jdum
c         call julian_day_hashtable(iyear, jdum, 2)
c         write(*,333) iyear, "m", jdum
c         write(*,*)
c      enddo
c 333  format (i4, 2x, a1, 2x, 12(i4))
c      end program

!       program test_sunrise
! c -------------------------------------------------------------------------------------------
! c     test wrapper to subroutine day_bounds
! c     compilation: ifort -o sunrise sunrise.f
! c     usage      : echo $day $month $year $localOffset $latitude $longitude $zenith     |sunrise
! c     example    : echo  25    6    1990     -4         40.9       -74.3   90.8333333   |./sunrise
! c -------------------------------------------------------------------------------------------
!       implicit none
!       integer :: day, month, year, localOffset
!       real    :: latitude, longitude, zenith
!       real    :: sunrise, sunset 
!       logical :: neverrise, neverset
!       read(5,*) day, month, year, localOffset,
!      +         latitude, longitude, zenith
! 
!       write(*,200) day, month, year,localOffset
!  200  format("day = ",i2,"    month = ",i2,"    year = ", i4,
!      +        "    UT offset(h) = ",i2)
!       write(*,205) latitude, longitude
!  205  format("latitude(deg) = ",f6.2,"    longitude(deg) = ",f6.2)
!       write(*,210) zenith
!  210  format("zenith(deg) for sun rise/setting = ",f8.4)
!  
!       call day_bounds(day, month, year, localOffset,
!      +                latitude, longitude, zenith, 
!      +                sunrise, sunset, neverrise, neverset)
! 
!       if (neverrise) write(*,*) "sun never rises at this day and place"
!       if (neverset)  write(*,*) "sun never sets at this day and place"
!       if (.not.(neverset.or.neverrise)) write(*,250) sunrise, sunset
!  250  format("sunrise(h) = ",f6.2,"    sunset(h) = ",f6.2)
! 
!       end program

      
!      program test_POSIX_time
!c -------------------------------------------------------------------------------------------
!c     test wrapper to subroutine calculating POSIX time 
!c     compilation: ifort -e90 -o POSIX time_tools.f -Llibtime -ltime77
!c     usage      : echo $year $month $day $hour $min $sec   |POSIX
!c     example    : echo 2009  2      13    23   31   30     |./POSIX
!c     answer     : 1234567890
!c -------------------------------------------------------------------------------------------
!      use time_tools      
!      implicit none
!      integer :: year, month, day, hour, min, sec 
!      integer :: secs
!      type(clock) :: clockin
!      
!      read(5,*) year, month, day, hour, min, sec
!      secs = hour*3600 + min*60 + sec
!      call set_clock(clockin,year,month, day, secs)
!      write(*,*) "Input clock : ", get_datetime(clockin)
!      write(*,*) "POSIX time  : ", get_POSIXtime(clockin)      
!      end program

