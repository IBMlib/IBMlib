ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     API for COARDS compliant 2+1D(time) netcdf data overlays
c
c     * internalize data grids which may differ from main data grid to allow multiple data grids
c     * target 2D netcdf variables stored as (x,y,t)
c     * maskerade if a float variable is stored as float or scaled integer
c     * homogenize grids on load, so that axes are in order x,y and coordinates ascenting     
c     * internalize netcdf handlers
c     * maintain IBMlib time-space interpolation as time sync first, then space interpolation      
c     * data are not necessary wet data, but e.g. atmospheric data, with grid specific validity mask
c     * adapt COARDS layout variations to IBMlib layout for easy access     
c     * follow https://ferret.pmel.noaa.gov/Ferret/documentation/coards-netcdf-conventions
c      
c     https://www.unidata.ucar.edu/software/netcdf/docs-fortran/f90-use-of-the-netcdf-library.html
c
c     In addition to modules imported, this module also applies string_tools.f and
c     libtime (via time_tools)
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module coards_netcdf_2D

      use netcdf       ! decorate and reexport (scope = public)
      use time_tools   ! type clock + clock methods
      use array_tools  ! search_sorted_list
      
      implicit none     
      private      ! default scope
      
      type simple_lonlat_grid
c     --------------------------------------------------------
c     Represent a regular lon-lat grid
c      
c     The grid can self initialize from a netcdf data variable
c     In principle swapxy belongs to the variable, because
c     axes ordering lies in variable definition
c     Support runtime remapping of invalid grid points
c     Make type attributes public for simplicity
c     --------------------------------------------------------
      integer         :: nx,ny
      logical,pointer :: mask(:,:)  ! valid data points (on homogenized grid, assumed static)
      integer,pointer :: remap(:,:) ! (4, nremaps) ; each column =  (ix, iy, use_this_ix, use_this_iy)
      real            :: lon1, lat1 ! SW grid corner point (ix,iy = 1)
      real            :: dlon, dlat ! dlon, dlat > 0
      logical         :: flipx      ! longitude axis descending in data set on file
      logical         :: flipy      ! latitude  axis descending in data set on file
      logical         :: swapxy     ! axis order is (latitude,longitude) of variable at load 
      end type

      
      type coards_xyt_variable
c     ------------------------------------------------------------------------
c     Represent gridded lon-lat real data
c      
c     The variable can load from an netcdf data set by attaching to a variable
c     Cover both static 2D data and time variable 2D data
c     Make type attributes public for simplicity
c     ------------------------------------------------------------------------      
      type(simple_lonlat_grid) :: grid          ! spatial grid of data
      
      integer                  :: ncid          ! attached netcdf set ID 
      integer                  :: varid         ! attached netcdf variable ID within set ncid
      character(len=999)       :: varname       ! store reference for easier tracking
      integer                  :: storage_type  ! int_type / float_type in netcdf data set
      real                     :: scale, offset ! for integer storage, x_real = scale*x_int + offset
      
      integer                  :: nt            ! number of time frames (len(time)). nt=0  for static data
      real, pointer            :: time(:)       ! time of frames (hours since offset) available from set ncid
      type(clock)              :: time_offset   !  offset for time
      
      type(clock)              :: now           ! time corresponding to data
      real, pointer            :: data(:,:)     ! nx,ny - data interpolated to now
      
      integer                  :: frame0        ! current left frame number in time
      real, pointer            :: databuf0(:,:) ! nx,ny - data corresponding to frame0
      integer                  :: frame1        ! current right frame number in time
      real, pointer            :: databuf1(:,:) ! nx,ny - data corresponding to frame1
      
      end type

      integer, parameter :: int_type   = 1
      integer, parameter :: float_type = 2
      

      interface flip_xaxis
        module procedure flip_xaxis_real
        module procedure flip_xaxis_logical
      end interface

      interface flip_yaxis
        module procedure flip_yaxis_real
        module procedure flip_yaxis_logical
      end interface

c     ------- define module public scope -------

      public ::  simple_lonlat_grid
      public ::  coards_xyt_variable
      
      public ::  load_grid                       ! args == (grid, varname, ncid, cast_mode, is_temporal)
      public ::  attach_variable                 ! args == (xyt, varname, ncid, cast_mode)
      public ::  sync_to_time                    ! args == (xyt, when)
      public ::  load_frame                      ! args == (xyt, iframe)
      public ::  interpolate_xyt_variable        ! args == (geo,xyt,deriv,result,status)
      public ::  write_valid_points              ! args == (grid,iunit)
      
      public ::  close_simple_lonlat_grid        ! args == (grid)
      public ::  close_coards_xyt_variable       ! args == (xyt)
         
c     =================================================================
                          contains
c     =================================================================

      subroutine resolve_time_parameters(ncid, time_vname, year, month,
     +                                   day, hour, minute, sec, tunit)
c     -----------------------------------------------------------------
c     Resolve time offset as and time unit in data set from time:unit attribute
c     
c     1) Resolve time offset as (year,month,day,hour,minute,sec)
c     2) Resolve time unit tunit as number of seconds per time unit, i.e.
c        for time unit = sec   tunit = 1
c        for time unit = hours tunit = 3600      
c        for time unit = days  tunit = 86400
c     Cloned from coards_netcdf.resolve_time_parameters 
c     -----------------------------------------------------------------
      integer, intent(in)         :: ncid
      character*(*), intent(in)   :: time_vname
      integer, intent(out)        :: year,month,day,hour,minute,sec
      integer, intent(out)        :: tunit ! number of seconds per time unit
      integer                     :: varid,ihyph,it
      character*999               :: toffstr, text, date, time, msg
      integer                     :: start(99), nwords
c     -----------------------------------------------------------------
c
c     resolve the time offset for the time variable in the data set
c     by parsing the attribute "time:units"
c
      call NetCDFcheck( nf90_inq_varid(ncid, time_vname, varid) )       ! fixed name
      call NetCDFcheck( nf90_get_att(ncid, varid, "units", toffstr) )   ! retrieve time:units attribute
      
c     --------- known variants ---------
c     1) "seconds since 1970-01-01T00:00:00Z"
c     2) "seconds since 1970-01-01 00:00:00"
c     3) "hours since 2014-04-01 01:00:00"
c     4) "seconds since 2014-4-1 00:00:00"
c     case 2-4 with no HH:MM:SS at end (implicitly 00:00:00)

      call tokenize(toffstr, start, nwords) ! in string_tools.f
      if (nwords == 3) then     ! assume variant 1) above
         it = index(toffstr, "T")
         if (it<1) then
            toffstr = trim(adjustl(toffstr)) // " 00:00:00" ! no T;  add default midnight
         else
            toffstr(it:it) = " "      ! remove "T" for easier parsing below
         endif
         call tokenize(toffstr, start, nwords) ! reparse, robust toward spaces
      endif

      
      if (nwords /= 4) then
         write(*,*) "resolve_time_parameters: unable to parse unit attr"
         write(*,*) "unit attr = >",trim(adjustl(toffstr)),"<"
      endif
c     --- capture unit name ---
      text = ""
      read(toffstr(start(1):start(2)-1), *) text
      if (trim(adjustl(text))     == "seconds") then
         tunit = 1
      elseif (trim(adjustl(text)) == "minutes") then
         tunit = 60
      elseif (trim(adjustl(text)) == "hours") then
         tunit = 3600
      elseif (trim(adjustl(text)) == "days") then
         tunit = 86400
      else
         write(*,*) "time unit = >",trim(adjustl(text)),"<"
         msg = "resolve_time_parameters: unable to parse time unit"
         call NetCDFstop(msg)
      endif
c     --- check "since" is next word ---
      text = ""
      read(toffstr(start(2):start(3)-1), *) text
      if (trim(adjustl(text)) /= "since") then
         write(*,*) "unit attr = >",trim(adjustl(toffstr)),"<"
         msg = "resolve_time_parameters: parse error of unit attr"
         call NetCDFstop(msg)
      endif
c     --- capture offset date ---
      date = toffstr(start(3):start(4)-1)   ! like 2014-04-01 or 2014-4-1 
      read(date, 708, err=832) year         !
      ihyph = index(date, "-", back=.true.) ! last "-" in date token
      read(date(6:ihyph-1), *, err=832) month 
      read(date(ihyph+1:),  *, err=832) day
c     --- capture offset time-in-day ---      
      time = toffstr(start(4):)
      read(time, 711, err=833) hour,minute,sec
 708  format(i4)             ! 4 digit year expected
 711  format(i2,1x,i2,1x,i2) ! match known time variants
      
      return
      
 832  write(*,*) "unit attr = >",trim(adjustl(toffstr)),"<"
      msg = "resolve_time_parameters: error parsing offset date"
      call NetCDFstop(msg)
         
 833  write(*,*) "unit attr = >",trim(adjustl(toffstr)),"<"
      msg = "resolve_time_parameters: error parsing offset time"
      call NetCDFstop(msg)
      
      end subroutine resolve_time_parameters

      

      subroutine load_grid(grid, varname, ncid, cast_mode, is_temporal)
c     -----------------------------------------------------------------
c     Initialize spatial grid corresponding to variable varname in set ncid
c     
c     Variable varname is an xy (2D static) or xyt (2D, time variable) type
c     One axis must have name "lon*", the other "lat*" (content after letter 3 not checked)
c     Facilitate normalization of data, so data is in (x,y) order, with ascending coordinates
c     Dimension name and variable name must be coincident. Lon/lat grids are assumed regular, with at least two points
c     Do not assess time aspect, if present. Time dimension must be last, if present
c      
c     cast_mode options:
c         nearest_valid: nearest wet points is used, if a grid point has invalid data
c         native : only valid points are used for interpolations
c
c     Diagnostic attributes to normalize data:      
c       flipx:    longitude axis descending in data set varname in ncid (stored )
c       flipy:    latitude  axis descending in data set varname in ncid
c       swapxy:   data order is (latitude,longitude) in data set varname in ncid
c       is_temporal: true, if variable has a time dimension      
c     -----------------------------------------------------------------
      type(simple_lonlat_grid),intent(out) :: grid
      character*(*),intent(in)             :: varname
      integer,intent(in)                   :: ncid
      character*(*),intent(in)             :: cast_mode
      logical,intent(out)                  :: is_temporal
c      
      integer              :: varid, ndims, dimids(99),xtype
      character*256        :: name1,name2
      integer              :: len1,len2, lonid, latid,istat
      integer              :: dry_int, nremap,ix,iy
      real                 :: dry_real
      real, allocatable    :: lon(:), lat(:), rbuf(:,:)
      logical, allocatable :: native_mask(:,:)
      integer, allocatable :: ibuf(:,:)
c     -----------------------------------------------------------------
      call NetCDFcheck(nf90_inq_varid(ncid, varname, varid)) 
      call NetCDFcheck(nf90_inquire_variable(ncid, varid,ndims=ndims,
     +     dimids=dimids))
c
c     ----- flag whether data is xy or xyt type
c
      if     (ndims == 3) then
         is_temporal = .true.     ! (2D, time variable)
      elseif (ndims == 2) then
         is_temporal = .false.    ! (2D static)
      else
         call NetCDFstop("variable rank /= 2,3")
      endif
c
c     ----- extract basic grid descriptors and orientations
c
c     (x,y) are first two axes - time is third (if present)      
      call NetCDFcheck( nf90_inquire_dimension(ncid,dimids(1),
     +                  name1,len1) )
      call NetCDFcheck( nf90_inquire_dimension(ncid,dimids(2),
     +     name2,len2) )
c     ------ variable is stored with axis order (lon,lat) ------      
      if     ((name1(1:3)=="lon").and.(name2(1:3)=="lat")) then
         grid%swapxy = .false.
         allocate( lon(len1) )
         call NetCDFcheck( nf90_inq_varid(ncid, name1, lonid) )
         call NetCDFcheck( nf90_get_var(ncid, lonid, lon) )
         allocate( lat(len2) )
         call NetCDFcheck( nf90_inq_varid(ncid, name2, latid) )
         call NetCDFcheck( nf90_get_var(ncid, latid, lat) )      
         grid%nx    = len1          
         grid%ny    = len2 
c     ------ variable is stored with opposite axis order (lat,lon) ------
      elseif ((name1(1:3)=="lat").and.(name2(1:3)=="lon")) then
         grid%swapxy = .true.   ! buffer need to be transposed after loading
         allocate( lon(len2) )
         call NetCDFcheck( nf90_inq_varid(ncid, name2, lonid) )
         call NetCDFcheck( nf90_get_var(ncid, lonid, lon) )
         allocate( lat(len1) )
         call NetCDFcheck( nf90_inq_varid(ncid, name1, latid) )
         call NetCDFcheck( nf90_get_var(ncid, latid, lat) )
         grid%nx    = len2
         grid%ny    = len1
      else
         call NetCDFstop("xy axes not recognized")
      endif
c      
c     ---- raw arrays (lon,lat) are now defined
c
c     ---- set longitude axis in grid   
         if (lon(2)>lon(1)) then ! ascending grid
            grid%lon1  = lon(1)
            grid%dlon  = lon(2)-lon(1)
            grid%flipx = .false.
         else                    ! descending grid   
            grid%lon1  = lon(size(lon)) ! last point
            grid%dlon  = lon(1)-lon(2)
            grid%flipx = .true.
         endif
c        ---- set latitude axis  in grid            
         if (lat(2)>lat(1)) then ! ascending grid
            grid%lat1  = lat(1)
            grid%dlat  = lat(2)-lat(1)
            grid%flipy = .false.
         else                    ! descending grid   
            grid%lat1  = lat(size(lat)) ! last point
            grid%dlat  = lat(1)-lat(2)
            grid%flipy = .true.
         endif
         
      deallocate(lon)
      deallocate(lat)
c
c     ----- assess the grid mask by probing for fill values / missing values
c           apply first frame for time dependent sets
c
      allocate( native_mask(len1,len2) )
      call NetCDFcheck( nf90_inquire_variable(ncid, varid,
     +     xtype=xtype) )
      
      if     ((xtype==NF90_SHORT).or.(xtype==NF90_INT)) then
c         ======== integer type storage ========
         istat = nf90_get_att(ncid, varid,
     +                        "_FillValue", dry_int) ! precedence
         if (istat /= nf90_noerr) then
            istat = nf90_get_att(ncid, varid,"missing_value",
     +                           dry_int) ! secondary choice
            if (istat /= nf90_noerr) then ! assume all points valid
               native_mask = .true.
               goto 561
            endif
         endif                  ! istat test 1

         allocate( ibuf(len1,len2) ) ! native shape
         if (is_temporal) then
            call NetCDFcheck( nf90_get_var(ncid, varid, ibuf,
     +           count=(/len1,len2,1/))) ! use first frame
         else
            call NetCDFcheck( nf90_get_var(ncid, varid, ibuf) )
         endif

         where (ibuf==dry_int)
            native_mask = .false.
         elsewhere
            native_mask = .true.
         end where
         
      elseif ((xtype==NF90_FLOAT).or.(xtype==NF90_DOUBLE)) then
c         ======== real type storage ========
         istat = nf90_get_att(ncid, varid,
     +                        "_FillValue", dry_real) ! precedence
         if (istat /= nf90_noerr) then
            istat = nf90_get_att(ncid, varid, "missing_value",
     +                           dry_real) ! secondary choice
            if (istat /= nf90_noerr) then ! assume all points valid
               native_mask = .true.
               goto 561
            endif
         endif                  ! istat test 1

         allocate( rbuf(len1,len2) ) ! native shape
         if (is_temporal) then
            call NetCDFcheck( nf90_get_var(ncid, varid, rbuf,
     +           count=(/len1,len2,1/))) ! use first frame
         else
            call NetCDFcheck( nf90_get_var(ncid, varid, rbuf) )
         endif

         where (rbuf==dry_real)
            native_mask = .false.
         elsewhere
            native_mask = .true.
         end where
         
      else   
         call NetCDFstop("load_grid: unhandled storage type")
      endif                     ! xtype == ...
      
 561  continue  ! 
c
c     normalize native_mask -> grid%mask
c      
      allocate( grid%mask(grid%nx, grid%ny) )
      if (grid%flipx) call flip_xaxis(native_mask)   ! inplace
      if (grid%flipy) call flip_yaxis(native_mask)   ! inplace
      if (grid%swapxy) then
         grid%mask = transpose(native_mask)
      else
         grid%mask = native_mask
      endif
c
c     post process normalized mask, generate/nullify remap
c
      if (cast_mode == "nearest_valid") then
         nremap = count(.not.grid%mask)  ! number of points to remap
         if (nremap == (grid%nx*grid%nx)) then
            call NetCDFstop("grid contains only invalid points") ! we can't fix this
         elseif (nremap == 0) then ! all points are valid, signal no remapping
            nullify (grid%remap)   ! implies grid%mask == .true. everywhere
         else                      ! in this case some point need remapping
            allocate( grid%remap(4,nremap) )
            call set_remap_nearest_valid(grid%mask, grid%remap)
            grid%mask = .true.     ! overwrite loaded mask
         endif
      elseif (cast_mode == "native") then
         nullify (grid%remap)   ! signal no remapping, keep loaded mask
      else
         call NetCDFstop("unknown cast_mode")
      endif
      deallocate( native_mask )
         
      end subroutine load_grid
      


      subroutine write_valid_points(grid,iunit)
c     ----------------------------------------------------
c     Simple diagnostic dump of valid grid points to iunit      
c     ----------------------------------------------------      
      type(simple_lonlat_grid),intent(in) :: grid
      integer,intent(in)                  :: iunit
      integer :: ix,iy
      real    :: x,y
c     ----------------------------------------------------            
      do iy=1,grid%ny
         y = grid%lat1 + (iy-1)*grid%dlat
         do ix=1,grid%nx
            x = grid%lon1 + (ix-1)*grid%dlon
            if (grid%mask(ix,iy)) write(iunit,*) x,y
         enddo
      enddo
      end subroutine write_valid_points
      

      
      subroutine attach_variable(xyt, varname, ncid, cast_mode)
c     -----------------------------------------------------------------
c     Initialize container xyt by attaching to the variable 
c     varname in necdf set ncid
c     Load time vector, but not data. Name "time" is mandatory
c     for time vector
c     -----------------------------------------------------------------
      type(coards_xyt_variable), intent(out)  :: xyt
      character*(*), intent(in)               :: varname
      integer, intent(in)                     :: ncid
      character*(*),intent(in)                :: cast_mode
c      
      integer                      :: xtype, istat, timeid,i
      integer                      :: dimids(99) ! allow graceful handing, if dim mismatch
      integer                      :: sec_of_day, tunit
      logical                      :: is_temporal
      integer                      :: year,month,day,hour,minute,sec
      character*256                :: msg,cdummy
c     -----------------------------------------------------------------
      xyt%ncid    = ncid
      xyt%varname = varname
      call load_grid(xyt%grid, varname, ncid, cast_mode, is_temporal)
      call NetCDFcheck( nf90_inq_varid(ncid, varname, xyt%varid))
      call NetCDFcheck( nf90_inquire_variable(ncid, xyt%varid, cdummy,
     +     xtype=xtype, dimids=dimids) )
c
c     -----  load time, if present   -----
c
      if (is_temporal) then
         call NetCDFcheck(nf90_inquire_dimension(ncid, dimids(3),
     +        len=xyt%nt))
         allocate( xyt%time(xyt%nt) )
         call NetCDFcheck( nf90_inq_varid(ncid, "time", timeid) )
         call NetCDFcheck( nf90_get_var(ncid, timeid, xyt%time) )         
         call resolve_time_parameters(ncid, "time", year, month,
     +        day, hour, minute, sec, tunit)
c
c        NB: paranthesis must be around (tunit/3600.) to avoid accuracy loss due to
c            left-to-right evaluation, when offset is in the far past
c 
         xyt%time = xyt%time*(tunit/3600.) ! convert time to hours since offset
         sec_of_day = hour*3600 + minute*60 + sec
         call set_clock(xyt%time_offset, year, month, day, sec_of_day)
      else
         xyt%nt = 0            ! static data set
         nullify( xyt%time )
      endif
c
c     ------- allocate remaining buffers + init handlers -------
c
c     do not set clock + buffers
      allocate( xyt%data(xyt%grid%nx, xyt%grid%ny) )
      xyt%frame0 = -1       ! signal unloaded
      allocate( xyt%databuf0(xyt%grid%nx, xyt%grid%ny) )
      xyt%frame1 = -1       ! signal unloaded
      allocate( xyt%databuf1(xyt%grid%nx, xyt%grid%ny) )     
c
c     ------- assess whether variable is stored as int/float -------
c      
      if     ((xtype==NF90_SHORT).or.(xtype==NF90_INT)) then
         xyt%storage_type = int_type
c 
         istat = nf90_get_att(ncid, xyt%varid,"scale_factor",
     +        xyt%scale)       ! does possibly not exist
         if (istat /= nf90_noerr) then
            msg = "attach_variable: scale_factor mandatory for"//
     +            "integer-compressed floats"
            call NetCDFstop(msg)
         endif
c 
         istat = nf90_get_att(ncid, xyt%varid,"add_offset",
     +        xyt%offset)       ! does possibly not exist
         if (istat /= nf90_noerr) then
             msg = "attach_variable: add_offset mandatory for"//
     +             "integer-compressed floats"
             call NetCDFstop(msg)
         endif          
c         
      elseif ((xtype==NF90_FLOAT).or.(xtype==NF90_DOUBLE)) then
         xyt%storage_type = float_type
      else
         write(*,*) "attach_variable: for name=", trim(adjustl(varname))
         write(msg,'(a,i)') "attach_variable: unhandled storage type",
     +                      xtype
         call NetCDFstop(msg)
      endif
      
      end subroutine attach_variable



      subroutine load_frame(xyt, iframe)
c     -----------------------------------------------------------------
c     Set data buffer in xyt to frame number iframe
c     -----------------------------------------------------------------      
      type(coards_xyt_variable), intent(inout) :: xyt
      integer,intent(in)                       :: iframe
      integer                                  :: hours, secs
      
c     -----------------------------------------------------------------
      xyt%frame0 = -1   ! signal unloaded
      xyt%frame1 = -1   ! signal unloaded
      call load_xy_frame_to_extbuf(xyt, xyt%data, iframe)
      call set_clock(xyt%now, xyt%time_offset)
      hours = nint( xyt%time(iframe) )
      secs  = 3600*nint( xyt%time(iframe) - hours )
      call add_hours_to_clock(xyt%now, hours)
      call add_seconds_to_clock(xyt%now, secs)
      end subroutine load_frame

      
      
      subroutine sync_to_time(xyt, when)
c     -----------------------------------------------------------------
c     Create linear interpolation in time corresponding to when
c     between neigboring frames      
c     -----------------------------------------------------------------      
      type(coards_xyt_variable), intent(inout) :: xyt
      type(clock),intent(in)                   :: when
c      
      integer          :: it0, it1, istat,i
      real             :: st
c     -----------------------------------------------------------------
      call get_time_coordinates(when, xyt%time, xyt%time_offset,
     +     it0, it1, st, istat)
      if (xyt%frame0 /= it0) then
         call load_xy_frame_to_extbuf(xyt, xyt%databuf0, it0)
         xyt%frame0 = it0
      endif
      if (xyt%frame1 /= it1) then
         call load_xy_frame_to_extbuf(xyt, xyt%databuf1, it1)
         xyt%frame1 = it1
      endif
      xyt%data = xyt%databuf0*(1-st) + xyt%databuf1*st  ! linear time interpolation 
      call set_clock(xyt%now, when)                     ! set time tag
      end subroutine sync_to_time


  
      subroutine get_time_coordinates(when, time, time_offset,
     +     it0, it1, s, istat)
c     -----------------------------------------------------------------
c     Compute interpolation coordinates of when in time array time
c      
c     time is hours since time_offset, corresponding to time frames
c     time does not need to be regular (equally spaced), but must not be degenerate (any dt=0) 
c     Interpolation coordinates is (it0, it1, st) so that
c
c         interpolation = data(it0)*(1-s) + data(it1)*s
c
c     It is guarantied that 1 <= it0 <= it1 <= len(time) and 0 <= s <= 1
c     so that extrapolation corresponds to end point mapping
c     istat = 0 for interior interpolation
c     istat = 1 for left extrapolation
c     istat = 2 for right extrapolation       
c     -----------------------------------------------------------------
      type(clock),intent(in)     :: when
      real,intent(in)            :: time(:)
      type(clock),intent(in)     :: time_offset
      integer,intent(out)        :: it0, it1 ! lower/upper frames for interpolation
      real,intent(out)           :: s
      integer,intent(out)        :: istat
      real                       :: thours, dt
      integer                    :: nt,i
c     -----------------------------------------------------------------
      nt = size(time)
      call get_period_length_hour(time_offset, when, thours)
      call search_sorted_list(thours, time, it0) ! 0 <= it0 <= nt
      it1 = it0 + 1
      it0 = min(max(1, it0), nt)  ! confine to valid range
      it1 = min(max(1, it1), nt)  ! confine to valid range
      if     (thours < time(1)) then
         istat = 1
      elseif (thours > time(nt)) then
         istat = 2
      else
         istat = 0  ! interior point
      endif
      if (it0 /= it1) then          ! assume time grid is not degenerate
         dt = time(it1)-time(it0)
         s  = (thours-time(it0))/dt
         s  = min(1.0, max(0.0, s)) ! enforce cap 0 <= s <= 1
         if (dt<1e-8) istat = 3     ! signal degenerate interval  
      else
         s = 0.5
      endif
      end subroutine get_time_coordinates


      
      subroutine load_xy_frame_to_extbuf(xyt, rbuffer, iframe)
c     -----------------------------------------------------------------
c     Load frame iframe from xyt into preallocated external rbuffer
c      
c     If iframe is absent, it is assumed to be a static xy variable, 
c     which is loaded correspondingly
c      
c     Float types relies on build-in cast
c     integer types are transformed as x_real = scale*x_int + offset
c     -----------------------------------------------------------------
      type(coards_xyt_variable), intent(in)  :: xyt
      real, intent(out)                      :: rbuffer(:,:)
      integer, intent(in), optional          :: iframe
c
      character*256                          :: msg
      integer, allocatable                   :: ibuffer(:,:)  ! buffer for integer-stored floats
      integer                                :: k,ix,iy,jx,jy
      integer                                :: start(3),count(3)
c     -----------------------------------------------------------------
      if ((size(rbuffer,1) /= xyt%grid%nx).or.
     +    (size(rbuffer,2) /= xyt%grid%ny)) then
         call NetCDFstop("rbuffer shape unsuitable")
      endif
      
      if (present(iframe)) then
         if ((iframe < 1).or.(iframe > xyt%nt))   
     +        call NetCDFstop("invalid frame requested")
         
         count = (/xyt%grid%nx, xyt%grid%ny, 1/)
         start = (/1,           1,           iframe/)
          if (xyt%storage_type == float_type) then ! just load, rely on build-in recast for different floats
             call NetCDFcheck( nf90_get_var(xyt%ncid, xyt%varid,
     +                     rbuffer, count=count, start=start) )
          else                  ! recast integer as x_real = scale*x_int + offset
             allocate( ibuffer(xyt%grid%nx, xyt%grid%ny) )
             call NetCDFcheck( nf90_get_var(xyt%ncid, xyt%varid,
     +                     ibuffer, count=count, start=start) )
          endif
       else     ! load static frame
          if (xyt%storage_type == float_type) then ! just load, rely on build-in recast for different floats
             call NetCDFcheck(nf90_get_var(xyt%ncid, xyt%varid,rbuffer))
          else                  
             allocate( ibuffer(xyt%grid%nx, xyt%grid%ny) )
             call NetCDFcheck(nf90_get_var(xyt%ncid, xyt%varid,ibuffer))
          endif
       endif

c      ------ recast integer as x_real = scale*x_int + offset, if loaded   
       if ( allocated(ibuffer) ) then
          rbuffer = xyt%scale*ibuffer + xyt%offset ! whole array transformation
          deallocate( ibuffer )
       endif
       
c      ------ flip axes (before remapping), if needed
       
       if (xyt%grid%flipx) call flip_xaxis(rbuffer) ! inplace flipping
       if (xyt%grid%flipy) call flip_yaxis(rbuffer) ! inplace flipping
       
c      ------ remap grid points, if feature set  ------  
       if ( associated(xyt%grid%remap) ) then
          do k = 1, size(xyt%grid%remap, 2)  ! shape = (4, nremaps) 
             ix = xyt%grid%remap(1,k) ! unpack invalid pt
             iy = xyt%grid%remap(2,k) ! unpack invalid pt
             jx = xyt%grid%remap(3,k) ! unpack source
             jy = xyt%grid%remap(4,k) ! unpack source
             rbuffer(ix,iy) = rbuffer(jx,jy)
          enddo
       endif

       
       end subroutine load_xy_frame_to_extbuf

      
      subroutine set_remap_nearest_valid(mask, remap)
c     -----------------------------------------------------------------
c     Given mask(nx,ny) generate a mapping (ix,iy) -> (jx,jy) for false points
c     in mask, where (jx,jy) represents nearest valid (true) for invalid point (ix,iy)
c      
c     This subroutine scans for nearest Euclidian point, assuming
c     isotropic grid. It is assumed that mask contains at least one valid point.
c     remap has been allocated prior to this subroutine
c     (min size = count(.not.mask))
c     Storage sequence: k remapping is remap(k,:) = (ix,iy, jx,jy)
c     -----------------------------------------------------------------
      logical, intent(in)  :: mask(:,:)    ! (nx,ny)
      integer, intent(out) :: remap(:,:)   ! (4,nremap)
      integer              :: nx,ny,ix,iy,jx,jy,ire
      real                 :: mindist2, dist2
c     -----------------------------------------------------------------
      nx  = size(mask,1)
      ny  = size(mask,2)
      ire = 1                   ! remap counter
      do iy=1,ny
         do ix=1,nx
            if (mask(ix,iy)) cycle   ! this points is valid, skip
            mindist2 = 1.0e12
            remap(1,ire) = ix        ! (ix,iy) is invalid
            remap(2,ire) = iy   ! (ix,iy) is invalid
c           ----- scan entire grid for nearest neighbor -----
            do jy=1,ny
               do jx=1,nx
                  if (.not.mask(jx,jy)) cycle     ! this includes (jx,jy) = (ix,iy)
                  dist2 = (ix-jx)**2 + (iy-jy)**2 ! squared Euclidian grid distance
                  if (dist2 < mindist2) then      ! store this candidate
                     mindist2     = dist2
                     remap(3,ire) = jx
                     remap(4,ire) = jy
                  endif
               enddo
            enddo
            ire = ire + 1 ! ready for next
         enddo
      enddo
      end subroutine set_remap_nearest_valid

      
      
      subroutine interpolate_xyt_variable(geo,xyt,deriv,result,status)
c     -----------------------------------------------------------------
c     Interpolate variable xyt at position geo
c      
c     This interpolation does not distinguish between wet/dry but only
c     considers the grid mask for including valid points
c     2D masked grid interpolation is based on weights ~ 1/dist1 for valid grid points
c     where dist1 is the 1-norm distance between interpolation point and grid point
c     Intracell coordinates satisfy 0 < (sx,sy) < 1 inside grid domain, but are unrestricted 
c     outside - this gives a more smooth extrapolation
c      
c     Applies the databuffer xyt%data, which is assumed syncronized to desired time,
c     for data interpolation
c     Interior grid range is considered to be (not same as wet definition)
c       0  < (x-lon1)/dlon < nx-1
c       0  < (y-lat1)/dlat < ny-1
c     i.e. an interpolation point must be surrounded by grid points (which
c     may have mask value false, though)
c
c     deriv  = 0: value interpolation in xyt%data
c     deriv  = 1: longitude derivative of value interpolation
c     deriv  = 2: latitude derivative of value interpolation
c      
c     status = 0: interior interpolation performed
c     status = 1: horizontal range violation, attempt extrapolation (superseedes status == 3 if both applies)
c     status = 3: rank deficit situation, set result = padval
c     -----------------------------------------------------------------     
      real, intent(in)                      :: geo(:)  ! only geo(1:2) is assessed
      type(coards_xyt_variable), intent(in) :: xyt
      integer, intent(in)                   :: deriv   ! 0=value; deriv = (1,2) gives derivative wrt (x,y) 
      real, intent(out)                     :: result
      integer, intent(out)                  :: status
      real, parameter                       :: padval = 0.0
      integer                               :: ix,iy,ix0,iy0,ix1,iy1
      real                                  :: x,y, sx,sy, w, d1
      real                                  :: val,wsum,dval,dwsum,dd1dg
c     -----------------------------------------------------------------      
      status = 0                     ! basic assertion
      x = (geo(1)-xyt%grid%lon1) / xyt%grid%dlon  ! unrestricted grid coordinate
      y = (geo(2)-xyt%grid%lat1) / xyt%grid%dlat  ! unrestricted grid coordinate
      sx  = x - int(x)  ! intra cell coordinate 0<sx<1
      sy  = y - int(y)  ! intra cell coordinate 0<sy<1
      if ((x<0).or.(x>(xyt%grid%nx-1))) status = 1 ! flag horizontal range violation, multiplicative
      if ((y<0).or.(y>(xyt%grid%ny-1))) status = 1 ! flag horizontal range violation, multiplicative
c
      ix  = 1 + int(x)           ! W grid point for interpolation
      iy  = 1 + int(y)           ! S grid point for interpolation
      ix0 = min(max(ix,    1),xyt%grid%nx)  ! resolve interpolation vertices
      ix1 = min(max(ix + 1,1),xyt%grid%nx)  ! resolve interpolation vertices
      iy0 = min(max(iy,    1),xyt%grid%ny)  ! resolve interpolation vertices
      iy1 = min(max(iy + 1,1),xyt%grid%ny) ! resolve interpolation vertices
      sx  = x - ix0 + 1  ! intra cell coordinate 0<sx<1 for interior, unrestricted outside
      sy  = y - iy0 + 1  ! intra cell coordinate 0<sy<1 for interior, unrestricted outside
c
c      write(*,*) iy0,iy1,sy
      
      val   = 0.0
      wsum  = 0.0
      dwsum = 0.0
      dval  = 0.0
      do iy = iy0, iy1
         do ix = ix0, ix1
            if ( xyt%grid%mask(ix,iy) ) then
               d1 = abs(sx-ix+ix0) + abs(sy-iy+iy0) + 1e-12 ! strictly d1 > 0
               if (deriv==1) dd1dg = sign(1.0/xyt%grid%dlon, sx-ix+ix0)  ! d(d1)/d(lon)
               if (deriv==2) dd1dg = sign(1.0/xyt%grid%dlat, sy-iy+iy0)  ! d(d1)/d(lat)
               w     = 1.0/d1
               val   = val + xyt%data(ix,iy)*w
               wsum  = wsum + w
               if (deriv>0) then
                 dval  = dval  - dd1dg * xyt%data(ix,iy) * w**2 
                 dwsum = dwsum - dd1dg * w**2
               endif
            endif
         enddo
      enddo
      
      if ((wsum==0.0).and.(status==0)) status = 3 ! do not superseede horizontal range violation flag
      
      if (deriv>0) then       ! derivative of interpolation
         if (wsum==0.0) then
            result = 0.0      ! derivative of padval
         else
            result = (dval*wsum - val*dwsum)/wsum**2
         endif
      else                    ! value interpolation
         if (wsum==0.0) then
            result = padval
         else
            result = val/wsum
         endif
      endif

      end subroutine interpolate_xyt_variable
     
      
      subroutine close_simple_lonlat_grid(grid)
      type(simple_lonlat_grid),intent(inout) :: grid
      if (associated(grid%mask))  deallocate(grid%mask)
      if (associated(grid%remap)) deallocate(grid%remap)
      end subroutine close_simple_lonlat_grid

      
      subroutine close_coards_xyt_variable(xyt)
      type(coards_xyt_variable),intent(inout) :: xyt
      call close_simple_lonlat_grid(xyt%grid)
      if (associated(xyt%time))      deallocate(xyt%time)
      if (associated(xyt%databuf0))  deallocate(xyt%databuf0)
      if (associated(xyt%databuf1))  deallocate(xyt%databuf1)
      end subroutine close_coards_xyt_variable


      subroutine flip_xaxis_real(buffer)
c     ------ flip xaxis inplace of 2D real buffer------
      real,intent(inout) :: buffer(:,:)
      real               :: tmp(size(buffer,2))
      integer            :: i,nx,ny
      nx = size(buffer,1)
      ny = size(buffer,2)
      do i = 1, nx/2            ! integer division; OK for nx=1
         tmp              = buffer(i,:)
         buffer(i,:)      = buffer(nx+1-i,:)
         buffer(nx+1-i,:) = tmp
      enddo
      end subroutine flip_xaxis_real

      
      subroutine flip_xaxis_logical(buffer)
c     ------ flip xaxis inplace of 2D logical buffer------
      logical,intent(inout) :: buffer(:,:)
      logical               :: tmp(size(buffer,2))
      integer            :: i,nx,ny
      nx = size(buffer,1)
      ny = size(buffer,2)
      do i = 1, nx/2            ! integer division; OK for nx=1
         tmp              = buffer(i,:)
         buffer(i,:)      = buffer(nx+1-i,:)
         buffer(nx+1-i,:) = tmp
      enddo
      end subroutine flip_xaxis_logical

      
      subroutine flip_yaxis_real(buffer)
c     ------ flip xaxis inplace of 2D real buffer------
      real,intent(inout) :: buffer(:,:)
      real               :: tmp(size(buffer,1))
      integer            :: i,nx,ny
      nx = size(buffer,1)
      ny = size(buffer,2)
      do i = 1, ny/2            ! integer division; OK for nx=1
         tmp              = buffer(:,i)
         buffer(:,i)      = buffer(:,ny+1-i)
         buffer(:,ny+1-i) = tmp
      enddo   
      end subroutine flip_yaxis_real

      
      subroutine flip_yaxis_logical(buffer)
c     ------ flip xaxis inplace of 2D logical buffer------
      logical,intent(inout) :: buffer(:,:)
      logical               :: tmp(size(buffer,1))
      integer            :: i,nx,ny
      nx = size(buffer,1)
      ny = size(buffer,2)
      do i = 1, ny/2            ! integer division; OK for nx=1
         tmp              = buffer(:,i)
         buffer(:,i)      = buffer(:,ny+1-i)
         buffer(:,ny+1-i) = tmp
      enddo         
      end subroutine flip_yaxis_logical


      
      subroutine NetCDFcheck(status)
c     -----------------------------------------------------------------
c     This subroutine supports the recommended reading 
c     style of NetCDF4:
c       call NetCDFcheck( nf90_get_var(ncid, varid, data_in) )
c     -----------------------------------------------------------------
      integer, intent (in) :: status
c     -----------------------------------------------------------------    
      if(status /= nf90_noerr) then 
         print *, trim(nf90_strerror(status))
         call NetCDFstop("NetCDFcheck:Stopped")
      end if
      end subroutine NetCDFcheck

      
      subroutine NetCDFstop(message)
c     -----------------------------------------------------------------
c     Fatal error condition encountered, stop netcdf handling
c     -----------------------------------------------------------------
      character*(*) :: message
c     -----------------------------------------------------------------    
      write(*,744) trim(message)
      stop "NetCDFstop:Stopped"
 744  format("Fatal netCDF condition:", a)
      end subroutine NetCDFstop

      
      end module 


      
c      program test
c      use netcdf             ! access raw API
c      use coards_netcdf_2D   ! add-on to netcdf
c      use time_tools   
c      use array_tools
c      implicit none 
c      character*(*), parameter  :: path = "ERA5_test_wave.nc"
c      integer                   :: ncid, status, ix,iy,i
c      real                      :: geo(2),res,v1,v0,dxy=0.001
c      type(coards_xyt_variable) :: xyt
c      type(clock)               :: when
cc     --------------------------------------------------------------------     
c      status = nf90_open(path, nf90_nowrite, ncid)  ! netcdf API
c      
cc     ---- instantiate grid + variable ----
c      call attach_variable(xyt,"shww",ncid,"native") ! args == (xyt, varname, ncid, cast_mode)
c      call attach_variable(xyt,"shww",ncid,"nearest_valid")  ! alternative cast
c      
cc     ---- specific data frame fetch ----      

c      call load_frame(xyt, iframe=7)
c      do ix = 1, xyt%grid%nx
c         do iy = 1, xyt%grid%ny
c            if (xyt%grid%mask(ix,iy)) write(21,*) ix,iy
c         enddo
c     enddo
cc     ---- interpolate data to a specific time ----
c      call set_clock(when, 2014, 2, 23, 3*3600)
c      call sync_to_time(xyt, when)
c      call write_valid_points(xyt%grid,23)
c      
cc     ---- do a line scan ----      
c      do i=1,1000
c         geo = (/4., 45.+i*0.02/)
c         geo = (/-5. + 0.015*i, 55./)   
c         call interpolate_xyt_variable(geo,xyt,0,res,status) ! args == (geo,xyt,deriv,result,status)
c         write(22,*) geo(1), res, status
c     enddo
c      
cc     ---- test derivative consistency ----            
c      geo = (/3.046, 55.127/)
c      call interpolate_xyt_variable(geo,xyt,1,res,status) ! args == (geo,xyt,deriv,result,status)
c      geo(1) = geo(1) + dxy
c      call interpolate_xyt_variable(geo,xyt,0,v1,status) 
c      geo(1) = geo(1) - 2*dxy
c      call interpolate_xyt_variable(geo,xyt,0,v0,status)
c      write(*,*) "dv_dx(analyt) =", res
c      write(*,*) "dv_dx(numeric)=", (v1-v0)/2/dxy
c      
c      geo = (/3.046, 55.127/)      
c      call interpolate_xyt_variable(geo,xyt,2,res,status) ! args == (geo,xyt,deriv,result,status)
c      geo(2) = geo(2) + dxy
c      call interpolate_xyt_variable(geo,xyt,0,v1,status) 
c      geo(2) = geo(2) - 2*dxy
c      call interpolate_xyt_variable(geo,xyt,0,v0,status)
c      write(*,*) "dv_dy(analyt) =", res
c      write(*,*) "dv_dy(numeric)=", (v1-v0)/2/dxy
      
c      call close_coards_xyt_variable(xyt)  ! includes grid
c      status = nf90_close(ncid) ! netcdf API
c      end program
