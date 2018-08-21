ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     COARDS like format with all frames in single file
c     ---------------------------------------------------
c
c     Currently
c         pick time frame closest to current time
c         data on same and static vertical z-grid layer (no time varying sea surface elevation added)
c         only physics
c         assume all data on same grid and at same times
c         only regular lon-lat is allowed horizontally
c       
c     NOTES: based on BIMS-ECO/POM_daily_mean.f
c     TODO: add a input handle to allow flipping w (if present) to
c           conform to IBMlib orientation (positive downward)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module physical_fields

      use mesh_grid   ! use cc current interpolation
      use horizontal_grid_transformations
      use geometry
      use array_tools
      use run_context, only: simulation_file
      use time_services           ! import clock type and time handling
      use input_parser
      use netcdf   ! use site installation

      implicit none
      private     


      public :: init_physical_fields     
      public :: close_physical_fields
      public :: get_master_clock ! from time_services
      public :: set_master_clock ! from time_services
      public :: update_physical_fields
      public :: interpolate_turbulence
      public :: interpolate_turbulence_deriv
      public :: interpolate_currents
      public :: interpolate_temp
      public :: interpolate_salty   

      public :: interpolate_wdepth
      public :: is_wet    
      public :: is_land
      public :: horizontal_range_check
      public :: coast_line_intersection
      public :: get_pbi_version

c     -------------------- module data --------------------  
      
      
c
c     ------ data frame handler ------
      
      character*999              :: hydroDBpath ! hydrographic data sets (proper trailing separator will be added at read, if not provided)
      type(clock)                :: time_offset ! 1970-01-01 00:00:00 - this is a particular to BIMS-ECO and should be removed from template ...
      
      
      character*(*), parameter   :: not_set_c = "none"     ! assume name conflict unlikely
      integer, parameter         :: not_set_i = -999999999 ! assume number conflict unlikely
                                                
      logical                    :: data_in_buffers     = .false. ! first time setting
      logical                    :: data_files_are_open = .false. ! first time setting   
      integer                    :: cur_frame    ! current frame number in buffer
      integer                    :: ncid         ! NetCDF file handler for currents (this is considered the root data set)
      integer                    :: ncid_salt    ! NetCDF file handler for salinity (possibly same as ncid )
      integer                    :: ncid_temp    ! NetCDF file handler for temperature (possibly same as ncid )        
      logical                    :: w_defined    ! vertical current available else set to zero
      integer                    :: varid_u, varid_v, varid_w ! NetCDF varIDs for current components
      integer                    :: varid_temp, varid_salt    ! NetCDF varIDs for temperature and salinity
                       
      real, parameter            :: molecular_diffusivity = 1.e-9    ! unit m2/s (lower limit on applied horizontal/vertical diffusivity)
      real                       :: constant_horiz_diffusivity       ! unit m2/s; 
      real                       :: constant_verti_diffusivity       ! unit m2/s; 

      real, allocatable          :: ccdepth0(:) ! static vertical grid
c
c     ---------- time grid supported is a regularly spaced sequence described by these three parameters:
c
      integer                    :: time1          ! seconds since time offset of first time frame
      integer                    :: dtime          ! time step between subsequent time frames
      integer                    :: nt             ! number of time frames in data set       

c     --- units    ---    
 

c     ===================================================
                            contains
c     ===================================================
  
      subroutine init_physical_fields(time)
c     ------------------------------------------
c     Do not trigger data load
c     primitive runtime linux/windows recognition, by allowing user to supply path separator 
c     ------------------------------------------
      type(clock), intent(in),optional :: time
      character(len=len(hydroDBpath))  :: path
      character                        :: lastchr
      integer                          :: lp, varid
      character*999                    :: toffstr
      integer                          :: year,month,day
      integer                          :: hour,minute,sec
      integer                          :: sec_of_day
c     ------------------------------------------
      if (present(time)) then
         call set_master_clock(time)
         write(*,*) "set simulation PBI clock = "
         call write_clock(time)
      endif
      write(*,*) trim(get_pbi_version()) 

      call read_control_data(simulation_file,"hydroDBpath", path)
      
c     --- provide primitive runtime linux/windows recognition, by allowing user to 
c         supply path separator to hydroDBpath

      lp      = len(path)
      path    = adjustr(path)
      lastchr = path(lp:lp)
      if     (lastchr == "/") then                  ! found explicit linux seperator
         hydroDBpath = trim(adjustl(path))          ! keep seperator as provided
      elseif (lastchr == "\") then                  ! found explicit windows seperator
         hydroDBpath = trim(adjustl(path))          ! keep seperator as provided
      else
         hydroDBpath = trim(adjustl(path)) // "/"   ! no separator, assume linux and prepend linux delimiter
      endif

      write(*,*) "init_physical_fields: hydrographic database path =", 
     +           trim(adjustl(hydroDBpath)) 
      

      call open_data_files()  ! initializes netCDF handlers
      call reset_frame_handler()
      call load_grid_desc()  ! invokes init_horiz_grid_transf
      call init_mesh_grid(init_biogeochem = .false.)  ! allocate core arrays, currently no biogeochem
      call init_topography() ! requires mesh_grid arrays are allocated


c
c     resolve the time offset for the time variable in the data set
c     by parsing the attribute "time:units"
c
      call NetCDFcheck( nf90_inq_varid(ncid, "time", varid) )   ! fixed name
      call NetCDFcheck( nf90_get_att(ncid, varid,                
     +     "units", toffstr) )                                  ! retrieve time:units attribute
      
c      0...0....1....1....2....2....3
c      1...5....0....5....0....5....0
c     "seconds since 1970-01-01T00:00:00Z"
c     "seconds since 1970-01-01 00:00:00"
 711  format(14x,i4,5(1x,i2))

      toffstr = adjustl(toffstr)
      read(toffstr,711) year,month,day,hour,minute,sec
      sec_of_day = 3600*hour + 60*minute + sec
      call set_clock(time_offset, year, month, day, sec_of_day) 
      write(*,*) "init_physical_fields: data time offset is"
      call write_clock(time_offset)
    

c
c     ---- currently no horizontal turbulent diffusivity in data set ----
c          allow user to set a static horizontal turbulent diffusivity for hdiffus  (optional)  
c          initialize it to molecular diffusivity lower limit, if no values 
c          are provided
c    
      if (count_tags(simulation_file, "horizontal_diffusivity")/=0) then
         call read_control_data(simulation_file,
     +                          "horizontal_diffusivity",
     +                          constant_horiz_diffusivity)
         
         write(*,563) "horizontal", constant_horiz_diffusivity
         hdiffus = constant_horiz_diffusivity
      else
         constant_horiz_diffusivity = -999.0  ! set handle to apply dynamic horiz_diffusivity, if available
         write(*,563) "horizontal", molecular_diffusivity
         hdiffus = molecular_diffusivity
      endif

c
c     ---- currently no vertical turbulent diffusivity in data set ----
c          allow user to set a static vertical turbulent diffusivity for vdiffus  (optional)  
c          initialize it to molecular diffusivity lower limit, if no values 
c          are provided
c    
      if (count_tags(simulation_file, "vertical_diffusivity")/=0) then
         call read_control_data(simulation_file,
     +                          "vertical_diffusivity",
     +                          constant_verti_diffusivity)
         
         write(*,563) "vertical", constant_verti_diffusivity
         vdiffus = constant_verti_diffusivity
      else
         constant_verti_diffusivity = -999.0  ! set handle to apply dynamic verti_diffusivity, if available
         write(*,563) "vertical", molecular_diffusivity
         vdiffus = molecular_diffusivity
      endif


  563 format("init_physical_fields: using constant ",a,
     +       " diffusivity =", e12.5," m2/s")


      
      end subroutine init_physical_fields

 

      character*100 function get_pbi_version()  
      get_pbi_version="COARDS simple offline pbi version: $Rev:  $"
      end function



      subroutine close_physical_fields()
c     ------------------------------------------ 
c     ------------------------------------------    
      call close_mesh_grid()
      call reset_frame_handler()
      call close_horiz_grid_transf()
c     --- local cleanup ---
      call close_data_files()
      deallocate( ccdepth0 )
c     ------------------------------------------
      end subroutine 



      subroutine update_physical_fields(time, dt)
c     ------------------------------------------  
      type(clock), intent(in),optional :: time
      integer, intent(in),optional     :: dt
      
      integer                          :: frame   ! time frame number, not sec since offset
      type(clock), pointer             :: aclock
c     ------------------------------------------  
      aclock => get_master_clock()
      if (present(time)) then
         call set_master_clock(time)
      elseif (present(dt)) then
         call add_seconds_to_clock(aclock, dt)
      endif
c
      call resolve_corresp_frame(aclock, frame)

      if ((.not.data_in_buffers).or.(frame /= cur_frame)) then
         call load_data_frame(frame)  ! updates state handlers
      endif
      
c     ------------------------------------------       
      end subroutine update_physical_fields



      
      subroutine load_grid_desc()
c     -------------------------------------------------------
c     Should only be invoked at start to pick up grid layout
c     assume grid is regular lon-lat oriented W->E, S->N 
c     invoke init_horiz_grid_transf
c     root netCDF data set must be opened
c     dimension and variable name of longitude and latitude must coinside
c     -------------------------------------------------------
      character*64         :: lon_name, lat_name ! assume long enough

      integer              :: dimid,varid,idum, gncid, ibuf(2)
      real                 :: lambda1, dlambda, phi1, dphi, rbuf(2) ! LOCAL DUMMIES
c     -------------------------------------------------------
      write(*,*) "load_grid_desc: loading grid from root data set"

      call read_control_data(simulation_file, "longitude_name", 
     +                       lon_name)
      call NetCDFcheck( nf90_inq_dimid(ncid, lon_name, dimid) )
      call NetCDFcheck( nf90_inquire_dimension(ncid, dimid, len=nx) )

      call read_control_data(simulation_file, "latitude_name", 
     +                       lat_name)
      call NetCDFcheck( nf90_inq_dimid(ncid, lat_name, dimid) )
      call NetCDFcheck( nf90_inquire_dimension(ncid, dimid, len=ny) )

c     ---  only vertical dim/var name "depth" is recognized
      call NetCDFcheck( nf90_inq_dimid(ncid, "depth", dimid) )
      call NetCDFcheck( nf90_inquire_dimension(ncid, dimid, len=nz) )
      
c      
c     --- probe horizontal grid descriptors from lon/lat arrays
c
      
      call NetCDFcheck( nf90_inq_varid(ncid, lon_name, varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, rbuf, count=(/2/)) )
      lambda1 = rbuf(1)
      dlambda = rbuf(2)-rbuf(1)   ! regular grid is assumed

      call NetCDFcheck( nf90_inq_varid(ncid, lat_name, varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, rbuf, count=(/2/)) )
      phi1    = rbuf(1)
      dphi    = rbuf(2)-rbuf(1)   ! regular grid is assumed
    
      write(*,231) nx, ny, nz
      write(*,232) lambda1, dlambda
      write(*,233) phi1,    dphi 

      write(*,*) "load_grid_desc: horizontal initialization"
      call init_horiz_grid_transf(lambda1, phi1, dlambda, dphi)

c      
c     --- assess vertical grid descriptor from depth arrays ---
c   
      allocate( ccdepth0(nz) )
      call NetCDFcheck( nf90_inq_varid(ncid, "depth", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, ccdepth0))
      write(*,241) "vertical", ccdepth0

  229 format(a,a)    
  231 format("load_grid_desc: 3d grid dim (nx,ny,nz) = ", i4,i4,i4)
  232 format("load_grid_desc: lambda1 = ",f12.7,
     +                      " dlambda = ",f12.7," degrees")  
  233 format("load_grid_desc: phi1    = ",f12.7,
     +                      " dphi    = ",f12.7," degrees")  
      
 241  format("load_grid_desc: static ", a, " grid = ", 999f9.3)

c      
c     --- assess time grid available (load it as integer)---
c  
      call NetCDFcheck( nf90_inq_dimid(ncid, "time", dimid) )   ! fixed name
      call NetCDFcheck( nf90_inquire_dimension(ncid, dimid, len=nt) )
      call NetCDFcheck( nf90_inq_varid(ncid, "time", varid) )   ! fixed name
      call NetCDFcheck( nf90_get_var(ncid, varid, ibuf, count=(/2/)) )
      time1 = ibuf(1)
      dtime = ibuf(2)-ibuf(1)   ! regular grid is assumed
      write(*,*) "load_grid_desc: data set contains", nt, "time frames"
      write(*,*) "load_grid_desc: time resolution = ", dtime, "sec"

      end subroutine load_grid_desc
      



      subroutine init_topography()
c     ------------------------------------------------------------
c     Assess topography by probing first temperature frame  
c     Set mesh_grid arrays (which are assumed initialized):
c        wetmask      (static)
c        bottom_layer (static)
c        wdepth       (static)
c        ccdepth      (static, no sea elevation)
c        acc_width    (static, no sea elevation)
c     If dry points cannot be detected, all points are assumed wet
c     ------------------------------------------------------------
      integer             :: dimid,varid,idum, gnci
      integer             :: ix,iy,iz,count4D(4),istat
      real,allocatable    :: buffer(:,:,:), acc_width0(:)
      character*64        :: temp_name
      real                :: dry_val, landfrac, wetfrac, layerwidth 
c     ------------------------------------------------------------
      allocate( buffer(nx,ny,nz) )
      count4D = (/nx,ny,nz,1/)
      call read_control_data(simulation_file, "temperature_name", 
     +                       temp_name)
      call NetCDFcheck( nf90_inq_varid(ncid_temp, temp_name, varid) )
      call NetCDFcheck( nf90_get_var(ncid_temp, varid, buffer,
     +                       count=count4D))

c
c     check if we can retrieve a missing value attribute, otherwise assume
c     nf90_fill_real is applied for dry points 
c
      istat = nf90_get_att(ncid_temp, varid,"missing_value", dry_val)  ! does possibly not exist
      if (istat /= nf90_noerr) then
          dry_val = nf90_fill_real  
          write(*,*) "init_topography: warning: no missing_value "//
     +         "in data set - assuming nf90_fill_real applies"
       else
          write(*,*) "init_topography: dry value identified"
      endif
      do ix = 1, nx
         do iy = 1, ny
            do iz = 1, nz
               if (abs(buffer(ix,iy,iz))>0.999*abs(dry_val)) exit
            enddo
            bottom_layer(ix,iy) = iz-1   ! last wet layer in this column
         enddo
      enddo

      where(bottom_layer > 0)
         wetmask      = 1       ! wet
      elsewhere
         wetmask      = 0       ! dry
      end where
      
      landfrac = 1.0 - 1.0*sum(wetmask)/nx/ny
      wetfrac  = 1.0*sum(bottom_layer)/nx/ny/nz
      write(*,*) "init_topography: land fraction        = ", landfrac
      write(*,*) "init_topography: wet volume fraction  = ", wetfrac
c
c     initialize dynamic 3D topographuc arrays ccdepth, acc_width and wdepth
c     which are static for this PBI
c  
      allocate ( acc_width0(nz+1) )
      acc_width0(1) = 0.0
      layerwidth    = 2*ccdepth0(1)
      acc_width0(2) = layerwidth
      do iz = 2, nz
         layerwidth = 2.0*(ccdepth0(iz)-ccdepth0(iz-1)) - layerwidth
         acc_width0(iz+1) = acc_width0(iz) + layerwidth
      enddo
      do ix = 1, nx
         do iy = 1, ny
            ccdepth(ix,iy,:)   = ccdepth0
            acc_width(ix,iy,:) = acc_width0
            wdepth(ix,iy)      = acc_width0(1 + bottom_layer(ix,iy)) ! OK for dry
         enddo
      enddo

      deallocate( buffer )
      deallocate( acc_width0 )

      end subroutine init_topography




      subroutine reset_frame_handler()
c     ------------------------------------------
      data_in_buffers = .false.
      cur_frame       = not_set_i
      end subroutine reset_frame_handler



      subroutine close_data_files()
c     ------------------------------------------
c     Only close, if data files are open
c     ------------------------------------------
      if (data_files_are_open) then
         call NetCDFcheck( nf90_close(ncid) )
         ncid = not_set_i          
         write(*,*) "close_data_files: closed open NetCDF set"
         data_files_are_open = .false.
      endif
      end subroutine close_data_files




      subroutine open_data_files()
c     ---------------------------------------------------
c     Explore data availability in hydroDBpath 
c     and open data files and set file handlers and prepare reading
c     Either there is one multi file corresponding to input tag (precedence)
c
c          hydrography_file
c
c     or three files corresponding to input tags
c
c          current_file
c          temperature_file
c          salinity_file
c
c     In any case, netcdf handlers (ncid, ncid_salt, ncid_temp) 
c     are defined, in the former case they are identical
c     ---------------------------------------------------
      character*999    :: fname ! assume long enough
      integer          :: timedimID, ntime, varid
      character*64     :: varname
c     ------------------------------------------
      call close_data_files() ! in case they are open ...

      if (count_tags(simulation_file, "hydrography_file")/=0) then
c     --- combifile input ---
           call read_control_data(simulation_file,"hydrography_file",
     +                       fname)
           fname = trim(hydroDBpath) // trim(adjustl(fname))  ! incl separator
           call NetCDFcheck( nf90_open(fname, NF90_NOWRITE, ncid))
           ncid_salt = ncid
           ncid_temp = ncid
           write(*,352) trim(fname)              
      else  
c     --- input split in three files ---
         call read_control_data(simulation_file,"current_file",
     +                       fname)
         fname = trim(hydroDBpath) // trim(adjustl(fname))  ! incl separator
         call NetCDFcheck( nf90_open(fname, NF90_NOWRITE, ncid))  
         write(*,353) "currents", trim(fname)
c
         call read_control_data(simulation_file,"temperature_file",
     +                       fname)
         fname = trim(hydroDBpath) // trim(adjustl(fname))  ! incl separator
         call NetCDFcheck( nf90_open(fname, NF90_NOWRITE, ncid_temp))  
         write(*,353) "temperature", trim(fname)
c
         call read_control_data(simulation_file,"salinity_file",
     +                       fname)
         fname = trim(hydroDBpath) // trim(adjustl(fname))  ! incl separator
         call NetCDFcheck( nf90_open(fname, NF90_NOWRITE, ncid_salt))  
         write(*,353) "salinity", trim(fname)
      endif

      data_files_are_open = .true.

 352  format("open_data_files: opened combi file: ", a)
 353  format("open_data_files: opened separate file (",a,"):",a) 

c
c     ----- resolve which variable names to look for in the netcdffile
c           and retrieve the corresponding variable IDs

      call read_control_data(simulation_file, "u_name", varname)
      call NetCDFcheck( nf90_inq_varid(ncid, varname, varid_u) )
      write(*,362) "u", trim(varname), varid_u   
c
      call read_control_data(simulation_file, "v_name", varname)
      call NetCDFcheck( nf90_inq_varid(ncid, varname, varid_v) )
      write(*,362) "v", trim(varname), varid_v
c     
      if (count_tags(simulation_file, "w_name")/=0) then
         call read_control_data(simulation_file, "w_name", varname)
         call NetCDFcheck( nf90_inq_varid(ncid, varname, varid_w) )
         write(*,362) "w", trim(varname), varid_v
         w_defined = .true.
      else
         write(*,363) 
         w_defined = .false.
      endif
c
      call read_control_data(simulation_file, "salinity_name", varname)
      call NetCDFcheck( nf90_inq_varid(ncid_salt, varname, varid_salt))
      write(*,362) "salinity", trim(varname), varid_salt
c     
      call read_control_data(simulation_file, "temperature_name", 
     +                       varname)
      call NetCDFcheck( nf90_inq_varid(ncid_temp, varname, varid_temp))
      write(*,362) "temperature", trim(varname), varid_temp
c  
 362  format("open_data_files: varname for ", a," = ", a, 
     +       " (varid=", i2, ")") 
 363  format("open_data_files: w not available")
 

      end subroutine open_data_files



      subroutine resolve_corresp_frame(aclock, frame)
c     ------------------------------------------------------------------------------------
c     Resolve which time frame (1 <= frame <= nt) to use in data set corresponding to aclock
c     based on time grid descriptors (time1, dtime, nt) resolved at initialization
c     In case of exterior time points, the first/last frame will be applied
c     ------------------------------------------------------------------------------------
      type(clock), intent(in) :: aclock
      integer, intent(out)    :: frame
      integer                 :: numsec, it
      real                    :: t
c     ------------------------------------------------------------------------------------
      call get_period_length_sec(time_offset, aclock, numsec)
      it     = nint(1.0*(numsec-time1)/dtime) ! time grid 
      if ((it < 1) .or. (it > nt)) then
         write(*,286) it
      endif
 286  format("resolve_corresp_frame: warning: time extrapolation (it=",
     +    i4,")")
      frame = min(nt, max(it,1))
      end subroutine resolve_corresp_frame



      subroutine load_data_frame(frame)
c     ------------------------------------------------------------------------------------
c     Verbose reading af NetCDF sets to grid buffers from multi frame data set
c     
c     It is assumed that grids do not change during a simulation (assertion not checked)
c     This version does not feature dynamic sea surface elevation, so vertical auxillary
c     arrays are static.
c     Variable handler are retrieved at initialization
c
c     fortran index order opposite CDL dump order
c     fortran: fastests index left most
c     ------------------------------------------------------------------------------------
      integer, intent(in) :: frame
      integer :: count4D(4), start4D(4)
c     ------------------------------------------------------------------------------------ 
      write(*,*) "load_data_frame: loading frame ", frame


c     -------------------- 3D frames ------------------------

      start4D = (/1,1,1,frame/) 
      count4D = (/nx,ny,nz, 1/) 
      
      
c     ---- load u current component ---- 
      call NetCDFcheck( nf90_get_var(ncid, varid_u, u,
     +                  count=count4D, start=start4D) )
      
c     ---- load v current component ---- 
      call NetCDFcheck( nf90_get_var(ncid, varid_v, v,
     +                  count=count4D, start=start4D) )

c     ---- load w current component (if present) ----
      if (w_defined) then
          call NetCDFcheck( nf90_get_var(ncid, varid_w, w,
     +        count=count4D, start=start4D) )
      else
          w = 0.  ! default
      endif

c     ---- load temperature ----
      call NetCDFcheck( nf90_get_var(ncid_temp, varid_temp, temp,
     +                  count=count4D, start=start4D) )

c     ---- load salinity ---- 
      call NetCDFcheck( nf90_get_var(ncid_salt, varid_salt, salinity,
     +                  count=count4D, start=start4D) )
      
c     ------------ update state handlers ------------
      data_in_buffers = .true.
      cur_frame       = frame

      end subroutine load_data_frame
      


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     NetCDF auxillaries  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine NetCDFcheck(status)
c     ------------------------------------------
c     This subroutine supports the recommended reading 
c     style of NetCDF4:
c       call NetCDFcheck( nf90_get_var(ncid, varid, data_in) )
c     Private to this module
c     ------------------------------------------
      integer, intent ( in) :: status
    
      if(status /= nf90_noerr) then 
         print *, trim(nf90_strerror(status))
         stop "NetCDFcheck:Stopped"
      end if
      end subroutine NetCDFcheck 



      end module

      
