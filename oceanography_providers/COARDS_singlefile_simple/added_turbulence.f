cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ------------------------------------------------------------------------------
c     COARDS like format with all frames in single file(s) with posteori turbulence
c     ------------------------------------------------------------------------------
c
c     Currently
c        * optional zero-equation posteori turbulence schemes
c        * Pick time frame closest to current time
c        * Grid is either in z-mode or sigma-mode   
c        * Data on same and static vertical z/sigma-grid layer (no time varying sea surface elevation added)
c        * Only physics
c        * Assume all data on same grid and at same times
c        * Only regular lon-lat is allowed horizontally
c             
c       
c     NOTES: based on coards_allframes_simple.f
c     TODO: add a input handle to allow flipping w (if present) to
c           conform to IBMlib orientation (positive downward)
c           add dynamic sea surface elevation, if present
c           slim load (e.g. only currents)
c           topography from currents rather than temperature
c           surface wind (uswind,vswind) not loaded       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module physical_fields

      use mesh_grid,
     +     redir_t  => interpolate_turbulence, ! shadow by redirection
     +     redir_dt => interpolate_turbulence_deriv ! shadow by redirection
      use water_density, only:
     +     rhow,
     +     init_water_density,
     +     update_water_density,
     +     close_water_density
      use horizontal_grid_transformations
      use geometry
      use array_tools
      use run_context, only: simulation_file
      use time_services           ! import clock type and time handling
      use input_parser
      use coards_netcdf   ! 
      use turbulence,
     +     ext_turb_intp   => interpolate_turbulence, ! hide by redirection
     +     ext_turbderiv_intp => interpolate_turbulence_deriv ! hide by redirection
      implicit none
      private     


      public :: init_physical_fields     
      public :: close_physical_fields
      public :: get_master_clock ! from time_services
      public :: set_master_clock ! from time_services
      public :: update_physical_fields
      public :: interpolate_turbulence       ! defined locally, not in mesh_grid or turbulence
      public :: interpolate_turbulence_deriv ! defined locally, not in mesh_grid or turbulence
      public :: interpolate_currents
      public :: interpolate_temp
      public :: interpolate_salty   

      public :: interpolate_wdepth
      public :: is_wet    
      public :: is_land
      public :: horizontal_range_check
      public :: coast_line_intersection
      public :: get_pbi_version

c      public :: dump_land_points ! non standard tmp
      
c     -------------------- module data --------------------  
      
      
c
c     ------ data frame handler ------
      
      character*999              :: hydroDBpath ! hydrographic data sets (proper trailing separator will be added at read, if not provided)
      type(clock)                :: time_offset ! 
      
      
      character*(*), parameter   :: not_set_c = "none"     ! assume name conflict unlikely
      integer, parameter         :: not_set_i = -999999999 ! assume number conflict unlikely
                                                
      logical                    :: data_in_buffers     = .false. ! first time setting
      logical                    :: data_files_are_open = .false. ! first time setting   
      integer                    :: cur_frame    ! current frame number in buffer
      integer                    :: ncid         ! NetCDF file handler for currents (this is considered the root data set)
      integer                    :: ncid_salt    ! NetCDF file handler for salinity (possibly same as ncid )
      integer                    :: ncid_temp    ! NetCDF file handler for temperature (possibly same as ncid )        
      logical                    :: w_defined    ! vertical current available else set to zero
      type(coards_xyzt_variable) :: u_var        ! proxy var for easy loading/unpacking
      type(coards_xyzt_variable) :: v_var        ! proxy var for easy loading/unpacking
      type(coards_xyzt_variable) :: w_var        ! proxy var for easy loading/unpacking  
      type(coards_xyzt_variable) :: t_var        ! proxy var for easy loading/unpacking 
      type(coards_xyzt_variable) :: s_var        ! proxy var for easy loading/unpacking
                       
      logical                    :: zmode        ! grid mode toggle: TRUE : z-mode ... FALSE : sigma-mode 
      real, allocatable          :: ccdepth0(:)  ! static vertical grid (z-mode)
      real, allocatable          :: sigma(:)     ! sigma grid           (sigma-mode)
      real, allocatable          :: wdepth0(:,:) ! static CC depths     (sigma-mode)

      real                       :: lambda1, dlambda, phi1, dphi ! local copies for convenience
      real, allocatable          :: uswind(:,:)  ! W:E surface wind [m/s] (currently unused)
      real, allocatable          :: vswind(:,:)  ! S:N surface wind [m/s]  (currently unused)
      real, allocatable          :: layer_width(:,:,:)  ! generated from acc_width
c
c     ---------- time grid supported is a regularly spaced sequence described by these three parameters:
c
      integer                    :: time1             ! time since offset of first time frame (in data set time unit) 
      integer                    :: dtime             ! time step between subsequent time frames (in data set time unit) 
      integer                    :: nt                ! number of time frames in data set       
      integer                    :: sec_per_time_unit ! seconds per time unit of data set
      
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

      call resolve_time_parameters(ncid, "time", year, month, day, 
     +                             hour, minute, sec, sec_per_time_unit)
      write(*,*) "init_physical_fields: sec_per_time_unit=", 
     +            sec_per_time_unit
      sec_of_day = 3600*hour + 60*minute + sec
      call set_clock(time_offset, year, month, day, sec_of_day) 
      write(*,*) "init_physical_fields: data time offset is"
      call write_clock(time_offset)


      call load_grid_desc()  ! invokes init_horiz_grid_transf, set zmode
      call init_mesh_grid(init_biogeochem = .false.) ! allocate core arrays, currently no biogeochem
      if (zmode) then
         call init_topography_z()     ! requires mesh_grid arrays are allocated
      else
         call init_topography_sigma() ! requires mesh_grid arrays are allocated
      endif

c
c     ---- redelegate turbulence to external module ----
c    
      allocate( uswind(nx,ny) )
      allocate( vswind(nx,ny) )
      allocate( layer_width(nx,ny,nz) )
      uswind = 0.0
      vswind = 0.0
      call init_water_density()       ! follows init_mesh_grid
      call init_turbulence(nx,ny,nz, lambda1, phi1, dlambda, dphi, 
     +                     bottom_layer)                  
 

      end subroutine init_physical_fields

 

      character*100 function get_pbi_version()  
      get_pbi_version="COARDS simple offline pbi version: $Rev:  $"
      end function



      subroutine close_physical_fields()
c     ------------------------------------------ 
c     ------------------------------------------    
      call close_mesh_grid()
      call close_water_density()
      call reset_frame_handler()
      call close_horiz_grid_transf()    
      call close_turbulence()
c     --- local cleanup ---
      call close_data_files()
      if (allocated(ccdepth0)) deallocate(ccdepth0)
      if (allocated(sigma))    deallocate(sigma)
      if (allocated(wdepth0))  deallocate(wdepth0)
      if (allocated(uswind))       deallocate(uswind)
      if (allocated(vswind))       deallocate(vswind)
      if (allocated(layer_width))  deallocate(layer_width)
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
         call load_data_frame(frame)  ! updates state handlers + secondary updates
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
c
c     determine whether data is z or sigma vertically
c     sigma grid is assumed, if a depth_file is found
c     -------------------------------------------------------
      character*64         :: lon_name, lat_name ! assume long enough
      character*999        :: fname ! assume long enough
      integer              :: dimid,varid,idum, gncid, ibuf(2)
      real                 :: wd, rbuf(2) 
      integer              :: iunit, ix, iy, nlines
c     -------------------------------------------------------
      write(*,*) "load_grid_desc: loading grid from root data set"

      call read_control_data(simulation_file, "longitude_name", 
     +                       lon_name)
      call NetCDFcheck( nf90_inq_dimid(ncid, lon_name, dimid) )
      call NetCDFcheck( nf90_inquire_dimension(ncid, dimid, len=nx) )  ! resolve nx

      call read_control_data(simulation_file, "latitude_name", 
     +                       lat_name)
      call NetCDFcheck( nf90_inq_dimid(ncid, lat_name, dimid) )
      call NetCDFcheck( nf90_inquire_dimension(ncid, dimid, len=ny) )  ! resolve ny

c     ---  only vertical dim/var name "depth" is recognized
      call NetCDFcheck( nf90_inq_dimid(ncid, "depth", dimid) )
      call NetCDFcheck( nf90_inquire_dimension(ncid, dimid, len=nz) )  ! resolve nz
      
c      
c     =======  probe horizontal grid descriptors from lon/lat arrays =======
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
c     ======= assess vertical grid descriptor from depth arrays =======
c   
c     ---- set zmode: determine wheter data is sigma or z type by 
c          asking for depth file ----

      if (count_tags(simulation_file, "depth_file")/=0) then  ! initialize as sigma

         zmode = .false.
         write(*,255) "sigma"  
c         
         allocate( sigma(nz) )
         call NetCDFcheck( nf90_inq_varid(ncid, "depth", varid) )
         call NetCDFcheck( nf90_get_var(ncid, varid, sigma)) ! depth interpreted as sigma
         write(*,241) "vertical sigma ", sigma
c
c     --- load depth file ---
c
         allocate( wdepth0(nx,ny) )
         call read_control_data(simulation_file, "depth_file",
     +                          fname)
         write(*,261) trim(adjustl(fname)) 
         call find_free_IO_unit(iunit)
         open(iunit, file=fname, status='old')
         nlines = 0
         wdepth0 = 0.0  ! default dry, if no entries in depth_file
         do 
            read(iunit,*,end=345) ix,iy,wd
            wdepth0(ix,iy) = wd
            nlines = nlines + 1 
         enddo
 345     write(*,265) nlines
         close(iunit)

      else   ! initialize as z mode

         zmode = .true.
         write(*,255) "z"      
         allocate( ccdepth0(nz) )
         call NetCDFcheck( nf90_inq_varid(ncid, "depth", varid) )
         call NetCDFcheck( nf90_get_var(ncid, varid, ccdepth0)) ! depth interpreted as ccdepth0
         write(*,241) "vertical z ", ccdepth0

      endif
      

  229 format(a,a)    
  231 format("load_grid_desc: 3d grid dim (nx,ny,nz) = ", i4,i4,i4)
  232 format("load_grid_desc: lambda1 = ",f12.7,
     +                      " dlambda = ",f12.7," degrees")  
  233 format("load_grid_desc: phi1    = ",f12.7,
     +                      " dphi    = ",f12.7," degrees")  
      
 241  format("load_grid_desc: static ", a, " grid = ", 999f9.3)
 255  format("load_grid_desc: grid type = ", a)
 261  format("load_grid_desc: reading depth data from ", a)
 265  format("load_grid_desc: found ", i, " depth points")
c      
c     ======= assess time grid available (load it as integer) =======
c  
      call NetCDFcheck( nf90_inq_dimid(ncid, "time", dimid) )   ! fixed name
      call NetCDFcheck( nf90_inquire_dimension(ncid, dimid, len=nt) )
      call NetCDFcheck( nf90_inq_varid(ncid, "time", varid) )   ! fixed name
      call NetCDFcheck( nf90_get_var(ncid, varid, ibuf, count=(/2/)) )
      time1 = ibuf(1)
      dtime = ibuf(2)-ibuf(1)   ! regular grid is assumed
      write(*,*) "load_grid_desc: data set contains", nt, "time frames"
c     --- sec_per_time_unit must be set ---
      write(*,*) "load_grid_desc: time resolution = ", 
     +            dtime*sec_per_time_unit, "sec"   

      end subroutine load_grid_desc
      



      subroutine init_topography_z()
c     ------------------------------------------------------------
c     z-mode topography initialization
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
      real,allocatable    :: acc_width0(:)
      character*64        :: temp_name
      real                :: dry_val, landfrac, wetfrac, layerwidth 
c     ------------------------------------------------------------
      call detect_topography(u_var, wetmask, bottom_layer)
      
      landfrac = 1.0 - 1.0*sum(wetmask)/nx/ny
      wetfrac  = 1.0*sum(bottom_layer)/nx/ny/nz
      write(*,*) "init_topography_z: land fraction        = ", landfrac
      write(*,*) "init_topography_z: wet volume fraction  = ", wetfrac
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

      deallocate( acc_width0 )

      end subroutine init_topography_z



      subroutine init_topography_sigma()
c     ------------------------------------------------------------
c     sigma-mode topography initialization
c     wdepth0 has already been loaded  
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
      real,allocatable    :: acc_width0(:)
      character*64        :: temp_name
      real                :: dry_val, landfrac, wetfrac, layerwidth 
c     ------------------------------------------------------------
      call detect_topography(u_var, wetmask, bottom_layer)
    
      landfrac = 1.0 - 1.0*sum(wetmask)/nx/ny
      wetfrac  = 1.0*sum(bottom_layer)/nx/ny/nz
      write(*,*) "init_topography_sigma: land fraction    = ", landfrac
      write(*,*) "init_topography_sigma: wet vol fraction = ", wetfrac
c
c     initialize dynamic 3D topographuc arrays ccdepth, acc_width and wdepth
c     which are static for this PBI
c  
      allocate ( acc_width0(nz+1) )
      acc_width0(1) = 0.0 
      layerwidth    = 2*sigma(1)
      acc_width0(2) = layerwidth    ! scaled depth
      do iz = 2, nz
         layerwidth = 2.0*(sigma(iz)-sigma(iz-1)) - layerwidth  ! scaled width
         acc_width0(iz+1) = acc_width0(iz) + layerwidth         ! scaled depth
      enddo
 
      do ix = 1, nx
         do iy = 1, ny
            ccdepth(ix,iy,:)   = sigma*wdepth0(ix,iy)
            acc_width(ix,iy,:) = acc_width0*wdepth0(ix,iy)
            wdepth(ix,iy)      = wdepth0(ix,iy) 
         enddo
      enddo

      deallocate( acc_width0 )

      end subroutine init_topography_sigma




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
      call attach_variable(u_var, varname, ncid)
      write(*,362) "u", trim(varname) 
c
      call read_control_data(simulation_file, "v_name", varname)
      call attach_variable(v_var, varname, ncid)
      write(*,362) "v", trim(varname)
c     
      if (count_tags(simulation_file, "w_name")/=0) then
         call read_control_data(simulation_file, "w_name", varname)
         call attach_variable(w_var, varname, ncid)
         write(*,362) "w", trim(varname)
         w_defined = .true.
      else
         write(*,363) 
         w_defined = .false.
      endif
c
      call read_control_data(simulation_file, "salinity_name", varname)
      call attach_variable(s_var, varname, ncid_salt)
      write(*,362) "salinity", trim(varname)
c     
      call read_control_data(simulation_file, "temperature_name", 
     +                       varname)
      call attach_variable(t_var, varname, ncid_temp)
      write(*,362) "temperature", trim(varname)
c  
 362  format("open_data_files: varname for ", a," = ", a) 
 363  format("open_data_files: w not available")
 

      end subroutine open_data_files



      subroutine resolve_corresp_frame(aclock, frame)
c     ------------------------------------------------------------------------------------
c     Resolve which time frame (1 <= frame <= nt) to use in data set corresponding to aclock
c     based on time grid descriptors (time1, dtime, nt) resolved at initialization.
c     Use cell centered time association, so nearest time point on grid is applied      
c     In case of exterior time points, the first/last frame will be applied
c     ------------------------------------------------------------------------------------
      type(clock), intent(in) :: aclock
      integer, intent(out)    :: frame
      integer                 :: numsec, it
      real                    :: t, numtics
c     ------------------------------------------------------------------------------------
      call get_period_length_sec(time_offset, aclock, numsec)
      numtics = 1.0*numsec/sec_per_time_unit    ! time distance to offset in data set time unit
      it      = 1 + nint((numtics-time1)/dtime) ! on time grid 
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
c     ------------------------------------------------------------------------------------ 
      write(*,*) "load_data_frame: loading frame ", frame

      call load_xyz_frame_asreal(u_var, u, frame)
      call load_xyz_frame_asreal(v_var, v, frame)

      if (w_defined) then
         call load_xyz_frame_asreal(w_var, w, frame)
      else
          w = 0.  ! default
      endif

      call load_xyz_frame_asreal(t_var, temp, frame)
      call load_xyz_frame_asreal(s_var, salinity, frame)

c     ------------ update state handlers ------------
      data_in_buffers = .true.
      cur_frame       = frame

c     ------------ secondary updates ------------

      layer_width = acc_width(:,:,2:nz+1)-acc_width(:,:,1:nz)
      call update_water_density()
      call update_turbulence(u,v,uswind,vswind,rhow,layer_width,wdepth)   
      
      end subroutine load_data_frame
      

      
      subroutine interpolate_turbulence(geo, r3, status)
c     ------------------------------------------
c     Relay to turbulence module
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
      real                 :: xyz(3)
c     ------------------------------------------ 
      call get_grid_coordinates(geo,xyz(1),xyz(2),xyz(3))
      call ext_turb_intp(xyz, r3)  ! redirected name - defined in module turbulence
      status = 0                   ! no error conditions are flagged
c     ------------------------------------------ 
      end subroutine interpolate_turbulence


      subroutine interpolate_turbulence_deriv(geo, r3, status)
c     ------------------------------------------ 
c     Cartesian derivative (diagonal) of vertical/horizontal turbulence
c     Relay to turbulence module
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
      real                 :: xyz(3), dxyz_dr(3,3), dk(3)
c     ------------------------------------------ 
      call get_grid_coordinates(geo,xyz(1),xyz(2),xyz(3),dxyz_dr)
      call ext_turbderiv_intp(xyz, dk)   ! redirected name - defined in module turbulence
      status = 0                         ! no error conditions are flagged
      r3     = matmul(dxyz_dr, dk)       ! apply chain rule to get Cartesian derivative
c     ------------------------------------------    
      end subroutine interpolate_turbulence_deriv

      


      subroutine dump_land_points(fname)
c     ------------------------------------
c     for debugging
c     ------------------------------------
      character*(*) :: fname
      integer       :: ix,iy
      open(82, file=fname)
      do ix=1,nx
         do iy=1,ny
            if (wetmask(ix,iy)==0) write(82,*) ix,iy
         enddo
      enddo
      close(82)
      end subroutine dump_land_points



      end module

      
