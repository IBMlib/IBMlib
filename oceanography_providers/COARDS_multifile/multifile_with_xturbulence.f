cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  ---------------------------------------------------------------------
c  COARDS like format with time-sliced files and a posteori turbulence
c  ---------------------------------------------------------------------
c
c     Currently
c        * optional zero-equation posteori turbulence schemes
c        * optional Stokes drift for near-surface drifters (need not be time synced with hydrography)
c        * Pick time frame closest to current time, assume time frames is regularly spaced
c        * Grid is either in z-mode or sigma-mode   
c        * Data on same and static vertical z/sigma-grid layer (no time varying sea surface elevation added)
c        * Only physics
c        * Assume all data on same grid and at same times
c        * Only regular lon-lat is allowed horizontally
c        * current, temperature and salinity time frames in same files (not possible to split over fields)
c        * Allow year bootstrapping
c        * As a special case, it should work with all hydrography/stokes frames in a single file
c        * time-to-file mapping should be defined in resolve_data_file (and resolve_stokes_file, if Stokes drift is applied)
c      
c     NOTES: based on added_turbulence.f
c     TODO: add a input handle to allow flipping w (if present) to
c           conform to IBMlib orientation (positive downward)
c           add dynamic sea surface elevation, if present
c           surface wind (uswind,vswind) not loaded       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module physical_fields

      use mesh_grid,
     +     redir_t  => interpolate_turbulence,                    ! shadow by redirection
     +     redir_dt => interpolate_turbulence_deriv,              ! shadow by redirection
     +     interpolate_currents_eulerian => interpolate_currents  ! overlay Stokes drift
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
      use coards_netcdf       ! 3D
      use coards_netcdf_2D,   ! coards_xyt_variable + methods
     +     attach_variable_2d =>  attach_variable             ! same name as coards_netcdf
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

      public :: write_topomask  ! from mesh_grid
      
c      public :: dump_land_points ! non standard tmp
      
c     -------------------- module data --------------------  
      
      integer, parameter         :: FNLEN = 999   ! file/path names
      integer, parameter         :: VNLEN = 256   ! variable names
c
c     ------ data frame handler ------
      
      character(len=FNLEN)       :: hydroDBpath ! hydrographic data sets (proper trailing separator will be added at read, if not provided)
      
      integer, parameter         :: no_bootstrap = -9999
      integer                    :: bootstrap_year ! allow year bootstrapping
      
      character(len=FNLEN), parameter :: not_set_c = "none"     ! assume name conflict unlikely
      integer, parameter              :: not_set_i = -999999999 ! assume number conflict unlikely
                                                  
      integer                    :: cur_frame    ! current frame number in buffer (implicit: not set = no data loaded)
      character(len=FNLEN)       :: cur_file     ! current file name associated with ncid (implicit: not set = no file open)
      integer                    :: ncid         ! NetCDF file handler for currents,temp and salt (this is considered the root data set)  
      logical                    :: w_defined    ! vertical current available else set to zero
      type(coards_xyzt_variable) :: u_var        ! proxy var for easy loading/unpacking; no allocation at instantiation
      type(coards_xyzt_variable) :: v_var        ! proxy var for easy loading/unpacking; no allocation at instantiation
      type(coards_xyzt_variable) :: w_var        ! proxy var for easy loading/unpacking; no allocation at instantiation
      type(coards_xyzt_variable) :: t_var        ! proxy var for easy loading/unpacking; no allocation at instantiation
      type(coards_xyzt_variable) :: s_var        ! proxy var for easy loading/unpacking; no allocation at instantiation
      character(len=VNLEN)       :: u_varname    ! varname to look for in netcdf file
      character(len=VNLEN)       :: v_varname    ! varname to look for in netcdf file
      character(len=VNLEN)       :: w_varname    ! varname to look for in netcdf file    
      character(len=VNLEN)       :: t_varname    ! varname to look for in netcdf file
      character(len=VNLEN)       :: s_varname    ! varname to look for in netcdf file
c
      logical                    :: zmode        ! grid mode toggle: TRUE : z-mode ... FALSE : sigma-mode 
      real, allocatable          :: ccdepth0(:)  ! static vertical grid (z-mode)
      real, allocatable          :: sigma(:)     ! sigma grid           (sigma-mode)
      real, allocatable          :: wdepth0(:,:) ! static CC depths     (sigma-mode)

      real                       :: lambda1, dlambda, phi1, dphi ! local copies for convenience
      real, allocatable          :: uswind(:,:)  ! W:E surface wind [m/s] (currently unused)
      real, allocatable          :: vswind(:,:)  ! S:N surface wind [m/s]  (currently unused)
      real, allocatable          :: layer_width(:,:,:) ! generated from acc_width

c
c     ---------- time grid supported is a regularly spaced sequence described by these three parameters:
c                updated each time a hydrography file is opened
c      
      type(clock)                :: time_offset       ! for netcdf time vector
      integer                    :: time1             ! time since offset of first time frame (in data set time unit) 
      integer                    :: dtime             ! time step between subsequent time frames (in data set time unit) 
      integer                    :: nt                ! number of time frames in data set (possibly only 1)      
      integer                    :: sec_per_time_unit ! seconds per time unit of data set
      
c
c     ---------- provision for optional Stokes drift (2D)
c
      type(coards_xyt_variable) :: ustokes            ! u_component_stokes_drift
      type(coards_xyt_variable) :: vstokes            ! v_component_stokes_drift
      real                      :: stokes_decay_depth ! for exponential extrapolation of Stokes drift from surface
      logical                   :: incl_stokes        ! logical switch to control whether Stokes drift is included
      integer                   :: ncid_stokes        ! data set handlers
      character(len=FNLEN)      :: cur_file_stokes    ! current file name associated with ncid_stokes
      character(len=VNLEN)      :: ust_varname        ! varname to look for in netcdf file
      character(len=VNLEN)      :: vst_varname        ! varname to look for in netcdf file


c     ===================================================
                            contains
c     ===================================================
  
      subroutine init_physical_fields(time)
c     ------------------------------------------
c     Do not trigger data load
c     primitive runtime linux/windows recognition, by allowing user to supply path separator
c     Either time must be provided as argument here or time set prior to calling 
c     init_physical_fields so that data files corresponding to time can be opened
c     and grids can be initialized according to first frames 
c     ------------------------------------------
      type(clock), intent(in),optional :: time
c
      character(len=FNLEN)             :: path, fname
      character                        :: lastchr
      integer                          :: lp, varid
      integer                          :: year,month,day,hour,minute,sec
      integer                          :: sec_of_day, start(99), nwords
      type(clock), pointer             :: aclock
      character(len=VNLEN)             :: cbuf
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

c 
c     --- enable time bootstrapping, if requested ---
c
      if (count_tags(simulation_file, "bootstrap_year")/=0) then
         call read_control_data(simulation_file, 
     +                     "bootstrap_year",bootstrap_year)
         write(*,*) "WARNING: bootstrap_year=",bootstrap_year
      else
         bootstrap_year = no_bootstrap
      endif

c      
c     ======  open data files, but do not load buffers ======
c
      call reset_frame_handlers()
c      
c     ----- resolve which variable names to look for in the netcdffile
c
      call read_control_data(simulation_file, "u_name", u_varname)
      write(*,362) "u", trim(u_varname) 
c
      call read_control_data(simulation_file, "v_name", v_varname)
      write(*,362) "v", trim(v_varname)
c     
      if (count_tags(simulation_file, "w_name")/=0) then
         call read_control_data(simulation_file, "w_name", w_varname)
         write(*,362) "w", trim(w_varname)
         w_defined = .true.
      else
         write(*,363) 
         w_defined = .false.
      endif
c
      call read_control_data(simulation_file, "salinity_name", 
     +                       s_varname)
      write(*,362) "salinity", trim(s_varname)
c     
      call read_control_data(simulation_file, "temperature_name", 
     +                       t_varname)
      write(*,362) "temperature", trim(t_varname)
c  
 362  format("init_physical_fields: varname for ", a," = ", a) 
 363  format("init_physical_fields: w not available")

      aclock => get_master_clock()  ! is set at this point
      call resolve_data_file(aclock, fname)
      call open_data_file(fname)   ! defines ncid

      
c     ----- optional stokes drift: resolve which variable names
c
      if (count_tags(simulation_file, "stokes_drift")/=0) then
         incl_stokes = .true.

         call read_control_data(simulation_file, "stokes_drift", cbuf) ! varname_u, varname_v
         call tokenize(cbuf, start, nwords) ! in string_tools.f
         if (nwords /= 2) then
            write(*,*) "init_physical_fields: "//
     +                 "expect 2 variable names for stokes_drift"
            stop 823
         endif
         read(cbuf(start(1):start(2)-1),'(a)') ust_varname
         read(cbuf(start(2):),          '(a)') vst_varname
      
         call read_control_data(simulation_file, "stokes_decay_depth",
     +     stokes_decay_depth)  !
      
         write(*,318) "u stokes drift",trim(adjustl(ust_varname))
         write(*,318) "v stokes drift",trim(adjustl(vst_varname))
         write(*,*) "stokes decay depth=",stokes_decay_depth,"m"

 318     format("varname for",1x,a,1x,"=",1x,a)

         call resolve_stokes_file(aclock, fname)
         call open_stokes_file(fname)
      else
         incl_stokes = .false.
         write(*,*)  "init_physical_fields: no Stokes drift"
      endif

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
      get_pbi_version="COARDS multifile offline pbi version: $Rev:  $"
      end function



      subroutine close_physical_fields()
c     --------------------------------------------------------- 
c     coards_xyzt_variable instances does not need deallocation
c     ---------------------------------------------------------    
      call close_mesh_grid()
      call close_water_density()
      call reset_frame_handlers()
      call close_horiz_grid_transf()    
      call close_turbulence()
      if (incl_stokes) then
         call close_coards_xyt_variable(ustokes) ! includes grid
         call close_coards_xyt_variable(vstokes) ! includes grid
      endif
c     --- local cleanup ---
      call close_data_file()
      call close_stokes_file()
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
c     Check whether new data fiels should be opened 
c     and whether new frames should be loaded
c     Only load, if needed 
c     ------------------------------------------  
      type(clock), intent(in),optional :: time
      integer, intent(in),optional     :: dt
      
      character(len=FNLEN)             :: filenam
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
      call resolve_data_file(aclock,  filenam)
      if ((cur_file == not_set_c).or.
     +    (cur_file /= filenam)) then
          call open_data_file(filenam)  ! triggers cur_frame = not_set_i
      endif

      call resolve_data_frame(aclock, frame)  ! resolution based on updated time grid
      if ((cur_frame == not_set_i).or.
     +    (cur_frame /= frame)) then
          call load_data_frame(frame)   ! updates cur_frame
      endif
         
c      
c     ---- stokes data does not necessary use same time grid
c
      if (incl_stokes) then
         aclock => get_master_clock()
         call resolve_stokes_file(aclock, filenam)
         if ((cur_file_stokes == not_set_c).or.
     +       (cur_file_stokes /= filenam)) then
             call open_stokes_file(filenam)   ! cause variable reattachment
         endif
         call sync_to_time(ustokes, aclock)   ! frame load externalized
         call sync_to_time(vstokes, aclock)   ! frame load externalized
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
      character(len=VNLEN) :: lon_name, lat_name ! assume long enough
      character(len=FNLEN) :: fname ! assume long enough
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
      integer              :: dimid,varid,idum, gnci
      integer              :: ix,iy,iz,count4D(4),istat
      real,allocatable     :: acc_width0(:)
      real                 :: dry_val, landfrac, wetfrac, layerwidth 
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
      integer              :: dimid,varid,idum, gnci
      integer              :: ix,iy,iz,count4D(4),istat
      real,allocatable     :: acc_width0(:)
      real                 :: dry_val, landfrac, wetfrac, layerwidth 
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




      subroutine reset_frame_handlers()
c     ------------------------------------------
c     Reset both hydrography / stokes drift data handlers
c     (even if stokes drift is not enabled)
c     ------------------------------------------
      cur_file         = not_set_c
      cur_frame        = not_set_i
      ncid             = not_set_i
      cur_file_stokes  = not_set_c
      ncid_stokes      = not_set_i
      end subroutine reset_frame_handlers



      subroutine close_data_file()
c     ------------------------------------------
c     Only close, if data files are open
c     ------------------------------------------
      if (cur_file == not_set_c) return

      call NetCDFcheck( nf90_close(ncid) )
      ncid     = not_set_i  
      cur_file = not_set_c
      write(*,*) "close_data_file: closed open NetCDF set"
      
      end subroutine close_data_file

      
      subroutine close_stokes_file()
c     ------------------------------------------
c     Only close, if data files are open
c     ------------------------------------------
      if (cur_file_stokes == not_set_c) return

      call NetCDFcheck( nf90_close(ncid_stokes) )
      ncid_stokes     = not_set_i  
      cur_file_stokes = not_set_c
      write(*,*) "close_stokes_file: closed open NetCDF set"
      
      end subroutine close_stokes_file



      subroutine open_data_file(fname)
c     ---------------------------------------------------
c     Open hydrography_file (currents, temperature and salinity) 
c     and attach variables to entries in netcdf file     
c     ---------------------------------------------------
      character(len=FNLEN),intent(in) :: fname     ! assume long enough
      character(len=FNLEN)            :: fullfname ! 
      integer                         :: year,month,day,hour,minute,sec
      integer                         :: sec_of_day, dimid, varid 
      integer                         :: ibuf(2)
c     ------------------------------------------
      call close_data_file() ! in case they are open ...
      cur_file  = fname      ! just local file name
      cur_frame = not_set_i  ! invalidate buffers

      fullfname = trim(hydroDBpath) // trim(adjustl(fname)) ! incl separator
      call NetCDFcheck( nf90_open(fullfname, NF90_NOWRITE, ncid))
      write(*,352) trim(fullfname) 
             
 352  format("open_data_files: opened hydrography file: ", a)

      call attach_variable(u_var, u_varname, ncid)
      call attach_variable(v_var, v_varname, ncid)
      if (w_defined) call attach_variable(w_var, w_varname, ncid)
      call attach_variable(s_var, s_varname, ncid)
      call attach_variable(t_var, t_varname, ncid)
c      
c     ======= assess time grid available (load it as integer) =======
c        
      call resolve_time_parameters(ncid, "time", year, month, day, 
     +                           hour, minute, sec, sec_per_time_unit)
      write(*,*) "open_data_files: sec_per_time_unit=", 
     +            sec_per_time_unit
      sec_of_day = 3600*hour + 60*minute + sec
      call set_clock(time_offset, year, month, day, sec_of_day) 
      write(*,*) "open_data_files: data time offset is"
      call write_clock(time_offset)

      call NetCDFcheck( nf90_inq_dimid(ncid, "time", dimid) )   ! fixed name
      call NetCDFcheck( nf90_inquire_dimension(ncid, dimid, len=nt) )
      call NetCDFcheck( nf90_inq_varid(ncid, "time", varid) ) ! fixed name
c     --- set time1 and time1, depending on nt
      if (nt>1) then
         call NetCDFcheck(nf90_get_var(ncid,varid,ibuf,count=(/2/)))
         time1 = ibuf(1)
         dtime = ibuf(2)-ibuf(1) ! regular grid is assumed
         write(*,*) "open_data_files: data set contains", nt, 
     +        "time frames"
c     --- sec_per_time_unit must be set ---
         write(*,*) "open_data_files: time resolution = ", 
     +        dtime*sec_per_time_unit, "sec" 
      else                     ! nt==1 : file only contains a single frame
         call NetCDFcheck( nf90_get_var(ncid, varid, time1) )
         dtime = 0 ! render defined 
         write(*,*) "open_data_files: data set contains 1 time frame"
      endif

      end subroutine open_data_file


      subroutine open_stokes_file(fname)
c     ---------------------------------------------------
c     Open file with 2D stokes drift fields
c     and attach variables to entries in netcdf file
c     This subroutie is only called, if Stokes drift is enabled
c     ---------------------------------------------------
      character(len=FNLEN),intent(in) :: fname     ! assume long enough
      character(len=FNLEN)            :: fullfname !    
      integer                         :: status
c     ---------------------------------------------------
      call close_stokes_file() ! in case they are open ..
      cur_file_stokes = fname  ! just local file name
c
c     ---------- optional Stokes drift ----------
c
      fullfname = trim(hydroDBpath)//trim(adjustl(fname))
      status = nf90_open(fullfname, nf90_nowrite, ncid_stokes)
      if (status /= nf90_noerr) then
         write(*,*) "error opening stokes_drift set",
     +                 trim(adjustl(fullfname))
         stop 22
      endif
      call attach_variable_2d(ustokes, ust_varname,
     +     ncid_stokes,"nearest_valid") 
      call attach_variable_2d(vstokes, vst_varname,
     +     ncid_stokes,"nearest_valid")
     
      write(*,*) "open_stokes_file: opened ", trim(fullfname)

      end subroutine open_stokes_file



      subroutine resolve_data_file(aclock, fname)
c     ----------------------------------------------------
c     Name template to be adapted 
c     Currently map to hydro_YYYY_MM.nc where data covers 
c     frames in month,year = MM,YYYY
c     ----------------------------------------------------
      type(clock), intent(in)           :: aclock
      character(len=FNLEN), intent(out) :: fname
      integer                           :: it,  iye, imo, idy
c     ----------------------------------------------------
      call get_date_from_clock(aclock, iye, imo, idy)
      write(fname, 311) iye, imo 
 311  format("hydro_", i4.4, "_", i2.2, ".nc")
      end subroutine resolve_data_file



      subroutine resolve_stokes_file(aclock, fname)
c     ----------------------------------------------------
c     Name template to be adapted
c     Currently map to stokes_YYYY_MM.nc where data covers 
c     frames in month,year = MM,YYYY
c     Only called if Stokes drift is enabled
c     ----------------------------------------------------
      type(clock), intent(in)           :: aclock
      character(len=FNLEN), intent(out) :: fname
      integer                           :: it,  iye, imo, idy
c     ----------------------------------------------------
      call get_date_from_clock(aclock, iye, imo, idy)
      write(fname, 312) iye, imo 
 312  format("stokes_", i4.4, "_", i2.2, ".nc")
      end subroutine resolve_stokes_file



      subroutine resolve_data_frame(aclock, frame)
c     ------------------------------------------------------------------------------------
c     Resolve which time frame (1 <= frame <= nt) to use in data set corresponding to aclock
c     based on current time grid descriptors (time1, dtime, nt) resolved at open_data_file
c     Use cell centered time association, so nearest time point on grid is applied      
c     In case of exterior time points, the first/last frame will be applied
c     Special case with one time frame in data set is handled
c
c     08 Nov 2021: replaced get_period_length_sec -> get_period_length_hour
c                  since we are hitting roof on period length > 69 years measured in seconds
c     16 Nov 2021: enable year bootstrapping, if set
c     ------------------------------------------------------------------------------------
      type(clock), intent(in) :: aclock
      type(clock)             :: clock_bootstrapped
      integer, intent(out)    :: frame
      integer                 :: it,  iye, imo, idy, isec
      real                    :: numhours, numtics
c     ------------------------------------------------------------------------------------
      if (bootstrap_year /= no_bootstrap) then
         call get_date_from_clock(aclock, iye, imo, idy)
         call get_second_in_day(aclock, isec)
         call set_clock(clock_bootstrapped, 
     +                  bootstrap_year, imo, idy, isec)
         call get_period_length_hour(time_offset, 
     +                  clock_bootstrapped, numhours)               ! apply bootstrapped clock
      else
         call get_period_length_hour(time_offset, aclock, numhours) ! no bootstrapping
      endif


      if (nt>1) then  ! multi frame data file
         numtics = 3600.0*numhours/sec_per_time_unit ! time distance to offset in data set time unit
         it      = 1 + nint((numtics-time1)/dtime) ! on time grid 
         if ((it < 1) .or. (it > nt)) then
            write(*,286) it
         endif
 286     format("resolve_data_frame::warning: time extrapolation (it=",
     +        i,")")
         frame = min(nt, max(it,1))
      else
         frame = 1    ! data file only contains a single frame
      endif


      end subroutine resolve_data_frame


      


      subroutine load_data_frame(frame)
c     ------------------------------------------------------------------------------------
c     Verbose reading af NetCDF sets to grid buffers from multi frame data set
c     
c     It is assumed that grids do not change during a simulation (assertion not checked)
c     This version does not feature dynamic sea surface elevation, so vertical auxillary
c     arrays are static.
c     Variable handler are retrieved at open_data_files
c     Trigger secondary update of managend properties dependent on primary hydrography (u,v,w,s,t)
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

      

      subroutine interpolate_currents(xyz, uvw, status)      
c     -------------------------------------------------------------- 
c     Add optional stokes drift to Eulerian currents
c     -------------------------------------------------------------- 
      real, intent(in)     :: xyz(:)  ! (lon,lat,depth)
      real, intent(out)    :: uvw(:)
      integer, intent(out) :: status

      integer              :: statu, statv, statw, ist1,ist2,ist
      real                 :: x,y,z, sx,sy,sz, uv(3)
c     -------------------------------------------------------------- 
      call interpolate_currents_eulerian(xyz, uvw, status) 
c
      if (incl_stokes) then
         call interpolate_xyt_variable(xyz,ustokes,0,uv(1),ist1)
         call interpolate_xyt_variable(xyz,vstokes,0,uv(2),ist2)
         ist   = max(ist1,ist2)
         uv(3) = 0                     ! no vertical component
         if (ist > 0) uv = 0.0         ! silently ignore problems
         uvw = uvw + uv * exp(-xyz(3)/stokes_decay_depth) ! vector operation
      endif
c     ------------------------------------------ 
      end subroutine interpolate_currents



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

      
