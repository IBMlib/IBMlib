
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     RCO-SCOBI format physics API implementation 
c     ---------------------------------------------------
c
c     Based on native RCO-SCOBI snap format (big_endian must be set)
c     snap files have one time frame per file; data is stored as
c     horizontal layers, where only wet points are stored
c     snap files consist of a header section and the a sequence of horizontal layers ("snaps")
c     content + order of snaps seem variable - this implementation does not assert a
c     specific ordering
c
c     API implementation based on pseudo fortran code fragments provided
c     by SMHI (see format_snap.txt)
c  
c     TODO:
c       interpolate_currents: w + RCO flux conservation scheme
c       convert potential temperature (provided) to physical temperature
c       load biogeochem
c       load diffusivity  
c 
c     NOTES: 
c     
c     LOG: 
c          
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module physical_fields

      use mesh_grid, 
     +          interpolate_currents_cc => interpolate_currents  ! use face-centered local

      use horizontal_grid_transformations
      use geometry
      use array_tools
      use run_context, only: simulation_file
      use time_services           ! import clock type and time handling
      use input_parser

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

c      public :: interpolate_wind    ! currently unused
c     public :: interpolate_zooplankton
      
c     ---- currently not implemented ----
      
      public :: interpolate_zooplankton
      public :: interpolate_oxygen
      public :: interpolate_nh4                       
      public :: interpolate_no3                       
      public :: interpolate_po4                        
      public :: interpolate_diatoms                   
      public :: interpolate_flagellates                 
      public :: interpolate_cyanobacteria              
      public :: interpolate_organic_detritus                   
      public :: interpolate_part_org_matter 
      public :: interpolate_DIC 
      public :: interpolate_alkalinity  
      public :: interpolate_DIN    
      public :: interpolate_chlorophyl  
c      

      

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
      
      integer, parameter         :: tag_lenght  = 10   ! 
      character*(*), parameter   :: not_set_c = "none" ! assume name conflict unlikely
      integer, parameter         :: not_set_i = -9999
     
      logical                    :: data_in_buffers     = .false.    ! first time setting
      logical                    :: data_files_are_open = .false.    ! first time setting
      character(len=tag_lenght)  :: cur_tag             = not_set_c  ! tag (YYYYMMDDHH) of current open set
      integer                    :: iunit_data          = -1         ! logical unit for data I/O    
      
      real, parameter            :: molecular_diffusivity = 1.e-9    ! unit m2/s
      real                       :: constant_horiz_diffusivity       ! unit m2/s; negative => load data from data set
      real                       :: constant_vertic_diffusivity      ! unit m2/s; negative => load data from data set
c     --- 3D grids ---
      real, allocatable          :: zoobuf(:,:,:)  ! buffer for adding micro+meso zooplankton 
      real,allocatable,target    :: ccdepth0(:,:,:)   ! reference cell center depth water below surface [m] (dslm=0)
      real,allocatable,target    :: acc_width0(:,:,:) ! accumulated water above this undisturbed layer [m] dim=nz+
c     --- 2D grids ---    
      real, allocatable          :: wdepth0(:,:)   ! undisturbed water depth, ncvariable h(x,y)
      integer, allocatable       :: kmu(:,:)       ! for unpacking u,v
c     --- 1D grids ---
  
c     --- units    ---    
 
      include "CERES_setup.h"
c     ===================================================
                            contains
c     ===================================================
  
      subroutine init_physical_fields(time)
c     ------------------------------------------
c     Do not trigger data load
c     primitive runtime linux/windows recognition, by allowing user to supply path separator 
c     ------------------------------------------
      type(clock), intent(in),optional :: time
      character(len=len(hydroDBpath))  :: sbuf
      character                        :: lastchr
      integer                          :: lp
      character(len=tag_lenght)        :: tag
      type(clock), pointer             :: aclock
c     ------------------------------------------
      if (present(time)) then
         call set_master_clock(time)
         write(*,*) "set simulation PBI clock = "
         call write_clock(time)
      endif
      write(*,*) trim(get_pbi_version()) 

      call read_control_data(simulation_file,"hydroDBpath", sbuf)

c     --- provide primitive runtime linux/windows recognition, by allowing user to 
c         supply path separator to hydroDBpath

      lp      = len(sbuf)
      sbuf    = adjustr(sbuf)
      lastchr = sbuf(lp:lp)
      if     (lastchr == "/") then                  ! found explicit linux seperator
         hydroDBpath = trim(adjustl(sbuf)) ! keep seperator as provided
      elseif (lastchr == "\") then                  ! found explicit windows seperator
         hydroDBpath = trim(adjustl(sbuf))          ! keep seperator as provided
      else
         hydroDBpath = trim(adjustl(sbuf)) // "/"   ! no separator, assume linux and prepend linux delimiter
      endif

      write(*,*) "init_physical_fields: hydrographic database path =", 
     +           trim(adjustl(hydroDBpath)) 

      if (present(time)) then
         call resolve_tag_from_time(time, tag)
      else
         aclock => get_master_clock()
         call resolve_tag_from_time(aclock, tag)
      endif

      call reset_frame_handler()
      call open_data_files(tag)
      call load_grid()  ! invokes init_horiz_grid_transf
      call init_mesh_grid(init_biogeochem = .false.)  ! allocate core arrays, currently no biogeochem
      call load_topography() ! requires mesh_grid arrays are allocated

      allocate( zoobuf(nx, ny, nz) )
     
c
c     ---- currently no horizontal turbulent diffusivity in data set ----
c          initialize it to molecular diffusivity lower limit, if no values 
c          are provided
c
c     resolved user provided constant horizontal_diffusivity (optional)  
c    
      if (count_tags(simulation_file, "horizontal_diffusivity")/=0) then
         call read_control_data(simulation_file,
     +                          "horizontal_diffusivity",
     +                          constant_horiz_diffusivity)
         
         write(*,563) constant_horiz_diffusivity
         hdiffus = constant_horiz_diffusivity
      else
         constant_horiz_diffusivity = -999.0  ! set handle to apply dynamic horiz_diffusivity, if available
         write(*,564) molecular_diffusivity
         hdiffus = molecular_diffusivity
      endif

  563 format("init_physical_fields: using constant horiz_diffusivity =",
     +       e12.5," m2/s")
  564 format("init_physical_fields: init default horiz_diffusivity   =",
     +       e12.5," m2/s")


c
c     ---- currently no vertical turbulent diffusivity in data set ----
c          initialize it to molecular diffusivity lower limit, if no values 
c          are provided
c
c     resolved user provided constant vertical_diffusivity (optional)  
c    
      if (count_tags(simulation_file, "vertical_diffusivity")/=0) then
         call read_control_data(simulation_file,
     +                          "vertical_diffusivity",
     +                          constant_vertic_diffusivity)
         
         write(*,573) constant_vertic_diffusivity
         vdiffus = constant_vertic_diffusivity
      else
         constant_vertic_diffusivity = -999.0  ! set handle to apply dynamic vertic_diffusivity, if available
         write(*,574) molecular_diffusivity
         vdiffus = molecular_diffusivity
      endif

  573 format("init_physical_fields: using constant vertic_diffusivity =",
     +       e12.5," m2/s")
  574 format("init_physical_fields: init default vertic_diffusivity   =",
     +       e12.5," m2/s")


c     ---- specials for this data set / version ----

      u_wind_stress = 0.0  ! no wind in this data set
      v_wind_stress = 0.0  ! no wind in this data set
      
      end subroutine init_physical_fields

 

      character*100 function get_pbi_version()  
      get_pbi_version="RCO-SCOBI offline pbi version: $Rev:  $"
      end function



      subroutine close_physical_fields()
c     ------------------------------------------ 
c     ------------------------------------------    
      call close_mesh_grid()
      call reset_frame_handler()
      call close_horiz_grid_transf()
c     --- local cleanup ---
      deallocate( zoobuf   )
      deallocate( wdepth0  )
      deallocate( kmu      )
c     ------------------------------------------
      end subroutine 


      subroutine update_physical_fields(time, dt)
c     ------------------------------------------  
      type(clock), intent(in),optional :: time
      integer, intent(in),optional     :: dt
      logical                          :: update
      character(len=tag_lenght)        :: tag
      type(clock), pointer             :: aclock
c     ------------------------------------------  
      aclock => get_master_clock()
      if (present(time)) then
         call set_master_clock(time)
      elseif (present(dt)) then
         call add_seconds_to_clock(aclock, dt)
      endif
c
      call resolve_tag_from_time(aclock, tag)
      call update_dataset(tag)
c     ------------------------------------------       
      end subroutine update_physical_fields

      

      subroutine resolve_tag_from_time(aclock, tag)
c     ------------------------------------------
c     Generate tag YYYYMMDDHH corresponding to aclock
c     SMHI data is provided as hourly frames representing
c     time stamp in file (mirrored in filename)
c
c     Add half hour to aclock and truncate corresponding time to
c     tag YYYYMMDDHH   (where HH = 00 ... 23)
c     File name format:
c     filename ~ d2004123109.snap1
c     tag      =  YYYYMMDDHH
c     ------------------------------------------      
      type(clock), intent(in)                :: aclock
      character(len=tag_lenght), intent(out) :: tag
c   
      type(clock)                            :: dummy_clock 
      integer                                :: isec,year,month,day,hour
c     ----------------------------------------------------------------------------+      
      call set_clock(dummy_clock, aclock)
      call add_seconds_to_clock(dummy_clock, 1800)  ! allow rounding by truncation
      call get_date_from_clock(dummy_clock, year, month, day)
      call get_second_in_day(dummy_clock, isec)
      hour = int(isec/3600) ! cast by truncation : 0 .. 23
      write(tag,433) year, month, day, hour
 433  format(i4.4, 3i2.2)
      end subroutine resolve_tag_from_time


      
      subroutine update_dataset(tag)
c     ------------------------------------------
c     Only open new data set file if necessary
c     In this setup there is only one frame in each file
c     ------------------------------------------
      character(len=tag_lenght), intent(inout) :: tag 
c     ------------------------------------------
      if ((.not.data_in_buffers).or.(tag /= cur_tag)) then
         call reset_frame_handler() ! close old
         call open_data_files(tag)
         call load_data_frame()   ! corresponding to tag
      endif
       
      end subroutine update_dataset




      
      subroutine load_grid()
c     -------------------------------------------------------
c     Should only be invoked at start
c     -------------------------------------------------------
      real*8  :: dx, dy, stlon, stlat, dxdeg, dydeg
      real*4  :: ird, ird2, ird3, ird4
      real    :: lambda1, dlambda, phi1, dphi ! LOCAL DUMMIES
      real*8  :: ird0,ird20,ird30,ird40
c     -------------------------------------------------------
      if (data_in_buffers) stop "load_grid:unexpected"
      
      write(*,*) "load_grid: loading grid"
      
c     --- fetch grid dimensions / spacings
      rewind(iunit_data)
      read(iunit_data) ird
      read(iunit_data) ird
      nz = ird                 ! Number of levels for 3D-fields - unsafe cast
      read(iunit_data) ird
      read(iunit_data) ird
      nx = ird                 ! Number of points along a parallel (x or i direction) - unsafe cast
      read(iunit_data) ird
      ny = ird                 ! Number of points along a meridian (y or j direction) - unsafe cast
      read(iunit_data) ird      
      read(iunit_data) ird      
      read(iunit_data) ird      
      read(iunit_data) ird             
      read(iunit_data) ird      
      read(iunit_data) ird     
      read(iunit_data) ird
      read(iunit_data) ird
      read(iunit_data) ird0,ird20,ird30

c     dxdeg   Step in longitude between data points in degrees = dlambda
c     dydeg   Step in latitude between data points in degrees  = dphi
      read(iunit_data) dx, dy, dxdeg, dydeg
      dlambda = dxdeg  ! real8->real
      dphi    = dydeg  ! real8->real
      
c     stlon   Longitude of RCO reference point
c     stlat   Latitude of RCO reference point
      read(iunit_data) stlon, stlat
      
      lambda1 = stlon + 0.5*dlambda
      phi1    = stlat + 0.5*dphi
      
      write(*,231) nx, ny, nz
      write(*,232) lambda1, dlambda
      write(*,233) phi1,    dphi 
    
      write(*,*) "load_grid: horizontal initialization"
      call init_horiz_grid_transf(lambda1, phi1, dlambda, dphi)

  231 format("load_grid: 3d grid dim (nx,ny,nz) = ", i4,i4,i4)
  232 format("load_grid: lambda1 = ",f12.7,
     +                 " dlambda = ",f12.7," degrees")  
  233 format("load_grid: phi1    = ",f12.7,
     +                 " dphi    = ",f12.7," degrees")  
      
      end subroutine load_grid
      


      subroutine load_topography()
c     ------------------------------------------------------------
c     Remaining grid/topography initialization (that requires
c     mesh_grid arrays are allocated by init_mesh_grid)  
c     Set mesh_grid arrays:
c        wetmask     
c        bottom_layer
c     ------------------------------------------------------------
      real*8  :: dx, dy, stlon, stlat, dxdeg, dydeg
      real*4  :: ird, ird2, ird3, ird4
      real*8  :: ird0,ird20,ird30,ird40

      integer              :: i,j,iz,nsnaps
      real*4, allocatable  :: rd1d(:,:), rd2d(:,:)
c     ------------------------------------------------------------
      if (data_in_buffers) stop "oad_topography:unexpected"
      if (nz /= nz_expect) then
         write(*,*) "load_topography: expected nz=", nz_expect
         write(*,*) "load_topography: found    nz=", nz
         stop "load_topography:unexpected"
      endif
      
      allocate( wdepth0(nx,ny) ) ! undisturbed water depth
      allocate( kmu(nx,ny) )     ! auxillary for unpacking u,v
      allocate( ccdepth0(nx,ny,nz)) 
      allocate( acc_width0(nx,ny,nz+1)  )
      
      write(*,*) "load_topography: loading topography"
      
      rewind(iunit_data)
      read(iunit_data) ird
      read(iunit_data) ird    
      read(iunit_data) ird
      read(iunit_data) ird      
      read(iunit_data) ird      
      read(iunit_data) ird      
      read(iunit_data) ird
      nsnaps = ird   ! Number of fields in this file where 3d-fields are counted as the number of levels they contain.
      read(iunit_data) ird      
      read(iunit_data) ird             
      read(iunit_data) ird      
      read(iunit_data) ird     
      read(iunit_data) ird
      read(iunit_data) ird
      read(iunit_data) ird0,ird20,ird30
      read(iunit_data) dx, dy, dxdeg, dydeg
      read(iunit_data) stlon, stlat
      
      allocate(rd1d(nsnaps,2))
      read(iunit_data) rd1d(1:nsnaps,1), rd1d(1:nsnaps,2)

      allocate(rd2d(nx,ny))
      read(iunit_data) rd2d(1:nx,1:ny)
      wetmask      = 0      ! default dry
      wdepth0      = -1.0   ! default dry
      bottom_layer =  0         ! default dry
      do iz=1,nz
         ccdepth0(:,:,iz)   = (iz-0.5)*layer_width
         acc_width0(:,:,iz) = (iz-1.0)*layer_width
      enddo
      acc_width0(:,:,nz+1) = nz*layer_width
      
      do i=1,nx
         do j=1,ny
            bottom_layer(i,j) = rd2d(i,j) ! Number of levels in each grid cell on t-grid.
            if (bottom_layer(i,j)>0) then
                wetmask(i,j) = 1    ! wet
                wdepth0(i,j) = acc_width0(1,1,bottom_layer(i,j)+1) 
            endif
         enddo
      enddo
c     .... initialize auxillary grid for unpacking u,v (kmt => bottom_layer)
      kmu = 0
      do i=1,nx-1
        do j=1,ny-1
           kmu(i,j) = min(bottom_layer(i,j),  bottom_layer(i+1,j),
     +                    bottom_layer(i,j+1),bottom_layer(i+1,j+1))
        enddo
      enddo


      ccdepth   = 0.0
      acc_width = 0.0

      deallocate(rd1d)
      deallocate(rd2d)
      
      end subroutine load_topography




      subroutine reset_frame_handler()
c     ------------------------------------------
c     no buffer deallocation
c     ------------------------------------------
      call close_data_files()
      data_in_buffers = .false.
      cur_tag         = not_set_c
      end subroutine reset_frame_handler



      subroutine close_data_files()
c     ------------------------------------------
c     Only close, if data files are open
c     ------------------------------------------
      if (data_files_are_open) then
         close(iunit_data)       
         write(*,*) "close_data_files: closed unit ", iunit_data
         iunit_data = -1
         data_files_are_open = .false.
      endif
      end subroutine close_data_files




      subroutine open_data_files(tag)
c     ---------------------------------------------------
c     Assemble file name from path and tag and open data set
c     file name template is in format 333
c     path: hydroDBpath /  d  + TAG  + .snap1
c     hydroDBpath includes trailing separator
c     ---------------------------------------------------
      character(len=tag_lenght),intent(in) :: tag
      character*999                        :: fname ! assume long enough
c     ------------------------------------------
      call close_data_files() ! in case they are open ...
      write(fname,333) trim(adjustl(hydroDBpath)), tag
 333  format(a,"d",a,".snap1")
c     --- we don't need to check whether iunit_data is still available
c     if it was used previously
      
      if (iunit_data<0) then                  ! starting condition
         call find_free_IO_unit(iunit_data)   ! runtime_tools
      endif
      open(unit=iunit_data, file=fname, status='old',  
     +     form='unformatted', action='read')
      
      write(*,*) "open_data_files: opened data set ", 
     +           trim(adjustl(fname))
c  
      data_files_are_open = .true.
      cur_tag = tag
      data_in_buffers     = .false. ! needs to be reloaded
      
      end subroutine open_data_files


      subroutine skip_header(nsnaps, nlen)
c     ------------------------------------------------------------------------------------
c     Advance open data file to point where level list is
c     print time stamp and return number of snaps (nsnaps)
c     and max buffer length (nlen)      
c     ------------------------------------------------------------------------------------
      integer, intent(out) :: nsnaps, nlen

      integer :: itt, km, nt, imt, jmt,i,j,iv,vlen
      integer :: year, month, day, hour, minute, second
      real*8  :: dx, dy, stlon, stlat, dtts, snapd, totsec
      real*8  :: dxdeg, dydeg
      real*4  :: ird, ird2, ird3, ird4
      real*8  :: ird0,ird20,ird30,ird40
      
      rewind(iunit_data)
      read(iunit_data) ird
      read(iunit_data) ird
      read(iunit_data) ird      
      read(iunit_data) ird      
      read(iunit_data) ird      
      read(iunit_data) ird
      NLEN = ird                ! Maximum number of active points in any field/level.            
      read(iunit_data) ird
      NSNAPS = ird              ! Number of fields in this file where 3d-fields are counted as       
      read(iunit_data) ird
      year = ird                ! year, month, day, hour, minute and second for which the fields
      read(iunit_data) ird              ! in this file are valid.
      month = ird
      read(iunit_data) ird
      day = ird
      read(iunit_data) ird
      hour = ird
      read(iunit_data) ird
      minute = ird
      read(iunit_data) ird
      second = ird

      read(iunit_data) ird0,ird20,ird30
      read(iunit_data) dx, dy, dxdeg, dydeg
      read(iunit_data) stlon, stlat

      write(*,411) year, month, day, hour, minute, second
      write(*,412) NSNAPS

 411  format("skip_header: reading header for time frame",i4,5(":",i2))
 412  format("skip_header: time frame contains",1x,i8,1x,"data layers")     
      end subroutine skip_header


      
      subroutine load_data_frame()
c     ------------------------------------------------------------------------------------
c     Verbose reading af current open data file
c         
c     ------------------------------------------------------------------------------------    
      integer :: nsnaps,nlen
      integer :: i,iv,vlen,ix,iy,iz
      real*4  :: ird
      integer, allocatable :: ispvar(:), isplev(:) 
      real*4, allocatable  :: rd1d(:,:), rd2d(:,:), snap1d(:)
      integer  :: check_dslm
      integer  :: check_u(nz_expect), check_v(nz_expect)
      integer  :: check_temp(nz_expect), check_salinity(nz_expect)
      real     :: dz
c     ------------------------------------------------------------------------------------
      call skip_header(nsnaps,nlen) ! after rewind
      write(*,*) "load_data_frame: loading frame ", cur_tag

c     ---- read table of contents ----
      allocate(rd1d(nsnaps,2))
      allocate(ispvar(nsnaps))
      allocate(isplev(nsnaps))
      read(iunit_data) rd1d(1:nsnaps,1), rd1d(1:nsnaps,2)
      do i=1,NSNAPS
         ispvar(i) = rd1d(i,1)  ! Array of parameter numbers.
         isplev(i) = rd1d(i,2)  ! Array of levels.
      enddo
      allocate(rd2d(nx,ny))     
      read(iunit_data) rd2d(1:nx,1:ny) ! assume topography unchanged
      allocate(snap1d(nlen))

      u        = 0
      v        = 0
      w        = 0
      temp     = 0
      salinity = 0
      dslm     = 0
c      
      check_dslm     = 0
      check_u        = 0
      check_v        = 0
      check_temp     = 0
      check_salinity = 0
c      
      !     Loop over the fields or until we found everything we want.
      do iv = 1, nsnaps
         ! Reading parameter ispvar(iv), level isplev(iv)
         read(iunit_data) ird
         vlen = ird             ! Number of active points in this field.
         if (vlen.gt.0) then
            read(iunit_data) snap1d(1:vlen)
c     write(*,332) ispvar(iv), isplev(iv), vlen
         else
            cycle 
         endif
         iz = isplev(iv)        ! 1 for 2D layers
         if ((iz<1).or.(iz>nz_expect)) then
            write(*,*) "load_data_frame: illegal iz=", iz
            stop 832
         endif
         if (ispvar(iv) == levnum_dslm) then
            if (iz /= 1) then
              write(*,*) "load_data_frame: unexpected iz=", iz
              stop 832
            endif
            call unpack_snap_data(snap1d(1:vlen),iz,bottom_layer,dslm)
            check_dslm = 1
         elseif (ispvar(iv) == levnum_t) then
            call unpack_snap_data(snap1d(1:vlen),iz,bottom_layer,
     +                            temp(:,:,iz))
            check_temp(iz) = 1
         elseif (ispvar(iv) == levnum_s) then
            call unpack_snap_data(snap1d(1:vlen),iz,bottom_layer,
     +                            salinity(:,:,iz))
            check_salinity(iz) = 1
         elseif (ispvar(iv) == levnum_u) then
            call unpack_snap_data(snap1d(1:vlen),iz,kmu,u(:,:,iz))
            check_u(iz) = 1
         elseif (ispvar(iv) == levnum_v) then
            call unpack_snap_data(snap1d(1:vlen),iz,kmu,v(:,:,iz))   
            check_v(iz) = 1
         endif
         
      enddo
c
      if (check_dslm /= 1) then
         write(*,*) "load_data_frame: dslm not found"
         stop 122
      endif
      if (sum(check_temp) /= nz_expect) then
         write(*,*) "load_data_frame: temp frames missing"
         write(*,*) "load_data_frame: check_temp = ", check_temp
         stop 123
      endif
      if (sum(check_salinity) /= nz_expect) then
         write(*,*) "load_data_frame: salinity frames missing"
         write(*,*) "load_data_frame: check_salinity = ", check_salinity
         stop 124
      endif
      if (sum(check_u) /= nz_expect) then
         write(*,*) "load_data_frame: u frames missing"
         write(*,*) "load_data_frame: check_u = ", check_u
         stop 125
      endif
      if (sum(check_v) /= nz_expect) then
         write(*,*) "load_data_frame: v frames missing"
         write(*,*) "load_data_frame: check_v = ", check_v
         stop 126
      endif
      
c     ------ generate auxillary wet point descriptors
c            ccdepth/acc_width already initialized
c 
      wdepth = 0.0 
      do ix=1,nx
         do iy=1,ny
            if (wetmask(ix,iy)==0) cycle ! nothing for dry points
            dz                  = dslm(ix,iy)
            wdepth(ix,iy)       = max(0.0, wdepth0(ix,iy) + dz)   
            acc_width(ix,iy,2:) = acc_width0(ix,iy,2:) + dz       
            ccdepth(ix,iy,1)    = ccdepth0(ix,iy,1)    + dz*0.5   
            ccdepth(ix,iy,2:)   = ccdepth0(ix,iy,2:)   + dz       
         enddo
      enddo
      
c
c     ---- enforce exact consistency between acc_width(:,:,nz+1) and wdepth 
c          (there are numerical diffs of order 10^-5)
c
      acc_width(:,:,nz+1) = wdepth 


c     ------------ update state handlers ------------
      
      data_in_buffers = .true.

c     ------------ clean up ------------
      
      deallocate( rd1d )
      deallocate( ispvar )
      deallocate( isplev )
      deallocate( rd2d )
      deallocate( snap1d )

      
c     ------------ misc debug ------------
c$$$      do ix=1,nx
c$$$         do iy=1,ny
c$$$            if (wetmask(ix,iy)==0) 
c$$$     +          write(77,*) (ix-1)*dlambda + lambda1,
c$$$     +                     (iy-1)*dphi    + phi1
c$$$            if (wetmask(ix,iy)==0) write(77,*) ix,iy
c$$$         enddo
c$$$      enddo
c$$$      stop 84465

      end subroutine load_data_frame

      
      
      subroutine unpack_snap_data(snap1d, lay, deepest, buf2D)
c     ---------------------------------------------------------------------------
c     Unpack data buffer snap1d on t/u grid 
c     into buf3D(1:nx,1:ny,lay); do not pad dry elements
c     For t-grid (cell centered) data, deepest should be bottum_layer
c     For u-grid (corner centered) data, deepest should be kmu
c     ---------------------------------------------------------------------------
      real*4, intent(in)  :: snap1d(:)
      integer, intent(in) :: lay
      integer, intent(in) :: deepest(:,:)  
      real, intent(out)   :: buf2D(:,:)
c 
      integer :: vlen, i, ix,iy
c     ---------------------------------------------------------------------------     
      vlen = size(snap1d)       ! Number of active points in this field.
      i    = 0
      do iy = 1, ny
         do ix = 1, nx
            if(deepest(ix,iy).ge.lay) then
               i=i+1
               if (i>vlen) then
                  write(*,*) "unpack_tgrid: buffer violation"
                  stop 233
               endif
               buf2D(ix,iy)=snap1d(i)
            endif
         enddo
      enddo
  
c     ---------------------------------------------------------------------------
      end subroutine unpack_snap_data




      subroutine interpolate_currents(xyz, uvw, status)
c     ---------------------------------------------------------------------------
c     Perform linear interpolation in current fields into vector uvw
c     using data on RCO u-grid (data on horizontal cell vertices)
c     
c     xyz is the interpolation position as (lon[degE],lat[degN],depth[m,positive down])
c     uvw are interpolated currents in m/s, shape == real(3+)
c     status is the status of the interpolation action
c       status = 0: all OK
c       status = 1: domain violation (but extrapolation provided)
c     status = 3: dry point
c
c     Notes: this is a simplified construction, based on simple horizontal transformation
c            of grid coordinates. Interpolation is unrestricted, implying currents = 0
c            at dry vertices
c
c     TODO
c       reconstruct w
c       check RCO flux conservation scheme in relation to interpolate_currents
c     --------------------------------------------------------------------------
      real, intent(in)     :: xyz(:)  ! (lon,lat,depth)
      real, intent(out)    :: uvw(:)
      integer, intent(out) :: status

      integer              :: statu, statv, statw
      integer              :: ix,iy,iz,ix0,iy0,iz0,ix1,iy1,iz1,ibot
      real                 :: x,y,z, sx,sy,sz,corners(8)
c     ------------------------------------------ 
      if (.not.is_wet(xyz)) then
         uvw    = 0.0
         status = 3  ! signal dry point
         return
      endif
      
c.....transform to continuous node-centered grid coordinates
c     and check ranges for this staggering
      call get_grid_coordinates(xyz,x,y,z)
      
      x = x - 0.5   ! transform for t to u-grid coordinates (data on )
      y = y - 0.5   ! transform for t to u-grid coordinates
      statu = 0
      statv = 0 
      statw = 0
      if ((x<0.5).or.(x>(nx-0.5))) statu = 1 ! flag x range violation
      if ((y<0.5).or.(y>(ny-0.5))) statv = 1 ! flag y range violation
      if ((z<0.5).or.(z>(nz+0.5))) statw = 1 ! flag z range violation
c
c     Currently ignore the mismatch between interpolation and simulation domain 
c
      status  = max(statu, statv, statw) 
    
c.....determine u-cell associations of point, constrained to 1 <= i <= n
      ix   = max(1, min(nint(x), nx)) 
      iy   = max(1, min(nint(y), ny)) 
      ibot = kmu(ix,iy) ! 1 <= ibot <= nz
      iz   = max(1, min(nint(z),ibot)) ! w at cell upper face
      
c.....determine vertices + relative intracell coorninate s, constrained to 0 < s < 1 
c     when point exceed coordinate bounds s is assigned 0 or 1
      ix0 =  max(1, min(int(x), nx))   ! truncate to get lower u-grid vertex
      iy0 =  max(1, min(int(y), ny))   ! truncate to get lower u-grid vertex
      iz0 =  max(1, min(int(z), ibot)) ! truncate to get lower u-grid vertex
      ix1 =  min(ix0+1, nx)   
      iy1 =  min(iy0+1, ny)   
      iz1 =  max(iz0+1, ibot)
      sx = max(0.d0, min(1.d0, x-real(ix0)))
      sy = max(0.d0, min(1.d0, y-real(iy0)))           
      sz = max(0.d0, min(1.d0, z-real(iz0)))  
c
      corners(1) = u(ix0,iy0,iz0)
      corners(2) = u(ix0,iy0,iz1)
      corners(3) = u(ix0,iy1,iz0)
      corners(4) = u(ix0,iy1,iz1)
      corners(5) = u(ix1,iy0,iz0)
      corners(6) = u(ix1,iy0,iz1)
      corners(7) = u(ix1,iy1,iz0)
      corners(8) = u(ix1,iy1,iz1)
      call interp_3Dbox_data(sx,sy,sz,corners,0,uvw(1))
c      
      corners(1) = v(ix0,iy0,iz0)
      corners(2) = v(ix0,iy0,iz1)
      corners(3) = v(ix0,iy1,iz0)
      corners(4) = v(ix0,iy1,iz1)
      corners(5) = v(ix1,iy0,iz0)
      corners(6) = v(ix1,iy0,iz1)
      corners(7) = v(ix1,iy1,iz0)
      corners(8) = v(ix1,iy1,iz1)
      call interp_3Dbox_data(sx,sy,sz,corners,0,uvw(2))
c     
      uvw(3) = 0.0   ! this implementation
c
c     ------------------------------------------ 
      end subroutine interpolate_currents



      end module

      
