ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     POLCOMS+ERSEM pbi for daily averaged data sets
c     ---------------------------------------------------
c     $Rev: 228 $
c     $LastChangedDate: 2011-01-26 02:52:59 +0100 (Wed, 26 Jan 2011) $
c     $LastChangedBy: asch $ 
c
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module physical_fields

      use mesh_grid
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

c      public :: interpolate_salty   ! currently unused
c      public :: interpolate_wind    ! currently unused
      public :: interpolate_zooplankton
      public :: interpolate_wdepth
      public :: is_wet    
      public :: is_land
      public :: horizontal_range_check
      public :: coast_line_intersection
      public :: get_pbi_version

c     -------------------- module data --------------------  
      
      
c
c     ------ data frame handler ------
      
      character*999              :: hydroDBpath ! hydrographic data sets

      integer, parameter         :: tag_lenght = 7 ! format == YYYY.MM
      character*(*), parameter   :: not_set_c = "not set"
      integer, parameter         :: not_set_i = -9999
     
      logical                    :: data_in_buffers     = .false.  ! first time setting
      logical                    :: data_files_are_open = .false.  ! first time setting
      character(len=tag_lenght)  :: cur_tag      ! tag YYYYMM of current open set
      integer                    :: cur_frame    ! current day-in-month time frame in buffer
      integer                    :: ncid         ! NetCDF file handler
      
c     --- 3D grids ---
      real, allocatable          :: zoobuf(:,:,:) 

c     --- 2D grids ---    
c     --- 1D grids ---
c     --- units    ---    
      real, parameter            :: mgCm3_to_kgDWm3 = 1.0e-6/0.32 ! conversion factor for mgC/m3 -> kgDW/m3

c     ===================================================
                            contains
c     ===================================================
  
      subroutine init_physical_fields(time)
c     ------------------------------------------
c     Do not trigger data load
c     ------------------------------------------
      type(clock), intent(in),optional :: time
      real,parameter :: very_deep = 1000.0
c     ------------------------------------------
      if (present(time)) call set_master_clock(time)
      write(*,*) trim(get_pbi_version()) 

      call read_control_data(simulation_file,"hydroDBpath",hydroDBpath)
      write(*,*) "init_physical_fields: hydrographic database path =", 
     +           trim(hydroDBpath) 
      call reset_frame_handler()
      call load_grid_desc()  ! invokes init_horiz_grid_transf
      call init_mesh_grid()  ! allocate core arrays  
      
      allocate( zoobuf(nx, ny, nz) )
c     ---- specials for this data set / version ----

      
      dslm             = 0.0        ! NA for sigma grids
      bottom_layer     = 0          ! default dry
      
      ccdepth(:,:,1)   = very_deep 
      ccdepth(:,:,2:)  = 2*very_deep
      acc_width(:,:,1) = 0.0        ! including surface layer (iz=1)
      acc_width(:,:,2) = 2*very_deep
      acc_width(:,:,3:)= 2*very_deep
      
      zoo  = 0.0 
      
      end subroutine init_physical_fields

 

      character*100 function get_pbi_version()  
      get_pbi_version="POLCOMS+ERSEM offline pbi version: $Rev: 228 $"
      end function



      subroutine close_physical_fields()
c     ------------------------------------------ 
c     ------------------------------------------    
      call close_mesh_grid()
      call reset_frame_handler()
      call close_horiz_grid_transf()
c     --- local cleanup ---
      deallocate( zoobuf )
c     ------------------------------------------
      end subroutine 


      subroutine update_physical_fields(time, dt)
c     ------------------------------------------  
      type(clock), intent(in),optional :: time
      integer, intent(in),optional     :: dt
      logical                          :: update
      character(len=tag_lenght)        :: tag
      integer                          :: frame ! number in that set
      type(clock), pointer             :: aclock
c     ------------------------------------------  
      aclock => get_master_clock()
      if (present(time)) then
         call set_master_clock(time)
      elseif (present(dt)) then
         call add_seconds_to_clock(aclock, dt)
      endif
c
      call resolve_corresp_dataset(aclock, tag, frame)
      call update_dataset(tag, frame)
c     ------------------------------------------       
      end subroutine update_physical_fields




      subroutine resolve_corresp_dataset(aclock, tag, frame)
c     ------------------------------------------
c     Resolve tag of necessary data set and frame number in that data set
c     corresponding to aclock
c
c     tag POLCOMS+ERSEM pbi for daily averaged data sets are YYYY.MM
c     and data frame is the day-in-month
c     ------------------------------------------
      type(clock), intent(in)                :: aclock
      character(len=tag_lenght), intent(out) :: tag
      integer, intent(out)                   :: frame

      integer                                :: year, month
c     ------------------------------------------
      call get_date_from_clock(aclock, year, month, frame)
      write(tag,455) year, month 
 455  format(i4.4,".",i2.2) !  YYYY.MM
      end subroutine resolve_corresp_dataset



      subroutine update_dataset(tag, frame)
c     ------------------------------------------
c     Only open new data set file if necessary
c     Only load new frames if necessary
c     ------------------------------------------
      character(len=tag_lenght), intent(inout) :: tag
      integer, intent(in)                      :: frame
      logical  :: needed
c     ------------------------------------------
      needed = .not.data_in_buffers
c
      if (tag /= cur_tag) then
         call reset_frame_handler() ! resets cur_frame+close old
         call open_data_files(tag)
         needed = .true. ! force update if accidentally frame == cur_frame
      endif
c      
      if ((frame /= cur_frame).or.needed) then
         call load_data_frames(frame)
      endif  
c 
      end subroutine update_dataset


      
      subroutine load_grid_desc()
c     ------------------------------------------
c     Should only be invoked at start
c     ------------------------------------------
      character(len=tag_lenght) :: tag
      integer                   :: dimid,varid,idum
      real,allocatable          :: lon(:),lat(:)
      type(clock), pointer      :: aclock
      real                      :: lambda1, dlambda, phi1, dphi ! LOCAL DUMMIES
c     ------------------------------------------
      if (data_in_buffers) stop "load_grid_desc:unexpected"
      aclock => get_master_clock()
      call resolve_corresp_dataset(aclock, tag, idum)
      call open_data_files(tag) ! -> ncid
c     --- fetch dimensions ---
      call NetCDFcheck( nf90_inq_dimid(ncid, "lon",  dimid) )
      call NetCDFcheck( nf90_inquire_dimension(ncid, dimid, len=nx) )
      call NetCDFcheck( nf90_inq_dimid(ncid, "lat",  dimid) )
      call NetCDFcheck( nf90_inquire_dimension(ncid, dimid, len=ny) )
      call NetCDFcheck( nf90_inq_dimid(ncid, "z",    dimid) )
      call NetCDFcheck( nf90_inquire_dimension(ncid, dimid, len=nz) )
      allocate(lon(nx))
      allocate(lat(ny))
      call NetCDFcheck( nf90_inq_varid(ncid, "lon", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, lon))
      call NetCDFcheck( nf90_inq_varid(ncid, "lat", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, lat))
c     --- probe grid dsecriptors from lon/lat arrays ---
      lambda1 = lon(1)
      dlambda = lon(2)-lon(1)
      phi1    = lat(1)
      dphi    = lat(2)-lat(1)
      deallocate(lon)
      deallocate(lat)

      write(*,231) nx, ny, nz
      write(*,232) lambda1, dlambda
      write(*,233) phi1,    dphi 
    
      write(*,*) "load_grid_desc: horizontal initialization"
      call init_horiz_grid_transf(lambda1, phi1, dlambda, dphi)

 229  format(a,a)    
 231  format("load_grid_desc: 3d grid dim (nx,ny,nz) = ", i4,i4,i4)
 232  format("load_grid_desc: lambda1 = ",f12.7,
     +                      " dlambda = ",f12.7," degrees")  
 233  format("load_grid_desc: phi1    = ",f12.7,
     +                      " dphi    = ",f12.7," degrees")  
      
      end subroutine load_grid_desc
      



      subroutine reset_frame_handler()
c     ------------------------------------------
      call close_data_files()
      data_in_buffers = .false.
      cur_tag         = not_set_c
      cur_frame       = not_set_i   
      end subroutine reset_frame_handler



      subroutine close_data_files()
c     ------------------------------------------
c     Only close, if data files are open
c     
c     ------------------------------------------
      if (data_files_are_open) then
         call NetCDFcheck( nf90_close(ncid) )
         ncid = not_set_i          
         write(*,*) "close_data_files: closed open NetCDF set"
         data_files_are_open = .false.
      endif
      end subroutine close_data_files


      subroutine open_data_files(tag)
c     ------------------------------------------
c     Open netcdf file set corresponding to 
c     tag == YYYY.MM 
c     
c     hydroDBpath/Daily.PolcomsErsem.YYYY.MM.SLAM.nc  -> ncid
c     
c     ------------------------------------------
      character(len=tag_lenght),intent(in) :: tag
      character*999                        :: fname
      integer :: timedimID, ntime, varid
c     ------------------------------------------
      call close_data_files() ! in case they are open ...
      write(fname,333) trim(adjustl(hydroDBpath)), tag
      call NetCDFcheck( nf90_open(fname, NF90_NOWRITE, ncid) )
 333  format(a,'/Daily.PolcomsErsem.',a,'.SLAM.nc')
 
      write(*,*) "open_data_files: opened NetCDF set ", 
     +           trim(adjustl(fname))
c  
      data_files_are_open = .true.
      cur_tag = tag

      end subroutine open_data_files



      subroutine load_data_frames(frame)
c     ------------------------------------------
c     Verbose reading af NetCDF sets consecutively to grid buffers 
c     from currently open NetCDF set corresponding to frame (day-in-month)
c
c     Syncronize auxillary grid fields after data load. 
c     It is assumed that grids do not
c     change during a simulation (assertion not checked)
c     
c     Currently load these fields (in native units/conventions):
c
c     float u(time, depth, lat, lon) ;  "m/s"  positive east
c     float v(time, depth, lat, lon) ;  "m/s"  positive north 
c     float w(time, depth, lat, lon) ;  "m/s", positive up 
c     float z(time, lat, lon) ;         "m"    positive down (!)
c     float t(time, depth, lat, lon) ;  "degC" 
c     float zoo(time, depth, lat, lon) ;"Zooplankton: kg DW/m3" 
c                
c     The is currently no horizontal diffusivity in the set (set it to zero)
c
c     fortran index order opposite CDL dump order
c     fortran: fastests index left most
c     ------------------------------------------
      integer, intent(in) :: frame ! pick frame, corresponding to frame hours since ref_clock 
      
      integer :: start2D(3),start3D(4),start4D(5),varid
      integer :: count2D(3),count3D(4),count4D(5)
      integer :: ix,iy,iz,no_fill
      real    :: dz
      logical :: not_ok
      real,pointer  :: cc(:), acc(:)
c     ------------------------------------------ 
      write(*,*) "load_data_frames: loading frame", frame
c     ---- 3D ----
      start3D(:) = 1
      start3D(4) = frame
c     ---- load u current component ---- 
      call NetCDFcheck( nf90_inq_varid(ncid, "ucurD", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, u, start=start3D) )
c     ---- load v current component ---- 
      call NetCDFcheck( nf90_inq_varid(ncid, "vcurD", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, v, start=start3D) )
c     ---- load w current component ---- 
c          (w positive downward in data set)
      call NetCDFcheck( nf90_inq_varid(ncid, "wcurD", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, w, start=start3D) )
c     ---- load temperature ----
      call NetCDFcheck( nf90_inq_varid(ncid, "ETWD", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, temp, start=start3D)) 
c     ---- load vertical turbulent diffusivity ----
c          notice that generic vertical allocation is nz+1, because
c          some data sets define data points at vertical faces
      call NetCDFcheck( nf90_inq_varid(ncid, "nuvD", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid,vdiffus(:,:,1:nz),
     +                  start=start3D)) ! wrong
      vdiffus(:,:,nz+1) = 1.e-9 ! render all elements defined
c     ---- load zooplankton ----
c          currently zooplankton is set to the sum of micro+meso  zooplankton    
      call NetCDFcheck( nf90_inq_varid(ncid, "Z4cD", varid) ) ! Mesozooplankton
      call NetCDFcheck( nf90_get_var(ncid, varid, zoo, start=start3D)) 
      call NetCDFcheck( nf90_inq_varid(ncid, "Z5cD", varid) ) ! Microzooplankton
      call NetCDFcheck( nf90_get_var(ncid,varid,zoobuf,start=start3D)) 
      zoo = zoo + zoobuf        ! add micro and mesozooplankton
      zoo = zoo*mgCm3_to_kgDWm3 ! convert biomass to unit kg DW / m3
c
c          NB: horizontal derivatives NOT implemented
c          when spatially non constant hdiffus are applied,
c          interpolate_turbulence_deriv must be updated
c
c     ---- currently no horizontal turbulent diffusivity: set to
c          molecular diffusivity lower limit ~ 1.e-9 m2/s   ---
      where (vdiffus<1.e-9)
         vdiffus = 1.e-9
      end where
      hdiffus = 1.e-9   
  
c     ---- load cell center depths (positive below water)---- 
c          depths are positive below water, so flip sign of input data, which
c          has opposite convention
      call NetCDFcheck( nf90_inq_varid(ncid, "depth", varid) )
      call NetCDFcheck( nf90_get_var(ncid,varid,ccdepth,start=start3D))
      ccdepth = -ccdepth  ! flip sign so wdepth is positive in water

c     ---- load vertical cell face depths (only used to defined total water depth) ----
c     ---- grid is sigma type so all vertical cells are wet, if any
c          depths are positive below water, so flip sign of input data, which
c          has opposite convention
c                               1       2    3   4   5
c          fortran-def: zbnd(n_zfaces, lon, lat, z, time) ;
c     
      start4D(:) = 1
      start4D(1) = 2 ! pick lower face
      start4D(4) = nz  ! bottom layer
      start4D(5) = frame
      count4D      = 1
      count4D(2:3) = shape(wdepth)
c
      call NetCDFcheck( nf90_inq_varid(ncid, "zbnd", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, wdepth, 
     +                  start=start4D, count=count4D) )
      do ix=1,nx
         do iy=1,ny
            if (abs(wdepth(ix,iy))<1.e-6) wdepth(ix,iy) = 1.0 ! we flip sign below
            if (is_NetCDF_fill_value_real(wdepth(ix,iy)))
     +          wdepth(ix,iy) = 1.0                           ! we flip sign below
         enddo
      enddo
      wdepth = -wdepth ! flip sign of whole array so wdepth is positive in water

c     --- define the wetmask auxillary  ---    
      where(wdepth>0)
         wetmask      = 1 ! wet
         bottom_layer = nz
      elsewhere
         wetmask      = 0 ! dry
         bottom_layer = 0
      end where 

c     ------ generate auxillary grid descriptors for wet points  
c            acc_width(ix,iy,1) = 0 (sea surface)  
      do ix=1,nx
         do iy=1,ny
            if (wetmask(ix,iy)>0) then
               cc  => ccdepth(ix,iy,:)
               acc => acc_width(ix,iy,:)
               call levels2accwidth(cc, acc)
            endif      
         enddo
      enddo
c
c     ---- enforce exact consistency between acc_width(:,:,nz+1) and wdepth 
c          (there are numerical diffs of order 10^-5)
c
      acc_width(:,:,nz+1) = wdepth 
           

      write(*,*) "load_data_frames: loaded frame:", frame
      data_in_buffers = .true.
      cur_frame = frame

c$$$      do ix=1,nx
c$$$         do iy=1,ny
c$$$            if (wetmask(ix,iy)==0) 
c$$$     +          write(77,*) (ix-1)*dlambda + lambda1,
c$$$     +                     (iy-1)*dphi    + phi1
c$$$            if (wetmask(ix,iy)==0) write(77,*) ix,iy
c$$$         enddo
c$$$      enddo
c$$$      stop 84465

      end subroutine load_data_frames
      


      subroutine levels2accwidth(levels, faces)
c     ------------------------------------------
c     Compute position of layer faces from center levels
c     ------------------------------------------
      real, intent(in)  :: levels(:)  ! 1:nz
      real, intent(out) :: faces(:)   ! 1:nz+1
      integer :: iz
      real    :: layer_width
c     ------------------------------------------ 
      faces(1) = 0
      do iz = 1,size(levels)
         layer_width   = 2.0*(levels(iz) - faces(iz))
         faces(iz+1)   = faces(iz) + layer_width
      enddo
      end subroutine levels2accwidth


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


      logical function is_NetCDF_fill_value_real(value)
c     ------------------------------------------------------
c     NetCDF auxillary query function (temporary implementation)
c     nf90_fill_real is defined in netcdf.mod 
c     ------------------------------------------------------
      real :: value
      is_NetCDF_fill_value_real = (value>0.999*nf90_fill_real)
      end function is_NetCDF_fill_value_real


      end module

      
