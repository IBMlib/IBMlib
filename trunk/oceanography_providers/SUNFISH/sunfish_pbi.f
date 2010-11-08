ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     PBI interface for SUNFISH cmod + NPZD data
c
c     Coast line reconstruction: this module features
c     a consistent definition of coast_line_intersection()
c     interpolate_wdepth(), and is_land() so that coast line is
c     characterized by wdepth=0 piecewise lines
c
c     Uses site installation of NetCDF
c
c     TODO:  
c       Currently turbulence + deriv is fixed to zero - update
c       test interpolations
c       validate w sign
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module physical_fields
      use time_tools           ! import clock type
      use constants
      use run_context, only: simulation_file
      use input_parser
      use netcdf   ! use site installation

      implicit none
      private     

      interface get_local_distance
         module procedure get_local_distance_scalar
         module procedure get_local_distance_vector      
      end interface


      interface is_NetCDF_fill_value
         module procedure is_NetCDF_fill_value_real     
      end interface


      public :: init_physical_fields     
      public :: close_physical_fields
      public :: get_master_clock
      public :: set_master_clock
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
      public :: get_horizontal_distance
      public :: get_local_distance
      public :: coast_line_intersection
      public :: d_cart2d_xy
      public :: d_xy2d_cart 
      public :: add_finite_step

c     -------------------- module data --------------------  
      
      type(clock), target :: master_clock

c     ------ grid dimensions:   ------

      integer :: nx,ny,nz
      real    :: lambda1, dlambda ! (lambda1,phi1) is the position of 
      real    :: phi1,    dphi    ! the first grid point (ix,iy)=(1,1)
c
c     ------ data frame handler ------
      
      character*999              :: hydroDBpath ! hydrographic data sets

      integer, parameter         :: tag_lenght = 10 ! YYYYMMDDHH
      character*(*), parameter   :: not_set_c = "not set"
      integer, parameter         :: not_set_i = -9999
      type(clock)                :: ref_clock  ! reference offset for hydrographic data

      logical                    :: data_in_buffers     = .false.  ! first time setting
      logical                    :: data_files_are_open = .false.  ! first time setting
      character(len=tag_lenght)  :: cur_tag      ! tag YYYYMMDDHH of current open set
      integer                    :: cur_h1900    ! current time frame in buffer
      integer,allocatable        :: h1900_map(:) ! map of time frames in current open set
      integer                    :: ncid_u       ! u NetCDF file handler
      integer                    :: ncid_v       ! v NetCDF handler
      integer                    :: ncid_w       ! w NetCDF handler
      integer                    :: ncid_dslm    ! sea leve  l NetCDF handler
      integer                    :: ncid_t       ! temperature NetCDF handler
      integer                    :: ncid_zoo     ! zooplankton NetCDF handler 
c      integer                    :: ncid_s     : salinity currently unused
c      integer                    :: ncid_light : light currently unused

c     --- 3D grids ---
             
      real,allocatable :: u(:,:,:)          ! u of current [m/s] (positive east)
      real,allocatable :: v(:,:,:)          ! v of current [m/s] (positive north)       
      real,allocatable :: w(:,:,:)          ! w of current [m/s] (positive down)   
      real,allocatable :: temp(:,:,:)       ! Water temp. [Celcius]
      real,allocatable :: vdiffus(:,:,:)    ! vertical   diffusivity [m**2/s]              
      real,allocatable :: hdiffus(:,:,:)    ! horizontal diffusivity [m**2/s]
      real,allocatable :: dslm(:,:)         ! current sea surface elevation over reference [m]
      real,allocatable :: zoo(:,:,:)        ! Zooplankton [10^-3 mol N/liter]
      real,allocatable,target :: ccdepth0(:,:,:)   ! reference cell center depth water below surface [m] (dslm=0)
      real,allocatable,target :: ccdepth(:,:,:)    ! cell center depth water below surface [m]; pointer
      real,allocatable,target :: acc_width0(:,:,:) ! accumulated water above this undisturbed layer [m] dim=nz+1   
      real,allocatable,target :: acc_width(:,:,:)  ! accumulated water above this layer [m] dim=nz+1  

c     --- 2D grids ---
      
      real,allocatable :: wdepth0(:,:)      ! reference depth at cell-center, positive  [m] - dry points negative
      real,allocatable :: wdepth (:,:)      ! current depth at cell-center, including dslm [m]
      real,allocatable :: wetmask(:,:)      ! auxillary to define coastal geometry (1: wet nodes, 0 dry)  
      integer, allocatable :: bottom_layer(:,:) ! last wet layer (0 for dry points) nx,ny

c     --- 1D grids ---

      real, allocatable :: levels(:)        ! center depth of layers for DSLM=0. [m] dim=nz 
      real, allocatable :: layer_width(:)   ! width for undisturbed layers [m] dim=nz 
      real, allocatable :: layer_faces(:)   ! accumulated width for undisturbed layers above this [m] dim=nz+1 
      
c     ===================================================
                            contains
c     ===================================================
  
      subroutine init_physical_fields(time)
c     ------------------------------------------
c     Do not trigger data load
c     ------------------------------------------
      type(clock), intent(in),optional :: time
c     ------------------------------------------
      if (present(time)) master_clock = time

      call read_control_data(simulation_file,"hydroDBpath",hydroDBpath)
      write(*,*) "init_physical_fields: hydrographic database path =", 
     +           trim(hydroDBpath)
    
      call read_grid_desc()
c.....Set reference offset for hydrographic data: 1900-01-01 00:00:00.
      call set_clock(ref_clock, 1900, 01, 01, 0)

      call reset_frame_handler()

      end subroutine init_physical_fields

 

   

      subroutine read_grid_desc()
c     ---------------------------------------------------
c     Read grid definition which is contained in file 
c     given as tag grid_desc in the simulation file simulation_file
c     Allocate grid arrays in this subroutine. 
c
c     grid dimensions:         nx,ny,nz
c     grid coordinate map:     lambda1, dlambda; phi1, dphi
c
c     grid point (ix,iy) = (1,1) is at (lambda1,phi1)
c     ---------------------------------------------------
      character*999        :: gd_fname
      type(control_file)   :: grid_ctrlfile
      integer              :: ix,iy,iz,ihit,idx,nwet,nacc,nbott
      real                 :: wd0
c     ---------------------------------------------------
      call read_control_data(simulation_file, "grid_desc", gd_fname)
      write(*,229) "read_grid_desc: reading grid description from: ", 
     +           trim(adjustl(gd_fname))

      call open_control_file(gd_fname, grid_ctrlfile)
    
c.....grid scale/dimensions
      
      call read_control_data(grid_ctrlfile,"lambda1", lambda1)
      call read_control_data(grid_ctrlfile,"dlambda", dlambda) 
      call read_control_data(grid_ctrlfile,"phi1",    phi1)
      call read_control_data(grid_ctrlfile,"dphi",    dphi) 
      call read_control_data(grid_ctrlfile,"nx",      nx) 
      call read_control_data(grid_ctrlfile,"ny",      ny) 
      call read_control_data(grid_ctrlfile,"nz",      nz) 
      write(*,231) nx, ny, nz
      write(*,232) lambda1, dlambda
      write(*,233) phi1,    dphi 
 229  format(a,a)    
 231  format("read_grid_desc: 3d grid dim (nx,ny,nz) = ", i4,i4,i4)
 232  format("read_grid_desc: lambda1 = ",f12.7,
     +                      " dlambda = ",f12.7," degrees")  
 233  format("read_grid_desc: phi1    = ",f12.7,
     +                      " dphi    = ",f12.7," degrees")  
     

c     ------ allocate grid arrays ------

      write(*,*) "read_grid_desc: allocate grid arrays: begin"

c     --- 3D grids ---
      allocate( u(nx,ny,nz)       )   
      allocate( v(nx,ny,nz)       )   
      allocate( w(nx,ny,nz)       )   
      allocate( temp(nx,ny,nz)    )  
      allocate( vdiffus(nx,ny,nz) )                      
      allocate( hdiffus(nx,ny,nz) ) 
      allocate( zoo(nx,ny,nz)     ) 
    
      allocate( ccdepth0(nx,ny,nz)) 
      allocate( ccdepth(nx,ny,nz) ) 
      allocate( acc_width0(nx,ny,nz+1)  )  
      allocate( acc_width(nx,ny,nz+1)   )

c     --- 2D grids ---

      allocate( dslm(nx,ny)       )       
      allocate( wdepth0(nx,ny)    ) 
      allocate( wdepth(nx,ny)     ) 
      allocate( wetmask(nx,ny)    )  ! not static
      allocate( bottom_layer(nx,ny) )

c     --- 1D grids ---
      
      allocate( levels(nz)        ) 
      allocate( layer_width(nz)   )  ! width for undisturbed layers in meters
      allocate( layer_faces(nz+1) )  ! meters of undisturbed water above layer
      
      write(*,*) "read_grid_desc: allocate grid arrays: OK"   


c.....vertical layer structure: levels -> layer_width0 + acc_width0

      call read_control_data(grid_ctrlfile,"levels", levels) 
      write(*,241) levels
 241  format("read_grid_desc: read levels = ",1000f6.1)
     
c     set acc_width0, layer_width0 from levels
c     last element acc_width0(nz+1) is total (dslm=0) max depth of grid
      layer_faces(1)    = 0.
      acc_width0(:,:,1) = 0. ! sea surface
      do iz=1,nz 
         layer_width(iz)      = 2.0*(levels(iz) - layer_faces(iz))
         layer_faces(iz+1)    = layer_faces(iz) + layer_width(iz)
         ccdepth0(:,:,iz)     = levels(iz)
         acc_width0(:,:,iz+1) = layer_faces(iz+1)
      enddo
      

      write(*,242) layer_width
      write(*,243) layer_faces
 242  format("read_grid_desc: calc ref layer_width = ",1000f6.1)     
 243  format("read_grid_desc: calc ref layer_faces = ",1000f6.1)
c
c.....read topography -> wdepth0 (only wety points provided, other assumed dry)
c
c     1) set bottom_layer
c     2) adjust ccdepth0, so that actual width of bottom layer is consistent 
c        with wdepth0 by stretching the last wet layer to the topographic bottom
c         
c
      write(*,*) "read_grid_desc: reading topography"
      wdepth0      = -1.0       ! default dry
      bottom_layer =  0         ! default dry
      ihit         =  1         ! start scanning for is_wet at ihit=1
      call read_control_data(grid_ctrlfile,"topography_points",nwet,
     +                       ihit)                     ! locate file position: ihit
      nacc   = 0
      nbott  = 0
      do idx = 1, nwet
         read(grid_ctrlfile%rawdata(ihit+idx),*) ix, iy, wd0
c        --- horizontal check ---
         if ((ix>nx).or.(ix<1)) then
            write(*,*) "topography point", idx, 
     +                 " exceeds x limits - point ignored"
            cycle
         endif
         if ((iy>ny).or.(iy<1)) then
            write(*,*) "topography point", idx, 
     +                 " exceeds y limits - point ignored"
            cycle
         endif
c
c        locate vertical cell iz containing wd0: 
c        layer_faces(iz) < wd0 < layer_faces(iz+1) 
c        iz=0    -> wd0 negative
c        iz=nz+1 -> wd0 below sea bed
c
         call search_sorted_list(wd0, layer_faces, iz) 

c        --- check of iz  ---
         if (iz<1) then
            write(*,*) "topography point", idx, 
     +                 "has negative depth  - point ignored"
            cycle
         endif
         if (iz>nz) then
            write(*,*) "topography point", idx, 
     +                 " exceeds grid z limit - bottum lowered"           
            iz     = nz
            nbott  = nbott + 1 
         endif
c
c        --- accept this point, if not bounced above ---
c            after here we know iz is appropriate: 1<=iz<=nz 
c            adjust grid bottum by stretching the last wet layer
         nacc                    = nacc + 1
         wdepth0(ix, iy)         = wd0
         bottom_layer(ix, iy)    = max(iz-1, 1)     
         acc_width0(ix, iy,iz+1) = wd0              
         ccdepth0(ix, iy, iz)    = 0.5*(acc_width0(ix, iy,iz) + 
     +                                  acc_width0(ix, iy,iz+1)) ! re-center mid point of bottom cell
c
      enddo ! idx = 1, nwet
      write(*,288) nwet, nacc
      write(*,289) nbott
 288  format("read_grid_desc: read",i6," wet points - ",i6, 
     +       " points accepted")
 289  format("read_grid_desc: grid bottum lowered in ",i6, " points")     
c
c     --- define the wetmask statically (no flooding/drying) --- 
c   
      where(wdepth0 > 0.0)
         wetmask = 1.0 ! wet
      elsewhere
         wetmask = 0.0 ! dry
      end where 

      write(*,*) "read_grid_desc: grid arrays initialized"

      call close_control_file(grid_ctrlfile)
     
c     ------------------------------------------       
      end subroutine read_grid_desc



      subroutine close_physical_fields()
c     ------------------------------------------ 
c     ------------------------------------------    
      if (allocated(h1900_map))    deallocate( h1900_map )  
      if (allocated(u))            deallocate( u )
      if (allocated(v))            deallocate( v )
      if (allocated(w))            deallocate( w )
      if (allocated(temp))         deallocate( temp )
      if (allocated(vdiffus))      deallocate( vdiffus )
      if (allocated(hdiffus))      deallocate( hdiffus )
      if (allocated(dslm))         deallocate( dslm )
      if (allocated(zoo))          deallocate( zoo )
      if (allocated(ccdepth0))     deallocate( ccdepth0 )
      if (allocated(ccdepth))      deallocate( ccdepth )
      if (allocated(acc_width0))   deallocate( acc_width0 )
      if (allocated(acc_width))    deallocate( acc_width )
      if (allocated(wdepth0))      deallocate( wdepth0 )
      if (allocated(wdepth))       deallocate( wdepth )
      if (allocated(wetmask))      deallocate( wetmask )
      if (allocated(bottom_layer)) deallocate( bottom_layer )
      if (allocated(levels))       deallocate( levels )
      if (allocated(layer_width))  deallocate( layer_width )
      if (allocated(layer_faces))  deallocate( layer_faces )
      call reset_frame_handler()
c     ------------------------------------------
      end subroutine 



      function   get_master_clock()
c     ------------------------------------------ 
      type(clock), pointer :: get_master_clock
c     ------------------------------------------ 
      get_master_clock => master_clock
c     ------------------------------------------ 
      end function 


      subroutine set_master_clock(time)
c     ------------------------------------------ 
      type(clock), intent(in) :: time
c     ------------------------------------------ 
      master_clock = time ! local copy
c     ------------------------------------------ 
      end subroutine 



      subroutine update_physical_fields(time, dt)
c     ------------------------------------------  
      type(clock), intent(in),optional :: time
      integer, intent(in),optional     :: dt
      logical                          :: update
      character(len=tag_lenght)        :: tag
      integer                          :: h1900
c     ------------------------------------------  
      if (present(time)) then
         call set_master_clock(time)
      elseif (present(dt)) then
         call add_seconds_to_clock(master_clock, dt)
      endif
c
      call resolve_corresp_dataset(master_clock, tag, h1900)
      call update_dataset(tag, h1900)
c     ------------------------------------------       
      end subroutine update_physical_fields


      subroutine resolve_corresp_dataset(aclock, tag, h1900)
c     ------------------------------------------
c     Resolve tag and hours since 1900-01-01 00:00:00" (== dataframe)
c     corresponding to aclock
c
c     Data frames are stored according to end point of
c     the hour interval, i.e time = 912385h represents
c     time interval ]912384h; 912385h] since 1900-01-01 00:00:00.
c     [ref == mam@dmu.dk; ma 11-10-2010 12:13]
c
c     Tag YYYYMMDDHH refers to the beginning of that day
c     ------------------------------------------
      type(clock), intent(in)                :: aclock
      character(len=tag_lenght), intent(out) :: tag
      integer, intent(out)                   :: h1900
      integer          :: year, month, day
      real             :: h1900_r
c     ------------------------------------------
      call get_date_from_clock(aclock, year, month, day)
      write(tag,455) year, month, day, 0   ! beginning of that day
      call get_period_length_hour(ref_clock, aclock, h1900_r)
      h1900 = nint(h1900_r + 0.5)     
 455  format(i4.4,3i2.2)
      end subroutine resolve_corresp_dataset




      subroutine update_dataset(tag, h1900)
c     ------------------------------------------
c     Only open new data set file if necessary
c     Only load new frames if necessary
c     ------------------------------------------
      character(len=tag_lenght), intent(inout) :: tag
      integer, intent(inout)                   :: h1900
      logical  :: needed
c     ------------------------------------------
      needed = .not.data_in_buffers
c
      if ((tag /= cur_tag).or.needed) then
         call reset_frame_handler() ! resets cur_h1900+close old
         call open_data_files(tag)
         cur_tag = tag
      endif
c      
      if ((h1900 /= cur_h1900).or.needed) then
         call load_data_frames(h1900)
      endif  
c 
      end subroutine update_dataset





      subroutine reset_frame_handler()
c     ------------------------------------------
      call close_data_files()
      data_in_buffers = .false.
      cur_tag         = not_set_c
      cur_h1900       = not_set_i   
      end subroutine reset_frame_handler



      subroutine close_data_files()
c     ------------------------------------------
c     Only close, if data files are open
c     
c     ------------------------------------------
      if (data_files_are_open) then
         call NetCDFcheck( nf90_close(ncid_u) )
         call NetCDFcheck( nf90_close(ncid_v) )
         call NetCDFcheck( nf90_close(ncid_w) )
         call NetCDFcheck( nf90_close(ncid_dslm) )
         call NetCDFcheck( nf90_close(ncid_t) )
         call NetCDFcheck( nf90_close(ncid_zoo) )
         ncid_u    = not_set_i  
         ncid_v    = not_set_i  
         ncid_w    = not_set_i  
         ncid_dslm = not_set_i  
         ncid_t    = not_set_i  
         ncid_zoo  = not_set_i          
         write(*,*) "close_data_files: closed open NetCDF set"
         data_files_are_open = .false.
      endif
      end subroutine close_data_files


      subroutine open_data_files(tag)
c     ------------------------------------------
c     Open netcdf file set corresponding to 
c     tag == YYYYMMDDHH and load the
c     time map h1900_map of data in this set
c     
c     tempdat.u_YYYYMMDDHH.00.nc  -> ncid_u
c     tempdat.v_YYYYMMDDHH.00.nc  -> ncid_v
c     tempdat.w_YYYYMMDDHH.00.nc  -> ncid_w
c     tempdat.z_YYYYMMDDHH.00.nc  -> ncid_dslm
c     tempdat.t_YYYYMMDDHH.00.nc  -> ncid_t
c     tempdat.s_YYYYMMDDHH.00.nc  -> (ncid_s     : currently unused)
c     biodat.light_YYYYMMDDHH.nc  -> (ncid_light : currently unused)
c     biodat.zoo_YYYYMMDDHH.nc    -> ncid_zoo
c     ------------------------------------------
      character(len=tag_lenght),intent(in) :: tag
      character*999                        :: fname
      integer :: timedimID, ntime, varid
c     ------------------------------------------
      call close_data_files() ! in case they are open ...
      write(fname,333) trim(adjustl(hydroDBpath)), "u", tag
      call NetCDFcheck( nf90_open(fname, NF90_NOWRITE, ncid_u) )
      write(fname,333) trim(adjustl(hydroDBpath)), "v", tag
      call NetCDFcheck( nf90_open(fname, NF90_NOWRITE, ncid_v) )
      write(fname,333) trim(adjustl(hydroDBpath)), "w", tag
      call NetCDFcheck( nf90_open(fname, NF90_NOWRITE, ncid_w) )
      write(fname,333) trim(adjustl(hydroDBpath)), "z", tag
      call NetCDFcheck( nf90_open(fname, NF90_NOWRITE, ncid_dslm) )
      write(fname,333) trim(adjustl(hydroDBpath)), "t", tag
      call NetCDFcheck( nf90_open(fname, NF90_NOWRITE, ncid_t) )
      write(fname,334) trim(adjustl(hydroDBpath)), "zoo", tag
      call NetCDFcheck( nf90_open(fname, NF90_NOWRITE, ncid_zoo) )
c
c     load  time map h1900_map of data (asserted identical for all in this set)
c
      if (allocated(h1900_map)) deallocate(h1900_map)
      call NetCDFcheck( nf90_inq_dimid(ncid_dslm, "time", timedimID) )
      call NetCDFcheck( nf90_inquire_dimension(ncid_dslm, timedimID, 
     +                  len = ntime) )
c     
      if (ntime>0) then
         allocate( h1900_map(ntime) )
         call NetCDFcheck( nf90_inq_varid(ncid_dslm, "time", varid) )
         call NetCDFcheck( nf90_get_var(ncid_dslm, varid, h1900_map))
      else
         write(*,*) "open_data_files error: NetCDF set ", tag,
     +              " contains ",ntime, "frames"
      endif
c
      write(*,*) "open_data_files: opened NetCDF set ", tag,
     +           " containing ",ntime, "frames"  
c  
      data_files_are_open = .true.
c     
 333  format(a,'/tempdat.',a,'_',a10,'.00.nc')
 334  format(a,'/biodat.',a,'_',a10,'.nc')
      end subroutine open_data_files



      subroutine load_data_frames(h1900)
c     ------------------------------------------
c     Verbose reading af NetCDF sets consecutively to grid buffers 
c     from open NetCDF sets (referenced by data frame handler)
c     corresponding to h1900 hours since 1900-01-01 00:00:00 (ref_clock)
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
c     float zoo(time, depth, lat, lon) ;"Zooplankton: 10^-3 mol N/liter" 
c                
c     The is currently no horizontal diffusivity in the set (set it to zero)
c
c     fortran index order opposite CDL dump order
c     fortran: fastests index left most
c     ------------------------------------------
      integer, intent(in) :: h1900 ! pick frame, corresponding to h1900 hours since ref_clock 
      integer :: iframe 
      integer :: start3D(4),start2D(3),varid
      integer :: ix,iy,iz
      real    :: dz
c     ------------------------------------------ 
c
c     1) locate frame to pick - time frame map h1900_map already loaded in open_data_files
c
      iframe = 1
      do while ((iframe <= size(h1900_map)).and.
     +         (h1900 /= h1900_map(iframe)))
         iframe = iframe + 1
      enddo
      if (iframe > size(h1900_map)) then
         write(*,*) "load_data_frames: unable to locate h1900 =",h1900
         write(*,*) "frames in present set =", h1900_map
         stop "load_data_frames: load error"
      endif
c
c     2) load iframe from set
c
c     ---- 2D ----
      start2D(:) = 1
      start2D(3) = iframe
      call NetCDFcheck( nf90_inq_varid(ncid_dslm,     "z",   varid) )
      call NetCDFcheck( nf90_get_var(ncid_dslm, varid, dslm, start2D) )
c     ---- 3D ----
      start3D(:) = 1
      start3D(4) = iframe
      call NetCDFcheck( nf90_inq_varid(ncid_u,     "u", varid) )
      call NetCDFcheck( nf90_get_var(ncid_u, varid, u, start3D))      
      call NetCDFcheck( nf90_inq_varid(ncid_v,     "v", varid) )
      call NetCDFcheck( nf90_get_var(ncid_v, varid, v, start3D))
      call NetCDFcheck( nf90_inq_varid(ncid_w,     "w", varid) )
      call NetCDFcheck( nf90_get_var(ncid_w, varid, w, start3D))
      call NetCDFcheck( nf90_inq_varid(ncid_t,     "t", varid) )
      call NetCDFcheck( nf90_get_var(ncid_t, varid, temp, start3D))
      call NetCDFcheck( nf90_inq_varid(ncid_zoo,     "zoo", varid) )
      call NetCDFcheck( nf90_get_var(ncid_zoo, varid, zoo, start3D))
c
c     3) postprocess data, suncronize auxillary fields, pad holes
c      
c     ---- currently no horizontal turbulent diffusivity: set to
c          molecular diffusivity lower limit ~ 1.e-9 m2/s   ---
c          NB: horizontal derivatives NOT implemented
c          when spatially non constant hdiffus are applied,
c          interpolate_turbulence_deriv must be updated
c
      vdiffus = 1.e-9
      hdiffus = 1.e-9   
c
c     ------ flip sign of w 
c     ------ physical_fields sign convention is positive down, cmod convention is positive up
c
      w = -w
c     ------ add DSLM to auxillary grid descriptors for wet points  
c            acc_width(ix,iy,1) = 0 (sea surface)  
      do ix=1,nx
         do iy=1,ny
            if (wetmask(ix,iy)<1.e-6) cycle ! nothing for dry points
            dz = dslm(ix,iy) ! change of sea level in this point, positive down
            wdepth(ix,iy)       = wdepth0(ix,iy)       - dz
            acc_width(ix,iy,2:) = acc_width0(ix,iy,2:) - dz
            ccdepth(ix,iy,1)    = ccdepth0(ix,iy,1)    - dz*0.5 
            ccdepth(ix,iy,2:)   = ccdepth0(ix,iy,2:)   - dz          
         enddo
      enddo
c      
      write(*,*) "load_data_frames: loaded frame:",
     +            iframe,"/", size(h1900_map)
      data_in_buffers = .true.
      cur_h1900 = h1900
      end subroutine load_data_frames


      subroutine search_sorted_list(v,vlist,i)
c     ------------------------------------------------------
c     Locate the position i in sorted list vlist where v should
c     appear so that  vlist(i) < v < vlist(i+1)
c     Especially return i=0, if  v < vlist(1) and
c     return i=n, if  v > vlist(n)
c     Search by bisection
c     ------------------------------------------------------
      real, intent(in)     :: v, vlist(:)
      integer, intent(out) :: i
      integer              :: ip,itest,n
c     ------------------------------------------------------
      n = size(vlist)
      if (v < vlist(1)) then
         i=0
         return
      endif
      if (v > vlist(n)) then
         i=n
         return
      endif
c     --- do loop: keep vlist(i) < v < vlist(ip) and i<ip
      i  = 1  ! we know v > vlist(1)
      ip = n  ! we know v < vlist(1)
      do while ((ip-i)>1)
         itest = (i+ip)/2 ! integer arithmetics: i<itest<ip
         if (vlist(itest)>v) then
            ip = itest
         else
            i  = itest
         endif
      enddo
c     ------------------------------------------------------
      end subroutine search_sorted_list



      subroutine interpolate_cc_3Dgrid_data(xyz,array,deriv,result,
     +                                      status)
c     -------------------------------------------------------------------- 
c     Interpolate on grid of corner centered data 3D array on point xyz.
c     Apply appropriate extrapolations near boundaries or return padvalue
c
c     deriv = 0 gives value, deriv = (1,2,3) gives derivative along (x,y,z)
c     Reserve deriv = 123 for overloaded full gradient ?
c     Currently, only deriv = (0,3) is implemented, until other derivatives 
c     are needed.
c
c     TODO:  multiply metric factor for deriv = (1,2), because
c
c     Return status:
c       status = 0: interior interpolation performed
c       status = 1: horizontal range violation, set result = padval
c       status = 2: vertical extension using boundary values
c                   for a dry point set result = padval
c       status = 3: rank deficit interpolation performed (this 
c                   is the case near the coast, where some data columns
c                   may be absent, or at land) status = 3 may include 
c                   vertical extension additionally. At interior land
c                   points, result = padval for deriv = 0
c
c     tested: deriv=0,3
c     --------------------------------------------------------------------
      real, intent(in)     :: xyz(:),array(:,:,:)
      integer, intent(in)  :: deriv  ! 0=value; (1,2,3) = along (x,y,z)
      real, intent(out)    :: result
      integer, intent(out) :: status
c      
      real, parameter     :: padval = 0. ! later move to argument list
c
      integer             :: ix,iy,iz,i,idum,ibot,ix0,ix1,iy0,iy1
      real                :: z,sx,sy,sz,depth,fbar,fup,flow
      integer             :: cx(4), cy(4), pad_mask(4)
      real                :: fz(4), dz(4), dfz(4), hweight(4)
      logical             :: valid(4)
      real, pointer       :: zgrid(:)
c     --------------------------------------------------------------------
      if (.not.horizontal_range_check(xyz)) then
          result = padval
          if (deriv>0) result = 0 
          status = 1
          return
      endif

      call interpolate_wdepth(xyz,depth,idum)
      if ((z<0).or.(z>depth)) then
         status = 2                ! signal vertical out-of-bound
         result = padval           ! derivative interpolation
         if (deriv>0) result = 0   ! value interpolation
         return
      endif
c
c     define corners clock wise
c         col1  (ix0,iy0)
c         col2  (ix0,iy1)
c         col3  (ix1,iy0)  
c         col4  (ix1,iy1)
c
      call get_horiz_grid_coordinates(xyz,ix,iy,sx,sy)
c
      ix0 = min(max(ix,  1),nx)  ! cover grid edges
      ix1 = min(max(ix+1,1),nx)  ! cover grid edges
      iy0 = min(max(iy,  1),ny)  ! cover grid edges
      iy1 = min(max(iy+1,1),ny)  ! cover grid edges
      cx(1) = ix0; cy(1) = iy0 
      cx(2) = ix0; cy(2) = iy1 
      cx(3) = ix1; cy(3) = iy0 
      cx(4) = ix1; cy(4) = iy1 
      hweight(1) = (1.0-sx)*(1.0-sy)
      hweight(2) = (1.0-sx)*(  sy  )
      hweight(3) = (  sx  )*(1.0-sy)
      hweight(4) = (  sx  )*(  s  y)

      z          = xyz(3)  ! short hand
      status     = 0       ! assume interpolation possible
      valid      = .false. ! default for dry/out-of-bounds pillars
      fz         = 0.      ! default for dry/out-of-bounds pillars
      dfz        = 0.      ! default for dry/out-of-bounds pillars
c
c     ------ setup 4 interpolation pillars ------
c
      do i = 1,4
         if (wetmask(cx(i),cy(i))>0) then
            valid(i) = .true.
            ibot     = bottom_layer(cx(i),cy(i)) ! ibot >= 1
            zgrid    => ccdepth(cx(i),cy(i),:)
            call search_sorted_list(z,zgrid,iz) ! 0<=iz<=nz+1
            ! ---- capture vertical extrapolation ----
            if ((iz<1).or.(iz>=ibot)) then 
               iz = min(max(iz,1),ibot) ! now iz = 1 or ibot
               flow    = array(cx(i),cy(i),iz)
               fup     = flow    ! => dfz = 0
               sz      = 0.5     ! set dummy
               status  = 2       ! flag vertical extrapolation
               dz(i)   = 1.0     ! avoid numerical problems, assign dummy
            ! ---- interior linear vertical interpolation ----
            else
               flow = array(cx(i),cy(i),iz)
               fup  = array(cx(i),cy(i),iz+1)
               dz(i)= zgrid(iz+1)-zgrid(iz) ! assumed > 0
               sz   = (z-zgrid(iz))/dz(i)
            endif
            fz(i)   = (1.0-sz)*flow + sz*fup
            dfz(i)  = fup-flow
         endif ! wetmask ...
      enddo
c
c     ------ fill in missing data where flagged ------
c
      if (all(valid)) then

         continue

      elseif (.not.any(valid)) then

         result = padval
         if (deriv>0) result = 0 
         status = 3  ! rank deficit exit
         return

      else ! we know at least one valid corner

         fbar = sum(fz)/count(valid)
         do i = 1,4
            if (.not.valid(i)) fz(i) = fbar ! relax to average
         enddo

      endif
c
c     ------ evaluate interpolations/derivatives ------
c
      result = 0
      if (deriv == 0) then     ! evaluate value
         do i = 1,4
            result = result + hweight(i)*fz(i) 
         enddo
      elseif (deriv == 3) then ! evaluate z derivative
         do i = 1,4
            result = result + hweight(i)*dfz(i)/dz(i) 
         enddo
      else
         write(*,*) "interpolate_cc_3Dgrid_data: deriv = ",
     +               deriv, "is not implemented"  
      endif

      end subroutine interpolate_cc_3Dgrid_data



      subroutine interpolate_cc_2Dgrid_data(xy,array,deriv,result,
     +                                      status)
c     --------------------------------------------------------------------
c     Interpolate grid of corner centered data 2D array on point xy
c     with data points at integer valued xy
c
c     deriv = 0 gives value, deriv = (1,2) gives derivative along (x,y)
c     Currently, only deriv = 0, until other derivatives 
c     are needed.
c
c     Return status:
c       status = 0: interior interpolation performed
c       status = 1: horizontal range violation, set result = padval
c     --------------------------------------------------------------------
      real, intent(in)     :: xy(:),array(:,:)
      integer, intent(in)  :: deriv  ! 0=value; deriv = (1,2) gives derivative along (x,y)
      real, intent(out)    :: result
      integer, intent(out) :: status
c      
      integer           :: ix,iy,ix0,ix1,iy0,iy1
      real              :: vc(4),sx,sy
      real, parameter   :: padval = 0. ! later move to argument list
c     --------------------------------------------------------------------
       if (.not.horizontal_range_check(xy)) then
          result = padval
          status = 1
          return
      endif
c
      call get_horiz_grid_coordinates(xy,ix,iy,sx,sy)
c      
      ix0 = min(max(ix,    1),nx)  ! cover boundary layers
      ix1 = min(max(ix + 1,1),nx)  ! cover boundary layers
c     
      iy0 = min(max(iy,    1),ny)  ! cover boundary layers
      iy1 = min(max(iy + 1,1),ny)  ! cover boundary layers
c
      vc(1) = array(ix0,iy0)
      vc(2) = array(ix0,iy1)
      vc(3) = array(ix1,iy0)
      vc(4) = array(ix1,iy1)

      if     (deriv == 0) then 
         call interp_2Dbox_data(sx,sy,vc,deriv,result)
      else
         stop "interpolate_cc_2Dgrid_data: unhandled deriv request"
      endif

      status = 0                ! signal interior interpolation performed
c     -------------------------------------------------------------------- 
      end subroutine interpolate_cc_2Dgrid_data




      subroutine interpolate_turbulence(xyz, r3, status)
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
      integer              :: statv,stath
c     ------------------------------------------ 
      call interpolate_cc_3Dgrid_data(xyz,hdiffus,0,r3(1),stath)
      call interpolate_cc_3Dgrid_data(xyz,vdiffus,0,r3(3),statv)
      r3(2)   = r3(1) ! horizontal isotropy
      status  = max(statv,stath) ! in the unextected case taht they differ ...
c     ------------------------------------------ 
      end subroutine

 

      subroutine interpolate_turbulence_deriv(xyz, r3, status)
c     ------------------------------------------ 
c     Currently do not support horizontal derivatives
c     ------------------------------------------ 
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
      integer              :: statv,stath
c     ------------------------------------------ 
      call interpolate_cc_3Dgrid_data(xyz,vdiffus,3,r3(3),statv)
      r3(1:2) = 0. ! Currently do not support horizontal derivatives
      status  = statv 
c     ------------------------------------------    
      end subroutine 




      subroutine interpolate_currents(xyz, uvw, status)
c     ------------------------------------------
c     Perform linear two-face interpolation in current fields into vector uvw
c     
c     xyz is the interpolation position as (lon[degE],lat[degN],depth[m,positive down])
c     uvw are interpolated currents in m/s, shape == real(3+)
c     status is the status of the interpolation action
c       status = 0: all OK
c       status = 1: domain violation (but extrapolation provided)
c     
c     ------------- data layout -----------------
c
c     u,v positive along lambda,phi in units m/s
c     w   positive downward         in units m/s
c    
c     u(ix,iy,iz) in grid position (ix+0.5, iy    , iz)     (i.e. eastern  cell face)
c     v(ix,iy,iz) in grid position (ix    , iy-0.5, iz)     (i.e. southern cell face)
c     w(ix,iy,iz) in grid position (ix    , iy,     iz-0.5) (i.e. upper    cell face)
c     The boundary condition:
c           w( ix, iy, bottom_layer(ix,iy)+0.5 ) = 0
c     is implicit 
c
c     Implied interpolation ranges in ncc coordinates:
c           1.5 < x < nx+0.5
c           0.5 < y < ny-0.5
c           0.5 < z < nz+0.5 (due to bottom BC)
c     ------------------------------------------ 
      real, intent(in)     :: xyz(:)  ! (lon,lat,depth)
      real, intent(out)    :: uvw(:)
      integer, intent(out) :: status

      integer              :: statu, statv, statw
      integer              :: ix,iy,iz,jwest,jnorth,jlow,ibot
      real                 :: x,y,z, sx,sy,sz
c     ------------------------------------------ 
c.....transform to continuous node-centered grid coordinates
c     and check ranges for this staggering
      call get_ncc_coordinates(xyz,x,y,z)
      statu = 0
      statv = 0 
      statw = 0
      if ((x<1.5).or.(x>(nx+0.5))) statu = 1 ! flag x range violation
      if ((y<0.5).or.(y>(ny-0.5))) statv = 1 ! flag x range violation
      if ((z<0.5).or.(z>(nz+0.5))) statw = 1 ! flag x range violation
      status  = max(statu, statv, statw) 

c.....determine cell associations of point, constrained to 1 <= i <= n
      ix   = max(1, min(nint(x), nx)) ! u at cell east  face
      iy   = max(1, min(nint(y), ny)) ! v at cell south face
      ibot = bottom_layer(ix,iy) ! 1 <= ibot <= nz
      iz   = max(1, min(nint(z),ibot)) ! w at cell upper face

c.....determine relative intracell coorninate s, constrained to 0 < s < 1 
c     when point exceed coordinate bounds s is assigned 0 or 1
      sx = max(0.d0, min(1.d0, x+0.5d0-real(ix)))  
      sy = max(0.d0, min(1.d0, y+0.5d0-real(iy)))           
      sz = max(0.d0, min(1.d0, z+0.5d0-real(iz)))  

c.....face-to-face interpolation 
c     cell indices satisfy 1 <= (ix,iy,iz) <= (nx,ny,ibot)
      jwest  = max(1, min(nint(x-1), nx)) ! cell west  face
      jnorth = max(1, min(nint(y+1), ny)) ! cell north face
      uvw(1) = sx*u(ix, iy, iz)    + (1.-sx)*u(jwest, iy, iz)
      uvw(2) = sy*v(ix, jnorth,iz) + (1.-sy)*v(ix, iy, iz)

      if (iz==ibot) then
         jlow = 0               ! for debugging
         uvw(3) = (1.-sz)*w(ix, iy, ibot) ! w=0 at sea bed, ibot corresponds to upper cell face
      else
         jlow = max(1, min(nint(z+1),ibot)) ! cell lower face
         uvw(3) = sz*w(ix, iy, jlow) + (1.-sz)*w(ix, iy, iz)
      endif
c     ------------------------------------------ 
      end subroutine interpolate_currents



      subroutine interpolate_temp (xyz, r, status) 
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
c     ------------------------------------------ 
      call interpolate_cc_3Dgrid_data(xyz,temp,0,r,status)
c     ------------------------------------------ 
      end subroutine 



      subroutine interpolate_zooplankton (xyz, r, status) 
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r(:)
      integer, intent(out) :: status
c     ------------------------------------------ 
      call interpolate_cc_3Dgrid_data(xyz,zoo,0,r(1),status)
c     ------------------------------------------ 
      end subroutine 



      subroutine interpolate_wdepth(xy, r, status) 
c     ------------------------------------------ 
c     Multiply wetmask to ensure piecewise linear coastlines
c     defined by wdepth=0 
c     ------------------------------------------ 
      real, intent(in)     :: xy(:) 
      real, intent(out)    :: r
      integer, intent(out) :: status
      integer           :: ixc,iyc
      real              :: xc,yc
c     ------------------------------------------ 
      call interpolate_cc_2Dgrid_data(xy,wdepth,0,r,status)
      call get_horiz_ncc(xy,ixc,iyc,xc,yc)
      r = r*wetmask(ixc,iyc)
c     ------------------------------------------ 
      end subroutine 




      LOGICAL function is_wet(xyz)
c     ------------------------------------------ 
c     Probe wdepth
c     return is_wet = .true. at horizontal range violation
c     ------------------------------------------ 
      real, intent(in) :: xyz(:) 
      real             :: wd
      integer          :: status
c     ------------------------------------------ 
      call interpolate_wdepth(xyz, wd, status) 
      if (status==1) then
         is_wet = .true.
         return
      elseif (status==0) then
         if ((xyz(3)>0).and.(xyz(3)<wd)) then
            is_wet = .true.
         else
            is_wet = .false.
         endif
      else
         stop "is_wet: unhandled return condition"
      endif
c     ------------------------------------------
      end function


 
      LOGICAL function is_land(xy)
c     ------------------------------------------ 
c     return is_land = .false. at horizontal range violation
c     ------------------------------------------ 
      real, intent(in) :: xy(:)
      real             :: wdepth
      integer          :: status
      integer          :: ixc,iyc
      real             :: xc,yc
c     ------------------------------------------ 
c      call interpolate_wdepth(xy, wdepth, status) 
  
      call get_horiz_ncc(xy,ixc,iyc,xc,yc)

c     avoid probing wetmask at range violation

      if (horizontal_range_check(xy)) then ! range OK
         call get_horiz_ncc(xy,ixc,iyc,xc,yc)
         is_land = abs(wetmask(ixc,iyc))<1.e-6
      else
         is_land = .false.
      endif    
c     ------------------------------------------
      end function





      LOGICAL function horizontal_range_check(xy)
c     ------------------------------------------
c     Return .true. if cell association of xy
c     0 < (ix,iy) < (nx,ny) is obeyed
c     upper limit (i<n): require data coverage on both
c     sides of an interior point
c     ------------------------------------------ 
      real, intent(in) :: xy(:)
      integer          :: ix,iy
      real             :: sx,sy
c     ------------------------------------------ 
      call get_horiz_grid_coordinates(xy,ix,iy,sx,sy)
      horizontal_range_check = (0<ix).and.(ix<nx)
     +                    .and.(0<iy).and.(iy<ny)
c     ------------------------------------------ 
      end function



      subroutine get_jacobian(xy,r3)
c     ------------------------------------------
c     Transform from coordinates space to Cartesian
c     dR = r3 * (dlon[deg], dlat[deg], dz[m]) [m**3]
c     ------------------------------------------
      real, intent(in)     :: xy(:) ! only 1:2 used
      real, intent(out)    :: r3(:) ! jacobian
      r3(1) = earth_radius*cos(xy(2)*deg2rad)*deg2rad
      r3(2) = earth_radius*deg2rad
      r3(3) = 1.0
      end subroutine



      subroutine get_horizontal_distance(xy1, xy2, r)
c     ------------------------------------------ 
      real, intent(in)     :: xy1(:), xy2(:)
      real, intent(out)    :: r
      real :: costheta, l1,l2,p1,p2, angle
c     ------------------------------------------ 
      p1 = xy1(1)*deg2rad
      p2 = xy2(1)*deg2rad
      l1 = xy1(2)*deg2rad
      l2 = xy2(2)*deg2rad
      costheta = cos(p1)*cos(p2) + sin(p1)*sin(p2)
      costheta = costheta*cos(l1)*cos(l2) + sin(l1)*sin(l2)
      angle    = acos(costheta) ! 0 < angle < pi
      r        = angle*earth_radius
      end subroutine 


 
      subroutine get_local_distance_scalar(xyz1, xyz2, r)
c     ------------------------------------------ 
      real, intent(in)     :: xyz1(:),xyz2(:)
      real, intent(out)    :: r
      real                 :: v(3)
c     ------------------------------------------ 
      call get_local_distance_vector(xyz1, xyz2, v)
      r = sqrt(sum(v*v))
      end subroutine 


      subroutine get_local_distance_vector(xyz1, xyz2, v)
c     ------------------------------------------
c     vector in meters (oriented from xyz1 to xyz2)
c     ------------------------------------------
      real, intent(in)     :: xyz1(:),xyz2(:)
      real, intent(out)    :: v(:)
      real                 :: xyzmid(3),jac(3)
      xyzmid = 0.5*(xyz1 + xyz2)
      call get_jacobian(xyzmid, jac)
      v = jac*(xyz2 - xyz1) ! element-by-element
      end subroutine 




    
      subroutine coast_line_intersection(xyz1, xyz2, anycross, xyzref, 
     +                                   xyzhit) 
c     ------------------------------------------ 
c     This function is a service function assisting in enforcing 
c     coastal boundary conditions on particle steps. xyz1 is the current 
c     valid (wet, inside domain) position. xyz2 is the next position that may or may not be in water.
c     The assumption is that the particle moves in a straight line 
c     between xyz1 and xyz2 in (lon,lat,z) coordinates.
c     Output: If the straight line xyz1 -> xyz2 has crossed a coast line, 
c     anycross is returned .true. else .false. (i.e. the line between xyz1 and xyz2 is all in water). 
c     If anycross is true, the following vectors are computed:
c     xyzref is the modified final position, if step from xyz1 to xyz2 is reflected
c     in the coast line. xyzhit is the first position where the line from xyz1 to xyz2 
c     crosses a coast line. 
c     If anycross is false xyzref and xyzhit is not computed but assigned an interface
c     dependent dummy value.
c     
c     In this implementation land crossing is flagged, if xyz2 is dry. 
c     If xyz1 is dry, a warning is issued and anycross is returned false
c     This misses potential rare edge cuttings (assuming a straight line of motion).
c     In this implementation the fine scale coast line geometry follows 
c     cell centered grid boxes. 
c     In this implementation does not hand multiple reflection effects
c     xyz1 and xyz2 must have same length (2 or 3) so they are vector addable
c     ------------------------------------------ 
      real, intent(in)     :: xyz1(:),xyz2(:)
      logical, intent(out) :: anycross
      real, intent(out)    :: xyzref(:), xyzhit(:)
c
      real                 :: c00(2), c01(2), c10(2), c11(2)
      real                 :: s, s1, s2, s3, s4
      logical              :: xyz1_is_land, xyz2_is_land
      logical              :: cross, any_solution
      character*1          :: reftype
      real                 :: direct(size(xyz1)), reflect(size(xyz1))
      integer              :: ixc, iyc
      real                 :: drycen(2)
c     ------------------------------------------ 
      xyz1_is_land = is_land(xyz1)
      xyz2_is_land = is_land(xyz2)
      if (xyz1_is_land) then
         anycross = .false.
         write(*,*) "warning: coast_line_intersection: xyz1 dry"
         return
      endif
c  
c     The most general condition is 
c     anycross = xyz1_is_land .XOR. xyz2_is_land
c
      anycross = xyz2_is_land  ! assumes .not.xyz1_is_land

c
c     If the coast line is crossed, compute xyzref, xyzhit
c     0 < s < 1 is the coordinate along the vector (xyz2-xyz1)
c 
      if (anycross) then
c
c        --- define node-centered grid cell containing xyz2 (c00, c01, c11, c10) --- 
c            The underlying assumption here is that (xyz1 -> xyz2) only
c            spans two adjecent cells
c
         call get_horiz_ncc_corners(xyz2, c00,c01,c10,c11)   ! in (lon,lat)    
c
c        --- Determine how the step enters the dry box (c00,c01,c10,c11))
c            Make a clockwise round tour: c00 -> c01 -> c11 -> c10 -> c00
c            determine intersection:          s1     s2     s3     s4
c            pick up the closests solution s = si
c
         any_solution = .false.
         s            = 1.0
         call cross_2Dline_segments(xyz1,xyz2,c00,c01,s1,cross)
         if (cross) then
            s            = min(s,s1)
            reftype      = "v"  ! vertical
            any_solution = .true.
         endif
         call cross_2Dline_segments(xyz1,xyz2,c01,c11,s2,cross)
         if (cross) then
            s            = min(s,s2)
            reftype      = "h"  ! horizontal
            any_solution = .true.
         endif
         call cross_2Dline_segments(xyz1,xyz2,c11,c10,s3,cross)
         if (cross) then
            s            = min(s,s3)
            reftype      = "v"  ! vertical
            any_solution = .true.
         endif
         call cross_2Dline_segments(xyz1,xyz2,c10,c00,s4,cross)
         if (cross) then
            s            = min(s,s4)
            reftype      = "h"  ! horizontal
            any_solution = .true.
         endif

c     --- there should be a solution, if we did end in this if clause
         if (.not.any_solution) 
     +       stop "coast_line_intersection: unexpected condition"
         
         direct  = xyz2-xyz1
         reflect = direct
         if     (reftype=="h") then
            reflect(2) = -reflect(2) ! horizontal box line reflection
         elseif (reftype=="v") then
            reflect(1) = -reflect(1) ! vertical   box line reflection
         else
            stop "coast_line_intersection: unhandled reflection type"
         endif
c        |direct| == |reflect|        
         xyzhit = xyz1   +     s*direct  ! position where coast line is hit
         xyzref = xyzhit + (1-s)*reflect ! position for reflection in coast line 
         
      endif
      end subroutine 

     

      subroutine get_horiz_grid_coordinates(xy,ix,iy,sx,sy)
c     ------------------------------------------------------
c     From the lon/lat vector xy (which may include the z component)
c     determine the cell association (ix,iy) where the node of the cell is the
c     lower left corner point, so that [ix:ix+1, iy:iy+1[ maps to cell (ix,iy)
c     i.e. a corner centered convention.
c     0 < (sx,sy) < 1 is the intra cell cooridnates. 
c     Grid origo (ix,iy) = (x,y) = (1,1) is associated with (lambda1,phi1)
c     Do not flag grid range violations, so that 
c     0 < (ix,iy) <= (nx,ny) is not enforced 
c     (horizontal_range_check uses this function to check grid range violations)
c     ------------------------------------------------------
      real, intent(in)     :: xy(:)
      integer, intent(out) :: ix,iy
      real, intent(out)    :: sx,sy
      real                 :: dx1,dy1
c     ------------------------------------------------------
      dx1 = (xy(1)-lambda1)/dlambda
      dy1 = (xy(2)-phi1)   /dphi
      ix  = 1 + int(dx1)    ! truncate decimals to get lon cell
      iy  = 1 + int(dy1)    ! truncate decimals to get lat cell
      sx  = dx1 - int(dx1)  ! intra cell coordinate 0<sx<1
      sy  = dy1 - int(dy1)  ! intra cell coordinate 0<sy<1
c     ------------------------------------------------------
      end subroutine get_horiz_grid_coordinates


      subroutine get_ncc_coordinates(xyz,x,y,z)
c     ------------------------------------------ 
c     Get continuous node-centered grid coordinates with grid points 
c     on integer values of (x,y,z) from xyz = (lon,lat,depth)
c     water surface is at z = 0.5, bottom at z = bottom_layer(...)+0.5
c     No range check imposed
c     ------------------------------------------ 
      real, intent(in)  :: xyz(:)
      real, intent(out) :: x,y,z

      integer           :: ixc,iyc,izc
      real, pointer     :: this_col(:)
      real              :: layerw
c     ------------------------------------------ 
      x   = 1.0 + (xyz(1)-lambda1)/dlambda
      y   = 1.0 + (xyz(2)-phi1)   /dphi
      ixc = nint(x)     
      iyc = nint(y)
c     --- locate vertical cell ---
      this_col => acc_width(ixc, iyc, :)
      call search_sorted_list(xyz(3), this_col, izc) 
      layerw = this_col(izc+1) -  this_col(izc)
      z      = 1.0 + (xyz(3)-ccdepth(ixc, iyc, izc))/layerw
c     ------------------------------------------ 
      end subroutine get_ncc_coordinates


      subroutine get_horiz_ncc(xyz,ixc,iyc,xc,yc)
c     ------------------------------------------ 
c     Resolve the node-centered grid cell containing xyz
c     (ixc,iyc) is grid indices of the closest node
c     (xc,yc)   is (lon, lat) coordinates of the closest node
c     no range check
c     ------------------------------------------ 
      real, intent(in)     :: xyz(:)
      integer, intent(out) :: ixc,iyc
      real, intent(out)    :: xc,yc
c     ------------------------------------------ 
      ixc = 1 + nint((xyz(1)-lambda1)/dlambda)     
      iyc = 1 + nint((xyz(2)-phi1)   /dphi)   
      xc  = lambda1 + (ixc-1)*dlambda
      yc  = phi1    + (iyc-1)*dphi
      end subroutine


      subroutine get_horiz_ncc_corners(xyz,c00,c01,c10,c11)
c     ------------------------------------------ 
c     Resolve corners (c00,c01,c10,c11) in (lon, lat) of the 
c     node-centered grid cell containing xyz
c     ------------------------------------------ 
      real, intent(in)  :: xyz(:)
      real, intent(out) :: c00(2),c01(2),c10(2),c11(2)
      integer           :: ixc,iyc
      real              :: xc,yc
c     ------------------------------------------ 
      call get_horiz_ncc(xyz,ixc,iyc,xc,yc)
      c00(1) = xc - 0.5*dlambda
      c01(1) = xc - 0.5*dlambda
      c10(1) = xc + 0.5*dlambda
      c11(1) = xc + 0.5*dlambda      
      c00(2) = yc - 0.5*dphi
      c01(2) = yc + 0.5*dphi
      c10(2) = yc - 0.5*dphi
      c11(2) = yc + 0.5*dphi
c     ------------------------------------------ 
      end subroutine



      subroutine d_cart2d_xy(xy, r2)
c     ------------------------------------------ 
      real, intent(in)     :: xy(:)
      real, intent(inout)  :: r2(:)
      real                 :: jac(3)
c     ------------------------------------------ 
      call get_jacobian(xy, jac)
      r2(1:2) = r2(1:2)/jac(1:2) ! element-by-element
      end subroutine 



      subroutine d_xy2d_cart(xy, r2) 
c     ------------------------------------------ 
      real, intent(in)     :: xy(:)
      real, intent(inout)  :: r2(:)
      real                 :: jac(3)
c     ------------------------------------------ 
      call get_jacobian(xy, jac)
      r2(1:2) = r2(1:2)*jac(1:2) ! element-by-element
      end subroutine 




      subroutine add_finite_step(xy, dR) 
c     ------------------------------------------ 
c     TODO: implementation deferred, use tmp tangent space 
c     implementation
c     ------------------------------------------ 
      real, intent(inout) :: xy(:)
      real, intent(in)    :: dR(:)
      real :: dxy(3) 
c     ------------------------------------------ 
      if (size(dR) == 3) then  
         dxy = dR
      elseif (size(dR) == 2) then 
         dxy(1:2) = dR
         dxy(3)   = 0.
      else
         write(*,*) "add_finite_step: array mismatch"
         stop
      endif
c
      call d_cart2d_xy(xy,dxy)   ! inplace transform of dxy
      xy = xy + dxy
      
      end subroutine 




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

     
