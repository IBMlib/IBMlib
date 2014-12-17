ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     BIMS-ECO pbi for daily averaged data sets
c     ---------------------------------------------------
c     $Rev: 228 $
c     $LastChangedDate: 2011-01-26 02:52:59 +0100 (Wed, 26 Jan 2011) $
c     $LastChangedBy: asch $ 
c
c     Based on PBI POLCOMS+ERSEM/offline
c    
c     grid def in black_sea.grid.nc
c     daily avg physics in restart.DDDD.nc       (where DDDD is days since 1990-01-01)
c     daily avg biology in black_sea.bio.DDDD.nc (where DDDD is days since 1990-01-01)
c
c     POM employs sigma coordinates, which are bottom-following coordinates that map the vertical coordinate from
c     the vertical coordinate from -H < z < eta onto -1 < sigma_POM < 0, where eta is sea level elevation wrt.
c     reference depth (ncvariable h(x,y) -> wdepth0)
c
c     TODO: add parallel load of biogeochemistry (optional load)
c           vertical/horizontal derivatives of diffusivity
c           smagorinsky horizontal diffusivity
c           add minor correction to convert potential temperature to physical temperature
c   
c     NOTES: black_sea.grid.4096.h (variables z and zz) suggests that nz is number of faces vertically
c            so that the number of wet cells is nz-1 vertically. Cell-centered arrays like t (temperature)
c            are padded with an arbitrary value in last element iz=35 (surface cell at iz=1)
c     
c     LOG: initial version without biogeochemistry tested 17 Dec 2014    
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
c      public :: interpolate_zooplankton

      public :: interpolate_rho ! extension of basic physics interface

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
      type(clock)                :: time_offset ! 1990-01-01 00:00:00
      
      
      integer, parameter         :: tag_lenght = 4     ! format == DDDD (days since 1990-01-01)
      character*(*), parameter   :: not_set_c = "none" ! assume name conflict unlikely
      integer, parameter         :: not_set_i = -9999
     
      logical                    :: data_in_buffers     = .false.  ! first time setting
      logical                    :: data_files_are_open = .false.  ! first time setting
      character(len=tag_lenght)  :: cur_tag      ! tag DDDD of current open set
      integer                    :: ncid         ! NetCDF file handler
      character(len=tag_lenght)  :: fixed_tag    ! test mode (one specified frame is recycled, no time variation)
      
      real, parameter            :: rhoref=1025.e0 ! reference water density kg/m3 (pom08.f + correspondance with Bettina)
      real, parameter            :: molecular_diffusivity = 1.e-9    ! unit m2/s
      real                       :: constant_horiz_diffusivity       ! unit m2/s; handle to override provided horiz_diffusivity


c     --- 3D grids ---
      real, allocatable          :: zoobuf(:,:,:)  ! buffer for adding micro+meso zooplankton 
      real, allocatable          :: rho(:,:,:)     ! water density kg/m3

c     --- 2D grids ---    
      real, allocatable          :: wdepth0(:,:)   ! undisturbed water depth, ncvariable h(x,y)

c     --- 1D grids ---
      real, allocatable          :: sigma_cc(:)    ! sigma of cell centre, 0(surf) < sigma < 1(bott) (NB: sign flipped, -1 < sigma_POM < 0)

c     --- units    ---    
 

c     ===================================================
                            contains
c     ===================================================
  
      subroutine init_physical_fields(time)
c     ------------------------------------------
c     Do not trigger data load
c     ------------------------------------------
      type(clock), intent(in),optional :: time
c     ------------------------------------------
      if (present(time)) then
         call set_master_clock(time)
         write(*,*) "set simulation PBI clock = "
         call write_clock(time)
      endif
      write(*,*) trim(get_pbi_version()) 

      call set_clock(time_offset,1990,1,1,0) 

      call read_control_data(simulation_file,"hydroDBpath",hydroDBpath)
      write(*,*) "init_physical_fields: hydrographic database path =", 
     +           trim(adjustl(hydroDBpath)) 
      call reset_frame_handler()

c     ---- test, if a fixed frame should be recycled (fast test mode)
      if (count_tags(simulation_file, "fixed_hydro_frame")/=0) then
         call read_control_data(simulation_file,
     +                          "fixed_hydro_frame", fixed_tag)
         write(*,*) "init_physical_fields: apply fixed hydroframe:", 
     +              fixed_tag
      else
         fixed_tag = not_set_c  ! apply time varying hydroframes
      endif


      call load_grid_desc()  ! invokes init_horiz_grid_transf
      call init_mesh_grid(init_biogeochem = .false.)  ! allocate core arrays, currently no biogeochem
      call init_topography() ! requires mesh_grid arrays are allocated

      allocate( zoobuf(nx, ny, nz) )
      allocate( rho(nx, ny, nz) )

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
 563  format("init_physical_fields: using constant horiz_diffusivity =",
     +       e12.5," m2/s")
 564  format("init_physical_fields: init default horiz_diffusivity   =",
     +       e12.5," m2/s")

c     ---- specials for this data set / version ----

      u_wind_stress = 0.0  ! no wind in this data set
      v_wind_stress = 0.0  ! no wind in this data set
      
      end subroutine init_physical_fields

 

      character*100 function get_pbi_version()  
      get_pbi_version="BINS-ECO offline pbi version: $Rev:  $"
      end function



      subroutine close_physical_fields()
c     ------------------------------------------ 
c     ------------------------------------------    
      call close_mesh_grid()
      call reset_frame_handler()
      call close_horiz_grid_transf()
c     --- local cleanup ---
      deallocate( zoobuf   )
      deallocate( rho   )
      deallocate( wdepth0  )
      deallocate( sigma_cc )
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
      call resolve_corresp_dataset(aclock, tag)
      call update_dataset(tag)
c     ------------------------------------------       
      end subroutine update_physical_fields




      subroutine resolve_corresp_dataset(aclock, tag)
c     ------------------------------------------
c     Resolve tag of data set corresponding to aclock
c
c     tag in BIMS-ECO PBI for daily averaged data sets is DDDD,
c     where DDDD is days since 1990-01-01
c     In test mode, a fixed frame can be recycled by setting
c     fixed_hydro_frame = DDDD in simulation file
c     ------------------------------------------
      type(clock), intent(in)                :: aclock
      character(len=tag_lenght), intent(out) :: tag
      integer                                :: nhours, ndays
c     ------------------------------------------
      if (fixed_tag == not_set_c) then
         call get_period_length_hour(time_offset, aclock, nhours) ! round to nearest hour
         ndays = int(nhours/24.0) 
         write(tag, 427) ndays 
 427     format(i4.4)           !  DDDD
      else
         tag = fixed_tag
         write(*,*) "resolve_corresp_dataset: applying fixed_tag:", tag
      endif

      end subroutine resolve_corresp_dataset




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
         call load_data_frame() 
      endif
       
      end subroutine update_dataset




      
      subroutine load_grid_desc()
c     -------------------------------------------------------
c     Should only be invoked at start
c     grid file = hydroDBpath/black_sea.grid.nc
c
c     assume grid is regular lon-lat oriented W->E, S->N 
c     -------------------------------------------------------
      character*999        :: fname ! assume long enough

      integer              :: dimid,varid,idum, gncid
      real,allocatable     :: gcoor(:,:)
      real                 :: lambda1, dlambda, phi1, dphi ! LOCAL DUMMIES
c     -------------------------------------------------------
      if (data_in_buffers) stop "load_grid_desc:unexpected"
      write(fname,336) trim(adjustl(hydroDBpath))
      call NetCDFcheck( nf90_open(fname, NF90_NOWRITE, gncid) )
 336  format(a,'/black_sea.grid.nc')
      write(*,*) "load_grid_desc: loading grid from", 
     +           trim(adjustl(fname))

c     --- fetch grid dimensions (NB: nz -> nz-1)---
      call NetCDFcheck( nf90_inq_dimid(gncid, "x",    dimid) )
      call NetCDFcheck( nf90_inquire_dimension(gncid, dimid, len=nx) )
      call NetCDFcheck( nf90_inq_dimid(gncid, "y",    dimid) )
      call NetCDFcheck( nf90_inquire_dimension(gncid, dimid, len=ny) )
      call NetCDFcheck( nf90_inq_dimid(gncid, "z",    dimid) )
      call NetCDFcheck( nf90_inquire_dimension(gncid, dimid, len=nz) )
c     --- In BIMS-ECO nz refers to number of faces including surcafe and bottom, in IBMlib number of cells
      nz = nz - 1   ! number of wet cells ( == IBMlib dimension)

c     --- probe grid descriptors from cell centered points east_e, north_e
      allocate(gcoor(nx,ny))
      call NetCDFcheck( nf90_inq_varid(gncid, "east_e", varid) )
      call NetCDFcheck( nf90_get_var(gncid, varid, gcoor))
      lambda1 = gcoor(1,1)
      dlambda = gcoor(2,1)-gcoor(1,1)
      call NetCDFcheck( nf90_inq_varid(gncid, "north_e", varid) )
      call NetCDFcheck( nf90_get_var(gncid, varid, gcoor))
c     --- probe grid dsecriptors from lon/lat arrays ---
      phi1    = gcoor(1,1)
      dphi    = gcoor(1,2)-gcoor(1,1)
      deallocate(gcoor)
     
      write(*,231) nx, ny, nz
      write(*,232) lambda1, dlambda
      write(*,233) phi1,    dphi 
    
      write(*,*) "load_grid_desc: horizontal initialization"
      call init_horiz_grid_transf(lambda1, phi1, dlambda, dphi)
c
c     allocate and load static layer and topographic descriptors
c
      allocate( wdepth0(nx,ny) )        ! undisturbed water depth
      call NetCDFcheck( nf90_inq_varid(gncid, "h", varid) )
      call NetCDFcheck( nf90_get_var(gncid, varid, wdepth0))
      
      allocate( sigma_cc(nz) )          ! sigma of cell centre
      call NetCDFcheck( nf90_inq_varid(gncid, "zz", varid) )
      call NetCDFcheck( nf90_get_var(gncid, varid, sigma_cc,
     +                  start=(/1/), count=(/nz/)) ) ! skip pad value at end
      sigma_cc = -sigma_cc    ! flip sign so: 0(sea surface) < sigma < 1 (bottom)

      write(*,*) "load_grid_desc: sigma_cc =", sigma_cc

      call NetCDFcheck( nf90_close(gncid) )


 229  format(a,a)    
 231  format("load_grid_desc: 3d grid dim (nx,ny,nz) = ", i4,i4,i4)
 232  format("load_grid_desc: lambda1 = ",f12.7,
     +                      " dlambda = ",f12.7," degrees")  
 233  format("load_grid_desc: phi1    = ",f12.7,
     +                      " dphi    = ",f12.7," degrees")  
      
      end subroutine load_grid_desc
      


      subroutine init_topography()
c     ------------------------------------------------------------
c     Remaining grid/topography initialization (that requires
c     mesh_grid arrays are allocated by init_mesh_grid)  
c     Set mesh_grid arrays:
c        wetmask     
c        bottom_layer
c     ------------------------------------------------------------
      character*999       :: fname ! assume long enough
      integer             :: dimid,varid,idum, gncid
      integer             :: ix,iy
      real,allocatable    :: fsm(:,:)
      real,parameter      :: very_deep = 1.0e5
c     ------------------------------------------------------------
      if (data_in_buffers) stop "init_topography:unexpected"
      write(fname,361) trim(adjustl(hydroDBpath))
      call NetCDFcheck( nf90_open(fname, NF90_NOWRITE, gncid) )
 361  format(a,'/black_sea.grid.nc')
      write(*,*) "init_topography: loading topography", 
     +           trim(adjustl(fname))

c     In this data set apparently the wet/dry condition of grid points is static
c     and accessible in black_sea.grid.h. Therefore load float variable fsm (0/1 free surface mask)
c     and set wetmask + bottom_layer correspondingly
c     
      call NetCDFcheck( nf90_inq_varid(gncid, "fsm", varid) )
      allocate( fsm(nx,ny) )  ! auxillary, typed real  
      call NetCDFcheck( nf90_get_var(gncid, varid, fsm))
      
      where(fsm > 0.0)
         wetmask      = 1 ! wet
         bottom_layer = nz
      elsewhere
         wetmask      = 0 ! dry
         bottom_layer = 0
         wdepth0      = 0 ! dry points in h apparently set to 1
      end where 
      deallocate( fsm )

      call NetCDFcheck( nf90_close(gncid) )

c     ---- render vertical grid tables well defined in all points, including dry
c          ccdepth/acc_width in wet points will be updated dynamically, 
c          when hydrographical frames are loaded
      
      ccdepth(:,:,1)   = very_deep 
      ccdepth(:,:,2:)  = 2*very_deep
      acc_width(:,:,1) = 0.0        ! including surface layer (iz=1)
      acc_width(:,:,2) = 2*very_deep
      acc_width(:,:,3:)= 2*very_deep


      end subroutine init_topography




      subroutine reset_frame_handler()
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
         call NetCDFcheck( nf90_close(ncid) )
         ncid = not_set_i          
         write(*,*) "close_data_files: closed open NetCDF set"
         data_files_are_open = .false.
      endif
      end subroutine close_data_files




      subroutine open_data_files(tag)
c     ------------------------------------------
c     Open netcdf file set corresponding to tag == DDDD: 
c     file = hydroDBpath/restart.DDDD.nc
c     (currently only load of physics)
c     ------------------------------------------
      character(len=tag_lenght),intent(in) :: tag
      character*999                        :: fname ! assume long enough
      integer :: timedimID, ntime, varid
c     ------------------------------------------
      call close_data_files() ! in case they are open ...
      write(fname,333) trim(adjustl(hydroDBpath)), tag

      call NetCDFcheck( nf90_open(fname, NF90_NOWRITE, ncid) )
 333  format(a,'/restart.',a,'.nc')
 
      write(*,*) "open_data_files: opened NetCDF set ", 
     +           trim(adjustl(fname))
c  
      data_files_are_open = .true.
      cur_tag = tag

      end subroutine open_data_files




      subroutine load_data_frame()
c     ------------------------------------------------------------------------------------
c     Verbose reading af NetCDF sets to grid buffers 
c     
c     This data set has only one frame per file; syncronize auxillary grid fields after data load. 
c     It is assumed that grids do not change during a simulation (assertion not checked)
c     
c     POM employs sigma coordinates, which are bottom-following coordinates that map the vertical coordinate from
c     the vertical coordinate from -H < z < eta onto -1 < sigma_POM < 0, where z is -depth below ref level (not sea surface)
c     H > 0 is the static variable h(x,y) (in grid file) loaded into wdepth0
c     eta is dynamic variable elb(x,y) (in physics file) )

c     Currently load these fields (in native units/conventions, c-declaration index order):
c
c     float elb(y, x)  = "surface elevation in external mode at -dt"   units = "metre" ;     staggering  = (east_e, north_e)
c     float u(z, y, x) = "x-velocity"                                  units = "metre/sec" ; staggering  = (east_u, north_u, zz) 
c     float v(z, y, x) = "y-velocity"                                  units = "metre/sec" ; staggering  = (east_v, north_v, zz) 
c     float w(z, y, x) = "sigma-velocity"                              units = "metre/sec" ; staggering  = (east_e, north_e, zz) ... or z ??
c     float t(z, y, x) = "potential temperature" ;                     units = "K" ;         staggering  = (east_e, north_e, zz) 
c     float s(z, y, x) = "salinity x rho / rhoref" ;                   units = "PSS" ;       staggering  = (east_e, north_e, zz) 
c     float kh(z, y, x) = "vertical diffusivity" ;                     units = "m^2/sec";    staggering  = (east_e, north_e, zz) 
c
c     fortran index order opposite CDL dump order
c     fortran: fastests index left most
c     ------------------------------------------------------------------------------------    
      integer :: varid
      integer :: ix,iy,iz,no_fill
      real    :: dz
      logical :: not_ok
      real,pointer  :: cc(:), acc(:)
      integer :: count3D(3)
c     ------------------------------------------------------------------------------------ 
      write(*,*) "load_data_frame: loading frame ", cur_tag


c     -------------------- dynamic topography --------------------

c     ---- load sea level elevation dslm ( <- elb)
      call NetCDFcheck( nf90_inq_varid(ncid, "elb", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, dslm) )
      

c     ------ generate auxillary wet point descriptors
c            ccdepth/acc_width already initialized
c 
      wdepth = 0.0 
      do ix=1,nx
         do iy=1,ny
            if (wetmask(ix,iy)>0) then ! wetmask set in load_grid_desc
               wdepth(ix,iy)    = dslm(ix,iy) + wdepth0(ix,iy)
               ccdepth(ix,iy,:) = sigma_cc*wdepth(ix,iy) ! vector*scalar multiplication
               cc               => ccdepth(ix,iy,:)
               acc              => acc_width(ix,iy,:)
               call levels2accwidth(cc, acc)
            endif      
         enddo
      enddo
c
c     ---- enforce exact consistency between acc_width(:,:,nz+1) and wdepth 
c          (there are numerical diffs of order 10^-5)
c
      acc_width(:,:,nz+1) = wdepth 


c     -------------------- 3D frames ------------------------

      count3D = (/nx,ny,nz/) ! 3D fields are stored with nz+1 elements vertically
      
c     ---- load u current component ---- 
      call NetCDFcheck( nf90_inq_varid(ncid, "u", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, u, count=count3D) )
      

c     ---- load v current component ---- 
      call NetCDFcheck( nf90_inq_varid(ncid, "v", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, v, count=count3D) )

c     ---- load w current component ---- 
      call NetCDFcheck( nf90_inq_varid(ncid, "w", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, w, count=count3D) )
      w = -w  ! w positive direction downward 


c     ---- load temperature (ignore difference beteen potential and physical temperature)----
c          unit seems to be Celcius, despite netcdf attribute t:units = "K" ;

      call NetCDFcheck( nf90_inq_varid(ncid, "t", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, temp, count=count3D) )
     
c     ---- load density ---- 
c     netcdf variable = density-1000)/rhoref
      call NetCDFcheck( nf90_inq_varid(ncid, "rho", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, rho, count=count3D) )
      rho = rho*rhoref + 1000.0 ! convert to physical density in kg/m3

c     ---- load salinity ---- 
c     netcdf variable = salinity * rho / rhoref
      call NetCDFcheck( nf90_inq_varid(ncid, "s", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, salinity, 
     +                  count=count3D) )
      where (rho > 1e-8)
         salinity = salinity*rhoref/rho ! convert to salinity in PSU 
      elsewhere
         salinity = 0.0 ! pad to avoid numerical exceptions
      end where  
               
c     ---- load vertical turbulent diffusivity ----
c          notice that generic vertical allocation is nz+1, because
c          some data sets define data points at vertical faces
      call NetCDFcheck( nf90_inq_varid(ncid, "kh", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, vdiffus)) 
      vdiffus(:,:,nz+1) = 1.e-9 ! set sub surface elements
      

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
      where (hdiffus<1.e-9)
         hdiffus = 1.e-9  
      end where     

c     ------------ update state handlers ------------
      data_in_buffers = .true.


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



      subroutine interpolate_rho(geo, r, status)
c     ------------------------------------------ 
c     extension of minimal physical suite
c     rho is available with this data set
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
c     ------------------------------------------ 
      call interpolate_cc_3Dgrid_data(geo, rho, 0, 0.0, r, status)
c     ------------------------------------------ 
      end subroutine interpolate_rho



      subroutine interpolate_currents(xyz, uvw, status)
c     ---------------------------------------------------------------------------
c     Perform linear two-face interpolation in current fields into vector uvw
c     
c     xyz is the interpolation position as (lon[degE],lat[degN],depth[m,positive down])
c     uvw are interpolated currents in m/s, shape == real(3+)
c     status is the status of the interpolation action
c       status = 0: all OK
c       status = 1: domain violation (but extrapolation provided)
c     
c     --------------------------------- data layout ------------------------------
c
c     u,v positive along lambda,phi in units m/s
c     w   positive downward         in units m/s
c    
c     u(ix,iy,iz) in grid position (ix-0.5, iy    , iz)     (i.e. western  cell face)
c     v(ix,iy,iz) in grid position (ix    , iy-0.5, iz)     (i.e. southern cell face)
c     w(ix,iy,iz) in grid position (ix    , iy,     iz-0.5) (i.e. upper    cell face)
c
c     The boundary condition:
c           w( ix, iy, nz+0.5 ) = 0
c     is implicit 
c
c     Implied interpolation ranges in ncc coordinates:
c           0.5 < x < nx-0.5
c           0.5 < y < ny-0.5
c           0.5 < z < nz+0.5 (due to bottom BC)
c
c     this means that the proper interpolation domain is different from the simulation domain
c     definition derived from the grid. Project these rim points onto the
c     interior domain and interpolate currents here and do not flag them as domain violation
c     --------------------------------------------------------------------------
      real, intent(in)     :: xyz(:)  ! (lon,lat,depth)
      real, intent(out)    :: uvw(:)
      integer, intent(out) :: status

      integer              :: statu, statv, statw
      integer              :: ix,iy,iz,ieast,inorth,jlow,ibot
      real                 :: x,y,z, sx,sy,sz
c     ------------------------------------------ 
      if (.not.is_wet(xyz)) then
         uvw    = 0.
         status = 3  ! signal dry point
         return
      endif
      
c.....transform to continuous node-centered grid coordinates
c     and check ranges for this staggering
      call get_grid_coordinates(xyz,x,y,z)
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
    
c.....determine cell associations of point, constrained to 1 <= i <= n
      ix   = max(1, min(nint(x), nx)) 
      iy   = max(1, min(nint(y), ny)) 
      ibot = bottom_layer(ix,iy) ! 1 <= ibot <= nz
      iz   = max(1, min(nint(z),ibot)) ! w at cell upper face
      
c.....determine relative intracell coorninate s, constrained to 0 < s < 1 
c     when point exceed coordinate bounds s is assigned 0 or 1
      sx = max(0.d0, min(1.d0, x+0.5d0-real(ix)))  
      sy = max(0.d0, min(1.d0, y+0.5d0-real(iy)))           
      sz = max(0.d0, min(1.d0, z+0.5d0-real(iz)))  

c.....face-to-face interpolation 
c     cell indices satisfy 1 <= (ix,iy,iz) <= (nx,ny,ibot)
      ieast  = max(1, min(ix+1, nx)) ! cell east  face
      inorth = max(1, min(iy+1, ny)) ! cell north face
      uvw(1) = sx*u(ieast, iy, iz) + (1.-sx)*u(ix, iy, iz)
      uvw(2) = sy*v(ix, inorth,iz) + (1.-sy)*v(ix, iy, iz)

      if (iz==ibot) then
         uvw(3) = (1.-sz)*w(ix, iy, ibot) ! w=0 at sea bed, ibot corresponds to upper cell face
      else
         jlow = max(1, min(iz+1,ibot)) ! cell lower face
         uvw(3) = sz*w(ix, iy, jlow) + (1.-sz)*w(ix, iy, iz)
      endif
c     ------------------------------------------ 
      end subroutine interpolate_currents


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

      
