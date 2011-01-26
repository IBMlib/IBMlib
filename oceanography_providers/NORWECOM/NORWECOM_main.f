ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     -----------------------------------------------------------
c     NORWEgian ECOlogical Model (NORWECOM) Oceanography Provider
c     -----------------------------------------------------------
c     $Rev: 231 $
c     $LastChangedDate: 2011-01-26 11:07:27 +0100 (Wed, 26 Jan 2011) $
c     $LastChangedBy: mpay $ 
c
c     This Oceanography Provider provides a general interface to the
c     outputs of the NORWECOM model. The provider is written to be as
c     generic as possible, making space for other alternative outputs
c     from this model source. The default configuration uses the entire
c     frame - however, configuration file options allow for only a 
c     subset of the total domain to be considered.
c
c     The provider has four component files:
c       NORWECOM_main.f          : The master module file
c       NORWECOM_coastline.f     : Deals with coastline interactions
c       NORWECOM_grid.f          : Grid related helper functions
c       NORWECOM_interpolators.f : Functions for interpolating 
c     
c    TODO:--------------------------------
c       * Check calculation of horizontal diffusivity
c       * Check coastline reflection
c    TODO later:--------------------------
c       * Check that velocity interpolation is appropriate
c         especially around the edge
c       * Test velocity vectors by comparison with ARIANE
c       * Improve interpolations around coastline
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module physical_fields
      use run_context, only: simulation_file
      use constants
      use geometry
      use netcdf
      use input_parser
      use time_tools
      implicit none
      private                     ! default visibility

c     -------------------- public interface --------------------  
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
      public :: coast_line_intersection
      public :: horizontal_range_check

c     -------------------- overloaded methods --------------------  

c     -------------------- module data --------------------  
      character*999 :: hydroDBpath ! hydrographic data sets
      character*1, parameter ::  path_separator = "/"   ! UNIX/Linux separator
      real, parameter :: NaN = 0.0 / 0.0
      integer, parameter :: verbose = 0   ! Zero for off

c     Time pieces
c     We are forced to make an ugly hack here, as our time library can't handle times before
c     1970, but the NORWECOM baseline is set to 1 Jan 1950. We therefore need to have a
c     constant integer that accounts for the difference between these two, NORWECOM_t_offset. The value of this integer is the time, in POSIX time, where NORWECOM time units are based
      integer, parameter     ::  NORWECOM_t_offset = 631148400 ! Number of seconds between base_time and NORWECOM time baseline
      type(clock), private, target :: master_clock
      type(clock)  ::  cur_daily_fields, cur_hourly_fields ! Time stamp of ocean frames currently loaded in memory

c     full grid dimensions, including grid specifier
      integer, parameter  ::  nxfull=160  
      integer, parameter  ::  nyfull=140
      integer, parameter  ::  nzfull=20
      real, parameter     ::  DX =10
      real, parameter     ::  DY =10
      real, parameter     ::  XPOLE=382
      real, parameter     ::  YPOLE=256
      real, parameter     ::  ALPHA=58
      real, parameter     ::  lat_true=60
      real  PHINUL

c     sub grid dimensions and offsets
      integer :: nx,ny,nz
      integer :: nxoffs,nyoffs   

c.....grid details
      real,    allocatable   ::  sigma_edges(:)      ! Dividing lines between layers, in sigma coordinates
      real,    allocatable   ::  layer_width(:,:,:)  ! width for layers incl. DSLM [=] m
      real,    allocatable   ::  acc_width(:,:,:)    ! cumulative width for layers above this incl. DSLM, [=] m
      integer, allocatable   :: bottom_layer(:,:) !Index of the very bottom layer

c.....lat/lon grid auxillaries      
      real, allocatable   :: xfacelen(:,:)  ! length of the southern cell face i.e. between cell(i,j-1) and cell(i,j) [m]
      real, allocatable   :: yfacelen(:,:)  ! length of the western  cell face i.e. between cell(i,j) and cell(i+1,j) [m]
      real, allocatable   :: cell_area(:,:) ! horizontal area of cell (i,j) in meters**2

c     --- 3D sub grids (nx,ny,nz) ---
      real,allocatable,target :: u(:,:,:)          ! u of current [m/s]   
      real,allocatable,target :: v(:,:,:)          ! v of current [m/s]        
      real,allocatable,target :: w(:,:,:)          ! w of current [m/s]      
      real,allocatable,target :: temp(:,:,:)       ! Water temp. [Celcius]
      real,allocatable,target :: salinity(:,:,:)      ! Salinity [PSU]
      real,allocatable,target :: vdiffus(:,:,:)    ! vertical diffusivity [m**2/s] 
      real,allocatable,target :: hdiffus(:,:,:)    ! horizontal diffusivity [m**2/s]

c     --- 2D sub grids (nx,ny) ---
      real,allocatable    :: bath(:,:)      ! positive below water surface [m] - dry points negative
      logical,allocatable :: land(:,:)       ! auxillary to define coastal geometry (F: wet nodes, T dry)

c     ===================================================
                         contains
c     ===================================================

c########################################################################
c#                           General Methods                            # 
c########################################################################

      subroutine init_physical_fields(time)
c     ---------------------------------------------------
c     Initialises the physical fields module, including
c     allocating the data arrays
c     ---------------------------------------------------
      character*999        :: bath_fname
      integer              :: ncid, varid,ix,iy,iz,topo_fill
      real                 :: d(3),pos1(3),pos2(3),geopos1(3)
      real                 :: geopos2(3),tmp
      integer,allocatable  :: intbath(:,:)
      type(clock), intent(in),optional :: time
c     ---------------------------------------------------
c.....Display version numbers for reference
      write(*,*) "NORWECOM Oceanography provider versions:"
      write(*,*) "Main comp     : $Rev: 231 $"
      call grid_component_version()
      call interpolator_component_version()
      call coastline_component_version()

c.....Set clock if present     
      if (present(time)) master_clock = time

c.....Read input path      
      call read_control_data(simulation_file,"hydroDBpath",hydroDBpath)
      write(*,*) "init_physical_fields: hydrographic database path =", 
     +           trim(hydroDBpath)
 
c.....Setup subgridding. Default is to use the entire grid, but later
c     versions should be capable of using a subgridding
      nx=nxfull
      ny=nyfull
      nz=nzfull
      nxoffs=1
      nyoffs=1

c     Setup PHINUL, the radian equivalent of lat_trueL
      PHINUL = lat_true * deg2rad        

c.....Read data from the bathmetry file
      allocate(intbath(nx,ny)) !Temporary integer bathymetry 
      allocate(sigma_edges(nz+1))
      call read_control_data(simulation_file,"bath_file", bath_fname) 
      write(*,*) "init_physical_fields: bathymetry file = ",
     +    trim(adjustl(bath_fname))
      call check_ncdf(nf90_open(trim(adjustl(bath_fname)),
     +        nf90_nowrite,ncid))
      call check_ncdf(nf90_inq_varid(ncid,"Topo",varid))
      call check_ncdf(nf90_get_var(ncid,varid,intbath))
      call check_ncdf(nf90_get_att(ncid,varid,"_FillValue",topo_fill))
      call check_ncdf(nf90_inq_varid(ncid,"Zedges",varid))
      call check_ncdf(nf90_get_var(ncid,varid,sigma_edges))
      call check_ncdf(nf90_close(ncid))

c.....setup bathymetry query arrays                
c     The NORWECOM grid is terrain following, so every point is wet, except those on land.
      allocate(bath(nx,ny))    
      allocate(land(nx,ny)) 
      allocate(bottom_layer(nx,ny))
      bottom_layer =nz
      land=.FALSE.
      bath = real(intbath)
      where(intbath==topo_fill )!Bathymetry is fill value for dry land. 
          land=.TRUE.      !Mark as dry 
          bottom_layer=0   !bottom layer is the 0th layer
          bath = NaN         !Set depth to NaN
      end where
      deallocate(intbath)

c.....setup vertical grid auxillary arrays
      allocate(layer_width(nx,ny,nz))   
      allocate(acc_width(nx,ny,nz+1))  
      acc_width(:,:,1)=0.
      do iz=1,nz
         layer_width(:,:,iz)=bath*(sigma_edges(iz+1)-sigma_edges(iz))
         acc_width(:,:,iz+1)=bath*sigma_edges(iz+1)
      enddo 

c.....setup cell physical dimensions
c     xfacelen(i,j) and yfacelen(i,j) give the length of the souuth and
c     west face of cell (i,j) respetively.
      write(*,*) "init_physical_fields: allocating physical dimensions"
      allocate(xfacelen(nx+1,ny+1))   !xfacelen(:,ny+1) are the lengths of the top edge of the grid
      allocate(yfacelen(nx+1,ny+1))   !yfacelen(nx+1,:) are the lengths of the right-hand edge of the grid
      allocate(cell_area(nx,ny)) 
      pos1=0.5             !Define the facelengths to be at the surface
      d=0
      do ix=1,nx+1
         do iy=1,ny+1
c           First the xfacelengths      
            pos1(1)=ix-1   !x-dir left hand side of cell
            pos1(2)=iy-1.0   !xfacelen is defined at the south side of cell
            pos2=pos1
            pos2(1)=ix   !x-dir right hand side of cell
            call xy2lonlat(pos1,geopos1)            
            call xy2lonlat(pos2,geopos2)            
            call get_horizontal_distance(geopos1,
     +               geopos2,xfacelen(ix,iy))
            
c           Then the yfacelengths 
            pos1(1)=ix-1.0   !yfacelen is defined at the west side of cell
            pos1(2)=iy-1.0   !y-dir bottom of cell
            pos2=pos1
            pos2(2)=iy   !y-dir top of cell
            call xy2lonlat(pos1,geopos1)            
            call xy2lonlat(pos2,geopos2)            
            call get_horizontal_distance(geopos1,
     +               geopos2,yfacelen(ix,iy))
          end do
      end do     
c     And now the cell area 
      do ix=1,nx
         do iy=1,ny
            cell_area(ix,iy)=0.5*(xfacelen(ix,iy)+xfacelen(ix,iy+1))*
     +            0.5*(yfacelen(ix,iy)+yfacelen(ix+1,iy))
         end do
      end do     

c.....allocate hydrographical arrays 
      write(*,*) "init_physical_fields: allocating physical data arrays"
      allocate( u(nx,ny,nz) ) 
      allocate( v(nx,ny,nz) ) 
      allocate( w(nx,ny,nz) ) 
      allocate( temp(nx,ny,nz)  )
      allocate( salinity(nx,ny,nz)  )       
      allocate( vdiffus(nx,ny,nz+1))
      allocate( hdiffus(nx,ny,nz))
      vdiffus=0.
      hdiffus=0.

      end subroutine init_physical_fields
      
      
      subroutine close_physical_fields()
c-------------------------------------------------------
c     Cleanup module
c-------------------------------------------------------
      write(*,*) "close_physical_fields: layer widths"
      if (allocated(layer_width))    deallocate(layer_width)
      if (allocated(acc_width))      deallocate(acc_width ) 
      if (allocated(bath))           deallocate(bath)
      if (allocated(bottom_layer))   deallocate(bottom_layer)

      write(*,*) "close_physical_fields: cell sizes"
      if (allocated(xfacelen))       deallocate(xfacelen)
      if (allocated(yfacelen))       deallocate(yfacelen)
      if (allocated(cell_area))      deallocate(cell_area)

      write(*,*) "close_physical_fields: currents"
      if (allocated(u))              deallocate(u) 
      if (allocated(v))              deallocate(v) 
      if (allocated(w))              deallocate(w) 
      write(*,*) "close_physical_fields: physical variables"
      if (allocated(temp))           deallocate(temp)
      if (allocated(salinity))       deallocate(salinity)

      write(*,*) "close_physical_fields: dellocate diffusivities etc"
      if (allocated(hdiffus))          deallocate(hdiffus)
      if (allocated(vdiffus))          deallocate(vdiffus)

      end subroutine close_physical_fields     

      
      function get_master_clock()
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


      subroutine update_physical_fields(force)
c     -----------------------------------------------------
c     Update physical_fields corresponding to match master clock
c     The decision to update  is based upon the distance between the 
c     current time and the current loaded fields - if this is greater 
c     than the time resolution  of the corresponding  fields, an 
c     update is required. An update can be forced by setting the
c     force parameter to TRUE
c     -----------------------------------------------------
      logical,  optional :: force
      integer       :: cyear, cmonth, cday, cisec
      integer       :: offset
      logical       :: updated
      type(clock), pointer :: ctime
      type(timer)   :: stopwatch
c     -----------------------------------------------------
      !Start stopwatch to time update
      call tic(stopwatch) 

      !Set pointer to clock
      ctime => get_master_clock()

c     First if we are forcing the update, lets jump straight there
      if(present(force))  then
      if(force) then
        write(*,*) "Forcing update..."
        call update_daily_physical_fields(ctime)
        call update_hourly_physical_fields(ctime)
        call update_derived_variables()
        call toc(stopwatch,"Physical_fields update complete :") 
        return
      endif
      endif 

c.....Decide whether an update of the daily fields is needed or not
      updated=.FALSE.
      call get_period_length_sec(ctime,cur_daily_fields,offset)
      if(abs(offset) > 43200)  then
        write(*,*) "Updating daily fields...."
        write(*,*) "Offset between daily fields and current time:",
     &       offset,"secs"
        call update_daily_physical_fields(ctime)
        updated = .TRUE.
      endif

c.....Now decide whether an update of the hourly fields is needed or not,
c     using similar logic
      call get_period_length_sec(ctime,cur_hourly_fields,offset)
      if(abs(offset) > 1800 ) then
        write(*,*) "Updating hourly fields...."
        write(*,*) "Offset between hourly fields and current time:",
     &    offset,"secs"
        call update_hourly_physical_fields(ctime)
        updated = .TRUE.
      endif

c.....If fields were updated then update dependent variables      
      if (updated) then
        call update_derived_variables()
        !How long did it all take?
        call toc(stopwatch,"Physical_fields update complete :") 
      endif
      
      end subroutine update_physical_fields


      subroutine update_daily_physical_fields(ctime)
c     -----------------------------------------------------
c     Update physical_fields corresponding to clock time ctime
c     We assume that the NORWECOM daily fields represent the
c     entirity of the given day. Thus, the first entry, which
c     is timestamped with 2006/01/01 12:00 covers everything
c     from 00:00 to 23:59
c     -----------------------------------------------------
      type(clock), intent(in)   :: ctime
      type(clock)   :: ctmp
      integer       :: cyear, cmonth, cday
      integer       :: jday,idx,leng
      real          :: actual_idx(1) !Index that we actually got
c     -----------------------------------------------------
c.....Figure out what the next data frame should be
      call get_date_from_clock(ctime,cyear,cmonth,cday)
      call get_julian_day(ctime,jday)  !The julian day is the index used to store the daily fields
      idx=jday    ! Index used to look up data

c.....Now load the new daily fields
      call read_NORWECOM_fields("STbio_daily","Temp",cyear,idx,temp)
      call land_fill(temp,NaN)  !Set land to NaN
      call read_NORWECOM_fields("STbio_daily","Salt",cyear,idx,salinity)
      call land_fill(salinity,NaN) !Set land to NaN

c.....Get the NORWECOM time that the array index corresponds to      
      call read_NORWECOM_fields("STbio_daily","T",cyear,idx,
     +         dat1d=actual_idx)

c.....And finish off the update by setting the time properly etc and
c     checking it agrees with result from NetCDF files
      call set_clock(ctmp,cyear,cmonth,cday,43200)
      leng = get_POSIXtime(ctmp)
      if (nint(actual_idx(1))*3600.ne.leng+NORWECOM_t_offset) then
          write(*,'(25("-"),a,25("-"))') "STOP"
          write(*,*)"Update_daily_physical_ fields: data frame" //
     +       " loaded does not match that desired."
          write(*,*) "Actual time : ",get_datetime(ctime)
          write(*,*) "Desired time frame : ",get_datetime(ctmp)
          write(*,*) "Index of loaded time frame : ",actual_idx  
          stop
      else
          call set_clock_from_clock(cur_daily_fields,ctmp)
      endif
c     Confirm update
      write(*,*) "New daily fields : ",get_datetime(cur_daily_fields)
 
      end subroutine update_daily_physical_fields


      subroutine update_hourly_physical_fields(ctime)
c     -----------------------------------------------------
c     Update physical_fields corresponding to clock time ctime
c     We assume that the NORWECOM fields are representative of 
c     the hour period in which they lie in the middle of  
c     eg the first field in the NCDF file,
c     which is time stamped 2006/01/01 01:00 represents the flow
c     fields between 00:30 and 01:30.
c     -----------------------------------------------------
      type(clock), intent(in)   :: ctime
      type(clock)   :: ctmp, cstart_of_year
      integer       :: cyear, cmonth, cday,csec
      integer       :: jday, idx, secon
      integer       :: i,j,leng, elapsed
      real          :: actual_idx(1)
c     -----------------------------------------------------
c.....Figure out what the next data frame should be
      call set_clock_from_clock(ctmp,ctime)
      call add_seconds_to_clock(ctmp,-1800) !Taking the clock back some 30min avoids problems at the very start of the year
      call get_date_from_clock(ctmp,cyear,cmonth,cday)
      call set_clock(cstart_of_year,cyear,1,1,0)  
      call get_period_length_sec(ctmp,cstart_of_year,elapsed) !Time elapsed since the start of the year (in secs)
      idx = ceiling(real(elapsed)/3600.)
      if (idx.eq.0) then 
        idx=1
      end if

c.....Now load the new hourly fields
      call read_NORWECOM_fields("velo_hourly","U",cyear,idx,u)
      call land_fill(u,0.0)
      call read_NORWECOM_fields("velo_hourly","V",cyear,idx,v)
      call land_fill(v,0.0)
      call read_NORWECOM_fields("velo_hourly","VertD",cyear,
     +        idx,vdiffus(:,:,1:nz))

c.....Polish up the vertical diffusivity - diffusivity at the boundaries
c     conditions should be zero. Otherwise, diffusivity should not
c     fall below the minimum background value of 1e-5.
      where (vdiffus(:,:,2:nz) < 1e-5)
          vdiffus(:,:,2:nz) = 1e-5
      end where
      vdiffus(:,:,1) =0
      vdiffus(:,:,nz+1) =0
      call land_fill(vdiffus(:,:,1:nz),0.0)

c.....Get the NORWECOM time that the array index corresponds to      
      call read_NORWECOM_fields("velo_hourly","T",cyear,idx,
     +         dat1d=actual_idx)

c.....And finish off the update by setting the time properly etc and
c     checking it agrees with result from NCDF files
      call set_clock_from_clock(ctmp,cstart_of_year)
      call add_seconds_to_clock(ctmp,idx*3600)
      leng= get_POSIXtime(ctmp)
      if (nint(actual_idx(1))*3600.ne.leng+NORWECOM_t_offset) then
          write(*,'(25("-"),a,25("-"))') "STOP"
          write(*,*)"Update_hourly_physical_ fields: data frame" //
     +       " loaded does not match that desired."
          write(*,*) "Actual time : ", get_datetime(ctime)
          write(*,*) "Desired time frame : ",get_datetime(ctmp)
          write(*,*) "Index of loaded time frame : ",actual_idx  
          stop
      else
          call set_clock_from_clock(cur_hourly_fields,ctmp)
      endif
c     Confirm update
      write(*,*) "New hourly fields : ",get_datetime(cur_hourly_fields)

      end subroutine update_hourly_physical_fields


      subroutine read_NORWECOM_fields(prefix,varName,year,idx,
     +                  dat3d,dat2d,dat1d)
c     -----------------------------------------------------
c     Reads data directly from the NORWECOM NetCDF database. The
c     subroutine accepts either a 3d, 2d or 1d output variable
c     which is then used to determine how the data is read
c     -----------------------------------------------------
      character(len=*),intent(in) ::prefix,varName
      integer,intent(in)          ::year,idx
      real,optional,intent(out)   ::dat3d(:,:,:)
      real,optional,intent(out)   ::dat2d(:,:)
      real,optional,intent(out)   ::dat1d(1)
      integer       :: ncid,varid
      character*999 :: ncfilename
      integer       :: status
      real          :: scf
c     -----------------------------------------------------
c     Setup filename
      write(ncfilename,801) trim(adjustl(hydroDBpath)),
     &        path_separator,prefix,year
 801  format(3a,i4,".nc")
      write(*,*) "Loading ",varName," data from: ",trim(ncfilename)
c     Now open the file
      call check_ncdf(nf90_open(ncfilename,nf90_nowrite,ncid))
      call check_ncdf(nf90_inq_varid(ncid,varName,varid))
c     Get scale factor if it exists 
      scf= 1.0
      status = nf90_inquire_attribute(ncid,varid,"scale_factor")
      if (status.eq. nf90_noerr) then
        call check_ncdf(nf90_get_att(ncid,varid,"scale_factor",scf))
      end if
c     Read the data
      if (present(dat3d)) then 
         call check_ncdf(nf90_get_var(ncid,varid,dat3d,
     &            start=(/nxoffs,nyoffs,1,idx/),count=(/nx,ny,nz,1/) )) 
         dat3d=dat3d*scf
      elseif (present(dat2d)) then
         call check_ncdf(nf90_get_var(ncid,varid,dat2d,
     &            start=(/nxoffs,nyoffs,idx/),count=(/nx,ny,1/) )) 
         dat2d=dat2d*scf
      elseif (present(dat1d)) then
         call check_ncdf(nf90_get_var(ncid,varid,dat1d,
     &            start=(/idx/),count=(/1/) )) 
         dat1d=dat1d*scf
      else
         write(*,*) "Arguments to read_NORWECOM_fields misspecified"
         stop
      endif

c     Tidy up
      call check_ncdf(nf90_close(ncid))


      end subroutine read_NORWECOM_fields


      subroutine land_fill(dat,val)
c     -----------------------------------------------------
c     Sets the data stored in the supplied matrix to a 
c     value of "val" at all points that are on land
      real, intent(inout) :: dat(:,:,:)
      real, intent(in)    :: val
      integer ix, iy
c     -----------------------------------------------------
      do ix=1,nx
      do iy=1,ny
       if(land(ix,iy)) dat(ix,iy,:) = val
      enddo
      enddo
      end subroutine land_fill


      subroutine check_ncdf(status)
c     -------------------------------------------------------
c     Checks that NetCDF operations worked
c     -------------------------------------------------------
      integer, intent ( in) :: status
    
      if(status /= nf90_noerr) then 
         print *, trim(nf90_strerror(status))
         call abort("Check_ncdf","Stopped")
      end if
      end subroutine check_ncdf


      subroutine update_derived_variables()
c     -----------------------------------------------------
c     Initiates a recalculation of derived variables the
c     flow fields have been updated
c     -----------------------------------------------------
      write(*,*) "Updating derived variables...."
      call update_w()
      call update_horizontal_turbulence()
      end subroutine update_derived_variables


      subroutine update_w()
c     -----------------------------------------------------
c     Update w(:,:,:) from u(:,:,:), v(:,:,:) using 
c     water incompressibility, starting at each bottom face
c     with w=0 and integrating upward
c     Notice that w(nx,:,:) and w(:,ny,:) can not be defined
c     by this procedure, due to open faces (these planes are 
c     assigned w=0)
c
c     w is midt point at cell upper face and positive upward 
c     so that sea bed boundary condition w=0 is implicit   
c     In principle, w(:,:,1) should equal (d/dt) dslm(:,:)
c     (assertion not checked)
c     -----------------------------------------------------
      integer :: ix, iy, iz
      real    :: dw, wlow, aoh, totwth
c     -----------------------------------------------------      
      w = 0.              ! pad value for undefined entries
      do ix=1,nx-1        ! w(nx,:,:)  undefined 
         do iy=1,ny-1     ! w(:,ny,:) undefined
            wlow = 0.     ! sea bed boundary condition
            do iz = bottom_layer(ix,iy), 1, -1  ! bottom_layer>= 1
               totwth = layer_width(ix,iy,iz)
               aoh = cell_area(ix,iy)/totwth
               dw  = yfacelen(ix,iy)  * u(ix,iy,iz) - 
     1               yfacelen(ix+1,iy)* u(ix+1,iy,iz) + 
     2               xfacelen(ix,iy)  * v(ix,iy,iz) - 
     3               xfacelen(ix,iy+1)* v(ix,iy+1,iz)
               w(ix,iy,iz) = wlow + dw/aoh
               wlow        = w(ix,iy,iz) ! is now next lower face
            enddo
         enddo
      enddo
  
      end subroutine update_w


      subroutine update_horizontal_turbulence()
c     -----------------------------------------------------
c     Smagorinsky parameterization of horizontal eddy diffusivity
c     Parameterisation employed here is based on that employed
c     natively in the NORWECOM code
c  
c     hdiffus is cell centered, i.e. hdiffus(ix,iy,iz) corresponds to (x,y,z) = (ix,iy,iz)
c     Evaluate shear by finite differencing (with one grid spacing from center position)
c     -----------------------------------------------------
c     Smagorinsky parameters for horizontal diffusion      
      real,parameter          :: smagorinsky_constant =0.03  !from M Skogen
c     local variables
      integer                 :: I,J,K, KB
c     ---------------------------------------------------------------------------
      write(*,*) "update_horizontal_turbulence: entered."
c     Setup KB counter to match NORWECOM code
      KB = NZ +1  !Number of sigma  layers +1 
c     Populate the hdiffus array
      DO K = 1,KB-1
        DO J = 2,NY-1
          DO I = 2,NX-1
            hdiffus(I,J,K)=smagorinsky_constant*DX*DY
     1       *SQRT( ((U(I+1,J,K)-U(I,J,K))/DX)**2
     2       +((V(I,J+1,K)-V(I,J,K))/DY)**2
     3       +.5*(.25*(U(I,J+1,K)+U(I+1,J+1,K)-U(I,J-1,K)-U(I+1,J-1,K))
     4        /DY   
     5       + .25*(V(I+1,J,K)+V(I+1,J+1,K)-V(I-1,J,K)-V(I-1,J+1,K))
     6       /DX)**2)
          ENDDO   !I looping in horizontal x
        ENDDO     !J looping in horizotnal y
      ENDDO       !K looping in vertical direction z

      !Deviation from NORWECOM source code
      !Unit conversions. DX and DY in our setup are in km, but we
      !want hdiffus in m^/2. Therefore multiple by 1000, twice
      hdiffus = hdiffus* 1000* 1000
      !/End deviation
C 
C     Extend the diffusivities to the boundary grid cells.
C
      DO  K = 1,KB-1
        DO  J = 2,NY-1
          hdiffus(1,J,K) = hdiffus(2,J,K)
          hdiffus(NX,J,K) = hdiffus(NX-1,J,K)
        ENDDO     !J looping in horizotnal y
      ENDDO       !K looping in vertical direction z
      DO  K = 1,KB-1
        DO  I = 2,NX-1
          hdiffus(I,1,K) = hdiffus(I,2,K)
          hdiffus(I,NY,K) = hdiffus(I,NY-1,K)
        ENDDO     !I looping in horizontal x
      ENDDO       !K looping in vertical direction z
      DO  K = 1,KB-1
        hdiffus(1,1,K) = hdiffus(2,2,K)
        hdiffus(1,NY,K) = hdiffus(2,NY-1,K)
        hdiffus(NX,NY,K) = hdiffus(NX-1,NY-1,K)
        hdiffus(NX,1,K) = hdiffus(NX-1,2,K)
      ENDDO       !K looping in vertical direction z
C
c     Original code has a 2D implementation here. We don't have anything
c     analogous in IBMlib, so ignore. However, have maintained code here
c     for consistency and reference
!!      DO 100 J=1,NY
!!        DO 100 I=1,NX
!!          AAM2D(I,J)=0.
!!100   CONTINUE
!!      DO 110 K = 1,KB-1
!!        DO 110 J = 1,NY
!!          DO 110 I = 1,NX
!!            AAM2D(I,J)=AAM2D(I,J)+AAM(I,J,K)*DZ(K)
!!110   CONTINUE
!!C
      end subroutine update_horizontal_turbulence



c########################################################################
c########################################################################
c#           Import other components                                    # 
c########################################################################
      include 'NORWECOM_grid.f'  !Grid relatied transformations and functions
      include 'NORWECOM_interpolators.f' !Interpolation functions and public accessors
      include 'NORWECOM_coastline.f' !Coastline intersection 

      end module
      

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$      PROGRAM test
c$$$C     ifort -e90 physical_fields.f run_context.o input_parser.o time_tools.o dmidata_interface.o libtime/libtime77.a gribex/libgribex.a
c$$$      use physical_fields
c$$$      integer :: i,j,k, ic,ix,iy,iz
c$$$      integer, parameter :: ns = 7
c$$$      real :: pos(3,ns**3), uvw(3,ns**3),x,y,z
c$$$c      real    :: spts(ns) = (/-0.7,  0.49, 0.51, 1.49,
c$$$c     +                         1.51, 2.49, 2.51, 8.0/)
c$$$  real    :: spts(ns) = (/1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1/)
c$$$
c$$$c     nx=2; ny=2; nz=2
c$$$      nx=4; ny=4; nz=4
c$$$
c$$$      allocate( u(nx,ny,nz)         ) 
c$$$      allocate( v(nx,ny,nz)         ) 
c$$$      allocate( w(nx,ny,nz)         ) 
c$$$      allocate( bottom_layer(nx,ny) ) 
c$$$
c$$$c$$$      u=0.
c$$$c$$$      v=0.
c$$$c$$$      w=0.
c$$$c     u=3x+1  v=4y+2  w=5z+7
c$$$      do ix=1,nx
c$$$         do iy=1,ny
c$$$            do iz=1,nz
c$$$               u(ix,iy,iz) = 3*(ix+0.5) + 1
c$$$               v(ix,iy,iz) = 4*(iy-0.5) + 2
c$$$               w(ix,iy,iz) = 5*(iz-0.5) + 7
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$      bottom_layer=nz
c$$$c      call random_number(pos)
c$$$c      pos(1,:) = pos(1,:)*1.5*nx - 1.0
c$$$c      pos(2,:) = pos(2,:)*1.5*ny - 1.0
c$$$c      pos(3,:) = pos(3,:)*1.5*nz - 1.0
c$$$      ic=1
c$$$      do i=1,ns
c$$$         do j=1,ns
c$$$            do k=1,ns
c$$$               pos(:,ic) = (/spts(i), spts(j), spts(k)/)
c$$$               ic = ic+1
c$$$            enddo
c$$$         enddo
c$$$      enddo
c$$$      
c$$$      call interpolate_currents(pos,uvw)
c$$$
c$$$      do ic=1,ns**3
c$$$         x = pos(1,ic)
c$$$         y = pos(2,ic)
c$$$         z = pos(3,ic)
c$$$         write(*,33) x,y,z, uvw(1,ic)-(3*x+1), 
c$$$     +                     uvw(2,ic)-(4*y+2), uvw(3,ic)-(5*z+7)
c$$$ 33      format(6f8.3)
c$$$      enddo
c$$$
c$$$      end PROGRAM 


