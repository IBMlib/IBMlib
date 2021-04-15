ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
c     Support COARDS compliant netcdf I/O for data on regular lon-lat grids
c     ---------------------------------------------------
c     Module manage grid- and netcdf-related issue
c     xy buffers are created and outside this module; this module just write xy frames
c     Currently only handles writing - currently just a single xy or xyt variable is allowed for
c           
c     Based on experimental CLAIM version of Jun 11, 2020
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module write_netcdf_lonlat_data
      use netcdf
     
      implicit none
      private                   ! default visibility

c --------------------------------------------------------
c ---------------- module data section -------------------
c --------------------------------------------------------

      type netcdf_grid_2Ddata_handler ! keep components public
c        --- grid ---
         integer :: nxs, nys     ! number of grid points
         integer :: nt           ! (current) number of time points
         real    :: lon1, lat1   ! SW corner of grid
         real    :: dlon, dlat   ! grid spacings
c        --- netcdf ---
         integer :: ncid         ! set to -999 when closed
         integer :: data_id, time_id
         character*999  :: fname ! file name asspciated with data
         character      :: mode  ! r(read)/w(write) - currently only writing supported
      end type
      public :: netcdf_grid_2Ddata_handler

      interface init_netcdf_dataset
c         module procedure init_netcdf_dataset_read   ! to be implemented
        module procedure init_netcdf_dataset_write
      end interface
      public :: init_netcdf_dataset
      public :: init_netcdf_dataset_write  ! must be added explicitly, not only indirectly through init_netcdf_dataset
      public :: close_netcdf_dataset

      interface write_netcdf_data
        module procedure write_xy             
        module procedure append_xy_timeframe
      end interface
      public :: write_netcdf_data

      public :: make_time_attribute
      public :: get_subgrid_coordinates
      public :: get_subgrid_indices
      public :: get_area_matrix

      real, parameter :: R_earth = 6371.0 ! km  
      real, parameter :: pi      = 3.1415926535897932385d0
      
c ========================================================
                      contains
c ========================================================
     
      subroutine make_time_attribute(timeattr, tunit, year, month,
     +                     day, hour, minute, second)
c     ------------------------------------------------------      
      character(len=*), intent(out) :: timeattr
      character*(*), intent(in)     :: tunit
      integer, intent(in)           :: year, month, day
      integer, intent(in), optional :: hour, minute, second
      integer                       :: ho, mi, se
c     ------------------------------------------------------
      ho = 0
      mi = 0
      se = 0
      if (present(hour))   ho = hour
      if (present(minute)) mi = minute
      if (present(second)) se = second
c      
      write(timeattr, 288) tunit, year, month, day, ho, mi, se
c      
 288  format(a, " since ", i4.4, 2("-",i2.2), " ", 2(i2.2,":"),i2.2)
      end subroutine make_time_attribute


      
      subroutine init_netcdf_dataset_write(nch, fname,  
     +                     nxs, nys, lon1, lat1, dlon, dlat,
     +                     type, toffs_attr, vname, vattr)
c     ---------------------------------------------------------
c     Currently limited to float data variables
c     type is 'xy' or 'xyt'
c     if 'xyt' time is bound to the unlimited dimension
c     
c     ---------------------------------------------------------
      type(netcdf_grid_2Ddata_handler), intent(out) :: nch
      character*(*), intent(in) :: fname
      integer, intent(in)       :: nxs, nys
      real, intent(in)          :: lon1, lat1, dlon, dlat
      character*(*), intent(in) :: type
      character*(*), intent(in), optional :: toffs_attr  ! ignored for type=xy
      character*(*), intent(in), optional :: vname       ! netcdf variable name (default 'data')
      character(len=*), intent(in), optional :: vattr(:,:)  ! variable metadata shape (2,:),  as (name,value) pairs, value only chars
c    
      integer :: nxs_id, nys_id, lon_id, lat_id, tdim_id
      integer :: ix,iy,iattr
      character*256 :: varname
c     ---------------------------------------------------------
      nch%nxs   = nxs
      nch%nys   = nys   
      nch%lon1  = lon1
      nch%lat1  = lat1                  
      nch%dlon  = dlon
      nch%dlat  = dlat
      nch%fname = trim(adjustl(fname))
      nch%mode  = 'w'
      if (present(vname)) then
         varname = trim(adjustl(vname))  ! assume not too long
      else
         varname = 'data'                ! default
      endif
c      
c     ---- create spatial dimensions and axes ----
c
      call nfcheck( nf90_create(fname, NF90_CLOBBER, nch%ncid) )
      call nfcheck( nf90_def_dim(nch%ncid, "lon", nxs, nxs_id))
      call nfcheck( nf90_def_dim(nch%ncid, "lat", nys, nys_id))
      call nfcheck( nf90_def_var(nch%ncid, "lon", NF90_FLOAT,
     +     (/nxs_id/), lon_id) )
      call nfcheck( nf90_put_att(nch%ncid, lon_id, "long_name",
     +     "longitude") )
      call nfcheck( nf90_put_att(nch%ncid, lon_id, "unit",
     +     "degrees East") )
      call nfcheck( nf90_def_var(nch%ncid, "lat", NF90_FLOAT,
     +     (/nys_id/), lat_id) )
      call nfcheck( nf90_put_att(nch%ncid, lat_id, "long_name",
     +     "latitude") )
      call nfcheck( nf90_put_att(nch%ncid, lat_id, "unit",
     +     "degrees North") )
c      
c     ---- create variable ----
c      
      if (type == 'xyt') then
         nch%nt = 0
         call nfcheck( nf90_def_dim(nch%ncid, "time", NF90_UNLIMITED, 
     +        tdim_id))
         call nfcheck( nf90_def_var(nch%ncid, "time", NF90_INT, 
     +        (/tdim_id/), nch%time_id))
         if (.not.present(toffs_attr)) then
            write(*,562) "time unit mandatory for type=xyt"
            stop
         endif
         call nfcheck( nf90_put_att(nch%ncid, nch%time_id, "unit", 
     +                 toffs_attr)) 
         call nfcheck( nf90_def_var(nch%ncid, varname, NF90_FLOAT,
     +        (/nxs_id, nys_id, tdim_id/), nch%data_id) )
      elseif (type == 'xy') then
         nch%nt = -1
         call nfcheck( nf90_def_var(nch%ncid, varname, NF90_FLOAT,
     +        (/nxs_id,nys_id/), nch%data_id) )
      else
         write(*,563) "unknown type", type
         stop 
      endif
c      
c     ---- add variable attributes, if provided ----
c
c          attribute data is provided as (name1,value1),  (name2,value2), ...
c      
      if ( present(vattr) ) then
         if (size(vattr,1) /= 2) then
            write(*,562) "attribute data must have shape (2,*)"
            stop
         endif
         do iattr = 1, size(vattr,2)
            call nfcheck( nf90_put_att(nch%ncid, nch%data_id,
     +           trim(adjustl(vattr(1,iattr))),
     +           trim(adjustl(vattr(2,iattr)))) )
         enddo
      endif       ! if ( present(vattr) )
      
c      
      call nfcheck( nf90_enddef(nch%ncid) )  ! leave define mode
c      
c     ---- write lon-lat axes ----
c
c     use implied do loops for dumpling lon,lat
c     items of loop must be cast into a vector by (/ ... /)
c      
      call nfcheck( nf90_put_var(nch%ncid, lon_id,
     +     (/ (lon1 + (ix-1)*dlon, ix=1,nxs)  /) ))
      call nfcheck( nf90_put_var(nch%ncid, lat_id,
     +     (/ (lat1 + (iy-1)*dlat, iy=1,nys)  /) ))
      
 562  format("init_netcdf_dataset_write:",a)
 563  format("init_netcdf_dataset_write:",a,":",a)
      
      end subroutine init_netcdf_dataset_write   
      

      
      subroutine close_netcdf_dataset(nch)
      type(netcdf_grid_2Ddata_handler), intent(inout) :: nch
      call nfcheck( nf90_close(nch%ncid) )
      nch%ncid = -999 ! trigger netcdf error if applied
      end subroutine close_netcdf_dataset

      

      subroutine write_xy(nch, buffer)
c     ---------------------------------------------------------
c     write a single time frame buffer at end corresponding
c     to time timeval
c     (added to CLAIM version)
c     ---------------------------------------------------------      
      type(netcdf_grid_2Ddata_handler), intent(inout) :: nch
      real, intent(in) :: buffer(:,:)
c     ---------------------------------------------------------
      call nfcheck( nf90_put_var(nch%ncid, nch%data_id, buffer))
      end subroutine write_xy


      
      subroutine append_xy_timeframe(nch, buffer, timeval)
c     ---------------------------------------------------------
c     write a single time frame buffer at end corresponding
c     to time timeval
c     (subroutine named write_xy_timeframe in CLAIM version)
c     ---------------------------------------------------------      
      type(netcdf_grid_2Ddata_handler), intent(inout) :: nch
      real, intent(in) :: buffer(:,:)
      integer, intent(in) :: timeval
c     ---------------------------------------------------------   
      if (nch%nt < 0) then
         write(*,*) "append_xy_timeframe: not in time append mode"
         stop
      endif
      nch%nt = nch%nt + 1 
      call nfcheck( nf90_put_var(nch%ncid, nch%data_id, buffer,
     +     start = (/1,1,nch%nt/)))
      call nfcheck( nf90_put_var(nch%ncid, nch%time_id, timeval,
     +     start = (/nch%nt/)))
      end subroutine append_xy_timeframe

      
      
      subroutine nfcheck(status)
c     -------------------------------------------------      
      integer, intent ( in) :: status
c     -------------------------------------------------      
      if(status /= nf90_noerr) then 
         write(0,*)  trim(nf90_strerror(status))
         stop  "Stopped"
      end if
      end subroutine nfcheck        


      subroutine get_subgrid_coordinates(nch, pos, sxy, inside)
c     -----------------------------------------------------------
c     inside includes rim [0.5;1] and [nx/y; nx/y+0.5]      
c     -----------------------------------------------------------      
      type(netcdf_grid_2Ddata_handler), intent(in) :: nch
      real, intent(in)     :: pos(:) ! query horizontal position as lon,lat
      real, intent(out)    :: sxy(:) ! grid coordinate corresponding to pos
      logical, intent(out) :: inside ! inside CC cells correspondint to grid
c     -----------------------------------------------------------
      sxy(1) = 1.0 + (pos(1) - nch%lon1) / nch%dlon
      sxy(2) = 1.0 + (pos(2) - nch%lat1) / nch%dlat
      inside = ((sxy(1)>0.5) .and. (sxy(1)<(nch%nxs + 0.5)) .and.
     +          (sxy(2)>0.5) .and. (sxy(2)<(nch%nys + 0.5)))
      end subroutine get_subgrid_coordinates

      

      subroutine get_subgrid_indices(nch, pos, ix, iy, inside)
c     -----------------------------------------------------------
c     Find cell association of position pos; (ix, iy) will only
c     be in valid range 1 <= (ix, iy) <= (nx, ny) if inside is true
c     inside includes rim [0.5;1] and [nx/y; nx/y+0.5]      
c     -----------------------------------------------------------      
      type(netcdf_grid_2Ddata_handler), intent(in) :: nch
      real, intent(in)     :: pos(:)    ! query horizontal position as lon,lat
      integer, intent(out) :: ix,iy  ! resolved cell indices
      logical, intent(out) :: inside ! inside CC cells correspondint to grid
      real                 :: sxy(2) 
c     -----------------------------------------------------------
      call get_subgrid_coordinates(nch, pos, sxy, inside)
      ix = nint(sxy(1))
      iy = nint(sxy(2))
      end subroutine get_subgrid_indices

      
      subroutine get_area_matrix(nch, area_km2)
c     -----------------------------------------------------------
c     Generate the approximative (tangent space) area per grid cell
c     for data normalization
c     area(lat/deg, dlon/deg, dlat/deg) = R*cos(pi*lat/180)*(pi*dlat/180) * R*(pi*dlon/180)
c                                       = (R*pi/180)**2 * cos(pi*lat/180) * dlon * dlat
c     where R is representative Earth radius in km
c     assume buffer area_km2(:,:) is at least (nxs, nys)
c     -----------------------------------------------------------          
      type(netcdf_grid_2Ddata_handler), intent(in) :: nch
      real, intent(out)                            :: area_km2(:,:)  ! nxs, nys
      integer    :: iy
      real       :: fac, lat
c     -----------------------------------------------------------     
      fac = (R_earth*pi/180.)**2 * nch%dlon * nch%dlat
      do iy=1, nch%nys
         lat = nch%lat1 + (iy-1)*nch%dlat
         area_km2(1:nch%nxs, iy) = fac*cos(pi*lat/180)
      enddo
      end subroutine get_area_matrix
      
      end module
