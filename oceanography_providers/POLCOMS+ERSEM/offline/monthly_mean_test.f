ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     POLCOMS + ERSEM pbi 
c     ---------------------------------------------------
c     $Rev: 228 $
c     $LastChangedDate: 2011-01-26 02:52:59 +0100 (Wed, 26 Jan 2011) $
c     $LastChangedBy: asch $ 
c
c     Test reading interface for offline POLCOMS+ERSEM data
c     
c     The data set is on a sigma type grid
c
c     Coast line reconstruction: this module features
c     a consistent definition of coast_line_intersection()
c     interpolate_wdepth(), and is_land() so that coast line is
c     characterized by wdepth=0 piecewise lines
c
c     Uses site installation of NetCDF
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module physical_fields
      use time_tools           ! import clock type
      use geometry
      use constants
      use run_context, only: simulation_file
      use input_parser
      use netcdf   ! use site installation

      implicit none
      private     


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
      public :: interpolate_wdepth
      public :: is_wet    
      public :: is_land
      public :: horizontal_range_check
      public :: coast_line_intersection
      public :: get_pbi_version

c     -------------------- module data --------------------  
      
      character*999 :: hydroDBpath ! hydrographic data sets

      type(clock), private, target :: master_clock

c     full grid dimensions: 
      integer :: nxfull,nyfull,nzfull

c     sub grid dimensions:  
      integer :: nx,ny,nz
      real    :: lambda1, dlambda 
      real    :: phi1, dphi     

c     sub grid offsets:     
      integer :: nxoffs,nyoffs
c
      logical :: data_in_buffers

c     --- full grid ---

      real,allocatable :: fullgrid4Dbuf(:,:,:,:)
      real,allocatable :: fullgrid3Dbuf(:,:,:) 
      real,allocatable :: fullgrid2Dbuf(:,:)     

c     --- 3D sub grids ---
             
      real,allocatable :: u(:,:,:)          ! u of current [m/s]   
      real,allocatable :: v(:,:,:)          ! v of current [m/s]        
      real,allocatable :: w(:,:,:)          ! w of current [m/s]      
      real,allocatable :: temp(:,:,:)       ! Water temp. [Celcius]
      real,allocatable :: vdiffus(:,:,:)    ! vertical   diffusivity [m**2/s]              
      real,allocatable :: hdiffus(:,:,:)    ! horizontal diffusivity [m**2/s]
      real,pointer     :: ccdepth(:,:,:)    ! cell center depth water below surface [m]; pointer: allow slicing

c     --- 2D sub grids ---
      real,allocatable :: wdepth(:,:)       ! positive below water surface [m] - dry points negative
      real,allocatable :: wetmask(:,:)      ! auxillary to define coastal geometry (1: wet nodes, 0 dry)
    
c     ===================================================
                            contains
c     ===================================================



  
      subroutine init_physical_fields(time) 
c     ------------------------------------------
      type(clock), intent(in),optional :: time
c     ------------------------------------------
      if (present(time)) master_clock = time
      write(*,*) trim(get_pbi_version()) 

      call read_control_data(simulation_file,"hydroDBpath",hydroDBpath)
      write(*,*) "init_physical_fields: hydrographic database path =", 
     +           trim(hydroDBpath)
    
      call read_grid_desc()
      call allocate_module_buffers()
      
      data_in_buffers = .false.

      end subroutine init_physical_fields

      character*100 function get_pbi_version()  
      get_pbi_version =  "POLCOMS + ERSEM pbi version: $Rev: 228 $"
      end function
    

      subroutine read_grid_desc()
c     ---------------------------------------------------
c     Read grid/subgrid definition which is contained in file 
c     given as grid_desc in the simulation file simulation_file
c     No array allocation in this subroutine. The sub grid is
c     coherent with full grid in overlap region.
c     Set the following module dimensions:
c
c     full grid dimensions:        nxfull,nyfull,nzfull
c     sub grid dimensions:         nx,ny,nz
c     sub grid horizontal offsets: nxoffs,nyoffs
c     sub grid coordinate map:     lambda1, dlambda; phi1, dphi
c
c     so that sub grid range in full grid arrays is
c     
c     array(nxoffs+1 : nxoffs+nx, 
c           nyoffs+1 : nyoffs+ny, 
c           1 : nz)
c
c     At some point, the full grid descriptors should be 
c     extracted directly from the self containerd data sets,
c     however this requires that module data allocation is 
c     deferred to first reading of a time frame. 
c     Grids are not allowed to change during a simulation
c     ---------------------------------------------------
      character*999        :: gd_fname
      type(control_file)   :: grctrlf
c
      real                 :: lambdaB, lambdaE
      real                 :: phiB, phiE
      real                 :: sublamB, sublamE
      real                 :: subphiB, subphiE
c     ---------------------------------------------------
      call read_control_data(simulation_file, "grid_desc", gd_fname)
      write(*,*) "read_grid_desc: reading grid description from: ", 
     +           trim(adjustl(gd_fname))

      call open_control_file(gd_fname, grctrlf)
      
c   
c.....read grid dimensions of full grid
c      
      call read_control_data(grctrlf,"lambda_start_fullgrid", lambdaB)
      call read_control_data(grctrlf,"lambda_end_fullgrid", lambdaE)
      call read_control_data(grctrlf,"dlambda_fullgrid", dlambda) ! set module data
      call read_control_data(grctrlf,"phi_start_fullgrid", phiB)
      call read_control_data(grctrlf,"phi_end_fullgrid", phiE)
      call read_control_data(grctrlf,"dphi_fullgrid", dphi)       ! set module data
      call read_control_data(grctrlf,"nz", nzfull)                ! set module data
c     
      nxfull = 1 + nint((lambdaE-lambdaB)/dlambda)
      nyfull = 1 + nint((phiE   -phiB)   /dphi)
      write(*,*) "read_grid_desc: nxfull =  ", nxfull
      write(*,*) "read_grid_desc: nyfull =  ", nyfull
      write(*,*) "read_grid_desc: nzfull =  ", nzfull
c      
c.....compute dimensions of sub grid (input need not be coherent with full grid)
c
      call read_control_data(grctrlf,"lambda_start_subgrid", sublamB) 
      call read_control_data(grctrlf,"lambda_end_subgrid", sublamE)
      call read_control_data(grctrlf,"phi_start_subgrid", subphiB)
      call read_control_data(grctrlf,"phi_end_subgrid", subphiE)
c      
      nx      = 1 + nint((sublamE-sublamB)/dlambda)
      ny      = 1 + nint((subphiE-subphiB)   /dphi)
      nz      = nzfull
      nxoffs  = nint((sublamB-lambdaB)/dlambda)
      nyoffs  = nint((subphiB-phiB)   /dphi)
      lambda1 = lambdaB + nxoffs*dlambda ! set module data
      phi1    = phiB    + nyoffs*dphi    ! set module data

      write(*,*) "read_grid_desc: nx     = ", nx
      write(*,*) "read_grid_desc: ny     = ", ny
      write(*,*) "read_grid_desc: nz     = ", nz
      write(*,*) "read_grid_desc: nxoffs = ", nxoffs
      write(*,*) "read_grid_desc: nyoffs = ", nyoffs
      write(*,*) "coherent lon grid = [",
     +                                lambda1,lambda1+(nx-1)*dlambda,"]"
      write(*,*) "coherent lat grid = [",
     +                                phi1,   phi1+(ny-1)*dphi,      "]"

      call close_control_file(grctrlf)
c     ------------------------------------------       
      end subroutine read_grid_desc

      

      subroutine allocate_module_buffers()
c     ------------------------------------------ 
c     Assume read_grid_desc has been invoked so that 
c     the module data dimensions have been set:
c
c     full grid dimensions:        nxfull,nyfull,nzfull
c     sub grid dimensions:         nx,ny,nz
c
c     Allocate one buffer to load the full grid
c     ------------------------------------------ 
      write(*,*) "allocate_module_buffers: begin"

      allocate( fullgrid4Dbuf(2,nxfull,nyfull,nzfull) ) ! for zbnd
      allocate( fullgrid3Dbuf(nxfull,nyfull,nzfull) )   
      allocate( fullgrid2Dbuf(nxfull,nyfull)        )

c     --- 3D sub grids ---
  
      allocate( u(nx,ny,nz)       )   
      allocate( v(nx,ny,nz)       )   
      allocate( w(nx,ny,nz)       )   
      allocate( temp(nx,ny,nz)    )  
      allocate( vdiffus(nx,ny,nz) )                      
      allocate( hdiffus(nx,ny,nz) ) 
      allocate( ccdepth(nx,ny,nz) ) 

c     --- 2D sub grids ---
      allocate( wdepth(nx,ny)     ) 
      allocate( wetmask(nx,ny)    ) 
      

      write(*,*) "allocate_module_buffers: end"
c     ------------------------------------------ 
      end subroutine allocate_module_buffers



      subroutine close_physical_fields()
c     ------------------------------------------ 
c     ------------------------------------------ 
      if (allocated(fullgrid4Dbuf)) deallocate(fullgrid4Dbuf) 
      if (allocated(fullgrid3Dbuf)) deallocate(fullgrid3Dbuf) 
      if (allocated(fullgrid2Dbuf)) deallocate(fullgrid2Dbuf)  
      if (allocated(vdiffus))       deallocate(vdiffus)                 
      if (allocated(hdiffus))       deallocate(hdiffus)     
      if (allocated(u))             deallocate(u) 
      if (allocated(v))             deallocate(v) 
      if (allocated(w))             deallocate(w) 
      if (allocated(temp))          deallocate(temp) 
      if (allocated(wdepth))        deallocate(wdepth) 
      if (allocated(wetmask))       deallocate(wetmask) 
      data_in_buffers = .false.
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
c     ------------------------------------------  
      if (present(time)) then
         call set_master_clock(time)
      elseif (present(dt)) then
         call add_seconds_to_clock(master_clock, dt)
      endif

      call check_need_for_update(update)
      if (update) call read_physical_fields()

c     ------------------------------------------       
      end subroutine 



      subroutine check_need_for_update(update)
c     ------------------------------------------ 
      logical, intent(out)    :: update
c     ------------------------------------------ 
      update = .not.data_in_buffers
c     ------------------------------------------ 
      end subroutine check_need_for_update


      

      subroutine read_physical_fields()
c     ------------------------------------------
c     Read NetCDF sets consecutively to full grid buffers  
c     and cast to sub grid buffers.
c     The POLCOMS+ERSEM data is on a sigma type grid
c     In this test mode, time argument is ignored until
c     a time sequence protocol is established
c     Neither synchronize module time with 
c     argument time, but set module temporary flag data_in_buffers
c     It is assumed that grids (sub and full) do not
c     change during a simulation (assertion not checked)
c     
c     Currently load these fields:
c
c       float ucurM(time, z, lat, lon); "m/s" 
c       float vcurM(time, z, lat, lon); "m/s" 
c       float wcurM(time, z, lat, lon); "m/s" downward component
c       float ETWM(time, z, lat, lon) ; "degrees Celsius" ;    
c       float nuvM(time, z, lat, lon) ; "Monthly Mean vertical turbulent diffusivity" ;"m^2/s" 
c       zbnd(time, z, lat, lon, n_zfaces) cell bounds (upper,lower) below surface "m"
c       depth(time, z, lat, lon)          depth of cell center below surface      "m"
c     
c     The is currently no horizontal diffusivity in the set (set it to zero)
c
c     fortran index order opposite CDL dump order
c     fortran: fastests index left most
c     ------------------------------------------
      integer :: ncid, varid    ! netCDF ID for the file and data variable.
      integer :: ix,iy
c     ------------------------------------------ 
c
c     When a time ordered data sets appear in multiple files,
c     a subroutine resolve_data_file_name should be invoked here
c     nf90_open appears fragile to long names or beginning/trailing white spaces
c     so make sure to truncate the file name string
c
      call NetCDFcheck( nf90_open(trim(adjustl(hydroDBpath)), 
     +                  NF90_NOWRITE, ncid) )
c     ---- load u current component ----
      call NetCDFcheck( nf90_inq_varid(ncid, "ucurM", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, fullgrid3Dbuf) )
      call sub_3Dgrid_cast(fullgrid3Dbuf, u)
c     ---- load v current component ----
      call NetCDFcheck( nf90_inq_varid(ncid, "vcurM", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, fullgrid3Dbuf) )
      call sub_3Dgrid_cast(fullgrid3Dbuf, v)
c     ---- load w current component ----
      call NetCDFcheck( nf90_inq_varid(ncid, "wcurM", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, fullgrid3Dbuf) )
      call sub_3Dgrid_cast(fullgrid3Dbuf, w)
c     ---- load temperature ----
      call NetCDFcheck( nf90_inq_varid(ncid, "ETWM", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, fullgrid3Dbuf) )
      call sub_3Dgrid_cast(fullgrid3Dbuf, temp)      
c     ---- load vertical turbulent diffusivity ----
      call NetCDFcheck( nf90_inq_varid(ncid, "nuvM", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, fullgrid3Dbuf) )
      call sub_3Dgrid_cast(fullgrid3Dbuf, vdiffus)    
 
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
      call NetCDFcheck( nf90_get_var(ncid, varid, fullgrid3Dbuf) )
      call sub_3Dgrid_cast(fullgrid3Dbuf, ccdepth) 
      ccdepth = -ccdepth  ! flip sign so wdepth is positive in water

c     ---- load vertical cell face depths (only used to defined total water depth) ----
c     ---- grid is sigma type so all vertical cells are wet, if any
c          depths are positive below water, so flip sign of input data, which
c          has opposite convention
      call NetCDFcheck( nf90_inq_varid(ncid, "zbnd", varid) )
      call NetCDFcheck( nf90_get_var(ncid, varid, fullgrid4Dbuf) )
      call sub_2Dgrid_cast(fullgrid4Dbuf(2,:,:,nz), wdepth) 
      do ix=1,nx
         do iy=1,ny
            if (abs(wdepth(ix,iy))<1.e-6) wdepth(ix,iy) = 1.e-3 ! we flip sign below
            if (is_NetCDF_fill_value_real(wdepth(ix,iy)))
     +          wdepth(ix,iy) = 1.e-3                           ! we flip sign below
         enddo
      enddo
      wdepth = -wdepth ! flip sign of whole array so wdepth is positive in water
c     --- define the wetmask auxillary  ---    
      where(wdepth>0)
         wetmask = 1.0 ! wet
      elsewhere
         wetmask = 0.0 ! dry
      end where 
      
c     Close the file, freeing all resources.
      call NetCDFcheck( nf90_close(ncid) )

      data_in_buffers = .true.

      end subroutine read_physical_fields


      subroutine sub_3Dgrid_cast(afull3Dgridbuff,asubgridarray)
c     ------------------------------------------------------ 
c     copy data, because full3Dgridbuff is reused for many fields
c     ------------------------------------------------------ 
      real, intent(in)     :: afull3Dgridbuff(:,:,:)
      real, intent(out)    :: asubgridarray(:,:,:)

      asubgridarray = afull3Dgridbuff(nxoffs+1 : nxoffs+nx, 
     +                              nyoffs+1 : nyoffs+ny, 
     +                              1 : nz)
c     ------------------------------------------------------ 
      end subroutine sub_3Dgrid_cast


      subroutine sub_2Dgrid_cast(afull2Dgridbuff,asubgridarray)
c     ------------------------------------------------------ 
c     copy data, because full2Dgridbuff is reused for many fields
c     ------------------------------------------------------ 
      real, intent(in)     :: afull2Dgridbuff(:,:)
      real, intent(out)    :: asubgridarray(:,:)

      asubgridarray = afull2Dgridbuff(nxoffs+1 : nxoffs+nx, 
     +                              nyoffs+1 : nyoffs+ny)
c     ------------------------------------------------------ 
      end subroutine sub_2Dgrid_cast


      
      subroutine get_horiz_grid_coordinates(xy,ix,iy,sx,sy)
c     ------------------------------------------------------
c     From the lon/lat vector xy (which may include the z component)
c     determine the cell association (ix,iy) and 
c     intra cell  0 < (sx,sy) < 1. The node of the cell is the
c     lower left corner point, so that [ix:ix+1, iy:iy+1[ maps to cell (ix,iy)
c     i.e. a corner centered convention
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



      subroutine get_vertical_octet(xyz,i38,status)
c     ------------------------------------------------------
c     Identify the grid indices i38of the eight corner points 
c     surrounding the position xyz. 
c     Corners are enumerated in the order 
c     (c000,c001, c010,c011, c100,c101, c110,c111) where indices refer
c     to xyz directions.
c     Horizontal resolution obtained from get_horiz_grid_coordinates 
c     and vertical resolution from the grid of cell center vertical 
c     positions ccdepth.
c     Return vertical indices  0 <= i38(3,:) <= nz 
c     i38(3,:) == 0    signals surface layer extrapolation is required
c     i38(3,:) == nz+1 signals bottom  layer extrapolation is required
c     upper/lower cell centers in that particular column
c     
c     Return status:
c       status = 0: OK, interior interpolation possible
c       status = 1: horizontal range violation, i38 set to -1
c     ------------------------------------------------------
      real, intent(in)     :: xyz(:)
      integer, intent(out) :: i38(3,8)  ! (xyz; c000,c001, ... c111)
      integer, intent(out) :: status 
      integer              :: ix,iy,ix0,iy0,ix1,iy1,icol,iz,i
      real                 :: sx,sy,z
      real, pointer        :: zgrid(:)
c     ------------------------------------------------------   
      if (.not.horizontal_range_check(xyz)) then
          i38 = -1
          status = 1
          return
      endif

      call get_horiz_grid_coordinates(xyz,ix,iy,sx,sy)
c
      ix0 = ix
      ix1 = ix + 1
      ix0 = min(max(ix0,1),nx)  ! cover boundary layers
      ix1 = min(max(ix1,1),nx)  ! cover boundary layers
c      
      iy0 = iy
      iy1 = iy + 1
      iy0 = min(max(iy0,1),ny)  ! cover boundary layers
      iy1 = min(max(iy1,1),ny)  ! cover boundary layers
c
c     load horizontal indices in order
c         c000,c001  (col1)  (ix0,iy0)
c         c010,c011  (col2)  (ix0,iy1)
c         c100,c101  (col3)  (ix1,iy0)  
c         c110,c111  (col4)  (ix1,iy1)
c     
      z        = xyz(3) ! short hand
      i38(3,8) = -1     ! signal extrapolation needed (if not set below)
c     ------ col1 ------
      zgrid => ccdepth(ix0,iy0,:)
      call search_sorted_list(z,zgrid,iz) ! 0<=iz<=nz
      i38(1,1:2) = ix0
      i38(2,1:2) = iy0      
      i38(3,1)   = iz   ! possibly i38==0
      i38(3,2)   = iz+1 ! possibly i38==nz+1 
c     ------ col2 ------
      zgrid => ccdepth(ix0,iy1,:)
      call search_sorted_list(z,zgrid,iz) ! 0<=iz<=nz
      i38(1,3:4) = ix0
      i38(2,3:4) = iy1      
      i38(3,3) = iz     ! possibly i38==0
      i38(3,4) = iz+1   ! possibly i38==nz+1 
c     ------ col3 ------
      zgrid => ccdepth(ix0,iy0,:)
      call search_sorted_list(z,zgrid,iz) ! 0<=iz<=nz
      i38(1,5:6) = ix1
      i38(2,5:6) = iy0      
      i38(3,5) = iz     ! possibly i38==0
      i38(3,6) = iz+1   ! possibly i38==nz+1 
c     ------ col4 ------
      zgrid => ccdepth(ix0,iy1,:)
      call search_sorted_list(z,zgrid,iz) ! 0<=iz<=nz
      i38(1,7:8) = ix1
      i38(2,7:8) = iy1      
      i38(3,7) = iz     ! possibly i38==0
      i38(3,8) = iz+1   ! possibly i38==nz+1 
c     ------------------
      status = 0 ! OK, interpolation possible

c     ------------------------------------------ 
      end subroutine get_vertical_octet



      subroutine interpolate_cc_3Dgrid_data(xyz,array,deriv,result,
     +                                      status)
c     -------------------------------------------------------------------- 
c     Interpolate grid of corner centered data 3D array on point xyz
c     with data points at integer valued xyz
c
c     deriv = 0 gives value, deriv = (1,2,3) gives derivative along (x,y,z)
c     Currently, only deriv = (0,3) is implemented, until other derivatives 
c     are needed.
 
c     Return status:
c       status = 0: interior interpolation performed
c       status = 1: horizontal range violation, set result = padval
c       status = 2: vertical   range violation, set result = padval
c     --------------------------------------------------------------------
      real, intent(in)     :: xyz(:),array(:,:,:)
      integer, intent(in)  :: deriv  ! 0=value; (1,2,3) = along (x,y,z)
      real, intent(out)    :: result
      integer, intent(out) :: status
c      
      integer             :: i38(3,8)  ! (xyz; c000,c001, ... c111)
      integer             :: istat, ix,iy,ic,idum1,idum2,idum3,iz
      real                :: z,sx,sy, depth,zc(8),vc(8), dpt2
      real, parameter     :: padval = 0. ! later move to argument list
c     --------------------------------------------------------------------
      call get_vertical_octet(xyz,i38,istat)  ! includes horizontal range check
c
      
      if (istat==1) then
         status = 1            ! signal horizontal range violation
         result = padval 
         return
      elseif (istat == 0) then ! horizontal range OK
         z = xyz(3)
         call interpolate_wdepth(xyz,depth,idum3)
         if ((z<0).or.(z>depth)) then
            status = 2         ! signal vertical range violation
            result = padval
         endif

c        --- pick corner data points -> (zc,vc)  
         do ic = 1,8
            ix = i38(1,ic)
            iy = i38(2,ic)
            iz = i38(3,ic)
c           ---- special case: point in upper half of surface layer             
            if (iz==0) then            
               zc(ic) = 0               ! project on surface
               vc(ic) = array(ix,iy,1)  ! reuse surface layer center value
c           ---- special case: point in lower half of bottom layer                   
            elseif (iz==nz+1) then 
               zc(ic) = wdepth(ix,iy)   ! project on bottom
               vc(ic) = array(ix,iy,nz) ! reuse bottom layer center value
c           ---- normal case: 1 <= iz <= nz
            else
               zc(ic) = ccdepth(ix,iy,iz)
               vc(ic) = array(ix,iy,iz)
            endif
         enddo

c        ---- pass data to actual interpolator ----
         call get_horiz_grid_coordinates(xyz,idum1,idum2,sx,sy) ! get sx,sy
         call interp_irregbox_data(sx,sy,z,zc,vc,deriv,result)         
c
         status = 0  ! signal interior interpolation performed
      else
         write(*,*) "interpolate_cc_3Dgrid_data: " // 
     +        "unhandled return consition of get_vertical_octet"
         stop 
      endif
c     -------------------------------------------------------------------- 
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
      ix0 = ix
      ix1 = ix + 1
      ix0 = min(max(ix0,1),nx)  ! cover boundary layers
      ix1 = min(max(ix1,1),nx)  ! cover boundary layers
c     
      iy0 = iy
      iy1 = iy + 1
      iy0 = min(max(iy0,1),ny)  ! cover boundary layers
      iy1 = min(max(iy1,1),ny)  ! cover boundary layers
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



      subroutine interpolate_currents(xyz, r3, status)
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
      integer              :: statu, statv, statw
c     ------------------------------------------
      call interpolate_cc_3Dgrid_data(xyz,u,0,r3(1),statu)
      call interpolate_cc_3Dgrid_data(xyz,v,0,r3(2),statv)
      call interpolate_cc_3Dgrid_data(xyz,w,0,r3(3),statw)
      status  = max(statu, statv, statw) ! in the unextected case taht they differ ...
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
      character*5          :: reftype
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
c        --- define node-centered grid cell containing xyz2 (c00, c01, c11, c10) --- 
c            The underlying assumption here is that (xyz1 -> xyz2) only
c            spans two adjecent cells
c
         call get_horiz_ncc_corners(xyz2, c00,c01,c10,c11)       
c
c        --- round tour: c00 -> c01 -> c11 -> c10 -> c00
c            intersection:   s1     s2     s3     s4
c            pick up the closests solution s=si
c
         any_solution = .false.
         s            = 1.0
         call cross_2Dline_segments(xyz1,xyz2,c00,c01,s1,cross)
         if (cross) then
            s            = min(s,s1)
            reftype      = "verti"
            any_solution = .true.
         endif
         call cross_2Dline_segments(xyz1,xyz2,c01,c11,s2,cross)
         if (cross) then
            s            = min(s,s2)
            reftype      = "horiz"
            any_solution = .true.
         endif
         call cross_2Dline_segments(xyz1,xyz2,c11,c10,s3,cross)
         if (cross) then
            s            = min(s,s3)
            reftype      = "verti"
            any_solution = .true.
         endif
         call cross_2Dline_segments(xyz1,xyz2,c10,c00,s4,cross)
         if (cross) then
            s            = min(s,s4)
            reftype      = "horiz"
            any_solution = .true.
         endif

c     --- there should be a solution, if we did end in this if clause
         if (.not.any_solution) 
     +       stop "coast_line_intersection: unexpected condition"
         
         direct  = xyz2-xyz1
         reflect = direct
         if     (reftype=="horiz") then
            reflect(2) = -reflect(2) ! horizontal box line reflection
         elseif (reftype=="verti") then
            reflect(1) = -reflect(1) ! vertical   box line reflection
         else
            stop "coast_line_intersection: unhandled reflection type"
         endif
c        |direct| == |reflect|        
         xyzhit = xyz1   +     s*direct  ! position where coast line is hit
         xyzref = xyzhit + (1-s)*reflect ! position for reflection in coast line 
         
      endif
      end subroutine 

      
      subroutine get_horiz_ncc(xyz,ixc,iyc,xc,yc)
c     ------------------------------------------ 
c     Resolve the node-centered grid cell containing xyz
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
c     Resolve corners (c00,c01,c10,c11) of the 
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
