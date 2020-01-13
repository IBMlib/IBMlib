      module regular_lonlat_grid   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This module provide a leight class for 2D/3D regular lonlat grids to allow mulitple grids
c     
c     Offer basic queries grid-related queries
c     Inlined subroutines from horizontal_grid_transformations_lonlatgrid (with grid argument added)
c      
c     Query return codes: see error codes/messages shared betweem grids/data below
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use geometry              ! cross_2Dlines
      use array_tools           ! search_sorted_list
      
      implicit none
      private           ! default scope
      
      type lonlat_grid
c     ----------------------------------------------------------------
c     Amalgation of 2D/3D grids
c     3D grids in 2D contexts will only apply 2D features
c     2D grids can not generally provide in 3D contexts
c     type attributes in public scope to allow clients (e.g. grid data types)
c     to access data
c     ----------------------------------------------------------------
        integer           :: nx, ny, nz             ! nz not accessed in 2D contexts
        real              :: lon1, lat1, dlon, dlat ! 2D/3D 
        logical           :: is3D ! True: 3D grid; False 2D grid
        integer, pointer  :: wetmask(:,:)      ! 2D/3D  1=wet(legal), 0=dry(illegal)
        integer, pointer  :: bottom_layer(:,:) ! 3D
        
c       ------ reference levels for grid (sea surface elevation = 0)
        real, pointer     :: ccdepth0(:,:,:)   ! 3D reference cell center depth water below surface [m] (dslm=0)
        real, pointer     :: acc_width0(:,:,:) ! 3D accumulated water above this undisturbed layer [m] dim=nz+1
        real, pointer     :: wdepth0(:,:)      ! 3D auxillary, for interpolate_wdepth_grid
c     ------ with added dynamic sea surface elevation
        integer(kind=selected_int_kind(16)) :: revtag ! dynamic revision number of topographic array
        real, pointer     :: ccdepth(:,:,:)    ! 3D actual cell center depth water below surface [m] 
        real, pointer     :: acc_width(:,:,:)  ! 3D actual accumulated water above this undisturbed layer [m] dim=nz+1
        real, pointer     :: wdepth(:,:)       ! 3D current depth at cell-center, including dslm [m]. wdepth = 0 at dry points
      end type
      public :: lonlat_grid

     

      public :: init_lonlat_grid
      public :: finalize_lonlat_grid
      
c     ------- inlined from horizontal_grid_transformations_lonlatgrid -------
      
      public :: horizontal_range_check_for_grid  ! grid specific, should not be confused with global version
      public :: is_land_for_grid                 ! grid specific, should not be confused with global version
      public :: is_wet_for_grid                  ! grid specific, should not be confused with global version
      public :: interpolate_wdepth_for_grid      ! grid specific, should not be confused with global version   
      public :: coastline_intersection_for_grid  ! (included code) grid specific, should not be confused with global version   
      
      public :: get_horiz_grid_coordinates
      public :: get_horiz_geo_coordinates
      public :: get_surrounding_box
      public :: get_horiz_ncc_index
      public :: get_grid_coordinates
      public :: update_topographic_arrays
      public :: dump_dry_or_wet_points   ! for debugging
      
      public :: err_reglonlat_success  
      public :: err_reglonlat_horiz_violation  
      public :: err_reglonlat_vertical_violation
      public :: err_reglonlat_dry_point         
      public :: err_reglonlat_runtimeerr       
      public :: err_reglonlat_inconsisty     
      public :: err_msg_reglonlat   

      
      interface get_horiz_grid_coordinates
         module procedure get_horiz_grid_coor_scalar
         module procedure get_horiz_grid_coor_vector
      end interface
      
      interface get_horiz_geo_coordinates
         module procedure get_horiz_geo_coor_scalar
         module procedure get_horiz_geo_coor_vector
      end interface

      
c     -------------------- module data --------------------  
      
      integer, parameter :: verbose = 0    ! debugging output control      
      real, parameter    :: htol   = 1.e-6 ! tolerance for surface/bottom
      real, parameter    :: padval = 0     ! numeric value applied where undefined

c     ------- error codes/messages shared betweem grids/data

      integer, parameter :: err_reglonlat_success            = 0
      integer, parameter :: err_reglonlat_horiz_violation    = 1
      integer, parameter :: err_reglonlat_vertical_violation = 2
      integer, parameter :: err_reglonlat_dry_point          = 3
      integer, parameter :: err_reglonlat_runtimeerr         = 4
      integer, parameter :: err_reglonlat_inconsisty         = 5
      
      character(len=256), parameter :: err_msg_reglonlat(0:5)= (/
     +     "successful query / operation                    ",   !  error 0
     +     "horizontal range violation (result = padval)    ",   !  error 1
     +     "vertical range violation (result = padval)      ",   !  error 2
     +     "dry point / rank deficit situation for operation",   !  error 3
     +     "runtime error (arrays not initialized)          ",   !  error 4
     +     "internal inconsistency                          "/)  !  error 5
      
c     ===================================================
      contains
c     ===================================================

                            
      subroutine init_lonlat_grid(this, nx, ny, 
     +     lon1, lat1, dlon, dlat, nz)
c     -----------------------------------------------------------------
c     Initialize a lonlat_grid instance
c
c     2D/3D initialization is signalled by absence/presense of nz argument      
c     topographic arrays must be set after invoking init_lonlat_data3D by service requester
c     As only array, wetmask is set to 1 by default at exit (all point wet/legal)
c     correspondingly bottom_layer = 1 will not work for both 2D/3D, since data buffers has different ranks
c     -----------------------------------------------------------------     
      type(lonlat_grid), intent(out)   :: this
      integer, intent(in)              :: nx, ny
      real, intent(in)                 :: lon1, lat1, dlon, dlat
      integer, intent(in), optional    :: nz  
c     -----------------------------------------------------------------
      if (present(nz)) then
         this%is3D = .true.
         this%nz   = nz
      else
         this%is3D = .false.
         this%nz   = -1 ! render defined
      endif
c     ---- 2D attributes and topographic arrays ----         
      this%nx   = nx
      this%ny   = ny
      this%lon1 = lon1
      this%lat1 = lat1
      this%dlon = dlon
      this%dlat = dlat
      allocate ( this%wetmask(nx,ny) )
      this%wetmask      = 1     ! default (all points wet/legal), may be overwritten by service requester 
c     ---- 3D topographic arrays ----
      if (this%is3D) then
         this%revtag = 0                           ! 3D: await update of topographic arrays
         allocate ( this%bottom_layer(nx,ny)  )    ! allocated, but unset
         
         allocate ( this%ccdepth0(nx,ny,nz)   )    ! allocated, but unset
         allocate ( this%acc_width0(nx,ny,nz+1) )  ! allocated, but unset
         allocate ( this%wdepth0(nx,ny)       )    ! allocated, but unset
         
         allocate ( this%ccdepth(nx,ny,nz)    )    ! allocated, but unset
         allocate ( this%acc_width(nx,ny,nz+1)  )  ! allocated, but unset
         allocate ( this%wdepth(nx,ny)       )     ! allocated, but unset
      else
         this%revtag = 1                ! 2D: consider ready with default wetmask = 1 
         nullify( this%bottom_layer )   ! render defined
         nullify( this%ccdepth0     )   ! render defined
         nullify( this%acc_width0   )   ! render defined
         nullify( this%wdepth0      )   ! render defined
         nullify( this%ccdepth      )   ! render defined
         nullify( this%acc_width    )   ! render defined
         nullify( this%wdepth       )   ! render defined
      endif
      end subroutine init_lonlat_grid

      
      subroutine finalize_lonlat_grid(this)
c     -----------------------------------------------------------------
c     finalize lonlat_data3D instance
c     -----------------------------------------------------------------     
      type(lonlat_grid), intent(out) :: this
c     -----------------------------------------------------------------
      this%revtag = 0 ! signal reset
      if (associated(this%wetmask))      deallocate(this%wetmask)
      if (associated(this%bottom_layer)) deallocate(this%bottom_layer)
      if (associated(this%ccdepth0))     deallocate(this%ccdepth0)
      if (associated(this%acc_width0))   deallocate(this%acc_width0)
      if (associated(this%ccdepth))      deallocate(this%ccdepth)
      if (associated(this%acc_width))    deallocate(this%acc_width)
      end subroutine finalize_lonlat_grid
                      
      
c=======================================================================      
c     from horizontal_grid_transformations_lonlatgrid 
c=======================================================================

      
      LOGICAL function horizontal_range_check_for_grid(this, geo)
c     ------------------------------------------
c     Query in relation to specific grid this (different from global query)      
c     Require grid coordinates 
c      (0.5 < x < nx+0.5) and (0.5 < y < ny+0.5)              
c     ------------------------------------------
      type(lonlat_grid), intent(in) :: this
      real, intent(in)              :: geo(:)
      integer                       :: ixc,iyc
c     ------------------------------------------
      call get_horiz_ncc_index(this,geo,ixc,iyc)
      horizontal_range_check_for_grid = (1<=ixc).and.(ixc<=this%nx)
     +                             .and.(1<=iyc).and.(iyc<=this%ny)
c     ------------------------------------------ 
      end function horizontal_range_check_for_grid


      
      LOGICAL function is_land_for_grid(this, geo)       
c     ------------------------------------------ 
c     return is_land = .false. at horizontal range violation
c     for this specific grid
c     ------------------------------------------
      type(lonlat_grid), intent(in) :: this
      real, intent(in)              :: geo(:)
      integer                        :: ixc,iyc
c     ------------------------------------------ 
      if (horizontal_range_check_for_grid(this, geo)) then ! range OK
         call get_horiz_ncc_index(this, geo,ixc,iyc)
         is_land_for_grid = (this%wetmask(ixc,iyc)==0)
      else
         is_land_for_grid = .false.
      endif    
c     ------------------------------------------
      end function


      
      LOGICAL function is_wet_for_grid(this, geo)
c     ------------------------------------------ 
c     Probe wdepth
c     return is_wet = .true. at horizontal range violation
c     accept points at sea surface/bottom as wet (htol)
c     ------------------------------------------
      type(lonlat_grid), intent(in) :: this      
      real, intent(in)              :: geo(:) 
      real                          :: wd
      integer                       :: status
c     ------------------------------------------ 
      call interpolate_wdepth_for_grid(this, geo, wd, status) ! can not signal vertical error 

      select case (status)
      case (err_reglonlat_success)  ! interior wet point
         if ((geo(3)>-htol).and.(geo(3)<wd+htol)) then
            is_wet_for_grid = .true.
         else
            is_wet_for_grid = .false.
         endif
      case (err_reglonlat_horiz_violation)  ! status = 1: horizontal range violation, set result = padval
         is_wet_for_grid = .true.
      case (err_reglonlat_dry_point)  ! status = 3: dry horizontal point
         is_wet_for_grid = .false. ! logical correct
      case (err_reglonlat_runtimeerr)  ! status = 4 : auxillary arrays for query not initialized
         write(*,'(a,a)') "is_wet_for_grid:",
     +        err_msg_reglonlat(err_reglonlat_runtimeerr)
         stop 322
      case default ! internal inconsistency, should not end here
         write(*,'(a,a)') "is_wet_for_grid:",
     +        err_msg_reglonlat(err_reglonlat_inconsisty)
         stop 323
      end select
c     ------------------------------------------
      end function


      
      subroutine interpolate_wdepth_for_grid(this, geo, result, status) 
c     -------------------------------------------------- 
c     Do not apply interpolate_cc_2Dgrid_data which uses
c     wet-point interpolation, but unrestricted interpolation (wdepth=0 at dry points)
c     (otherwise problems in coastal regions)
c     Any horizontal point that tests wet should be assigned a water depth
c     error 2 can not be raised since this is a horizontal query
c     --------------------------------------------------
      type(lonlat_grid), intent(in) :: this        
      real, intent(in)              :: geo(:) 
      real, intent(out)             :: result
      integer, intent(out)          :: status   
c
      integer           :: ix,iy,ix0,iy0,ix1,iy1
      real              :: sx,sy,vc(4)    
c     ------------------------------------------ 
      if (.not.horizontal_range_check_for_grid(this, geo)) then
          result = padval
          status = err_reglonlat_horiz_violation
          return
      endif
      if (is_land_for_grid(this, geo)) then
          result = 0.0 ! fixed value for dry points
          status = err_reglonlat_dry_point 
          return
       endif
       if (this%revtag .le. 0) then
          result = 0.0 ! fixed value for dry points
          status = err_reglonlat_runtimeerr   ! uninitialized dynamic topographic arrays
          return
       endif
c     --- delegate normal interior interpolation to interp_2Dbox_data
      call get_surrounding_box(this,geo,ix,iy,sx,sy)
      ix0 = min(max(ix,    1), this%nx)  ! cover boundary layers
      ix1 = min(max(ix + 1,1), this%nx)  ! cover boundary layers
      iy0 = min(max(iy,    1), this%ny)  ! cover boundary layers
      iy1 = min(max(iy + 1,1), this%ny)  ! cover boundary layers
      vc(1:2) = this%wdepth(ix0, iy0:iy1)
      vc(3:4) = this%wdepth(ix1, iy0:iy1)
      call interp_2Dbox_data(sx, sy, vc, 0, result)   ! defined in grid_interpolations.f
      status = err_reglonlat_success
c     ------------------------------------------ 
      end subroutine interpolate_wdepth_for_grid

      
      

      subroutine get_horiz_grid_coor_scalar(this,geo,x,y)
c     -------------------------------------------------------- 
c     Compute continuous horizontal grid coordinates x,y
c     of lon-lat position geo (which may include the z components etc)
c     Mesh points are on integer values of (x,y) and the
c     first mesh point is at (x,y) = (1,1), i.e. fortran offset
c     no range check     
c     --------------------------------------------------------
      type(lonlat_grid), intent(in) :: this
      real, intent(in)   :: geo(:)
      real, intent(out)  :: x,y
c     -------------------------------------------------------- 
      x  = 1.0 + (geo(1)-this%lon1)/this%dlon
      y  = 1.0 + (geo(2)-this%lat1)/this%dlat
      end subroutine get_horiz_grid_coor_scalar

      
      subroutine get_horiz_grid_coor_vector(this,geo,xy)  
c     -------------------------------------------------------- 
c     Compute continuous horizontal grid coordinates xy
c     of lon-lat position geo (which may include the z components etc)
c     Copy geo(3:max) to xy, if both buffers are larger than 2
c     no horizontal range check     
c     --------------------------------------------------------
      type(lonlat_grid), intent(in) :: this
      real, intent(in)   :: geo(:)
      real, intent(out)  :: xy(:)
      integer            :: n
c     -------------------------------------------------------- 
      n = min(size(geo), size(xy))
      if (n>2) xy(3:n) = geo(3:n)  ! copy rest of geo buffer to xy, if any
      call get_horiz_grid_coor_scalar(this,geo,xy(1),xy(2))
      end subroutine get_horiz_grid_coor_vector



      subroutine get_horiz_geo_coor_scalar(this,x,y,geo)
c     -------------------------------------------------------- 
c     Compute lon-lat position geo corresponding to 
c     continuous horizontal grid coordinates (x,y)
c     Buffer geo may contain z value 
c     --------------------------------------------------------
      type(lonlat_grid), intent(in) :: this
      real, intent(in)    :: x,y
      real, intent(out)   :: geo(:)
c     -------------------------------------------------------- 
      geo(1) = this%lon1 + (x-1.0)*this%dlon 
      geo(2) = this%lat1 + (y-1.0)*this%dlat 
      end subroutine get_horiz_geo_coor_scalar

      
      subroutine get_horiz_geo_coor_vector(this,xy,geo)
c     -------------------------------------------------------- 
c     Compute lon-lat position geo corresponding to 
c     continuous horizontal grid coordinates xy
c     Buffer geo may contain z value 
c     Copy xy(3:max) to geo, if both buffers are larger than 2
c     --------------------------------------------------------
      type(lonlat_grid), intent(in) :: this
      real, intent(in)    :: xy(:)
      real, intent(out)   :: geo(:)
      integer             :: n
c     -------------------------------------------------------- 
      n = min(size(geo), size(xy))
      if (n>2) geo(3:n) = xy(3:n) ! copy rest of xy buffer to geo, if any
      call get_horiz_geo_coor_scalar(this,xy(1),xy(2),geo)
      end subroutine get_horiz_geo_coor_vector



      subroutine get_surrounding_box(this,geo,ix,iy,sx,sy) ! NB previously: get_horiz_grid_coordinates
c     ------------------------------------------------------
c     From the lon/lat vector geo (which may include the z component)
c     determine the cell association (ix,iy) where the node of the cell is the
c     lower left corner point, so that [ix:ix+1, iy:iy+1[ maps to cell (ix,iy)
c     i.e. a corner centered convention.
c     0 < (sx,sy) < 1 is the intra cell cooridnates. 
c     Grid origo (ix,iy) = (x,y) = (1,1) is associated with (lambda1,phi1)
c     Do not flag grid range violations, so that 
c     0 < (ix,iy) <= (nx,ny) is not enforced 
c     (horizontal_range_check uses this function to check grid range violations)
c     NB: previously, this was subroutine called get_horiz_grid_coordinates
c     ------------------------------------------------------
      type(lonlat_grid), intent(in) :: this
      real, intent(in)     :: geo(:)
      integer, intent(out) :: ix,iy
      real, intent(out)    :: sx,sy
      real                 :: dx1,dy1
c     ------------------------------------------------------
      dx1 = (geo(1)-this%lon1)/this%dlon
      dy1 = (geo(2)-this%lat1)/this%dlat
      ix  = 1 + int(dx1)    ! truncate decimals to get lon cell
      iy  = 1 + int(dy1)    ! truncate decimals to get lat cell
      sx  = dx1 - int(dx1)  ! intra cell coordinate 0<sx<1
      sy  = dy1 - int(dy1)  ! intra cell coordinate 0<sy<1
c     ------------------------------------------------------
      end subroutine get_surrounding_box



      subroutine get_horiz_ncc_index(this,geo,ixc,iyc) !,xc,yc)
c     ------------------------------------------ 
c     Resolve the node-centered grid cell containing geo
c     (ixc,iyc) is grid indices of the closest node
c     (xc,yc)   is (lon, lat) coordinates of the closest node
c     Notice that (ixc,iyc) is not coinciding with (ix,iy)
c     in get_horiz_grid_coordinates(geo,ix,iy,sx,sy) above
c     no range check, so mapping NOT restricted to grid limits
c     1 <= ixc,iyc <= nx,ny
c     ------------------------------------------
      type(lonlat_grid), intent(in) :: this
      real, intent(in)     :: geo(:)
      integer, intent(out) :: ixc,iyc
c     ------------------------------------------ 
      ixc = 1 + nint((geo(1)-this%lon1)/this%dlon)     
      iyc = 1 + nint((geo(2)-this%lat1)/this%dlat)   
      end subroutine get_horiz_ncc_index


      
      subroutine get_grid_coordinates(this,geo,x,y,z)  ! formerly named get_ncc_coordinates
c     -------------------------------------------------------------------------------- 
c     Get continuous node-centered grid coordinates with grid points 
c     on integer values of (x,y,z) from geo = (lon,lat,depth)
c     water surface is at z = 0.5, sea bed at z = bottum_layer+0.5
c     It is not checked that z is above the sea bed. The inter grid
c     range is 0.0 <= z <= nz+0.5.
c     If vertical range is exceeded the first/last layer, respectively,
c     is used to extrapolate a vertical grid coordinate smoothly
c     (no extrapolation is flagged) 
c     Include intracell interpolation of layer spacings
c     --------------------------------------------------------------------------------
      type(lonlat_grid), intent(in) :: this
      real, intent(in)              :: geo(:)
      real, intent(out)             :: x,y,z
c     --------------------------------------------------------------------------------
      integer           :: ix,iy,ix0,iy0,ix1,iy1,iz
      real, pointer     :: z00(:), z01(:)
      real, pointer     :: z10(:), z11(:)
      real, allocatable :: z0(:), z1(:), acclay(:)
      real              :: layerw,sx,sy,cwd
c     ------------------------------------------
      if (.not.this%is3D) then
         write(*,*) "get_grid_coordinates: grid is not 3D"
         stop 854
      endif
      allocate( z0(this%nz+1)     )
      allocate( z1(this%nz+1)     )
      allocate( acclay(this%nz+1) )
      call get_horiz_grid_coordinates(this,geo,x,y)   ! define x,y
c
c     ----- interpolate local layer spacings from acc_width: acclay -----
c
      call get_surrounding_box(this,geo,ix,iy,sx,sy)
c     
      ix0 = min(max(ix,  1),this%nx)  ! cover grid edges
      ix1 = min(max(ix+1,1),this%nx)  ! cover grid edges
      iy0 = min(max(iy,  1),this%ny)  ! cover grid edges
      iy1 = min(max(iy+1,1),this%ny)  ! cover grid edges
c     --- trilin interpolation  ---
      z00 => this%acc_width(ix0, iy0, :)  ! range = 1:nz+1
      z01 => this%acc_width(ix0, iy1, :)  ! range = 1:nz+1
      z10 => this%acc_width(ix1, iy0, :)  ! range = 1:nz+1
      z11 => this%acc_width(ix1, iy1, :)  ! range = 1:nz+1
      z0  = z00 + sx*(z10 - z00)   ! vector operation (south face)
      z1  = z01 + sx*(z11 - z01)   ! vector operation (north face)
      acclay = z0 + sy*(z1-z0)     ! vector operation (north-south)  
c
c     ----- interpolate z from local layer spacings -----
c
c     search_sorted_list: result = iz (0 <= iz <= nz+1)
c                         acclay(iz) < geo(3) < acclay(iz+1)   (iz<=nz)
c                         acclay(1) = 0
c
      call search_sorted_list(geo(3), acclay, iz)  ! iz is layer numer, where geo(3) belongs
      iz = min(max(iz,1),this%nz) ! capture vertical range excess 
      layerw = acclay(iz+1) -  acclay(iz) 
      z      = (iz - 0.5) + (geo(3) - acclay(iz))/layerw ! provides smooth extrapolation at range excess 
c
      deallocate( z0, z1, acclay  )
      nullify( z00, z01, z10, z11 )
c     ------------------------------------------ 
      end subroutine get_grid_coordinates

      
      
      subroutine update_topographic_arrays(this, ss_elev, revtag)
c     ---------------------------------------------------------
c     Update topographic arrays of grid in relation to dynamic
c     sea surface elevation ss_elev (positive up)
c     Optionally set revision handle to specific (positive) revtag
c     otherwise just increment current value (default)
c     ---------------------------------------------------------      
      type(lonlat_grid), intent(inout) :: this
      real, intent(in)                 :: ss_elev(:,:) ! this%nx,this%ny
      integer(kind=selected_int_kind(16)), optional :: revtag
      
      integer :: ix,iy
      real    :: dz
c     ---------------------------------------------------------  
      do ix=1,this%nx
         do iy=1,this%ny
            if (this%wetmask(ix,iy)==0) cycle ! nothing for dry points
            dz                      = ss_elev(ix,iy)
            this%wdepth(ix,iy)      =max(0.0, this%wdepth0(ix,iy) + dz)   ! NB: flipped sign
            this%acc_width(ix,iy,2:)=this%acc_width0(ix,iy,2:) + dz       ! NB: flipped sign
            this%ccdepth(ix,iy,1)   =this%ccdepth0(ix,iy,1)    + dz*0.5   ! NB: flipped sign
            this%ccdepth(ix,iy,2:)  =this%ccdepth0(ix,iy,2:)   + dz       ! NB: flipped sign
         enddo
      enddo
      if (present(revtag)) then
         if (this%revtag .le. 0) then
            write(*,*) "update_topographic_arrays: received revtag<0:",
     +           this%revtag
            stop 381
         else
            this%revtag = revtag   ! accept user provided revtag
         endif
      else                                            ! default: just increment counter
         if (this%revtag == huge(this%revtag)) then   ! HUGE(X) The largest positive number
            this%revtag = 1                           ! restart counter for increments
         else
            this%revtag = this%revtag + 1             ! 
         endif
      endif   ! present(revtag)
      
      end subroutine update_topographic_arrays


      subroutine dump_dry_or_wet_points(this, iunit, print_dry)
c     ----------------------------------------------------------------------
c     write coordinates of dry/wet points of this grid to logical unit iunit
c     For debugging
c     ----------------------------------------------------------------------      
      type(lonlat_grid), intent(in) :: this   
      integer, intent(in)           :: iunit
      logical, intent(in), optional :: print_dry ! default .true.
      integer :: ix,iy
      real    :: xy(2)
      logical :: do_print_dry
c     ----------------------------------------------------------------------   
      do_print_dry = .true. ! default 
      if (present(print_dry)) do_print_dry = print_dry

      do ix=1,this%nx 
            do iy=1,this%ny
               call get_horiz_geo_coor_scalar(this, 1.0*ix, 1.0*iy, xy)
               if (do_print_dry) then
                  if (this%wetmask(ix,iy)==0) write(iunit,*) xy  ! print dry points
               else
                  if (this%wetmask(ix,iy)==1) write(iunit,*) xy  ! print wet points
               endif
         enddo
      enddo
      
      end subroutine dump_dry_or_wet_points
      

c     ----------------------------------------------      
      include "coastline_intersection_for_grid.f"   ! provide subroutine coastline_intersection_for_grid
c     ----------------------------------------------     
      
      end module

      
      
