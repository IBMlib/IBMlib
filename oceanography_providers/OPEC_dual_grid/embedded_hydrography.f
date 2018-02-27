      module embedded_hydrography
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Experimental module to represent embedded hydrography (on a regular lonlat grid)
c     within a base grid as defined by e.g. mesh_grid
c
c     Try to keep relation between embedded_hydrography and associating context similar to mesh_grid
c     but inline necessary pieces from  horizontal_representation + horizontal_grid_transformations,
c     as mesh_grid occupies use association (horizontal_representation + horizontal_grid_transformations
c     are singletons). 
c     Currently only represent currents (and topography objects). interpolate_currents
c     assumes a particular staggering
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
      use array_tools                 ! search_sorted

      implicit none
      private           ! default scope

      public :: init_embedded_hydrography   ! init  this module
      public :: close_embedded_hydrography  ! close this module
      
      public :: interpolate_currents    ! corresponding to UV staggering as below
      public :: interpolate_wdepth     
      public :: is_wet    
c      public :: is_land                 
      public :: horizontal_range_check  
     

c     -------------------- module data --------------------  

      integer, parameter :: verbose = 0  ! debugging output control
      real, parameter    :: htol = 1.e-6 ! tolerance for surface/bottom

c
c     ------ grid dimensions/descriptors:   ------
c     
      integer, public        :: nx,ny,nz        
      real                   :: lambda1, dlambda
      real                   :: phi1,    dphi

c     --- 3D grids ---
             
      real,allocatable,public :: u(:,:,:)          ! u of current [m/s] (positive east), staggered on eastern cell faces
      real,allocatable,public :: v(:,:,:)          ! v of current [m/s] (positive north), staggered on southern cell face       
      real,allocatable,public :: w(:,:,:) ! w of current [m/s] (positive down), staggered on upper    cell face
      
      
      real,allocatable,target,public :: ccdepth(:,:,:)    ! cell center depth water below surface [m]; pointer 
      real,allocatable,target,public :: acc_width(:,:,:)  ! accumulated water above this layer [m] dim=nz+1  
   
      real, parameter   :: padval = 0.0      
      
c     --- 2D grids ---
      
      real,allocatable,public     :: wdepth(:,:)       ! current depth at cell-center, including dslm [m]. wdepth = 0 at dry points
      integer, allocatable,public :: bottom_layer(:,:) ! last wet layer (0 for dry points) nx,ny
      integer, allocatable,public :: wetmask(:,:)      ! relocated to this module
      real,allocatable,public     :: dslm(:,:)         ! current sea surface elevation over reference [m] (positive up)

c     ===================================================
                            contains
c     ===================================================

                      
      subroutine init_embedded_hydrography(mx,my,mz,l1,p1,dl,dp)
c     ------------------------------------------------------
c     Initialize module
c     ------------------------------------------------------
      integer, intent(in) ::  mx,my,mz     ! embedded grid dimensions
      real, intent(in)    :: l1,p1,dl,dp   ! embedded grid descriptors
c     ------------------------------------------------------
      write(*,*) "init_embedded_hydrography: "//
     +           "allocate grid arrays: begin" 

c     ---- set module data ----
      nx      = mx
      ny      = my
      nz      = mz
      lambda1 = l1
      dlambda = dl
      phi1    = p1
      dphi    = dp
      
c     ----- physics -----
      allocate( u(nx,ny,nz)       )   
      allocate( v(nx,ny,nz)       )   
      allocate( w(nx,ny,nz)       )    
      allocate( ccdepth(nx,ny,nz) )  
      allocate( acc_width(nx,ny,nz+1)   )

c     --- 2D grids ---
   
      allocate( wdepth(nx,ny)     ) 
      allocate( bottom_layer(nx,ny) )
      allocate( wetmask(nx,ny) )
      allocate( dslm(nx,ny) )
      
      write(*,*) "init_embedded_hydrography: allocate grid arrays: OK"

      end subroutine init_embedded_hydrography
       



      subroutine close_embedded_hydrography()
c     ------------------------------------------------------
c     ------------------------------------------------------
      if (allocated(u))            deallocate( u )
      if (allocated(v))            deallocate( v )
      if (allocated(w))            deallocate( w )
      if (allocated(ccdepth))      deallocate( ccdepth )     
      if (allocated(acc_width))    deallocate( acc_width )     
      if (allocated(wdepth))       deallocate( wdepth )
      if (allocated(bottom_layer)) deallocate( bottom_layer )
      if (allocated(wetmask))      deallocate( wetmask )
      if (allocated(dslm))         deallocate( dslm )

      end subroutine close_embedded_hydrography



      
      subroutine interpolate_wdepth(geo, result, status) 
c     -------------------------------------------------- 
c     Do not apply interpolate_cc_2Dgrid_data which uses
c     wet-point interpolation, but unrestricted interpolation (wdepth=0 at dry points)
c     (otherwise problems in coastal regions)
c     -------------------------------------------------- 
      real, intent(in)     :: geo(:) 
      real, intent(out)    :: result
      integer, intent(out) :: status   
c
      integer           :: ix,iy,ix0,iy0,ix1,iy1
      real              :: sx,sy,vc(4)    
c     ------------------------------------------ 
      if (.not.horizontal_range_check(geo)) then
          result = padval
          status = 1
          return
      endif
      if (is_land(geo)) then
          result = 0.0 ! fixed value for dry points
          status = 3
          return
      endif
c     --- delegate normal interior interpolation to interp_2Dbox_data
      call get_surrounding_box(geo,ix,iy,sx,sy)
      ix0 = min(max(ix,    1),nx)  ! cover boundary layers
      ix1 = min(max(ix + 1,1),nx)  ! cover boundary layers
      iy0 = min(max(iy,    1),ny)  ! cover boundary layers
      iy1 = min(max(iy + 1,1),ny)  ! cover boundary layers
      vc(1:2) = wdepth(ix0, iy0:iy1)
      vc(3:4) = wdepth(ix1, iy0:iy1)
      call interp_2Dbox_data(sx, sy, vc, 0, result)
      status = 0
c     ------------------------------------------ 
      end subroutine 


      LOGICAL function is_wet(geo)
c     ------------------------------------------ 
c     Probe wdepth
c     return is_wet = .true. at horizontal range violation
c     accept points at sea surface/bottom as wet (htol)
c     ------------------------------------------ 
      real, intent(in) :: geo(:) 
      real             :: wd
      integer          :: status
c     ------------------------------------------ 
      call interpolate_wdepth(geo, wd, status) 
      if (status==1) then     ! status = 1: horizontal range violation, set result = padval
         is_wet = .true.
         return
      elseif (status==0) then ! status = 0: interior interpolation performed
         if ((geo(3)>-htol).and.(geo(3)<wd+htol)) then
            is_wet = .true.
         else
            is_wet = .false.
         endif
      else
         stop "is_wet: unhandled return condition"
      endif
c     ------------------------------------------
      end function


      
      subroutine get_grid_coordinates(geo,x,y,z)  
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
      real, intent(in)  :: geo(:)
      real, intent(out) :: x,y,z

      integer           :: ix,iy,ix0,iy0,ix1,iy1,iz
      
      real, pointer     :: z00(:), z01(:)
      real, pointer     :: z10(:), z11(:)
      real, allocatable :: z0(:), z1(:), acclay(:)
      real              :: layerw,sx,sy,cwd
c     ------------------------------------------ 
      allocate( z0(nz+1)     )
      allocate( z1(nz+1)     )
      allocate( acclay(nz+1) )
      call get_horiz_grid_coor_scalar(geo,x,y)   ! define x,y
c
c     ----- interpolate local layer spacings from acc_width: acclay -----
c
      call get_surrounding_box(geo,ix,iy,sx,sy)
c     
      ix0 = min(max(ix,  1),nx)  ! cover grid edges
      ix1 = min(max(ix+1,1),nx)  ! cover grid edges
      iy0 = min(max(iy,  1),ny)  ! cover grid edges
      iy1 = min(max(iy+1,1),ny)  ! cover grid edges
c     --- trilin interpolation  ---
      z00 => acc_width(ix0, iy0, :)  ! range = 1:nz+1
      z01 => acc_width(ix0, iy1, :)  ! range = 1:nz+1
      z10 => acc_width(ix1, iy0, :)  ! range = 1:nz+1
      z11 => acc_width(ix1, iy1, :)  ! range = 1:nz+1
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
      iz = min(max(iz,1),nz) ! capture vertical range excess 
      layerw = acclay(iz+1) -  acclay(iz) 
      z      = (iz - 0.5) + (geo(3) - acclay(iz))/layerw ! provides smooth extrapolation at range excess 
c
      deallocate( z0, z1, acclay  )
      nullify( z00, z01, z10, z11 )
c     ------------------------------------------ 
      end subroutine get_grid_coordinates

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     inlined subroutine from horizontal_grid_transformations
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine get_surrounding_box(geo,ix,iy,sx,sy) 
c     ------------------------------------------------------
c     From the lon/lat vector geo (which may include the z component)
c     determine the cell association (ix,iy) where the node of the cell is the
c     lower left corner point, so that [ix:ix+1, iy:iy+1[ maps to cell (ix,iy)
c     i.e. a corner centered convention.
c     0 < (sx,sy) < 1 is the intra cell cooridnates. 
c     Grid origo (ix,iy) = (x,y) = (1,1) is associated with (lambda1,phi1)
c     Do not flag grid range violations, so that 
c     0 < (ix,iy) <= (nx,ny) is not enforced 
c     Inlined get_surrounding_box from horizontal_grid_transformations     
c     ------------------------------------------------------
      real, intent(in)     :: geo(:)
      integer, intent(out) :: ix,iy
      real, intent(out)    :: sx,sy
      real                 :: dx1,dy1
c     ------------------------------------------------------
      dx1 = (geo(1)-lambda1)/dlambda
      dy1 = (geo(2)-phi1)   /dphi
      ix  = 1 + int(dx1)    ! truncate decimals to get lon cell
      iy  = 1 + int(dy1)    ! truncate decimals to get lat cell
      sx  = dx1 - int(dx1)  ! intra cell coordinate 0<sx<1
      sy  = dy1 - int(dy1)  ! intra cell coordinate 0<sy<1
c     ------------------------------------------------------
      end subroutine get_surrounding_box

      
      subroutine get_horiz_grid_coor_scalar(geo,x,y)
c     -------------------------------------------------------- 
c     Compute continuous horizontal grid coordinates x,y
c     of lon-lat position geo (which may include the z components etc)
c     Mesh points are on integer values of (x,y) and the
c     first mesh point is at (x,y) = (1,1), i.e. fortran offset
c     no range check
c     Inlined get_horiz_grid_coor_scalar from horizontal_grid_transformations     
c     -------------------------------------------------------- 
      real, intent(in)   :: geo(:)
      real, intent(out)  :: x,y
c     -------------------------------------------------------- 
      x  = 1.0 + (geo(1)-lambda1)/dlambda 
      y  = 1.0 + (geo(2)-phi1)/dphi         
      end subroutine get_horiz_grid_coor_scalar    

      
      LOGICAL function horizontal_range_check(geo)
c     ------------------------------------------
c     Require grid coordinates 
c      (0.5 < x < nx+0.5) and (0.5 < y < ny+0.5)              
c     ------------------------------------------ 
      real, intent(in) :: geo(:)
      integer          :: ixc,iyc
c     ------------------------------------------ 
      call get_horiz_ncc_index(geo,ixc,iyc) 
      horizontal_range_check = (1<=ixc).and.(ixc<=nx)
     +                    .and.(1<=iyc).and.(iyc<=ny)
c     ------------------------------------------ 
      end function horizontal_range_check


      
      subroutine get_horiz_ncc_index(geo,ixc,iyc) !,xc,yc)
c     ------------------------------------------ 
c     Resolve the node-centered grid cell containing geo
c     (ixc,iyc) is grid indices of the closest node
c     (xc,yc)   is (lon, lat) coordinates of the closest node
c     Notice that (ixc,iyc) is not coinciding with (ix,iy)
c     in get_horiz_grid_coordinates(geo,ix,iy,sx,sy) above
c     no range check
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      integer, intent(out) :: ixc,iyc
c     ------------------------------------------ 
      ixc = 1 + nint((geo(1)-lambda1)/dlambda)     
      iyc = 1 + nint((geo(2)-phi1)   /dphi)   
      end subroutine get_horiz_ncc_index

      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     inlined from horizontal_representation
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      LOGICAL function is_land(geo)        ! COPY -> horizontal_representation_stay_wet.f
c     ------------------------------------------
c     Inlined function from horizontal_representation
c     return is_land = .false. at horizontal range violation
c     ------------------------------------------ 
      real, intent(in) :: geo(:)
      integer          :: ixc,iyc
c     ------------------------------------------ 
      if (horizontal_range_check(geo)) then ! range OK
         call get_horiz_ncc_index(geo,ixc,iyc)
         is_land = (wetmask(ixc,iyc)==0)
      else
         is_land = .false.
      endif    
c     ------------------------------------------
      end function

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     inlined from opec.f
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      
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
c
c     The boundary condition:
c           w( ix, iy, bottom_layer(ix,iy)+0.5 ) = 0
c     is implicit 
c
c     Implied interpolation ranges in ncc coordinates:
c           1.5 < x < nx+0.5
c           0.5 < y < ny-0.5
c           0.5 < z < nz+0.5 (due to bottom BC)
c
c     this means that the proper interpolation domain is different from the simulation domain
c     definition derived from the grid. Project these rim points onto the
c     interior domain and interpolate currents here and do not flag them as domain violation
c
c     Inlined subroutine interpolate_currents(xyz, uvw, status) from opec.f
c     ------------------------------------------ 
      real, intent(in)     :: xyz(:)  ! (lon,lat,depth)
      real, intent(out)    :: uvw(:)
      integer, intent(out) :: status

      integer              :: statu, statv, statw
      integer              :: ix,iy,iz,jwest,jnorth,jlow,ibot
      real                 :: x,y,z, sx,sy,sz
c     ------------------------------------------ 
      if (.not.horizontal_range_check(xyz)) then
         uvw    = 0.
         status = 1  ! signal range violation
         return
      endif
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
      if ((x<1.5).or.(x>(nx+0.5))) statu = 1 ! flag x range violation
      if ((y<0.5).or.(y>(ny-0.5))) statv = 1 ! flag y range violation
      if ((z<0.5).or.(z>(nz+0.5))) statw = 1 ! flag z range violation
c
c     Currently ignore the mismatch between interpolation and simulation domain 
c     original test: status  = max(statu, statv, statw) 
c
      status = 0 ! signal all OK, after range/wet check

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
      jnorth = max(1, min(nint(y+1), ny)) ! cell north fce
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

      
      end module
