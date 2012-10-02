ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     -----------------------------------------------------------
c     NORWEgian ECOlogical Model (NORWECOM) Oceanography Provider
c     Grid modules
c     -----------------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
c     This component of the NORWECOM oceanography provider deals with
c     the grid and grid-related functions. It is in included as a module by 
c     the NORWECOM_main.f and by the horizontal_representation module. 
c     Brief details of the grid are found below - more detailed 
c     information can be found in the docs dirctory
c
c     Grid type: Rotated polar-stereographic grid, with sigma
c                coordinates in the vertical. Arakawa C staggering
c     Define grid according to NORWECOM
c       u is midt point at cell west   face and positive eastward   unit == m/s
c       v is midt point at cell south  face and positive northward  unit == m/s
c       w is midt point at cell upper  face and positive upward     unit == m/s
c     Coordinate definitions:
c       See sph2gr.f for description of coordinate system. Generally speaking, the
c       (i,j)th cell has coordinates (i-0.5,j-0.5). The SW corner of the
c       first cell therefore has the coordinates (0,0).
c     Orientations:
c       Velocities u and v are oriented parallel to the NORWECOM grid axes
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      module horizontal_grid_transformations
      use constants
      use geometry
      implicit none
      private

c     ----------Standard Interface ----------------
      public :: init_horiz_grid_transf   ! canonical name too long ...
      public :: close_horiz_grid_transf  ! canonical name too long ...
      public :: horizontal_range_check
      public :: get_horiz_grid_coordinates
      public :: get_horiz_geo_coordinates
      public :: get_surrounding_box
      public :: get_horiz_ncc_index
      public :: get_horiz_ncc_corners
      interface get_horiz_grid_coordinates
         module procedure get_horiz_grid_coor_scalar
         module procedure get_horiz_grid_coor_vector
      end interface
      interface get_horiz_geo_coordinates
         module procedure get_horiz_geo_coor_scalar
         module procedure get_horiz_geo_coor_vector
      end interface

c     ----------Extras interface ----------------
      public ::  get_jacobian
      public ::  get_inv_jacobian
      public ::  resolve_vector
      public ::  grid_range_check
      public ::  get_cell_index_from_xy
      public ::  get_CC_cube

c     full grid dimensions, including grid specifier
      integer, public,parameter::  nx=160
      integer, public,parameter::  ny=140
      real, public,parameter   ::  DX =10
      real, public,parameter   ::  DY =10
      real, parameter          ::  XPOLE=382
      real, parameter          ::  YPOLE=256
      real, parameter          ::  ALPHA=58
      real, parameter          ::  lat_true=60
      real  PHINUL
      
      real, parameter          :: horz_buf=0.5    !Tolerance, in grid units, around practical boundary

c.....lat/lon grid auxillaries      
      real, allocatable,public   :: xfacelen(:,:)  ! length of the southern face of cell  [m]
      real, allocatable,public   :: yfacelen(:,:)  ! length of the western  face of cell  [m]
      real, allocatable,public   :: cell_area(:,:) ! horizontal area of cell (i,j) in meters**2

c     ----------------
      contains
c     ----------------

c     Import NORWECOM grid <-> lat,lon routines
      include "gr2sph.f"
      include "sph2gr.f"


      subroutine init_horiz_grid_transf() 
c     ------------------------------------------
c     Sets up key aspects of the horizontal gri
      real                 :: d(3),pos1(3),pos2(3),geopos1(3)
      real                 :: geopos2(3),tmp
      integer              :: ix,iy
c     ------------------------------------------
c     Setup PHINUL, the radian equivalent of lat_trueL
      PHINUL = lat_true * deg2rad        

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
            call get_horiz_geo_coordinates(pos1,geopos1)            
            call get_horiz_geo_coordinates(pos2,geopos2)            
            call get_horizontal_distance(geopos1,
     +               geopos2,xfacelen(ix,iy))
            
c           Then the yfacelengths 
            pos1(1)=ix-1.0   !yfacelen is defined at the west side of cell
            pos1(2)=iy-1.0   !y-dir bottom of cell
            pos2=pos1
            pos2(2)=iy   !y-dir top of cell
            call get_horiz_geo_coordinates(pos1,geopos1)            
            call get_horiz_geo_coordinates(pos2,geopos2)            
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
      end subroutine init_horiz_grid_transf
      


      subroutine close_horiz_grid_transf() 
c     ------------------------------------------
c     Empty for this grid type
c     ------------------------------------------
      end subroutine close_horiz_grid_transf


      subroutine get_horiz_grid_coor_vector(geo,xy)  
c     -------------------------------------------------------
c     Horizontal vectorized coordinate transformation (longitude/latitude -> xy) 
c     Acts as a wrapper to the sph2gr function direct 
c     from NORWECOM 
c     Validity of the xy and geo range is not checked
c     -------------------------------------------------------
      real, intent(in)  :: geo(:)  ! assumed shape (2 <= 
      real, intent(out) :: xy(:)  ! assumed shape (do not check consistency with lp)
      integer           :: fail
c     --------------------------------------------------
      call sph2gr(geo(1),geo(2), XPOLE, YPOLE, DX, DY, ALPHA,
     >                  xy(1), xy(2), fail)
      if(fail .NE. 0) then
         write(*,*) "Illegal coordinate conversion, sph2gr()"
         write(*,*) "geo: ", geo
         write(*,*) "xy : ", xy
         write(*,*) "failure code = ", fail
         stop
      endif
      end subroutine get_horiz_grid_coor_vector


      subroutine get_horiz_grid_coor_scalar(geo,x,y)
c     -------------------------------------------------------- 
c     Compute continuous horizontal grid coordinates (x,y)
c     of lon-lat position xy (which may include the z component)
c     Mesh points are on integer values of (x,y) and the
c     first mesh point is at (x,y) = (1,1), i.e. fortran offset
c     no range check     
c     -------------------------------------------------------- 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: x,y
      real    :: xy(size(geo))
c     -------------------------------------------------------- 
      call get_horiz_grid_coor_vector(geo,xy)
      x=xy(1)
      y=xy(2)
      end subroutine get_horiz_grid_coor_scalar


      subroutine get_horiz_geo_coor_vector(xy,geo)
c     --------------------------------------------------------  
c     Horizontal single coordinate transformation (xy -> longitude/latitude)
c     Wraps the gr2sph function direct from NORWECOM
c     Validity of the xy and geo range is not checked
c     --------------------------------------------------------  
      real, intent(in)    :: xy(:)
      real, intent(out)   :: geo(:)
      integer           :: fail
c     --------------------------------------------------
      call gr2sph(xy(1), xy(2), XPOLE, YPOLE, DX, DY, ALPHA,
     >                  geo(1), geo(2), fail)
      if(fail .NE. 0) then
         write(*,*) "Illegal coordinate conversion, gr2sph()"
         write(*,*) "xy  : ", xy
         write(*,*) "geo : ", geo
         write(*,*) "failure code = ", fail
         stop
      endif
      end subroutine get_horiz_geo_coor_vector


      subroutine get_horiz_geo_coor_scalar(x,y,geo)
c     -------------------------------------------------------- 
c     Compute lon-lat position geo corresponding to 
c     continuous horizontal grid coordinates (x,y)
c     Buffer geo may contain z value 
c     -------------------------------------------------------- 
      real, intent(in)    :: x,y
      real, intent(out)   :: geo(:)
c     -------------------------------------------------------- 
      call get_horiz_geo_coor_vector((/x,y/),geo)
      end subroutine get_horiz_geo_coor_scalar


      logical function horizontal_range_check(geo)
c     ------------------------------------------------------
c     Public interface that checks whether a geo position (lon, lat, possibly z)
c     is within the range of the grid or not
c     Check is farmed out to grid_range_check
c     ------------------------------------------------------
      real,    intent(in)  :: geo(:)
      real    :: xyz(size(geo))
c     ------------------------------------------------------
      !Don't bother checking if its really obviously wrong
      if(geo(2)<-90 .or. geo(2)>90) then
         horizontal_range_check = .FALSE.
         return         
      endif
      call get_horiz_grid_coordinates(geo,xyz)
      horizontal_range_check = grid_range_check(xyz)
      end function horizontal_range_check


      logical function grid_range_check(xy)
c     ------------------------------------------------------
c     Check xyz wrt. interpolation ranges in grid space:
c           0.5 < x < nx-1
c           0.5 < y < ny-1
c     which is the relevant range for interpolation of all parameters
c     Within this range position xy is always bracketed in each
c     direction, allowing interpolation
c     return .true. if xy is interior to above ranges, else .false.
c     For safety, we define on the line as out
c     ------------------------------------------------------
      real,    intent(in)  :: xy(:)
      integer :: cross(2)
c     ------------------------------------------------------
      cross  = 0
      if  (xy(1) <= 0.5 + horz_buf)    cross(1) = -1
      if  (xy(1) >= nx-1-horz_buf) cross(1) =  1
      if  (xy(2) <= 0.5+horz_buf)    cross(2) = -1
      if  (xy(2) >= ny-1-horz_buf) cross(2) =  1
      grid_range_check = .not.any(cross.ne.0)
      end function grid_range_check


      subroutine get_surrounding_box(geo,ix,iy,sx,sy)
c     ------------------------------------------------------
c     For any given point in grid space (x,y), there will be four nodes that immediately
c     surround it. This function determines the cell index, (ix,iy) corresponding to 
c     the lower-left of these nodes (the others can then simply be calculated by 
c     adding 1 in each direction). We then return the relative position 
c     in the x and y directions  between the corresponding notes (sx,sy) 
c     ie the relative coordinates of (x,y) in the  box formed by these four nodes
c     0 < (sx,sy) < 1 
c     Do not flag grid range violations
c     NB: previously, this was subroutine called get_horiz_grid_coordinates
c         and get_location_in_mesh
c     ------------------------------------------------------
      real, intent(in)     :: geo(:)
      integer, intent(out) :: ix,iy
      real, intent(out)    :: sx,sy
      real                 :: dx1,dy,x,y
c     ------------------------------------------------------
      call get_horiz_grid_coordinates(geo,x,y)
      ix = nint(x)
      iy = nint(y)
      sx  = x - real(ix) +0.5   ! intra cell coordinate 0<sx<1
      sy  = y - real(iy) +0.5   ! intra cell coordinate 0<sy<1
c     ------------------------------------------------------
      end subroutine get_surrounding_box


      subroutine get_CC_cube(xyz,ixyz0,ixyz1,sxyz)
c     ----------------------------------------------------------------
c     Returns the indices of the corners corners of a 1x1x1 cube in 
c     grid space that encloses the point xyz. The corner indices can
c     thus be used in a lookup function to get data for 
c     interpolation. Also return sxyz, the relative position
c     of the point xyz within that cube (0 <= s <= 1)
c     Both land and surface/bottom breaking are handled elsewhere
c     ----------------------------------------------------------------
      real, intent(in)    :: xyz(:)
      integer, intent(out) :: ixyz0(:), ixyz1(:)
      real, intent(out) :: sxyz(:)
c     ----------------------------------------------------------------
c     determine corners around which to interpolate
      ixyz0  = nint(xyz)   !SW corner ie (0,0,0)
      ixyz1   = ixyz0 + 1  !NE corner ie (1,1,1)
c     determine locatin relative to corners (0 <= s <= 1)
      sxyz    = xyz-real(ixyz0)+0.5
      end subroutine get_CC_cube


      subroutine get_horiz_ncc_index(geo,ixc,iyc) 
c     ------------------------------------------ 
c     Resolve the indices of the node-centered grid cell containing geo
c     Actual calculation is farmed out
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      integer, intent(out) :: ixc,iyc
      real  :: xy(2)
      integer :: ixy(2)
c     ------------------------------------------ 
      call get_horiz_grid_coordinates(geo,xy)
      call get_cell_index_from_xy(xy,ixy)
      ixc=ixy(1)
      iyc=ixy(2)
      end subroutine get_horiz_ncc_index


      subroutine get_cell_index_from_xy(xy,ixy)
c     ------------------------------------------ 
c     Resolve the indices of the node-centered grid cell 
c     containing the grid coordinates, xy
c     Notice that (ixc,iyc) is not the same as (ix,iy)
c     in get_surrounding_box(geo,ix,iy,sx,sy) above
c     no range check
c     ------------------------------------------ 
      real, intent(in)     :: xy(:)
      integer, intent(out) :: ixy(:)
c     ------------------------------------------ 
      ixy(1)=ceiling(xy(1))
      ixy(2)=ceiling(xy(2))
      end subroutine get_cell_index_from_xy


      subroutine get_horiz_ncc_corners(ixc,iyc,c00,c01,c10,c11)
c     ------------------------------------------ 
c     Resolve corners (c00,c01,c10,c11) in grid units of the 
c     node-centered horizontal grid cell with indices (ixc,iyc)
c     Previously, this function returned lon-lat of corners
c     ------------------------------------------ 
      integer, intent(in) :: ixc,iyc
      real, intent(out)   :: c00(2),c01(2),c10(2),c11(2)
c     ------------------------------------------     
      c00(1) = ixc - 1.0
      c00(2) = iyc - 1.0
      c01(1) = ixc - 1.0
      c01(2) = iyc
      c10(1) = ixc 
      c10(2) = iyc - 1.0
      c11(1) = ixc       
      c11(2) = iyc 
c     ------------------------------------------ 
      end subroutine get_horiz_ncc_corners


      subroutine get_jacobian(xy, jacob)
c------------------------------------------------------- 
c     2D jacobian for the transformation from grid space to
c     lon-lat-vertical representation 
c     A displacement in grid coordinates dxy from position xyz
c     corresponds to the displacement 
c       dll = jacob %*% dxy
c     in lon-lat coordinate representation
c     Expressions for the individual elements of the jacobian
c     are given as Equation 17 in the NORWECOM oceanography
c     provider documentation
c     Validity of the xyz range is not checked
c------------------------------------------------------- 
      real, intent(in)   :: xy(:) ! shape (2+) 
      real, intent(out)  :: jacob(:,:)  ! shape (2+,2+)
      real               :: x,y      ! Positions in xyz grid space
      integer            :: ix,iy    ! Rounded positions in xyz grid space
      real               :: d, g     ! Temporary variables used in calculations
      real               :: lon,lat     ! lon-lat at position
c     --------------------------------------------------
c     Shorthand
      x = xy(1)
      y = xy(2)

c     First the constants 
      d =sqrt(DX**2*(x-XPOLE)**2 + DX**2 * (y-YPOLE)**2)
      g = earth_radius/1000*(1+sin(PHINUL))  !earth_radius is in m, whereas DX, DY are in km
c     Changes in latitude      
      jacob(1,1:2)= DX*DY/(d**2)
      jacob(1,1) = -jacob(1,1)*(y-YPOLE)
      jacob(1,2) = jacob(1,2)*(x-XPOLE)
c     Then the changes in latitude      
      jacob(2,1:2) = -2 /(g*d*(1+(d/g)**2))
      jacob(2,1) = jacob(2,1)*DX**2*(x-XPOLE)
      jacob(2,2) = jacob(2,2)*DY**2*(y-YPOLE)
      end subroutine get_jacobian


      subroutine get_inv_jacobian(xy, invjacob)
c------------------------------------------------------- 
c     2D transformation from lon-lat-vertical representation to xyz grid-space
c     A displacement in spherical coordinates dll from grid position xyz
c     corresponds to the displacement 
c       dxy = invjacob %*% dll
c     in xyz grid coordinate representation
c     Code uses get_jacobian function as its basis, evaluating the
c     jacobian at position xyz, and then inverting the resultant matrix
c     In this way, changes in get_jacobian propigate directly here
c     Validity of the xyz range is not checked
c------------------------------------------------------- 
      real, intent(in)   :: xy(:) ! shape (2+) 
      real, intent(out)  :: invjacob(2,2)  ! shape (2+,2+)
      real               :: jacob(2,2),a,b,c,d
c     --------------------------------------------------
      call get_jacobian(xy,jacob)
      a=jacob(1,1)
      b=jacob(1,2)
      c=jacob(2,1)
      d=jacob(2,2)
      invjacob(1,:) =(/d,-b/)
      invjacob(2,:) =(/-c,a/)
      invjacob = invjacob /(a*d-b*c)
      end subroutine get_inv_jacobian


      subroutine resolve_vector(xy,N,G)
c------------------------------------------------------- 
c     Handles the transformation from NORWECOM grid oriented
c     displacement vectors to geo grid oriented displacement vectors
c     Input: displacement vector N located at NORWECOM grid position xyz,
c       and with components N(1), N(2) and N(3) oriented along
c       the NORWECOM axes x,y,z
c     Output: displacement vector G of the same magnitude as N, but represented
c       instead in terms of the components G(1), G(2) and G(3) that
c       are oriented parallel to the lon, lat, and depth axes
c     The transformation between the two is given by 
c       G = A %*% N
c       A = J %*% B
c     where J is the jacobiani, A is the transformation matrix and 
c     B is a matrix that normalises the Jacobian. For more details, 
c     see the NORWECOM oceanography provider documentation.
c     Note that as the vertical axes in both NORWECOM and lon-lat-depth
c     are coincident, this axis is essentially unchanged and thus we only
c     need to deal with the horizontal axes.
c     Validity of the xyz range is not checked
c------------------------------------------------------- 
      real, intent(in)   :: xy(:), N(:) ! shape (3+) 
      real, intent(out)  :: G(:)  ! shape (3+)
      real               :: J(2,2), C(2,2),B(2,2), A(2,2),ll(2),phi
c     --------------------------------------------------
      call get_jacobian(xy,J)
c     Convert NORWECOM position  to lon lat
      call get_horiz_geo_coordinates(xy,ll)
      phi=ll(2)*deg2rad 
c     Calcualte the sphereical scaling matrix, B
      B=0.0
      B(1,1)=cos(phi)
      B(2,2)=1.0
c     Calculate the normalisation matrix, C
      C=0.0
      C(1,1) = 1/sqrt(J(1,1)**2*cos(phi)**2 + J(2,1)**2)
      C(2,2) = 1/sqrt(J(1,2)**2*cos(phi)**2 + J(2,2)**2)
c     Calculate the transformation matrix, A
      A = matmul(J,C)
      A = matmul(B,A)
c     Now do the transformation
      G(1:2) = matmul(A,N(1:2))
      G(3) = N(3)   !No transform required 
      end subroutine resolve_vector




      end module

