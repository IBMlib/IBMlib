ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     -----------------------------------------------------------
c     NORWEgian ECOlogical Model (NORWECOM) Oceanography Provider
c     Grid components
c     -----------------------------------------------------------
c     $Rev: 218 $
c     $LastChangedDate: 2011-01-25 15:09:13 +0100 (Tue, 25 Jan 2011) $
c     $LastChangedBy: mpay $ 
c
c     This component of the NORWECOM oceanography provider deals with
c     the grid and grid-related functions. It is in included by 
c     the NORWECOM_main.f file. Brief details of the grid are found
c     below - more detailed information can be found in the docs 
c     dirctory
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
c       Unsure - how are the velocities stored? Check
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      character*100 function grid_component_version() 
c     Writes the SVN revision number out for logging purposes
      grid_component_version =  "Grid comp     : $Rev: 218 $"
      end function

c########################################################################
c########################################################################
c#           Grid related transformations and functions                 # 
c########################################################################
c     Import NORWECOM grid <-> lat,lon routines
      include "gr2sph.f"
      include "sph2gr.f"

      subroutine cell_index(xyz,ix,iy,iz)
c     -------------------------------------------------------
c     Returns the cell index for a given set of grid coordinates
c     Note that this function is completely unchecked for the 
c     validity of the range - this is the responsbility of
c     data_grid_check()
c     -------------------------------------------------------
      real, intent(in)  :: xyz(:)  
      integer, intent(out) :: ix,iy
      integer, intent(out),optional :: iz
c     --------------------------------------------------
      ix   = ceiling(xyz(1))
      iy   = ceiling(xyz(2))
      if(present(iz)) then
        if(xyz(3)==0) then
          iz=1    !Ceiling doesn't quite work properly for z=0
        else 
          iz   = ceiling(xyz(3))
        endif
      endif
      end subroutine cell_index


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

      
      subroutine get_horiz_corners(ix,iy,c00,c01,c10,c11) 
c     ----------------------------------------------------------------
c     Returns the indices of the corners corners of the cell ix, iy 
c     in grid space
c     ----------------------------------------------------------------
      integer, intent(in) :: ix,iy
      real, intent(out)   :: c00(2),c01(2),c10(2),c11(2)
c     ----------------------------------------------------------------
      c00(1) = ix - 1.0
      c00(2) = iy - 1.0
      c01(1) = ix - 1.0
      c01(2) = iy
      c10(1) = ix 
      c10(2) = iy - 1.0
      c11(1) = ix       
      c11(2) = iy 
      end subroutine get_horiz_corners 
      

      logical function grid_range_check(xyz)
c     ------------------------------------------------------
c     Check xyz wrt. interpolation ranges in grid space:
c           0.5 <= x < nx-1
c           0.5 <= y < ny-1
c           0   <= z <= nz (due to bottom BC)
c     which is the relevant range for interpolation of all parameters
c     Within this range position xyz is always bracketed in each
c     direction, allowing interpolation
c     return .true. if xyz is interior to above ranges, else .false.
c     For safety, we define on the line as out
c     ------------------------------------------------------
      real,    intent(in)  :: xyz(:)
      integer :: cross(3)
c     ------------------------------------------------------
      cross  = 0
      if  (xyz(1) < 0.5)    cross(1) = -1
      if  (xyz(1) >= nx-1.0) cross(1) =  1
      if  (xyz(2) < 0.5)    cross(2) = -1
      if  (xyz(2) >= ny-1.0) cross(2) =  1
      if(size(xyz)>2) then
        if  (xyz(3) < 0.0)   cross(3) = -1
        if  (xyz(3) > nz)     cross(3) =  1
      endif
      grid_range_check = .not.any(cross.ne.0)
      end function grid_range_check


      logical function geo_range_check(geo)
c     ------------------------------------------------------
c     Checks whether a geo position (lon, lat, depth)
c     is within the range of the grid or not. Checks for
c     positions above the surface and below the bottom
c     in addition to the more standard horizontal checks
c     The actual check is farmed out to grid_range_check
c     ------------------------------------------------------
      real,    intent(in)  :: geo(3)  !3+
      real    :: xyz(3)
c     ------------------------------------------------------
      call lonlat2xy(geo(1:2),xyz(1:2))
      geo_range_check = grid_range_check(xyz(1:2))
      if(geo_range_check) then  !Only check the vertical if the horizontal is ok
         call depth2z(geo(3),xyz(1:2),xyz(3))
         geo_range_check = grid_range_check(xyz)
      endif
      end function geo_range_check


      logical function horizontal_range_check(ll)
c     ------------------------------------------------------
c     Public interface that checks whether a geo position (lon, lat)
c     is within the range of the grid or not
c     The actual check is farmed out to grid_range_check
c     ------------------------------------------------------
      real,    intent(in)  :: ll(:)
      real    :: xy(size(ll))
c     ------------------------------------------------------
      call lonlat2xy(ll,xy)
      horizontal_range_check = grid_range_check(xy(1:2))
      end function horizontal_range_check


      logical function is_wet(ll)
c     ------------------------------------------------------
c     Check whether geo (lon,lat,depth) corresponds to a wet position or not by
c     looking up in the wet array
c     Assert validity the ll position. Check should be handled elsewhere 
c     ------------------------------------------------------
      real,    intent(in) :: ll(:)
      real     :: xy(size(ll))
      integer  :: ix, iy
c     ------------------------------------------------------
      call lonlat2xy(ll,xy)
      call cell_index(xy(1:2),ix,iy)    
      is_wet = .not. land(ix,iy)
      end function is_wet


      logical function is_land(geo)
c     ------------------------------------------------------
c     Check whether geo (lon,lat,depth) corresponds to a wet position or not by
c     looking up in the wet array
c     ------------------------------------------------------
      real,    intent(in) :: geo(:)
c     ------------------------------------------------------
      is_land = .not. is_wet(geo)
      end function is_land


      subroutine lonlat2xy(ll, xy)
c     -------------------------------------------------------
c     Horizontal vectorized coordinate transformation (longitude/latitude -> xy) 
c     Acts as a wrapper to the sph2gr function direct 
c     from NORWECOM 
c     Validity of the xy and ll range is not checked
c     -------------------------------------------------------
      real, intent(in)  :: ll(:)  ! assumed shape (2 <= 
      real, intent(out) :: xy(:)  ! assumed shape (do not check consistency with lp)
      integer           :: fail
c     --------------------------------------------------
      call sph2gr(ll(1),ll(2), XPOLE, YPOLE, DX, DY, ALPHA,
     >                  xy(1), xy(2), fail)
      if(fail .NE. 0) then
         write(*,*) "Illegal coordinate conversion, sph2gr()"
         write(*,932) "lon   ", ll(1)
         write(*,932) "lat   ", ll(2)
         write(*,932) "x     ", xy(1)
         write(*,932) "y     ", xy(2)
         write(*,*) "failure code = ", fail
932      format(a," = ",f12.7)
         stop
      endif
      end subroutine lonlat2xy


      subroutine xy2lonlat(xy,ll)
c     -------------------------------------------------------
c     Horizontal single coordinate transformation (xy -> longitude/latitude)
c     Wraps the gr2sph function direct from NORWECOM
c     Validity of the xy and ll range is not checked
c     -------------------------------------------------------
      real, intent(in)  :: xy(:) 
      real, intent(out) :: ll(:)
      integer           :: fail
c     --------------------------------------------------
      call gr2sph(xy(1), xy(2), XPOLE, YPOLE, DX, DY, ALPHA,
     >                  ll(1), ll(2), fail)
      if(fail .NE. 0) then
         write(*,*) "Illegal coordinate conversion, gr2sph()"
         write(*,931) "x     ", xy(1)
         write(*,931) "y     ", xy(2)
         write(*,931) "lon   ", ll(1)
         write(*,931) "lat   ", ll(2)
         write(*,*) "failure code = ", fail
931      format(a," = ",f12.7)
         stop
      endif
      end subroutine xy2lonlat


      subroutine depth2z(depth,xy,z)
c-------------------------------------------------------
c     geographical single coordinate transformation at
c     a given xy position (usually obtained first from
c     lonlat2xy) ie
c     depth | xy -> z 
c     Assert that xy is a valid range
c     Interpolation /extrapolation functionality is maintained
c     for historical consistency. However, this may be removed
c     in the future
c-------------------------------------------------------
      real, intent(in)  :: depth,xy(:)
      real, intent(out) :: z
      integer           :: ix, iy, iz
      real              :: sz
c--------------------------------------------------
c     Get parent cells
      call cell_index(xy(1:2),ix,iy)

c     Are we on land? If so, then that's it - set the z dim
c     to a very small but positive number  and exit. 
      if(land(ix,iy)) then
        z = 1e-16
        return
      endif

c     Else transform the vertical coordinates
      if (depth <= 0)  then                                  ! above surface (z=0.0)
         z =  depth/layer_width(ix,iy,1)               ! z < 0.0
      elseif (depth >= bath(ix,iy)) then           ! below grid bottom
         sz = (depth - bath(ix,iy))/layer_width(ix,iy,nz)
         z = sz + nz
      else                                      ! within normal range: 0.0 < z < nz
         do iz=1,nz ! cell association -> iz
            if (acc_width(ix,iy,iz+1)>depth) exit
         enddo
         sz = (depth - acc_width(ix,iy,iz))/layer_width(ix,iy,iz) ! relative cell completion
         z = real(iz) -1 + sz
      endif
      end subroutine depth2z  


      subroutine z2depth(z,xy,depth)
c-------------------------------------------------------
c     geographical single coordinate transformation at
c     a given xy position (usually obtained first from
c     lonlat2xy) ie
c     z | xy -> depth
c     Assert that xy is a valid range
c     Interpolation /extrapolation functionality is maintained
c     for historical consistency. However, this may be removed
c     in the future
c-------------------------------------------------------
      real, intent(in)  :: z,xy
      real, intent(out) :: depth
      real              ::  sz 
      integer           :: ix, iy, iz
      logical           :: rangeOK
c--------------------------------------------------
c     Now interpolate/extrapolate as necessary
      if (z<=0.0) then      ! above surface (z<0.0), extrapolate from first layer
         depth = z*layer_width(ix,iy,1)  ! depth < 0 as z<0
      elseif (z>=nz) then ! below grid bottom
         sz    =  z - real(nz)
         depth = bath(ix,iy) + sz*layer_width(ix,iy,nz)
      else     ! z is inside valid range
         sz    = z-floor(z) 
         depth = acc_width(ix,iy,iz) + sz*layer_width(ix,iy,iz)
      endif
      end subroutine z2depth

      
      subroutine get_jacobian(xyz, jacob)
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
      real, intent(in)   :: xyz(:) ! shape (2+) 
      real, intent(out)  :: jacob(:,:)  ! shape (2+,2+)
      real               :: x,y      ! Positions in xyz grid space
      integer            :: ix,iy    ! Rounded positions in xyz grid space
      real               :: d, g     ! Temporary variables used in calculations
      real               :: lon,lat     ! lon-lat at position
c     --------------------------------------------------
      x = xyz(1)
      y = xyz(2)

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


      subroutine get_inv_jacobian(xyz, invjacob)
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
      real, intent(in)   :: xyz(:) ! shape (2+) 
      real, intent(out)  :: invjacob(2,2)  ! shape (2+,2+)
      real               :: jacob(2,2),a,b,c,d
c     --------------------------------------------------
      call get_jacobian(xyz,jacob)
      a=jacob(1,1)
      b=jacob(1,2)
      c=jacob(2,1)
      d=jacob(2,2)
      invjacob(1,:) =(/d,-b/)
      invjacob(2,:) =(/-c,a/)
      invjacob = invjacob /(a*d-b*c)
      end subroutine get_inv_jacobian


      subroutine resolve_vector(xyz,N,G)
c------------------------------------------------------- 
c     Handles the transformation from NORWECOM grid oriented
c     vectors to geo grid oriented vectors
c     Input: vector N located at NORWECOM grid position xyz,
c       and with components N(1), N(2) and N(3) oriented along
c       the NORWECOM axes x,y,z
c     Output: vector G of the same magnitude as N, but represented
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
      real, intent(in)   :: xyz(:), N(:) ! shape (3+) 
      real, intent(out)  :: G(:)  ! shape (3+)
      real               :: J(2,2), B(2,2), A(2,2)
c     --------------------------------------------------
      call get_jacobian(xyz,J)
c     Calculate the normalisation matrix, B
      B=0.0
      B(1,1) = 1/sqrt(J(1,1)**2 + J(2,1)**2)
      B(2,2) = 1/sqrt(J(1,2)**2 + J(2,2)**2)
c     Calculate the transformation matrix, A
      A = matmul(J,B)
c     Now do the transformation
      G(1:2) = matmul(A,N(1:2))
      G(3) = N(3)   !No transform required 
      end subroutine resolve_vector



