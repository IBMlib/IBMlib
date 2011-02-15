      module horizontal_grid_transformations
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Provider of horizontal grid transformation interface
c     for regular lon-lat grid 
c     
c     $Rev: $
c     $LastChangedDate:  $
c     $LastChangedBy: $ 
c      
c     Coordinate conventions for grid space: 
c       * mesh points at integer value grid coordinates, 
c         first mesh point at (1,1) in grid space
c       * corner centered cells are those given by the mesh
c       * node centered cells has centers at mesh points
c         and cell faces at half-valued lines in grid space
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      private
      integer, public :: nx   ! mesh generic, set directly by user
      integer, public :: ny   ! mesh generic, set directly by user

      real :: lambda1         ! set via init_horiz_grid_transf
      real :: dlambda         ! set via init_horiz_grid_transf
      real :: phi1            ! set via init_horiz_grid_transf
      real :: dphi            ! set via init_horiz_grid_transf

      public :: init_horiz_grid_transf   ! canonical name too long ...
      public :: close_horiz_grid_transf  ! canonical name too long ...
      public :: horizontal_range_check
      public :: get_horiz_grid_coordinates
      public :: get_horiz_geo_coordinates
      public :: get_location_in_mesh
      public :: get_horiz_ncc
      public :: get_horiz_ncc_corners

      interface get_horiz_grid_coordinates
         module procedure get_horiz_grid_coor_scalar
         module procedure get_horiz_grid_coor_vector
      end interface
      
      interface get_horiz_geo_coordinates
         module procedure get_horiz_geo_coor_scalar
         module procedure get_horiz_geo_coor_vector
      end interface

c     ----------------
      contains
c     ----------------


      subroutine init_horiz_grid_transf(l1,p1,dl,dp) 
c     ------------------------------------------
      real, intent(in) :: l1,p1,dl,dp
c     ------------------------------------------
      write(*,*) "initialize interface horizontal_grid_transformations"
      write(*,*) "grid type = regular lon-lat"
      write(*,*) "implementation = "//
     +           "horizontal_grid_transformations_lonlatgrid"
      write(*,232) lambda1, dlambda
      write(*,233) phi1,    dphi    
 232  format("init_horiz_grid_transf: lambda1 = ",f12.7,
     +                      " dlambda = ",f12.7," degrees")  
 233  format("init_horiz_grid_transf: phi1    = ",f12.7,
     +                      " dphi    = ",f12.7," degrees") 
 
      lambda1 = l1 ! set module data
      phi1    = p1 ! set module data
      dlambda = dl ! set module data
      dphi    = dp ! set module data
      
      end subroutine init_horiz_grid_transf
      


      subroutine close_horiz_grid_transf() 
c     ------------------------------------------
c     Empty for this grid type
c     ------------------------------------------
      end subroutine close_horiz_grid_transf



      LOGICAL function horizontal_range_check(geo)
c     ------------------------------------------
c     Require grid coordinates 
c      (0.5 < x < nx+0.5) and (0.5 < y < ny+0.5)              
c     ------------------------------------------ 
      real, intent(in) :: geo(:)
      integer          :: ixc,iyc
c     ------------------------------------------ 
      call get_horiz_ncc(geo,ixc,iyc) 
      horizontal_range_check = (1<=ixc).and.(ixc<=nx)
     +                    .and.(1<=iyc).and.(iyc<=ny)
c     ------------------------------------------ 
      end function horizontal_range_check



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
c     -------------------------------------------------------- 
      x  = 1.0 + (geo(1)-lambda1)/dlambda 
      y  = 1.0 + (geo(2)-phi1)/dphi         
      end subroutine get_horiz_grid_coor_scalar

      
      subroutine get_horiz_grid_coor_vector(geo,xy)  
c     -------------------------------------------------------- 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: xy(:)
c     -------------------------------------------------------- 
      call get_horiz_grid_coor_scalar(geo,xy(1),xy(2))
      end subroutine get_horiz_grid_coor_vector



      subroutine get_horiz_geo_coor_scalar(x,y,geo)
c     -------------------------------------------------------- 
c     Compute lon-lat position geo corresponding to 
c     continuous horizontal grid coordinates (x,y)
c     Buffer geo may contain z value 
c     -------------------------------------------------------- 
      real, intent(in)    :: x,y
      real, intent(out)   :: geo(:)
c     -------------------------------------------------------- 
      geo(1) = lambda1 + (x-1.0)*dlambda 
      geo(2) = phi1    + (y-1.0)*dphi 
      end subroutine get_horiz_geo_coor_scalar

      
      subroutine get_horiz_geo_coor_vector(xy,geo)
c     --------------------------------------------------------  
      real, intent(in)    :: xy(:)
      real, intent(out)   :: geo(:)
c     -------------------------------------------------------- 
      call get_horiz_geo_coor_scalar(xy(1),xy(2),geo)
      end subroutine get_horiz_geo_coor_vector




      subroutine get_location_in_mesh(geo,ix,iy,sx,sy) ! NB previously: get_horiz_grid_coordinates
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
      end subroutine get_location_in_mesh



      subroutine get_horiz_ncc(geo,ixc,iyc) !,xc,yc)
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
      end subroutine get_horiz_ncc



      subroutine get_horiz_ncc_corners(ixc,iyc,c00,c01,c10,c11)
c     ------------------------------------------ 
c     Resolve corners (c00,c01,c10,c11) in grid units of the 
c     node-centered horizontal grid cell with indices (ixc,iyc)
c     Previously, this function returned lon-lat of corners
c     ------------------------------------------ 
      integer, intent(in) :: ixc,iyc
      real, intent(out)   :: c00(2),c01(2),c10(2),c11(2)
c     ------------------------------------------     
      c00(1) = ixc - 0.5
      c01(1) = ixc - 0.5
      c10(1) = ixc + 0.5
      c11(1) = ixc + 0.5   
      c00(2) = iyc - 0.5
      c01(2) = iyc + 0.5
      c10(2) = iyc - 0.5
      c11(2) = iyc + 0.5
c     ------------------------------------------ 
      end subroutine get_horiz_ncc_corners


      end module



