ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Geometry module
c     ---------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
c     Provides basic geography functions and transformations
c
c     TODO:  fix implementation of add_finite_step  
c   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      module geometry
      use constants

      implicit none
      private

      public ::  cross_2Dline_segments
      public ::  cross_2Dlines
      public ::  dcart2dxy
      public ::  dxy2dcart
      public ::  get_local_distance
      public ::  get_horizontal_distance
      public ::  add_finite_step   

      interface get_local_distance
         module procedure get_local_distance_scalar
         module procedure get_local_distance_vector      
      end interface


c     ============================================================
                          contains
c     ============================================================

      subroutine cross_2Dline_segments(x0,x1,y0,y1,s,cross)
c     ----------------------------------------------------------
c     Calculate whether the line from x0 to x1 crosses
c     the line from y0 to y1. s is the coordinate along the
c     vector from x0 to x1, and 0<s<1 if cross == .true.
c     If cross == .false. lines do not cross between (x0 to x1)
c     and between (y0 to y1) or lines are parallel.
c
c     Wrapper variant of cross_2Dlines to ensure only interior solution
c     ----------------------------------------------------------
      implicit none
      real, intent(in)     :: x0(*),x1(*),y0(*),y1(*)
      real, intent(out)    :: s
      logical, intent(out) :: cross
      real                 :: vx(2),vy(2),lvx2,lvy2,det,t
      real,parameter       :: s_parallel = 1.e20
      real,parameter       :: s_undef    = 2.e20
c     ----------------------------------------------------------
      call cross_2Dlines(x0,x1,y0,y1,s,t,cross)
      cross = cross .and. in_princip_range(s) .and. in_princip_range(t)

      contains
 
      logical function in_princip_range(x)
      real, intent(in)     :: x
      in_princip_range = ((x >= 0.0).and.(x <= 1.0))
      end function
c     ----------------------------------------------------------
      end subroutine cross_2Dline_segments


      subroutine cross_2Dlines(x0,x1,y0,y1,s,t,cross)
c     ----------------------------------------------------------
c     Calculate whether and where the line through (x0, x1) crosses
c     the line through (y0, y1). s is the coordinate along the
c     vector from x0 to x1, and t is the coordinate along the
c     vector from y0 to y1. The solution (s,t) may be any real values 
c     if it exists (cross == .true.). If cross == .false. lines 
c     are parallel or either (x0==x1) or (y0==y1). 
c     ----------------------------------------------------------
      implicit none
      real, intent(in)     :: x0(*),x1(*),y0(*),y1(*)
      real, intent(out)    :: s,t
      logical, intent(out) :: cross
      real                 :: vx(2),vy(2),sum_vx,sum_vy,det
      real,parameter       :: s_parallel = 1.e20
      real,parameter       :: s_undef    = 2.e20
c     ----------------------------------------------------------
      vx     = x1(1:2)-x0(1:2)
      vy     = y1(1:2)-y0(1:2)
      sum_vx = sum(abs(vx))
      sum_vy = sum(abs(vy))
      if ((sum_vx < 1.e-12).or.(sum_vx < 1.e-12)) then ! angle undef
         s = s_undef
         t = s_undef
         cross=.false.
         return
      endif   
      det = vx(2)*vy(1) - vx(1)*vy(2)     
      if (abs(det)< 1.e-12) then ! lines parallel
         s = s_parallel
         t = s_parallel
         cross=.false.
         return
      endif
c
c     --- now we know lines are not parallel or degenerate ---
c
c     Simplify[Solve[{x01 + vx1 s == y01 + vy1 t, x02 + vx2 s == y02 + vy2 t}, {s,t}]]//FortranForm
c
c              vy2 x01 - vy1 x02 - vy2 y01 + vy1 y02
c         s -> -------------------------------------, 
c                       vx2 vy1 - vx1 vy2
c     
c              vx2 x01 - vx1 x02 - vx2 y01 + vx1 y02
c         t -> -------------------------------------
c                       vx2 vy1 - vx1 vy2
c
c     condition for interior crossing: 0 < s,t < 1
c

      cross = .true.
      s = (vy(2)*x0(1) - vy(1)*x0(2) - vy(2)*y0(1) + vy(1)*y0(2))/det
      t = (vx(2)*x0(1) - vx(1)*x0(2) - vx(2)*y0(1) + vx(1)*y0(2))/det
c     ----------------------------------------------------------
      end subroutine cross_2Dlines



      subroutine get_jacobian(xy,r3) ! consider as public
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



      subroutine dcart2dxy(xy, r2)
c     ------------------------------------------ 
c     Convert inplace a Cartesian displacement 
c     vector r2 (meters) at xy  to a displacement 
c     vector in (longitude,latitude,depth), with
c     Cartesian space oriented along (longitude,latitude)
c     tangent space
c     ------------------------------------------ 
      real, intent(in)     :: xy(:)  !dimension 2+
      real, intent(inout)  :: r2(:)  !dimension 2+
      real                 :: jac(3)
c     ------------------------------------------  
      call get_jacobian(xy, jac)
      r2(1:2) = r2(1:2)/jac(1:2) ! element-by-element
      end subroutine 


      subroutine dxy2dcart(xy, r2)
c     ------------------------------------------ 
c     Convert inplace a (longitude,latitude,depth)
c     displacement into a Cartesian displacement 
c     r2 (meters), with Cartesian space oriented 
c     along (longitude,latitude) tangent space
c     ------------------------------------------ 
      real, intent(in)     :: xy(:)  !dimension 2+
      real, intent(inout)  :: r2(:)  !dimension 2+
      real                 :: jac(3)
c     ------------------------------------------ 
      call get_jacobian(xy, jac)
      r2(1:2) = r2(1:2)*jac(1:2) ! element-by-element
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


      subroutine get_horizontal_distance(geo1,geo2,d)
c-------------------------------------------------------
c     Calculate distance between two points on the surface 
c     of earth using the great circle approximation.
c     Input is lon/lat positions (possibly depth, but is not required)
c     Return value, d, is distance in metres
c     Earths radius used here is from constants.f
c     Note that as the distances involved are relatively
c     small, and therefore only a fraction of a radian of
c     the earths circumference, a high degree of accuracy is
c     required for dealing with the cosines. Hence, the
c     variables used here are real(16), as opposed to the
c     bog standard real(4) variables
c-------------------------------------------------------
      real, intent(in)     :: geo1(:),geo2(:)
      real, intent(out)    :: d
      real(8)             :: cosdR1,cosdR2,lat1,lat2,lon1,lon2
c     --------------------------------------------------
c     Extract values and convert to radians
      lon1=geo1(1)*deg2Rad
      lat1=geo1(2)*deg2Rad
      lon2=geo2(1)*deg2Rad
      lat2=geo2(2)*deg2Rad
c     Calculate distance      
      cosdR1 =sin(lat1)*sin(lat2)
      cosdR2 =cos(lat1)*cos(lat2)*cos(lon1-lon2)
      d=earth_radius*acos(cosdR1+cosdR2)
      end subroutine get_horizontal_distance


c      subroutine get_horizontal_distance(xy1, xy2, r)
cc     ------------------------------------------
c      Alternative formula - check close range performance
c      i.e. xy1 ~ xy2
cc     ------------------------------------------
c      real, intent(in)     :: xy1(:), xy2(:)
c      real, intent(out)    :: r
c      real :: costheta, l1,l2,p1,p2, angle
cc     ------------------------------------------ 
c      p1 = xy1(1)*deg2rad
c      p2 = xy2(1)*deg2rad
c      l1 = xy1(2)*deg2rad
c      l2 = xy2(2)*deg2rad
c      costheta = cos(p1)*cos(p2) + sin(p1)*sin(p2)
c      costheta = costheta*cos(l1)*cos(l2) + sin(l1)*sin(l2)
c      angle    = acos(costheta) ! 0 < angle < pi
c      r        = angle*earth_radius
c      end subroutine 



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
      call dcart2dxy(xy,dxy)   ! inplace transform of dxy
      xy = xy + dxy
      
      end subroutine 


      end module
