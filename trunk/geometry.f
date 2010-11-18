c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Geometry module
c     ---------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
c     Provides basic geography functions and transformations
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      module geometry
      use constants

      implicit none
      private

      public ::  cross_2Dline_segments
      public ::  dcart2dxy
      public ::  get_horizontal_distance

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
c     Corrected 28 Oct 2010 to ensure only interior solution
c     ----------------------------------------------------------
      implicit none
      real, intent(in)     :: x0(*),x1(*),y0(*),y1(*)
      real, intent(out)    :: s
      logical, intent(out) :: cross
      real                 :: vx(2),vy(2),lvx2,lvy2,det,t
      real,parameter       :: s_parallel = 1.e20
      real,parameter       :: s_undef    = 2.e20
c     ----------------------------------------------------------
      vx    = x1(1:2)-x0(1:2)
      vy    = y1(1:2)-y0(1:2)
      lvx2  = sum(vx*vx)
      lvy2  = sum(vy*vy)
      if ((lvx2 < 1.e-12).or.(lvy2 < 1.e-12)) then ! angle undef
         s = s_undef
         cross=.false.
         return
      endif   
      det = vx(2)*vy(1) - vx(1)*vy(2)     
      if (abs(det)< 1.e-8) then ! lines parallel
         s = s_parallel
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

      s = (vy(2)*x0(1) - vy(1)*x0(2) - vy(2)*y0(1) + vy(1)*y0(2))/det
      t = (vx(2)*x0(1) - vx(1)*x0(2) - vx(2)*y0(1) + vx(1)*y0(2))/det

      if (in_princip_range(s).and.in_princip_range(t)) then
         cross = .true.
      else
         cross = .false.
      endif
      
      contains

      logical function in_princip_range(x)
      real, intent(in)     :: x
      in_princip_range = ((x >= 0.0).and.(x <= 1.0))
      end function
c     ----------------------------------------------------------
      end subroutine cross_2Dline_segments


      subroutine dcart2dxy(xy, r2)
c     ------------------------------------------ 
c     Convert inplace a Cartesian displacement 
c     vector r2 (meters) at xy  to a displacement 
c     vector in (longitude,latitude,depth)
c     ------------------------------------------ 
      real, intent(in)     :: xy(:)  !dimension 2+
      real, intent(inout)  :: r2(:)  !dimension 2+
      real                 :: jac(2)
c     ------------------------------------------ 
      jac(1) = earth_radius*cos(xy(2)*deg2rad)*deg2rad
      jac(2) = earth_radius*deg2rad
      r2(1:2) = r2(1:2)/jac(1:2) !element-by-element
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
      real(16)             :: cosdR1,cosdR2,lat1,lat2,lon1,lon2
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


      end module
