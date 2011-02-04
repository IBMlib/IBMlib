      module horizontal_grid_transformations
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Provider of horizontal grid transformations 
c     for any grid which is a mesh, which maps onto a
c     square lattice grid space
c     Coordinate conventions for grid space: 
c       * mesh points at integer value grid coordinates, 
c         first mesh point at (1,1) in grid space
c       * corner centered cells are those given by the mesh
c       * node centered cells has centers at mesh points
c         and cell faces at half-valued lines in grid space
c
c     Below is the regular lon-lat grid implementation
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      private
      integer, parameter,public :: nx = 20,ny = 20
      real, parameter    :: lambda1 = -4.0
      real, parameter    :: dlambda =  1.0
      real, parameter    :: phi1    = 49.0
      real, parameter    :: dphi    =  1.0

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


      end module horizontal_grid_transformations


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     paste-in auxillaries
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module geometry
      
      contains

      subroutine cross_2Dlines(x0,x1,y0,y1,s,t,cross) ! paste-in dummy from geometry.f
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
      real                 :: vx(2),vy(2),lvx2,lvy2,det
      real,parameter       :: s_parallel = 1.e20
      real,parameter       :: s_undef    = 2.e20
c     ----------------------------------------------------------
      vx    = x1(1:2)-x0(1:2)
      vy    = y1(1:2)-y0(1:2)
      lvx2  = sum(vx*vx)
      lvy2  = sum(vy*vy)
      if ((lvx2 < 1.e-12).or.(lvy2 < 1.e-12)) then ! angle undef
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

      end module geometry
      

      include "../generic_elements/coastal_handler_hashtables"

      module physical_fields
      use horizontal_grid_transformations
      use coastal_handler
      
      contains

      LOGICAL function is_land(geo)        ! COPY2 -> regular_lonlat_grid.f
c     ------------------------------------------ 
c     return is_land = .false. at horizontal range violation
c     ------------------------------------------ 
      real, intent(in) :: geo(:)
      real             :: xy(2)
      integer          :: ixc,iyc
c     ------------------------------------------ 
c     avoid probing wetmask at range violation

      if (horizontal_range_check(geo)) then ! range OK
         call get_horiz_ncc(geo,ixc,iyc)
         call get_horiz_grid_coordinates(geo,xy)
         is_land = ((wetmask(ixc,iyc)==0).or.at_coast_line(xy))
      else
         is_land = .false.
      endif    
c     ------------------------------------------
      end function

      end module physical_fields



      module particle_tracking
c           

      use physical_fields
      integer,parameter,private  :: verbose = 0

      contains 
c     >>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<
c     >>>>>>>>>> particle_tracking.f:begin    <<<<<<<<<<<<
c     >>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<


      subroutine multiple_reflection_path(s0, s1, anycross, sref, shit1)  ! COPY2 ->  particle_tracking.f 
c-------------------------------------------------------------------------
c     Compute key points of multiple horizontal coastal reflection path
c     by iterative application of coast_line_intersection primitive
c     when trying to move from (valid, wet) s0 to new position s1
c     If a coast line is crossed retrun anycross == .true.
c     Then sref is the (multiple) horizontal coastal reflected end point
c     which is valid and wet. shit1 is the point of (first) coastal intersection
c     If anycross == .false. (sref, shit1) is undefined
c-------------------------------------------------------------------------
      real,intent(in)      :: s0(3), s1(3)
      logical,intent(out)  :: anycross
      real,intent(out)     :: sref(3), shit1(3)

      real                 :: start(3), sfin(3), shit(3)
      integer              :: iter
      integer, parameter   :: max_reflections = 5 ! otherwise time step too long ...
c---------------------------------------------------------------
      if (verbose>0) write(*,*) " multiple reflection analysis begin"

      iter  = 0
      call coast_line_intersection(s0, s1, anycross, sref, shit1)

      if (verbose>0) then
         if (anycross) then 
            write(*,422) iter, s0(1:2), s1(1:2), sref(1:2), shit1(1:2)  
            write(*,*) "is_land(reflected point) = ", is_land(sref)
            write(*,*) "continue multiple reflection analysis"
         else
            write(*,423) iter, s0(1:2), s1(1:2)
            write(*,*)"wet reflection: multiple reflection analysis end"
         endif
      endif ! verbose

      if (.not.anycross) return
      
c
c     --- we hit the coast line, start multiple reflection analysis
c         (this means (sref, shit1) are defined and anycross == .true.)
c
      iter  = 1
      shit  = shit1  ! save first coastal hit
c     note: is_land = .false. at horizontal range violation
c     The loop below computes the final reflection, sref
c     In most cases we will not enter the while loop, because sref is a wet point
c
      do while (anycross .and. 
     +          is_land(sref).and.
     +          (iter <= max_reflections))
         start = shit
         sfin  = sref
         call coast_line_intersection(start, sfin, anycross, sref, shit)
         if (.not.anycross) sref=sfin ! roll back, in case sref is overwritten
         if (verbose>0) then
            if (anycross) then 
               write(*,422) iter, start(1:2), sfin(1:2), 
     +                      sref(1:2), shit(1:2)
               write(*,*) "is_land(reflected point) = ", is_land(sref)
            else
               write(*,423) iter, start(1:2), sfin(1:2)
            endif
         endif ! verbose

         iter  = iter + 1
      enddo
      if (iter > max_reflections) then
         write(*,*) "multiple_reflection_path: max_reflections exceeded"
         stop
      endif
 
      if (verbose>0) write(*,*) " multiple reflection analysis end"

 422  format("step ", i3, " :", 2f7.2, " ->", 2f7.2, " : ref=",
     +        2f7.2, " hit=", 2f7.2)    
 423  format("step ", i3, " :", 2f7.2, " ->", 2f7.2, " : no cross")
c-----------------------------------------------------------
      end subroutine

c     >>>>>>>>>> particle_tracking.f:end      <<<<<<<<<<<<    
      end module particle_tracking


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      program test_it
c     ------------------------------------------------------------------------
c     Compile & run:  ifort -e90 coast_line_intersection_hash_tables_general_mesh.f; a.out
c     ------------------------------------------------------------------------
c     Plot results in R       
c     dat <- read.table("fort.10");
c     plot(NA,asp=1,xlim=range(dat[,c(1,3,5,7)]),ylim=range(dat[,c(2,4,6,8)]));
c     segments(dat$V1,dat$V2,dat$V3,dat$V4);
c     segments(dat$V5,dat$V6,dat$V7,dat$V8);
c     abline(v=-180:180+0.5,h=-90:90+0.5,col="lightgrey",lty=3)
c     ----------------------------------------------------------------
      use physical_fields
      use particle_tracking
      implicit none
      real    :: geo1(3), geo2(3), georef(3), geohit(3),deg2rad
      real    :: dgeo(3),x,y,s
      logical :: anycross,ldum
      real    :: angle,r2(2)
      integer :: ixc,iyc,ix,iy,istep,ityp
      real,parameter :: rjump = 0.5  ! Bee jump range
      real    :: stability
c     ----------------------------------------------------------------

c$$$c     ------------------------ 
c$$$c     ---- MPA test block ----
c$$$c     ------------------------ 
c$$$      geo1 = (/2.0, 53.0, 0.0/)
c$$$      call get_horiz_ncc(geo1,ixc,iyc)     
c$$$      !Surround the point with dry cells
c$$$      wetmask   = 0 ! 0 dry, 1 wet
c$$$      wetmask(ixc,iyc) = 1  ! 
c$$$      !Test multiple reflection code
c$$$      deg2rad = 3.14159/180
c$$$      angle = 0 
c$$$      geo2 = geo1 + 0.5*
c$$$     +             (/cos(angle*deg2rad), sin(angle*deg2rad),0.0/)
c$$$      call multiple_reflection_path(geo1,geo2,anycross,georef,geohit)

c$$$c     ------------------------------- 
c$$$c     ---- cross coast line test ----
c$$$c     ------------------------------- 
c$$$      ixc = 10
c$$$      iyc = 10 q
c$$$      wetmask  = 0 ! 0 dry, 1 wet
c$$$      wetmask(ixc,iyc)     = 1  ! 0 dry, 1 wet 
c$$$c      do s = ixc-0.501, ixc-0.499, 1.0e-6
c$$$      do s = ixc+0.499, ixc+0.501, 1.0e-6
c$$$         x = ixc
c$$$         y = s
c$$$         geo1 = (/x2lon(x), y2lat(y), 0.0/)
c$$$         geo2 = geo1
c$$$         call repel_from_coast_line(geo2)
c$$$         write(*,*) s, sum(geo2-geo1)
c$$$      enddo
      
c     ------------------ 
c     ---- Bee test ----
c     ------------------
      geo1 = (/2.0, 53.0, 0.0/)
      stability = 6   ! log10 of number of steps
      call get_horiz_ncc(geo1,ixc,iyc)     
      !Surround the point with dry cells
      wetmask   = 0 ! 0 dry, 1 wet
      wetmask(ixc,iyc)     = 1  ! 
      wetmask(ixc+1,iyc)   = 1  ! 
      wetmask(ixc+1,iyc+1) = 1  ! 
      call configure_coastal_hash_tables()


      do istep = 1, nint(10**stability)
         call random_number(r2)
         angle = 8*atan(1.0)*r2(2)
         dgeo  = r2(1)*rjump*(/cos(angle), sin(angle), 0.0/)
         geo2  = geo1+dgeo
         call multiple_reflection_path(geo1,geo2,anycross,georef,geohit)
         if (anycross) geo2 = georef
c
c        -------- write step --------
c
c$$$         if (anycross) then
c$$$            write(*,112) istep, geo1(1:2),geo2(1:2),
c$$$     +                   geo1(1:2)+dgeo(1:2),georef(1:2),geohit(1:2)
c$$$         else
c$$$            write(*,113) istep, geo1(1:2), geo2(1:2)
c$$$         endif
c$$$ 112     format(i8," g1=", 2f10.6,"  g2=",2f10.6,
c$$$     +             " (g1+dg=", 2f10.6,
c$$$     +             " ref=",    2f10.6,
c$$$     +             " hit=",    2f10.6, ")")
c$$$ 113     format(i8," g1=", 2f10.6,"  g2=",2f10.6, " (no reflections)")

c
c        -------- check step --------
c
         if (is_land(geo2)) then
            if (anycross) then
               write(*,*) "geo2 is dry (crossed coast)"
            else
               write(*,*) "geo2 is dry (no coast cross)"
            endif
            stop
         endif

         write(65,*) geo1(1:2)
c         write(65,*) geo1(1:2),istep
c         if (is_land(geo1)) then
c            write(*,*) "geo1 dry at step",istep
c            write(*,*) "geo1 = ", geo1
c            write(*,*) geo1(1:2)
c            stop
c         endif

c
c        -------- roll state --------
c
         geo1 = geo2 ! accept step
         
      enddo



      end program

      
