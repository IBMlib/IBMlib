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




      end module





cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Implementation of coast_line_intersection based on 
c     hash tables for a static coastal geometry.
c     Appears to be somewhat faster than other version which
c     features runtime analysis of local coastal topology
c     ASC, 14Dec2010
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module test_coast_line_intersection
      use horizontal_grid_transformations
      implicit none
c      private     
    
      integer            :: wetmask(nx,ny) = 0 ! (1: wet nodes, 0 dry)  
  
      real, parameter    :: facial_tol(2) = 1.0e-5 ! note: in grid units lon part should be lat dependent
      integer,parameter  :: verbose = 0

c     ----- static coast line auxillary tables for fast classification -----
      real    :: vertical_coastlines(nx+1,ny+1)   ! longitude in degrees
      real    :: horizontal_coastlines(nx+1,ny+1) ! latitude in degrees
      real    :: corners_coastlines(nx+1,ny+1,2)  ! (longitude, latitude) in degrees
      logical :: is_coastal(nx+1,ny+1,3)          ! ityp: 1=west, 2=south,3=sw corner  


      contains 

      
      subroutine configure_coastal_hash_tables()  ! COPY2 -> regular_lonlat_grid.f
c     ------------------------------------------ 
c     Analyze coastal configuration based on wetmask and
c     generate static coastal hash tables
c     Cell index association convention:
c       parent cells run over (1,1) <= (ix,iy) <= (nx+1,ny+1)  
c       ityp refers to: 1=west, 2=south,3=sw corner of parent cell  
c     Parent cells are node centered    
c     ------------------------------------------ 
      integer :: ix, iy, fac, fac2(2), wsum, ixf, iyf
      real    :: x, y, xy(2)
      logical :: picked_a_cell
c     ------------------------------------------ 
c     ... set pad/default values ...
      vertical_coastlines   = -999.
      horizontal_coastlines = -999.
      corners_coastlines    = -999.
      is_coastal            = .false.
c
c     Scan over inner vertical faces (type 1). 
c     Set is_coastal(2:nx, 1:ny, 1) and vertical_coastlines(2:nx,1:ny)
c      
      do ix=2,nx
         x = ix-0.5             ! grid coordinate of west face of this cell 
         do iy=1,ny
            if (sum(wetmask(ix-1:ix,iy))==1) then  ! wet/dry junction -> wsum == 1
               is_coastal(ix,iy,1) = .true.
               if (wetmask(ix,iy)==1) then
                  fac =  2
               else
                  fac = -2
               endif
               vertical_coastlines(ix,iy) = x + fac*facial_tol(1)
            endif
         enddo
      enddo
c
c     Scan over inner horizontal faces (type 2). 
c     Set is_coastal(1:nx, 2:ny, 2) and horizontal_coastlines(1:nx, 2:ny)
c  
      do iy=2,ny
         y = iy-0.5     ! grid coordinate of south face of this cell
         do ix=1,nx
            if (sum(wetmask(ix,iy-1:iy))==1) then ! wet/dry junction -> wsum == 1
               is_coastal(ix,iy,2) = .true.
               if (wetmask(ix,iy)==1) then
                  fac =  2
               else
                  fac = -2
               endif
               horizontal_coastlines(ix,iy) = y + fac*facial_tol(2)
            endif
         enddo
      enddo
c
c     Scan over corners (type 3). 
c     Set is_coastal(2:nx, 2:ny, 3) and corners_coastlines(2:nx, 2:ny)
c  
      do iy=2,ny
         xy(2) = iy-0.5     ! grid coordinate of south face of this cell
         do ix=2,nx
            xy(1) = ix-0.5  ! grid coordinate of west face of this cell
            wsum = sum(wetmask(ix-1:ix,iy-1:iy))
            if ((wsum>0).and.(wsum<4)) then ! coastal condition: mix of wet/dry in cluster
               is_coastal(ix,iy,3) = .true. 
c              --- just pick first encountered wet cell
               picked_a_cell = .false.
               do ixf=ix-1,ix
                  do iyf=iy-1,iy
                     if (wetmask(ixf,iyf)==1) then
                        fac2(1) = 4*(ixf-xy(1))
                        fac2(2) = 4*(iyf-xy(2))
                        picked_a_cell = .true.
                     endif
                  enddo
               enddo
               if (.not.picked_a_cell) then
                  write(*,*) "configure_coastal_hash_tables: unexpected"
                  stop 
               endif
               corners_coastlines(ix,iy,:) = xy + fac2*facial_tol
            endif ! wsum
         enddo
      enddo
c
c     Set vertical domain bounds vertical faces (type 1). 
c  
      is_coastal(1,    1:ny, 1) = .true. ! flag left  domain bound as coastal
      is_coastal(nx+1, 1:ny, 1) = .true. ! flag right domain bound as coastal
      vertical_coastlines(1,    1:ny) =      0.5 + facial_tol(1)
      vertical_coastlines(nx+1, 1:ny) = nx + 0.5 - facial_tol(1)
c     --- grid nodes ----
      is_coastal(1,    1:ny, 3) = .true. ! flag left  domain bound as coastal
      is_coastal(nx+1, 1:ny, 3) = .true. ! flag right domain bound as coastal
      do iy=2,ny
c        --- left side ---
         xy(2) = iy-0.5     ! grid coordinate of south face of this cell
         xy(1) = 0.5  
         fac2(1) = 2 ! point to right side
         fac2(2) = 0 ! default
         if      (wetmask(1,iy)==1) then
            fac2(2) =  2
         elseif  (wetmask(1,iy)==1) then
            fac2(2) = -2
         endif ! no else clause
         corners_coastlines(1,iy,:) = xy + fac2*facial_tol
c        --- right side ---
         xy(1) = nx + 0.5 ! same xy(2) 
         fac2(1) = -2     ! point to left side
         fac2(2) =  0     ! default
         if      (wetmask(nx,iy)==1) then
            fac2(2) =  2
         elseif  (wetmask(nx,iy)==1) then
            fac2(2) = -2
         endif ! no else clause
         corners_coastlines(nx+1,iy,:) = xy + fac2*facial_tol
c
      enddo
c
c     Set horizontal domain bounds vertical faces (type 2). 
c  
      is_coastal(1:nx, 1,    2) = .true. ! flag lower domain bound as coastal
      is_coastal(1:nx, ny+1, 2) = .true. ! flag upper domain bound as coastal
      horizontal_coastlines(1:nx, 1   ) =    0.5+facial_tol(2)
      horizontal_coastlines(1:nx, ny+1) = ny+0.5-facial_tol(2) 
c     --- grid nodes ----
      is_coastal(1:nx,    1, 3) = .true. ! flag left  domain bound as coastal
      is_coastal(1:nx, ny+1, 3) = .true. ! flag right domain bound as coastal
      do ix=2,nx
c        --- lower side ---
         xy(1) = ix-0.5     ! grid coordinate of south face of this cell
         xy(2) = 0.5  
         fac2(2) = 2 ! point to upper side
         fac2(1) = 0 ! default
         if      (wetmask(ix,1)==1) then
            fac2(1) =  2
         elseif  (wetmask(ix,1)==1) then
            fac2(1) = -2
         endif ! no else clause
         corners_coastlines(ix,1,:) = xy + fac2*facial_tol
c        --- upper side ---
         xy(2) = ny + 0.5 ! same xy(1) 
         fac2(2) = -2     ! point to lower side
         fac2(1) =  0     ! default
         if      (wetmask(ix,ny)==1) then
            fac2(1) =  2
         elseif  (wetmask(ix,ny)==1) then
            fac2(1) = -2
         endif ! no else clause
         corners_coastlines(ix,ny+1,:) = xy + fac2*facial_tol
c
      enddo
c
c     Finally set domain bound corners (type 3). 
c    
      is_coastal(1, 1, 3)    = .true. ! sw corner
      fac2 = (/2,   2/)
      xy   = (/0.5, 0.5/)
      corners_coastlines(1,1,:)       = xy + fac2*facial_tol
c
      is_coastal(nx+1, 1, 3) = .true. ! se corner
      fac2 = (/-2,     2/)
      xy   = (/nx+0.5, 0.5/)
      corners_coastlines(nx+1, 1,:)   = xy + fac2*facial_tol
c
      is_coastal(1, ny+1, 3) = .true. ! nw corner
      fac2 = (/2,   -2/)
      xy   = (/0.5, ny+0.5/)
      corners_coastlines(1,ny+1,:)    = xy + fac2*facial_tol
c
      is_coastal(nx+1, ny+1, 3) = .true. ! ne corner
      fac2 = (/-2,   -2/)
      xy   = (/nx+0.5, ny+0.5/)
      corners_coastlines(nx+1,ny+1,:) = xy + fac2*facial_tol

      end subroutine configure_coastal_hash_tables





      subroutine coast_line_intersection(geo1, geo2,cross_coast,georef,  ! COPY2 -> regular_lonlat_grid.f
     +                                   geohit) 
c     ------------------------------------------ 
c     Trajectory tracking variant
c
c     This function is a service function assisting in enforcing 
c     coastal boundary conditions on particle steps. geo1 is the current 
c     valid (wet, inside domain) position. geo2 is the next position that may or may not be in water.
c     The assumption is that the particle moves in a straight line 
c     between geo1 and geo2 in (lon,lat,z) coordinates. This algorithm 
c     solves the problems exactly for the any coastal geometry represented
c     as regular cell centered squares (without previous assumptions about
c     small steps geo1 -> geo2, so that any number of neighbor cells can be crossed)

c     Output: 
c     * anycross: If the straight line geo1 -> geo2 has crossed a coast line, 
c       anycross is returned .true. otherwise .false. (i.e. the line between geo1 and 
c        geo2 is all in water). 
c
c     If anycross is true, the following vectors are computed:
c
c     * georef is the modified final position, if step from geo1 to geo2 is reflected
c       in the coast line
c     * geohit is the first position where the line from geo1 to geo2 
c       crosses a coast line. 
c     Both (geohit,georef) are repelled from the coastline, if the they end there
c     to within a numerical tolerance, by a small displacement
c     Therefore geohit should test wet
c
c     If anycross is false neither georef and geohit are assigned.
c     If geo1 is dry, a warning is issued and anycross is returned false
c     In this implementation does not handle multiple reflection effects,
c     i.e. it does not test that georef is a wet point (the step steps geo1 -> geo2
c     may become reflected onto land in concave coastal geometries for large steps)
c     geo1 and geo2 must have same length (2 or 3) so they are vector addable
c     Currently it is being asserted that geo2 does not violate the outer domain (not checked)
c     
c     Conceptual revision ASC/MPA Nov 24, 2010 
c     Infinitesimal coast line repulsion added ASC Dec 06, 2010 
c     Resolve association conflict at cell transition analysis, ASC Jan 18, 2011 
c     Internal grid space implementation, ASC Feb 02, 2011  
c     ---------------------------------------------------- 
      real, intent(in)     :: geo1(:),geo2(:)
      logical, intent(out) :: cross_coast
      real, intent(out)    :: georef(:), geohit(:)
c   
      real                 :: xy1(size(geo1)),xy2(size(geo2))         ! grid coor set
      real                 :: xyref(size(georef)),xyhit(size(geohit)) ! grid coor set

      real                 :: s
      logical              :: leaving,loopflag,ldum
      character*1          :: reftype, xface
      real                 :: direct(size(geo1)), reflect(size(geo1))
      integer              :: ix,iy,ix2,iy2,istep,maxsteps
c     ---------------------------------------------------- 
      call get_horiz_ncc(geo1,ix,iy)    ! resolve start point cell
      call get_horiz_ncc(geo2,ix2,iy2)  ! resolve end point cell
c
c     Test assertion that geo1 is wet
c
      if (wetmask(ix,iy)==0) then
         write(*,347) geo1(1:2)        
 347     format("coast_line_intersection: warning: geo1=",
     +          2f10.6, " is a dry point")

         cross_coast = .false.
         return
      endif
c
c     Test if geo1 and geo2 are in same cell => no coast line crossing
c     most cases terminate here
c
      if ((ix==ix2).and.(iy==iy2)) then
         if (verbose>0) write(*,*) "geo1 -> geo2 in same cell"
         call get_horiz_grid_coordinates(geo2,xy2)
         cross_coast = check_coastal_osculation(xy2,xyref,xyhit)
         call get_horiz_geo_coordinates(xyref,georef) ! return in geo coordinates
         call get_horiz_geo_coordinates(xyhit,geohit) ! return in geo coordinates
         if ((verbose>0).and.cross_coast) 
     +       write(*,*) "coastal osculation of geo2 handled"
         return
      endif
      if (verbose>0) write(*,*) "trajectory analysis: begin"

c     Now we know that we are leaving the cell. But where? 
c     Follow sequence moving from one cell to the next until we either
c     reach the end of the vector (s=1), or encounter land
c     Interior part here is in grid coordinates

      call get_horiz_grid_coordinates(geo1,xy1)
      call get_horiz_grid_coordinates(geo2,xy2)
      cross_coast = .false. ! default unless detected below
      loopflag    = .true.
      istep       = 1
      maxsteps    = max(1,abs(ix-ix2))*max(1,abs(iy-iy2)) ! safety plug

      do while(loopflag)
         call assess_cell_crossing(xy1,xy2,ix,iy,s,xface)
         if (verbose>0) write(*,421) ix,iy,xface
 421     format("leaving cell ",2i4," via ", a," face") 
         !Move into the next cell
         select case(xface)
         case ("N")
            iy = iy +1   
         case ("S")
            iy = iy -1
         case ("E")
            ix = ix +1
         case ("W")
            ix = ix -1
         case default
            write (*,*) "coast_line_intersection: Exit face error"
            stop
         end select
         !Is the next cell wet or dry? If its dry, then we've hit the coast.  
         !Move into the next cell
         if (verbose>0) then
            if(wetmask(ix,iy)>0) then
               write(*,422) ix,iy,"wet"
            else
               write(*,422) ix,iy,"dry" 
            endif
422        format("new cell : ",2i4," is ", a) 
         endif

         if(wetmask(ix,iy) == 0) then     ! have crossed onto land
            cross_coast = .true.          ! so setup for exit
            loopflag    = .false.        
         elseif ((ix==ix2).and.(iy==iy2)) then ! we are in final cell without crossing the coast  
            loopflag    = .false.         ! so setup for exit, still cross_coast == .false.
         elseif (istep>maxsteps) then
            write(*,*) "coast_line_intersection: maxsteps exceeded)"
            stop
         endif

         istep = istep + 1

      enddo
c  
c     if (cross_coast == .true.) continue and resolve xyref and xyhit
c     else let xyref and xyhit remain undefined
c
      if (.not.cross_coast) then
         cross_coast = check_coastal_osculation(xy2,xyref,xyhit) ! xy2 resolved at loop entry
         call get_horiz_geo_coordinates(xyref,georef) ! return in geo coordinates
         call get_horiz_geo_coordinates(xyhit,geohit) ! return in geo coordinates
         if ((verbose>0).and.cross_coast) 
     +       write(*,*) "coastal osculation of geo2 handled"
         return
      endif
c
c     If the coast line is crossed, compute xyref, xyhit
c     0 < s < 1 is the coordinate along the vector (xy2-xy1) 
c     corresponding to the point xyhit
c     
      direct  = xy2-xy1
      reflect = direct
      if     (xface=="N" .or. xface =="S") then
         reflect(2) = -reflect(2) ! horizontal box line reflection
      elseif (xface=="E" .or. xface =="W") then
         reflect(1) = -reflect(1) ! vertical   box line reflection
      else
         stop "coast_line_intersection: unhandled reflection type"
      endif
c     |direct| == |reflect|        
      xyhit = xy1   +     s*direct ! position where coast line is hit
      xyref = xyhit + (1-s)*reflect ! position for reflection in coast line
 
c     Repel both (xyhit,xyref) from the coastline, if the they end there
c     to within a numerical tolerance, by a small displacement

      ldum = repel_from_coast_line(xyhit)   ! do not parse result
      ldum = repel_from_coast_line(xyref)   ! do not parse result
      call get_horiz_geo_coordinates(xyref,georef) ! return in geo coordinates
      call get_horiz_geo_coordinates(xyhit,geohit) ! return in geo coordinates

      return   ! from coast_line_intersection

c     -----------------------------
      contains  ! local subroutines
c     -----------------------------

      subroutine assess_cell_crossing (xy1,xy2,ix,iy,s,exitface) ! internal to coast_line_intersection
c     --------------------------------------------------------------
c     Determine where the direct line from xy1 through xy2 
c     leaves the cell with center (ix,iy)
c     Return 0 < s  which is the coordinate along the vector (xy2-xy1) 
c     where the vector exits the cell. Notice s may exceed 1 
c     (which is the case, if cell (ix,iy) contains xy2)
c     Also identify the exit face as one of "N","S","E","W" 
c     so the proper reflection can be calculated
c
c     Modified from assess_cell_leaving to avoid decision conflicts ASC/18Jan2011
c     --------------------------------------------------------------
      real, intent(in)         :: xy1(:),xy2(:)
      integer, intent(in)      :: ix,iy
      real, intent(out)        :: s
      character*1, intent(out) :: exitface
c
      logical                  :: yes
      real                     :: stest,tdum,dxy(2)
      real                     :: c00(2),c01(2),c10(2),c11(2)
c     --------------------------------------------------------------
      call get_horiz_ncc_corners(ix,iy,c00,c01,c10,c11) ! resolve faces          
c
c     Determine how the step leaves the box defined by (c00,c01,c10,c11))
c     We can make efficiency gains by taking advantage of the direction
c     that the particle is travelling in - a particle travelling NE
c     can only leave by the North and East faces etc
c
      exitface = "@"    ! we should not pick this one up now, but ...
      s        = 1.0e12 ! must be set before first comparison
      dxy     = xy2(1:2)-xy1(1:2)
c
c     First the east-west direction
c
      if(dxy(1)>0) then    !we're moving to the east

         call cross_2Dlines(xy1,xy2,c10,c11,stest,tdum,yes)
         if (yes.and.(stest<s)) then ! rely on right-to-left evaluation
            s         = stest
            exitface  = "E"   
         endif
      else   !we're moving to the west
         call cross_2Dlines(xy1,xy2,c00,c01,stest,tdum,yes)
         if (yes.and.(stest<s)) then ! rely on right-to-left evaluation
            s         = stest
            exitface  = "W"   
         endif
      endif
c
c     Then the north-south direction
c
      if(dxy(2) >0) then !we're moving to the North
         call cross_2Dlines(xy1,xy2,c01,c11,stest,tdum,yes)
         if (yes.and.(stest<s)) then ! rely on right-to-left evaluation
            s         = stest
            exitface  = "N"   
         endif
      else  !we're moving to the south
         call cross_2Dlines(xy1,xy2,c00,c10,stest,tdum,yes)
         if (yes.and.(stest<s)) then ! rely on right-to-left evaluation
            s         = stest
            exitface  = "S"   
         endif
      endif
c
c     check exit conditions
c
      if (s<0) then
         write(*,*) "assess_cell_crossing: unexpected s<0"
         write(*,*) "found s =",s
         write(*,*) "xy1    =",xy1
         write(*,*) "xy2    =",xy2
         write(*,*) "(ix,iy) =",ix,iy
         stop       
      endif
      end subroutine assess_cell_crossing     ! local subroutine

      end subroutine coast_line_intersection ! embedding subroutine




      logical function check_coastal_proximity(xy,ixf,iyf,ityp) ! COPY2 -> regular_lonlat_grid.f
c----------------------------------------------
c     Test whether xy is within the coastal perimeter
c     Do not modify xy
c
c     Output: check_coastal_proximity = .false.
c                -> xy is not within the coastal eprimeter
c             ixf,iyf is assigned, but not meaningful
c             ityp is not assigned
c
c             check_coastal_proximity = .true.
c                -> xy is within the coastal eprimeter
c             
c             ixf,iyf is the center indices of the parent cell of the detected coastal face/corner
c             ityp = 1: western  face of the parent cell (x = ixf - 0.5)
c             ityp = 2: southern face of the parent cell (y = iyf - 0.5)
c             ityp = 3: south-western corner of the parent cell 
c                       (x,y) = (ixf - 0.5, iyf - 0.5)       
c----------------------------------------------
      real, intent(in)      :: xy(:)
      integer, intent(out)  :: ixf,iyf,ityp   
      logical               :: at_x_face, at_y_face, at_corner
      real                  :: xh, yh, dxh, dyh
      logical               :: potentially_at_coast
c
c     1) test whether xy is on a cell boundary (and potentially coastal)
c  
      xh = xy(1)+0.5
      yh = xy(2)+0.5
      ixf = nint(xh)     ! parent face cell
      iyf = nint(yh)     ! parent face cell
      dxh = abs(xh-ixf)
      dyh = abs(yh-iyf)
      at_x_face = .false.    ! default case
      at_y_face = .false.    ! default case
      at_corner   = .false.    ! default case
      if (dxh < facial_tol(1))     at_x_face   = .true. ! should be latitude dependent ...
      if (dyh < facial_tol(2))     at_y_face   = .true.
      if (at_x_face.and.at_y_face) at_corner = .true.
      potentially_at_coast = at_x_face.or.at_y_face
c
c     2) Lookup in hash table is_coastal whether this element is coastal
c 
      if (potentially_at_coast) then
         if (at_corner) then  ! both at_x_face/at_y_face remain .true.
            ityp = 3          ! (ixf,iyf) gives correct look-up
         elseif (at_y_face) then
            ityp = 2
            ixf  = nint(xy(1))    ! enforce correct x axis look-up
         elseif (at_x_face) then
            ityp = 1
            iyf  = nint(xy(2))    ! enforce correct y axis look-up
         else
            stop "check_coastal_proximity: unexpected"
         endif
         check_coastal_proximity = is_coastal(ixf,iyf,ityp)   
      else
         check_coastal_proximity = .false.
         return
      endif
      end function  check_coastal_proximity
         
     

      logical function at_coast_line(xy)             ! COPY2 -> regular_lonlat_grid.f
c----------------------------------------------
c     Just return the main result whether xy
c     xy is numerically on a coast line and
c     waste detailed results from check_coastal_proximity
c     Do not modify xy 
c----------------------------------------------
      real, intent(in)    :: xy(:)
      integer             :: ixf,iyf,ityp       ! dummies for check_coastal_proximity
c----------------------------------------------      
      at_coast_line = check_coastal_proximity(xy,ixf,iyf,ityp)      
      end function at_coast_line


      
      logical function repel_from_coast_line(xy)   !  COPY2 -> regular_lonlat_grid.f
c----------------------------------------------
c     Internal function that places xy(1:2) outside 
c     the coastal perimeter, if tested inside. 
c----------------------------------------------
      real, intent(inout) :: xy(:)
      integer             :: ixf,iyf,ityp
c----------------------------------------------
      repel_from_coast_line = check_coastal_proximity(xy,ixf,iyf,ityp) 

      if (.not.repel_from_coast_line) return
c
c     At this point we know that xy has been flagged as 
c     coastal point; repel it
c      
      if     (ityp==1) then
         xy(1) = vertical_coastlines(ixf,iyf)
      elseif (ityp==2) then
         xy(2) = horizontal_coastlines(ixf,iyf)
      elseif (ityp==3) then
         xy(1:2) = corners_coastlines(ixf,iyf,:)
      else
         stop "repel_from_coast_line: unexpected"
      endif

      end function repel_from_coast_line
      


      
      logical function check_coastal_osculation(xy,xyref,xyhit) ! COPY2 -> regular_lonlat_grid.f
c     ----------------------------------------------------------- 
c     Auxillary function that chacks whether xy is in the coastal proximity
c     i.e. at the coastal line within a predefined tolerence
c
c     Output:   check_coastal_osculation = .false.
c                  xy is not in the coastal proximity
c                  xy unchanged
c                  xyref = xy
c                  xyhit = xy
c
c               check_coastal_osculation = .true.
c                  xy is in the coastal proximity
c                  xy unchanged
c                  xyref: is the position xy displaced infinitesimally
c                          outside the coastal proximity detection
c                  xyhit: a copy of input xy (reciding in the coastal proximity)
c     ----------------------------------------------------------- 
      real, intent(in)  :: xy(:)
      real, intent(out) :: xyref(:), xyhit(:)    
c     -----------------------------------------------------------
      xyhit = xy
      xyref = xy
      check_coastal_osculation = repel_from_coast_line(xyref) 
      end function check_coastal_osculation


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



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     paste-in auxillaries
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      

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


      end module

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
      use horizontal_grid_transformations
      use test_coast_line_intersection
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
c$$$      iyc = 10 
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

      
