cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Implementation of coast_line_intersection based on 
c     hash tables for a static coastal geometry.
c     Appears to be somewhat faster than other version which
c     features runtime analysis of local coastal topology
c     ASC, 14Dec2010
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module test_coast_line_intersection
      
      implicit none

      integer, parameter :: nx = 20,ny = 20
      integer            :: wetmask(nx,ny) = 0 ! (1: wet nodes, 0 dry)  
      real, parameter    :: lambda1 = -4.0
      real, parameter    :: dlambda =  1.0
      real, parameter    :: phi1    = 49.0
      real, parameter    :: dphi    =  1.0
      real, parameter    :: facial_tol(2) = 1.0e-4 ! note: lon part should be lat dependent
      integer,parameter  :: verbose = 0

c     ----- static coast line auxillary tables for fast classification -----
      real    :: vertical_coastlines(nx+1,ny+1)   ! longitude in degrees
      real    :: horizontal_coastlines(nx+1,ny+1) ! latitude in degrees
      real    :: corners_coastlines(nx+1,ny+1,2)  ! (longitude, latitude) in degrees
      logical :: is_coastal(nx+1,ny+1,3)          ! ityp: 1=west, 2=south,3=sw corner  

      contains 

      
      subroutine configure_coastal_hash_tables()
c     ------------------------------------------ 
c     Analyze coastal configuration based on wetmask and
c     generate static coastal hash tables
c     Cell index association convention:
c       parent cells run over (1,1) <= (ix,iy) <= (nx+1,ny+1)  
c       ityp refers to: 1=west, 2=south,3=sw corner of parent cell      
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
               vertical_coastlines(ix,iy) = x2lon(x + fac*facial_tol(1))
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
               horizontal_coastlines(ix,iy) = 
     +                            y2lat(y + fac*facial_tol(2))
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
               corners_coastlines(ix,iy,:) = 
     +                   xy2lonlat(xy + fac2*facial_tol)
            endif ! wsum
         enddo
      enddo
c
c     Set vertical domain bounds vertical faces (type 1). 
c  
      is_coastal(1,    1:ny, 1) = .true. ! flag left  domain bound as coastal
      is_coastal(nx+1, 1:ny, 1) = .true. ! flag right domain bound as coastal
      vertical_coastlines(1,    1:ny) = x2lon(     0.5 + facial_tol(1))
      vertical_coastlines(nx+1, 1:ny) = x2lon(nx + 0.5 - facial_tol(1))
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
         corners_coastlines(1,iy,:) = xy2lonlat(xy + fac2*facial_tol)
c        --- right side ---
         xy(1) = nx + 0.5 ! same xy(2) 
         fac2(1) = -2     ! point to left side
         fac2(2) =  0     ! default
         if      (wetmask(nx,iy)==1) then
            fac2(2) =  2
         elseif  (wetmask(nx,iy)==1) then
            fac2(2) = -2
         endif ! no else clause
         corners_coastlines(nx+1,iy,:) = xy2lonlat(xy + fac2*facial_tol)
c
      enddo
c
c     Set horizontal domain bounds vertical faces (type 2). 
c  
      is_coastal(1:nx, 1,    2) = .true. ! flag lower domain bound as coastal
      is_coastal(1:nx, ny+1, 2) = .true. ! flag upper domain bound as coastal
      horizontal_coastlines(1:nx, 1   ) = y2lat(   0.5+facial_tol(2))
      horizontal_coastlines(1:nx, ny+1) = y2lat(ny+0.5-facial_tol(2)) 
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
         corners_coastlines(ix,1,:) = xy2lonlat(xy + fac2*facial_tol)
c        --- upper side ---
         xy(2) = ny + 0.5 ! same xy(1) 
         fac2(2) = -2     ! point to lower side
         fac2(1) =  0     ! default
         if      (wetmask(ix,ny)==1) then
            fac2(1) =  2
         elseif  (wetmask(ix,ny)==1) then
            fac2(1) = -2
         endif ! no else clause
         corners_coastlines(ix,ny+1,:) = xy2lonlat(xy + fac2*facial_tol)
c
      enddo
c
c     Finally set domain bound corners (type 3). 
c    
      is_coastal(1, 1, 3)    = .true. ! sw corner
      fac2 = (/2,   2/)
      xy   = (/0.5, 0.5/)
      corners_coastlines(1,1,:)       = xy2lonlat(xy + fac2*facial_tol)
c
      is_coastal(nx+1, 1, 3) = .true. ! se corner
      fac2 = (/-2,     2/)
      xy   = (/nx+0.5, 0.5/)
      corners_coastlines(nx+1, 1,:)   = xy2lonlat(xy + fac2*facial_tol)
c
      is_coastal(1, ny+1, 3) = .true. ! nw corner
      fac2 = (/2,   -2/)
      xy   = (/0.5, ny+0.5/)
      corners_coastlines(1,ny+1,:)    = xy2lonlat(xy + fac2*facial_tol)
c
      is_coastal(nx+1, ny+1, 3) = .true. ! ne corner
      fac2 = (/-2,   -2/)
      xy   = (/nx+0.5, ny+0.5/)
      corners_coastlines(nx+1,ny+1,:) = xy2lonlat(xy + fac2*facial_tol)
c     -----------------------------
      contains  ! local subroutines
c     -----------------------------

      real function x2lon(x)    
c     ------------------------------------------------------
      real, intent(in)  :: x
      x2lon = lambda1 + (x-1.0)*dlambda
      end function x2lon

      real function y2lat(y) 
c     ------------------------------------------------------
      real, intent(in)  :: y
      y2lat = phi1 + (y - 1.0)*dphi 
      end function y2lat

      function xy2lonlat(xy)
c     ------------------------------------------------------
      real, intent(in)  :: xy(:)
      real              :: xy2lonlat(2)
      xy2lonlat(1) = x2lon(xy(1))
      xy2lonlat(2) = y2lat(xy(2))
      end function xy2lonlat
c
      end subroutine configure_coastal_hash_tables





      subroutine coast_line_intersection(geo1, geo2,cross_coast,georef,  ! -> regular_lonlat_grid.f
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
c     ------------------------------------------ 
      real, intent(in)     :: geo1(:),geo2(:)
      logical, intent(out) :: cross_coast
      real, intent(out)    :: georef(:), geohit(:)
c   
      real                 :: s, s_min
     
      logical              :: leaving,loopflag,ldum
      character*1          :: reftype, xface
      real                 :: direct(size(geo1)), reflect(size(geo1))
      integer              :: ix,iy,ix2,iy2
c     ------------------------------------------ 
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
         cross_coast = check_coastal_osculation(geo2,georef,geohit)
         if ((verbose>0).and.cross_coast) 
     +       write(*,*) "coastal osculation of geo2 handled"
         return
      endif
      if (verbose>0) write(*,*) "trajectory analysis: begin"

c     Now we know that we are leaving the cell. But where? 
      call assess_cell_leaving(geo1,geo2,ix,iy,leaving,s,xface) 

c     Follow sequence moving from one cell to the next until we either
c     reach the end of the vector, or encounter land
      s_min = 1.0           ! upper bound
      cross_coast = .false. ! default unless detected below
      loopflag  = .true.
      do while(loopflag)
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
         !Is the next cell wet or dry? If its dry, then we've hit the
         !coast.  
         !Move into the next cell
         if (verbose>0) then
            if(wetmask(ix,iy)>0) then
               write(*,422) ix,iy,"wet"
            else
               write(*,422) ix,iy,"dry" 
            endif
422        format("new cell : ",2i4," is ", a) 
         endif

         if(wetmask(ix,iy) == 0) then  !have crossed onto land
           cross_coast = .true.        !so setup for exit
           loopflag    = .false.        
         else   !do we leave this new cell?
           call assess_cell_leaving(geo1,geo2,ix,iy,leaving,s,xface) 
           loopflag = leaving     !loop if we leave this new cell
         endif
      enddo
c  
c     if (cross_coast == .true.) continue and resolve georef and geohit
c     else let georef and geohit remain undefined
c
      if (.not.cross_coast) then
         cross_coast = check_coastal_osculation(geo2,georef,geohit)
         if ((verbose>0).and.cross_coast) 
     +       write(*,*) "coastal osculation of geo2 handled"
         return
      endif
c
c     If the coast line is crossed, compute georef, geohit
c     0 < s_min < 1 is the coordinate along the vector (geo2-geo1)
c          
      direct  = geo2-geo1
      reflect = direct
      if     (xface=="N" .or. xface =="S") then
         reflect(2) = -reflect(2) ! horizontal box line reflection
      elseif (xface=="E" .or. xface =="W") then
         reflect(1) = -reflect(1) ! vertical   box line reflection
      else
         stop "coast_line_intersection: unhandled reflection type"
      endif
c     |direct| == |reflect|        
      geohit = geo1   +     s*direct ! position where coast line is hit
      georef = geohit + (1-s)*reflect ! position for reflection in coast line
 
c     Repel both (geohit,georef) from the coastline, if the they end there
c     to within a numerical tolerance, by a small displacement

      ldum = repel_from_coast_line(geohit)   ! do not parse result
      ldum = repel_from_coast_line(georef)   ! do not parse result

      return   ! from coast_line_intersection

c     -----------------------------
      contains  ! local subroutines
c     -----------------------------

      subroutine assess_cell_leaving(geo1,geo2,ix,iy,leave,s,exitface) ! internal to coast_line_intersection
c     --------------------------------------------------------------
c     Determine whether and where the direct line from geo1 to geo2 
c     leaves the cell with center (ix,iy)
c     If so return leave=.true. (else .false.)
c     If so return leave=.true., determine 0 < s < 1 which is the 
c     coordinate along the vector (geo2-geo1) where the vector exits the
c     cell. Also identify the exit face as one of "N","S","E","W" 
c     so the proper reflection can be calculated
c     If leave=.false. (s,exitface) are undefined
c     --------------------------------------------------------------
      real, intent(in)         :: geo1(:),geo2(:)
      integer, intent(in)      :: ix,iy
      logical, intent(out)     :: leave
      real, intent(out)        :: s
      character*1, intent(out) :: exitface
c
      logical                  :: yes
      real                     :: stest,dgeo(2)
      real                     :: c00(2),c01(2),c10(2),c11(2)
c     --------------------------------------------------------------
      call get_horiz_ncc_corners(ix,iy,c00,c01,c10,c11) ! resolve faces          
c
c     Determine how the step leaves the box defined by (c00,c01,c10,c11))
c     We can make efficiency gains by taking advantage of the direction
c     that the particle is travelling in - a particle travelling NE
c     can only leave by the North and East faces etc
c
      leave = .false.
      exitface = "@"
      s     = 1.0 
      dgeo = geo2(1:2)-geo1(1:2)
c     First the east-west direction
      if(dgeo(1)>0) then    !we're moving to the east
         call cross_2Dline_segments(geo1,geo2,c10,c11,stest,yes)
         if (yes.and.(stest<=s)) then ! rely on right-to-left evaluation
            leave = .true.
            s     = stest
            exitface  = "E"   
         endif
      else   !we're moving to the west
         call cross_2Dline_segments(geo1,geo2,c00,c01,stest,yes)
         if (yes.and.(stest<=s)) then ! rely on right-to-left evaluation
            leave = .true.
            s     = stest
            exitface  = "W"   
         endif
      endif
      if(dgeo(2) >0) then !we're moving to the North
         call cross_2Dline_segments(geo1,geo2,c01,c11,stest,yes)
         if (yes.and.(stest<=s)) then ! rely on right-to-left evaluation
            leave = .true.
            s     = stest
            exitface  = "N"   
         endif
      else  !we're moving to the south
         call cross_2Dline_segments(geo1,geo2,c00,c10,stest,yes)
         if (yes.and.(stest<=s)) then ! rely on right-to-left evaluation
            leave = .true.
            s     = stest
            exitface  = "S"   
         endif
      endif
      end subroutine assess_cell_leaving     ! local subroutine
      end subroutine coast_line_intersection ! embedding subroutine




      logical function check_coastal_proximity(geo,ixf,iyf,ityp)
c----------------------------------------------
c     Test whether geo is within the coastal perimeter
c     Do not modify geo
c
c     Output: check_coastal_proximity = .false.
c                -> geo is not within the coastal eprimeter
c             ixf,iyf is assigned, but not meaningful
c             ityp is not assigned
c
c             check_coastal_proximity = .true.
c                -> geo is within the coastal eprimeter
c             
c             ixf,iyf is the center indices of the parent cell of the detected coastal face/corner
c             ityp = 1: western  face of the parent cell (x = ixf - 0.5)
c             ityp = 2: southern face of the parent cell (y = iyf - 0.5)
c             ityp = 3: south-western corner of the parent cell 
c                       (x,y) = (ixf - 0.5, iyf - 0.5)       
c----------------------------------------------
      real, intent(in)      :: geo(:)
      integer, intent(out)  :: ixf,iyf,ityp   
      logical               :: at_lon_face, at_lat_face, at_corner
      real                  :: x, y, xh, yh, dxh, dyh
      logical               :: potentially_at_coast
c
c     1) test whether geo is on a cell boundary (and potentially coastal)
c 
      x  = 1.0 + (geo(1)-lambda1)/dlambda   ! grid coordinate
      y  = 1.0 + (geo(2)-phi1)/dphi         ! grid coordinate
      xh = x+0.5
      yh = y+0.5
      ixf = nint(xh)     ! parent face cell
      iyf = nint(yh)     ! parent face cell
      dxh = abs(xh-ixf)
      dyh = abs(yh-iyf)
      at_lon_face = .false.    ! default case
      at_lat_face = .false.    ! default case
      at_corner   = .false.    ! default case
      if (dxh < facial_tol(1))     at_lon_face   = .true. ! should be latitude dependent ...
      if (dyh < facial_tol(2))     at_lat_face   = .true.
      if (at_lon_face.and.at_lat_face) at_corner = .true.
      potentially_at_coast = at_lon_face.or.at_lat_face
c
c     2) Lookup in hash table is_coastal whether this element is coastal
c 
      if (potentially_at_coast) then
         if (at_corner) then  ! both at_lon_face/at_lat_face remain .true.
            ityp = 3          ! (ixf,iyf) gives correct look-up
         elseif (at_lat_face) then
            ityp = 2
            ixf  = nint(x)    ! enforce correct x axis look-up
         elseif (at_lon_face) then
            ityp = 1
            iyf  = nint(y)    ! enforce correct y axis look-up
         else
            stop "check_coastal_proximity: unexpected"
         endif
         check_coastal_proximity = is_coastal(ixf,iyf,ityp)   
      else
         check_coastal_proximity = .false.
         return
      endif
      end function  check_coastal_proximity
         
     

      logical function at_coast_line(geo)             ! -> regular_lonlat_grid.f
c----------------------------------------------
c     Just return the main result whether geo
c     geo is numerically on a coast line and
c     waste detailed results from check_coastal_proximity
c     Do not modify geo 
c----------------------------------------------
      real, intent(in)    :: geo(:)
      integer             :: ixf,iyf,ityp       ! dummies for check_coastal_proximity
c----------------------------------------------      
      at_coast_line = check_coastal_proximity(geo,ixf,iyf,ityp)      
      end function at_coast_line


      
      logical function repel_from_coast_line(geo)   ! -> regular_lonlat_grid.f
c----------------------------------------------
c     Internal function that places geo(1:2) outside 
c     the coastal perimeter, if tested inside. 
c----------------------------------------------
      real, intent(inout) :: geo(:)
      integer             :: ixf,iyf,ityp

c----------------------------------------------
      repel_from_coast_line = check_coastal_proximity(geo,ixf,iyf,ityp) 

      if (.not.repel_from_coast_line) return
c
c     At this point we know that geo has been flagged as 
c     coastal point; repel it
c      
      if     (ityp==1) then
         geo(1) = vertical_coastlines(ixf,iyf)
      elseif (ityp==2) then
         geo(2) = horizontal_coastlines(ixf,iyf)
      elseif (ityp==3) then
         geo(1:2) = corners_coastlines(ixf,iyf,:)
      else
         stop "repel_from_coast_line: unexpected"
      endif
      

      end function repel_from_coast_line
      


      
      logical function check_coastal_osculation(geo,georef,geohit) ! -> regular_lonlat_grid.f
c     ----------------------------------------------------------- 
c     Auxillary function that chacks whether geo is in the coastal proximity
c     i.e. at the coastal line within a predefined tolerence
c
c     Output:   check_coastal_osculation = .false.
c                  geo is not in the coastal proximity
c                  geo unchanged
c                  georef = geo
c                  geohit = geo
c
c               check_coastal_osculation = .true.
c                  geo is in the coastal proximity
c                  geo unchanged
c                  georef: is the position geo displaced infinitesimally
c                          outside the coastal proximity detection
c                  geohit: a copy of input geo (reciding in the coastal proximity)
c     ----------------------------------------------------------- 
      real, intent(in)  :: geo(:)
      real, intent(out) :: georef(:), geohit(:)    
c     -----------------------------------------------------------
      geohit = geo
      georef = geo
      check_coastal_osculation = repel_from_coast_line(georef) 
      end function check_coastal_osculation






      logical function is_land(geo)                  ! update -> regular_lonlat_grid.f
c     ------------------------------------------------------
c     Check whether geo (lon,lat,depth) corresponds to a wet position or not by
c     looking up in the wet array
c     Dec 08, 2010: include coast line in land definition
c     ASC: fixed wetmask type + grid indices lookup
c     ------------------------------------------------------
      real,    intent(in) :: geo(:)
      integer             :: ixc,iyc
c     ------------------------------------------------------
      call get_horiz_ncc(geo,ixc,iyc)
      is_land = .false. ! default wet
      if ((ixc<1).or.(iyc<1).or.(ixc>nx).or.(iyc>ny)) return 
      is_land = ((wetmask(ixc,iyc) == 0).or.at_coast_line(geo)) ! valid array access
      end function is_land


c     >>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<
c     >>>>>>>>>> particle_tracking.f:begin    <<<<<<<<<<<<
c     >>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<


      subroutine multiple_reflection_path(s0, s1, anycross, sref, shit1)  ! update to particle_tracking.f 
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



      subroutine get_horiz_ncc(geo,ixc,iyc) !,xc,yc)   ! paste-in dummy
c     ------------------------------------------ 
c     Resolve the node-centered grid cell containing xyz
c     (ixc,iyc) is grid indices of the closest node
c     (xc,yc)   is (lon, lat) coordinates of the closest node
c     no range check
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      integer, intent(out) :: ixc,iyc
c      real, intent(out)    :: xc,yc
c     ------------------------------------------ 
      ixc = 1 + nint((geo(1)-lambda1)/dlambda)     
      iyc = 1 + nint((geo(2)-phi1)   /dphi)   
c      xc  = lambda1 + (ixc-1)*dlambda  ! unused
c      yc  = phi1    + (iyc-1)*dphi     ! unused
      end subroutine get_horiz_ncc



      subroutine get_horiz_ncc_corners(ixc,iyc,c00,c01,c10,c11)   ! paste-in dummy
c     ------------------------------------------ 
c     Resolve corners (c00,c01,c10,c11) in (lon, lat) of the 
c     node-centered horizontal grid cell containing with 
c     cell center (ixc,iyc)
c     ------------------------------------------ 
      integer, intent(in) :: ixc,iyc
      real, intent(out)   :: c00(2),c01(2),c10(2),c11(2)
      real                :: xc,yc
c     ------------------------------------------ 
      xc  = lambda1 + (ixc-1)*dlambda  
      yc  = phi1    + (iyc-1)*dphi     
      c00(1) = xc - 0.5*dlambda
      c01(1) = xc - 0.5*dlambda
      c10(1) = xc + 0.5*dlambda
      c11(1) = xc + 0.5*dlambda      
      c00(2) = yc - 0.5*dphi
      c01(2) = yc + 0.5*dphi
      c10(2) = yc - 0.5*dphi
      c11(2) = yc + 0.5*dphi
c     ------------------------------------------ 
      end subroutine get_horiz_ncc_corners


      
      subroutine cross_2Dline_segments(x0,x1,y0,y1,s,cross)    ! paste-in dummy
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

      end module

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      program test_it
c     ------------------------------------------------------------------------
c     Compile & run:  ifort -e90 coast_line_intersection_hash_tables.f; a.out
c     ------------------------------------------------------------------------
c     Plot results in R       
c     dat <- read.table("fort.10");
c     plot(NA,asp=1,xlim=range(dat[,c(1,3,5,7)]),ylim=range(dat[,c(2,4,6,8)]));
c     segments(dat$V1,dat$V2,dat$V3,dat$V4);
c     segments(dat$V5,dat$V6,dat$V7,dat$V8);
c     abline(v=-180:180+0.5,h=-90:90+0.5,col="lightgrey",lty=3)
c     ----------------------------------------------------------------
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
      stability = 8   ! log10 of number of steps
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

c         write(65,*) geo1(1:2)
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

      
