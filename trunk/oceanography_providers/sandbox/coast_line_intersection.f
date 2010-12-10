c     ------------------------------------------
c     TODO  
c     2) multi reflect update ??
c     ---------------------------------------
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


      contains 


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
         write(*,*) "coast_line_intersection: warning: geo1=",
     +              geo1, "is a dry point"
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


      
      logical function repel_from_coast_line(geo)   ! -> regular_lonlat_grid.f
c----------------------------------------------
c     Internal function that repels geo(1:2) inplace
c     from the coast line. 
c     Includes test of whether geo(1:2) is numerically 
c     on a coast line. Enforce outter domain range: 
c     0.5 < geo(1:2) < (nx,ny) + 0.5
c     For corner points, coordinate is displaced to 
c     a (strictly) random valid adjecent wet cell to avoid bias
c----------------------------------------------
      real, intent(inout) :: geo(:)
      logical             :: east_west,north_south,corner   ! for check_coastal_proximity
      integer             :: iwe, iso, dir(2,4), ndir       ! for check_coastal_proximity
      integer             :: pick
      real                :: rnum
c----------------------------------------------
      repel_from_coast_line = check_coastal_proximity(geo, 
     +                           east_west, north_south, corner,  
     +                           iwe, iso, dir, ndir)

      if (.not.repel_from_coast_line) return
c
c     At this point we know that geo has been flagged as 
c     coastal point; repel it
c      
      if (ndir == 0) stop "repel_from_coast_line: no wet cells"
        
      if (east_west) then   ! face coordinate == iwe + 0.5

         if (ndir /= 1) stop "repel_from_coast_line: unexpected e/w"
         geo(1) = x2lon(iwe + 0.5) + dir(1,1)*facial_tol(1) ! only modify lon

      elseif (north_south) then ! face coordinate == iso + 0.5
      
         if (ndir /= 1) stop "repel_from_coast_line: unexpected n/s"
         geo(2) = y2lat(iso + 0.5) + dir(2,1)*facial_tol(2) ! only modify lat

      elseif (corner) then

         if (ndir > 3) stop "repel_from_coast_line: unexpected corner"
         call random_number(rnum)  
         pick = 1 + max(int(ndir*rnum),ndir-1) ! 1 <= pick <= ndir
         geo(1) = x2lon(iwe + 0.5) + dir(1,pick)*facial_tol(1)
         geo(2) = y2lat(iso + 0.5) + dir(2,pick)*facial_tol(2) 

      else
         stop "repel_from_coast_line: unexpected"
      endif

      end function repel_from_coast_line
      




      logical function at_coast_line(geo)             ! -> regular_lonlat_grid.f
c----------------------------------------------
c     Just return the main result whether geo
c     geo is numerically on a coast line and
c     waste detailed results from check_coastal_proximity
c     Do not modify geo 
c----------------------------------------------
      real, intent(in)    :: geo(:)
      logical             :: east_west,north_south,corner   ! dummies for check_coastal_proximity
      integer             :: iwe, iso, dir(2,4), ndir       ! dummies for check_coastal_proximity
c----------------------------------------------      
      at_coast_line = check_coastal_proximity(geo,    
     +                   east_west, north_south, corner,  
     +                   iwe, iso, dir, ndir)      
      end function at_coast_line





      logical function check_coastal_proximity(geo,     ! -> regular_lonlat_grid.f
     +                       east_west, north_south, corner,  
     +                       iwe, iso, dir, ndir)
c---------------------------------------------------------------------
c     Flag a coastal point (logical work horse in this setup). 
c
c     The test consists of two steps:
c
c       1) the point straddles one (or two) cell faces withing numerical limit 
c          set by at_lon_face/at_lat_face 
c       2) The straddled face (or corner) separates land and water
c
c     Input:  horizontal position geo (lon,lat)
c     Output: return value :        geo is/is not a coastal point
c                                   if geo exceeds numerical coast range of
c                                   outer domain edges, it will be flagged as coastal
c                                   at nearest outer domain face  
c
c             logical east_west:    geo is on an east/west face   OR
c             logical north_south:  geo is on an north/south face OR
c             logical corner        geo is in a  corner
c                                   (at maximum one of ast_west/north_south/corner is flagged)
c
c             integer iwe    left side cell of east/west face (only) resolved for
c                            east_west/corner situation. The face has 
c                            The face has x grid coordinates iwe+0.5
c                            Return range (if set) 0 <= iwe <= nx
c                            iwe = 0/nx corresponds to west/east edge of outer domain
c
c             integer iso    lower side cell of north/south face (only) resolved for
c                            north_south/corner situation
c                            The face has y grid coordinates iso+0.5
c                            Return range (if set) 0 <= iso <= ny
c                            iso = 0/ny corresponds to south/north edge of outer domain
c
c             integer dir(2,4) Non normalized direction vector (geodir, dir_number) indicating direction
c                              indicating direction of a valid wet cell starting to offending point
c             integer ndir     number of resolving directions in dir
c                              As a special case, ndir==0 means that no directions lead to
c                              a wet (valid) cell - this exception is left for external handling
c
c---------------------------------------------------------------------
      real,    intent(in)  :: geo(:)    ! accesses (1:2)
      logical, intent(out) :: east_west
      logical, intent(out) :: north_south
      logical, intent(out) :: corner   
      integer, intent(out) :: iwe    
      integer, intent(out) :: iso  
      integer, intent(out) :: dir(:,:) ! writes dimension (1:2,1:4) 
      integer, intent(out) :: ndir
c
      integer              :: ixc,iyc
      logical              :: wet_n, wet_s, wet_e, wet_w
c     -------------------------------------------------------  
c     Part 1): test point straddles one (or two) cell faces
c     ------------------------------------------------------- 
      check_coastal_proximity = .false.   ! default value
      east_west   = at_lon_face(geo, iwe) ! iwe = western  side
      north_south = at_lat_face(geo, iso) ! iwe = southern side
      corner      = east_west.and.north_south
c     
c     test if quick exit is possible
c
      if (.not.(east_west.or.north_south)) return  ! quick exit - do not resolve further descrioptors
         
c     .... have only one situation east_west/north_south/corner flagged ....
      if (corner) then  
         east_west   = .false.
         north_south = .false.
      endif
      
c     -------------------------------------------------------
c     Part 2): test whether straddled face (or corner) separates land and water
c              and identify directions (dir) pointing to adjecent wet cells
c     -------------------------------------------------------
c     now either east_west/north_south/corner and iwe and iso 
c     Proceed to resolve resolving directions
c     For east_west/north_south one shold just spot the wet cell
c 
      call get_horiz_ncc(geo,ixc,iyc) ! resolve cell coordinate of non-facial direction
c
      dir  = 0    ! default value dir(2,4) 
      ndir = 0    ! default value
c     ------------------------------------
      if (east_west) then
c     ---------------- handle left domain bound separately
         if (iwe==0) then   
            if (wetmask(1,iyc)>0) then
               dir(1,1) = 1     ! we can only move right
               ndir     = 1
            else   
               continue  ! only diagnose exception here (dir = ndir = 0)
            endif
            return       ! no further processing
c        ------------ handle right domain bound separately
         elseif (iwe==nx) then
            if (wetmask(nx,iyc)>0) then
               dir(1,1) = -1    ! we can only move left
               ndir     = 1
            else   
               continue  ! only diagnose exception here (dir = ndir = 0)
            endif
            return       ! no further processing
c        ------------ interior cells 
         else
            wet_w = (wetmask(iwe,iyc)   == 1) ! array access valid
            wet_e = (wetmask(iwe+1,iyc) == 1) ! array access valid
c
            if (wet_w.and.wet_e) return ! this is a wet-wet face, exit
            if (.not.(wet_w.or.wet_e)) then ! both cells dry
               return ! only diagnose exception here (dir = ndir = 0)
            endif
c           ... now we know either wet_w OR wet_e          
            check_coastal_proximity = .true.  ! this is a coastal case
            if     (wet_w) then
               dir(1,1) = -1     ! move right
               ndir     = 1
            elseif (wet_e) then
               dir(1,1) =  1     ! move right
               ndir     = 1
            else
               stop "check_coastal_proximity: unexpected w/e"
            endif
         endif ! if (iwe==
c     ------------------------------------          
      elseif (north_south) then
c        ---------------- handle lower domain bound separately
         if (iso==0) then   
            if (wetmask(ixc,1)>0) then
               dir(2,1) = 1     ! we can only move up
               ndir     = 1
            else   
               continue  ! only diagnose exception here (dir = ndir = 0)
            endif
            return       ! no further processing
c        ------------ handle upper domain bound separately
         elseif (iso==ny) then
            if (wetmask(ixc,ny)>0) then
               dir(2,1) = -1    ! we can only move down
               ndir     = 1
            else   
               continue  ! only diagnose exception here (dir = ndir = 0)
            endif
            return       ! no further processing
c        ------------ interior cells 
         else
            wet_s = (wetmask(ixc,iso)   == 1) ! array access valid
            wet_n = (wetmask(ixc,iso+1) == 1) ! array access valid
c
            if (wet_s.and.wet_n) return ! this is a wet-wet face, exit
            if (.not.(wet_s.or.wet_n)) then ! both cells dry
               return ! only diagnose exception here (dir = ndir = 0)
            endif
c           ... now we know either wet_s OR wet_n          
            check_coastal_proximity = .true.  ! this is a coastal case
            if     (wet_s) then
               dir(2,1) = -1     ! move down
               ndir = 1
            elseif (wet_n) then
               dir(2,1) =  1     ! move up
               ndir = 1
            else
               stop "check_coastal_proximity: unexpected n/s"
            endif
         endif ! if iso==

c     ------------------------------------ 
      elseif (corner) then
c   
c     loop over potential (interior) cells
c         
         do ixc = iwe, iwe+1
            if ((ixc<1).or.(ixc>nx)) cycle      ! avoid exterior cells
            do iyc = iso, iso+1
               if ((iyc<1).or.(iyc>ny)) cycle   ! avoid exterior cells
               if (wetmask(ixc,iyc) == 0) then  ! valid access
                  cycle                         ! dry cell, proceed to next
               else                             
                  ndir = ndir + 1               ! wet cell, accept it
                  dir(1,ndir) = 2*(ixc-iwe) - 1
                  dir(2,ndir) = 2*(iyc-iso) - 1
               endif
            enddo
         enddo
c        ---- ndir == 4 means all wet ----
         if (ndir < 4) check_coastal_proximity = .true. ! conclude corner analysis
c         
      else
         stop "repel_from_coast_line: unexpected (overall case)"
      endif

      end function check_coastal_proximity



      
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




      logical function at_lon_face(geo, iwe)              ! -> regular_lonlat_grid.f
c--------------------------------------------------------------------------------
c     Output: at_lon_face = .true. if geo(1) is at a east/west face (else .false.)
c
c             iwe: left cell of the face; the face has x grid coordinates iwe+0.5
c                  iwe is only meaningful, if at_lon_face = .true.
c
c     Use half the coastal push length (facial_tol) to avoid placing particles at detection edge 
c     If geo exceeds east/west domain faces, flag geo as coastal and range truncate
c     range truncate iwe to 0 <= iwe <= nx 
c--------------------------------------------------------------------------------
      real,    intent(in)  :: geo(:)    ! accesses (1:2)
      integer, intent(out) :: iwe
      real                 :: x,xh,dxh
c------------------------------------
      at_lon_face = .false.    ! default case
      call lon2x(geo, x, iwe) ! iwe is left cell, if geo on a east/west cell face
      xh  = x-0.5
      dxh = abs(xh-nint(xh))
      if (dxh < 0.5*facial_tol(1)) at_lon_face = .true.  ! close to interior coastal face, iwe set
      if ((iwe<0).or.(iwe>nx)) then              ! out-of-domian case
         at_lon_face = .true.
         iwe = max(0, min(iwe,nx))
      endif
      end function at_lon_face
         
     

      logical function at_lat_face(geo, iso)              ! -> regular_lonlat_grid.f
c--------------------------------------------------------------------------------
c     Output: at_lat_face = .true. if geo(2) is at a north/south face (else .false.)
c
c             iso: south cell of the face; the face has y grid coordinates iso+0.5
c                  iso is only meaningful, if at_lat_face = .true.
c
c     Use half the coastal push length (facial_tol) to avoid placing particles at detection edge 
c     If geo exceeds north/south domain faces, flag geo as coastal and range truncate
c     range truncate iso to 0 <= iso <= ny 
c--------------------------------------------------------------------------------
      real,    intent(in)  :: geo(:)    ! accesses (1:2)
      integer, intent(out) :: iso
      real                 :: y,yh,dyh
c------------------------------------
      at_lat_face = .false.    ! default case
      call lat2y(geo, y, iso) ! iso is left cell, if geo on a east/west cell face
      yh  = y-0.5
      dyh = abs(yh-nint(yh))
      if (dyh < 0.5*facial_tol(2)) at_lat_face = .true.  ! close to interior coastal face, iso set
      if ((iso<0).or.(iso>ny)) then              ! out-of-domian case
         at_lat_face = .true.
         iso = max(0, min(iso,ny))
      endif
      end function at_lat_face



      subroutine lon2x(geo, x, iwe)              ! -> regular_lonlat_grid.f
c--------------------------------------------------------------------------------
c     Raw transformation to virtual grid coordinates, without range check
c     Input:  geo(1:2) is a position in (lon,lat)
c     Output: x:  x grid coordinate, cell centers at intergers, faces at half integers
c             iwe: left cell of the face, if geo(1) is on a east-west face; the face has 
c             x grid coordinates iwe+0.5
c--------------------------------------------------------------------------------
      real, intent(in)     :: geo(:)
      real, intent(out)    :: x
      integer, intent(out) :: iwe
      x  = 1.0 + (geo(1)-lambda1)/dlambda   ! virtual coordinate
      iwe = nint(x-0.5)                    ! left cell
      end subroutine lon2x


      subroutine lat2y(geo, y, iso)              ! -> regular_lonlat_grid.f
c--------------------------------------------------------------------------------
c     Raw transformation to virtual grid coordinates, without range check
c     Input:  geo(1:2) is a position in (lon,lat)
c     Output: y:  y grid coordinate, cell centers at intergers, faces at half integers
c             iso: south cell of the face, if geo(2) is on a north-south face; the face has 
c             y grid coordinates iso+0.5
c--------------------------------------------------------------------------------
      real, intent(in)     :: geo(:)
      real, intent(out)    :: y
      integer, intent(out) :: iso
      y  = 1.0 + (geo(2)-phi1)/dphi   ! virtual coordinate
      iso = nint(y-0.5)              ! left cell
      end subroutine lat2y



      real function x2lon(x)  ! -> regular_lonlat_grid.f
c     ------------------------------------------------------
      real, intent(in)  :: x
      x2lon = lambda1 + (x-1.0)*dlambda
      end function x2lon



      real function y2lat(y)  ! -> regular_lonlat_grid.f
c     ------------------------------------------------------
      real, intent(in)  :: y
      y2lat = phi1 + (y - 1.0)*dphi 
      end function y2lat



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


      subroutine multiple_reflection_path(s0, s1, anycross, sref, shit1)  ! paste-in from particle_tracking.f ??
c-------------------------------------------------------------------------
c     Compute key points of multiple horizontal coastal reflection path
c     by iterative application of coast_line_intersection primitive
c     when trying to move from (valid, wet) s0 to new position s1
c     If a coast line is crossed retrun anycross == .true.
c     Then sref is the (multiple) horizontal coastal reflected end point
c     which is valid and wet. shit1 is the point of first coastal intersection
c     If anycross == .false. (sref, shit) is undefined
c-------------------------------------------------------------------------
      real,intent(in)      :: s0(3), s1(3)
      logical,intent(out)  :: anycross
      real,intent(out)     :: sref(3), shit1(3)

      real                 :: start(3), sfin(3), shit(3)
      integer              :: iter
      integer, parameter   :: max_reflections = 5 ! otherwise time step too long ...
c---------------------------------------------------------------
      if (verbose>0) write(*,*) " multiple reflection analysis begin"

      call coast_line_intersection(s0, s1, anycross, sref, shit1)

      if (verbose>0) then
         if (anycross) then 
            write(*,422) s0, s1, sref, shit1  
            write(*,*) "is_land(reflected point) = ", is_land(sref)
         else
            write(*,423) s0, s1
         endif
      endif ! verbose

      if (.not.anycross) return
      
c
c     --- we hit the coast line, start multiple reflection analysis
c         (this means (sref, shit1) are defined and anycross == .true.)
c
      iter  = 1
      shit  = shit1
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
               write(*,422) start, sfin, sref, shit
               write(*,*) "is_land(reflected point) = ", is_land(sref)
            else
               write(*,423) start, sfin
            endif
         endif ! verbose

         iter  = iter + 1
      enddo
      if (iter > max_reflections) then
         write(*,*) "multiple_reflection_path: max_reflections exceeded"
         stop
      endif
 
 422  format(3f7.2, " ->", 3f7.2, ":ref=", 3f7.2, " hit=", 3f7.2)    
 423  format(3f7.2, " ->", 3f7.2, ":no cross")
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
c     -------------------------------------------------
c     Compile & run:  ifort -e90 coast_line_intersection.f; a.out
       
c     Plot results in R       
c     dat <- read.table("fort.10");
c     plot(NA,asp=1,xlim=range(dat[,c(1,3,5,7)]),ylim=range(dat[,c(2,4,6,8)]));
c     segments(dat$V1,dat$V2,dat$V3,dat$V4);
c     segments(dat$V5,dat$V6,dat$V7,dat$V8);
c     abline(v=-180:180+0.5,h=-90:90+0.5,col="lightgrey",lty=3)
c     -------------------------------------------------
      use test_coast_line_intersection
      implicit none
      real    :: geo1(3), geo2(3), georef(3), geohit(3),deg2rad
      real    :: dgeo(3),x,y,s
      logical :: anycross
      real    :: angle,r2(2)
      integer :: ixc,iyc,ix,iy,istep
      real,parameter :: rjump = 0.5  ! Bee jump range
      real    :: stability
c     -------------------------------------------------

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
      stability = 9   ! log10 of number of steps
      call get_horiz_ncc(geo1,ixc,iyc)     
      !Surround the point with dry cells
      wetmask   = 0 ! 0 dry, 1 wet
      wetmask(ixc,iyc)     = 1  ! 
      wetmask(ixc+1,iyc)   = 1  ! 
      wetmask(ixc+1,iyc+1) = 1  ! 
      do istep = 1, nint(10**stability)
         call random_number(r2)
         angle = 8*atan(1.0)*r2(2)
         dgeo  = r2(1)*rjump*(/cos(angle), sin(angle), 0.0/)
         geo2  = geo1+dgeo
         call multiple_reflection_path(geo1,geo2,anycross,georef,geohit)
         if (anycross) geo2 = georef

         if (is_land(geo2)) then
            if (anycross) then
               write(*,522) istep, geo1(1:2),geo1(1:2)+dgeo(1:2),
     +                    georef(1:2),geohit(1:2)
            else
               write(*,523) istep, geo1(1:2), geo2(1:2)
            endif
 522        format(i8,"    cross g1=", 2f12.6," g1+dg=",
     +             2f12.6," ref=",2f12.6," hit=",2f12.6)
 523        format(i8," no cross g1=", 2f12.6," g1+dg=",2f12.6)
         else
            geo1 = geo2 ! accept step
         endif

c         write(*,*) geo1(1:2)
c         if (is_land(geo1)) then
c            write(*,*) "geo1 dry at step",istep
c            write(*,*) "geo1 = ", geo1
c            write(*,*) geo1(1:2)
c            stop
c         endif
      enddo



      end program

      
