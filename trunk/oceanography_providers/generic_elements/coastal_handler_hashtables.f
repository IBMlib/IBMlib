ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     -----------------------------------------------------------
c     Hashtable-base coastal handler
c     -----------------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
c     Module based implementation of coast_line_intersection based on 
c     hash tables for a static coastal geometry.
c     Appears to be somewhat faster than other version which
c     features runtime analysis of local coastal topology
c     Based on original implemention by ASC, 14Dec2010
c     Extracted into module by MPA, 4Feb2011
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module coastal_handler
      use horizontal_grid_transformations
      use geometry
      implicit none
c      private     
    
      integer,public     :: wetmask(nx,ny) = 0 ! (1: wet nodes, 0 dry)  
  
      real, parameter    :: facial_tol(2) = 1.0e-5 ! note: in grid units lon part should be lat dependent
      integer,parameter,private  :: verbose = 0

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

      end subroutine coast_line_intersection ! embedding subroutine


      subroutine assess_cell_crossing (xy1,xy2,ix,iy,s,exitface) 
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
      end subroutine assess_cell_crossing     


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

      end module
