      module horizontal_representation
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Generic module for static horizontal issues
c     (i.e. not wdepth) including the coast line. 
c     In principle, it should be possible to change the coast line representation
c     without having to change modules above this. Therefore exported 
c     array wetmask should only be used to signal, whether data points on grid
c     are wet/dry      
c     
c     $Rev: $
c     $LastChangedDate:  $
c     $LastChangedBy: $ 
c
c     Coast line topography: 
c         the node centered cells of the lon-lat grid may be wet/dry
c         i.e. the coast line is constituted by NS/EW conected line segments
c         with resolution given by the underying lon-lat grid 
c         subrouitne coast_line_intersection/is_land/is_wet strictly 
c         enforces this coast line 
c
c     Bottom topography:
c         water depth is mapped on the lon-lat grid points. In between
c         lon-lat grid points bottom topography follows bilinear interpolation
c         and is thus continuous everywhere, except at the coast line,
c         where is jumps to zero.
c
c     This version of coast_line_intersection "staywet" abandon the 
c     coast line cushin approach, but rather focuses 
c     on the consistency in decisions, essentially.
c     The "traj" tag indicates that this is the trajectory resolving
c     variant for intermediate cells between initial and final cells.
c
c       1) is_land is the authoritative 
c       2) if coast_line_intersection finds a crossing, then
c          hit point must test wet by island (but reflection point may be wet/dry)
c       3) if multiple_reflection_path finds a crossing, then
c          both first hit point and reflection point must be wet
c     
c     Generalized from regular_lonlat_grid.f:  ASC  03 Feb, 2011
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use horizontal_grid_transformations   ! injects nx,ny
      use geometry                          ! cross_2Dlines
     
      implicit none
      private           ! default scope

      public :: init_horizontal_representation   ! init  this module
      public :: close_horizontal_representation  ! close this module
   
      public :: horizontal_range_check  ! reexport from horizontal_grid_transformations
      public :: is_land
      public :: coast_line_intersection      

c     -------------------- module data --------------------  
      
      integer, parameter :: verbose = 0  ! debugging output control
c
c     ------ grid dimensions:   ------
c     
      public      :: nx,ny ! reexport from horizontal_grid_transformations 
      
c      
c     ----- wetmask tells wheter a grid point is wet/dry -----
c     
      integer, allocatable,public :: wetmask(:,:) ! auxillary to define coastal geometry (1: wet nodes, 0 dry)  


c     ===================================================
                            contains
c     ===================================================

                      
      subroutine init_horizontal_representation()
c     ------------------------------------------------------
c     Assumes (nx,ny) has been set by client module
c     NB: no initialization of sub module horizontal_grid_transformations
c     because argument sequence to init_horiz_grid_transf 
c     is specific to the particular horizontal  
c     No parts of this initialization depends on whether 
c     module horizontal_grid_transformations has been initialized (check)
c     ------------------------------------------------------
      write(*,*) "init_horizontal_representation: "//
     +           "coast_line_intersection algorithm = staywet"
      write(*,*) "init_horizontal_representation: allocating wetmask"

      allocate( wetmask(nx,ny)    )  

      end subroutine init_horizontal_representation
       



      subroutine close_horizontal_representation()
c     ------------------------------------------------------
c     ------------------------------------------------------    
      if (allocated(wetmask)) deallocate( wetmask )
      end subroutine close_horizontal_representation


      

      subroutine coast_line_intersection(geo1, geo2,
     +                                   cross_coast,
     +                                   georef, geohit) ! COPY -> horizontal_representation_stay_wet.f
c     ------------------------------------------ 
c     Trajectory tracking variant
c
c     This elementary function is a service function assisting in enforcing 
c     coastal boundary conditions on particle steps. geo1 is the current 
c     valid (wet, inside domain) position. geo2 is the next position that may or may not be in water.
c     The assumption is that the particle moves in a straight line 
c     between geo1 and geo2 in grid space. This algorithm 
c     solves the problems exactly (within numerical precision) 
c     for the any coastal geometry represented as regular cell centered squares 
c     in grid space (without previous assumptions about
c     small steps geo1 -> geo2, so that any number of neighbor cells can be crossed)

c     Output: 
c     * anycross: If the straight line geo1 -> geo2 has crossed a coast line, 
c       anycross is returned .true. otherwise .false. (i.e. the line between geo1 and 
c        geo2 is all in water). 
c
c     If anycross is .true., the following vectors are computed:
c
c     * georef is the modified final position, if step from geo1 to geo2 is reflected
c       in the coast line. georef may be both wet/dry
c     * geohit is the first position where the line from geo1 to geo2 
c       crosses a coast line. geohit is guarentied a wet point
c
c     If anycross is .false. neither georef and geohit are assigned.
c
c     If geo1 is dry, a warning is issued and anycross is returned false
c     (and neither georef/geohit computed)
c     This implementation does not handle multiple reflection effects, so the
c     reflection may cross the coast line additionally (the step steps geo1 -> geo2
c     may become reflected onto land in concave coastal geometries for large steps)
c     geo1 and geo2 must have same vector length (2 or 3) so they are vector addable
c     Currently it is being asserted that geo2 does not violate the outer domain (not checked)
c     It is also assumed that wet points exist arbitrarely close to geo1 in the
c     direction of geo2.
c
c     Numerical robustness: if geo2 is dry and geo1 wet, this algo will always signal 
c     cross_coast == .true.
c
c     Conceptual revision ASC/MPA Nov 24, 2010 
c     Infinitesimal coast line repulsion added ASC Dec 06, 2010 
c     Resolve association conflict at cell transition analysis, ASC Jan 18, 2011 
c     Internal grid space implementation, ASC Feb 02, 2011
c     Backtracking algortihm to assure that geohit is always wet, if computed, ASC Feb 11, 2011
c     Improved handling of extremely rare cases, ASC Feb 11, 2011   
c     ---------------------------------------------------- 
      real, intent(in)     :: geo1(:),geo2(:)
      logical, intent(out) :: cross_coast
      real, intent(out)    :: georef(:), geohit(:)
c   
      real                 :: xy1(size(geo1)),xy2(size(geo2))         ! grid coor set
      real                 :: xyref(size(georef)),xyhit(size(geohit)) ! grid coor set
      real                 :: assured_wet(size(geo2))
      real                 :: s
      logical              :: leaving,loopflag,ldum
      character*1          :: reftype, xface
      real                 :: direct(size(geo1)), reflect(size(geo1))
      integer              :: ix,iy,ix1,iy1,ix2,iy2,istep,maxsteps
c     ---------------------------------------------------- 
      call get_horiz_ncc_index(geo1,ix1,iy1)  ! resolve start point cell
      call get_horiz_ncc_index(geo2,ix2,iy2)  ! resolve end point cell
c
c     Test assertion that geo1 is wet. 
c     Do not resolve geohit/georef if geo1 is dry
c
      if (wetmask(ix1,iy1)==0) then
         write(*,347) geo1(1:2)        
 347     format("coast_line_intersection: warning: geo1=",
     +          2f10.6, " is a dry point - abandon analysis")
         cross_coast = .false.
         if (verbose>0) call write_summary() ! only if verbose
         return ! assertion failed: geo1 dry
      endif
c
c     Test if geo1 and geo2 are in same cell => no coast line crossing
c     most cases terminate here
c
      if ((ix1==ix2).and.(iy1==iy2)) then
         if (verbose>0) write(*,*) "geo1 -> geo2: step in same wet cell"
c        --- handle coastal osculation, if detected
         cross_coast = point_on_land(geo2,geo1,assured_wet)
         if (cross_coast) then
            geohit = assured_wet
            georef = assured_wet
         endif
         if (verbose>0) then
            if (cross_coast) then
               write(*,*) "geo2 coastal osculation handled"
            else
               write(*,*) "geo2 is clear of coast"
            endif
            call write_summary() ! only if verbose
         endif
         return ! geo1 and geo2 in same box - most cases terminate here
      endif

c     Now we know that we are leaving the cell. But where? 
c     Follow sequence moving from one cell to the next until we either
c     reach the end of the vector (s=1), or encounter land
c     Interior part here is in grid coordinates

      if (verbose>0) then
         write(*,*) "geo1 -> geo2 is an intercell step"
         write(*,*) "trajectory analysis:begin"
      endif

      call get_horiz_grid_coordinates(geo1,xy1)
      call get_horiz_grid_coordinates(geo2,xy2)
      cross_coast = .false. ! default unless detected below
      loopflag    = .true.
      istep       = 1
      maxsteps    = max(1,abs(ix1-ix2))*max(1,abs(iy1-iy2)) ! safety plug

      ix = ix1   ! start in geo1 cell
      iy = iy1   ! start in geo1 cell
      do while(loopflag)
         call assess_cell_crossing(xy1,xy2,ix,iy,s,xface)
         if (verbose>0) write(*,421) istep,ix,iy,xface
 421     format("substep",i3,": leaving cell ",2i4," via ",a," face") 
         !Move into the next cell
         select case(xface)
         case ("N")
            iy = iy + 1   
         case ("S")
            iy = iy - 1
         case ("E")
            ix = ix + 1
         case ("W")
            ix = ix - 1
c        
c        We end in the default case step when xy1 -> xy2 is essentially a zero step
c        in which case tracking equations in assess_cell_crossing are 
c        ill conditioned. As final strategy, pick any potential intercell step
c
         case default  
            if     (ix2 > ix) then ! E step potential, do that
               xface = "E"
               ix    = ix + 1   
            elseif (ix2 < ix) then ! W step potential, do that
               xface = "W"
               ix    = ix - 1
            elseif (iy2 > iy) then ! N step potential, do that
               xface = "N"
               iy    = iy + 1
            elseif (iy2 < iy) then ! S step potential, do that
               xface = "S"
               iy    = iy - 1
            endif
            s = 0.5 ! render s defined
         end select
c
c        Is the next cell wet or dry? If its dry, then we've hit the coast.  
c        Move into the next cell
c
         if (verbose>0) then
            if(wetmask(ix,iy)>0) then
               write(*,422) ix,iy,"wet"
            else
               write(*,422) ix,iy,"dry" 
            endif
422        format("new cell : ",2i4," is ", a) 
         endif

         ! --- check if we should abort trajectory analysis (set loopflag = .false.):
         !     either coast is crossed / geo1 -> geo2 all in water / maxsteps exceeded
         if(wetmask(ix,iy) == 0) then          ! 1) have crossed onto land
            cross_coast = .true.               !    so setup for exit
            loopflag    = .false.  
            if (verbose>0) then
               write(*,*) "coast line crossed, resolve geohit/georef" 
            endif
         elseif ((ix==ix2).and.(iy==iy2)) then ! 2) we are in final cell without crossing the coast  
            loopflag    = .false.              ! so setup for exit, still cross_coast == .false.
            if (verbose>0) then
               write(*,*) "new cell is final, geo1 -> geo2 all in water" 
            endif
         elseif (istep>maxsteps) then          ! 3) loop exhausted, we should not end here, no plan B
            write(*,*) "coast_line_intersection: maxsteps exceeded)"
            stop
         endif

         istep = istep + 1

      enddo
c  
c     if (cross_coast == .false.) only set geohit and georef if
c     geo2 osculates the coast, otherwise let geohit and georef 
c     remain undefined
c
      if (.not.cross_coast) then
c        --- handle coastal osculation, if detected
         cross_coast = point_on_land(geo2,geo1,assured_wet)
         if (cross_coast) then
            geohit = assured_wet
            georef = assured_wet
         endif
         if (verbose>0) then
            if (cross_coast) then
               write(*,*) "geo2 coastal osculation handled"
            else
               write(*,*) "geo2 is clear of coast"
            endif
            call write_summary() ! only if verbose
         endif
         return  ! no coast line crossed (or coastal osculation)
      endif

c
c     The coast line was crossed geo1->geo2, resolve xyref and xyhit
c     0 < s < 1 is the coordinate along the vector (xy2-xy1) 
c     corresponding to the point xyhit
c     
      direct  = xy2-xy1
      reflect = direct
      if     (xface=="N" .or. xface =="S") then
         if (verbose>0) write(*,*) "compute N/S reflection"
         reflect(2) = -reflect(2) ! horizontal box line reflection
      elseif (xface=="E" .or. xface =="W") then
         if (verbose>0) write(*,*) "compute E/W reflection"
         reflect(1) = -reflect(1) ! vertical   box line reflection
      else
         stop "coast_line_intersection: unhandled reflection type"
      endif
c     |direct| == |reflect|        
      xyhit = xy1   +     s*direct  ! position where coast line is hit
      xyref = xyhit + (1-s)*reflect ! position for reflection in coast line
      call get_horiz_geo_coordinates(xyref,georef) ! return to geo coordinates
      call get_horiz_geo_coordinates(xyhit,geohit) ! return to geo coordinates
c
c     Finally, when coast line was crossed, make sure geohit tests wet
c     do not test/modify georef
c
      ldum = point_on_land(geohit,geo1,assured_wet) ! waste logical result
      if (ldum) geohit = assured_wet

      if (verbose>0) then
         if (ldum) then
            write(*,*) "geohit coastal osculation handled"
         else
            write(*,*) "geohit is clear of coast"
         endif
         call write_summary() ! only if verbose
      endif  
       
      return   ! ! coast line crossed, georef/geohit computed

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
      if (dxy(1) > 0) then    !we're moving to the east

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
      if (dxy(2) > 0) then !we're moving to the North
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



      logical function point_on_land(maybewet, wetpt, assured_wetpt)
c     --------------------------------------------------------------
c     Ensure that maybewet is wet (i.e. is_land(maybewet) == .false.) 
c     by nudging it toward wetpt in very small steps to weed
c     out uncertainty in numerical solution of coast line intersection
c     equations. Perform algortihm in geo coordinates. Do not affect vertical
c     component of maybewet (if provided)
c     
c     Input:  maybewet      : may or may not be dry
c             wetpt         : assumed wet point
c
c     Output: assured_wetpt : a guarantied wet point if maybewet is dry
c                             undefined if maybewet tests wet          
c             point_on_land:  return .true. when maybewet was dry 
c                             and nudged into assured_wetpt    
c                             return .false. when maybewet was wet 
c                             and return assured_wetpt == maybewet 
c         
c     Assumptions: there exists wet points arbitrary close to wetpt
c                  is_land is the authoritative function for horizontal 
c                  wet/dry condition
c
c     ASC/09Feb2011
c     --------------------------------------------------------------
      real, intent(in)    :: maybewet(:), wetpt(:)
      real, intent(out)   :: assured_wetpt(:)      ! assumed same length as maybewet
      real                :: dg(2),scale,dgmax
      integer             :: ticks,istep
      logical             :: keep_going
      real,parameter      :: resol = 1.0e-6   ! approx last sign digit
      real,parameter      :: scale_default = 0.99    ! scale_min < scale_default < scale_max 
      real,parameter      :: scale_max     = 0.999 ! 0 < scale_min < scale_max < 1
      real,parameter      :: scale_min     = 0.5     ! 0 < scale_min < scale_max < 1
      integer,parameter   :: max_steps     = 1000
c     ------------------------------------------------------------      
      point_on_land = is_land(maybewet) ! .true. : needs nudging

      if (.not.point_on_land) return    ! no further processing
c     
c     resolve scaling factor (0 < scale_min < scale < scale_max < 1)
c     
      dg    = maybewet(1:2) - wetpt(1:2)
      dgmax = maxval(abs(dg))
      ticks = nint(dgmax/resol)
      if (ticks>0) then
         scale = abs(1.0*(ticks-1.0)/ticks)
         scale = min(max(scale, scale_min), scale_max)
      else
         scale = scale_default
      endif
c     
c     nudging loop: set assured_wetpt
c            
      istep         = 1
      assured_wetpt = maybewet ! load vertical component, if present
      keep_going    = .true. ! we know we need at least one step
      do while (keep_going)
         dg = dg*scale  ! scale down difference vector iteratively
         assured_wetpt(1:2) = wetpt(1:2) + dg ! 
         keep_going = is_land(assured_wetpt)
         istep = istep + 1
         if (istep > max_steps) then
            write(*,*) "point_on_land: max_steps exceeded"
            write(*,*) "maybewet=", maybewet
            write(*,*) "wetpt   =", wetpt
            write(*,*) "scale   =", scale
            stop ! fatal, no plan B available
         endif
      enddo
      
      end function point_on_land


      
      subroutine write_summary()
c     -------------------------------------------
c     internal subroutine for debugging
c     use scope of coast_line_intersection 
c     for in/out variables
c     -------------------------------------------
      write(*,*) "summary of coast_line_intersection"
      if (cross_coast) then 
         write(*,622) geo1(1:2), geo2(1:2), georef(1:2), geohit(1:2)
         write(*,*) "is_land(geo1)   = ", is_land(geo1)
         write(*,*) "is_land(geo2)   = ", is_land(geo2)
         write(*,*) "is_land(georef) = ", is_land(georef)
         write(*,*) "is_land(geohit) = ", is_land(geohit)
      else
         write(*,623) geo1(1:2), geo2(1:2)
         write(*,*) "is_land(geo1)   = ", is_land(geo1)
         write(*,*) "is_land(geo2)   = ", is_land(geo2)
         write(*,*) "is_land(georef) = undefined"
         write(*,*) "is_land(geohit) = undefined"
      endif
 622  format("summary of coast_line_intersection: ",
     +        2f7.2, " ->", 2f7.2, " : crossed coast :" 
     +       "ref=", 2f7.2, " hit=", 2f7.2)    
 623  format("summary of coast_line_intersection: ",
     +        2f7.2, " ->", 2f7.2, " : no coast crossed")              
      end subroutine write_summary

      end subroutine coast_line_intersection ! embedding subroutine




      LOGICAL function is_land(geo)        ! COPY -> horizontal_representation_stay_wet.f
c     ------------------------------------------ 
c     return is_land = .false. at horizontal range violation
c     ------------------------------------------ 
      real, intent(in) :: geo(:)
      integer          :: ixc,iyc
c     ------------------------------------------ 
      if (horizontal_range_check(geo)) then ! range OK
         call get_horiz_ncc_index(geo,ixc,iyc)
         is_land = (wetmask(ixc,iyc)==0)
      else
         is_land = .false.
      endif    
c     ------------------------------------------
      end function




      end module
