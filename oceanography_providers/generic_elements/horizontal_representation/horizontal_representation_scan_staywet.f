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
c     The "scan" tag indicates that this is the brute scan variant
c     over intermediate cells between initial and final cells (rather
c     than resolving an intercell trajectory, which may fail in very
c     rare cases)
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
c     Brute scan variant
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
c
c     Output: 
c     * anycross: If the straight line (linear in grid coordinates) from geo1 -> geo2 
c       has crossed a coast line, anycross is returned .true. otherwise .false. 
c       (i.e. the line between geo1 and geo2 is all in water). 
c
c     If anycross is .true., the following vectors are computed:
c
c     * georef is the (single) reflected final position, if step from geo1 to geo2 is reflected
c       in the coast line. georef may be both wet/dry
c     * geohit is the first position where the line from geo1 to geo2 
c       crosses a coast line. geohit is guarentied a wet point
c     Both georef and geohit will be assumed to have length as geo1/geo2 if assigned
c
c     If anycross is .false. neither georef nor geohit are assigned.
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
c     Numerical robustness: 
c        * If geo2 is dry and geo1 wet, this algo will always signal cross_coast == .true.
c        * If geo1 -> geo2 crosses a coastline, but step is smaller than step_resol_limit,
c          then geo1 = georef = geohit is returned and cross_coast = .true.
c
c     Conceptual revision ASC/MPA Nov 24, 2010 
c     Infinitesimal coast line repulsion added ASC Dec 06, 2010 
c     Resolve association conflict at cell transition analysis, ASC Jan 18, 2011 
c     Internal grid space implementation, ASC Feb 02, 2011
c     Backtracking algortihm to assure that geohit is always wet, if computed, ASC Feb 11, 2011
c     Improved handling of extremely rare cases, ASC Feb 11, 2011  
c     Go back to brute scan over potential intermediate cells (no trajectory analysis), ASC Feb 25, 2011   
c     Add lower step size limit for resolving coastal crossings
c     ---------------------------------------------------- 
      real, intent(in)     :: geo1(:),geo2(:)
      logical, intent(out) :: cross_coast
      real, intent(out)    :: georef(:), geohit(:)
c     --- inherit vector length from geo1 to accept both length=2/3
      real                 :: direct(size(geo1)) 
      real                 :: reflect(size(geo1))       
      real                 :: assured_wet(size(geo1))                 
      real                 :: xy1(size(geo1))
      real                 :: xy2(size(geo1))         
      real                 :: xyref(size(geo1))
      real                 :: xyhit(size(geo1)) 
      
      real                 :: s, stest
      logical              :: ldum,all_wet,did_cross,final_cell_is_wet
      character*1          :: reftype, face,facetest
   
      integer              :: ix,iy,ix1,iy1,ix2,iy2
      integer              :: ixmin,ixmax,iymin,iymax

      real, parameter      :: step_resol_limit = 5.e-5 ! limit for solving crossing equations in grid coordinates 
c     ---------------------------------------------------- 
      call get_horiz_ncc_index(geo1,ix1,iy1)  ! resolve start point cell
      call get_horiz_ncc_index(geo2,ix2,iy2)  ! resolve end point cell
      cross_coast = .false.                   ! default, unless set .true. below
c
c     Test assertion that geo1 is wet. 
c     Do not resolve geohit/georef if geo1 is dry
c
      if (wetmask(ix1,iy1)==0) then
         write(*,347) geo1(1:2)        
 347     format("coast_line_intersection: warning: geo1=",
     +          2f10.6, " is a dry point - abandon analysis")
         if (verbose>0) call write_summary() ! only if verbose
         return ! assertion failed: geo1 dry
      endif
c
c     Test if:
c       1) geo1 and geo2 are in same cell => no coast line crossing
c       2) all potential interdemiate cell are wet => no coast line crossing
c     most cases terminate here
c
      
      if ((ix1==ix2).and.(iy1==iy2)) then 
         all_wet = .true.         ! implied since (ix1,iy1) is wet
      else                        ! split logical expression to 
         ixmin   = min(ix1,ix2)   ! avoid the min/max operation in most cases
         ixmax   = max(ix1,ix2)
         iymin   = min(iy1,iy2)
         iymax   = max(iy1,iy2)
         all_wet = all(wetmask(ixmin:ixmax,iymin:iymax) > 0) ! indices in correct order
      endif
      
      
      if (all_wet) then
         if (verbose>0) then
                write(*,*) "coast_line_intersection: "//
     +          "geo1 -> geo2: step in same cell/all wet cells"
                call write_summary() ! only if verbose
         endif
         return ! geo1 and geo2 in same box/wet region - most cases terminate here
      endif

c     Now we know that either: 
c       1) geo1 -> geo2 is at least one intercell step and possibly
c       2) the final cell (of geo2) is dry. 
c     The continued analysis beyond this point is in grid coordinates

      if (verbose>0) then
         write(*,'(a)') "coast_line_intersection: geo1 -> geo2 is an "//
     +              "intercell step scan intermediate cells: begin"
      endif
      
      call get_horiz_grid_coordinates(geo1,xy1)
      call get_horiz_grid_coordinates(geo2,xy2)
c
c     Infinitesimal step exit point: capture very small intercell steps into a dry box
c     at the numerical resolution limit, where intercell transition
c     points can not be resolved (consider this as a coastal osculation case)     
c     Infinitesimal steps are into side/corner neighbor cells of starting cell

      if ((wetmask(ix2,iy2) == 0).and.  ! step into dry box
     +    (sum(abs(xy2(1:2)-xy1(1:2))) < step_resol_limit)) then
         cross_coast = .true.
         geohit      = geo1   ! geo1 is tested wet
         georef      = geo1   ! geo1 is tested wet
         if (verbose>0) then
            write(*,*) "geo1 -> geo2 is below resolution limit"//
     +              "for intercell transition dxy=",xy2-xy1
            call write_summary() ! only if verbose
         endif
         return ! geo1 -> geo2 cross coast, but step below resolution limit
      endif    
c     
c     ---------- intermediate cell loop ----------
c
      face        = "?"     ! set as undefined
      s           = 1.0e12  ! set to very large
c
      do ix = ix1,ix2       ! start in geo1 cell
         do iy = iy1,iy2    ! start in geo1 cell

            if ((ix==ix1).and.(iy==iy1)) then
               if (verbose>0) write(*,422) ix,iy," : omit initial cell"
               cycle ! omit initial cell
            endif
            if ((ix==ix2).and.(iy==iy2)) then
               if (verbose>0) write(*,422) ix,iy," : omit final cell"
               cycle ! omit final cell
            endif
           
            if (wetmask(ix,iy)==0) then ! only analyze dry cells
               if (verbose>0) write(*,422) ix,iy," : is dry (analyze)" 
               call test_cell_crossing(xy1, xy2, ix, iy,
     +                                 stest, facetest, did_cross)
               if (did_cross) then   
                  cross_coast = .true. ! we found at least one solution
                  if (stest<s) then    ! store this solution, if it is closer
                     s    = stest
                     face = facetest
                  endif
               endif 
            else
               if (verbose>0) write(*,422) ix,iy," : is wet (continue)" 
            endif

         enddo ! iy = iy1,iy2 
      enddo    ! ix = ix1,ix2 
 422  format("considering cell ",2i4," is ", a) 

c     
c     ---------- final cell analysis ----------
c        
c     Only perform final cell analysis, if final cell is dry and 
c     no intermediate crossings where detyected. An intermediate 
c     crossing is logically ranked before a final cell.
c     If final cell analysis is invoked, insist on a solution
c     because final cell is dry (i.e. handle potential numerical problems robustly) 
c     If final cell is wet and no intermediate crossing
c     conclude no coastal crossings 1->2 

      if ((.not.cross_coast).and.(wetmask(ix2,iy2)==0)) then
         cross_coast = .true. ! insist on a solution  
         call locate_cell_entrance(xy1,ix1,iy1, xy2,ix2,iy2, s,face)
      endif  ! otherwise cross_coast is still .false.

c     
c     ---------- conclude inter cell analysis ----------
c 
c     if still (cross_coast == .false.) only set geohit and georef if
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
            write(*,*) "concluding no coastal crossing geo1->geo2"
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
c     The coast line was crossed geo1->geo2, now resolve xyref and xyhit
c     0 < s < 1 is the coordinate along the vector (xy2-xy1) 
c     corresponding to the point xyhit
c     
      direct  = xy2-xy1
      reflect = direct
      if     (face=="N" .or. face =="S") then
         if (verbose>0) write(*,*) "compute N/S reflection"
         reflect(2) = -reflect(2) ! horizontal box line reflection
      elseif (face=="E" .or. face =="W") then
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
      ldum = point_on_land(geohit,geo1,assured_wet) 
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


      subroutine test_cell_crossing(xy1, xy2, ix, iy,
     +                              s, face, did_cross) ! internal to coast_line_intersection
c     --------------------------------------------------------------
c     Determine if the direct line from xy1 through xy2 
c     enters the cell with center (ix,iy): did_cross
c     If (did_cross == .true.) then compute 0 < s  which is the coordinate 
c     along the vector (xy2-xy1) where the vector enters the cell. 
c     Also identify the entrence face as one of "N","S","E","W" 
c     so the proper reflection can be calculated
c     If (did_cross == .false.) (i.e. xy1 -> xy2 is not crossing this cell)
c     then (s, face) is set to (1.0e12, "%")
c
c     Modified from assess_cell_crossing to avoid decision conflicts ASC/25Feb2011
c     Also uses cross_2Dline_segments rather than cross_2Dlines
c     --------------------------------------------------------------
      real, intent(in)         :: xy1(:),xy2(:)
      integer, intent(in)      :: ix,iy
      real, intent(out)        :: s
      character*1, intent(out) :: face
      logical, intent(out)     :: did_cross
c
      logical                  :: yes
      real                     :: stest,tdum,dxy(2)
      real                     :: c00(2),c01(2),c10(2),c11(2)
c     --------------------------------------------------------------
      call get_horiz_ncc_corners(ix,iy,c00,c01,c10,c11) ! resolve faces          
c
c     Determine is the path xy1 -> xy2 enters a box defined by (c00,c01,c10,c11))
c     We can make efficiency gains by taking advantage of the facts that
c     the direction xy1 -> xy2 shadows two (or three) faces on a box, so
c     we only need to consider two (or one) face of the box
c
c       West  face: c00,c01
c       East  face: c10,c11
c       North face: c01,c11
c       South face: c00,c10
c
      did_cross = .false.
      face      = "%"    ! default is no crossings
      s         = 1.0e12 ! must be set before first comparison
      dxy       = xy2(1:2)-xy1(1:2) 
c
c     ---------- First the east-west direction ----------
c
      if    (dxy(1) > 0) then  ! moving to the east: consider west face impact

         call cross_2Dline_segments(xy1,xy2,c00,c01,stest,yes)
         if (yes.and.(stest<s)) then ! rely on right-to-left evaluation
            s         = stest
            face      = "W" 
            did_cross = .true.
         endif

      elseif (dxy(1) < 0) then ! moving to the west: consider east face impact

         call cross_2Dline_segments(xy1,xy2,c10,c11,stest,yes)
         if (yes.and.(stest<s)) then ! rely on right-to-left evaluation
            s         = stest
            face      = "E" 
            did_cross = .true.
         endif

      endif  ! no action for dxy(1) = 0

c
c     ---------- Then the north-south direction ----------
c
      if     (dxy(2) > 0) then ! moving to the north: consider south face impact

         call cross_2Dline_segments(xy1,xy2,c00,c10,stest,yes)
         if (yes.and.(stest<s)) then ! rely on right-to-left evaluation
            s         = stest
            face      = "S"  
            did_cross = .true.
         endif
         
      elseif (dxy(2) < 0) then ! moving to the south: consider north face impact

         call cross_2Dline_segments(xy1,xy2,c01,c11,stest,yes)
         if (yes.and.(stest<s)) then ! rely on right-to-left evaluation
            s         = stest
            face      = "N" 
            did_cross = .true.
         endif
         
      endif  ! no action for dxy(2) = 0

      end subroutine test_cell_crossing     ! local subroutine



      subroutine locate_cell_entrance(xy1,ix1,iy1, 
     +                              xy2,ix2,iy2, 
     +                              s,face)    ! internal to coast_line_intersection
c     --------------------------------------------------------------
c     Determine where direct line from xy1 through xy2 
c     enters the cell with center (ix2,iy2). xy1 is situated in cell (ix1,iy1)
c     This subroutine should always generate one solution (s,face)   
c     s > 0  which is the coordinate along the vector (xy2-xy1) where the vector 
c     enters the cell at one of the faces: "N","S","E","W"  
c
c     This subroutine picks which entrance faces to consider based on comparing
c     (ix1,iy1) and (ix2,iy2) to avoid problems due to finite numerical precision
c     This subroutine is hardened against finite numerical precision in 
c     resolving line crossings
c
c     Assume (ix1,iy1) assignment consistent with xy1 (not checked for speed)
c     Assume (ix2,iy2) assignment consistent with xy2 (not checked for speed)
c     --------------------------------------------------------------
      real, intent(in)         :: xy1(:),  xy2(:)
      integer, intent(in)      :: ix1,iy1, ix2,iy2
      real, intent(out)        :: s
      character*1, intent(out) :: face
      
      logical                  :: yes
      real                     :: ssol(2),ss,tt,pen(2),dxy(2)
      character*1              :: cside(2)
      real                     :: c00(2),c01(2),c10(2),c11(2)
      integer                  :: dix,diy,nsol
c     --------------------------------------------------------------
c     
c     -------- 1) Build stack of solutions --------
c
c          entrance at         consider when:
c       West  face: c00,c01      dix > 0 
c       East  face: c10,c11      dix < 0 
c       North face: c01,c11      diy < 0
c       South face: c00,c10      diy > 0
c
       
      dix  = ix2 - ix1   ! steps in x direction
      diy  = iy2 - iy1   ! steps in y direction 
      nsol = 0           ! number of potential solutions
      call get_horiz_ncc_corners(ix2,iy2,c00,c01,c10,c11) ! resolve faces if (ix2,iy2) 
c
      if (dix > 0) then         ! consider west face entrance
         call cross_2Dlines(xy1,xy2,c00,c01,ss,tt,yes)
         if (yes) then ! stack and rank this solution
            nsol        = nsol + 1
            ssol(nsol)  = ss
            pen(nsol)   = penalty(tt) ! ranking
            cside(nsol) = "W"
         endif   
      endif
      if (dix < 0) then         ! consider east face entrance
         call cross_2Dlines(xy1,xy2,c10,c11,ss,tt,yes)
         if (yes) then ! stack and rank this solution
            nsol        = nsol + 1
            ssol(nsol)  = ss
            pen(nsol)   = penalty(tt) ! ranking
            cside(nsol) = "E"
         endif   
      endif
      if (diy < 0) then         ! consider north face entrance
         call cross_2Dlines(xy1,xy2,c01,c11,ss,tt,yes)
         if (yes) then ! stack and rank this solution
            nsol        = nsol + 1
            ssol(nsol)  = ss
            pen(nsol)   = penalty(tt) ! ranking
            cside(nsol) = "N"
         endif   
      endif
      if (diy > 0) then         ! consider south face entrance
         call cross_2Dlines(xy1,xy2,c00,c10,ss,tt,yes)
         if (yes) then ! stack and rank this solution
            nsol        = nsol + 1
            ssol(nsol)  = ss
            pen(nsol)   = penalty(tt) ! ranking
            cside(nsol) = "S"
         endif   
      endif
c     
c     -------- 2) Select the best solution -> (s,face) --------
c      
c     if we end up in a badly ill conditioned case (nsol != 1/2)
c     this is probably due to a very small step attempt.
c     just assign something to (s,face) is used as last resort in this case 
c
      if     (nsol == 2) then ! pick the one with lowest penalty
         if (pen(1) < pen(2)) then
            s    = ssol(1)
            face = cside(1)
         else
            s    = ssol(2)
            face = cside(2)
         endif
      elseif (nsol == 1) then ! just take that one
         s    = ssol(1)
         face = cside(1)
      else                    ! bail out, we should not end here at all ...
         dxy = xy2(1:2)-xy1(1:2)
         write(*,263) sum(abs(dxy)),dix,diy
 263     format("warning: locate_cell_entrance: ",
     +          "bail out by spooky reflection dxy =",
     +          e8.3," di(x,y) =",2i3)
         s    = 0.0  
         if (abs(dix)==1) then
            face = "W"
         else
            face = "N"  ! trashbin for other cases
         endif   
      endif  
      end subroutine locate_cell_entrance



      real function penalty(tval)   ! internal to coast_line_intersection
c     -----------------------------------------------
c     return a penalty >= 0, when (0 < tval < 1) is not satisfied 
c     -----------------------------------------------
      real, intent(in)    :: tval
      if     (tval > 1.0) then 
         penalty =  tval - 1.0
      elseif (tval < 0.0) then
         penalty = -tval
      else
         penalty = 0.0  ! (0 < tval < 1) is true
      endif
      end function penalty




      logical function point_on_land(maybewet, wetpt, assured_wetpt)
c     --------------------------------------------------------------
c     Ensure that maybewet is wet (i.e. is_land(maybewet) == .false.) 
c     by nudging it toward wetpt in very small steps to weed
c     out uncertainty in numerical solution of coast line intersection
c     equations. Perform algortihm in geo coordinates. Do not affect vertical
c     component of maybewet (if provided)
c 
c     The idea is to apply a scaling function scale_factor(i=1..max_steps)
c     where 
c         1) scale_factor(0)         = 1
c         2) scale_factor(max_steps) = 0
c         3) d(scale_factor(0))/di ~ -numerical_resolution 
c
c     Apply the parameterization
c   
c         scale_factor(i) = 1.0 - i*sigma + alpha*i**2
c        
c     where alpha = (sigma*max_steps - 1.0)/max_steps**2 implies 2) is satisfied
c
c  
c     Input:  maybewet      : may or may not be dry (accept both length = 2/3)
c             wetpt         : assumed wet point     (accept both length = 2/3)
c
c     Output: assured_wetpt : a guarantied wet point if maybewet is dry
c                             undefined if maybewet tests wet 
c                             ONLY the horizontal components of maybewet are manipulated         
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
      real, intent(in)    :: maybewet(:), wetpt(:) ! accept both length = 2/3
      real, intent(out)   :: assured_wetpt(:)      ! assumed same length as maybewet
      real                :: dg(2),dg_scaled(2)    ! only manipulate horizontal components
      real                :: scale_factor
      integer             :: istep
      logical             :: keep_going
      
      real, parameter    :: resol         = 1.0e-6   ! relative numerical resultion     
      integer,parameter  :: max_steps     = 1000
      real, parameter    :: alpha = (max_steps*resol-1.0)/max_steps**2    
c     ------------------------------------------------------------      
      point_on_land = is_land(maybewet) ! .true. : needs nudging

      if (.not.point_on_land) return    ! no further processing
c     
c     nudging loop: set assured_wetpt by scaling down difference vector iteratively
c   
      dg            = maybewet(1:2) - wetpt(1:2) ! horizontal vector only
      istep         = 1
      assured_wetpt = maybewet ! transfer vertical component, if present
      keep_going    = .true.   ! we know we need at least one step
      do while (keep_going)
         scale_factor = max(0.0, 1.0 - istep*resol + alpha*(istep**2))
         scale_factor = min(scale_factor, 1.0)
         dg_scaled           = dg*scale_factor
         assured_wetpt(1:2)  = wetpt(1:2) + dg_scaled ! only manipulate horizontal component
         keep_going          = is_land(assured_wetpt)
         if (istep > max_steps) then
            write(*,*) "point_on_land: max_steps exceeded"
            write(*,*) "maybewet         =", maybewet
            write(*,*) "wetpt            =", wetpt
            write(*,*) "is_land(maybewet)=", is_land(maybewet)
            write(*,*) "is_land(wetpt)   =", is_land(wetpt)
            stop ! fatal, no plan B available
         endif
         istep = istep + 1
      enddo
c      write(91,*) istep-1
      

      end function point_on_land


      

      subroutine write_summary()
c     -------------------------------------------
c     internal subroutine for debugging
c     use scope of coast_line_intersection 
c     for in/out variables
c     Can be applied when output variables have 
c     a consistent setting
c     -------------------------------------------
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
      write(*,*) "end of summary"
      write(*,*)  ! put empty line for visual clarity
 622  format("summary of coast_line_intersection: ",
     +        2f8.3, " ->", 2f8.3, " : crossed coast :" 
     +       "ref=", 2f8.3, " hit=", 2f8.3)    
 623  format("summary of coast_line_intersection: ",
     +        2f8.3, " ->", 2f8.3, " : no coast crossed")              
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
