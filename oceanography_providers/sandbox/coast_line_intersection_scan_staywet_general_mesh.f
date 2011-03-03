cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Implementation horizontal_representation pieces
c      
c     Stabilize the "staywet" algorithmic pairs and
c     study the stability of the interplay between 
c     coast_line_intersection and multiple_reflection_path.
c
c     This version abandon the coast line cushin approach, but rather focuses 
c     on the consistency in decisions, essentially
c
c       1) is_land is the authoritative 
c       2) if coast_line_intersection finds a crossing, then
c          hit point must test wet by island (but reflection point may be wet/dry)
c       3) if multiple_reflection_path finds a crossing, then
c          both first hit point and reflection point must be wet
c
c     There was very rare problems with the version (coast_line_intersection_staywet_general_mesh.f)
c     that resolves the intercell trajectory between initail and final cell.
c     This version rather scans ALL potential intermediate cell to avoid
c     decision conflict problems due to final numerical accuracy. The problems are 
c     related to very small steps and limit cases. This "scan" algorith should be 
c     rather efficient though, if a reasonable step size is applied.
c     
c     Needs: 
c          horizontal_grid_transformations.mod
c          geometry.mod
c          constants.mod (in geometry)
c 
c     Prepare:
c        ln -s  ../../constants.f 
c        ln -s  ../../geometry.f 
c        ln -s ../generic_elements/horizontal_grid_transformations/horizontal_grid_transformations_X.f horizontal_grid_transformations.f
c        ifort -c -e90 constants.f geometry.f horizontal_grid_transformations.f 
c     where X is the grid type (e.g. X==lonlat)
c
c     Compile:
c        ifort -e90 coast_line_intersection_scan_staywet_general_mesh.f constants.o geometry.o horizontal_grid_transformations.o
c       
c     ASC, 25Feb2011
c     UPDATED:
c           coast_line_intersection  ->  horizontal_representation_scan_stay_wet.f
c           multiple_reflection_path ->  particle_tracking.f 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module horiz_repres_pieces            ! implements sub module of horizontal_representation
      use horizontal_grid_transformations   ! apply full module
      use geometry                          ! cross_2Dlines
      implicit none

      private           ! default scope

      public :: nx,ny                            ! reexport from horizontal_grid_transformations
      public :: init_horizontal_representation   ! init  this module
      public :: set_verbose_horiz_repres
      public :: horizontal_range_check           ! reexport from horizontal_grid_transformations
      public :: is_land                          
      public :: coast_line_intersection 
      
c  
      integer, allocatable,public :: wetmask(:,:) ! auxillary to define coastal geometry (1: wet nodes, 0 dry) 
      integer         :: verbose = 0

      contains 
      

      subroutine init_horizontal_representation()
c     ------------------------------------------------------
      allocate(wetmask(nx,ny))  
      end subroutine init_horizontal_representation
c      


      subroutine set_verbose_horiz_repres(ival)     ! COPY2 ->  particle_tracking.f 
      integer,intent(in) :: ival
      verbose = ival
      end subroutine


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
         all_wet = .true.         ! implied since (ix1,ix2) is wet
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
c     Infinitesimal step exit point: capture very small intercell steps
c     at the numerical resolution limit, where intercell transition
c     points can not be resolved (consider this as a coastal osculation case)     
c   

      if ((wetmask(ix2,iy2) == 0).and.
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
c     Only perform final cell analysis if final cell is dry and 
c     no intermediate crossings. An intermediate crossing is 
c     logically ranked before a final cell.
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
c     The coast line was crossed geo1->geo2, resolve xyref and xyhit
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
      real             :: xy(2)
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


      end module ! horiz_repres_pieces



c############################################################################################      


      module particle_tracking_pieces
      use horiz_repres_pieces 
      implicit none

      integer,private :: verbose=0

      contains
      

      subroutine set_verbose_particle_tracking(ival)                      ! COPY2 ->  particle_tracking.f 
      integer,intent(in) :: ival
      verbose = ival
      end subroutine


      subroutine multiple_reflection_path(s0, s1, anycross, sref, shit0)  ! COPY2 ->  particle_tracking.f 
c-------------------------------------------------------------------------
c     Compute key points of multiple horizontal coastal reflection path
c     by iterative application of coast_line_intersection primitive
c     when trying to move from (valid, wet) s0 to new position s1
c     If a coast line is crossed return anycross == .true.
c     Then sref is the (multiple) horizontal coastal reflected end point
c     which is guarentied wet. shit0 is the point of (first) coastal intersection
c     and shit0 is guarentied wet.
c     If anycross == .false. (sref, shit0) are left undefined (unaltered)
c     multiple_reflection_path has same interface as coast_line_intersection
c   
c     ASC 10Feb2011: improved verbose logging + race condition trap
c                    fixed exit condition bug for higher reflections
c                    added trace to enable debugging + transparency
c-------------------------------------------------------------------------
      real,intent(in)      :: s0(:), s1(:)
      logical,intent(out)  :: anycross
      real,intent(out)     :: sref(:)   ! buffer assumed same size as s0/s1
      real,intent(out)     :: shit0(:)  ! buffer assumed same size as s0/s1

      integer, parameter   :: max_reflections = 20 ! otherwise time step too long ...
      real                 :: refs(size(s0),max_reflections) ! first index = fastest 
      real                 :: hits(size(s0),max_reflections) ! first index = fastest 
      logical              :: isla, latercross
      real                 :: ds  
      integer              :: i,k
      real, parameter      :: race_limit = 1.0e-5  ! step limit for flagging race condition
c---------------------------------------------------------------
      if (verbose>0) write(*,*) " multiple reflection analysis begin"

      i  = 1
      call coast_line_intersection(s0,s1,anycross,refs(:,i),hits(:,i))

      if (verbose>0) then
         if (anycross) then 
            write(*,422) i,s0(1:2),s1(1:2),refs(1:2,i),hits(1:2,i) 
            write(*,*) "coastal reflection detected"
            isla = is_land(refs(1:2,i))
            if (isla) then
               write(*,*) "reflected point point dry: "//
     +                    "continue multiple reflection analysis"
            else
               write(*,*) "reflected point point wet: "//
     +                    "check reflection path is wet"
            endif ! isla
         else
            write(*,423) i, s0(1:2), s1(1:2)
            write(*,*)"no coastal crossing: "//
     +                "multiple reflection analysis end"
         endif
      endif ! verbose

      if (.not.anycross) return  ! (sref, shit0) undefined
      
c
c     --- we hit the coast line, start multiple reflection analysis
c         (this means (sref, shit0) are defined and anycross == .true.)
c         
c
      latercross = .true.    ! enter multi loop and keep anycross == .true. as exit condition
      shit0      = hits(:,1) ! first hit defined when anycross == .true.
c
c     Recursively apply coast_line_intersection to (hit,reflect) pairs to 
c     end up with a sub path that does not cross the coast line (which
c     implies that reflect is wet, because hit is wet). Note that it is not 
c     sufficient to check that reflect is wet, because step may jump 
c     over land (wet-to-wet)
c
      do while (latercross .and. (i < max_reflections)) ! left-to_right evaluation
         i  = i + 1       
         call coast_line_intersection(hits(:,i-1),refs(:,i-1),
     +                                latercross, refs(:,i),hits(:,i)) 

         if (verbose>0) then
            if (latercross) then 
               write(*,422) i, hits(1:2,i-1), refs(1:2,i-1), 
     +                      refs(1:2,i), hits(1:2,i)
               write(*,*) "is_land(ref point) = ", is_land(refs(1:2,i))
            else
               write(*,423) i, hits(1:2,i-1), refs(1:2,i-1)
            endif
         endif ! verbose
      enddo
c
c     If multi reflection loop was timed out, check for race condition
c     If race condition is detected set sref = shit (last wet point on trajectory)
c     When timed out the set (hits(1:max_reflections), refs(1:max_reflections))
c     are defined, because step = max_reflections resulted in latercross = .true.
c
      if (i >= max_reflections) then 
         sref = hits(:,max_reflections) ! last wet point in analysis
         k    = max_reflections
         ds   = sum(abs(hits(1:2,k)-hits(1:2,k-1))) 
         if (verbose>0) then
            write(*,*) "multiple reflection analysis: "//
     +                 "race condition trapped"
            write(*,*) "multiple reflection analysis: "//
     +                 "return sref = last coastal hit = ", sref
            write(*,*) "race condition parameter ds = ", ds
         endif
         return  ! max_reflections exceeded exit point
      endif
c
c     If latercross == .false. last tested sub path (hits(i-1), refs(i-1))
c     was all in water and therefore refs(i-1) is the final reflection
c
      if (.not.latercross) then ! 
         sref = refs(:,i-1)
      else
         write(*,*) "multiple_reflection_path: unexpected"
         stop
      endif

      if (is_land(sref)) then ! sref defined
         write(*,*) "multiple_reflection_path: assertion failed"
         write(*,*) "unexpected dry test: island(sref) = ",is_land(sref)
         stop ! no plan B
      endif

      if (verbose>0) write(*,*) " multiple reflection analysis end"

 422  format("multiref step ", i3, " :", 2f11.6, " ->", 2f11.6, 
     +       " : crossed coast   ref=", 2f11.6, " hit=", 2f11.6)    
 423  format("multiref step ", i3, " :", 2f11.6, " ->", 2f11.6, 
     +       " : no cross")

c     ------------------------
      contains
c     ------------------------

      subroutine write_trace(last_valid)
c     ---- uses scope of multiple_reflection_path ----
      integer, intent(in) :: last_valid
      integer             :: j
      write(*,*) "multiple_reflection_path: trace"
      write(*,*) "starting point s0 = ", s0, " is_land=", is_land(s0)
      write(*,*) "ending   point s1 = ", s1, " is_land=", is_land(s1)
      write(*,*) "any coastal crossing on path s0->s1:", anycross
      write(*,*) "step       hit_point        dry?"//
     +               "       ref_point        dry?"
      do j=1, last_valid
         write(*,432) j, hits(1:2,j), is_land(hits(1:2,j)), 
     +                   refs(1:2,j), is_land(refs(1:2,j))
      enddo
 432  format(i3, " :", 2f11.6, l4, 2x, 2f11.6, l4)
           

      end subroutine write_trace ! internal to multiple_reflection_path
      end subroutine multiple_reflection_path


      end module ! particle_tracking_pieces


c############################################################################################


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
      use horiz_repres_pieces 
      use particle_tracking_pieces
      use geometry
      use constants
      implicit none
      real    :: geo1(3), geo2(3), georef(3), geohit(3)
      real    :: dgeo(2),x,y,s
      logical :: anycross,ldum
      real    :: angle,r2(2),psc,dpsc
      integer :: ixc,iyc,ix,iy,ityp
      real,parameter :: rjump = 0.5  ! Bee jump range (lon,lat) in degrees
      real    :: stability
      integer(kind=8) :: number_of_steps, istep, nextstop,logging_step ! allow very large number of steps
      integer :: verbose = 0  ! for main program
c     ----------------------------------------------------------------
      real, parameter    :: lambda1    = -4.0 ! -4.0 for backward comparisons
      real, parameter    :: lambda_max = 10.0 !  only for determining nx
      real, parameter    :: dlambda    =  1.0 !  1.0 for backward comparisons
      real, parameter    :: phi1       = 49.0 ! 49.0 for backward comparisons
      real, parameter    :: phi_max    = 60.0 !  only for determining ny
      real, parameter    :: dphi       =  1.0 !  1.0 for backward comparisons
c     ----------------------------------------------------------------
      nx = 1 + int((lambda_max - lambda1)/dlambda) ! set directly in horizontal_grid_transformations
      ny = 1 + int((phi_max    - phi1   )/dphi)    ! set directly in horizontal_grid_transformations
      call init_horiz_grid_transf(lambda1,phi1,dlambda,dphi) 
      call init_horizontal_representation()

c$$$c     ------------------------ 
c$$$c     ---- MPA test block ----
c$$$c     ------------------------ 
c$$$      geo1 = (/2.0, 53.0, 0.0/)
c$$$      call get_horiz_ncc_index(geo1,ixc,iyc)     
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
      call get_horiz_ncc_index(geo1,ixc,iyc)     
      !Surround the point with dry cells
      wetmask   = 0 ! 0 dry, 1 wet
      wetmask(ixc,iyc)     = 1  ! 
      wetmask(ixc+1,iyc)   = 1  ! 
      wetmask(ixc+1,iyc+1) = 1  ! 

c     -- number of steps / print frequency --
      psc      = 1
      dpsc     = 0.3
      nextstop = nint(10**psc, kind=8)  ! -1 for no print
      number_of_steps = nint(10**stability, kind=8)
      
      do istep = 1, number_of_steps

         if (verbose>0) write(*,33) istep
 33      format(25("="),1x,i8,1x,25("="))

c         if (mod(istep,100)==0) write(*,*) istep
c         if (istep==737300) call set_verbose_particle_tracking(1)

         if (istep == nextstop) then
            write(*,*) istep
            psc      = psc + dpsc
            nextstop = nint(10**psc, kind=8) 
         endif
         
c        -------- generate a RW step --------       
         call random_number(r2)
         angle = 8*atan(1.0)*r2(2) ! 0 < angle < 2*Pi
         dgeo  = r2(1)*rjump*(/cos(angle), sin(angle)/)
         geo2(1:2)  = geo1(1:2) + dgeo
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
         
c         write(65,*) geo1(3)
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
         
c$$$c        --- dump from a specific step ---
c$$$         logging_step = 1088048
c$$$         if (istep==(logging_step-1)) then
c$$$            call set_verbose_particle_tracking(1) 
c$$$            call set_verbose_horiz_repres(1) 
c$$$            verbose = 1
c$$$         endif
c$$$         if (istep==1088049) stop 9876
         

      enddo



      end program

      
