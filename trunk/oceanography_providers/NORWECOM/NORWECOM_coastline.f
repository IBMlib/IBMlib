ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     -----------------------------------------------------------
c     NORWEgian ECOlogical Model (NORWECOM) Oceanography Provider
c     Coastline intersection routine
c     -----------------------------------------------------------
c     $Rev: 212 $
c     $LastChangedDate: 2011-01-24 12:08:32 +0100 (Mon, 24 Jan 2011) $
c     $LastChangedBy: mpay $ 
c
c     This component of the NORWECOM oceanography provider deals with
c     the coastline intersection function. This routine is relatively
c     complex and subject to multiple versions, and thus has been split
c     off from the rest of the grid-related functions. 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine coastline_component_version() 
c     Writes the SVN revision number out for logging purposes
      write(*,*) "Coastline comp: $Rev: 212 $"
      end subroutine


      subroutine coast_line_intersection(geo1, geo2, cross_coast,
     +                                                 georef, geohit) 
c     ------------------------------------------ 
c     Trajectory tracking variant
c
c     This function is a service function assisting in enforcing 
c     coastal boundary conditions on particle steps. geo1 is the current 
c     valid (wet, inside domain) position. geo2 is the next position thati 
c     may or may not be in water.
c     The assumption is that the particle moves in a straight line 
c     between geo1 and geo2 in NORWECOM (x,y,z) coordinates i.e. between
c     the points xyz1 and xyz2. Simlarly, reflection also takes place in 
c     NORWECOM space rather than real space. However, for short displacements, 
c     this approximation is probably ok
c
c     This algorithm solves the problems exactly for the any coastal geometry 
c     represented as regular cell centered squares (without previous 
c     assumptions about  small steps xyz1 -> xyz2, so that any number of 
c     neighbor cells can be crossed)
c 
c     Output: 
c     * anycross: If the straight line xyz1 -> xyz2 has crossed a coast line, 
c       anycross is returned .true. otherwise .false. (i.e. the line 
c       between xyz1 and  xyz2 is all in water). 
c
c     If anycross is true, the following vectors are computed:
c
c     * xyzref is the modified final position, if step from xyz1 to 
c       xyz2 is reflected in the coast line
c     * xyzhit is the first position where the line from xyz1 to xyz2 
c       crosses a coast line. 
c
c     If anycross is false xyzref and xyzhit are assigned NaNs.
c     If xyz1 is dry, a warning is issued and anycross is returned false
c     In this implementation does not handle multiple reflection effects,
c     i.e. it does not test that xyzref is a wet point (the step steps 
c     xyz1 -> xyz2 may become reflected onto land in concave coastal 
c     geometries for large steps). This should be handled by the calling
c     function
c     xyz1 and xyz2 must have same length (2 or 3) so they are vector addable
c     
c     Modified to fit NORWECOM from conceptual revision by ASC/MPA Nov 24, 2010 
c     ------------------------------------------ 
      real, intent(in)     :: geo1(:),geo2(:)
      logical, intent(out) :: cross_coast
      real, intent(out)    :: georef(:), geohit(:)
c     -----locals-----------
      logical              :: leaving,loopflag
      character*1          :: reftype, xface
      real                 :: direct(2), reflect(2)
      real                 :: xyz1(2), xyz2(2),s,xyzhit(2),xyzref(2)
      integer              :: ix,iy,ix2,iy2,refdim
c     ------------------------------------------ 
c     Switch to grid units, xy1, xy2 and get cell indices 
      call lonlat2xy(geo1(1:2),xyz1)
      call lonlat2xy(geo2(1:2),xyz2)
      call cell_index(xyz1,ix,iy)
      call cell_index(xyz2,ix2,iy2)

      !Wet assert xy1 is wet - check it
      if(verbose>0) write (*,*) "xyz1 is dry?",land(ix,iy)
c     Test if xyz1 and xyz2 are in same cell => no coast line crossing
c     most cases terminate here
c
      if(verbose>0) write(*,*) "xyz1, ix ,iy  :", xyz1,ix, iy
      if(verbose>0) write(*,*) "xyz2, ix2,iy2 :", xyz2,ix2, iy2
      if ((ix==ix2).and.(iy==iy2)) then
         cross_coast = .false.
         if (verbose>0) write(*,*) "xyz1 -> xyz2 in same cell"
         return
      endif
      if (verbose>0) write(*,*) "trajectory analysis: begin"

c     Now we know that we are leaving the cell. But where? 
      call assess_cell_leaving(xyz1,xyz2,ix,iy,leaving,s,xface) 

c     Follow sequence moving from one cell to the next until we either
c     reach the end of the vector, or encounter land
      cross_coast = .false. ! default
      loopflag  = .true.
      do while(loopflag)
         if(verbose>0)   write(*,421) ix,iy,xface
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
            write (*,*) "coast_line_intersection: Unknown exit face :",
     +           xface
            stop
         end select
         !Is the next cell wet or dry? If its dry, then we've hit the
         !coast.  
         if (verbose>0) then
            if(land(ix,iy)) then
               write(*,422) ix,iy,"dry"
            else
               write(*,422) ix,iy,"wet" 
            endif
         endif
         if(land(ix,iy)) then   !have crossed onto land
           cross_coast = .true. !so setup for exit
           loopflag    = .false.        
         else   !do we leave this new cell?
           call assess_cell_leaving(xyz1,xyz2,ix,iy,leaving,s,xface) 
           loopflag = leaving     !loop if we leave this new cell,otherwise done
         endif
      enddo
421   format("leaving  ",2i4," via ", a," face") 
422   format("entering ",2i4," which is ", a) 
c  
c     if (cross_coast == .true.) continue and resolve xyzref and xyzhit
c     else let xyzref and xyzhit remain undefined
      if (.not.cross_coast) then
         geohit = NaN
         georef = NaN
         return
      endif
c
c     If the coast line is crossed, compute xyzref, xyzhit
c     0 < s_min < 1 is the coordinate along the vector (xyz2-xyz1)
c          
      direct  = xyz2-xyz1
      reflect = direct
      if     (xface=="N" .or. xface =="S") then
         refdim = 2 ! mirror the y dimension 
      elseif (xface=="E" .or. xface =="W") then
         refdim = 1 ! mirror the x dimension
      else
         stop "coast_line_intersection: unhandled reflection type"
      endif
      reflect(refdim) = -reflect(refdim) 
c     |direct| == |reflect|        
      xyzhit = xyz1   +     s*direct ! position where coast line is hit
      xyzref = xyzhit + (1-s)*reflect ! position for reflection in coast line 

c     We need to be aware of issues caused by numerical precision
c     If the endpoint is very close to a boundary, we can potentially
c     end up in a grey zone of unpredictable behaviour caused by
c     numerical precision - especially here in NORWECOM where we 
c     make a series of transforms between units. Thus, force the distance
c     moved after reflection to more than ~1m (1e-4 grid units) away from the
c     axis about which we are reflecting.
      if(abs((1-s)*reflect(refdim)) < 1e-4) then  ! we're too close
         if(verbose>0) then 
           write(*,*) "Particle too close to coast"
           write(*,*) "Dist       : ",abs((1-s)*reflect(refdim))  
           write(*,*) "Reflect    : ", reflect
           write(*,*) "xyzhit     : ", xyzhit
           write(*,*) "xyzref     : ", xyzref
         endif
         xyzref(refdim) = xyzhit(refdim)  
     +                     - 1e-4*direct(refdim)/abs(direct(refdim))
         if(verbose>0)  write(*,*) "xyzref new : ", xyzref
      endif

c     Switch to back to geo units, geohit, georef --- 
      call xy2lonlat(xyzhit,geohit(1:2))
      call xy2lonlat(xyzref,georef(1:2))

c     Finally, calculate the vertical aspects
c     The reflected point georef is not not effeted by reflections in
c     the 2D plane. Geohit(3) is the fraction of the distance moved
c     in the vertical plane
      georef(3) = geo2(3)      !Not effected by 2d reflections
      geohit(3) = (geo2(3)-geo1(3))*s + geo1(3) !Hits part of the way along movement

      end subroutine 


      subroutine assess_cell_leaving(xyz1,xyz2,ix,iy,leave,s,exitface) 
c     --------------------------------------------------------------
c     Determine whether and where the direct line from xyz1 to xyz2 
c     leaves the cell with center (ix,iy)
c     If so return leave=.true. (else .false.)
c     If so return leave=.true., determine 0 < s < 1 which is the 
c     coordinate along the vector (xyz2-xyz1) where the vector exits the
c     cell. Also identify the exit face as one of "N","S","E","W" 
c     so the proper reflection can be calculated
c     If leave=.false. (s,exitface) are set to the default values, 1.0
c     and "@"
c     --------------------------------------------------------------
      real, intent(in)         :: xyz1(:),xyz2(:)
      integer, intent(in)      :: ix,iy
      logical, intent(out)     :: leave
      real, intent(out)        :: s
      character*1, intent(out) :: exitface
c
      logical                  :: yes
      real                     :: stest,dxyz(2)
      real                     :: c00(2),c01(2),c10(2),c11(2)
c     --------------------------------------------------------------
      call get_horiz_corners(ix,iy,c00,c01,c10,c11) ! resolve faces          
c
c     Determine how the step leaves the box defined by (c00,c01,c10,c11))
c     We can make efficiency gains by taking advantage of the direction
c     that the particle is travelling in - a particle travelling NE
c     can only leave by the North and East faces etc
c
      leave = .false.
      exitface = "@"
      s     = 1.0 
      dxyz = xyz2(1:2)-xyz1(1:2)
c     First the east-west direction
      if(dxyz(1)>0) then    !we're moving to the east
         call cross_2Dline_segments(xyz1,xyz2,c10,c11,stest,yes)
         if (yes.and.(stest<=s)) then ! rely on right-to-left evaluation
            leave = .true.
            s     = stest
            exitface  = "E"   
         endif
      else   !we're moving to the west
         call cross_2Dline_segments(xyz1,xyz2,c00,c01,stest,yes)
         if (yes.and.(stest<=s)) then ! rely on right-to-left evaluation
            leave = .true.
            s     = stest
            exitface  = "W"   
         endif
      endif
      if(dxyz(2) >0) then !we're moving to the North
         call cross_2Dline_segments(xyz1,xyz2,c01,c11,stest,yes)
         if (yes.and.(stest<=s)) then ! rely on right-to-left evaluation
            leave = .true.
            s     = stest
            exitface  = "N"   
         endif
      else  !we're moving to the south
         call cross_2Dline_segments(xyz1,xyz2,c00,c10,stest,yes)
         if (yes.and.(stest<=s)) then ! rely on right-to-left evaluation
            leave = .true.
            s     = stest
            exitface  = "S"   
         endif
      endif
      end subroutine assess_cell_leaving     ! local subroutine

