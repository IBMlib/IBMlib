cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     Polygons module
c     ---------------------------------------------------
c     $Rev: 274 $
c     $LastChangedDate: 2011-02-10 16:41:38 +0100 (Thu, 10 Feb 2011) $
c     $LastChangedBy: asch $ 
c
c     Provides basic functionality about polygons
c     
c     Enable speed-up in processing by classifying polygons at instantiation
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module polygons
      use geometry
      implicit none
      private

      public :: lonlat_polygon  
      public :: init_lonlat_polygon   ! generic constructor
      public :: init_alligned_rectangle ! variant constructor for alligned rectangles, allow speed-up
      public :: delete_lonlat_polygon
      public :: print_lonlat_polygon
      public :: get_polygon_name
      public :: get_bounding_box
      public :: is_inside_polygon
      public :: polygon_member_mask
      public :: get_ICESrectangle_sw_corner     ! TAG     -> LON,LAT
C     public :: get_surrounding_ICESrectangle   ! LON,LAT -> TAG
     
c     --------------------------------
c     Polytype polygon
c     --------------------------------
      type lonlat_polygon     
      private
         logical       :: is_alligned_rectangle ! polygon is flagged as rectangle alligned with lon,lat axes
         real, pointer :: nodes(:,:)        ! input (2,nnode) implicit closure between first and last point  
         character*256 :: name              ! input: auxillary name   
         real          :: bounding_box(2,2) ! (SWlon,SWlat,NElon,NElat) of surrounding bounding box; inferred from nodes
         real          :: extpt(2)          ! auxillary external point; inferred from nodes
      end type


c     ============================================================
                          contains
c     ============================================================


      subroutine get_ICESrectangle_sw_corner(icesrec, lon_sw, lat_sw)
c     ------------------------------------------------------------------
c     Resolve the SW corner of an ICES statistical rectangle from tag LLPP
c
c     LL = 09 -> lat =  40       latitude  step / rectangle = 0.5 degrees
c     PP = B0 -> lon = -40       longitude step / rectangle = 1.0 degrees
c     
c
c     Notes: * PP = A? is iregular so it is treated as special case   
c            * Letter in PP must be upper case
c            * This implementation makes analytical continuation onto land
c              of rectangle designation, so it will not complain about 
c              fiktive rectangles on land. All wet rectangles are correct.
c     ------------------------------------------------------------------
      character(len=4), intent(in) :: icesrec
      real, intent(out)            :: lon_sw, lat_sw
      character                    :: lonhead 
      integer                      :: latix, lonminor
c     ----------------------------------------------------------------
      read(icesrec(1:2),*) latix 
      lonhead = icesrec(3:3)
      read(icesrec(4:4),*) lonminor
      lat_sw = 40.0 + (latix-9)*0.5   ! LL = 09 -> lat(sw) = 40 ; lat step = 0.5  
      if (lonhead == "A") then        ! A? is iregular - only covers A0,A1,A2,A3
         lon_sw = -36 + 1.0*lonminor  ! lon step = 1.0
      else
         lon_sw = -40 + 10*(ichar(lonhead)-ichar("B")) + 1.0*lonminor ! lon step = 1.0
      endif
      end subroutine get_ICESrectangle_sw_corner

      
      
      subroutine init_lonlat_polygon(poly, nodes, name)
c     ----------------------------------------------------
c     Generic polygon constructor
c
c     TODO: add a numerically robust algorithm to spot and flag alligned rectangles
c     ----------------------------------------------------
      type(lonlat_polygon), intent(inout) :: poly
      real, intent(in)                    :: nodes(:,:) ! (2+,nnode)  ignore vertical coordinate, if present
      character*(*), intent(in), optional :: name 
c
      integer                             :: n1,n2,icut
c     ----------------------------------------------------
      n1 = size(nodes, 1)
      n2 = size(nodes, 2)
      if (n1<2) then
         write(*,*) "init_lonlat_polygon: space dim < 2 ", n1
         stop
      endif
      if (n2<3) then
         write(*,*) "init_lonlat_polygon: too few nodes ", n2
         stop
      endif
      
      allocate(poly%nodes(2,n2))
      poly%nodes = nodes(1:2,:)  ! copy nodes; ignore vertical coordinate, if present
c     ------ polygon classification flags ------
      poly%is_alligned_rectangle = .false. ! only for processing speed-up
c     ------ infer bounding box ------
      poly%bounding_box(1,1) = minval(poly%nodes(1,:)) ! SWlon
      poly%bounding_box(2,1) = minval(poly%nodes(2,:)) ! SWlat
      poly%bounding_box(1,2) = maxval(poly%nodes(1,:)) ! NElon
      poly%bounding_box(2,2) = maxval(poly%nodes(2,:)) ! NElat
c     ------ infer auxillary exterior point ------
      poly%extpt = poly%bounding_box(:,2) + 1.0  ! add arbitrary safety dist to NE corner to get an external point
c     ------ set name
      if (present(name)) then
         icut = min(len(poly%name), len(name))
         poly%name = name(1:icut) 
         poly%name = adjustl(poly%name)   
      else
         poly%name = "no_name" 
      endif 
c     ----------------------------------------------------      
      end subroutine init_lonlat_polygon
      



      subroutine init_alligned_rectangle(poly, corners, name)
c     ----------------------------------------------------
c     Variant constructor for alligned rectangles, allow speed-up
c
c     Opposed to init_lonlat_polygon, only a set of opposite corners
c     should be provided at instantiation
c
c     TODO: make a variant that accepts vector corners = (SWlon, SWlat, NElon, NElat)
c           along with an interface
c     ----------------------------------------------------
      type(lonlat_polygon), intent(inout) :: poly
      real, intent(in)                    :: corners(:,:) ! (2+,2)  ignore vertical coordinate, if present
      character*(*), intent(in), optional :: name 
c
      integer                             :: n1,n2,icut
c     ----------------------------------------------------
      n1 = size(corners, 1)
      n2 = size(corners, 2)
      if (n1<2) then
         write(*,*) "init_alligned_rectangle: space dim < 2 ", n1
         stop
      endif
      if (n2/=2) then
         write(*,*) "init_alligned_rectangle: nodes /= corners ", n2
         stop
      endif
c     ------ infer bounding box from corners ------
      poly%bounding_box(1,1) = minval(corners(1,:)) ! SWlon
      poly%bounding_box(2,1) = minval(corners(2,:)) ! SWlat
      poly%bounding_box(1,2) = maxval(corners(1,:)) ! NElon
      poly%bounding_box(2,2) = maxval(corners(2,:)) ! NElat
c     ------ assign nodes counter clockwise from bounding box ------ 
c            ignore vertical coordinates, if present
      allocate(poly%nodes(2,4))
      poly%nodes(:,1) =  poly%bounding_box(:,1) ! SW corner
      poly%nodes(1,2) =  poly%bounding_box(1,2) ! SE corner, lon
      poly%nodes(2,2) =  poly%bounding_box(2,1) ! SE corner, lat
      poly%nodes(:,3) =  poly%bounding_box(:,2) ! NE corner
      poly%nodes(1,4) =  poly%bounding_box(1,1) ! NW corner, lon
      poly%nodes(2,4) =  poly%bounding_box(2,2) ! NW corner, lat
c     ------ polygon classification flags ------
      poly%is_alligned_rectangle = .true. ! only for processing speed-up
c     ------ infer auxillary exterior point ------
      poly%extpt = poly%bounding_box(:,2) + 1.0  ! add arbitrary safety dist to NE corner to get an external point
c     ------ set name
      if (present(name)) then
         icut = min(len(poly%name), len(name))
         poly%name = name(1:icut) 
         poly%name = adjustl(poly%name)   
      else
         poly%name = "no_name" 
      endif 
c     ----------------------------------------------------      
      end subroutine init_alligned_rectangle




      subroutine delete_lonlat_polygon(poly)
c     ----------------------------------------------------
      type(lonlat_polygon), intent(inout) :: poly
c     ----------------------------------------------------
      if (associated(poly%nodes)) deallocate(poly%nodes) 
      poly%bounding_box = -999.0
      poly%extpt        = -999.0
      poly%name         = "deleted"
c     ----------------------------------------------------
      end subroutine delete_lonlat_polygon



      subroutine print_lonlat_polygon(poly,iunit)
c     ----------------------------------------------------
c     Print lonlat_polygon instance to logical unit iunit
c     ----------------------------------------------------
      type(lonlat_polygon), intent(in) :: poly
      integer, intent(in)              :: iunit
      integer                          :: ino, nn
c     ----------------------------------------------------
      nn = size(poly%nodes, 2)
      write(iunit,533) trim(adjustl(poly%name)), nn
      write(iunit,535) poly%is_alligned_rectangle
      write(iunit,540) "SW", poly%bounding_box(:,1)
      write(iunit,540) "NE", poly%bounding_box(:,2)
      write(iunit,542) poly%extpt
      write(iunit,550) "node", "degE", "degN"
      do ino = 1, nn
         write(iunit,551) ino, poly%nodes(:,ino)
      enddo
c     -------------
 533  format("lonlat_polygon instance ", a, " with ", i5, " nodes")
 535  format("is_alligned_rectangle flag = ", l8)
 540  format("bounding box ",a, " corner = ", 
     +       f10.4," degE ", f10.4," degN")  
 542  format("auxillary external point", 
     +       f10.4," degE ", f10.4," degN")    
 550  format(3a10)
 551  format(i10,2f10.4)
c     ----------------------------------------------------
      end subroutine print_lonlat_polygon




      subroutine get_polygon_name(poly, polyname)
c     ----------------------------------------------------
c     Print lonlat_polygon instance to logical unit iunit
c     ----------------------------------------------------
      type(lonlat_polygon), intent(in) :: poly
      character*(*), intent(out)       :: polyname
      integer                          :: buflen, namlen
c     ----------------------------------------------------
      namlen = len(poly%name)
      buflen = len(polyname)
      if (namlen>buflen) then
         polyname = poly%name(1:buflen)
      else
         polyname = poly%name
      endif
c     ----------------------------------------------------
      end subroutine get_polygon_name

      


      subroutine get_bounding_box(poly,BB)
c     ------------------------------------------------------------
c     return bounding_box as BB = SWlon, SWlat, NElon, NElat
c     ------------------------------------------------------------
      type(lonlat_polygon), intent(in) :: poly
      real, intent(out)                :: BB(:)
      BB(1:4) = reshape(poly%bounding_box, (/4/))
      end subroutine get_bounding_box


      
      logical function is_inside_bounding_box(poly, xy)
c     ----------------------------------------------------
c     Test whether xy is inside bounding box of polygon:
c
c        (SWlon < x < NElon) .and. (SWlat  < y < NElat)
c
c     This function could be moved to the public interface, 
c     if useful outside
c     ----------------------------------------------------
      type(lonlat_polygon), intent(in) :: poly
      real, intent(in)                 :: xy(:)  ! size = 2+;  ignore vertical coordinate, if present
c
      is_inside_bounding_box = ((poly%bounding_box(1,1) < xy(1)).and.
     +                          (poly%bounding_box(1,2) > xy(1)).and.
     +                          (poly%bounding_box(2,1) < xy(2)).and.
     +                          (poly%bounding_box(2,2) > xy(2)))
      end function is_inside_bounding_box



      logical function is_inside_polygon(poly, xy)
c     ----------------------------------------------------
c     Test whether xy is inside polygon or not, using the ray tracing algortihm  
c     ray tracing algortihm: odd number of crossings of polygon side going from xy to exterior point => xy inside 
c   
c     If flag is_alligned_rectangle is .true. redirect to is_inside_bounding_box which is faster   
c     TODO: update scan over line segments with faster algorithm which 
c           uses a smarter external point - the numerical robustness
c           needs to be assessed though
c     ----------------------------------------------------
      type(lonlat_polygon), intent(in) :: poly
      real, intent(in)                 :: xy(:)  ! size = 2+;  ignore vertical coordinate, if present
c
      integer                          :: ncross, nnodes, inode
      real                             :: s
      logical                          :: cross
c     ----------------------------------------------------
      if (poly%is_alligned_rectangle) then
         is_inside_polygon = is_inside_bounding_box(poly, xy)  ! faster if applicable
         return
      endif
c     ------- apply generic ray tracing algortihm for inside test 
      ncross = 0
      nnodes = size(poly%nodes, 2)
      do inode = 1, nnodes-1
         call cross_2Dline_segments(poly%nodes(:,inode),
     +                              poly%nodes(:,inode+1),
     +                              xy, poly%extpt, s, cross)
         if (cross) ncross = ncross + 1
      enddo   
c     finally include implicit closure between final and first point
      call cross_2Dline_segments(poly%nodes(:,nnodes),
     +                           poly%nodes(:,1),
     +                           xy, poly%extpt, s, cross)
      if (cross) ncross = ncross + 1
      is_inside_polygon = (mod(ncross,2) == 1)
c     ----------------------------------------------------    
      end function is_inside_polygon



      subroutine polygon_member_mask(poly, xyset, mask)
c     ----------------------------------------------------
c     Set mask(i) on wheter xyset(:,i) is inside poly
c     Vectorized variant of is_inside_polygon
c     ----------------------------------------------------
      type(lonlat_polygon), intent(in) :: poly
      real, intent(in)                 :: xyset(:,:)  ! size = (2+,npt);  ignore vertical coordinate, if present
      logical, intent(out)             :: mask(:)     ! assume sufficiently large, not checked
      integer                          :: ipt,npt
c     ----------------------------------------------------
      npt = size(xyset, 2)
      if (size(mask)<npt) then
         write(*,*) "polygon_member_mask: mask size insufficient"
         stop
      endif
      do ipt = 1, npt
         mask(ipt) = is_inside_polygon(poly, xyset(:,ipt))
      enddo
      
c     ----------------------------------------------------
      end subroutine polygon_member_mask

      end module


c$$$      program test
c$$$c     ----------------------------------------------------
c$$$c     ifort geometry.o polygons.f 
c$$$c     ----------------------------------------------------      
c$$$      use polygons
c$$$      real :: nodes(3,4)
c$$$      type(lonlat_polygon) :: apoly
c$$$      nodes(3,:)   = 2.0
c$$$      nodes(1:2,1) = (/0.0, 0.0/)
c$$$      nodes(1:2,2) = (/2.0, 0.0/)
c$$$      nodes(1:2,3) = (/1.0, 1.1/)
c$$$      nodes(1:2,4) = (/0.0, 0.8/)
c$$$      call init_lonlat_polygon(apoly, nodes, "sldifh")
c$$$      write(*,*) is_inside_polygon(apoly, (/ 0.5,  0.5/))
c$$$      write(*,*) is_inside_polygon(apoly, (/ 5.0,  0.5/))
c$$$      write(*,*) is_inside_polygon(apoly, (/-5.0, -3.0/))
c$$$      end program

