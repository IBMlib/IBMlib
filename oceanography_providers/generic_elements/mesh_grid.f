      module regular_lonlat_grid
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This module holds core arrays and data transformations
c     representing horizontal regular lon-lat grids and various 
c     vertical grids, e.g. sigma and z-grid with open/closed surface
c     
c     This module provides the interface physical_fields except:
c
c       init_physical_fields(time)
c       close_physical_fields()
c       update_physical_fields(time, dt)
c
c     The module arrays are allocation/deallocation) by:
c
c       init_regular_lonlat_grid
c       close_regular_lonlat_grid
c
c     The module arrays/grid descriptors have public scope.
c     The module arrays are supposed to be updated by the client
c     module, when update_physical_fields are invoked in the client
c     module.

      update doc
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use time_tools           ! import clock type
      use constants
      use geometry
      use array_tools

      implicit none
      private           ! default scope

      public :: init_regular_lonlat_grid   ! init  this module
      public :: close_regular_lonlat_grid  ! close this module

     

      public :: get_master_clock
      public :: set_master_clock
      
      public :: interpolate_turbulence
      public :: interpolate_turbulence_deriv
      public :: interpolate_currents
      public :: interpolate_temp

c      public :: interpolate_salty   ! currently unused
c      public :: interpolate_wind    ! currently unused
      public :: interpolate_zooplankton
      public :: interpolate_wdepth
      public :: is_wet    
      public :: is_land
      public :: horizontal_range_check
      public :: coast_line_intersection
      
c     ---------- other exports ----------

      public :: get_ncc_coordinates

c     -------------------- module data --------------------  
      
      type(clock), target,public :: master_clock

      integer, parameter :: verbose = 0  ! debugging output control
      real, parameter    :: htol = 1.e-6 ! tolerance for surface/bottom


c     ------ grid dimensions:   ------
c     
      integer,public :: nz         ! basic grid dimensions - (nx,ny) from horiz grid provider)


c     --- 3D grids ---
             
      real,allocatable,public :: u(:,:,:)          ! u of current [m/s] (positive east)
      real,allocatable,public :: v(:,:,:)          ! v of current [m/s] (positive north)       
      real,allocatable,public :: w(:,:,:)          ! w of current [m/s] (positive down)   
      real,allocatable,public :: temp(:,:,:)       ! Water temp. [Celcius]
      real,allocatable,public :: vdiffus(:,:,:)    ! vertical   diffusivity [m**2/s]              
      real,allocatable,public :: hdiffus(:,:,:)    ! horizontal diffusivity [m**2/s]
      real,allocatable,public :: dslm(:,:)         ! current sea surface elevation over reference [m]
      real,allocatable,public :: zoo(:,:,:)        ! Zooplankton [10^-3 mol N/liter]
      real,allocatable,target,public :: ccdepth(:,:,:)    ! cell center depth water below surface [m]; pointer 
      real,allocatable,target,public :: acc_width(:,:,:)  ! accumulated water above this layer [m] dim=nz+1  

c     --- 2D grids ---
      
      real,allocatable,public     :: wdepth (:,:)      ! current depth at cell-center, including dslm [m]
      integer, allocatable,public :: wetmask(:,:)      ! auxillary to define coastal geometry (1: wet nodes, 0 dry)  
      integer, allocatable,public :: bottom_layer(:,:) ! last wet layer (0 for dry points) nx,ny


c     ----- static coast line auxillary tables/parameters for fast classification (private) -----

      real, parameter      :: facial_tol(2) = 1.0e-4     ! note: lon part should be lat dependent

      real, allocatable    :: vertical_coastlines(:,:)   ! longitude in degrees
      real, allocatable    :: horizontal_coastlines(:,:) ! latitude in degrees
      real, allocatable    :: corners_coastlines(:,:,:)  ! (longitude, latitude) in degrees
      logical, allocatable :: is_coastal(:,:,:)          ! ityp: 1=west, 2=south,3=sw corner  



c     ===================================================
                            contains
c     ===================================================

                      
      subroutine init_regular_lonlat_grid()
c     ------------------------------------------------------
c     Assumes (nx,ny,nz) has been set by client module
c     ------------------------------------------------------
      write(*,*) "init_regular_lonlat_grid: allocate grid arrays: begin" 

      allocate( u(nx,ny,nz)       )   
      allocate( v(nx,ny,nz)       )   
      allocate( w(nx,ny,nz)       )   
      allocate( temp(nx,ny,nz)    )  
      allocate( vdiffus(nx,ny,nz) )                      
      allocate( hdiffus(nx,ny,nz) ) 
      allocate( zoo(nx,ny,nz)     ) 
      allocate( ccdepth(nx,ny,nz) )  
      allocate( acc_width(nx,ny,nz+1)   )

c     --- 2D grids ---

      allocate( dslm(nx,ny)       )       
      allocate( wdepth(nx,ny)     ) 
      allocate( wetmask(nx,ny)    )  ! not static
      allocate( bottom_layer(nx,ny) )
      write(*,*) "init_regular_lonlat_grid:: allocate grid arrays: OK"

c     -----  coast line auxillary tables -----

      allocate( vertical_coastlines(nx+1,ny+1))   
      allocate( horizontal_coastlines(nx+1,ny+1)) 
      allocate( corners_coastlines(nx+1,ny+1,2)) 
      allocate( is_coastal(nx+1,ny+1,3))   
       
      call configure_coastal_hash_tables()

      end subroutine init_regular_lonlat_grid
       


      subroutine close_regular_lonlat_grid()
c     ------------------------------------------------------
c     ------------------------------------------------------
      if (allocated(u))            deallocate( u )
      if (allocated(v))            deallocate( v )
      if (allocated(w))            deallocate( w )
      if (allocated(temp))         deallocate( temp )
      if (allocated(vdiffus))      deallocate( vdiffus )
      if (allocated(hdiffus))      deallocate( hdiffus )
      if (allocated(dslm))         deallocate( dslm )
      if (allocated(zoo))          deallocate( zoo )
      if (allocated(ccdepth))      deallocate( ccdepth )     
      if (allocated(acc_width))    deallocate( acc_width )     
      if (allocated(wdepth))       deallocate( wdepth )
      if (allocated(wetmask))      deallocate( wetmask )
      if (allocated(bottom_layer)) deallocate( bottom_layer )
      
      if (allocated(vertical_coastlines))   
     +    deallocate( vertical_coastlines )
      if (allocated(horizontal_coastlines)) 
     +    deallocate( horizontal_coastlines )
      if (allocated(corners_coastlines))    
     +    deallocate( corners_coastlines )
      if (allocated(is_coastal ))           
     +    deallocate( is_coastal )

      end subroutine close_regular_lonlat_grid


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




      function   get_master_clock()
c     ------------------------------------------ 
      type(clock), pointer :: get_master_clock
c     ------------------------------------------ 
      get_master_clock => master_clock
c     ------------------------------------------ 
      end function 

      subroutine set_master_clock(time)
c     ------------------------------------------ 
      type(clock), intent(in) :: time
c     ------------------------------------------ 
      master_clock = time ! local copy
c     ------------------------------------------ 
      end subroutine 




      subroutine interpolate_cc_3Dgrid_data(xyz,array,deriv,result,
     +                                      status)
c     -------------------------------------------------------------------- 
c     Interpolate on grid of corner centered data 3D array on point xyz.
c     Apply appropriate extrapolations near boundaries or return padvalue
c     For points 0 < z <ccdepth(., ., 1) and ccdepth(., ., ibot) < z < wd
c     the surface/bottom layer value is used, i.e. the value is assumed 
c     constant in upper half of first layer/lower half of bottom layer
c     
c     deriv = 0 gives value, deriv = (1,2,3) gives derivative along (x,y,z)
c     Reserve deriv = 123 for overloaded full gradient ?
c     Currently, only deriv = (0,3) is implemented, until other derivatives 
c     are needed.
c
c     TODO:  multiply metric factor for deriv = (1,2), because
c
c     Return status:
c       status = 0: interior interpolation performed
c       status = 1: horizontal range violation, set result = padval
c       status = 2: vertical extension using boundary values
c                   for a dry point set result = padval
c       status = 3: dry point / rank deficit situation not permitting interpolation
c                   Return result = padval for deriv = 0  
c    
c     tested: deriv=0,3
c     --------------------------------------------------------------------
      real, intent(in)     :: xyz(:),array(:,:,:)
      integer, intent(in)  :: deriv  ! 0=value; (1,2,3) = along (x,y,z)
      real, intent(out)    :: result
      integer, intent(out) :: status
c      
      real, parameter     :: padval = 0. ! later move to argument list
c
      integer             :: ix,iy,iz,i,idum,ibot,ix0,ix1,iy0,iy1
      real                :: z,sx,sy,sz,depth,fbar,fup,flow
      integer             :: cx(4), cy(4), pad_mask(4)
      real                :: fz(4), dz(4), dfz(4), hweight(4)
      logical             :: valid(4)
      real, pointer       :: zgrid(:)
c     --------------------------------------------------------------------
      if (.not.horizontal_range_check(xyz)) then
          result = padval
          if (deriv>0) result = 0 
          status = 1
          return
      endif

      call interpolate_wdepth(xyz,depth,idum)
      if ((xyz(3)<-htol).or.(xyz(3)>depth+htol)) then
         status = 2                ! signal vertical out-of-bound
         result = padval           ! derivative interpolation
         if (deriv>0) result = 0   ! value interpolation
         return
      endif
c
c     define corners in this order
c         col1  (ix0,iy0)
c         col2  (ix0,iy1)
c         col3  (ix1,iy0)  
c         col4  (ix1,iy1)
c
      call get_location_in_mesh(xyz,ix,iy,sx,sy)
c
      ix0 = min(max(ix,  1),nx)  ! cover grid edges
      ix1 = min(max(ix+1,1),nx)  ! cover grid edges
      iy0 = min(max(iy,  1),ny)  ! cover grid edges
      iy1 = min(max(iy+1,1),ny)  ! cover grid edges
      cx(1) = ix0; cy(1) = iy0 
      cx(2) = ix0; cy(2) = iy1 
      cx(3) = ix1; cy(3) = iy0 
      cx(4) = ix1; cy(4) = iy1 
      hweight(1) = (1.0-sx)*(1.0-sy)
      hweight(2) = (1.0-sx)*(  sy  )
      hweight(3) = (  sx  )*(1.0-sy)
      hweight(4) = (  sx  )*(  sy  )

      z          = xyz(3)  ! short hand
      status     = 0       ! assume interpolation possible
      valid      = .false. ! default for dry/out-of-bounds pillars
      fz         = 0.      ! default for dry/out-of-bounds pillars
      dfz        = 0.      ! default for dry/out-of-bounds pillars
c
c     ------ setup 4 interpolation pillars ------
c            project z onto these pillars and interpolate
c
      do i = 1,4
         if (wetmask(cx(i),cy(i))>0) then
            valid(i) = .true.
            ibot     = bottom_layer(cx(i),cy(i)) ! ibot >= 1
            zgrid    => ccdepth(cx(i),cy(i),:)   ! range = 1:nz
            call search_sorted_list(z,zgrid,iz) ! 0<=iz<=nz
            ! ---- handle projection extrapolation 
            !      but do not flag flag vertical extrapolation
            !      since 0<z<depth
            if ((iz<1).or.(iz>=ibot)) then 
               iz = min(max(iz,1),ibot) ! now iz = 1 or ibot
               flow    = array(cx(i),cy(i),iz)
               fup     = flow    ! => dfz = 0
               sz      = 0.5     ! set dummy
c               status  = 2       ! flag vertical extrapolation
               dz(i)   = 1.0     ! avoid numerical problems, assign dummy
            ! ---- interior linear vertical interpolation ----
            else
               flow = array(cx(i),cy(i),iz)
               fup  = array(cx(i),cy(i),iz+1)
               dz(i)= zgrid(iz+1)-zgrid(iz) ! assumed > 0
               sz   = (z-zgrid(iz))/dz(i)
            endif
            fz(i)   = (1.0-sz)*flow + sz*fup
            dfz(i)  = fup-flow
         endif ! wetmask ...
      enddo
c
c     ------ fill in missing data where flagged ------
c
      if (all(valid)) then

         continue

      elseif (.not.any(valid)) then

         result = padval
         if (deriv>0) result = 0 
         status = 3  ! rank deficit exit
         return

      else ! we know at least one valid corner

         fbar = sum(fz)/count(valid)
         do i = 1,4
            if (.not.valid(i)) fz(i) = fbar ! relax to average
         enddo

      endif
c
c     ------ evaluate interpolations/derivatives ------
c
      result = 0
      if (deriv == 0) then     ! evaluate value
         do i = 1,4
            result = result + hweight(i)*fz(i) 
         enddo
      elseif (deriv == 3) then ! evaluate z derivative
         do i = 1,4
            result = result + hweight(i)*dfz(i)/dz(i) 
         enddo
      else
         write(*,*) "interpolate_cc_3Dgrid_data: deriv = ",
     +               deriv, "is not implemented"  
      endif

      end subroutine interpolate_cc_3Dgrid_data




      subroutine interpolate_cc_2Dgrid_data(xy,array,deriv,result,
     +                                      status)
c     --------------------------------------------------------------------
c     Interpolate grid of corner centered data 2D array on point xy
c     with data points at integer valued xy
c
c     deriv = 0 gives value, deriv = (1,2) gives derivative along (x,y)
c     Currently, only deriv = 0, until other derivatives 
c     are needed.
c
c     Return status:
c       status = 0: interior interpolation performed
c       status = 1: horizontal range violation, set result = padval
c     --------------------------------------------------------------------
      real, intent(in)     :: xy(:),array(:,:)
      integer, intent(in)  :: deriv  ! 0=value; deriv = (1,2) gives derivative along (x,y)
      real, intent(out)    :: result
      integer, intent(out) :: status
c      
      integer           :: ix,iy,ix0,ix1,iy0,iy1
      real              :: vc(4),sx,sy
      real, parameter   :: padval = 0. ! later move to argument list
c     --------------------------------------------------------------------
       if (.not.horizontal_range_check(xy)) then
          result = padval
          status = 1
          return
      endif
c
      call get_location_in_mesh(xy,ix,iy,sx,sy)
c      
      ix0 = min(max(ix,    1),nx)  ! cover boundary layers
      ix1 = min(max(ix + 1,1),nx)  ! cover boundary layers
c     
      iy0 = min(max(iy,    1),ny)  ! cover boundary layers
      iy1 = min(max(iy + 1,1),ny)  ! cover boundary layers
c
      vc(1) = array(ix0,iy0)
      vc(2) = array(ix0,iy1)
      vc(3) = array(ix1,iy0)
      vc(4) = array(ix1,iy1)

      if     (deriv == 0) then 
         call interp_2Dbox_data(sx,sy,vc,deriv,result)
      else
         stop "interpolate_cc_2Dgrid_data: unhandled deriv request"
      endif

      status = 0                ! signal interior interpolation performed
c     -------------------------------------------------------------------- 
      end subroutine interpolate_cc_2Dgrid_data



      subroutine interpolate_currents(xyz, r3, status)
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
      integer              :: statu, statv, statw
c     ------------------------------------------
      call interpolate_cc_3Dgrid_data(xyz,u,0,r3(1),statu)
      call interpolate_cc_3Dgrid_data(xyz,v,0,r3(2),statv)
      call interpolate_cc_3Dgrid_data(xyz,w,0,r3(3),statw)
      status  = max(statu, statv, statw) ! in the unextected case taht they differ ...
c     ------------------------------------------ 
      end subroutine interpolate_currents


      subroutine interpolate_turbulence(xyz, r3, status)
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
      integer              :: statv,stath
c     ------------------------------------------ 
      call interpolate_cc_3Dgrid_data(xyz,hdiffus,0,r3(1),stath)
      call interpolate_cc_3Dgrid_data(xyz,vdiffus,0,r3(3),statv)
      r3(2)   = r3(1) ! horizontal isotropy
      status  = max(statv,stath) ! in the unextected case taht they differ ...
c     ------------------------------------------ 
      end subroutine

 

      subroutine interpolate_turbulence_deriv(xyz, r3, status)
c     ------------------------------------------ 
c     Currently do not support horizontal derivatives
c     ------------------------------------------ 
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
      integer              :: statv,stath
      
c     ------------------------------------------ 
      call interpolate_cc_3Dgrid_data(xyz,vdiffus,3,r3(3),statv)
      r3(1:2) = 0. ! Currently do not support horizontal derivatives
      status  = statv 
c     ------------------------------------------    
      end subroutine 






      subroutine interpolate_temp (xyz, r, status) 
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
c     ------------------------------------------ 
      call interpolate_cc_3Dgrid_data(xyz,temp,0,r,status)
c     ------------------------------------------ 
      end subroutine 



      subroutine interpolate_zooplankton (xyz, r, status) 
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r(:)
      integer, intent(out) :: status
c     ------------------------------------------ 
      call interpolate_cc_3Dgrid_data(xyz,zoo,0,r(1),status)
c     ------------------------------------------ 
      end subroutine 



      subroutine interpolate_wdepth(xy, r, status) 
c     ------------------------------------------ 
c     Multiply wetmask to ensure piecewise linear coastlines
c     defined by wdepth=0 
c     ------------------------------------------ 
      real, intent(in)     :: xy(:) 
      real, intent(out)    :: r
      integer, intent(out) :: status
      integer           :: ixc,iyc
c     ------------------------------------------ 
      call interpolate_cc_2Dgrid_data(xy,wdepth,0,r,status)
      call get_horiz_ncc(xy,ixc,iyc)
      r = r*wetmask(ixc,iyc)
c     ------------------------------------------ 
      end subroutine 





      LOGICAL function is_wet(xyz)
c     ------------------------------------------ 
c     Probe wdepth
c     return is_wet = .true. at horizontal range violation
c     accept points at sea surface/bottom as wet (htol)
c     ------------------------------------------ 
      real, intent(in) :: xyz(:) 
      real             :: wd
      integer          :: status
c     ------------------------------------------ 
      call interpolate_wdepth(xyz, wd, status) 
      if (status==1) then
         is_wet = .true.
         return
      elseif (status==0) then
         if ((xyz(3)>-htol).and.(xyz(3)<wd+htol)) then
            is_wet = .true.
         else
            is_wet = .false.
         endif
      else
         stop "is_wet: unhandled return condition"
      endif
c     ------------------------------------------
      end function



      LOGICAL function is_land(xy)
c     ------------------------------------------ 
c     return is_land = .false. at horizontal range violation
c     ------------------------------------------ 
      real, intent(in) :: xy(:)
      real             :: wdepth
      integer          :: status
      integer          :: ixc,iyc
c     ------------------------------------------ 
c     avoid probing wetmask at range violation

      if (horizontal_range_check(xy)) then ! range OK
         call get_horiz_ncc(xy,ixc,iyc)
         is_land = ((wetmask(ixc,iyc)==0).or.at_coast_line(xy))
      else
         is_land = .false.
      endif    
c     ------------------------------------------
      end function



      LOGICAL function horizontal_range_check(xy)
c     ------------------------------------------
c     Return .true. if cell association of xy
c     0 < (ix,iy) < (nx,ny) is obeyed
c     upper limit (i<n): require data coverage on both
c     sides of an interior point
c     ------------------------------------------ 
      real, intent(in) :: xy(:)
      integer          :: ix,iy
      real             :: sx,sy
c     ------------------------------------------ 
      call get_location_in_mesh(xy,ix,iy,sx,sy)
      horizontal_range_check = (0<ix).and.(ix<nx)
     +                    .and.(0<iy).and.(iy<ny)
c     ------------------------------------------ 
      end function





      subroutine coast_line_intersection(geo1, geo2,cross_coast,georef,  ! -> regular_lonlat_grid.f
     +                                   geohit) 
c     ------------------------------------------ 
c     Trajectory tracking variant 
c
c     (sandbox: oceanography_providers/sandbox/coast_line_intersection_hash_tables.f)
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
c     ---------------------------------------------------- 
      real, intent(in)     :: geo1(:),geo2(:)
      logical, intent(out) :: cross_coast
      real, intent(out)    :: georef(:), geohit(:)
c   
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
         cross_coast = check_coastal_osculation(geo2,georef,geohit)
         if ((verbose>0).and.cross_coast) 
     +       write(*,*) "coastal osculation of geo2 handled"
         return
      endif
      if (verbose>0) write(*,*) "trajectory analysis: begin"

c     Now we know that we are leaving the cell. But where? 
c     Follow sequence moving from one cell to the next until we either
c     reach the end of the vector (s=1), or encounter land

      cross_coast = .false. ! default unless detected below
      loopflag    = .true.
      istep       = 1
      maxsteps    = max(1,abs(ix-ix2))*max(1,abs(iy-iy2)) ! safety plug

      do while(loopflag)
         call assess_cell_crossing(geo1,geo2,ix,iy,s,xface)
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
c     0 < s < 1 is the coordinate along the vector (geo2-geo1) 
c     corresponding to the point geohit
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

      subroutine assess_cell_crossing (geo1,geo2,ix,iy,s,exitface) ! internal to coast_line_intersection
c     --------------------------------------------------------------
c     Determine where the direct line from geo1 through geo2 
c     leaves the cell with center (ix,iy)
c     Return 0 < s  which is the coordinate along the vector (geo2-geo1) 
c     where the vector exits the cell. Notice s may exceed 1 
c     (which is the case, if cell (ix,iy) contains geo2)
c     Also identify the exit face as one of "N","S","E","W" 
c     so the proper reflection can be calculated
c
c     Modified from assess_cell_leaving to avoid decision conflicts ASC/18Jan2011
c     --------------------------------------------------------------
      real, intent(in)         :: geo1(:),geo2(:)
      integer, intent(in)      :: ix,iy
      real, intent(out)        :: s
      character*1, intent(out) :: exitface
c
      logical                  :: yes
      real                     :: stest,tdum,dgeo(2)
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
      dgeo     = geo2(1:2)-geo1(1:2)
c
c     First the east-west direction
c
      if(dgeo(1)>0) then    !we're moving to the east

         call cross_2Dlines(geo1,geo2,c10,c11,stest,tdum,yes)
         if (yes.and.(stest<s)) then ! rely on right-to-left evaluation
            s         = stest
            exitface  = "E"   
         endif
      else   !we're moving to the west
         call cross_2Dlines(geo1,geo2,c00,c01,stest,tdum,yes)
         if (yes.and.(stest<s)) then ! rely on right-to-left evaluation
            s         = stest
            exitface  = "W"   
         endif
      endif
c
c     Then the north-south direction
c
      if(dgeo(2) >0) then !we're moving to the North
         call cross_2Dlines(geo1,geo2,c01,c11,stest,tdum,yes)
         if (yes.and.(stest<s)) then ! rely on right-to-left evaluation
            s         = stest
            exitface  = "N"   
         endif
      else  !we're moving to the south
         call cross_2Dlines(geo1,geo2,c00,c10,stest,tdum,yes)
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
         write(*,*) "geo1    =",geo1
         write(*,*) "geo2    =",geo2
         write(*,*) "(ix,iy) =",ix,iy
         stop       
      endif
      end subroutine assess_cell_crossing     ! local subroutine

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


     

      subroutine get_location_in_mesh(xy,ix,iy,sx,sy)
c     ------------------------------------------------------
c     From the lon/lat vector xy (which may include the z component)
c     determine the cell association (ix,iy) where the node of the cell is the
c     lower left corner point, so that [ix:ix+1, iy:iy+1[ maps to cell (ix,iy)
c     i.e. a corner centered convention.
c     0 < (sx,sy) < 1 is the intra cell cooridnates. 
c     Grid origo (ix,iy) = (x,y) = (1,1) is associated with (lambda1,phi1)
c     Do not flag grid range violations, so that 
c     0 < (ix,iy) <= (nx,ny) is not enforced 
c     (horizontal_range_check uses this function to check grid range violations)
c     Previously, this subroutine was called get_horiz_grid_coordinates
c     ------------------------------------------------------
      real, intent(in)     :: xy(:)
      integer, intent(out) :: ix,iy
      real, intent(out)    :: sx,sy
      real                 :: dx1,dy1
c     ------------------------------------------------------
      dx1 = (xy(1)-lambda1)/dlambda
      dy1 = (xy(2)-phi1)   /dphi
      ix  = 1 + int(dx1)    ! truncate decimals to get lon cell
      iy  = 1 + int(dy1)    ! truncate decimals to get lat cell
      sx  = dx1 - int(dx1)  ! intra cell coordinate 0<sx<1
      sy  = dy1 - int(dy1)  ! intra cell coordinate 0<sy<1
c     ------------------------------------------------------
      end subroutine get_location_in_mesh


      subroutine get_ncc_coordinates(xyz,x,y,z)
c     ------------------------------------------ 
c     Get continuous node-centered grid coordinates with grid points 
c     on integer values of (x,y,z) from xyz = (lon,lat,depth)
c     water surface is at z = 0.5, sea bed at z = bottum_layer+0.5
c     It is not checked that z is above the sea bed. The inter grid
c     range is 0.0 <= z <= nz+0.5.
c     If vertical range is exceeded the first/last layer, respectively,
c     is used to extrapolate a vertical grid coordinate 
c     (no extrapolation is flagged) 
c     ------------------------------------------ 
      real, intent(in)  :: xyz(:)
      real, intent(out) :: x,y,z

      integer           :: ixc,iyc,izc
      real, pointer     :: this_col(:)
      real              :: layerw
c     ------------------------------------------ 
      x   = 1.0 + (xyz(1)-lambda1)/dlambda
      y   = 1.0 + (xyz(2)-phi1)   /dphi
      ixc = nint(x)     
      iyc = nint(y)
c     --- locate vertical cell ---
      this_col => acc_width(ixc, iyc, :) ! range = 1:nz+1
      call search_sorted_list(xyz(3), this_col, izc) ! 0<=izc<=nz+1
      izc = min(max(izc,1),nz) ! capture vertical range excess
      layerw = this_col(izc+1) -  this_col(izc) 
      z      = (izc - 0.5) + (xyz(3)-this_col(izc))/layerw
c     ------------------------------------------ 
      end subroutine get_ncc_coordinates


      subroutine get_horiz_ncc(xyz,ixc,iyc) !,xc,yc)
c     ------------------------------------------ 
c     Resolve the node-centered grid cell containing xyz
c     (ixc,iyc) is grid indices of the closest node
c     (xc,yc)   is (lon, lat) coordinates of the closest node
c     no range check
c     ------------------------------------------ 
      real, intent(in)     :: xyz(:)
      integer, intent(out) :: ixc,iyc
c      real, intent(out)    :: xc,yc
c     ------------------------------------------ 
      ixc = 1 + nint((xyz(1)-lambda1)/dlambda)     
      iyc = 1 + nint((xyz(2)-phi1)   /dphi)   
c      xc  = lambda1 + (ixc-1)*dlambda  ! unused
c      yc  = phi1    + (iyc-1)*dphi     ! unused
      end subroutine get_horiz_ncc



      subroutine get_horiz_ncc_corners(ixc,iyc,c00,c01,c10,c11)
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




      end module
