      module mesh_grid
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This module represents data on a 3D array where the vertical
c     mesh is perpendicular to the horizontal and the
c     horizontal grid are a mesh (without further restrictions)
c     This includes sigma and z-grid with open/closed surface
c     
c     This module hosts data arrays, add vertical issues to horizontaal to get 3D
c     and proivide 2D/3D interpolation
c
c     * provide interpolate_X
c     * provide time services
c     * provide is_wet
c     * host nz
c     * host 3D data arrays
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use horizontal_representation   ! injects nx,ny
      use time_tools                  ! import clock type
      use constants
      use array_tools                 ! search_sorted

      implicit none
      private           ! default scope

      public :: init_mesh_grid   ! init  this module
      public :: close_mesh_grid  ! close this module


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
      public :: is_land                 ! reexport from horizontal_representation
      public :: horizontal_range_check  ! reexport from horizontal_grid_transformations
      public :: coast_line_intersection ! reexport from horizontal_representation
      
c     ---------- other exports ----------

      public :: get_grid_coordinates     ! formerly named get_ncc_coordinates

c     -------------------- module data --------------------  
      
      type(clock), target,public :: master_clock

      integer, parameter :: verbose = 0  ! debugging output control
      real, parameter    :: htol = 1.e-6 ! tolerance for surface/bottom

c
c     ------ grid dimensions:   ------
c     
      public         :: nx,ny ! reexport from horizontal_grid_transformations 
      integer,public :: nz    ! basic vertical grid dimension 
       

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
      integer, allocatable,public :: bottom_layer(:,:) ! last wet layer (0 for dry points) nx,ny



c     ===================================================
                            contains
c     ===================================================

                      
      subroutine init_mesh_grid()
c     ------------------------------------------------------
c     Assumes (nx,ny,nz) has been set by client module
c     ------------------------------------------------------
      write(*,*) "init_mesh_grid: allocate grid arrays: begin" 

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
      allocate( bottom_layer(nx,ny) )
      write(*,*) "init_mesh_grid: allocate grid arrays: OK"

      call init_horizontal_representation()

      end subroutine init_mesh_grid
       



      subroutine close_mesh_grid()
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
      if (allocated(bottom_layer)) deallocate( bottom_layer )
      
      call close_horizontal_representation()

      end subroutine close_mesh_grid



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




      subroutine interpolate_cc_3Dgrid_data(geo,array,deriv,result,
     +                                      status)
c     -------------------------------------------------------------------- 
c     Interpolate on grid of corner centered data 3D array on point geo.
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
      real, intent(in)     :: geo(:),array(:,:,:)
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
      if (.not.horizontal_range_check(geo)) then
          result = padval
          if (deriv>0) result = 0 
          status = 1
          return
      endif

      call interpolate_wdepth(geo,depth,idum)
      if ((geo(3)<-htol).or.(geo(3)>depth+htol)) then
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
      call get_location_in_mesh(geo,ix,iy,sx,sy)
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

      z          = geo(3)  ! short hand
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




      subroutine interpolate_cc_2Dgrid_data(geo,array,deriv,result,
     +                                      status)
c     --------------------------------------------------------------------
c     Interpolate grid of corner centered data 2D array on point geo
c     with data points at integer valued geo
c
c     deriv = 0 gives value, deriv = (1,2) gives derivative along (x,y)
c     Currently, only deriv = 0, until other derivatives 
c     are needed.
c
c     Return status:
c       status = 0: interior interpolation performed
c       status = 1: horizontal range violation, set result = padval
c     --------------------------------------------------------------------
      real, intent(in)     :: geo(:),array(:,:)
      integer, intent(in)  :: deriv  ! 0=value; deriv = (1,2) gives derivative along (x,y)
      real, intent(out)    :: result
      integer, intent(out) :: status
c      
      integer           :: ix,iy,ix0,ix1,iy0,iy1
      real              :: vc(4),sx,sy
      real, parameter   :: padval = 0. ! later move to argument list
c     --------------------------------------------------------------------
       if (.not.horizontal_range_check(geo)) then
          result = padval
          status = 1
          return
      endif
c
      call get_location_in_mesh(geo,ix,iy,sx,sy)
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



      subroutine interpolate_currents(geo, r3, status)
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
      integer              :: statu, statv, statw
c     ------------------------------------------
      call interpolate_cc_3Dgrid_data(geo,u,0,r3(1),statu)
      call interpolate_cc_3Dgrid_data(geo,v,0,r3(2),statv)
      call interpolate_cc_3Dgrid_data(geo,w,0,r3(3),statw)
      status  = max(statu, statv, statw) ! in the unextected case taht they differ ...
c     ------------------------------------------ 
      end subroutine interpolate_currents


      subroutine interpolate_turbulence(geo, r3, status)
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
      integer              :: statv,stath
c     ------------------------------------------ 
      call interpolate_cc_3Dgrid_data(geo,hdiffus,0,r3(1),stath)
      call interpolate_cc_3Dgrid_data(geo,vdiffus,0,r3(3),statv)
      r3(2)   = r3(1) ! horizontal isotropy
      status  = max(statv,stath) ! in the unextected case taht they differ ...
c     ------------------------------------------ 
      end subroutine

 

      subroutine interpolate_turbulence_deriv(geo, r3, status)
c     ------------------------------------------ 
c     Currently do not support horizontal derivatives
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
      integer              :: statv,stath
      
c     ------------------------------------------ 
      call interpolate_cc_3Dgrid_data(geo,vdiffus,3,r3(3),statv)
      r3(1:2) = 0. ! Currently do not support horizontal derivatives
      status  = statv 
c     ------------------------------------------    
      end subroutine 






      subroutine interpolate_temp (geo, r, status) 
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
c     ------------------------------------------ 
      call interpolate_cc_3Dgrid_data(geo,temp,0,r,status)
c     ------------------------------------------ 
      end subroutine 



      subroutine interpolate_zooplankton (geo, r, status) 
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r(:)
      integer, intent(out) :: status
c     ------------------------------------------ 
      call interpolate_cc_3Dgrid_data(geo,zoo,0,r(1),status)
c     ------------------------------------------ 
      end subroutine 



      subroutine interpolate_wdepth(geo, r, status) 
c     ------------------------------------------ 
c     Multiply by is_land to ensure piecewise linear coastlines
c     defined by wdepth=0 
c     ------------------------------------------ 
      real, intent(in)     :: geo(:) 
      real, intent(out)    :: r
      integer, intent(out) :: status    
c     ------------------------------------------ 
      call interpolate_cc_2Dgrid_data(geo,wdepth,0,r,status)
      if (is_land(geo)) r=0.0
c     ------------------------------------------ 
      end subroutine 



      LOGICAL function is_wet(geo)
c     ------------------------------------------ 
c     Probe wdepth
c     return is_wet = .true. at horizontal range violation
c     accept points at sea surface/bottom as wet (htol)
c     ------------------------------------------ 
      real, intent(in) :: geo(:) 
      real             :: wd
      integer          :: status
c     ------------------------------------------ 
      call interpolate_wdepth(geo, wd, status) 
      if (status==1) then
         is_wet = .true.
         return
      elseif (status==0) then
         if ((geo(3)>-htol).and.(geo(3)<wd+htol)) then
            is_wet = .true.
         else
            is_wet = .false.
         endif
      else
         stop "is_wet: unhandled return condition"
      endif
c     ------------------------------------------
      end function



      subroutine get_grid_coordinates(geo,x,y,z)  ! formerly named get_ncc_coordinates
c     ------------------------------------------ 
c     Get continuous node-centered grid coordinates with grid points 
c     on integer values of (x,y,z) from geo = (lon,lat,depth)
c     water surface is at z = 0.5, sea bed at z = bottum_layer+0.5
c     It is not checked that z is above the sea bed. The inter grid
c     range is 0.0 <= z <= nz+0.5.
c     If vertical range is exceeded the first/last layer, respectively,
c     is used to extrapolate a vertical grid coordinate 
c     (no extrapolation is flagged) 
c     ------------------------------------------ 
      real, intent(in)  :: geo(:)
      real, intent(out) :: x,y,z

      integer           :: ixc,iyc,izc
      real, pointer     :: this_col(:)
      real              :: layerw
c     ------------------------------------------ 
      call get_horiz_grid_coordinates(geo,x,y)
      call get_horiz_ncc(geo,ixc,iyc) 
c     --- locate vertical cell ---
      this_col => acc_width(ixc, iyc, :) ! range = 1:nz+1
      call search_sorted_list(geo(3), this_col, izc) ! 0<=izc<=nz+1
      izc = min(max(izc,1),nz) ! capture vertical range excess
      layerw = this_col(izc+1) -  this_col(izc) 
      z      = (izc - 0.5) + (geo(3)-this_col(izc))/layerw
c     ------------------------------------------ 
      end subroutine get_grid_coordinates


      end module
