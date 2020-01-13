      module regular_lonlat_data   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This module provide a leight class for 2D/3D regular lonlat data to allow mulitple grids
C     
c     Offer basic interpolation + data-related queries
c     Inlined subroutines from horizontal_grid_transformations_lonlatgrid (with grid argument added)
c     Do not reexport any imports from regular_lonlat_grid
c
c     Query return codes from regular_lonlat_grid
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use regular_lonlat_grid
      use array_tools
      
      implicit none
      private           ! default scope

      type lonlat_data2D
        type(lonlat_grid), pointer :: grid            ! instance does not own grid; allow lonlat_data*D instances to share a common lonlat_grid
        real, pointer              :: buf(:,:)        ! instance owns its data buffer
        integer(kind=selected_int_kind(16)) :: revtag ! dynamic revision number of data array; user must update it
      end type
      public :: lonlat_data2D 
      
      type lonlat_data3D
        type(lonlat_grid), pointer :: grid            ! instance does not own grid; allow lonlat_data*D instances to share a common lonlat_grid
        real, pointer              :: buf(:,:,:)      ! instance owns its data buffer
        integer(kind=selected_int_kind(16)) :: revtag ! dynamic revision number of data array; user must update it
      end type
      public :: lonlat_data3D

      public :: init_lonlat_data
      public :: finalize_lonlat_data
      public :: interpolate_cc_data
      public :: update_revtag_data   
      public :: err_msg_reglonlat ! reexport from regular_lonlat_grid

      interface init_lonlat_data
         module procedure init_lonlat_data2D
         module procedure init_lonlat_data3D 
      end interface   

      interface finalize_lonlat_data
         module procedure finalize_lonlat_data2D
         module procedure finalize_lonlat_data3D
      end interface

      interface interpolate_cc_data 
         module procedure interpolate_cc_data2D
         module procedure interpolate_cc_data3D
      end interface  

      interface update_revtag_data
         module procedure update_revtag_data2D
         module procedure update_revtag_data3D
      end interface  
      
c     -------------------- module data --------------------  
      
      integer, parameter :: verbose = 0    ! debugging output control
      real, parameter    :: htol   = 1.e-6 ! tolerance for surface/bottom
      real, parameter    :: padval = 0     ! numeric value applied where undefined

      
c     ===================================================
      contains
c     ===================================================

      subroutine init_lonlat_data2D(this, grid)
      type(lonlat_data2D), intent(inout)    :: this
      type(lonlat_grid), intent(in), target :: grid
      this%grid => grid         ! instance does not own grid instance
      allocate ( this%buf(grid%nx, grid%ny) )
      this%revtag = 0           ! signal unset data buffer
      end subroutine init_lonlat_data2D
      
      subroutine finalize_lonlat_data2D(this)
      type(lonlat_data2D), intent(inout)    :: this
      nullify(this%grid )  ! instance does not own grid instance
      if (associated(this%buf))  deallocate ( this%buf )
      this%revtag = -1           ! signal void data buffer
      end subroutine finalize_lonlat_data2D                
      
      subroutine init_lonlat_data3D(this, grid)
      type(lonlat_data3D), intent(inout)    :: this
      type(lonlat_grid), intent(in), target :: grid
      this%grid => grid         ! instance does not own grid instance
      if (.not.grid%is3D) then
         write(*,*) "init_lonlat_data3D: grid is only 2D"
         stop 822
      endif
      allocate ( this%buf(grid%nx, grid%ny, grid%nz) )
      this%revtag = 0           ! signal unset data buffer
      end subroutine init_lonlat_data3D
      
      subroutine finalize_lonlat_data3D(this)
      type(lonlat_data3D), intent(inout)    :: this
      nullify(this%grid )  ! instance does not own grid instance
      if (associated(this%buf))  deallocate ( this%buf )
      this%revtag = -1 
      end subroutine finalize_lonlat_data3D



      

      subroutine interpolate_cc_data2D(this, geo, deriv, padval,
     +                                      result, status)
c     --------------------------------------------------------------------
c     Interpolate grid of cell centered data 2D array on point geo
c     with data points at integer valued ngeo
c
c     deriv = 0 gives value, deriv = (1,2) gives derivative along (x,y)
c     Currently, only deriv = 0, until other derivatives 
c     are needed.
c     Do not check vertival position for horizontal interpolations
c
c     Return status: error codes from regular_lonlat_grid
c     Return result = padval for deriv = 0
c      
c     ASCH 12 Dec 2019: transcribed from mesh_grid::interpolate_cc_2Dgrid_data (array = this%buf)    
c     --------------------------------------------------------------------
      type(lonlat_data2D), intent(in) :: this
      real, intent(in)                :: geo(:)
      integer, intent(in)             :: deriv  ! 0=value; deriv = (1,2) gives derivative along (x,y)
      real, intent(in)                :: padval
      real, intent(out)               :: result
      integer, intent(out)            :: status
c      
      integer           :: i,ix,iy,ix0,ix1,iy0,iy1
      real              :: vc(4),sx,sy,vcbar
      integer           :: cx(4), cy(4)
      logical           :: valid(4)
c     --------------------------------------------------------------------
c     ... set defaults ...
      if (deriv>0) then
         result = 0             ! derivative interpolation
      else
         result = padval        ! value interpolation
      endif
      
      if (.not.horizontal_range_check_for_grid(this%grid, geo)) then
          status = err_reglonlat_horiz_violation
          return
      endif
c     
      if (is_land_for_grid(this%grid, geo)) then
          status = err_reglonlat_dry_point
          return
      endif
c
      if (this%revtag .le. 0) then
          status = err_reglonlat_runtimeerr  ! uninitialized data array
          return
       endif
c
      call get_surrounding_box(this%grid, geo,ix,iy,sx,sy)

      ix0 = min(max(ix,    1), this%grid%nx)  ! cover boundary layers
      ix1 = min(max(ix + 1,1), this%grid%nx)  ! cover boundary layers
      iy0 = min(max(iy,    1), this%grid%ny)  ! cover boundary layers
      iy1 = min(max(iy + 1,1), this%grid%ny)  ! cover boundary layers
c
c     define corners in this order
c         col1  (ix0,iy0)
c         col2  (ix0,iy1)
c         col3  (ix1,iy0)  
c         col4  (ix1,iy1)
c
      cx(1) = ix0; cy(1) = iy0 
      cx(2) = ix0; cy(2) = iy1 
      cx(3) = ix1; cy(3) = iy0 
      cx(4) = ix1; cy(4) = iy1 

    
      status     = err_reglonlat_success ! assume interpolation possible
      valid      = .false.               ! default for dry/out-of-bounds pillars
      vc         = 0.                    ! default for dry/out-of-bounds pillars (NOT padval!)
c
c     ------ setup 4 interpolation corners ------
c
      do i = 1,4
         if (this%grid%wetmask(cx(i),cy(i))>0) then
            valid(i) = .true.
            vc(i)    = this%buf(cx(i),cy(i))
         endif ! wetmask ...
      enddo
c
c     ------ fill in missing data where flagged ------
c
      if (all(valid)) then

         continue

      elseif (.not.any(valid)) then

         status = err_reglonlat_dry_point  ! rank deficit exit
         return

      else ! we know at least one valid corner

         vcbar     = sum(vc)/count(valid)    ! invalid => vc=0    
         do i = 1,4
            if (.not.valid(i)) vc(i) = vcbar
         enddo

      endif
c
c     ------ delegate to horizontal interpolation ------
c 
      if     (deriv == 0) then 
         call interp_2Dbox_data(sx,sy,vc,0,result) ! interpolate value from vc
      else
         stop "interpolate_cc_2Dgrid_data: unhandled deriv request"
      endif

c     -------------------------------------------------------------------- 
      end subroutine interpolate_cc_data2D


      
      
      subroutine interpolate_cc_data3D(this, geo, deriv, padval,
     +                                      result, status)
c     --------------------------------------------------------------------
c     
c     Interpolate on grid of cell centered data 3D array on point geo.
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
c     Return status: return codes from regular_lonlat_grid
c       err_reglonlat_success:            interior interpolation performed
c       err_reglonlat_horiz_violation:    set result = padval
c       err_reglonlat_vertical_violation: set result = padval
c       err_reglonlat_dry_point:          set result = padval for deriv = 0, 0 for deriv >0
c      
c     tested: deriv=0,3
c
c     ASCH 12 Dec 2019: transcribed from mesh_grid::interpolate_cc_3Dgrid_data (aŕray = this%buf)
c     --------------------------------------------------------------------
      type(lonlat_data3D), intent(in) :: this
      real, intent(in)                :: geo(:)
      integer, intent(in)             :: deriv  ! 0=value; (1,2,3) = along (x,y,z)
      real, intent(in)                :: padval
      real, intent(out)               :: result
      integer, intent(out)            :: status
c      
      integer             :: ix,iy,iz,i,idum,ibot,ix0,ix1,iy0,iy1
      real                :: z,sx,sy,sz,depth,fbar,fup,flow,dfzdzbar
      integer             :: cx(4), cy(4)
      real                :: fz(4), dz, dfzdz(4)
      logical             :: valid(4)
      real, pointer       :: zgrid(:)
c     --------------------------------------------------------------------
c     ... set defaults ...
      if (deriv>0) then
         result = 0             ! derivative interpolation
      else
         result = padval        ! value interpolation
      endif
      
      if (.not.horizontal_range_check_for_grid(this%grid, geo)) then
          status = err_reglonlat_horiz_violation
          return
      endif
  
      if (is_land_for_grid(this%grid, geo)) then
          status = err_reglonlat_dry_point
          return
      endif

      if (this%revtag .le. 0) then
          status = err_reglonlat_runtimeerr   ! uninitialized data array
          return
       endif
       
c     --- accept points at bottum and on water surface --- 
c     --- within a numerical tolerance htol            ---

      call interpolate_wdepth_for_grid(this%grid,geo,depth,idum)
      if ((geo(3)<-htol).or.(geo(3)>depth+htol)) then
         status = err_reglonlat_vertical_violation  ! signal vertical out-of-bound
         return
      endif
c
c     define corners in this order
c         col1  (ix0,iy0)
c         col2  (ix0,iy1)
c         col3  (ix1,iy0)  
c         col4  (ix1,iy1)
c
      call get_surrounding_box(this%grid,geo,ix,iy,sx,sy)
    
      ix0 = min(max(ix,  1), this%grid%nx)  ! cover grid edges
      ix1 = min(max(ix+1,1), this%grid%nx)  ! cover grid edges
      iy0 = min(max(iy,  1), this%grid%ny)  ! cover grid edges
      iy1 = min(max(iy+1,1), this%grid%ny)  ! cover grid edges
      cx(1) = ix0; cy(1) = iy0 
      cx(2) = ix0; cy(2) = iy1 
      cx(3) = ix1; cy(3) = iy0 
      cx(4) = ix1; cy(4) = iy1 

      z          = geo(3)  ! short hand
      status     = err_reglonlat_success ! assume interpolation possible
      valid      = .false. ! default for dry/out-of-bounds pillars
      fz         = 0.      ! default for dry/out-of-bounds pillars (NOT padval!)
      dfzdz      = 0.      ! default for dry/out-of-bounds pillars 
c
c     ------ setup 4 interpolation pillars ------
c            project z onto these pillars and interpolate
c
      do i = 1,4
         if (this%grid%wetmask(cx(i),cy(i))>0) then
            valid(i) = .true.
            ibot     =  this%grid%bottom_layer(cx(i),cy(i)) ! ibot >= 1
            zgrid    => this%grid%ccdepth(cx(i),cy(i),1:ibot)   ! range = 1:ibot<=nz
            call search_sorted_list(z,zgrid,iz) ! 0<=iz<=nz
            ! ---- handle projection extrapolation 
            !      but do not flag flag vertical extrapolation
            !      since 0<z<depth
            if ((iz<1).or.(iz>=ibot)) then 
               iz = min(max(iz,1),ibot) ! now iz = 1 or ibot
               flow    = this%buf(cx(i),cy(i),iz)
               fup     = flow    ! => dfzdz = 0
               sz      = 0.5     ! set dummy
c               status  = err_reglonlat_vertical_violation       ! flag vertical extrapolation
               dz      = 1.0     ! avoid numerical problems, assign dummy
            ! ---- interior linear vertical interpolation ----
            else
               flow = this%buf(cx(i),cy(i),iz)
               fup  = this%buf(cx(i),cy(i),iz+1)
               dz   = zgrid(iz+1)-zgrid(iz) ! assumed > 0
               sz   = (z-zgrid(iz))/dz
            endif
            fz(i)    = (1.0-sz)*flow + sz*fup
            dfzdz(i) = (fup-flow)/dz
         endif ! wetmask ...
      enddo
c
c     ------ fill in missing data where flagged ------
c
      if (all(valid)) then

         continue

      elseif (.not.any(valid)) then

         status = err_reglonlat_dry_point  ! rank deficit exit
         return

      else ! we know at least one valid corner

         fbar     = sum(fz)/count(valid)    ! invalid => fz=0
         dfzdzbar = sum(dfzdz)/count(valid) 
         do i = 1,4
            if (.not.valid(i)) then
               fz(i)    = fbar     ! relax to average
               dfzdz(i) = dfzdzbar ! relax to average
            endif
         enddo

      endif
c
c     ------ delegate to horizontal interpolation ------
c   
      if     (deriv == 0) then     ! interpolate value from fz
         call interp_2Dbox_data(sx,sy,fz,0,result)
      elseif (deriv == 3) then     ! interpolate value from dfzdz
         call interp_2Dbox_data(sx,sy,dfzdz,0,result)
      else
         write(*,*) "interpolate_cc_data3D: deriv = ",
     +               deriv, "is not implemented"  
      endif

      end subroutine interpolate_cc_data3D


      
      subroutine update_revtag_data2D(this, revno)
c     --------------------------------------------------------
c     update/increment revtag
c     if optional revno is absent, increment current revtag
c     --------------------------------------------------------      
      type(lonlat_data2D), intent(out) :: this
      integer(kind=selected_int_kind(16)), optional :: revno
c     --------------------------------------------------------      
      if (.not.associated(this%buf)) then
         write(*,*) "update_revtag_data2D: buffer not associated"
         stop 733
      endif
c
      if (present(revno)) then
         this%revtag = revno                          ! apply user provided revtag
      else                                            ! default: just increment counter
         if (this%revtag == huge(this%revtag)) then   ! HUGE(X) The largest positive number
            this%revtag = 1                           ! restart counter for increments
         else
            this%revtag = this%revtag + 1             ! 
         endif
      endif   ! present(revno)
      end subroutine update_revtag_data2D


      
      subroutine update_revtag_data3D(this, revno)
c     --------------------------------------------------------
c     update/increment revtag
c     if optional revno is absent, increment current revtag
c     --------------------------------------------------------      
      type(lonlat_data3D), intent(out) :: this
      integer(kind=selected_int_kind(16)), optional :: revno
c     --------------------------------------------------------      
      if (.not.associated(this%buf)) then
         write(*,*) "update_revtag_data3D: buffer not associated"
         stop 733
      endif
c
      if (present(revno)) then
         this%revtag = revno                          ! apply user provided revtag
      else                                            ! default: just increment counter
         if (this%revtag == huge(this%revtag)) then   ! HUGE(X) The largest positive number
            this%revtag = 1                           ! restart counter for increments
         else
            this%revtag = this%revtag + 1             ! 
         endif
      endif   ! present(revno)
      end subroutine update_revtag_data3D


      
      end module
      
