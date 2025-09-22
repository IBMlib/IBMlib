      module mesh_grid_env
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Study variant of interpolate_wdepth that converges to wdepth=0
c     at numerical shorelines (at faces to dry cells, when depth is at cell centers)
c     in a minimal 2D subset of mesh_grid context; elements from horizontal_representation,
c     from horizontal_grid_transformations are included here
c
c     ifort -I/usr/local_intel/include coastal_depth_interpolation.f -L/usr/local_intel/lib  -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz  -o test_coast_intp
c  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      private           ! default scope

      public :: init_mesh_grid   ! init  this module   
      public :: interpolate_wdepth

c      
      public :: is_land                 ! reexport from horizontal_representation
      public :: horizontal_range_check  ! reexport from horizontal_grid_transformations

      
c     -------------------- module data --------------------  

      real,public :: lambda1         ! inlined from horizontal_grid_transformations_lonlatgrid.f
      real,public :: dlambda         ! inlined from horizontal_grid_transformations_lonlatgrid.f
      real,public :: phi1            ! inlined from horizontal_grid_transformations_lonlatgrid.f
      real,public :: dphi            ! inlined from horizontal_grid_transformations_lonlatgrid.f
      real,public :: lambda_max      ! aux
      real,public :: phi_max         ! aux
c
c     ------ grid dimensions:   ------
c     
      integer,public      :: nx,ny ! reexport from horizontal_grid_transformations 
    
c     --- 2D grids ---
      
      real,allocatable    :: wdepth(:,:)      ! current depth at cell-center, including dslm [m]. wdepth = 0 at dry points
      integer,allocatable :: wetmask(:,:)           
      real,parameter      :: padval_wdepth   = 0
      logical, parameter  :: dirichlet_coastal_wdepth_BC = .false. ! if true, wdepth=0 at numerical shoreline NEW
      
c     ===================================================
                            contains
c     ===================================================

      subroutine init_mesh_grid()
c     ------------------------------------------------------
c     setup testing environment 
c     ------------------------------------------------------
      integer, parameter :: nlan = 8
      integer            :: ixlist(nlan) = (/1,3,4,3,4,5,6,6/)
      integer            :: iylist(nlan) = (/1,3,3,4,4,4,5,2/)
      integer            :: i, ix, iy, nse
      integer            :: seed(3) = (/1,2,3/)
c     ------------------------------------------------------      
      nx = 8
      ny = 7
      allocate( wdepth(nx,ny)  )
      allocate( wetmask(nx,ny) )
      nse = size(seed) 
      CALL RANDOM_SEED( )  ! Processor reinitializes the seed
      CALL RANDOM_SEED(SIZE = nse)  
      CALL RANDOM_SEED(PUT=seed)
      CALL RANDOM_NUMBER(wdepth) ! 0 < wdepth < 1
      wdepth = wdepth*20         ! 0 < wdepth < 20
      wetmask = 1                ! default wet
c     --- add island in middle ---
      do i=1, nlan
         ix = ixlist(i)
         iy = iylist(i)
         wdepth(ix,iy) = 0
         wetmask(ix,iy) = 0 ! dry
      enddo
c     --- lon-lat grid descriptors ---
      lambda1    = 5.0
      dlambda    = 0.08      
      phi1       = 55.0      
      dphi       = 0.1
      lambda_max = lambda1 + (nx-1)*dlambda ! last center 
      phi_max    = phi1    + (ny-1)*dphi    ! last center 
      end subroutine init_mesh_grid


      
      subroutine get_surrounding_box(geo,ix,iy,sx,sy)
c     ------------------------------------------------------
c     >>>>>>>>>>>>>> inlined from horizontal_grid_transformations_lonlatgrid.f
c      
c     From the lon/lat vector geo (which may include the z component)
c     determine the cell association (ix,iy) where the node of the cell is the
c     lower left corner point, so that [ix:ix+1, iy:iy+1[ maps to cell (ix,iy)
c     i.e. a corner centered convention.
c     0 < (sx,sy) < 1 is the intra cell cooridnates. 
c     Grid origo (ix,iy) = (x,y) = (1,1) is associated with (lambda1,phi1)
c     Do not flag grid range violations, so that 
c     0 < (ix,iy) <= (nx,ny) is not enforced 
c     (horizontal_range_check uses this function to check grid range violations)
c     NB: previously, this was subroutine called get_horiz_grid_coordinates
c     ------------------------------------------------------
      real, intent(in)     :: geo(:)
      integer, intent(out) :: ix,iy
      real, intent(out)    :: sx,sy
      real                 :: dx1,dy1
c     ------------------------------------------------------
      dx1 = (geo(1)-lambda1)/dlambda
      dy1 = (geo(2)-phi1)   /dphi
      ix  = 1 + floor(dx1)    ! truncate decimals to get lon cell
      iy  = 1 + floor(dy1)    ! truncate decimals to get lat cell
      sx  = dx1 - floor(dx1)  ! intra cell coordinate 0<sx<1
      sy  = dy1 - floor(dy1)  ! intra cell coordinate 0<sy<1
c     ------------------------------------------------------
      end subroutine get_surrounding_box

      

      LOGICAL function is_land(geo)        ! COPY -> horizontal_representation_stay_wet.f
c     ------------------------------------------
c     >>>>>>>>>>>>>> inlined from horizontal_representation_traj_staywet.f
c
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

      
      
      LOGICAL function horizontal_range_check(geo)
c     ------------------------------------------
c     >>>>>>>>>>>>>> inlined from horizontal_grid_transformations_lonlatgrid.f

c     Require grid coordinates 
c      (0.5 < x < nx+0.5) and (0.5 < y < ny+0.5)              
c     ------------------------------------------ 
      real, intent(in) :: geo(:)
      integer          :: ixc,iyc
c     ------------------------------------------ 
      call get_horiz_ncc_index(geo,ixc,iyc) 
      horizontal_range_check = (1<=ixc).and.(ixc<=nx)
     +                    .and.(1<=iyc).and.(iyc<=ny)
c     ------------------------------------------ 
      end function horizontal_range_check


      

      subroutine get_horiz_ncc_index(geo,ixc,iyc) !,xc,yc)
c     ------------------------------------------
c     >>>>>>>>>>>>>> inlined from horizontal_grid_transformations_lonlatgrid.f
c
c     Resolve the node-centered grid cell containing geo
c     (ixc,iyc) is grid indices of the closest node
c     (xc,yc)   is (lon, lat) coordinates of the closest node
c     Notice that (ixc,iyc) is not coinciding with (ix,iy)
c     in get_horiz_grid_coordinates(geo,ix,iy,sx,sy) above
c     no range check
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      integer, intent(out) :: ixc,iyc
c     ------------------------------------------ 
      ixc = 1 + nint((geo(1)-lambda1)/dlambda)     
      iyc = 1 + nint((geo(2)-phi1)   /dphi)   
      end subroutine get_horiz_ncc_index

      

      subroutine interpolate_wdepth(geo, result, status) 
c     -------------------------------------------------- 
c     Interpolate water depth from 2D grid of depths
c     Do not apply interpolate_cc_2Dgrid_data which uses wet-point interpolation
c     This subroutine has two modes, controlled by module switch dirichlet_coastal_wdepth_BC :
c      
c       dirichlet_coastal_wdepth_BC = .true. : enforce wdepth goes linearly to zero
c         toward numerical shore lines
c         This is done inserting zero clamp points on a subgrid at numerical shore line 
c      
c       dirichlet_coastal_wdepth_BC = .false.: use unrestricted interpolation using wdepth=0
c         at dry points. This implies coastal cliffs along numerical shore line,
c         where wdepth in general is nonzero
c      
c     The interpolation subgrids nodes are indexed as follows:
c     vc = (2,4)         vcx = (2,7,4)
c          (1,3)               (6,5,9)
c                              (1,8,3)
c     so corners of vcx corresponds to vc. Nodes (5,6,7,8,9) are possible clamp points
c     If either (1,2,3,4) are dry, 5 is always clamped to zero. For wet x/y directed pairs
c     of (1,2,3,4) corresponding nodes in (6,7,8,9) are interpolated linearly.    
c     If all nodes (1,2,3,4) are wet, result is the same irrespoctive setting of
c     dirichlet_coastal_wdepth_BC. Previous version of interpolate_wdepth corresponded to
c     dirichlet_coastal_wdepth_BC = .false.
c      
c     LOG: 22 Sep 2025: added dirichlet_coastal_wdepth_BC = .true. variant 
c     -------------------------------------------------- 
      real, intent(in)     :: geo(:) 
      real, intent(out)    :: result
      integer, intent(out) :: status   
c
      integer           :: ix,iy,ix0,iy0,ix1,iy1
      real              :: sx,sy,vc(4),vcsub(9), vcx(4)
      logical           :: any_dry
c     ------------------------------------------ 
      if (.not.horizontal_range_check(geo)) then
          result = padval_wdepth
          status = 1
          return
      endif
      if (is_land(geo)) then
          result = 0.0 ! fixed value for dry points
          status = 3
          return
      endif
c     --- delegate normal interior interpolation to interp_2Dbox_data
      call get_surrounding_box(geo,ix,iy,sx,sy)
      ix0 = min(max(ix,    1),nx)  ! cover boundary layers
      ix1 = min(max(ix + 1,1),nx)  ! cover boundary layers
      iy0 = min(max(iy,    1),ny)  ! cover boundary layers
      iy1 = min(max(iy + 1,1),ny) ! cover boundary layers
      any_dry = any(wetmask(ix0:ix1, iy0:iy1)==0)
      vc(1) = wdepth(ix0, iy0)
      vc(2) = wdepth(ix0, iy1)
      vc(3) = wdepth(ix1, iy0)
      vc(4) = wdepth(ix1, iy1)
c      
      if (dirichlet_coastal_wdepth_BC .and. any_dry) then
         vcsub      = 0 ! center point 5 of subgrid, always 0 when any corner is dry
         vcsub(1:4) = vc
         if (wetmask(ix0,iy0)*wetmask(ix0,iy1) == 1) ! both nodes 1,2 wet -> subnode nonzero
         
     +        vcsub(6) = 0.5*(vc(1)+vc(2))
         
         if (wetmask(ix0,iy1)*wetmask(ix1,iy1) == 1) ! both nodes 2,4 wet -> subnode nonzero

     +        vcsub(7) = 0.5*(vc(2)+vc(4))
         
         if (wetmask(ix0,iy0)*wetmask(ix1,iy0) == 1) ! both nodes 1,3 wet -> subnode nonzero
         
     +        vcsub(8) = 0.5*(vc(1)+vc(3))
         
         if (wetmask(ix1,iy0)*wetmask(ix1,iy1) == 1) ! both nodes 3,4 wet -> subnode nonzero
     +        vcsub(9) = 0.5*(vc(3)+vc(4))
         
         sx = 2*sx  ! NB: now 0 < sx < 2
         sy = 2*sy  ! NB: now 0 < sy < 2
c        
         if     ((sx <= 1).and.(sy <= 1)) then ! interpolation point in sector 1
            
            vcx = (/vcsub(1), vcsub(6), vcsub(8), vcsub(5)/)
            call interp_2Dbox_data(sx, sy, vcx, 0, result)

         elseif ((sx <= 1).and.(sy  > 1)) then ! interpolation point in sector 2
            
            vcx = (/vcsub(6), vcsub(2), vcsub(5), vcsub(7)/)
            call interp_2Dbox_data(sx, sy-1, vcx, 0, result)

         elseif ((sx  > 1).and.(sy <= 1)) then ! interpolation point in sector 3
            
            vcx = (/vcsub(8), vcsub(5), vcsub(3), vcsub(9)/)
            call interp_2Dbox_data(sx-1, sy, vcx, 0, result)

         elseif ((sx  > 1).and.(sy  > 1)) then ! interpolation point in sector 4
           
            vcx = (/vcsub(5), vcsub(7), vcsub(9), vcsub(4)/)
            call interp_2Dbox_data(sx-1, sy-1, vcx, 0, result)

         else
            write(*,*) "interpolate_wdepth: unhandled case sx,sy=",sx,sy
            stop 32          
         endif
         
      else  ! all corners of interpolation rectangle wet and/or no dirichlet BC
         call interp_2Dbox_data(sx, sy, vc, 0, result)
         status = 0
      endif
c     ------------------------------------------ 
      end subroutine interpolate_wdepth
      end module

 

      subroutine interp_2Dbox_data(sx,sy,vc,deriv,result)
c     ------------------------------------------------------------
c     >>>>>>>>>>>>>> inlined from  grid_interpolations.f, outside module
c      
c     Basic bilinear interpolation of 2D grid data
c     deriv = 0 gives value, deriv = (1,2) gives derivative (wrt. sx,sy) along (x,y)
c     Currently, only deriv = 0 is implemented, until other derivatives 
c     are needed.
c     ------------------------------------------------------------
      implicit none
      real,intent(in)  :: sx,sy
      real,intent(in)  :: vc(*) ! 4 corners values : (v00,v01,v10,v11)
      integer,intent(in) :: deriv
      real,intent(out) :: result
      real             :: v00,v01,v10,v11, vl,vu
c     ------------------------------------------------------------
      v00 = vc(1)
      v01 = vc(2)
      v10 = vc(3)
      v11 = vc(4)

      vl  = (1-sx)*v00 + sx*v10
      vu  = (1-sx)*v01 + sx*v11

      if     (deriv == 0) then
         result = (1-sy)*vl + sy*vu
      else
         write(*,*) "interp_2Dbox_data: unhandled request deriv=",deriv
         stop 
      endif
c     ------------------------------------------------------------
      end subroutine interp_2Dbox_data

      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     test code
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program test_coast_intp
      use mesh_grid_env
      use netcdf
      implicit none
      integer            :: i,j
      integer, parameter :: mx = 1000 ! plotting grid
      integer, parameter :: my = 1000 ! plotting grid
      real               :: lon(mx), lat(my), wdepth(mx,my),dx,dy
      integer            :: ncid, mx_dimid, my_dimid,dimids(2)
      integer            :: lon_varid,lat_varid,varid,ist
      CHARACTER(LEN=50)  :: outfile = "wdepth.nc"
c     --------------------------------------------------------------------------------      
      call init_mesh_grid()
c     ------ setup outputfile ------      
      CALL check(nf90_create(outfile, NF90_CLOBBER, ncid))
      CALL check(nf90_def_dim(ncid,"mx",mx,mx_dimid))
      CALL check(nf90_def_dim(ncid,"my",my,my_dimid))
      CALL check(nf90_def_var(ncid,"lon",NF90_REAL,mx_dimid,lon_varid))
      CALL check(nf90_def_var(ncid,"lat",NF90_REAL,my_dimid,lat_varid))
      dimids = (/ mx_dimid, my_dimid /)
      CALL check(nf90_def_var(ncid, "wdepth",NF90_FLOAT,dimids,varid))
      CALL check(nf90_enddef(ncid))
c     ------ interpolate data ------
      dx = (lambda_max - lambda1) / (mx-1)         
      do i = 1, mx
         lon(i) = lambda1 + (i-1)*dx 
      enddo
c      
      dy = (phi_max    - phi1)    / (my-1)     
      do j = 1, my
         lat(j) = phi1 + (j-1)*dy 
      enddo
c
      do i = 1, mx
         do j = 1, my
            CALL interpolate_wdepth((/lon(i),lat(j)/),wdepth(i,j),ist)
         enddo
      enddo
c     ------ write data ------      
      CALL check(nf90_put_var(ncid, lon_varid, lon))
      CALL check(nf90_put_var(ncid, lat_varid, lat))
      CALL check(nf90_put_var(ncid, varid, wdepth))
      CALL check(nf90_close(ncid))
c     ---------- 
      contains
c     ---------- 
      subroutine check(anint)
      integer,intent(in) :: anint
      if (anint /= 0) then
         write(*,*) "check - ierr=", anint
      endif
      end subroutine check
c      
      end program test_coast_intp
