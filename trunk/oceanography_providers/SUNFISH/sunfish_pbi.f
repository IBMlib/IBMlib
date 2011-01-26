ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     SUNFISH pbi 
c     ---------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
c     PBI interface for SUNFISH cmod + NPZD data
c
c     Coast line reconstruction: this module features
c     a consistent definition of coast_line_intersection()
c     interpolate_wdepth(), and is_land() so that coast line is
c     characterized by wdepth=0 piecewise lines
c
c     Uses site installation of NetCDF
c
c     TODO:  
c     Currently turbulence + deriv is fixed to zero - update
c       test interpolations
c       validate w sign
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module physical_fields

      use regular_lonlat_grid, 
     +          interpolate_currents_cc => interpolate_currents  ! use face-centered local

      use geometry
      use run_context, only: simulation_file
      use time_tools           ! import clock type
      use input_parser
      use netcdf   ! use site installation

      implicit none
      private     


      public :: init_physical_fields     
      public :: close_physical_fields
      public :: get_master_clock
      public :: set_master_clock
      public :: update_physical_fields
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
      public :: get_pbi_version

c     -------------------- module data --------------------  
      
      
c
c     ------ data frame handler ------
      
      character*999              :: hydroDBpath ! hydrographic data sets

      integer, parameter         :: tag_lenght = 10 ! YYYYMMDDHH
      character*(*), parameter   :: not_set_c = "not set"
      integer, parameter         :: not_set_i = -9999
      type(clock)                :: ref_clock  ! reference offset for hydrographic data

      logical                    :: data_in_buffers     = .false.  ! first time setting
      logical                    :: data_files_are_open = .false.  ! first time setting
      character(len=tag_lenght)  :: cur_tag      ! tag YYYYMMDDHH of current open set
      integer                    :: cur_h1900    ! current time frame in buffer
      integer,allocatable        :: h1900_map(:) ! map of time frames in current open set
      integer                    :: ncid_u       ! u NetCDF file handler
      integer                    :: ncid_v       ! v NetCDF handler
      integer                    :: ncid_w       ! w NetCDF handler
      integer                    :: ncid_dslm    ! sea leve  l NetCDF handler
      integer                    :: ncid_t       ! temperature NetCDF handler
      integer                    :: ncid_zoo     ! zooplankton NetCDF handler 
c      integer                    :: ncid_s     : salinity currently unused
c      integer                    :: ncid_light : light currently unused

c     --- 3D grids ---
                
      real,allocatable,target :: ccdepth0(:,:,:)   ! reference cell center depth water below surface [m] (dslm=0)
      real,allocatable,target :: acc_width0(:,:,:) ! accumulated water above this undisturbed layer [m] dim=nz+1   
      
c     --- 2D grids ---
      
      real,allocatable :: wdepth0(:,:)      ! reference depth at cell-center, positive  [m] - dry points negative
 

c     --- 1D grids ---

      real, allocatable :: levels(:)        ! center depth of layers for DSLM=0. [m] dim=nz 
      real, allocatable :: layer_width(:)   ! width for undisturbed layers [m] dim=nz 
      real, allocatable :: layer_faces(:)   ! accumulated width for undisturbed layers above this [m] dim=nz+1 
      
c     ===================================================
                            contains
c     ===================================================
  
      subroutine init_physical_fields(time)
c     ------------------------------------------
c     Do not trigger data load
c     ------------------------------------------
      type(clock), intent(in),optional :: time
c     ------------------------------------------
      if (present(time)) master_clock = time
      write(*,*) trim(get_pbi_version()) 

      call read_control_data(simulation_file,"hydroDBpath",hydroDBpath)
      write(*,*) "init_physical_fields: hydrographic database path =", 
     +           trim(hydroDBpath)
    
      call read_grid_desc()
c.....Set reference offset for hydrographic data: 1900-01-01 00:00:00.
      call set_clock(ref_clock, 1900, 01, 01, 0)

      call reset_frame_handler()

      end subroutine init_physical_fields

 

      character*100 function get_pbi_version()  
      get_pbi_version =  "SUNFISH pbi version: $Rev$"
      end function

   

      subroutine read_grid_desc()
c     ---------------------------------------------------
c     Read grid definition which is contained in file 
c     given as tag grid_desc in the simulation file simulation_file
c     Allocate grid arrays in this subroutine. 
c
c     grid dimensions:         nx,ny,nz
c     grid coordinate map:     lambda1, dlambda; phi1, dphi
c
c     grid point (ix,iy) = (1,1) is at (lambda1,phi1)
c     ---------------------------------------------------
      character*999        :: gd_fname
      type(control_file)   :: grid_ctrlfile
      integer              :: ix,iy,iz,ihit,idx,nwet,nacc,nbott
      real                 :: wd0
c     ---------------------------------------------------
      call read_control_data(simulation_file, "grid_desc", gd_fname)
      write(*,229) "read_grid_desc: reading grid description from: ", 
     +           trim(adjustl(gd_fname))

      call open_control_file(gd_fname, grid_ctrlfile)
    
c.....grid scale/dimensions (held globally in module regular_lonlat_grid)
      
      call read_control_data(grid_ctrlfile,"lambda1", lambda1)
      call read_control_data(grid_ctrlfile,"dlambda", dlambda) 
      call read_control_data(grid_ctrlfile,"phi1",    phi1)
      call read_control_data(grid_ctrlfile,"dphi",    dphi) 
      call read_control_data(grid_ctrlfile,"nx",      nx) 
      call read_control_data(grid_ctrlfile,"ny",      ny) 
      call read_control_data(grid_ctrlfile,"nz",      nz) 
      write(*,231) nx, ny, nz
      write(*,232) lambda1, dlambda
      write(*,233) phi1,    dphi 
    
 229  format(a,a)    
 231  format("read_grid_desc: 3d grid dim (nx,ny,nz) = ", i4,i4,i4)
 232  format("read_grid_desc: lambda1 = ",f12.7,
     +                      " dlambda = ",f12.7," degrees")  
 233  format("read_grid_desc: phi1    = ",f12.7,
     +                      " dphi    = ",f12.7," degrees")  
     

      write(*,*) "read_grid_desc: allocate grid arrays: begin"

      call init_regular_lonlat_grid()  ! 

c     --- allocate specific auxillary arrays ---      

c     --- 3D grids ---
      
      allocate( ccdepth0(nx,ny,nz)) 
      allocate( acc_width0(nx,ny,nz+1)  )  
     
c     --- 2D grids ---
     
      allocate( wdepth0(nx,ny)    ) 
     
c     --- 1D grids ---
      
      allocate( levels(nz)        ) 
      allocate( layer_width(nz)   )  ! width for undisturbed layers in meters
      allocate( layer_faces(nz+1) )  ! meters of undisturbed water above layer
      
      write(*,*) "read_grid_desc: allocate grid arrays: OK"   


c.....vertical layer structure: levels -> layer_width0 + acc_width0

      call read_control_data(grid_ctrlfile,"levels", levels) 
      write(*,241) levels
 241  format("read_grid_desc: read levels = ",1000f6.1)
     
c     set acc_width0, layer_width0 from levels
c     last element acc_width0(nz+1) is total (dslm=0) max depth of grid
      layer_faces(1)    = 0.
      acc_width0(:,:,1) = 0. ! sea surface (DSLM=0)
      do iz=1,nz 
         layer_width(iz)      = 2.0*(levels(iz) - layer_faces(iz))
         layer_faces(iz+1)    = layer_faces(iz) + layer_width(iz)
         ccdepth0(:,:,iz)     = levels(iz)
         acc_width0(:,:,iz+1) = layer_faces(iz+1)
      enddo
      acc_width = acc_width0 ! allow certain grid functions to work on land
      ccdepth   = ccdepth0   ! allow certain grid functions to work on land

      write(*,242) layer_width
      write(*,243) layer_faces
 242  format("read_grid_desc: calc ref layer_width = ",1000f6.1)     
 243  format("read_grid_desc: calc ref layer_faces = ",1000f6.1)
c
c.....read topography -> wdepth0 (only wety points provided, other assumed dry)
c
c     1) set bottom_layer
c     2) adjust ccdepth0, so that actual width of bottom layer is consistent 
c        with wdepth0 by stretching the last wet layer to the topographic bottom
c         
c
      write(*,*) "read_grid_desc: reading topography"
      wdepth0      = -1.0       ! default dry
      bottom_layer =  0         ! default dry
      ihit         =  1         ! start scanning for is_wet at ihit=1
      call read_control_data(grid_ctrlfile,"topography_points",nwet,
     +                       ihit)                     ! locate file position: ihit
      nacc   = 0
      nbott  = 0
      do idx = 1, nwet
         read(grid_ctrlfile%rawdata(ihit+idx),*) ix, iy, wd0
c        --- horizontal check ---
         if ((ix>nx).or.(ix<1)) then
            write(*,*) "topography point", idx, 
     +                 " exceeds x limits - point ignored"
            cycle
         endif
         if ((iy>ny).or.(iy<1)) then
            write(*,*) "topography point", idx, 
     +                 " exceeds y limits - point ignored"
            cycle
         endif
c
c        locate vertical cell iz containing wd0: 
c        layer_faces(iz) < wd0 < layer_faces(iz+1) 
c        iz=0    -> wd0 negative
c        iz=nz+1 -> wd0 below sea bed
c
         call search_sorted_list(wd0, layer_faces, iz) 

c        --- check of iz  ---
         if (iz<1) then
            write(*,*) "topography point", idx, 
     +                 "has negative depth  - point ignored"
            cycle
         endif
         if (iz>nz) then
            write(*,*) "topography point", idx, 
     +                 " exceeds grid z limit - bottum lowered"           
            iz     = nz
            nbott  = nbott + 1 
         endif
c
c        --- accept this point, if not bounced above ---
c            after here we know iz is appropriate: 1<=iz<=nz 
c            adjust grid bottum by stretching the last wet layer
         nacc                    = nacc + 1
         wdepth0(ix, iy)         = wd0
         bottom_layer(ix, iy)    = max(iz-1, 1)     
         acc_width0(ix, iy,iz+1) = wd0              
         ccdepth0(ix, iy, iz)    = 0.5*(acc_width0(ix, iy,iz) + 
     +                                  acc_width0(ix, iy,iz+1)) ! re-center mid point of bottom cell
c
      enddo ! idx = 1, nwet

      write(*,288) nwet, nacc
      write(*,289) nbott
 288  format("read_grid_desc: read",i6," wet points - ",i6, 
     +       " points accepted")
 289  format("read_grid_desc: grid bottum lowered in ",i6, " points")     
c
c     --- define the wetmask statically (no flooding/drying) --- 
c   
      where(wdepth0 > 0.0)
         wetmask = 1 ! wet
      elsewhere
         wetmask = 0 ! dry
      end where 

      write(*,*) "read_grid_desc: grid arrays initialized"

      call close_control_file(grid_ctrlfile)
     
c     ------------------------------------------       
      end subroutine read_grid_desc



      subroutine close_physical_fields()
c     ------------------------------------------ 
c     ------------------------------------------    
      if (allocated(h1900_map))    deallocate( h1900_map )  
            
      if (allocated(ccdepth0))     deallocate( ccdepth0 )
      if (allocated(acc_width0))   deallocate( acc_width0 )
      if (allocated(acc_width))    deallocate( acc_width )
      if (allocated(wdepth0))      deallocate( wdepth0 )
 
      if (allocated(levels))       deallocate( levels )
      if (allocated(layer_width))  deallocate( layer_width )
      if (allocated(layer_faces))  deallocate( layer_faces )

      call close_regular_lonlat_grid()
      call reset_frame_handler()
c     ------------------------------------------
      end subroutine 


      subroutine update_physical_fields(time, dt)
c     ------------------------------------------  
      type(clock), intent(in),optional :: time
      integer, intent(in),optional     :: dt
      logical                          :: update
      character(len=tag_lenght)        :: tag
      integer                          :: h1900
c     ------------------------------------------  
      if (present(time)) then
         call set_master_clock(time)
      elseif (present(dt)) then
         call add_seconds_to_clock(master_clock, dt)
      endif
c
      call resolve_corresp_dataset(master_clock, tag, h1900)
      call update_dataset(tag, h1900)
c     ------------------------------------------       
      end subroutine update_physical_fields




      subroutine resolve_corresp_dataset(aclock, tag, h1900)
c     ------------------------------------------
c     Resolve tag and hours since 1900-01-01 00:00:00" (== dataframe)
c     corresponding to aclock
c
c     Data frames are stored according to end point of
c     the hour interval, i.e time = 912385h represents
c     time interval ]912384h; 912385h] since 1900-01-01 00:00:00.
c     [ref == mam@dmu.dk; ma 11-10-2010 12:13]
c
c     In the SUNFISH set, data are bundled by months
c     YYYYMM0100 refers to the beginning of the month corresponding to the data
c     so that YYYYMM0100 is the tag corresponding to time YYYYMMDDHH
c     ------------------------------------------
      type(clock), intent(in)                :: aclock
      character(len=tag_lenght), intent(out) :: tag
      integer, intent(out)                   :: h1900
      integer          :: year, month, day
      real             :: h1900_r
c     ------------------------------------------
      call get_date_from_clock(aclock, year, month, day)
      write(tag,455) year, month, 1, 0   ! beginning of that month
      call get_period_length_hour(ref_clock, aclock, h1900_r)
      h1900 = nint(h1900_r + 0.5)     
 455  format(i4.4,3i2.2)
      end subroutine resolve_corresp_dataset




      subroutine update_dataset(tag, h1900)
c     ------------------------------------------
c     Only open new data set file if necessary
c     Only load new frames if necessary
c     ------------------------------------------
      character(len=tag_lenght), intent(inout) :: tag
      integer, intent(inout)                   :: h1900
      logical  :: needed
c     ------------------------------------------
      needed = .not.data_in_buffers
c
      if ((tag /= cur_tag).or.needed) then
         call reset_frame_handler() ! resets cur_h1900+close old
         call open_data_files(tag)
         cur_tag = tag
         needed = .true. ! force update if accidentally h1900 == cur_h1900
      endif
c      
      if ((h1900 /= cur_h1900).or.needed) then
         call load_data_frames(h1900)
      endif  
c 
      end subroutine update_dataset





      subroutine reset_frame_handler()
c     ------------------------------------------
      call close_data_files()
      data_in_buffers = .false.
      cur_tag         = not_set_c
      cur_h1900       = not_set_i   
      end subroutine reset_frame_handler



      subroutine close_data_files()
c     ------------------------------------------
c     Only close, if data files are open
c     
c     ------------------------------------------
      if (data_files_are_open) then
         call NetCDFcheck( nf90_close(ncid_u) )
         call NetCDFcheck( nf90_close(ncid_v) )
         call NetCDFcheck( nf90_close(ncid_w) )
         call NetCDFcheck( nf90_close(ncid_dslm) )
         call NetCDFcheck( nf90_close(ncid_t) )
         call NetCDFcheck( nf90_close(ncid_zoo) )
         ncid_u    = not_set_i  
         ncid_v    = not_set_i  
         ncid_w    = not_set_i  
         ncid_dslm = not_set_i  
         ncid_t    = not_set_i  
         ncid_zoo  = not_set_i          
         write(*,*) "close_data_files: closed open NetCDF set"
         data_files_are_open = .false.
      endif
      end subroutine close_data_files


      subroutine open_data_files(tag)
c     ------------------------------------------
c     Open netcdf file set corresponding to 
c     tag == YYYYMMDDHH and load the
c     time map h1900_map of data in this set
c     
c     tempdat.u_YYYYMMDDHH.00.nc  -> ncid_u
c     tempdat.v_YYYYMMDDHH.00.nc  -> ncid_v
c     tempdat.w_YYYYMMDDHH.00.nc  -> ncid_w
c     tempdat.z_YYYYMMDDHH.00.nc  -> ncid_dslm
c     tempdat.t_YYYYMMDDHH.00.nc  -> ncid_t
c     tempdat.s_YYYYMMDDHH.00.nc  -> (ncid_s     : currently unused)
c     biodat.light_YYYYMMDDHH.nc  -> (ncid_light : currently unused)
c     biodat.zoo_YYYYMMDDHH.nc    -> ncid_zoo
c     ------------------------------------------
      character(len=tag_lenght),intent(in) :: tag
      character*999                        :: fname
      integer :: timedimID, ntime, varid
c     ------------------------------------------
      call close_data_files() ! in case they are open ...
      write(fname,333) trim(adjustl(hydroDBpath)), "u", tag
      call NetCDFcheck( nf90_open(fname, NF90_NOWRITE, ncid_u) )
      write(fname,333) trim(adjustl(hydroDBpath)), "v", tag
      call NetCDFcheck( nf90_open(fname, NF90_NOWRITE, ncid_v) )
      write(fname,333) trim(adjustl(hydroDBpath)), "w", tag
      call NetCDFcheck( nf90_open(fname, NF90_NOWRITE, ncid_w) )
      write(fname,333) trim(adjustl(hydroDBpath)), "z", tag
      call NetCDFcheck( nf90_open(fname, NF90_NOWRITE, ncid_dslm) )
      write(fname,333) trim(adjustl(hydroDBpath)), "t", tag
      call NetCDFcheck( nf90_open(fname, NF90_NOWRITE, ncid_t) )
      write(fname,334) trim(adjustl(hydroDBpath)), "zoo", tag
      call NetCDFcheck( nf90_open(fname, NF90_NOWRITE, ncid_zoo) )
c
c     load  time map h1900_map of data (asserted identical for all in this set)
c
      if (allocated(h1900_map)) deallocate(h1900_map)
      call NetCDFcheck( nf90_inq_dimid(ncid_dslm, "time", timedimID) )
      call NetCDFcheck( nf90_inquire_dimension(ncid_dslm, timedimID, 
     +                  len = ntime) )
c     
      if (ntime>0) then
         allocate( h1900_map(ntime) )
         call NetCDFcheck( nf90_inq_varid(ncid_dslm, "time", varid) )
         call NetCDFcheck( nf90_get_var(ncid_dslm, varid, h1900_map))
      else
         write(*,*) "open_data_files error: NetCDF set ", tag,
     +              " contains ",ntime, "frames"
      endif
c
      write(*,*) "open_data_files: opened NetCDF set ", tag,
     +           " containing ",ntime, "frames"  
c  
      data_files_are_open = .true.
c     
 333  format(a,'/tempdat.',a,'_',a10,'.00.nc')
 334  format(a,'/biodat.',a,'_',a10,'.nc')
      end subroutine open_data_files



      subroutine load_data_frames(h1900)
c     ------------------------------------------
c     Verbose reading af NetCDF sets consecutively to grid buffers 
c     from open NetCDF sets (referenced by data frame handler)
c     corresponding to h1900 hours since 1900-01-01 00:00:00 (ref_clock)
c     Syncronize auxillary grid fields after data load. 
c     It is assumed that grids do not
c     change during a simulation (assertion not checked)
c     
c     Currently load these fields (in native units/conventions):
c
c     float u(time, depth, lat, lon) ;  "m/s"  positive east
c     float v(time, depth, lat, lon) ;  "m/s"  positive north 
c     float w(time, depth, lat, lon) ;  "m/s", positive up 
c     float z(time, lat, lon) ;         "m"    positive down (!)
c     float t(time, depth, lat, lon) ;  "degC" 
c     float zoo(time, depth, lat, lon) ;"Zooplankton: 10^-3 mol N/liter" 
c                
c     The is currently no horizontal diffusivity in the set (set it to zero)
c
c     fortran index order opposite CDL dump order
c     fortran: fastests index left most
c     ------------------------------------------
      integer, intent(in) :: h1900 ! pick frame, corresponding to h1900 hours since ref_clock 
      integer :: iframe 
      integer :: start3D(4),start2D(3),varid
      integer :: ix,iy,iz,no_fill
      real    :: dz
      logical :: not_ok
      real    :: u_fill,v_fill,w_fill
c     ------------------------------------------ 
c
c     1) locate frame to pick - time frame map h1900_map already loaded in open_data_files
c
      iframe = 1
      do while ((iframe <= size(h1900_map)).and.
     +         (h1900 /= h1900_map(iframe)))
         iframe = iframe + 1
      enddo
      if (iframe > size(h1900_map)) then
         write(*,*) "load_data_frames: unable to locate h1900 =",h1900
         write(*,*) "frames in present set =", h1900_map
         stop "load_data_frames: load error"
      endif
c
c     2) load iframe from set 
c
c     ---- 2D ----
      start2D(:) = 1
      start2D(3) = iframe
      call NetCDFcheck( nf90_inq_varid(ncid_dslm,     "z",   varid) )
      call NetCDFcheck( nf90_get_var(ncid_dslm, varid, dslm, start2D) )
      
c     ---- 3D ----
      start3D(:) = 1
      start3D(4) = iframe
      call NetCDFcheck( nf90_inq_varid(ncid_u,     "u", varid) )
      call NetCDFcheck( nf90_get_var(ncid_u, varid, u, start3D))   
c      
      call NetCDFcheck( nf90_inq_varid(ncid_v,     "v", varid) )
      call NetCDFcheck( nf90_get_var(ncid_v, varid, v, start3D))
c
      call NetCDFcheck( nf90_inq_varid(ncid_w,     "w", varid) )
      call NetCDFcheck( nf90_get_var(ncid_w, varid, w, start3D))
c
      call NetCDFcheck( nf90_inq_varid(ncid_t,     "t", varid) )
      call NetCDFcheck( nf90_get_var(ncid_t, varid, temp, start3D))   
c
      call NetCDFcheck( nf90_inq_varid(ncid_zoo,     "zoo", varid) )
      call NetCDFcheck( nf90_get_var(ncid_zoo, varid, zoo, start3D))   

c

c     3) postprocess data, suncronize auxillary fields, pad holes
c      
c     ---- currently no horizontal turbulent diffusivity: set to
c          molecular diffusivity lower limit ~ 1.e-9 m2/s   ---
c          NB: horizontal derivatives NOT implemented
c          when spatially non constant hdiffus are applied,
c          interpolate_turbulence_deriv must be updated
c
      vdiffus = 1.e-9
      hdiffus = 1.e-9   
c
c     4) replace fill values for (u,v,w), because boundary faces are given the fill value
c        and boundary faces are needed for interpolation. Implicit boundary condition is 
c        (u,v,w)=0. Unfortunately the nf90 fill query appears broken so apply fill=100
c        as detection limit (hack - unable to do better right now)
c
      where(abs(u)>100.0) u=0.0
      where(abs(v)>100.0) v=0.0
      where(abs(w)>100.0) w=0.0
c
c     ------ flip sign of w 
c     ------ physical_fields sign convention is positive down, cmod convention is positive up
c
      w = -w
c     ------ add DSLM to auxillary grid descriptors for wet points  
c            acc_width(ix,iy,1) = 0 (sea surface)  
      do ix=1,nx
         do iy=1,ny
            if (wetmask(ix,iy)==0) cycle ! nothing for dry points
            dz = dslm(ix,iy) ! change of sea level in this point, positive down
            wdepth(ix,iy)       = wdepth0(ix,iy)       - dz
            acc_width(ix,iy,2:) = acc_width0(ix,iy,2:) - dz
            ccdepth(ix,iy,1)    = ccdepth0(ix,iy,1)    - dz*0.5 
            ccdepth(ix,iy,2:)   = ccdepth0(ix,iy,2:)   - dz          
         enddo
      enddo
c      
      write(*,*) "load_data_frames: loaded frame:",
     +            iframe,"/", size(h1900_map)
      data_in_buffers = .true.
      cur_h1900 = h1900


      end subroutine load_data_frames



      subroutine interpolate_currents(xyz, uvw, status)
c     ------------------------------------------
c     Perform linear two-face interpolation in current fields into vector uvw
c     
c     xyz is the interpolation position as (lon[degE],lat[degN],depth[m,positive down])
c     uvw are interpolated currents in m/s, shape == real(3+)
c     status is the status of the interpolation action
c       status = 0: all OK
c       status = 1: domain violation (but extrapolation provided)
c     
c     ------------- data layout -----------------
c
c     u,v positive along lambda,phi in units m/s
c     w   positive downward         in units m/s
c    
c     u(ix,iy,iz) in grid position (ix+0.5, iy    , iz)     (i.e. eastern  cell face)
c     v(ix,iy,iz) in grid position (ix    , iy-0.5, iz)     (i.e. southern cell face)
c     w(ix,iy,iz) in grid position (ix    , iy,     iz-0.5) (i.e. upper    cell face)
c
c     The boundary condition:
c           w( ix, iy, bottom_layer(ix,iy)+0.5 ) = 0
c     is implicit 
c
c     Implied interpolation ranges in ncc coordinates:
c           1.5 < x < nx+0.5
c           0.5 < y < ny-0.5
c           0.5 < z < nz+0.5 (due to bottom BC)
c
c     this means that the proper interpolation domain is different from the simulation domain
c     definition derived from the grid. Project these rim points onto the
c     interior domain and interpolate currents here and do not flag them as domain violation
c
c     ------------------------------------------ 
      real, intent(in)     :: xyz(:)  ! (lon,lat,depth)
      real, intent(out)    :: uvw(:)
      integer, intent(out) :: status

      integer              :: statu, statv, statw
      integer              :: ix,iy,iz,jwest,jnorth,jlow,ibot
      real                 :: x,y,z, sx,sy,sz
c     ------------------------------------------ 
      if (.not.horizontal_range_check(xyz)) then
         uvw    = 0.
         status = 1  ! signal range violation
         return
      endif
      if (.not.is_wet(xyz)) then
         uvw    = 0.
         status = 3  ! signal dry point
         return
      endif
      
c.....transform to continuous node-centered grid coordinates
c     and check ranges for this staggering
      call get_ncc_coordinates(xyz,x,y,z)
      statu = 0
      statv = 0 
      statw = 0
      if ((x<1.5).or.(x>(nx+0.5))) statu = 1 ! flag x range violation
      if ((y<0.5).or.(y>(ny-0.5))) statv = 1 ! flag y range violation
      if ((z<0.5).or.(z>(nz+0.5))) statw = 1 ! flag z range violation
c
c     Currently ignore the mismatch between interpolation and simulation domain 
c     original test: status  = max(statu, statv, statw) 
c
      status = 0 ! signal all OK, after range/wet check

c.....determine cell associations of point, constrained to 1 <= i <= n
      ix   = max(1, min(nint(x), nx)) ! u at cell east  face
      iy   = max(1, min(nint(y), ny)) ! v at cell south face
      ibot = bottom_layer(ix,iy) ! 1 <= ibot <= nz
      iz   = max(1, min(nint(z),ibot)) ! w at cell upper face

c.....determine relative intracell coorninate s, constrained to 0 < s < 1 
c     when point exceed coordinate bounds s is assigned 0 or 1
      sx = max(0.d0, min(1.d0, x+0.5d0-real(ix)))  
      sy = max(0.d0, min(1.d0, y+0.5d0-real(iy)))           
      sz = max(0.d0, min(1.d0, z+0.5d0-real(iz)))  

c.....face-to-face interpolation 
c     cell indices satisfy 1 <= (ix,iy,iz) <= (nx,ny,ibot)
      jwest  = max(1, min(nint(x-1), nx)) ! cell west  face
      jnorth = max(1, min(nint(y+1), ny)) ! cell north fce
      uvw(1) = sx*u(ix, iy, iz)    + (1.-sx)*u(jwest, iy, iz)
      uvw(2) = sy*v(ix, jnorth,iz) + (1.-sy)*v(ix, iy, iz)

      if (iz==ibot) then
         jlow = 0               ! for debugging
         uvw(3) = (1.-sz)*w(ix, iy, ibot) ! w=0 at sea bed, ibot corresponds to upper cell face
      else
         jlow = max(1, min(nint(z+1),ibot)) ! cell lower face
         uvw(3) = sz*w(ix, iy, jlow) + (1.-sz)*w(ix, iy, iz)
      endif
c     ------------------------------------------ 
      end subroutine interpolate_currents


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     NetCDF auxillaries  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine NetCDFcheck(status)
c     ------------------------------------------
c     This subroutine supports the recommended reading 
c     style of NetCDF4:
c       call NetCDFcheck( nf90_get_var(ncid, varid, data_in) )
c     Private to this module
c     ------------------------------------------
      integer, intent ( in) :: status
    
      if(status /= nf90_noerr) then 
         print *, trim(nf90_strerror(status))
         stop "NetCDFcheck:Stopped"
      end if
      end subroutine NetCDFcheck 


      logical function is_NetCDF_fill_value_real(value)
c     ------------------------------------------------------
c     NetCDF auxillary query function (temporary implementation)
c     nf90_fill_real is defined in netcdf.mod 
c     ------------------------------------------------------
      real :: value
      is_NetCDF_fill_value_real = (value>0.999*nf90_fill_real)
      end function is_NetCDF_fill_value_real


      end module

      
