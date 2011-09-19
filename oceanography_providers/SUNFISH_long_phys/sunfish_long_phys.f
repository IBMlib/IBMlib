ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     SUNFISH (long physics time series) pbi 
c     ---------------------------------------------------
c     $Rev: 353 $
c     $LastChangedDate: 2011-06-16 23:29:14 +0200 (Thu, 16 Jun 2011) $
c     $LastChangedBy: asch $ 
c
c     The data set contains no NPZD data
c    
c     TODO:  
c     Currently turbulence + deriv is fixed to zero - update
c       test interpolationes
c       validate w sign
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module physical_fields

c     plugin for reading native compressed cmod output (rename data, avoid mixing types)
      use read_cmod, u_cmod=>u,v_cmod=>v,w_cmod=>w,z_cmod=>z,s_cmod=>s,
     +               t_cmod=>t,wu_cmod=>wu,wv_cmod=>wv,
     +               m1_cmod=>m1,hz_cmod=>hz,cellw_cmod=>cellw 

      use mesh_grid, 
     +          interpolate_currents_cc => interpolate_currents  ! use face-centered local
      use horizontal_grid_transformations

      use geometry
      use array_tools
      use run_context, only: simulation_file
      use time_services           ! import clock type and time handling
      use input_parser
     

      implicit none
      private     


      public :: init_physical_fields     
      public :: close_physical_fields
      public :: get_master_clock   ! from time_services
      public :: set_master_clock   ! from time_services
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
      
      real,parameter :: molecular_diffusivity = 1.e-9 ! unit m2/s
c
c     ------ data frame handler ------
      
      character*999           :: hydroDBpath      ! path for hydrographic data sets     
      integer                 :: data_set_handler ! which cmod data map will be used

c     --- 3D grids ---
                
      real,allocatable,target :: ccdepth0(:,:,:)   ! reference cell center depth water below surface [m] (dslm=0)
      real,allocatable,target :: acc_width0(:,:,:) ! accumulated water above this undisturbed layer [m] dim=nz+1   
      
c     --- 2D grids ---
      
      real,allocatable :: wdepth0(:,:)      ! reference depth at cell-center, positive  [m] - dry points negative
 
c     --- 1D grids ---

      
c     ===================================================
                            contains
c     ===================================================
  
      subroutine init_physical_fields(time)
c     ------------------------------------------
c     Do not trigger data load
c     ------------------------------------------
      type(clock), intent(in),optional :: time
      real                             :: rdum
c     ------------------------------------------
      if (present(time)) call set_master_clock(time)
      write(*,*) trim(get_pbi_version()) 

      call read_control_data(simulation_file,"hydroDBpath",hydroDBpath)
      write(*,*) "init_physical_fields: hydrographic database path =", 
     +           trim(hydroDBpath)
    
      call read_grid_desc()   ! incl allocation of 3D arrays and init_horiz_grid_transf
c
c     ---- currently no horizontal turbulent diffusivity in data set ----
c          set it to molecular diffusivity lower limit, if no values 
c          are provided
c
c     resolved user provided constant horizontal_diffusivity (optional)  
c    
      if (count_tags(simulation_file, "horizontal_diffusivity")/=0) then
         call read_control_data(simulation_file,
     +                          "horizontal_diffusivity",rdum)
         rdum = max(molecular_diffusivity, rdum) ! never exceed lover limit
         write(*,563) "horizontal_diffusivity",rdum
      else
         rdum = molecular_diffusivity
         write(*,564) "horizontal_diffusivity",rdum
      endif
      hdiffus = rdum
c 
c     resolved user provided constant vertical_diffusivity (optional)  
c    
      if (count_tags(simulation_file, "vertical_diffusivity")/=0) then
         call read_control_data(simulation_file,
     +                          "vertical_diffusivity",rdum)
         rdum = max(molecular_diffusivity, rdum) ! never exceed lover limit
         write(*,563) "vertical_diffusivity",rdum
      else
         rdum = molecular_diffusivity
         write(*,564) "vertical_diffusivity",rdum
      endif
      vdiffus = rdum
c
 563  format("init_physical_fields: read const  ",a,"=", e12.5," m2/s")
 564  format("init_physical_fields: using const ",a,"=", e12.5," m2/s")

      end subroutine init_physical_fields

 

      character*100 function get_pbi_version()  
      get_pbi_version =  "SUNFISH (long physics)"//
     +                   "pbi version: $Rev: 353 $"
      end function

   

      subroutine read_grid_desc()
c     ---------------------------------------------------
c     Grid definitions is $(hydroDBpath)/cfg.nml (which points to auxillary
c     grid descriptors).
c     Allocate grid arrays in this subroutine. 
c
c     Set principal grid parameters:
c
c       horizontal grid dimensions: nx,ny (exported by module horizontal_representation)
c       vertical grid dimensions:   nz    (exported by module mesh_grid)
c       grid coordinate map:        lambda1, dlambda; phi1, dphi (passed to horizontal_grid_transformations)
c
c       grid point (ix,iy) = (1,1) is at (lambda1,phi1)
c     ---------------------------------------------------
      character*99         :: data_set_id
      integer              :: ix,iy,iz,idx1d,j
      
      real                 :: lambda1, dlambda, phi1, dphi ! LOCAL DUMMIES
c     ---------------------------------------------------
      call init_read_cmod(hydroDBpath) ! defined in module read_cmod

      call read_control_data(simulation_file,"cmod_data_set",
     +                       data_set_id)
      write(*,*) "Using cmod data set:", adjustl(trim(data_set_id))

c.....set grid scale/dimensions (held globally in module regular_lonlat_grid)
c     get_grid_descriptors is defined in module read_cmod

      call get_grid_descriptors(data_set_id, data_set_handler,     
     +                          lambda1,dlambda, phi1,dphi, nx,ny,nz) 

      write(*,231) nx, ny, nz
      write(*,232) lambda1, dlambda
      write(*,233) phi1,    dphi 
   
 231  format("read_grid_desc: 3d grid dim (nx,ny,nz) = ", i4,i4,i4)
 232  format("read_grid_desc: lambda1 = ",f12.7,
     +                      " dlambda = ",f12.7," degrees")  
 233  format("read_grid_desc: phi1    = ",f12.7,
     +                      " dphi    = ",f12.7," degrees")  
     

      write(*,*) "read_grid_desc: horizontal initialization"
      call init_horiz_grid_transf(lambda1, phi1, dlambda, dphi)


      write(*,*) "read_grid_desc: allocate grid arrays: begin"
      call init_mesh_grid()  ! incl allocation of 3D arrays
      
c     --- allocate specific auxillary arrays ---      

c     --- 3D grids ---
      
      allocate( ccdepth0(nx,ny,nz)) 
      allocate( acc_width0(nx,ny,nz+1)  )  
     
c     --- 2D grids ---
     
      allocate( wdepth0(nx,ny)    ) 
 
      write(*,*) "read_grid_desc: allocate grid arrays: OK"   

c     ----------------------------------------------------------------
c     set auxillary grid arrays: ccdepth0, acc_width0, bottom_layer, wetmask,wdepth0
c     ----------------------------------------------------------------   
      wdepth0      = -1.0   ! default dry
      bottom_layer =  0     ! default dry
c     do not assign a value to ccdepth0/acc_width0 to be able to spot unknown dry point violations

      acc_width0(:,:,1) = 0. ! sea surface (DSLM=0)    
      do ix=1,nx
         do iy=1,ny
            if (m1_cmod(data_set_handler)%p(ny+1-iy, ix, 1) > 0) then  ! wet/dry longitunal point test  
c             
c              Vertical loop for wet points. At loop exit iz will point to first dry layer
c              F90 std condition for full do-loop termination is iz=nz+1
c            
               do iz=1,nz 
                  idx1d = m1_cmod(data_set_handler)%p(ny+1-iy, ix, iz)
                  if (idx1d > 0) then                             ! wet/dry signal for cmod   
                     ccdepth0(ix,iy,iz) = 
     +                           hz_cmod(data_set_handler)%p(idx1d)
                     acc_width0(ix,iy,iz+1) = acc_width0(ix,iy,iz) + 
     +                           cellw_cmod(data_set_handler)%p(idx1d)
                  else
                     exit 
                  endif   
               enddo
               bottom_layer(ix,iy) = iz-1 ! last wet layer
               wdepth0(ix,iy) = acc_width0(ix,iy,1+bottom_layer(ix,iy))

            endif ! wet/dry longitunal point
         enddo
      enddo

      write(*,*) "read_grid_desc: transferred cmod bathymetry"    
c
c     --- define the wetmask statically (no flooding/drying) --- 
c   
      where(wdepth0 > 0.0)
         wetmask = 1 ! wet
      elsewhere
         wetmask = 0 ! dry
      end where 
c     ------------------------------------------       
      end subroutine read_grid_desc



      subroutine close_physical_fields()
c     ------------------------------------------ 
c     ------------------------------------------    
            
      if (allocated(ccdepth0))     deallocate( ccdepth0 )
      if (allocated(acc_width0))   deallocate( acc_width0 )
      if (allocated(wdepth0))      deallocate( wdepth0 )
 
      call close_read_cmod()
      call close_mesh_grid()
      call close_horiz_grid_transf()
c     ------------------------------------------
      end subroutine 


      subroutine update_physical_fields(time, dt)
c     -----------------------------------------------------------  
c     Map updated time to a file name and a frame in file
c     Files are currently bundled in 12 frames labelled by last frame, e.g. 
c
c       archive.2001031600   contains  2001.03.15 13 to 2001.03.16 00
c       archive.2001031612   contains  2001.03.16 01 to 2001.03.16 12
c
c     file name / frame map tested
c     -----------------------------------------------------------  
      type(clock), intent(in),optional :: time
      integer, intent(in),optional     :: dt
      character(len=18)                :: fname ! archive.2001031600 
      character(len=256)               :: full_fname ! incl path
      type(clock), pointer             :: aclock
      type(clock)                      :: dummy_clock 
      integer                          :: year,month,day,hour,isec
c     ------------------------------------------
      aclock => get_master_clock()
      if (present(time)) then
         call set_master_clock(time)
      elseif (present(dt)) then
         call add_seconds_to_clock(aclock, dt)
      endif    
c 
c     resolve file name
c
      call set_clock(dummy_clock, aclock)
      call add_seconds_to_clock(dummy_clock, -1800)
      call get_date_from_clock(dummy_clock, year, month, day)
      call get_second_in_day(dummy_clock, isec)
      if (isec >= 43200) then  ! must be >=
         write(fname, 824) year, month, day, 0
      else ! get next day by adding 86400 to day  start
         write(fname, 824) year, month, day, 12
      endif
 824  format("archive.",i4.4,i2.2,i2.2,i2.2)
      full_fname = adjustl(trim(hydroDBpath))// adjustl(trim(fname))
 
c
c     resolve time frame to pick from file
c     round toward nearest hour:  hour = intdiv((isec+1800)/3600)
c
      call set_clock(dummy_clock, aclock)
      call add_seconds_to_clock(dummy_clock, 1800)
      call get_date_from_clock(dummy_clock, year, month, day)
      call get_second_in_day(dummy_clock, isec)
      hour = isec/3600 ! integer division
c     
      call update_buffers(full_fname,year,month,day,hour)
      call syncronize_data()
c     ------------------------------------------       
      end subroutine update_physical_fields




      subroutine syncronize_data()
c     ------------------------------------------
c     Syncronize auxillary grid fields after data load. 
c     It is assumed that grids do not
c     change during a simulation (assertion not checked)
c     
c     load all fields in data (except wind):
c         u_cmod(area)%p(:)     :    u current
c         v_cmod(area)%p(:)     :    v current
c         w_cmod(area)%p(:)     :    w current
c         z_cmod(area)%p(:)     :    sea surf elev
c         s_cmod(area)%p(:)     :    salinity
c         t_cmod(area)%p(:)     :    temperature
c
c     The is currently no horizontal diffusivity in the set (set it to zero)
c     ------------------------------------------
      integer :: ix,iy,iz,j
      real    :: dz
      logical :: not_ok
      real    :: u_fill,v_fill,w_fill
c     ------------------------------------------ 
c
c     transfer raw data from cmod compact to mesh_grid 
c
      do ix=1,nx
         do iy=1,ny
            if (wetmask(ix,iy)==0) cycle ! nothing for dry points
            do iz=1,bottom_layer(ix,iy)
               j = m1_cmod(data_set_handler)%p(ny+1-iy, ix, iz)
               if (j==0) then
                  write(*,*) "grid inconsistency: ",ix,iy,iz
                  stop
               endif
               u(ix,iy,iz)        = u_cmod(data_set_handler)%p(j)
               v(ix,iy,iz)        = v_cmod(data_set_handler)%p(j)
               w(ix,iy,iz)        = w_cmod(data_set_handler)%p(j)
               temp(ix,iy,iz)     = t_cmod(data_set_handler)%p(j)
               salinity(ix,iy,iz) = s_cmod(data_set_handler)%p(j)
             enddo    
         enddo
      enddo
c
c     postprocess data, suncronize auxillary fields, pad holes
c      
c     ---- currently no turbulent diffusivity in data set: 
c          for now, this is read from input file in init_physical_fields
c          and vdiffus/hdiffus is set to this value
c
c          NB: horizontal derivatives NOT implemented
c          when spatially non constant hdiffus are applied,
c          interpolate_turbulence_deriv must be updated
c
c  

c
c     ------ flip sign of w 
c     ------ physical_fields sign convention is positive down, cmod convention is positive up
c
      w = -w
c     ------ add DSLM to auxillary grid descriptors for wet points  
c            acc_width(ix,iy,1) = 0 (sea surface)  
c            Jun 16, 2011: capture rare cases where dz creates negative depths

      do ix=1,nx
         do iy=1,ny
            if (wetmask(ix,iy)==0) cycle ! nothing for dry points
            j  = m1_cmod(data_set_handler)%p(ny+1-iy, ix, 1)
            if (j==0) then
               write(*,*) "grid inconsistency: ",ix,iy,iz
               stop
            else
               dz =  z_cmod(data_set_handler)%p(j)
            endif
            wdepth(ix,iy)       = max(0.0, wdepth0(ix,iy) - dz)
            acc_width(ix,iy,2:) = acc_width0(ix,iy,2:) - dz
            ccdepth(ix,iy,1)    = ccdepth0(ix,iy,1)    - dz*0.5 
            ccdepth(ix,iy,2:)   = ccdepth0(ix,iy,2:)   - dz     
         enddo
      enddo
    
      end subroutine syncronize_data




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
      call get_grid_coordinates(xyz,x,y,z)
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



      end module

      
