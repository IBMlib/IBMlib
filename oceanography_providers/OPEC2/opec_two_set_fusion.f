ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     OPEC cmod+ERGOM pbi (experimental)
c     ---------------------------------------------------
c     $Rev: 353 $
c     $LastChangedDate: 2011-06-16 23:29:14 +0200 (Thu, 16 Jun 2011) $
c     $LastChangedBy: asch $ 
c    
c     Data fusion version, where user define a custom grid that extends
c     a high resolution grid to larger domain, with data (and topography) interpolated from
c     embedding coarse grids outside the high resolution grid.
c     Mesh grid module is configured for this extended
c     high resolution grid, which is treated as an original source
c
c     v0: canonical hierachical interpolation
c     v1: There is a minor topographic issue: the NW corner of the fine grid extends to the
c         North Sea, but is treated dry for Eulerian reasons; in Lagrangian contexts, this
c         area then appears dry by the hierachical interpolation approach. This area is manually excluded
c         by the hierachical interpolation approach with the above_jutland_shoulder() query function
c      
c     currently omit wind stress as this is not used
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module physical_fields

c     plugin for reading native compressed cmod output (rename data, avoid mixing types)
c     test_undef: safe value > 0 to detect undef values by real comparisons, like: if (myval < unsetlim) < actions >
      
      use read_cmod_ergom, 
     +             unsetlim   => test_undef,
     +             m1_cmod    => m1,
     +             hz_cmod    => hz,
     +             cellw_cmod => cellw

      use regular_lonlat_grid  
      use regular_lonlat_data

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
      public :: interpolate_salty 
      public :: interpolate_wdepth
c
c      public :: interpolate_zooplankton
c      public :: interpolate_oxygen
c      public :: interpolate_nh4                       
c      public :: interpolate_no3                       
c      public :: interpolate_po4                        
c      public :: interpolate_diatoms                   
c      public :: interpolate_flagellates                 
c      public :: interpolate_cyanobacteria              
c      public :: interpolate_organic_detritus                   
c      public :: interpolate_part_org_matter 
c      public :: interpolate_DIC 
c      public :: interpolate_alkalinity  
c      public :: interpolate_DIN    
c      public :: interpolate_chlorophyl  
c      
      public :: is_wet    
      public :: is_land
      public :: horizontal_range_check
      public :: coast_line_intersection
      public :: get_pbi_version

c     -------------------- module data --------------------  
      
      real,parameter :: molecular_diffusivity = 1.e-9    ! unit m2/s
      real,parameter :: mmolN2kgDW = 14./1000/0.074/1000 ! conversion factor mmol N/m3 -> kg DW/m3
      real           :: hdiffus                          ! constant horizontal diffusivity (m2/s)
      real           :: vdiffus                          ! constant vertical diffusivity  (m2/s)
      real,parameter :: padval = 0                       ! module level placeholder for invalid queries
c
c     ------ data frame handler ------
   
      character*999      :: hydroDBpath       ! path for hydrographic data sets     
      logical            :: include_bio      = .true. ! handle to control whether biodata is read (physics is always read)

c     --- simulation domain < outer domain
      
      type(lonlat_grid) :: custom_grid     ! defines domain + topography, only horizontal 
      
c     --- native grid buffers ---

      integer, parameter         :: ngrids = 2      ! inner, outer
      integer                    :: set_handler(ngrids)
      type(lonlat_grid),target   :: grid(ngrids)    ! in descenting resolution order (first = prefer)
      type(lonlat_data3D)        :: u(ngrids) ! x currents; in descenting resolution order (first = prefer)
      type(lonlat_data3D)        :: v(ngrids) ! y currents; in descenting resolution order (first = prefer)
      type(lonlat_data3D)        :: w(ngrids) ! z currents; in descenting resolution order (first = prefer)
      type(lonlat_data2D)        :: z(ngrids) ! sea surface elevation;    NB: 2D not 3D
      type(lonlat_data3D)        :: s(ngrids) ! salinity; in descenting resolution order (first = prefer)
      type(lonlat_data3D)        :: t(ngrids) ! temperature; in descenting resolution order (first = prefer)
      
c     ===================================================
                            contains
c     ===================================================
  
      subroutine init_physical_fields(time)
c     ------------------------------------------
c     Do not trigger data load
c     ------------------------------------------
      type(clock), intent(in),optional :: time
c     ------------------------------------------
      if (present(time)) call set_master_clock(time)
      write(*,*) trim(get_pbi_version()) 

      call read_control_data(simulation_file,"hydroDBpath",hydroDBpath)
      write(*,*) "init_physical_fields: hydrographic database path =", 
     +           trim(adjustl(hydroDBpath))
    
      call setup_grids()        ! incl allocation of 3D arrays
c
c      call dump_dry_or_wet_points(custom_grid, 33)
c
c     ---- currently no horizontal turbulent diffusivity in data set ----
c          set it to molecular diffusivity lower limit, if no values 
c          are provided
c
c     resolved user provided constant horizontal_diffusivity (optional)  
c    
      if (count_tags(simulation_file, "horizontal_diffusivity")/=0) then
         call read_control_data(simulation_file,
     +                          "horizontal_diffusivity",hdiffus)
         hdiffus = max(molecular_diffusivity, hdiffus) ! never exceed lover limit
         write(*,563) "horizontal_diffusivity",hdiffus
      else
         hdiffus = molecular_diffusivity
         write(*,564) "horizontal_diffusivity",hdiffus
      endif
c 
c     resolved user provided constant vertical_diffusivity (optional)  
c    
      if (count_tags(simulation_file, "vertical_diffusivity")/=0) then
         call read_control_data(simulation_file,
     +        "vertical_diffusivity",vdiffus)
         vdiffus = max(molecular_diffusivity, vdiffus) ! never exceed lover limit
         write(*,563) "vertical_diffusivity",vdiffus
      else
         vdiffus = molecular_diffusivity
         write(*,564) "vertical_diffusivity",vdiffus
      endif
c

 563  format("init_physical_fields: read const  ",a,"=", e12.5," m2/s")
 564  format("init_physical_fields: using const ",a,"=", e12.5," m2/s")

      end subroutine init_physical_fields

 

      character*100 function get_pbi_version()  
      get_pbi_version =  "OPEC (grid fusion)"//
     +                   "pbi version: $Rev: 353 $"
      end function

   

      subroutine setup_grids()
c     ---------------------------------------------------
c     Grid definitions is $(hydroDBpath)/cfg.nml (which points to auxillary
c     grid descriptors).
c     Allocate grid arrays in this subroutine. 
c     grid point (ix,iy) = (1,1) is at (lambda1,phi1)
c     ---------------------------------------------------
      integer              :: nx, ny, mx, my, mz ! LOCAL DUMMIES
      integer              :: ix, iy, iz, jx, jy, ig, iset
      integer              :: idx1d, ibot,ixsw,iysw,ixne,iyne
      integer              :: ixi,iyi
      real                 :: lon1, lat1, dlon, dlat ! LOCAL DUMMIES
      real                 :: lambda1,phi1,dlambda,dphi
      real                 :: lonmax_outer,lonmax_custom
      real                 :: latmax_outer,latmax_custom
      real                 :: bbx(4), geo(2), xy(2), x, y
      type(lonlat_grid), pointer :: inner_grid, outer_grid, g
      character*99         :: set_id
      logical              :: inside_inner
c     ---------------------------------------------------
c
c     --- first resolve whether to read biogeochemistry: optional tag include_biogeochem
c   
      if (count_tags(simulation_file, "include_biogeochem") /= 0) then
         call read_control_data(simulation_file,"include_biogeochem",
     +                       include_bio)
      else 
         include_bio = .true. ! set default
      endif

c     ---- remove this plug when biogeochemistry is implemented ---
      if (include_bio) then
         write(*,*) "include_bio = True currently not implemented"
         stop 622
      endif
      
      call init_read_cmod_ergom(hydroDBpath, include_bio) ! defined in module read_cmod

      
c     ======= read custom configuration  =======
c
c     inner_grid: embedded high resolution grid    -> inner_set_handler
c     custom_grid:  horizontal simulation zone (synthesized high resolution 2D topographic description)      
c     outer_grid: embedding coarse resolution grid -> outer_set_handler
c
c
      call read_control_data(simulation_file,"inner_domain", set_id)        ! string identifier
      write(*,*) "Inner cmod set: ", trim(adjustl(set_id))

      call get_grid_descriptors(set_id, set_handler(1),
     +     lon1, dlon, lat1, dlat, mx, my, mz) 

      call init_lonlat_grid(grid(1),mx,my,lon1,lat1,dlon,dlat,mz)
      inner_grid => grid(1)
      
      call read_control_data(simulation_file,"outer_domain",set_id)        ! string identifier
    
      write(*,*) "Outer cmod set: ", trim(adjustl(set_id))

      call get_grid_descriptors(set_id, set_handler(2),
     +     lon1, dlon, lat1, dlat, mx, my, mz)  
      
      call init_lonlat_grid(grid(2),mx,my,lon1,lat1,dlon,dlat,mz)
      outer_grid => grid(2)
      
      call read_control_data(simulation_file,"custom_domain", bbx) ! lonmin,latmin,lonmax,latmax
      if ((bbx(1)>bbx(3)).or.(bbx(2)>bbx(4))) then
         write(*,*) "invalid custom_domain: ", bbx
         stop 654
      endif
      
c     --- compute and syncronize custom_domain (2D) ---

      dlambda = inner_grid%dlon ! inherit inner longitude step
      dphi    = inner_grid%dlat ! inherit inner latitude step
      call get_horiz_ncc_index(inner_grid,bbx(1:2),ix,iy)  ! custom_domain SW cell indices on inner_grid
      call get_horiz_ncc_index(inner_grid,bbx(3:4),jx,jy)  ! custom_domain NE cell indices on inner_grid
      call get_horiz_geo_coordinates(inner_grid, 1.0*ix, 1.0*iy, geo) ! custom SW corner shifted to inner grid
      lambda1 = geo(1)
      phi1    = geo(2)
      nx      = 1 + jx-ix                      ! custom longitude dimension
      ny      = 1 + jy-iy                      ! custom latitude dimension

c     --- check  (custom_domain < outer_domain) ---
      
      if (lambda1 < outer_grid%lon1) then
         write(*,*) "custom_domain exceeds outer domain (x-)"
         stop 115
      endif
      lonmax_outer = outer_grid%lon1 + 
     +               (outer_grid%nx-1) * outer_grid%dlon
      lonmax_custom = lambda1 + (nx-1)*dlambda
      if (lonmax_custom > lonmax_outer) then
         write(*,*) "custom_domain exceeds outer domain (x+)"
         stop 116
      endif
      
      if (phi1 < outer_grid%lat1) then
         write(*,*) "custom_domain exceeds outer domain (y-)"
         stop 117
      endif
      latmax_outer = outer_grid%lat1 +
     +               (outer_grid%ny-1) * outer_grid%dlat
      latmax_custom = phi1 + (ny-1)*dphi
      if (latmax_custom > latmax_outer) then
         write(*,*) "custom_domain exceeds outer domain (y+)"
         stop 118
      endif

      call init_lonlat_grid(custom_grid, nx, ny, 
     +     lambda1, phi1, dlambda, dphi)   !   2D grid initialized
      
      write(*,231) nx, ny
      write(*,232) lambda1, dlambda
      write(*,233) phi1,    dphi 
   
 231  format("setup_grids: custom grid: 3d grid dim (nx,ny) = ",
     + i4,i4)
 232  format("setup_grids: custom grid: lambda1 = ",f12.7,
     +     " dlambda = ",f12.7," degrees")  
 233  format("setup_grids: custom grid: phi1    = ",f12.7,
     +     " dphi    = ",f12.7," degrees")
      
c     ================ initialize native data buffers ================ 

c     ---------------- physics ----------------
      do ig=1,ngrids
        call init_lonlat_data(u(ig), grid(ig))
        call init_lonlat_data(v(ig), grid(ig))
        call init_lonlat_data(w(ig), grid(ig))
        call init_lonlat_data(z(ig), grid(ig))    ! 2D not 3D
        call init_lonlat_data(s(ig), grid(ig))
        call init_lonlat_data(t(ig), grid(ig))
      enddo
      
      
c     ---------------- biogeochemistry ----------------            
      
      
c     ================ initialize topographic arrays  ================
c 
c     set auxillary grid arrays: ccdepth0, acc_width0, bottom_layer, wetmask,wdepth0
c     
      do ig=1,ngrids            ! inner, outer
         iset             =  set_handler(ig)
         g                => grid(ig)
         g%wdepth0        = -1.0   ! default dry
         g%bottom_layer   =  0     ! default dry
c        do not assign a value to ccdepth0/acc_width0 to be able to spot unknown dry point violations

         g%acc_width0(:,:,1)  = 0.   ! sea surface (DSLM=0) 
         g%acc_width0(:,:,2:) = 1.e6 ! acc_width0 must be inceasing
         do ix=1,g%nx
            do iy=1,g%ny
               if (m1_cmod(iset)%p(g%ny+1-iy,ix,1)>0) then  ! wet/dry longitunal point test  
c             
c              Vertical loop for wet points. At loop exit iz will point to first dry layer
c              F90 std condition for full do-loop termination is iz=g%nz+1
c            
                  do iz=1,g%nz 
                     idx1d = m1_cmod(iset)%p(g%ny+1-iy,ix,iz)
                     if (idx1d > 0) then ! wet/dry signal for cmod   
                        g%ccdepth0(ix,iy,iz) = hz_cmod(iset)%p(idx1d)
                        g%acc_width0(ix,iy,iz+1) = 
     +                       g%acc_width0(ix,iy,iz) + 
     +                       cellw_cmod(iset)%p(idx1d)
                     else
                        exit 
                     endif   
                  enddo         ! iz loop

                  ibot                         = iz-1 ! last wet layer
                  g%bottom_layer(ix,iy) = ibot
                  g%wdepth0(ix,iy)      =
     +                 g%acc_width0(ix,iy,1+ibot)                    
c              

               endif            ! wet/dry longitunal point
            enddo               ! iy loop
         enddo                  ! ix loop

         write(*,*) "read_grid_desc: transferred cmod bathymetry", ig    
c
c     --- define the wetmask statically (no flooding/drying) --- 
c   
         where(g%wdepth0 > 0.0)
            g%wetmask = 1         ! wet
         elsewhere
            g%wetmask = 0         ! dry
         end where
         
      enddo                       ! ig loop

c
c === set custom_grid%wetmask used for global topography 
c
c     ---- compute window of inner_grid in custom_grid : (ixsw:ixne, iysw:iyne) ----
c          inner_grid is coherent with custom_grid        
c
      geo(1) = inner_grid%lon1
      geo(2) = inner_grid%lat1
      call get_horiz_ncc_index(custom_grid,geo,ixsw,iysw)   
      ixne = ixsw + inner_grid%nx - 1
      iyne = iysw + inner_grid%ny - 1
c      
c     ---- map from wetmasks of inner_grid / outer_grid ----
c          but discard NW corner of fine grid (ig==1)
      custom_grid%wetmask = 1 
      do ix=1,nx
         do iy=1,ny
            inside_inner = (ix.ge.ixsw).and.(ix.le.ixne).and.
     +           (iy.ge.iysw).and.(iy.le.iyne)
            xy(1) = 1.0*ix
            xy(2) = 1.0*iy
            call get_horiz_geo_coordinates(custom_grid,xy,geo)
            if (above_jutland_shoulder(geo)) inside_inner = .false.  ! manually discard NW corner of fine grid
            if (inside_inner) then
               ixi = ix-ixsw+1  ! map ix to this index on inner grid
               iyi = iy-iysw+1  ! map iy to this index on inner grid
               custom_grid%wetmask(ix,iy) = inner_grid%wetmask(ixi,iyi)
            else                ! cast outer_grid topography onto custom_grid
               if (is_land_for_grid(outer_grid, geo))
     +              custom_grid%wetmask(ix,iy) = 0
            endif  ! inside_inner
         enddo     ! iy = ..
      enddo        ! ix = ..
c     ------------------------------------------       
      end subroutine setup_grids



      subroutine close_physical_fields()
c     ------------------------------------------
      integer :: ig
c     ------------------------------------------
      do ig = 1,ngrids
         call finalize_lonlat_data(u(ig))
         call finalize_lonlat_data(v(ig))
         call finalize_lonlat_data(w(ig))
         call finalize_lonlat_data(z(ig))
         call finalize_lonlat_data(s(ig))
         call finalize_lonlat_data(t(ig))
         call finalize_lonlat_grid(grid(ig)) ! last action for this grid
      enddo    
      call close_read_cmod_ergom()   
c     ------------------------------------------
      end subroutine close_physical_fields


      
      subroutine update_physical_fields(time, dt)
c     -----------------------------------------------------------  
c     File mapping / frame buffering delegated to read_cmod_ergom.update_buffers
c     Only invoke syncronize_data, if new data is loaded
c     -----------------------------------------------------------  
      type(clock), intent(in),optional :: time
      integer, intent(in),optional     :: dt
      type(clock), pointer             :: aclock
      logical                          :: newdata
c     ------------------------------------------
      aclock => get_master_clock()
      if (present(time)) then
         call set_master_clock(time)
      elseif (present(dt)) then
         call add_seconds_to_clock(aclock, dt)
      endif     
c
      newdata = update_buffers(aclock)
      if (newdata) call syncronize_data() !
c     ------------------------------------------       
      end subroutine update_physical_fields




      subroutine syncronize_data()
c     ------------------------------------------
c     Syncronize auxillary grid fields after data load (recent load flagged by update_buffers()). 
c     It is assumed that grids do not
c     change during a simulation (assertion not checked)
c  
c     The is currently no horizontal diffusivity in the set (set it to zero)
c     ------------------------------------------
      integer :: ix,iy,iz
      real    :: dz
      integer :: iset, ig
      type(lonlat_grid), pointer :: g
c     ------------------------------------------ 
c     
c     arrays are padded with read_cmod_ergom.dryval for dry points
c     by get_X()
c
      do ig=1,ngrids
         
         iset =  set_handler(ig)
         g                => grid(ig)
         
         call get_u(u(ig)%buf, iset)
         call get_v(v(ig)%buf, iset)
         call get_w(w(ig)%buf, iset)
         call get_z(z(ig)%buf, iset)
         call get_s(s(ig)%buf, iset)
         call get_t(t(ig)%buf, iset)

c      if (include_bio) then
c         call get_eco_3D(oxygen, iset,                "oxy")
c         call get_eco_3D(zoo, iset,                   "zoo") ! sum of zooplankton classes
c         call get_eco_3D(nh4, iset,                   "nh4") 
c         call get_eco_3D(no3, iset,                   "no3")                      
c         call get_eco_3D(po4, iset,                   "po4")                      
c         call get_eco_3D(diatoms, iset,               "dia")            
c         call get_eco_3D(flagellates, iset,           "fla") 
c         call get_eco_3D(cyanobacteria, iset,         "cya")         
c         call get_eco_3D(organic_detritus, iset,      "odt")                
c         call get_eco_3D(part_org_matter, iset,       "pom") 
c         call get_eco_3D(dissolv_inorg_carbon, iset,  "dic")    
c         call get_eco_3D(alkalinity, iset,            "alk")
c         call get_eco_3D(dissolv_inorg_nitrogen, iset,"din") 
c         call get_eco_3D(chlorophyl, iset,            "chl") 
c      endif 
c
c     avoid tricky current interpolations where the number of wet layers change
c
         where (u(ig)%buf > unsetlim) 
            u(ig)%buf=0
         end where
         where (v(ig)%buf > unsetlim) 
            v(ig)%buf=0
         end where
         where (w(ig)%buf > unsetlim) 
            w(ig)%buf=0
         end where    
c
c     ==== postprocess data, suncronize auxillary fields, pad holes, change units ====
c      
c     ---- currently no turbulent diffusivity in data set: 
c          for now, this is read from input file in init_physical_fields
c          and vdiffus/hdiffus is set to this value
c
c          NB: horizontal derivatives NOT implemented
c          when spatially non constant hdiffus are applied,
c          interpolate_turbulence_deriv must be updated
c
c     ------ flip sign of w:  physical_fields sign convention is positive down, cmod convention is positive up
c
         where (w(ig)%buf < unsetlim) 
            w(ig)%buf = -w(ig)%buf
         end where
c
c     ------ convert zooplankton concentration from unit mmol N/m3 -> kg DW/m3
c
c      if (include_bio) then ! zoo only allocated if true
c         where (zoo < unsetlim) 
c            zoo = zoo*mmolN2kgDW
c         end where 
c      endif
c
c     ------ add DSLM to auxillary grid descriptors for wet points and 
c            update revision tags
         
         call update_topographic_arrays(g, z(ig)%buf) ! includes revision tag update
         call update_revtag_data(u(ig), g%revtag)     ! sync to grid revision tag
         call update_revtag_data(v(ig), g%revtag)     ! sync to grid revision tag
         call update_revtag_data(w(ig), g%revtag)     ! sync to grid revision tag
         call update_revtag_data(z(ig), g%revtag)     ! sync to grid revision tag
         call update_revtag_data(s(ig), g%revtag)     ! sync to grid revision tag
         call update_revtag_data(t(ig), g%revtag)     ! sync to grid revision tag

      enddo  ! ig=1,ngrids   
         
      end subroutine syncronize_data


      
      logical function above_jutland_shoulder(geo)
c     --------------------------------------------------------------------
c     Test whether horizontal position geo is above the Jutland shoulder,
c     defined by the line from pS(2) to pN(2)
c     For manual fixing of topographic issue, see top of module
c     --------------------------------------------------------------------
      real, intent(in) :: geo(:)
      real, parameter :: pN(2) = (/10.4167, 57.6/)
      real, parameter :: pS(2) = (/10.2500, 57.5/)
      real, parameter :: orthovec(2) = (/pS(2)-pN(2),  pN(1)-pS(1)/) ! (-y,x)
      real            :: dxy(2), project
c     --------------------------------------------------------------------
      dxy = geo(1:2)-pS
      if (sum(dxy*orthovec)>0) then
         above_jutland_shoulder = .true.
      else
         above_jutland_shoulder = .false.
      endif
      end function above_jutland_shoulder


      
      subroutine select_interpolation_grid(xyz, ig, status)
c     ------------------------------------------------------------
c     Pick first matching grid, since grids are ordered by decreasing resolution,
c     i.e. interpolation is preferred on grid with highest local resolution
c     ------------------------------------------------------------
      real, intent(in)     :: xyz(:) ! interpolation point
      integer, intent(out) :: ig     ! preferred grid
      integer, intent(out) :: status ! interpolation flag
      logical              :: inside
c     ------------------------------------------------------------     
      ig = -1  ! invalid
      inside = horizontal_range_check_for_grid(custom_grid, xyz)
      if (.not.inside) then
         status = err_reglonlat_horiz_violation             ! signal horizontal violation
      else
         status = err_reglonlat_success 
         do ig = 1,ngrids
            if ((ig==1).and.above_jutland_shoulder(xyz)) cycle  ! manually skip NW corner of fine grid (ig==1)
            if (horizontal_range_check_for_grid(grid(ig), xyz)) return
         enddo
         write(*,*) "select_interpolation_grid: internal error" ! we should not end here
         stop 292
      endif
      end subroutine select_interpolation_grid

      

      
      subroutine interpolate_3Dbuffer_set(buf, xyz, result, status)
c     ------------------------------------------------------------
c     prototype value interpolator for standard-ordered buffer set      
c     -----------------------------------------------------------
      type(lonlat_data3D), intent(in) :: buf(:)  ! interpolation buffers
      real, intent(in)                :: xyz(:)  ! interpolation point
      real, intent(out)               :: result  ! interpolation result
      integer, intent(out)            :: status  ! interpolation flag
      integer :: ig
c     -----------------------------------------------------------      
      call select_interpolation_grid(xyz, ig, status)
      if (status == err_reglonlat_success) then ! interpolate value at xyz in grid ig 
         call interpolate_cc_data(buf(ig), xyz, 0, padval,
     +                            result, status)   ! pass on status from interpolate_cc_data3D
      else  ! keep status > 0
         result = padval
      endif
      end subroutine interpolate_3Dbuffer_set

      
      
      subroutine interpolate_turbulence(xyz, r3, status)
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
      status = err_reglonlat_success ! skip range check 
      r3(1:2) = hdiffus
      r3(3)   = vdiffus
      end subroutine 

      
      subroutine interpolate_turbulence_deriv(xyz, r3, status)
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
      status = err_reglonlat_success ! skip range check
      r3 = 0.0
      end subroutine 

      
      subroutine interpolate_temp (xyz, r, status) 
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
      call interpolate_3Dbuffer_set(t, xyz, r, status)
      end subroutine 

      
      subroutine interpolate_salty(xyz, r, status)
      real, intent(in)     :: xyz(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
      call interpolate_3Dbuffer_set(s, xyz, r, status)
      end subroutine 


      subroutine interpolate_wdepth(xy, result, status)
c     -----------------------------------------------------------
c     Pass depth interpolation onto best grid
c     gate keeper function for is_wet
c     -----------------------------------------------------------            
      real, intent(in)     :: xy(:)
      real, intent(out)    :: result
      integer, intent(out) :: status
      integer :: ig
      logical :: inside
c     -----------------------------------------------------------
      result = padval
      inside = horizontal_range_check_for_grid(custom_grid, xy)
      if (.not.inside) then
         status = err_reglonlat_horiz_violation       ! signal horizontal violation
      elseif (is_land_for_grid(custom_grid, xy)) then 
         status = err_reglonlat_dry_point             ! signal dry point
      else  ! inside and wet
         call select_interpolation_grid(xy, ig, status)
         call interpolate_wdepth_for_grid(grid(ig), xy, result, status)  ! determines final status
      endif
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
      real, parameter  :: htol = 1.e-6 ! tolerance for surface/bottom test
c     ------------------------------------------ 
      call interpolate_wdepth(geo, wd, status)

      select case (status)
      case (err_reglonlat_success)  ! interior wet point
         if ((geo(3)>-htol).and.(geo(3)<wd+htol)) then
            is_wet = .true.
         else
            is_wet = .false.
         endif
      case (err_reglonlat_horiz_violation)    ! horizontal range violation, set result = padval
         is_wet = .true.
      case (err_reglonlat_dry_point)          ! dry horizontal point
         is_wet = .false. ! logical correct
      case (err_reglonlat_runtimeerr)         ! auxillary arrays for query not initialized
         write(*,'(a,a)') "is_wet:",
     +        err_msg_reglonlat(err_reglonlat_runtimeerr)
         stop 328
      case default ! other internal problem, should not end here
         write(*,'(a)') "is_wet: unhandled error"
         stop 329
      end select
c     ------------------------------------------
      end function


      
      LOGICAL function is_land(xy)
c     ---------------------------------------------------------------
c     pass query onto custom_grid
c     --------------------------------------------------------------           
      real, intent(in)     :: xy(:)
      is_land = is_land_for_grid(custom_grid, xy)
      end function

      
      LOGICAL function horizontal_range_check(xy)
c     ---------------------------------------------------------------
c     pass query onto custom_grid
c     --------------------------------------------------------------      
      real, intent(in)     :: xy(:)
      horizontal_range_check =
     +     horizontal_range_check_for_grid(custom_grid, xy)
      end function

      

      subroutine coast_line_intersection(geo1, geo2, cross_coast,
     +     georef, geohit)
c     --------------------------------------------------------------
c     pass query onto custom_grid
c     --------------------------------------------------------------
      real, intent(in)              :: geo1(:),geo2(:)
      logical, intent(out)          :: cross_coast
      real, intent(out)             :: georef(:), geohit(:)
c     --------------------------------------------------------------
      call coastline_intersection_for_grid(custom_grid,
     +     geo1, geo2, cross_coast,
     +     georef, geohit)
      
      end subroutine coast_line_intersection
      

            
      subroutine interpolate_currents(xyz, uvw, status)
c     ------------------------------------------
c     Perform linear two-face interpolation in current fields into vector uvw, scanning
c     a set of nested grids for optimal resolution. Scanning procedure is consistent with
c     select_interpolation_grid()
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
c     definition derived from the grid. Notice that staggering implies that in rim cases
c     currents may be interpolated from a grid different than that used for cell-centered data
c     Project these rim points onto the interior domain, if last grid is selected
c     and do not flag them as domain violation
c     ------------------------------------------ 
      real, intent(in)     :: xyz(:)  ! (lon,lat,depth)
      real, intent(out)    :: uvw(:)
      integer, intent(out) :: status

      integer              :: statu, statv, statw
      integer              :: ix,iy,iz,jwest,jnorth,jlow,ibot,ig
      real                 :: x,y,z, sx,sy,sz
      logical              :: inside
      type(lonlat_grid), pointer :: g
      real, pointer              :: uu(:,:,:), vv(:,:,:), ww(:,:,:)
c     -----------------------------------------------------------
c     custom_grid defines overall valid domain 
      inside = horizontal_range_check_for_grid(custom_grid, xyz)
      if (.not.inside) then
         uvw    = 0.
         status = err_reglonlat_horiz_violation  ! signal range violation
         return
      endif
      
c.... identify grid (ig) in sequence that should be used for current interpolation at xyz
c     if we end at last grid (ig==ngrids), apply this for extrapolation or padding, if needed,
c     when overall horizontal_range_check for custom_grid above has been passed
c
      do ig = 1,ngrids
         g  => grid(ig)
         if ((ig==1).and.above_jutland_shoulder(xyz)) cycle  ! manually skip NW corner of fine grid (ig==1)
         call get_grid_coordinates(g,xyz,x,y,z)
         if (((x<1.5).or.(x>(g%nx+0.5))).or.    ! x range violation for this grid
     +       ((y<0.5).or.(y>(g%ny-0.5))).or.    ! y range violation for this grid
     +       ((z<0.5).or.(z>(g%nz+0.5)))) then  ! z range violation for this grid
            continue              ! not interior for this grid
         else
            exit                  ! interior for this grid, keep ig
         endif
      enddo            ! do ig=1,ngrids
      ig = min(ig,ngrids) ! ig==ngrids+1 at exit, if unsuccesful
c      
      g  => grid(ig)      ! create short hand
      uu => u(ig)%buf     ! create short hand
      vv => v(ig)%buf     ! create short hand 
      ww => w(ig)%buf     ! create short hand
       
      if (.not.is_wet_for_grid(g, xyz)) then
         uvw    = 0.
         status = err_reglonlat_dry_point      ! signal dry point
         return
      endif
      status = err_reglonlat_success             ! signal all OK, after range/wet checks passed


c.....determine cell associations of point, constrained to 1 <= i <= n
      ix   = max(1, min(nint(x), g%nx)) ! u at cell east  face
      iy   = max(1, min(nint(y), g%ny)) ! v at cell south face
      ibot = g%bottom_layer(ix,iy)      ! 1 <= ibot <= nz
      iz   = max(1, min(nint(z),ibot))  ! w at cell upper face

c.....determine relative intracell coorninate s, constrained to 0 < s < 1 
c     when point exceed coordinate bounds s is assigned 0 or 1
      sx = max(0.d0, min(1.d0, x+0.5d0-real(ix)))  
      sy = max(0.d0, min(1.d0, y+0.5d0-real(iy)))           
      sz = max(0.d0, min(1.d0, z+0.5d0-real(iz)))  

c.....face-to-face interpolation 
c     cell indices satisfy 1 <= (ix,iy,iz) <= (nx,ny,ibot)
      jwest  = max(1, min(nint(x-1), g%nx)) ! cell west  face
      jnorth = max(1, min(nint(y+1), g%ny)) ! cell north fce
      uvw(1) = sx*uu(ix, iy, iz)    + (1.-sx)*uu(jwest, iy, iz)
      uvw(2) = sy*vv(ix, jnorth,iz) + (1.-sy)*vv(ix, iy, iz)

      if (iz==ibot) then
         jlow = 0               ! for debugging
         uvw(3) = (1.-sz)*ww(ix, iy, ibot) ! w=0 at sea bed, ibot corresponds to upper cell face
      else
         jlow = max(1, min(nint(z+1),ibot)) ! cell lower face
         uvw(3) = sz*ww(ix, iy, jlow) + (1.-sz)*ww(ix, iy, iz)
      endif
      
c.....deassociate pointers      
      
      nullify( g  )
      nullify( uu )
      nullify( vv )
      nullify( ww )
      
c     ------------------------------------------ 
      end subroutine interpolate_currents



      end module

      
