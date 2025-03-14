      module mesh_grid
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This module represents data on a 3D array where the vertical
c     mesh is perpendicular to the horizontal and the
c     horizontal grid are a mesh (without further restrictions)
c     This includes sigma and z-grid with open/closed surface
c     Serves both the physical and biogeochemical interface.
c     By default, only the physical interface is initialized by init_mesh_grid 
c
c     $Rev: $
c     $LastChangedDate:  $
c     $LastChangedBy: $ 
c            
c     This module hosts data arrays, add vertical issues to horizontaal to get 3D
c     and proivide 2D/3D interpolation
c
c     * provide interpolate_X
c     * provide is_wet
c     * host nz
c     * host 3D data arrays
c     
c     LOG: 
c     * stripped out time components (deal only with space/data) ASC Feb 16, 2011
c     * made biogeochemical components optional                  ASC May 2013
c     * added optional data revision logging, to allow plugins 
c       to minimize update overhead                              ASC Jun 2013
c
c     TODO:
c       * introduce an "only" list for fields (applied at init time), to reduce memory usage
c         as number of water quality parameters is steadily growing
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use horizontal_representation       ! injects nx,ny
      use horizontal_grid_transformations 
   
      use constants
      use array_tools                 ! search_sorted

      implicit none
      private           ! default scope

      public :: init_mesh_grid   ! init  this module
      public :: close_mesh_grid  ! close this module
      
      public :: interpolate_turbulence
      public :: interpolate_turbulence_deriv
      public :: interpolate_currents
      public :: interpolate_temp
      public :: interpolate_salty  
      public :: interpolate_wind_stress    
      public :: interpolate_wdepth
c
      public :: interpolate_zooplankton
      public :: interpolate_oxygen
      public :: interpolate_nh4                       
      public :: interpolate_no3                       
      public :: interpolate_po4                        
      public :: interpolate_diatoms                   
      public :: interpolate_flagellates                 
      public :: interpolate_cyanobacteria              
      public :: interpolate_organic_detritus                   
      public :: interpolate_part_org_matter 
      public :: interpolate_DIC 
      public :: interpolate_alkalinity  
      public :: interpolate_DIN    
      public :: interpolate_chlorophyl   
c
      public :: is_wet    
      public :: is_land                 ! reexport from horizontal_representation
      public :: horizontal_range_check  ! reexport from horizontal_grid_transformations
      public :: coast_line_intersection ! reexport from horizontal_representation
    
c     ---------- data revision monitoring ----------
c
c     allows independent services to defer updating data dependent on data 
c     in mesh_grid until it is required. Currently revision is a integer*16 (data_revision_tag) 
c     making it extremely unlikely to accidentially miss an update (due to the fact that
c     revision cycles after a finite number of updates
c  
      public :: get_data_rev_tag_mesh_grid    ! poll current data revision
      public :: update_data_rev_tag_mesh_grid ! update current data revision (after modifications)
      
c     ---------- other exports ----------

      public :: get_grid_coordinates       ! formerly named get_ncc_coordinates
      public :: interpolate_cc_3Dgrid_data ! generic interpolator for 3D CC data
      public :: interpolate_cc_2Dgrid_data ! generic interpolator for 2D CC data

c     ---------- diagnostic exports ----------
      
      public :: write_topomask
      
c     -------------------- module data --------------------  

      integer, parameter :: verbose = 0  ! debugging output control
      real, parameter    :: htol = 1.e-6 ! tolerance for surface/bottom
      real, parameter    :: default_padding = 0.0  ! used, if not specified by user
      integer(kind=selected_int_kind(16)) :: data_revision_tag = -1 ! allow plugins to defer updates (-1 = unset data)
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
      real,allocatable,public :: salinity(:,:,:)   ! Salinity. [psu]
      real,allocatable,public :: vdiffus(:,:,:)    ! vertical   diffusivity [m**2/s]              
      real,allocatable,public :: hdiffus(:,:,:)    ! horizontal diffusivity [m**2/s]
      real,allocatable,public :: dslm(:,:)         ! current sea surface elevation over reference [m] (positive up)
      real,allocatable,target,public :: ccdepth(:,:,:)    ! cell center depth water below surface [m]; pointer 
      real,allocatable,target,public :: acc_width(:,:,:)  ! accumulated water above this layer [m] dim=nz+1  

c     ---- optional biogeochemistry
      logical                 :: incl_biogeochem = .false.     ! default: inactive
      real,allocatable,public :: zoo(:,:,:)                    ! Zooplankton [kg DW/m3]
      real,allocatable,public :: oxygen(:,:,:)                 ! O2 concentration [mmol/m3]
      real,allocatable,public :: nh4(:,:,:)                    ! [mmol/m3]          
      real,allocatable,public :: no3(:,:,:)                    ! [mmol/m3]                
      real,allocatable,public :: po4(:,:,:)                    ! [mmol/m3]                 
      real,allocatable,public :: diatoms(:,:,:)                ! [mmol/m3]                
      real,allocatable,public :: flagellates(:,:,:)            ! [mmol/m3]                 
      real,allocatable,public :: cyanobacteria(:,:,:)          ! [mmol/m3]                 
      real,allocatable,public :: organic_detritus(:,:,:)       ! [mmol/m3]                         
      real,allocatable,public :: part_org_matter(:,:,:)        ! particular organic matter [mmol/m3]       
      real,allocatable,public :: dissolv_inorg_carbon(:,:,:)   ! dissolved inorganic carbon [mmol/m3]       
      real,allocatable,public :: alkalinity(:,:,:)             ! [mmol/m3]       
      real,allocatable,public :: dissolv_inorg_nitrogen(:,:,:) ! dissolved inorganic nitrogen [mmol/m3]          
      real,allocatable,public :: chlorophyl(:,:,:)             ! [mmol/m3]       

      real   :: padval_u        = default_padding
      real   :: padval_v        = default_padding          
      real   :: padval_w        = default_padding         
      real   :: padval_temp     = default_padding    
      real   :: padval_salinity = default_padding
      real   :: padval_vdiffus  = default_padding             
      real   :: padval_hdiffus  = default_padding 
      real   :: padval_dslm     = default_padding
      real   :: padval_uwind    = default_padding
      real   :: padval_vwind    = default_padding
      real   :: padval_bgc      = default_padding      ! apply to all biogeochemistry

c     --- 2D grids ---
      
      real,allocatable,public     :: wdepth(:,:)      ! current depth at cell-center, including dslm [m]. wdepth = 0 at dry points
      integer, allocatable,public :: bottom_layer(:,:) ! last wet layer (0 for dry points) nx,ny
      real,allocatable,public     :: u_wind_stress(:,:) ! u of wind stress[N/m2] (positive east)
      real,allocatable,public     :: v_wind_stress(:,:) !  v of wind stress[N/m2] (positive north)   
      public                      :: wetmask           ! hosted by horizontal_representation

      real   :: padval_wdepth   = default_padding 

c     ===================================================
                            contains
c     ===================================================

                      
      subroutine init_mesh_grid(padding,
     +                     pad_u, pad_v, pad_w, pad_temp, pad_salinity, 
     +                     pad_vdiffus, pad_hdiffus, pad_dslm, pad_bgc, 
     +                     pad_wdepth, init_biogeochem)
c     ------------------------------------------------------
c     Assumes (nx,ny,nz) has been set by client module
c     Notice that module horizontal_grid_transformations must
c     be initialized directly from the top level physical_fields
c     by calling init_horiz_grid_transf(<args>), since the
c     specific arguments args depends directly on the grid type
c     The top level physical_fields should also directly import this module
c     By default (init_biogeochem == .false.), only the physical interface is 
c     initialized by init_mesh_grid 
c     
c     Pad value setup for grid interpolations at invalid points:
c       1) If no pad values are provided, the compiled-in value default_padding 
c       2) If the overall pad value padding is specified, this will overwrite 1)
c       3) If specific pad values are provided, they will overwrite 1) and 2)
c     
c     Allocate vdiffus(nx,ny,nz+1) to handle that many data sets define
c     data points at vertical faces (i.e. not cell-centered data and interpolation)
c     ------------------------------------------------------
      real, intent(in), optional :: padding ! overall value
      real, intent(in), optional :: pad_u 
      real, intent(in), optional :: pad_v            
      real, intent(in), optional :: pad_w          
      real, intent(in), optional :: pad_temp     
      real, intent(in), optional :: pad_salinity 
      real, intent(in), optional :: pad_vdiffus               
      real, intent(in), optional :: pad_hdiffus  
      real, intent(in), optional :: pad_dslm    
      real, intent(in), optional :: pad_bgc
      real, intent(in), optional :: pad_wdepth
      logical, intent(in), optional :: init_biogeochem
c     ------------------------------------------------------
      write(*,*) "init_mesh_grid: allocate grid arrays: begin" 

c     ----- physics -----
      allocate( u(nx,ny,nz)       )   
      allocate( v(nx,ny,nz)       )   
      allocate( w(nx,ny,nz)       )   
      allocate( temp(nx,ny,nz)    )  
      allocate( salinity(nx,ny,nz))  
      allocate( vdiffus(nx,ny,nz+1) ) ! enable data points at vertical faces, default is cell-centered data
      allocate( hdiffus(nx,ny,nz) ) 
      allocate( ccdepth(nx,ny,nz) )  
      allocate( acc_width(nx,ny,nz+1)   )

c     ----- biogeochemistry -----
      if (present(init_biogeochem)) incl_biogeochem = init_biogeochem ! else compile time value apply
      if (incl_biogeochem) then
         write(*,*) "init_mesh_grid: allocating optional biogeochem" 
         allocate( zoo(nx,ny,nz)     ) 
         allocate( oxygen(nx,ny,nz)  ) 
         allocate( nh4(nx,ny,nz)  )     
         allocate( no3(nx,ny,nz)  )       
         allocate( po4(nx,ny,nz)  )                
         allocate( diatoms(nx,ny,nz)  )           
         allocate( flagellates(nx,ny,nz)  )             
         allocate( cyanobacteria(nx,ny,nz)  )              
         allocate( organic_detritus(nx,ny,nz)  )                        
         allocate( part_org_matter(nx,ny,nz)  ) 
         allocate( dissolv_inorg_carbon (nx,ny,nz)  )      
         allocate( alkalinity(nx,ny,nz)  ) 
         allocate( dissolv_inorg_nitrogen(nx,ny,nz)  )   
         allocate( chlorophyl(nx,ny,nz)  ) 
      else
         write(*,*) "init_mesh_grid: biogeochem not allocated" 
      endif

c     --- 2D grids ---

      allocate( dslm(nx,ny)       )       
      allocate( wdepth(nx,ny)     ) 
      allocate( bottom_layer(nx,ny) )
      allocate( u_wind_stress(nx,ny) )
      allocate( v_wind_stress(nx,ny) )
      write(*,*) "init_mesh_grid: allocate grid arrays: OK"

      call init_horizontal_representation()

c
c     --- set default padvalue, if present ----
c     
      if(present(padding)) then
         padval_u        = padding
         padval_v        = padding    
         padval_w        = padding    
         padval_temp     = padding 
         padval_salinity = padding
         padval_vdiffus  = padding             
         padval_hdiffus  = padding 
         padval_dslm     = padding 
         padval_bgc      = padding
         padval_wdepth   = padding 
         write(*,*) "init_mesh_grid: default padding = ",padding
      else
         write(*,*) "init_mesh_grid: ",
     +       "no default padding provided, using", default_padding
      endif
c
c     --- overwrite with specific values, if present ----
c          
      if(present(pad_u))        padval_u        = pad_u 
      if(present(pad_v))        padval_v        = pad_v            
      if(present(pad_w))        padval_w        = pad_w          
      if(present(pad_temp))     padval_temp     = pad_temp     
      if(present(pad_salinity)) padval_salinity = pad_salinity 
      if(present(pad_vdiffus))  padval_vdiffus  = pad_vdiffus                  
      if(present(pad_hdiffus))  padval_hdiffus  = pad_hdiffus  
      if(present(pad_dslm))     padval_dslm     = pad_dslm   
      if(present(pad_bgc))      padval_bgc      = pad_bgc
      if(present(pad_wdepth))   padval_wdepth   = pad_wdepth
      
      write(*,*) "init_mesh_grid: resolved padding values:"
      write(*,*) "padval_u       :", padval_u
      write(*,*) "padval_v       :", padval_v
      write(*,*) "padval_w       :", padval_w
      write(*,*) "padval_temp    :", padval_temp
      write(*,*) "padval_salinity:", padval_salinity
      write(*,*) "padval_vdiffus :", padval_vdiffus
      write(*,*) "padval_hdiffus :", padval_hdiffus
      write(*,*) "padval_dslm    :", padval_dslm
      write(*,*) "padval_bgc     :", padval_bgc

      data_revision_tag = -1      !  -1 = unset data

      end subroutine init_mesh_grid
       



      subroutine close_mesh_grid()
c     ------------------------------------------------------
c     ------------------------------------------------------
      if (allocated(u))            deallocate( u )
      if (allocated(v))            deallocate( v )
      if (allocated(w))            deallocate( w )
      if (allocated(temp))         deallocate( temp )
      if (allocated(salinity))     deallocate( salinity )
      if (allocated(vdiffus))      deallocate( vdiffus )
      if (allocated(hdiffus))      deallocate( hdiffus )
      if (allocated(dslm))         deallocate( dslm )
      if (allocated(ccdepth))      deallocate( ccdepth )     
      if (allocated(acc_width))    deallocate( acc_width )     
      if (allocated(wdepth))       deallocate( wdepth )
      if (allocated(bottom_layer)) deallocate( bottom_layer )
      if (allocated(u_wind_stress)) deallocate(u_wind_stress )
      if (allocated(v_wind_stress)) deallocate(v_wind_stress )
c
      if (allocated( zoo )) deallocate( zoo )
      if (allocated( oxygen )) deallocate( oxygen )
      if (allocated( nh4 )) deallocate( nh4 )
      if (allocated( no3 )) deallocate( no3 )
      if (allocated( po4 )) deallocate( po4 )         
      if (allocated( diatoms )) deallocate( diatoms )     
      if (allocated( flagellates )) deallocate( flagellates )     
      if (allocated( cyanobacteria )) deallocate(cyanobacteria)       
      if (allocated( organic_detritus )) 
     +                        deallocate( organic_detritus )        
      if (allocated( part_org_matter)) 
     +                        deallocate( part_org_matter )
      if (allocated( dissolv_inorg_carbon )) 
     +                        deallocate( dissolv_inorg_carbon )
      if (allocated( alkalinity )) deallocate( alkalinity )
      if (allocated( dissolv_inorg_nitrogen )) 
     +                        deallocate( dissolv_inorg_nitrogen ) 
      if (allocated( chlorophyl )) deallocate( chlorophyl )

      data_revision_tag = -1      !  -1 = unset data

      call close_horizontal_representation()

      end subroutine close_mesh_grid

c
c     data revision monitoring tools
c

      subroutine get_data_rev_tag_mesh_grid(i16)
c     -------------------------------------------------------
c     Return current revision in i16 (with representation KIND = selected_int_kind(16))
c     i16 = -1 means unset data (or that mesh_grid clients
c     has not invoked update_data_rev_tag_mesh_grid) 
c     -------------------------------------------------------
      integer(kind=selected_int_kind(16)), intent(out) :: i16
      i16 = data_revision_tag ! 
      end subroutine get_data_rev_tag_mesh_grid


      subroutine update_data_rev_tag_mesh_grid()
c     -------------------------------------------------------
c     mesh_grid clients should invoke this method once
c     they are done with modifying mesh_grid data
c     This allows other mesh_grid clients only perform 
c     necessary updates.
c     Currently counter data_revision_tag is incremented
c     - we could also set a random number
c     -------------------------------------------------------
      if (data_revision_tag<0) data_revision_tag = 0
      if (data_revision_tag == huge(data_revision_tag)) 
     +   data_revision_tag = 0    ! HUGE(X) The largest positive number
      data_revision_tag = data_revision_tag + 1
      end subroutine update_data_rev_tag_mesh_grid
      



      subroutine interpolate_cc_3Dgrid_data(geo,array,deriv,padval,
     +                                      result,status)
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
c     asc@23Sep2011: added slicing 1:ibot for "zgrid => ccdepth ..." to avoid dependence
c                    on that ccdepth is initialized below bottom. search_sorted assumes an increasing array
c     --------------------------------------------------------------------
      real, intent(in)     :: geo(:),array(:,:,:)
      integer, intent(in)  :: deriv  ! 0=value; (1,2,3) = along (x,y,z)
      real, intent(in)     :: padval
      real, intent(out)    :: result
      integer, intent(out) :: status
c      
      integer             :: ix,iy,iz,i,idum,ibot,ix0,ix1,iy0,iy1
      real                :: z,sx,sy,sz,depth,fbar,fup,flow,dfzdzbar
      integer             :: cx(4), cy(4)
      real                :: fz(4), dz, dfzdz(4)
      logical             :: valid(4)
      real, pointer       :: zgrid(:)
c     --------------------------------------------------------------------
      if (.not.horizontal_range_check(geo)) then
          result = padval          ! value interpolation 
          if (deriv>0) result = 0  ! derivative interpolation
          status = 1
          return
      endif

      if (is_land(geo)) then
          result = padval          ! value interpolation 
          if (deriv>0) result = 0  ! derivative interpolation
          status = 3
          return
      endif

c     --- accept points at bottum and on water surface --- 
c     --- within a numerical tolerance htol            ---

      call interpolate_wdepth(geo,depth,idum)
      if ((geo(3)<-htol).or.(geo(3)>depth+htol)) then
         status = 2                ! signal vertical out-of-bound
         result = padval           ! value interpolation
         if (deriv>0) result = 0   ! derivative interpolation
         return
      endif
c
c     define corners in this order
c         col1  (ix0,iy0)
c         col2  (ix0,iy1)
c         col3  (ix1,iy0)  
c         col4  (ix1,iy1)
c
      call get_surrounding_box(geo,ix,iy,sx,sy)
c
      ix0 = min(max(ix,  1),nx)  ! cover grid edges
      ix1 = min(max(ix+1,1),nx)  ! cover grid edges
      iy0 = min(max(iy,  1),ny)  ! cover grid edges
      iy1 = min(max(iy+1,1),ny)  ! cover grid edges
      cx(1) = ix0; cy(1) = iy0 
      cx(2) = ix0; cy(2) = iy1 
      cx(3) = ix1; cy(3) = iy0 
      cx(4) = ix1; cy(4) = iy1 

      z          = geo(3)  ! short hand
      status     = 0       ! assume interpolation possible
      valid      = .false. ! default for dry/out-of-bounds pillars
      fz         = 0.      ! default for dry/out-of-bounds pillars (NOT padval!)
      dfzdz      = 0.      ! default for dry/out-of-bounds pillars 
c
c     ------ setup 4 interpolation pillars ------
c            project z onto these pillars and interpolate
c
      do i = 1,4
         if (wetmask(cx(i),cy(i))>0) then
            valid(i) = .true.
            ibot     = bottom_layer(cx(i),cy(i)) ! ibot >= 1
            zgrid    => ccdepth(cx(i),cy(i),1:ibot)   ! range = 1:ibot<=nz
            call search_sorted_list(z,zgrid,iz) ! 0<=iz<=nz
            ! ---- handle projection extrapolation 
            !      but do not flag flag vertical extrapolation
            !      since 0<z<depth
            if ((iz<1).or.(iz>=ibot)) then 
               iz = min(max(iz,1),ibot) ! now iz = 1 or ibot
               flow    = array(cx(i),cy(i),iz)
               fup     = flow    ! => dfzdz = 0
               sz      = 0.5     ! set dummy
c               status  = 2       ! flag vertical extrapolation
               dz      = 1.0     ! avoid numerical problems, assign dummy
            ! ---- interior linear vertical interpolation ----
            else
               flow = array(cx(i),cy(i),iz)
               fup  = array(cx(i),cy(i),iz+1)
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

         result = padval
         if (deriv>0) result = 0 
         status = 3  ! rank deficit exit
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
         write(*,*) "interpolate_cc_3Dgrid_data: deriv = ",
     +               deriv, "is not implemented"  
      endif

      end subroutine interpolate_cc_3Dgrid_data




      subroutine interpolate_cc_2Dgrid_data(geo,array,deriv,padval,
     +                                      result,status)
c     --------------------------------------------------------------------
c     Interpolate grid of corner centered data 2D array on point geo
c     with data points at integer valued geo
c
c     deriv = 0 gives value, deriv = (1,2) gives derivative along (x,y)
c     Currently, only deriv = 0, until other derivatives 
c     are needed.
c     Do not check vertival position for horizontal interpolations
c
c     Return status:
c       status = 0: interior interpolation performed
c       status = 1: horizontal range violation, set result = padval
c       status = 3: dry point / rank deficit situation not permitting interpolation
c                   Return result = padval for deriv = 0  
c     --------------------------------------------------------------------
      real, intent(in)     :: geo(:),array(:,:)
      integer, intent(in)  :: deriv  ! 0=value; deriv = (1,2) gives derivative along (x,y)
      real, intent(in)     :: padval
      real, intent(out)    :: result
      integer, intent(out) :: status
c      
      integer           :: i,ix,iy,ix0,ix1,iy0,iy1
      real              :: vc(4),sx,sy,vcbar
      integer           :: cx(4), cy(4)
      logical           :: valid(4)
c     --------------------------------------------------------------------
      if (.not.horizontal_range_check(geo)) then
          result = padval
          status = 1
          return
      endif
c     
      if (is_land(geo)) then
          result = padval          ! value interpolation 
          if (deriv>0) result = 0  ! derivative interpolation
          status = 3
          return
      endif

c
      call get_surrounding_box(geo,ix,iy,sx,sy)

      ix0 = min(max(ix,    1),nx)  ! cover boundary layers
      ix1 = min(max(ix + 1,1),nx)  ! cover boundary layers
      iy0 = min(max(iy,    1),ny)  ! cover boundary layers
      iy1 = min(max(iy + 1,1),ny)  ! cover boundary layers
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

    
      status     = 0       ! assume interpolation possible
      valid      = .false. ! default for dry/out-of-bounds pillars
      vc         = 0.      ! default for dry/out-of-bounds pillars (NOT padval!)
c
c     ------ setup 4 interpolation corners ------
c
      do i = 1,4
         if (wetmask(cx(i),cy(i))>0) then
            valid(i) = .true.
            vc(i)    = array(cx(i),cy(i))
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
      end subroutine interpolate_cc_2Dgrid_data



      subroutine interpolate_currents(geo, r3, status)
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
      integer              :: statu, statv, statw
c     ------------------------------------------
      call interpolate_cc_3Dgrid_data(geo,u,0,padval_u,r3(1),statu)
      call interpolate_cc_3Dgrid_data(geo,v,0,padval_v,r3(2),statv)
      call interpolate_cc_3Dgrid_data(geo,w,0,padval_w,r3(3),statw)
      status  = max(statu, statv, statw) ! in the unextected case taht they differ ...
c     ------------------------------------------ 
      end subroutine interpolate_currents


      subroutine interpolate_turbulence(geo, r3, status)
c     ------------------------------------------
c     Default interpolation assumes vdiffus is cell-centered
c     otherwise a customized interpolater must be applied
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r3(:)
      integer, intent(out) :: status
      integer              :: statv,stath
c     ------------------------------------------ 
      call interpolate_cc_3Dgrid_data(geo,hdiffus,0,padval_hdiffus,
     +                                r3(1),stath)
      call interpolate_cc_3Dgrid_data(geo,vdiffus,0,padval_vdiffus,
     +                                r3(3),statv)
      r3(2)   = r3(1) ! horizontal isotropy
      status  = max(statv,stath) ! in the unextected case that they differ ...
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
c     padval = 0 is hard coded for derivatives
      call interpolate_cc_3Dgrid_data(geo,vdiffus,3,0.,r3(3),statv)
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
      call interpolate_cc_3Dgrid_data(geo,temp,0,padval_temp,r,status)
c     ------------------------------------------ 
      end subroutine 


      subroutine interpolate_salty (geo, r, status) 
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
c     ------------------------------------------ 
      call interpolate_cc_3Dgrid_data(geo,salinity,0,padval_salinity, 
     +                                r,status)
c     ------------------------------------------ 
      end subroutine 


      subroutine interpolate_wdepth(geo, result, status) 
c     -------------------------------------------------- 
c     Do not apply interpolate_cc_2Dgrid_data which uses
c     wet-point interpolation, but unrestricted interpolation (wdepth=0 at dry points)
c     (otherwise problems in coastal regions)
c     -------------------------------------------------- 
      real, intent(in)     :: geo(:) 
      real, intent(out)    :: result
      integer, intent(out) :: status   
c
      integer           :: ix,iy,ix0,iy0,ix1,iy1
      real              :: sx,sy,vc(4)    
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
      iy1 = min(max(iy + 1,1),ny)  ! cover boundary layers
      vc(1) = wdepth(ix0, iy0)
      vc(2) = wdepth(ix0, iy1)
      vc(3) = wdepth(ix1, iy0)
      vc(4) = wdepth(ix1, iy1)

      call interp_2Dbox_data(sx, sy, vc, 0, result)
      status = 0
c     ------------------------------------------ 
      end subroutine 


      subroutine interpolate_wind_stress(geo, r2, status)
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r2(:)
      integer, intent(out) :: status
      integer              :: statu, statv
c     ------------------------------------------
      call interpolate_cc_2Dgrid_data(geo, u_wind_stress, 0, 
     +                                padval_uwind, r2(1), statu)
      call interpolate_cc_2Dgrid_data(geo, v_wind_stress, 0, 
     +                                padval_vwind, r2(2), statv)
      status  = max(statu, statv) ! in the unextected case that they differ ...
c     ------------------------------------------ 
      end subroutine interpolate_wind_stress



      subroutine interpolate_zooplankton (geo, r, status) 
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r(:)
      integer, intent(out) :: status
c     ------------------------------------------ 
      if (incl_biogeochem) then
         call interpolate_cc_3Dgrid_data(geo,zoo,0,
     +        padval_bgc,r(1),status)
      else
         write(*,*) "interpolate_zooplankton: "//
     +              "biogeochem is not initialized"
         stop
      endif
c     ------------------------------------------ 
      end subroutine 


      subroutine interpolate_oxygen (geo, r, status) 
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
c     ------------------------------------------ 
      if (incl_biogeochem) then
         call interpolate_cc_3Dgrid_data(geo,oxygen,0,
     +                 padval_bgc,r,status)
      else
         write(*,*) "interpolate_oxygen: "//
     +              "biogeochem is not initialized"
         stop
      endif   
c     ------------------------------------------ 
      end subroutine 

      subroutine interpolate_nh4 (geo, r, status)  
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
c     ------------------------------------------ 
      if (incl_biogeochem) then
         call interpolate_cc_3Dgrid_data(geo,nh4,0,
     +                 padval_bgc,r,status)
      else
         write(*,*) "interpolate_nh4: "//
     +              "biogeochem is not initialized"
         stop
      endif   
c     ------------------------------------------ 
      end subroutine 


      subroutine interpolate_no3(geo, r, status) 
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
c     ------------------------------------------ 
      if (incl_biogeochem) then
         call interpolate_cc_3Dgrid_data(geo,no3,0,
     +                 padval_bgc,r,status)
      else
         write(*,*) "interpolate_no3: "//
     +              "biogeochem is not initialized"
         stop
      endif   
c     ------------------------------------------ 
      end subroutine   

           
      subroutine interpolate_po4(geo, r, status)  
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
c     ------------------------------------------ 
      if (incl_biogeochem) then
         call interpolate_cc_3Dgrid_data(geo,po4,0,
     +                 padval_bgc,r,status)
      else
         write(*,*) "interpolate_po4: "//
     +              "biogeochem is not initialized"
         stop
      endif   
c     ------------------------------------------ 
      end subroutine    

              
      subroutine interpolate_diatoms(geo, r, status) 
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
c     ------------------------------------------ 
      if (incl_biogeochem) then
         call interpolate_cc_3Dgrid_data(geo,diatoms,0,
     +                 padval_bgc,r,status)
      else
         write(*,*) "interpolate_diatoms: "//
     +              "biogeochem is not initialized"
         stop
      endif   
c     ------------------------------------------ 
      end subroutine        

      
      subroutine interpolate_flagellates(geo, r, status) 
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
c     ------------------------------------------ 
      if (incl_biogeochem) then
         call interpolate_cc_3Dgrid_data(geo,flagellates,0,
     +                 padval_bgc,r,status)
      else
         write(*,*) "interpolate_flagellates: "//
     +              "biogeochem is not initialized"
         stop
      endif   
c     ------------------------------------------ 
      end subroutine     

      
      subroutine interpolate_cyanobacteria (geo, r, status) 
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
c     ------------------------------------------ 
      if (incl_biogeochem) then
         call interpolate_cc_3Dgrid_data(geo,cyanobacteria,0,
     +                 padval_bgc,r,status)
      else
         write(*,*) "interpolate_cyanobacteria: "//
     +              "biogeochem is not initialized"
         stop
      endif   
c     ------------------------------------------ 
      end subroutine  

     
      subroutine interpolate_organic_detritus(geo, r, status) 
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
c     ------------------------------------------ 
      if (incl_biogeochem) then
         call interpolate_cc_3Dgrid_data(geo,organic_detritus,0,
     +                 padval_bgc,r,status)
      else
         write(*,*) "interpolate_organic_detritus: "//
     +              "biogeochem is not initialized"
         stop
      endif   
c     ------------------------------------------ 
      end subroutine   

           
      subroutine interpolate_part_org_matter(geo, r, status) 
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
c     ------------------------------------------ 
      if (incl_biogeochem) then
         call interpolate_cc_3Dgrid_data(geo,part_org_matter,0,
     +                 padval_bgc,r,status)
      else
         write(*,*) "interpolate_part_org_matter: "//
     +              "biogeochem is not initialized"
         stop
      endif   
c     ------------------------------------------ 
      end subroutine 


      subroutine interpolate_DIC(geo, r, status) 
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
c     ------------------------------------------ 
      if (incl_biogeochem) then
         call interpolate_cc_3Dgrid_data(geo,dissolv_inorg_carbon,0,
     +                 padval_bgc,r,status)
      else
         write(*,*) "interpolate_DIC: "//
     +              "biogeochem is not initialized"
         stop
      endif   
c     ------------------------------------------ 
      end subroutine 


      subroutine interpolate_alkalinity(geo, r, status) 
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
c     ------------------------------------------ 
      if (incl_biogeochem) then
         call interpolate_cc_3Dgrid_data(geo,alkalinity,0,
     +                 padval_bgc,r,status)
      else
         write(*,*) "interpolate_alkalinity: "//
     +              "biogeochem is not initialized"
         stop
      endif   
c     ------------------------------------------ 
      end subroutine 


      subroutine interpolate_DIN(geo, r, status) 
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
c     ------------------------------------------ 
      if (incl_biogeochem) then
         call interpolate_cc_3Dgrid_data(geo,dissolv_inorg_nitrogen,0,
     +                 padval_bgc,r,status)
      else
         write(*,*) "interpolate_DIN: "//
     +              "biogeochem is not initialized"
         stop
      endif   
c     ------------------------------------------ 
      end subroutine 


      subroutine interpolate_chlorophyl(geo, r, status) 
c     ------------------------------------------ 
c     ------------------------------------------ 
      real, intent(in)     :: geo(:)
      real, intent(out)    :: r
      integer, intent(out) :: status
c     ------------------------------------------ 
      if (incl_biogeochem) then
         call interpolate_cc_3Dgrid_data(geo,chlorophyl,0,
     +                 padval_bgc,r,status)
      else
         write(*,*) "interpolate_chlorophyl: "//
     +              "biogeochem is not initialized"
         stop
      endif   
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
      if (status==1) then     ! status = 1: horizontal range violation, set result = padval
         is_wet = .true.
         return
      elseif (status==0) then ! status = 0: interior interpolation performed
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


      subroutine get_grid_coordinates(geo,x,y,z,dxyz_dr)  ! formerly named get_ncc_coordinates
c     -------------------------------------------------------------------------------- 
c     Get continuous node-centered grid coordinates with grid points 
c     on integer values of (x,y,z) from geo = (lon,lat,depth)
c     water surface is at z = 0.5, sea bed at z = bottum_layer+0.5
c     layer boundaries at half-integers inbetween.
c     It is not checked that z is above the sea bed. The inter grid
c     range is 0.5 <= z <= nz+0.5.
c     If vertical range is exceeded the first/last layer, respectively,
c     is used to extrapolate a vertical grid coordinate smoothly
c     (no extrapolation is flagged) 
c     Include intracell interpolation of layer spacings
c     Optionally compute grid coordinate derivative wrt. Cartesian coordinates 
c     if argument dxyz_dr is provided. dxyz_dr(i,j) is ds_j/dr_i with s = (x,y,z) so that
c     dfdr = matmul(dxyz_dr, df_dxyz)
c     -------------------------------------------------------------------------------- 
      real, intent(in)  :: geo(:)
      real, intent(out) :: x,y,z
      real, intent(out), optional :: dxyz_dr(:,:) ! required shape at least (3,3)
      
      integer           :: ix,iy,ix0,iy0,ix1,iy1,iz
      
      real, pointer     :: z00(:), z01(:)
      real, pointer     :: z10(:), z11(:)
      real, allocatable :: z0(:), z1(:), acclay(:)
      real, allocatable :: dacc_dX(:), dacc_dY(:)
      real              :: layerw,sx,sy,cwd
      real              :: dlw_dX,dlw_dY,facx,facy
c     ------------------------------------------
      allocate( z0(nz+1)     )
      allocate( z1(nz+1)     )
      allocate( acclay(nz+1) )
      allocate( dacc_dX(nz+1)) ! Cartesian derivative of acclay, for dxyz_dr
      allocate( dacc_dY(nz+1)) ! Cartesian derivative of acclay, for dxyz_dr
      call get_horiz_grid_coordinates(geo,x,y)   ! define x,y
c
c     ----- interpolate local layer spacings from acc_width: acclay -----
c
      call get_surrounding_box(geo,ix,iy,sx,sy)
c     
      ix0 = min(max(ix,  1),nx)  ! cover grid edges
      ix1 = min(max(ix+1,1),nx)  ! cover grid edges
      iy0 = min(max(iy,  1),ny)  ! cover grid edges
      iy1 = min(max(iy+1,1),ny)  ! cover grid edges
c     --- trilin interpolation  ---
      z00 => acc_width(ix0, iy0, :)  ! range = 1:nz+1
      z01 => acc_width(ix0, iy1, :)  ! range = 1:nz+1
      z10 => acc_width(ix1, iy0, :)  ! range = 1:nz+1
      z11 => acc_width(ix1, iy1, :)  ! range = 1:nz+1
      z0  = z00 + sx*(z10 - z00)   ! vector operation (south face)
      z1  = z01 + sx*(z11 - z01)   ! vector operation (north face)
      acclay = z0 + sy*(z1-z0)  ! vector operation (north-south)
c
c     ----- interpolate z from local layer spacings -----
c
c     search_sorted_list: result = iz (0 <= iz <= nz+1)
c                         acclay(iz) < geo(3) < acclay(iz+1)   (iz<=nz)
c                         acclay(1) = 0        
c
      call search_sorted_list(geo(3), acclay, iz)  ! iz is layer numer, where geo(3) belongs
      iz = min(max(iz,1),nz) ! capture vertical range excess 
      layerw = acclay(iz+1) -  acclay(iz)
      layerw = max(layerw, 1e-6)  ! avoid numerical exceptions for źero width layers
      z      = (iz - 0.5) + (geo(3) - acclay(iz))/layerw ! provides smooth extrapolation at range excess
c
c     --- optionally evaluate Cartesian derivative (X,Y,Z) of (x,y,z) ---
c
      if ( present(dxyz_dr) ) then
         dxyz_dr(1:3,1:3) = 0.0
         call get_horiz_grid_coor_deriv(geo, dxyz_dr(1:2,1:2))
         dacc_dX = (z10-z00)*dxyz_dr(1,1) + (z01-z00)*dxyz_dr(1,2) +
     +           (z11+z00-z01-z10)*(sx*dxyz_dr(1,2) + dxyz_dr(1,1)*sy)
         dacc_dY = (z10-z00)*dxyz_dr(2,1) + (z01-z00)*dxyz_dr(2,2) +
     +           (z11+z00-z01-z10)*(sx*dxyz_dr(2,2) + dxyz_dr(2,1)*sy)
         dlw_dX       = dacc_dX(iz+1) - dacc_dX(iz)
         dlw_dY       = dacc_dY(iz+1) - dacc_dY(iz) 
         dxyz_dr(3,3) = 1/layerw ! dz/dZ
         dxyz_dr(1,3) = -dacc_dX(iz)*layerw - (geo(3)-acclay(iz))*dlw_dX ! dz/dX
         dxyz_dr(1,3) = dxyz_dr(1,3)/layerw**2
         dxyz_dr(2,3) = -dacc_dY(iz)*layerw - (geo(3)-acclay(iz))*dlw_dY ! dz/dY
         dxyz_dr(2,3) = dxyz_dr(2,3)/layerw**2
      endif      
      deallocate( z0, z1, acclay, dacc_dX, dacc_dY)
      nullify( z00, z01, z10, z11 )
c     ------------------------------------------ 
      end subroutine get_grid_coordinates


cccccccccccccccccccc     diagnostic   subroutines   ccccccccccccccccccc

      subroutine write_topomask(iunit, dry)
c     -----------------------------------------------------------------
c     print wet/dry grid points to logical unit iunit
c     dry = .TRUE.  -> print dry (i.e. land) grid points (default)
c     dry = .FALSE. -> print wet grid points (+ water depth at grid point)
c     wet/dry condition accessed from wetmask (must have been initialized)
c     water depth accessed with interpolate_wdepth
c     -----------------------------------------------------------------      
      integer, intent(in)           :: iunit
      logical, intent(in),optional  :: dry
      logical  :: printdry
      real     :: geo(2),wd
      integer  :: ix,iy,idum
c     -----------------------------------------------------------------      
      if (present(dry).and.(dry==.false.)) then
         printdry = .false.
      else                     
         printdry = .true.
      endif
      do ix=1, nx
         do iy=1,ny
            call get_horiz_geo_coordinates(1.0*ix, 1.0*iy, geo)
            if (printdry) then
               if (wetmask(ix,iy)==0) write(iunit,*) geo
            else
               if (wetmask(ix,iy)==1) then
                  call interpolate_wdepth(geo,wd,idum)
                  write(iunit,*) geo, wd
               endif   
            endif
         enddo
      enddo
      end subroutine write_topomask

      end module
