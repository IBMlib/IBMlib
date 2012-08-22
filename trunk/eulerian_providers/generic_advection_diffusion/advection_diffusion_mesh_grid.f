      module advection_diffusion
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     This module represents implements a preliminary interface for 
c     Eulerian simulation of production combined advection-diffusion
c
c     The module provides a class eulerian_field representing this
c     Several independent instances can be propagated by update_eulerian_field
c     Each  Eulerian field is associated with a horizontal sub grid which is a 
c     dense sub set of the parent grid. Sub grid is selected at instantiation
c     No vertical sub gridding is currently supported
c 
c     This implementation is based on the following staggering:
c
c        u,v positive along lambda,phi in units m/s
c        w   positive downward         in units m/s  
c        u(ix,iy,iz) in grid position (ix+0.5, iy    , iz)     (i.e. eastern  cell face)
c        v(ix,iy,iz) in grid position (ix    , iy-0.5, iz)     (i.e. southern cell face)
c        w(ix,iy,iz) in grid position (ix    , iy,     iz-0.5) (i.e. upper    cell face)
c
c     Cell geometry: it is currently assumed that coordinate vectors are
c                    perdendicular when generating cell area/cell volume arrays 
c                    - this assumption can later on be relaxed
c
c     Eulerian dynamics is performed on a sub grid of the parent grid. Multiple sub grids 
c     (for different Eulerian fields) are allowed; sub grid may be parent domain
c     Vertically, each sub grid must start at the surface, but needs not go to the bottom.
c
c     Divergence correction (input tag eulerian_div_corr, default value .false.). 
c     When .true. a term is added at forward integration step so that a homogeneous field
c     stay homogeneous when water fluxes across cell faces are not summing to exactly zero.
c     When .false. forward integration step exactly conserves mass (NB water elevation fluctuations may
c     interfere with mass conservation between updates of hydrographic frames, to be checked further)
c 
c     F90 notice: if an array is allocated with offset /= 1 like allocate(data(0:n))
c                 a dummy argument preserving the offset /= 1 can be declared as real :: dummy(0:)
c    
c     $Rev: $
c     $LastChangedDate:  $
c     $LastChangedBy: $ 
c     
c     FIXES/AMENDMENTS:
c         * May 23 2012: made cell face flux evaluation symmetric so that it also works for negative currents
c         * May 24 2012: optional divergence correction in update_eulerian_field (corrects for water fluxes
c                        into a control volume not being exactly zero)
c     TODO:
c         cleanup before commit 
c         non conservation correction
c         Eulerian subgrid class
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use physical_fields          ! injects nx,ny
      use mesh_grid
      use horizontal_grid_transformations, only:
     +           get_horiz_geo_coordinates, get_horiz_grid_coordinates
      use geometry,      only: get_horizontal_distance
      use eulerian_production_rate ! -> point_production_rate
      use run_context, only: simulation_file
      use input_parser


      implicit none
      private           ! default scope

c     ---------- set scope ----------

c     --- generic interface ---
      public :: init_advection_diffusion   ! init  this module
      public :: close_advection_diffusion  ! close this module
      
      public :: init_eulerian_field
      public :: close_eulerian_field
      public :: interpolate_eulerian_field
      public :: update_eulerian_field
      public :: update_eulerian_auxillaries

c     --- specific interface (for this implementation) ---
c         other implementations of module advection_diffusion may not provide these 
c
      public :: dump_2Dframe   ! dump a data layer to file
      public :: get_data       ! export pointer to data

c     ------------------------------------------------
c     The interior domain of the field is data(1:nxs, 1:nys, 1:nzs)  
c     elements with 0 and/or n*s+1 are rim points (needed for making field gradients)
c     point (i,j,k) on sub grid corresponds to (i+ixoff, j+iyoff, k) on full grid so 
c     that (ixoff, iyoff) == (0,0) means coincident indexing with parent grid
c     No z offset is provided (always start include the surface)
c     Notice that no assumptions are made about the unit X of the quantity 
c     described by the field, but data(:,:,:) have implied unit X/m3 
c     ------------------------------------------------
      type eulerian_field
c      private
        real,pointer    :: data(:,:,:)  ! concentration as X/m3 in sub grid cell (i,j,k)
        integer         :: ixoff,iyoff  ! subgrid offsets in full grid (nx,ny,nz) (zero means coincident)
        real,pointer    :: dataBC(:)    ! constant applied as BC (vertically layered)
c       ------ cell falgs to make it easier to perform Eulerian dynamics ------
        logical,pointer :: fluxE(:,:,:)   ! whether flux over East side of cell
        logical,pointer :: fluxS(:,:,:)   ! whether flux over South side of cell
        logical,pointer :: fluxR(:,:,:)   ! whether flux over cell roof
        logical,pointer :: dyncell(:,:,:) ! whether cell concentration may change
c       other handles
      end type
      public :: eulerian_field


c     -------------------- module data --------------------  
      
      real, allocatable :: facelen_E(:,:) ! length of E face of cell (i,j)     (static)
      real, allocatable :: facelen_S(:,:) ! length of S face of cell (i,j)     (static)
      real, allocatable :: area_E(:,:,:)  ! area of E face of cell (i,j,k)     (dynamic)
      real, allocatable :: area_S(:,:,:)  ! area of S face of cell (i,j,k)     (dynamic)
      real, allocatable :: area_R(:,:)    ! area of upper face of cell (i,j)   (static)
      real, allocatable :: dx(:,:)        ! distance [m] between cell centers (i,j,k) and (i+1,j,k) (static)
      real, allocatable :: dy(:,:)        ! distance [m] between cell centers (i,j,k) and (i,j-1,k) (static)
      real, allocatable :: dz(:,:,:)      ! distance [m] between cell centers (i,j,k) and (i,j,k-1) (dynamic)
      real, allocatable :: cellvol(:,:,:) ! volume of cell (i,j,k) in [m3]                          (dynamic)
      real, allocatable :: ku(:,:,:)      ! diffusivity [m2/s] at E face of cell (i,j,k)            (dynamic)
      real, allocatable :: kv(:,:,:)      ! diffusivity [m2/s] at S face of cell (i,j,k)            (dynamic)
      real, allocatable :: kw(:,:,:)      ! diffusivity [m2/s] at upper face of cell (i,j,k)        (dynamic)
     
c     -------------------- module control handlers --------------------
      logical  :: div_corr ! correct for sum of water fluxes across cell faces are not exactly zero

c     ===================================================
                            contains
c     ===================================================

                      
      subroutine init_advection_diffusion()
c     ------------------------------------------------------
c     Module initialization
c     Allocate auxillary arrays shared by different Eulerian fields
c     Auxillary arrays cover full grid area, since the set
c     of Eulerian fields handled by this module is case specific and not necessary static 
c     This subroutine must be invoked after init_physical_fields
c     so that geometric aspects are initialized
c     Do not assume that update_physical_fields has been invoked
c     ------------------------------------------------------
      integer :: i,j,k
      real    :: pos0(2), posE(2), posS(2), pos1(2), lenx, leny
c     ------------------------------------------------------
      allocate( facelen_E(nx,ny) )  ! length of E face of cell (i,j)     (static)  
      allocate( facelen_S(nx,ny) )  ! length of S face of cell (i,j)     (static)   
      allocate( area_E(nx,ny,nz) )  ! area of E face of cell (i,j,k)     (dynamic)
      allocate( area_S(nx,ny,nz) )  ! area of S face of cell (i,j,k)     (dynamic)
      allocate( area_R(nx,ny) )     ! area of upper face of cell (i,j)   (static)
      allocate( dx(nx,ny) )         ! distance [m] between cell centers (i,j) and (i+1,j)     (static)
      allocate( dy(nx,ny) )         ! distance [m] between cell centers (i,j) and (i,j-1)     (static)
      allocate( dz(nx,ny,nz) )      ! distance [m] between cell centers (i,j,k) and (i,j,k-1) (dynamic)
      allocate( cellvol(nx,ny,nz) ) ! volume of cell (i,j,k) in [m3]                          (dynamic)
      allocate( ku(nx,ny,nz) )      ! diffusivity [m2/s] at E face of cell (i,j,k)            (dynamic)
      allocate( kv(nx,ny,nz) )      ! diffusivity [m2/s] at S face of cell (i,j,k)            (dynamic)
      allocate( kw(nx,ny,nz) )      ! diffusivity [m2/s] at upper face of cell (i,j,k)        (dynamic)
c
c     ------ setup dx,dy: currently base all horizontal geometry on this
c                         only set for center connections inside grid
c
      dx = 0.0
      dy = 0.0
      do j=2,ny
         do i=1,nx-1
            call get_horiz_geo_coordinates(1.0*i,     1.0*j,     pos0)
            call get_horiz_geo_coordinates(1.0*(i+1), 1.0*j,     posE)
            call get_horiz_geo_coordinates(1.0*i,     1.0*(j-1), posS)
            call get_horizontal_distance(pos0, posE, dx(i,j))   ! module geometry
            call get_horizontal_distance(pos0, posS, dy(i,j))   ! module geometry
         enddo
      enddo
c     ------ setup area_R ------
      area_R = 0.0
      do j=1,ny
         do i=1,nx
c           ------ probe East face length ------
            call get_horiz_geo_coordinates(i+0.5, j+0.5, pos0)
            call get_horiz_geo_coordinates(i+0.5, j-0.5, pos1)
            call get_horizontal_distance(pos0, pos1, lenx)   ! module geometry
            facelen_E(i,j) = lenx
c           ------ probe South face length ------            
            call get_horiz_geo_coordinates(i+0.5, j-0.5, pos0)
            call get_horiz_geo_coordinates(i-0.5, j-0.5, pos1)
            call get_horizontal_distance(pos0, pos1, leny)   ! module geometry
            facelen_S(i,j) = leny
c           ------ assume orthogonality of dx and dy - can be generalized later ------
            area_R(i,j) = lenx*leny
         enddo
      enddo
c
c     ---- load module options ----
c
      if (count_tags(simulation_file, "eulerian_div_corr")>0) then
         call read_control_data(simulation_file, 
     +                          "eulerian_div_corr", div_corr)
      else 
         div_corr = .false. ! set default: no correction applied
      endif
      if (div_corr) then
         write(*,*) "init_advection_diffusion: "//
     +              "apply Eulerian divergence correction"
      else
         write(*,*) "init_advection_diffusion: "//
     +              "no Eulerian divergence correction applied"
      endif

c
c     ---- pass init to sub modules ----
c

      call init_eulerian_production_rate()

c
      end subroutine init_advection_diffusion 
       



      subroutine close_advection_diffusion()
c     ------------------------------------------------------
c     Module close-down
c     Deallocate modula date
c     Instantiated eulerian_field objects must be deallocated in external context
c     ------------------------------------------------------
      deallocate( facelen_E )
      deallocate( facelen_S )
      deallocate( area_E )  
      deallocate( area_S )  
      deallocate( area_R )     
      deallocate( dx )        
      deallocate( dy )      
      deallocate( dz )      
      deallocate( cellvol ) 
      deallocate( ku )      
      deallocate( kv )     
      deallocate( kw )       
c
c     ---- pass close to sub modules ----
c
      call close_eulerian_production_rate()
c
      end subroutine close_advection_diffusion



      subroutine init_eulerian_field(eulerf, geoSW, geoNE, 
     +                               data0, dataBC)
c     ------------------------------------------------------
c     Eulerian field constructor
c     
c     geoSW, geoNE are corners of Eulerian sub grid
c     geoSW, geoNE are projected onto parent grid, so that
c     sub grid is a dense sub set of parent grid
c
c     Currently, only the full water coloumn (within horizontal domain) can be selected 
c     It is somewhat tricky to pick a part of the water coloumn, because the sea surface
c     height is (usually) variable and the bottom topography is (usually) variable over 
c     the horizontal sub grid domain. However, accept vectors geoSW, geoNE as length 2+
c     for future implementation
c
c     Index in data0/dataBC refers to vertical layer. If length of is 1, data
c     is applied to all layers. If length of data is > 1, data is applied layer wise vertically
c     If length < nzs, last layer is padded downward
c     ------------------------------------------------------
      type(eulerian_field),intent(inout) :: eulerf    ! instance to be initialized
      real, intent(in)                   :: geoSW(:)  ! SW corner of Eulerian grid as (lon,lat,depth) in degrees
      real, intent(in)                   :: geoNE(:)  ! NE corner of Eulerian grid as (lon,lat,depth) in degrees 
      real, intent(in)                   :: data0(:)  ! initialization value of Eulerian field (layer wise)
      real, intent(in)                   :: dataBC(:) ! BC appplied to boundary layer after each time step (layer wise)
c
      real    :: xSW, ySW, xNE, yNE, d
      integer :: nxs, nys, nzs, i, j, k, kmax
      integer :: ixoff, iyoff, ixend, iyend
c     ------------------------------------------------------
c 
c     resolve requested sub grid domain
c
      call get_horiz_grid_coordinates(geoSW, xSW, ySW)
      call get_horiz_grid_coordinates(geoNE, xNE, yNE)
      ixoff = int(xSW + 0.5) - 1  ! so that SW corner corresponds to ixoff + 1 on full grid
      iyoff = int(ySW + 0.5) - 1  ! so that SW corner corresponds to iyoff + 1 on full grid
      ixend = int(xNE + 0.5)      ! full grid index of NE corner
      iyend = int(yNE + 0.5)      ! full grid index of NE corner
      if ((ixoff < 1).or.(iyoff < 1)) then  ! NB: need rim at parent grid as cushin
         write(*,*) "init_eulerian_field: SW corner outside parent grid", geoSW
         stop
      endif
      if ((ixend >= nx).or.(iyend >= ny)) then  ! NB: need rim at parent grid as cushin
         write(*,*) "init_eulerian_field: NE corner outside parent grid", geoNE
         stop
      endif
      nxs = ixend - ixoff 
      nys = iyend - iyoff 
      if ((nxs<1).or.(nys<1)) then
         write(*,*) "init_eulerian_field: invalid sub grid", nxs,nys
         write(*,*) "SW corner = ", geoSW
         write(*,*) "NE corner = ", geoNE
         stop
      endif    

c     Currently, only the full water coloumn (within the interior horizontal domain 
c     corresponding to the sub grid) can be selected 
      nzs = maxval(bottom_layer(ixoff+1:ixend, iyoff+1:iyend))  ! nzs <= nz
      if ((nzs<1).or.(nzs>nz)) then
         write(*,*) "init_eulerian_field: invalid nzs for sub grid", nzs
         stop
      endif   
c     
c     report resolved subgrid
c 
      write(*,*) "init_eulerian_field: sub grid dims = ", nxs,nys,nzs
      write(*,*) "init_eulerian_field: interior x index range = ", 
     +           ixoff+1, ixend
      write(*,*) "init_eulerian_field: interior y index range = ", 
     +           iyoff+1, iyend
      write(*,*) "init_eulerian_field: interior z index range = ", 
     +           1,nzs
      write(*,*) "init_eulerian_field: field initial value = ", data0
      write(*,*) "init_eulerian_field: Dirichlet BC        = ", dataBC
c 
c     initialize instance
c
      allocate ( eulerf%data(0:nxs+1, 0:nys+1, 1:nzs) ) ! no vertical cushin layers
      allocate ( eulerf%dataBC(nzs) ) 
      eulerf%ixoff  = ixoff
      eulerf%iyoff  = iyoff
c     --- set initial condition --- 
      kmax        = size(data0)    ! kmax >= 1
      eulerf%data = data0(kmax)    ! padding value
      do k = 1, min(nzs, kmax-1)   ! possibly void loop if kmax < 2
         eulerf%data(:,:,k) = data0(k)
      enddo
c     --- parse and set boundary condition ---
      kmax          = size(dataBC) ! kmax >= 1
      eulerf%dataBC = dataBC(kmax) ! padding value
      do k = 1, min(nzs, kmax-1)   ! possibly void loop if kmax < 2
         eulerf%dataBC(k) = dataBC(k)
      enddo
      call impose_boundary_conditions(eulerf)
c 
c     initialize topographic handles 
c
      allocate ( eulerf%fluxE(  0:nxs+1, 0:nys+1, 1:nzs) )
      allocate ( eulerf%fluxS(  0:nxs+1, 0:nys+1, 1:nzs) )
      allocate ( eulerf%fluxR(  0:nxs+1, 0:nys+1, 1:nzs+1) )
      allocate ( eulerf%dyncell(0:nxs+1, 0:nys+1, 1:nzs) )
      eulerf%fluxE   = .FALSE.  ! default
      eulerf%fluxS   = .FALSE.  ! default
      eulerf%fluxR   = .FALSE.  ! default
      eulerf%dyncell = .FALSE.  ! default 
c
c     wet_cell() is a local function
c
      do k=1,nzs
         do j=1,nys
            do i=0,nxs
               if (wet_cell(i,j,k).and.wet_cell(i+1,j,k)) then
                  eulerf%fluxE(i,j,k) = .TRUE.
               endif
            enddo
         enddo
      enddo
c
      do k=1,nzs
         do j=1,nys+1
            do i=1,nxs
               if (wet_cell(i,j,k).and.wet_cell(i,j-1,k)) then
                  eulerf%fluxS(i,j,k) = .TRUE.
               endif
            enddo
         enddo
      enddo
c
      do k=2,nzs    ! no flux through sea surface (k=1) of tracer
         do j=1,nys
            do i=1,nxs
               if (wet_cell(i,j,k).and.wet_cell(i,j,k-1)) then
                  eulerf%fluxR(i,j,k) = .TRUE.
               endif
            enddo
         enddo
      enddo      
c
      do k=1,nzs
         do j=1,nys
            do i=1,nxs
               if (wet_cell(i,j,k)) eulerf%dyncell(i,j,k) = .TRUE.
            enddo
         enddo
      enddo      
c
      contains  ! -------- local functions/subroutines -------- 
c
      logical function wet_cell(ixs,iys,izs)
c     --------------------------------------------
c     Campact test of whether sub grid cell (ixs,iys,izs) is wet/dry
c     Local scope access to offsets ixoff,iyoff
c     --------------------------------------------
      integer, intent(in) :: ixs,iys,izs
c     --------------------------------------------
      if (wetmask(ixs+ixoff,iys+iyoff) == 1) then
         if ((izs > bottom_layer(ixs+ixoff,iys+iyoff)).or.(izs<1)) then
            wet_cell = .FALSE.
         else
            wet_cell = .TRUE.
         endif   
      else
         wet_cell = .FALSE.
      endif
      end function wet_cell
c
      end subroutine init_eulerian_field

      


      subroutine close_eulerian_field(eulerf)
c     ------------------------------------------------------
c     Eulerian field destructor
c     ------------------------------------------------------
      type(eulerian_field),intent(inout) :: eulerf
c     ------------------------------------------------------
      if (associated(eulerf%data))    deallocate (eulerf%data)
      if (associated(eulerf%dataBC))  deallocate (eulerf%dataBC)
      if (associated(eulerf%fluxE))   deallocate (eulerf%fluxE)
      if (associated(eulerf%fluxS))   deallocate (eulerf%fluxS)
      if (associated(eulerf%fluxR))   deallocate (eulerf%fluxR)
      if (associated(eulerf%dyncell)) deallocate (eulerf%dyncell)
      end subroutine close_eulerian_field



      subroutine resolve_inner_domain_dims(eulerf,nxs,nys,nzs)
c     ------------------------------------------------------
c     The interior domain of the field is data(1:nxs, 1:nys, 1:nzs)  
c     field is allocated with index ranges (0:nxs+1, 0:nys+1, 1:nzs) 
c     ------------------------------------------------------
      type(eulerian_field),intent(in) :: eulerf
      integer,intent(out)             :: nxs,nys,nzs
c     ------------------------------------------------------
      nxs = size(eulerf%data, 1) - 2  ! interior domain 1:nxs
      nys = size(eulerf%data, 2) - 2  ! interior domain 1:nys
      nzs = size(eulerf%data, 3)      ! no vertical cushin layers
      end subroutine resolve_inner_domain_dims



      subroutine interpolate_eulerian_field(geo, eulerf, result, status)
c     ------------------------------------------------------------------
c     Interpolate Eulerian field eulerf -> result at position geo
c     Currently no interpolation between neighboring points, just 
c     return field value of the cell geo belongs to
c
c     Return status:
c       status = 0: interior interpolation performed
c       status = 1: horizontal range violation
c       status = 2: vertical range violation
c       status = 3: dry point / rank deficit situation not permitting interpolat
c
c     For return status /= 0, set status = 0.0 (no optinal pad value currently)
c     ------------------------------------------------------------------
      type(eulerian_field),intent(in) :: eulerf
      real,intent(in)                 :: geo(:)
      real,intent(out)                :: result
      integer,intent(out)             :: status
      real    :: x,y,z,xs,ys
      integer :: ixs,iys,izs,nxs, nys, nzs, ixoff,iyoff
c     ------------------------------------------------------------------
      status = 0    ! set default
      result = 0.0  ! set default

      if (is_land(geo)) then
         status = 3     ! flag dry point 
         return
      endif

      call get_grid_coordinates(geo,x,y,z)               !  on parent grid
      call resolve_inner_domain_dims(eulerf,nxs,nys,nzs) ! of sub grid
c
      ixoff = eulerf%ixoff 
      iyoff = eulerf%iyoff
      xs    = x - ixoff            ! convert xy to actual sub grid coordinate
      ys    = y - iyoff            ! convert xy to actual sub grid coordinate
      
      ixs   = int(xs + 0.5)              ! sub grid cell association
      iys   = int(ys + 0.5)              ! sub grid cell association
      izs    = int(z + 0.5)  
      if ((ixs>nxs).or.(ixs<1).or.(iys>nys).or.(iys<1)) then
         status = 1     ! flag horizontal range violation
         return
      endif
      if ((izs > bottom_layer(ixs+ixoff,iys+iyoff)).or.
     +    (izs > nzs).or.(izs<1)) then
         status = 2     ! flag vertical range violation 
         return
      endif
c
      result = eulerf%data(ixs,iys,izs) ! regular look-up
c
      end subroutine interpolate_eulerian_field



      subroutine update_eulerian_auxillaries()
c     ------------------------------------------------------------------
c     Update auxillary (full parent grid) fields needed for Eulerian dynamics
c     Must be invoked once for every time step (after update_physical_fileds), 
c     before update_eulerian_field are applied to active eulerian_field instances
c
c     Update all variable arrays derived from geometry/oceanography
c     which are shared between tracers: cell areas/volumes, diffusivity
c     Assume hdiffus is cell-centered
c     Assume vdiffus is cell-centered
c     Assume coordinate vectors (WE, NS) are locally orthogonal
c     ------------------------------------------------------------------
      integer   :: i,j,k, iea,jso
      real      :: h0, hE, hS
c     ------------------------------------------------------------------

c      u=1.0
c      v=1.0
c      w=0
c      write(*,*) "update_eulerian_auxillaries: WARNING - booted uvw"

      dz = 0.0
      do k=2,nz ! k=1 not def
         do j=1,ny
            do i=1,nx
               dz(i,j,k) = ccdepth(i,j,k)-ccdepth(i,j,k-1)
            enddo
         enddo
      enddo

      area_E  = 0.0
      area_S  = 0.0
      cellvol = 1.0 ! avoid divergence at division
      do k=1,nz
         do j=1,ny
            jso = max(1,j-1)      ! southern cell index (if any)
            do i=1,nx
               iea = min(nx,i+1)  ! eastern cell index (if any)
               h0             = acc_width(i,  j,k+1)-acc_width(i,  j,k)   ! width at cell center
               hE             = acc_width(iea,j,k+1)-acc_width(iea,j,k)   ! width at East  face
               hS             = acc_width(i,jso,k+1)-acc_width(i,jso,k)   ! width at South face
               area_E(i,j,k)  = facelen_E(i,j)*(h0+hE)/2.0
               area_S(i,j,k)  = facelen_S(i,j)*(h0+hS)/2.0
               cellvol(i,j,k) = facelen_E(i,j)*facelen_S(i,j)*h0 ! need angle, if coordinate vectors non orthogonal
            enddo
         enddo
      enddo
c 
c     Update cell face diffusivities from cell-centered horizontal/vertical diffusivities:
c     Assume hdiffus/vdiffus are cell centered
c       ku(i,j,k) : diffusivity [m2/s] at E face of cell (i,j,k)      <- hdiffus     
c       kv(i,j,k) : diffusivity [m2/s] at S face of cell (i,j,k)      <- hdiffus     
c       kw(i,j,k) : diffusivity [m2/s] at upper face of cell (i,j,k)  <- vdiffus      
c   
      do k=1,nz
         do j=1,ny
            do i=1,nx
               ku(i,j,k) = 0.5 * (hdiffus(i,j,k) +
     +                            hdiffus(min(nx,i+1),j,k))
               kv(i,j,k) = 0.5 * (hdiffus(i,j,k) + 
     +                            hdiffus(i,max(1,j-1),k))
               kw(i,j,k) = 0.5 * (vdiffus(i,j,k) + 
     +                            vdiffus(i,j,max(1,k-1)))
            enddo
         enddo
      enddo
c
      call update_eulerian_production_rate()  ! production_rate auxillaries
c
      end subroutine update_eulerian_auxillaries



      subroutine setup_fcc_gradients(eulerf, dcdx, dcdy, dcdz)
c     ------------------------------------------------------------------
c     Evaluate gradients [X/m4] at cell faces of Eulerian field with unit X
c     gradients staggered as current vectors
c     For simplicity, do not consider topography, but derive gradients for
c     entire sub grid domain, since Eulerian field is defined all over
c     the sub grid (possibly a pad value)
c     Notice declaration of dummies dcdx, dcdy, dcdz, which conserves index offset == 0
c     ------------------------------------------------------------------
      type(eulerian_field),intent(in) :: eulerf
      real,intent(out)                :: dcdx(0:,0:,0:)     ! x conc gradient, staggered as u
      real,intent(out)                :: dcdy(0:,0:,0:)     ! y conc gradient, staggered as v
      real,intent(out)                :: dcdz(0:,0:,0:)     ! z conc gradient, staggered as w
      integer                         :: i,j,k, nxs,nys,nzs, ix0,iy0
c     ------------------------------------------------------------------
      call resolve_inner_domain_dims(eulerf,nxs,nys,nzs)
      ix0 = eulerf%ixoff 
      iy0 = eulerf%iyoff
c
c     ------ x gradient at cell East face ------
c
      dcdx = 0.0
      do k=1,nzs
         do j=1,nys
            do i=0,nxs
               dcdx(i,j,k) = (eulerf%data(i+1,j,k)-eulerf%data(i,j,k))
     +                        / dx(i+ix0, j+iy0)
            enddo
         enddo
      enddo
c
c     ------ y gradient at cell South face ------
c
      dcdy = 0.0
      do k=1,nzs
         do j=1,nys+1
            do i=1,nxs
               dcdy(i,j,k) = (eulerf%data(i,j,k)-eulerf%data(i,j-1,k))
     +                        / dy(i+ix0, j+iy0)
            enddo
         enddo
      enddo      
c
c     ------ z gradient at cell roof ------
c
      dcdz = 0.0
      do k=2,nzs   ! no surface/bottom flux
         do j=1,nys
            do i=1,nxs
               dcdz(i,j,k) = (eulerf%data(i,j,k)-eulerf%data(i,j,k-1))
     +                        / dz(i+ix0, j+iy0, k)
            enddo
         enddo
      enddo
c
      end subroutine setup_fcc_gradients




      subroutine evaluate_production_rate(eulerf, prorate)
c     ------------------------------------------------------------------
c     Evaluate production rate of Eulerian field [X/m3/s] at cell centers
c     of Eulerian field with unit X 
c     Only evaluate production in wet cells, since value of physical
c     fields in dry cells need not have a meaningfull value (arbitrary pad value)
c     prorate is only defined in the interior part of sub grid (1:nxs, 1:nys, 1:nzs)
c     Need to generalize of temperature,salinity,depth etc dependence
c     Module eulerian_production_rate  defines production_rate
c     ------------------------------------------------------------------
      type(eulerian_field),intent(in) :: eulerf
      real,intent(out)                :: prorate(:,:,:)  ! indices on sub grid, range = (1:nxs, 1:nys, 1:nzs)

      integer                         :: nxs,nys,nzs,i,j,k
      integer                         :: ixoff,iyoff
      real, pointer                   :: data_interior(:,:,:)
      logical, pointer                :: incl_interior(:,:,:)
c     ------------------------------------------------------------------
      prorate = 0.0 ! pad value for dry/boundary cells
      call resolve_inner_domain_dims(eulerf,nxs,nys,nzs)
      ixoff = eulerf%ixoff 
      iyoff = eulerf%iyoff
      data_interior => eulerf%data(1:nxs, 1:nys, 1:nzs)    ! exclude sub grid boundary
      incl_interior => eulerf%dyncell(1:nxs, 1:nys, 1:nzs) ! exclude sub grid boundary
c  
      call point_production_rate(prorate, data_interior, incl_interior,
     +                           ixoff, iyoff, nxs, nys, nzs,
     +                           area_R, cellvol)


      end subroutine evaluate_production_rate



      subroutine impose_boundary_conditions(eulerf)
c     ------------------------------------------------------------------
c     Impose BC on cushin layers
c     currently only Dirichlet BC are implemented
c     ------------------------------------------------------------------
      type(eulerian_field),intent(inout) :: eulerf
      integer :: nxs,nys,nzs,i,j,k,ix0,iy0
c     ------------------------------------------------------------------
      ix0 = eulerf%ixoff 
      iy0 = eulerf%iyoff
      call resolve_inner_domain_dims(eulerf,nxs,nys,nzs) ! of sub grid
c     ---- apply BC to horizontal cushins ----
      do k=1,nzs
         eulerf%data(0,:,k)     = eulerf%dataBC(k)
         eulerf%data(nxs+1,:,k) = eulerf%dataBC(k)
         eulerf%data(:,0,k)     = eulerf%dataBC(k)
         eulerf%data(:,nys+1,k) = eulerf%dataBC(k)
      enddo
c     ---- apply BC below sea bed ----
      do j=1,nys
         do i=1,nxs
            do k = bottom_layer(i+ix0,j+iy0)+1, nzs
               eulerf%data(i,j,k) = eulerf%dataBC(k)
            enddo 
         enddo
      enddo    
c
      end subroutine impose_boundary_conditions



      subroutine update_eulerian_field(eulerf, dt, status)
c     ------------------------------------------------------------------
c     Propagate forward Eulerian field eulerf by time step dt(seconds)
c
c     Do not advance master clock, as several tasks may run
c     in parallel, e.g. multiple Eulerian fields.
c     update_cell_characteristics must have been invoked prior to calling this subroutine
c     Currently no error conditions are flagged (return status = 0)
c
c     This implementation is based on the following staggering:
c
c        u,v positive along lambda,phi in units m/s
c        w   positive downward         in units m/s  
c        u(ix,iy,iz) in grid position (ix+0.5, iy    , iz)     (i.e. eastern  cell face)
c        v(ix,iy,iz) in grid position (ix    , iy-0.5, iz)     (i.e. southern cell face)
c        w(ix,iy,iz) in grid position (ix    , iy,     iz-0.5) (i.e. upper    cell face)
c    
c     May 24 2012: optional divergence correction, correcting for 
c                  loaded water fluxes into control volumes not being exactly zero
c     ------------------------------------------------------------------
      type(eulerian_field),intent(inout) :: eulerf
      real,intent(in)                    :: dt
      integer, intent(out)               :: status
      integer              :: nxs,nys,nzs             ! subgrid shape
      integer              :: ix0,iy0,iz0,ix1,iy1,iz1 ! subgrid delimiters
      integer              :: i,j,k
      real                 :: mflux, trvec, curr, wflux, c0

c     .... temporary arrays ....
      
      real,allocatable     :: prorate(:,:,:)  ! cell production rate, cell centered
      real,allocatable     :: fu(:,:,:)       ! x face mass flux, staggered as u
      real,allocatable     :: fv(:,:,:)       ! y face mass flux, staggered as v
      real,allocatable     :: fw(:,:,:)       ! z face mass flux, staggered as w
      real,allocatable     :: fu0(:,:,:)      ! x face water flux, staggered as u
      real,allocatable     :: fv0(:,:,:)      ! y face water flux, staggered as v
      real,allocatable     :: fw0(:,:,:)      ! z face water flux, staggered as w
      real,allocatable     :: dcdx(:,:,:)     ! x conc gradient, staggered as u
      real,allocatable     :: dcdy(:,:,:)     ! y conc gradient, staggered as v
      real,allocatable     :: dcdz(:,:,:)     ! z conc gradient, staggered as w
      real,allocatable     :: dcdt(:,:,:)     ! conc increment rate
c     ------------------------------------------------------------------
c
c     ------ resolve sub grid size / grid delimiters ------
c
      call resolve_inner_domain_dims(eulerf,nxs,nys,nzs)
      ix0 = eulerf%ixoff 
      iy0 = eulerf%iyoff
      iz0 = 0     ! currrently no vertical offset allowed
c
c     ------ allocate/init auxillary grids ------
c
      allocate ( prorate(nxs, nys, nzs) )       ! only interior domain
      allocate ( dcdt   (nxs, nys, nzs) )       ! only interior domain
      allocate ( fu     (0:nxs+1, 0:nys+1, 0:nzs+1) )
      allocate ( fv     (0:nxs+1, 0:nys+1, 0:nzs+1) )
      allocate ( fw     (0:nxs+1, 0:nys+1, 0:nzs+1) )
      allocate ( fu0    (0:nxs+1, 0:nys+1, 0:nzs+1) )
      allocate ( fv0    (0:nxs+1, 0:nys+1, 0:nzs+1) )
      allocate ( fw0    (0:nxs+1, 0:nys+1, 0:nzs+1) )
      allocate ( dcdx   (0:nxs+1, 0:nys+1, 0:nzs+1) )
      allocate ( dcdy   (0:nxs+1, 0:nys+1, 0:nzs+1) )
      allocate ( dcdz   (0:nxs+1, 0:nys+1, 0:nzs+1) )
      fu = 0.0
      fv = 0.0
      fw = 0.0
      fu0 = 0.0   ! keeps this value if div_corr is .false.
      fv0 = 0.0   ! keeps this value if div_corr is .false.
      fw0 = 0.0   ! keeps this value if div_corr is .false.
      dcdx = 0.0
      dcdy = 0.0
      dcdz = 0.0
c
c     ------------  setup cell face fluxes ------------
c
c     assign fu(0:nxs, 1:nys,   1:nzs)   mass flux over EW faces 
c            fv(1:nxs, 1:nys+1, 1:nzs)   mass flux over NS faces 
c            fw(1:nxs, 1:nys,   1:nzs+1) mass flux over upper vertical faces 
c     
      call setup_fcc_gradients(eulerf, dcdx, dcdy, dcdz)
      call evaluate_production_rate(eulerf, prorate)
c  
      do k=1,nzs
         do j=1,nys
            do i=0,nxs
               if (eulerf%fluxE(i,j,k)) then
                  curr = u(i+ix0, j+iy0, k+iz0)
                  if (curr > 0) then
                     trvec = curr*eulerf%data(i,j,k) ! upwind
                  else
                     trvec = curr*eulerf%data(i+1,j,k) ! upwind
                  endif
                  fu(i,j,k) = area_E(i+ix0, j+iy0, k+iz0)*
     +                 (trvec - ku(i+ix0, j+iy0, k+iz0) * dcdx(i,j,k)) ! mass flux over EW faces 
                  if (div_corr) then
                     fu0(i,j,k) = curr * area_E(i+ix0, j+iy0, k+iz0)   ! water flux over EW faces 
                  endif
               endif
            enddo
         enddo
      enddo

      do k=1,nzs
         do j=1,nys+1
            do i=1,nxs
               if (eulerf%fluxS(i,j,k)) then
                  curr = v(i+ix0, j+iy0, k+iz0)
                  if (curr > 0) then
                     trvec = curr*eulerf%data(i,j-1,k) ! upwind
                  else
                     trvec = curr*eulerf%data(i,j,k)   ! upwind
                  endif
                  fv(i,j,k) = area_S(i+ix0, j+iy0, k+iz0)*
     +                 (trvec - kv(i+ix0, j+iy0, k+iz0) * dcdy(i,j,k)) ! mass flux over NS faces 
                  if (div_corr) then
                     fv0(i,j,k) = curr * area_S(i+ix0, j+iy0, k+iz0)   ! water flux over NS faces 
                  endif
               endif
            enddo
         enddo
      enddo
c
c     k loop delimiters: no vertical flux through surface/grid bottom
c     notice that grid bottom is (usually) below the sea bed
c
      do k=2,nzs  ! render fw(:,:,1) = fw(:,:,nzs+1) = 0
         do j=1,nys
            do i=1,nxs
               if (eulerf%fluxR(i,j,k)) then
                  curr = w(i+ix0, j+iy0, k+iz0)
                  if (curr > 0) then
                     trvec = curr*eulerf%data(i,j,k-1) ! upwind
                  else
                     trvec = curr*eulerf%data(i,j,k)   ! upwind
                  endif
                  fw(i,j,k) = area_R(i+ix0, j+iy0)*
     +                 (trvec - kw(i+ix0, j+iy0, k+iz0) * dcdz(i,j,k)) !  mass flux over vertical faces 
                  if (div_corr) then
                     fw0(i,j,k) = curr * area_R(i+ix0, j+iy0)          ! water flux over vertical faces 
                  endif
               endif
            enddo
         enddo
      enddo
c
c     ------ remove vertical fluxes through sea bed
c
      do j=1,nys
         do i=1,nxs
            if (wetmask(i+ix0,j+iy0) == 1) then
               fw(i, j, bottom_layer(i+ix0,j+iy0)+1) = 0 ! no vertical flux through bottom
               if (div_corr) then
                     fw0(i,j,bottom_layer(i+ix0,j+iy0)+1) = 0 ! no vertical flux through bottom 
               endif
            endif
         enddo
      enddo         
c
c     ------------  assemble concentration rate and update fields ------------
c
      wflux = 0.0 ! default, if div_corr is .false.
      c0    = 0.0
      do k=1,nzs
         do j=1,nys
            do i=1,nxs
               mflux = - fu(i,j,k) + fu(i-1,j,k)
     +                 + fv(i,j,k) - fv(i,j+1,k)
     +                 + fw(i,j,k) - fv(i,j,k+1)
               if (div_corr) then  ! 
                  wflux = - fu0(i,j,k) + fu0(i-1,j,k)
     +                    + fv0(i,j,k) - fv0(i,j+1,k)
     +                    + fw0(i,j,k) - fv0(i,j,k+1)
                  c0          = eulerf%data(i,j,k)
               endif
               
               dcdt(i,j,k) = (mflux - c0*wflux) / 
     +                       cellvol(i+ix0, j+iy0, k+iz0) +
     +                       prorate(i,j,k)
            enddo
         enddo
      enddo

c     -------- debug section: begin --------
c      i=1; j=1; k=1
c      write(*,*) ">"
c      do i=1,nxs
c         do j=1,nys
c            do k=1,bottom_layer(i+ix0,j+iy0)
c               if (cellvol(i+ix0, j+iy0, k+iz0)<1) then
c                  write(*,*) i,j,k,cellvol(i+ix0, j+iy0, k+iz0)
c               endif
c            enddo
c         enddo
c      enddo
c      write(*,*) "<"
c      write(*,*) bottom_layer(i+ix0,j+iy0)
c      write(*,*) (cellvol(i+ix0, j+iy0, k+iz0), k=1, nzs)
c
c      stop  888
c     -------- debug section: end   --------

c
c     ------- update  inner domain and impose boundary conditions -------
c
      eulerf%data(1:nxs,1:nys,1:nzs) = eulerf%data(1:nxs,1:nys,1:nzs) + 
     +                                 dt*dcdt  ! vector assignment for inner domain
      call impose_boundary_conditions(eulerf)
c
c     ------ deallocate auxillary grids ------
c
      deallocate ( prorate )
      deallocate ( dcdt    )
      deallocate ( fu      )
      deallocate ( fv      )
      deallocate ( fw      )
      deallocate ( dcdx    ) 
      deallocate ( dcdy    )
      deallocate ( dcdz    )
c
      status = 0   ! no conditions to flag
c
      end subroutine update_eulerian_field

ccccccccccccccccc specific interface (for this implementation) ccccccccccccccc

      
      subroutine dump_2Dframe(eulerf, ilay, fname)
c     ------------------------------------------------------------------
c     Dump layer number ilay of eulerf%data to file fname
c     do not check for wet/dry
c     Mainly for debugging
c     ------------------------------------------------------------------
      type(eulerian_field),intent(in) :: eulerf
      integer,intent(in)              :: ilay
      character*(*),intent(in)        :: fname   

      integer           :: nxs,nys,nzs
      integer           :: ix,iy,iunit
      character*(*),parameter :: fmt = '(9999(f8.3,1x))'
c     ------------------------------------------------------------------
      call find_free_IO_unit(iunit) ! runtime_tools.o
      open(iunit,file=fname) ! assume unused
      call resolve_inner_domain_dims(eulerf,nxs,nys,nzs)
c     --- skip boundary layers ---
      do iy=1,nys
         write(iunit,fmt) (eulerf%data(ix,iy,ilay), ix=1,nxs)
      enddo
      close(iunit)
      end subroutine dump_2Dframe


      function get_data(eulerf) 
c     -----------------------------------------------------
c     index offsets of eulerf%data are passed to get_data
c     -----------------------------------------------------
      type(eulerian_field),intent(inout) :: eulerf
      real,pointer :: get_data(:,:,:)
      get_data => eulerf%data
      end function get_data


      end module
