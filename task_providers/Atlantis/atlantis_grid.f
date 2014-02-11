      module atlantis_grid
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Represent an Atlantis model grid and provide functionality for extracting 
c     hydrodynamic data for Atlantis input
c     
c     Nearly minimal implementation.
c     An Atlantis model grid is horizontally is a set of polygons + adjoining faces 
c     Vertically, layers are vertically alligned between boxes (so the Nth layer
c     from the top has same width in all boxes), but layer widths are flexible.
c     Especially, the lowest layer adjecent to the sea bed can not have a flexible width in each box, 
c     but must follow the common vertical structure.
c     The number of wet layers per box is flexible, allowing for a rough representation of actual topography.
c     Atlantis layers are indexed 0,1,... from the bottom and up. In the internal representation in this module
c     we flip indexing, so it is 1,2,... starting from the surface and then at output time
c     flip to Atlantis indexing.
c     Positive flux direction is up (?).
c     The module can represent only one Atlantis grid (a singleton)
c     NOTE: the BGM file must be preprocessed to places all whitespaces (except "\n") with " "
c           using e.g. the script regularize_blanks.py
c     
c     assumptions: land mask fixed (i.e. no fluctuating coast line)
c                  allow SSH fluctuatrions
c     prepare: flux residual minimization
c     Atlantis boundary box number == 0 
c
c     Basic validation: horizontal/vertical fluxes:       OK, 30 Jan 2014
c                       box avg temperatures/salinities:  OK, 11 Feb 2014
c
c     TODO: latitude dependent normal vectors along faces in type AtlantisFace
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use input_parser
      use geometry              ! get_horizontal_distance
      use constants             ! deg2rad
      use polygons              ! lonlat_polygon + methods
      use time_tools
      use physical_fields
      use netcdf

      implicit none
      private

      type AtlantisBox
c     -------------------------------------------------------------------
c     Represent a vertical column of Atlantis polygons, given by its 
c     horizontal coordinates. Each polygon may have an arbitrary number of vertices.
c
c     Provide integral sampling meshes + weights
c     Horizontal integral sampling mesh is derived from a regular lon-lat grid and truncated by the polygon
c     and only hold points inside the polygon
c     Support surface integrals (vertical flux) and volume integrals (e.g. T,S)
c     -------------------------------------------------------------------
      type(lonlat_polygon) :: shape             ! geometric information
      integer, pointer     :: faces(:)          ! surrounding (dynamic) face numbers, referring to list Faces
      integer, pointer     :: neighbors(:)      ! neighboring (dynamic) box numbers , referring to list Boxes
      integer              :: nwet              ! number of wet polygons in this column. Layers numbered 1 ... from surface and down
      real, pointer        :: xy(:,:)           ! input (2, 1:nintgpts) wet horizontal integral mesh inside polygon
      real, pointer        :: da(:)             ! (1:nintgpts) area associated each point in xy
      end type
      


      type AtlantisFace
c     -------------------------------------------------------------------
c     Represent a vertical column linear oriented faces, which divides polygons in
c     the Atlantis grid. Only faces between dynamic polygons are represented as faces
c     
c     Support surface integrals and surface flux integrals
c     Hold topological information passed in bgm file
c     Provide integral sampling meshes + weights
c     If a more accurate representation is needed (or the face is very long latitude wise)
c     a normal vector should be associated with each xy point
c     -------------------------------------------------------------------
      real             :: nodes(2,2)        ! first and last point horizontal point of face
      real             :: normal(3)         ! defines positive flux direction through this face (have 3 elemnts to allow vector operations)
      integer          :: left_right(2)     ! neighboring box numbers, referring to list Boxes
      real, pointer    :: xy(:,:)           ! input (2, 1:nintgpts) wet horizontal integral mesh inside polygon
      real, pointer    :: dl(:)             ! (1:nintgpts) length associated with each point in xy
      end type


c     ------------------------------------------------------------------------------------------ 
c                                    Module data section
c     ------------------------------------------------------------------------------------------ 
c     Boxes and Faces are referred to index 0 to allign with Atlantis indexing
c     only dynamics boxes (with model state variables) are listed in Boxes(:)
c     and only faces between dynamics boxes are listed in Faces(:). i.e. faces to
c     model domain are not listed here.
c     Currently do not load boundary vertices ("bnd_vert"), as we do not need them in this context
c     
      type(AtlantisBox),  allocatable, target :: Boxes(:) ! 0: horizontal dynamic polygons in the Atlantis model setup
      type(AtlantisFace), allocatable, target :: Faces(:) ! 0: horizontal dynamic faces in the Atlantis model setup

c     ------------------------------------------------------------------------------------------ 
c     In the Atlantis grid, the vertical structure is identical between horizontal polygons 
c     Internally in this module layers are numbered 1,2 ... from the surface and down; at output
c     time, layer indexing is cast into Atlantis ordering, i.e. 0,1,2 ... from the bottom and up;
c     ------------------------------------------------------------------------------------------ 

      real, allocatable :: zgrid_bounds(:)  ! layer bounds for layer m is zgrid_bounds(m, m+1) (1 : nlayers+1)
      real, allocatable :: zgrid_CC(:)      ! depth of cell center in each layer               (1 : nlayers)
c
c     ------------ intra layer sampling: a linear grid in each layer  --------------
c                  intra layer sampling point m in layer n has depth
c                  z = zoffset(n) + (m-1)*dz(n)    m = 1 ... nzpts(n), n = 1 ... nlayers
c 
      real, allocatable    :: dz(:)         ! vertical intra layer spacings in each layer         (1 : nlayers)
      integer, allocatable :: nzpts(:)      ! number of intra layer sampling points in each layer (1 : nlayers)
      real, allocatable    :: zoffset(:)    ! start of intra layer sampling for each layer        (1 : nlayers)  

c     -------------------- misc module data ---------------------------------------------------
      integer, parameter :: taglength = 64    ! string length for identifiers
      integer            :: verbose = 1

c     -------------------- NetCDF aux data ---------------------------------------------------
c     1) hydro data set    
c
      integer            :: ncid_hydro        ! netCDF ID of hydro data set    
      integer            :: t_dimid_hydro     ! dimid for dimension t (time) in hydro data set        
      integer            :: b_dimid_hydro     ! dimid for dimension b (boxes) in hydro data set    
      integer            :: z_dimid_hydro     ! dimid for dimension z (layers) in hydro data set    
      integer            :: dest_dimid_hydro  ! dimid for dimension dest (flux destinations) in hydro data set     
      integer            :: t_id_hydro        ! varid for variable t (time) in hydro data set     
      integer            :: exchange_id_hydro ! varid for variable exchange in hydro data set     
      integer            :: dest_b_id_hydro   ! varid for variable dest_b (box dest map) in hydro data set     
      integer            :: dest_k_id_hydro   ! varid for variable dest_k (layer dest map) in hydro data set  
      integer, allocatable :: dest_b(:,:,:)       ! box destination map, static during run
      integer, allocatable :: dest_k(:,:,:)       ! layer destination map, static during run
      integer, allocatable :: h_destmap(:,:,:)    ! (3,nlayers,nfaces) destination lookup for hflux
      integer, allocatable :: v_destmap(:,:,:)    ! (3,nlayers,nboxes) destination lookup for vflux
      integer              :: dest            ! actual maximum number of connections for this topology
      real(kind=8), parameter  :: hydro_fill = 0.d0  ! define as kind=8, otherwise writing problems with other vars ...
      integer, parameter       :: dest_fill  = -1
c
c     2) temperature data set   
c
      integer            :: ncid_temp        ! netCDF ID of temp data set    
      integer            :: t_dimid_temp     ! dimid for dimension t (time) in temp data set        
      integer            :: b_dimid_temp     ! dimid for dimension b (boxes) in temp data set    
      integer            :: z_dimid_temp     ! dimid for dimension z (layers) in temp data set    
      integer            :: t_id_temp        ! varid for variable t (time) in temp data set  
      integer            :: temperature_id_temp ! varid for variable temperature in temp data set  
      real, allocatable  :: temperature(:,:,:)  ! output buffer for temperature
      real(kind=8), parameter :: temp_fill = 0.d0  ! define as kind=8, otherwise writing problems with other vars ...
c
c     3) salinity data set   
c
      integer            :: ncid_salt        ! netCDF ID of salt data set    
      integer            :: t_dimid_salt     ! dimid for dimension t (time) in salt data set        
      integer            :: b_dimid_salt     ! dimid for dimension b (boxes) in salt data set    
      integer            :: z_dimid_salt     ! dimid for dimension z (layers) in salt data set    
      integer            :: t_id_salt        ! varid for variable t (time) in salt data set  
      integer            :: salinity_id_salt ! varid for variable salinity in salt data set  
      real, allocatable  :: salinity(:,:,:)  ! output buffer for salinity
      real(kind=8), parameter  :: salt_fill = 0.d0  ! define as kind=8, otherwise writing problems with other vars ...


c     ---------------- define the public scope -------------------- 

      public :: initialize_atlantis_grid
      public :: close_atlantis_grid  
      public :: make_input_files

c     =========================================================================================
                                        contains
c     =========================================================================================
      


      subroutine initialize_atlantis_grid(bgmfname, lw, numpts, drhoz)
c     --------------------------------------------------------------------
c     Initialize Atlantis grid representaion in this module from a
c     BGM formatted file bgmfname
c     Integer vector lw is a vector of layer widths in meters, counted from the surface and down
c     numpts is the number of samplings per layer for vertical sampling grid
c     drhoz is the length scale for horizontal sampling grids
c     --------------------------------------------------------------------
      character*(*), intent(in) :: bgmfname
      real, intent(in)          :: lw(:) ! meters
      integer, intent(in)       :: numpts
      real, intent(in)          :: drhoz ! meters
      integer                   :: nbox, nface, ibox, iface
      type(control_file)        :: ctrlfile
c     --------------------------------------------------------------------
c
c     -------- initialize vertical structures  --------
c
      call initialize_vertical_grid(lw, numpts)
c
c     -------- initialize horizontal structures --------
c

      call open_control_file(bgmfname, ctrlfile, 
     +                       separator_token = " ", 
     +                       comment_token = "/")        ! more precisely "//" but "/" can not appear in a name/value
      call read_control_data(ctrlfile, "nbox",  nbox)
      call read_control_data(ctrlfile, "nface", nface)
      
c     -------- initialize grid boxes from BGM file --------   
      allocate( Boxes(0 : nbox-1) )  ! apply offset 0 to conform with Atlantis convention
      do ibox = 0, nbox-1
         call load_box_from_bgm_file(ctrlfile, ibox, drhoz)
      enddo
      if (verbose > 0) then
         write(*,467) "initialize_atlantis_grid:", nbox, "boxes"
      endif
c     --------  initialize box faces from BGM file --------    
      allocate( Faces(0 : nface-1) ) ! apply offset 0 to conform with Atlantis convention
      do iface = 0, nface-1
         call load_face_from_bgm_file(ctrlfile, iface, drhoz)
      enddo
      if (verbose > 0) then
         write(*,467) "initialize_atlantis_grid:",nface,"faces"
      endif
 467  format(a,1x,i8,1x,a)

      call close_control_file(ctrlfile)
 
      end subroutine initialize_atlantis_grid
    


      subroutine initialize_vertical_grid(lw, numpts)
c     --------------------------------------------------------------------
c     Initialize vertical grid auxillaries from input.
c     
c     In the Atlantis grid, layers are vertically alligned between boxes (so the Nth layer
c     from the top has same width in all boxes).
c
c     Input
c         lw:     vector of layer widths in meters, counted from the surface and down - this 
c                 also implies the max number wet of layers per box     
c         numpts: samplings per layer (for intra layer integrals) - currently constant
c
c     Set module data:
c         zgrid_bounds(:)  ! layer bounds for layer m is zgrid_bounds(m, m+1) (1 : nlayers+1)
c         zgrid_CC(:)      ! depth of cell center in each layer               (1 : nlayers)
c         dz(:)            ! vertical intra layer spacings in each layer
c         nzpts(:)         ! number of intra layer sampling points in each layer ()
c         zoffset(:)       ! start of intra layer sampling for each layer
c     --------------------------------------------------------------------
      real, intent(in)    :: lw(:)   ! layer widths in meters
      integer, intent(in) :: numpts  ! samplings per layer, currently constant
      integer             :: nlayers, ilay
c     --------------------------------------------------------------------
      if (verbose > 0) then
         write(*,332) "initialize_vertical_grid: layers     =", lw
         write(*,333) "initialize_vertical_grid: sub sampl. =", numpts
 332     format(a,99f8.2)
 333     format(a,i6)
      endif
      nlayers = size(lw)
      allocate( zgrid_bounds(nlayers+1) )
      allocate( zgrid_CC(nlayers)       )
      allocate( dz(nlayers)             )
      allocate( nzpts(nlayers)          )     
      allocate( zoffset(nlayers)        )
c
      zgrid_CC(1) = 0.5*lw(1)
      do ilay = 2, nlayers
         zgrid_CC(ilay) = zgrid_CC(ilay-1) + 0.5*(lw(ilay-1) + lw(ilay))
      enddo
c
      zgrid_bounds(1) = 0.0    ! sea surface = zero
      do ilay = 1, nlayers
         zgrid_bounds(ilay+1) = zgrid_bounds(ilay) + lw(ilay)
      enddo
c
c     ------------ intra layer sampling: a linear grid in each layer  --------------
c                  intra layer sampling point m in layer n has depth
c                  z = zoffset(n) + (m-1)*dz(n)    m = 1 ... nzpts(n)
c 
      dz    = lw/numpts  ! vector assignment
      nzpts = numpts     ! currently same across layers
      do ilay = 1, nlayers
         zoffset(ilay) = zgrid_bounds(ilay) + 0.5*dz(ilay) 
      enddo
      
      end subroutine initialize_vertical_grid




      subroutine make_prop_tag(tag, typ, i, prop)
c     ------------------------------------------------------------------------
c     Generate property tag like typ<i>.prop (e.g. "box2.vertmix")
c     for BGM file reading 
c     ------------------------------------------------------------------------
      character(len=taglength),intent(out) :: tag   ! assume output buffer sufficient
      integer,intent(in)                   :: i     
      character*(*),intent(in)             :: typ, prop 
      character(len=taglength)             :: i_as_str
c
      write(i_as_str, *) i
      tag = trim(adjustl(typ)) // trim(adjustl(i_as_str)) // 
     +                     "." // trim(adjustl(prop))
      end subroutine make_prop_tag
      




      subroutine load_box_from_bgm_file(ctrlfile, ibox, drhoz)
c     --------------------------------------------------------------------
c     Initialize box ibox from parameters in ctrlfile
c
c     drhoz[meters] is the isotropic horizontal scale used to initialize integral sampling grid
c     Assume Boxes has been allocated. 
c     Only retain a necessary sub set of the informatoin provided in BGM file 
c     Use parameter boxN.botz to statically select number of wet layers in each box
c     Parameter format in BGM file (example for box 2)
c         # Box number 2
c         box2.label      Box2
c         box2.inside     11.173199238623864 56.30559340280887
c         box2.nconn      3
c         box2.iface      7 8 99 
c         box2.ibox       27 3 4 
c         box2.botz       -33.0
c         box2.area       1.261383763075706
c         box2.vertmix    9.99999997475E-7
c         box2.horizmix   1.0
c         box2.vert 11.24336228544731 55.80020684195673
c         box2.vert 12.127242435108569 55.66716261584352
c         box2.vert 10.428889554898678 57.586861347572196
c         box2.vert 10.219156023058297 56.75538033731939
c         box2.vert 11.24336228544731 55.80020684195673
c     --------------------------------------------------------------------
      type(control_file), intent(in)  :: ctrlfile
      integer, intent(in)             :: ibox
      real, intent(in)                :: drhoz  ! unit = meters

      character(len=taglength)        :: proptag
      character(len=12)               :: boxname
      real                            :: botz, geo(2), bbox(4)
      real                            :: nw(2), se(2), sw(2), jacob(3)
      real                            :: drlon, drlat, dlon, dlat
      integer                         :: nnodes, inode, inext 
      integer                         :: nxsub, nysub, naccept, ix, iy
      integer                         :: nconn, ilay, maxlayers, i
      real, allocatable               :: nodebuf(:,:), xybuf(:,:)
c     --------------------------------------------------------------------
      call make_prop_tag(proptag, "box", ibox, "nconn")
      call read_control_data(ctrlfile, proptag, nconn)
c
      allocate(Boxes(ibox)%neighbors(nconn))
      call make_prop_tag(proptag, "box", ibox, "ibox")
      call read_control_data(ctrlfile, proptag, Boxes(ibox)%neighbors)
c
      allocate(Boxes(ibox)%faces(nconn))
      call make_prop_tag(proptag, "box", ibox, "iface")
      call read_control_data(ctrlfile, proptag, Boxes(ibox)%faces)
c
c     Use parameter boxN.botz to statically select number of wet layers in each box
c     layer m: zgrid_bounds(m, m+1)
c
      call make_prop_tag(proptag, "box", ibox, "botz")
      call read_control_data(ctrlfile, proptag, botz)
      botz      = -botz    ! now botz > 0
      maxlayers = size(zgrid_CC)
      do ilay = 1, maxlayers
         if (zgrid_bounds(ilay+1) > botz) exit
      enddo
      if (ilay > maxlayers) then ! standard exit condition at loop bound
         write(*,*) "load_box_from_bgm_file: botz > zgrid_bounds"
         write(*,*) "box number = ", ibox
         write(*,*) "botz       = ", botz
         write(*,*) zgrid_bounds(1), " < zgrid < ", 
     +              zgrid_bounds(maxlayers+1)
         stop 
      else ! loop terminated by break
         Boxes(ibox)%nwet = ilay
      endif
c  
c     load box perimeter
c
      call make_prop_tag(proptag, "box", ibox, "vert")
      nnodes = count_tags(ctrlfile, proptag) - 1    ! first point repeated at last
      if (nnodes < 3) then ! triangle limit 
         write(*,*) "load_box_from_bgm_file: nnodes =", nnodes, "< 3"
         write(*,*) "box number = ", ibox
         stop 
      else ! loop terminated by break
         allocate ( nodebuf(2,nnodes) )
         inext = 1
         do inode = 1, nnodes ! skip repeated last node    
            call read_control_data(ctrlfile, proptag, 
     +                             nodebuf(:,inode), inext) ! stepping through with inext 
            inext = inext + 1
         enddo
         write(boxname,*) ibox 
         boxname =  "box" // trim(adjustl(boxname))
         call init_lonlat_polygon(Boxes(ibox)%shape, nodebuf, boxname) ! copies nodebuf
         deallocate ( nodebuf )
      endif
c     
c     Use bounding box to define integral sampling mesh (lon-lat sub grid)
c
      call get_bounding_box(Boxes(ibox)%shape, bbox)   ! bbox  = (SWlon,SWlat,NElon,NElat)
      sw    = bbox(1:2)
      nw(1) = bbox(1)
      nw(2) = bbox(4)
      se(1) = bbox(3)
      se(2) = bbox(2)    
      call get_horizontal_distance(sw, se, drlon) ! width of bounding box in meters
      call get_horizontal_distance(sw, nw, drlat) ! height of bounding box in meters
      nxsub = 1 + int(drlon/drhoz)
      nysub = 1 + int(drlat/drhoz)
      dlon = (bbox(3)-bbox(1))/nxsub ! longitute sampling step
      dlat = (bbox(4)-bbox(2))/nysub ! latitude sampling step      
      allocate( xybuf(2,nxsub*nysub) )
      naccept = 0
      do ix = 1, nxsub
         geo(1) = bbox(1) + (ix-0.5)*dlon
         do iy = 1, nysub
            geo(2) = bbox(2) + (iy-0.5)*dlat
            if (is_inside_polygon(Boxes(ibox)%shape, geo)) then
               naccept = naccept + 1
               xybuf(:,naccept) = geo           
            endif        
         enddo
      enddo
      
c     -------- transfer sampling points from buffer and set area element -------- 
      if (naccept > 0) then
         allocate( Boxes(ibox)%xy(2,naccept) )
         allocate( Boxes(ibox)%da(naccept)   )
         Boxes(ibox)%xy = xybuf(:,1:naccept) ! copy buffer
         do i = 1, naccept
            call get_xyz2cart_jacobian(xybuf(:,i),jacob)
            Boxes(ibox)%da(i) = jacob(1)*dlon*jacob(2)*dlat ! apply tangent space area
         enddo      
      else
         write(*,*) "load_box_from_bgm_file: no sub grid"
         write(*,*) "box number = ", ibox
         stop 
      endif
      deallocate( xybuf )
      
      end subroutine load_box_from_bgm_file





      subroutine load_face_from_bgm_file(ctrlfile, iface, drhoz)
c     --------------------------------------------------------------------
c     Initialize face iface from parameters in ctrlfile
c
c     drhoz[meters] is the isotropic horizontal scale used to initialize integral sampling grid
c     Assume Faces has been allocated. 
c     Assume Boxes has been allocated and initialized (access number of wet layers). 
c     Only retain a necessary sub set of the informatoin provided in BGM file 
c     Parameter format in BGM file (example for box 2)
c         # Face number 4
c         face4.p1 22.59791297912599 58.23918838500983
c         face4.p2 22.133067554897845 58.0461552937827
c         face4.length 0.5033319409043573
c         face4.cs 0.9235365103055799 0.3835105137184525
c         face4.lr 8 28                          // 8 to left, 28 to right when looking from p1 to p2
c
c     According to ... positive flux direction is right-to-left across 
c     a face, i.e. vector from p1 to p2 rotated 90 degrees counter clockwise CHECK
c     
c     --------------------------------------------------------------------
      type(control_file), intent(in)  :: ctrlfile
      integer, intent(in)             :: iface
      real, intent(in)                :: drhoz  ! unit = meters
      character(len=taglength)        :: proptag
      real                            :: nvec(3), dist, p1p2(2), p1(2)
      real                            :: p2(2), dlonlat(2), jacob(3)
      real                            :: dr(2), tmp, phimid
      integer                         :: nsmpl, ipt
c     --------------------------------------------------------------------
      call make_prop_tag(proptag, "face", iface, "p1")
      call read_control_data(ctrlfile, proptag, p1)
      call make_prop_tag(proptag, "face", iface, "p2")
      call read_control_data(ctrlfile, proptag, p2)
      call make_prop_tag(proptag, "face", iface, "lr")
      call read_control_data(ctrlfile, proptag, Faces(iface)%left_right)
      Faces(iface)%nodes(:,1) = p1
      Faces(iface)%nodes(:,2) = p2
      p1p2                    = p2-p1
      phimid                  = 0.5*(p1(2)+p2(2)) ! latitude mid point on face
      nvec(1:2) = p1p2
      nvec(1)   = nvec(1)*cos(deg2rad*phimid)     ! just multiply metric lon/lat ratio, earth_radius normalized out
c     -- rotate nvec(1:2) 90 degrees counter clockwise
      tmp       = nvec(2)
      nvec(2)   = nvec(1)
      nvec(1)   = -tmp
      nvec(3)   = 0.0                                   ! normal vector is horizontal
      Faces(iface)%normal = nvec/sqrt(sum(nvec*nvec))   ! normalize vector
c
      call get_horizontal_distance(Faces(iface)%nodes(:,1), 
     +                             Faces(iface)%nodes(:,2), dist) ! width of bounding box in meters
      nsmpl = 1 + int(dist/drhoz)   
      allocate( Faces(iface)%xy(2,nsmpl) )
      allocate( Faces(iface)%dl(nsmpl)   )
      dlonlat = p1p2/nsmpl
c     ---- setup linear sampling mesh
      do ipt = 1, nsmpl
         Faces(iface)%xy(:,ipt) = p1 + (ipt-0.5)*dlonlat
         call get_xyz2cart_jacobian(Faces(iface)%xy(:,ipt),jacob)
         dr = jacob(1:2)*dlonlat
         Faces(iface)%dl(ipt) = sqrt(sum(dr*dr)) ! apply length of tangent space vector
      enddo
c
      end subroutine load_face_from_bgm_file
 

c     ==========================  Box  analyzers  ========================================
              
      subroutine get_box_phys_averages(temp, salt)
c     -------------------------------------------------------
c     Evaluate temperature/salinity averages in Atlantis boxes 
c     at the current time into buffers temp, salt.
c     Use same layout for output buffers (temp, salt) as Atlantis requires for input
c     Layerwise, box quantities are stored a little peculiar:
c     If nwet = 4, nlayers = 7 then
c       Atlantis layer indexing: 0, 1, 2, 3, _, _,  _, _, d
c       Internal layer indexing: 4, 3, 2, 1, _, _,  _, _, d         
c     where d is the demersal layer. In terms of physics, we
c     probe to quantity at the bottom and associate this with 
c     layer d. "_" is the fill value. Notice that layers written is
c     nlayers + 1, the last being sediment
c     -------------------------------------------------------
      real, intent(out)     :: temp(0:,0:) ! [deg C] shape = (0:nlaýers, 0:nboxes-1) = Atlantis layout
      real, intent(out)     :: salt(0:,0:) ! [PSU]   shape = (0:nlaýers, 0:nboxes-1) = Atlantis layout
      integer               :: ibox, ilay, atlay, istatt, istats
      integer               :: nboxes, nlayers, iz, ixy, istat
      real                  :: t_acc, s_acc, vol_acc, area_acc
      real                  :: geo(3), t,s,v,a,cwd
      type(AtlantisBox), pointer  :: box
c     -------------------------------------------------------
      nboxes  = size(Boxes)
      nlayers = size(zgrid_CC) 
      temp    = temp_fill   ! pad value
      salt    = salt_fill   ! pad value
      do ibox = 0, nboxes-1
         box => Boxes(ibox)
c           ---- average over each layer for this box ---- 
         do ilay = 1, box%nwet       ! ilay refers to internal layer indexing in this module          
            atlay = box%nwet - ilay ! corresponding Atlantis layer index
            t_acc   = 0.0
            s_acc   = 0.0
            vol_acc = 0.0  ! accepted volume
            do iz = 1, nzpts(ilay) ! intra layer loop
               geo(3) = zoffset(ilay) + (iz-1)*dz(ilay)    
               do ixy = 1, size(box%xy, 2)  ! horizontal mesh loop
                  geo(1:2) = box%xy(:,ixy)
                  call interpolate_temp(geo, t, istatt)
                  call interpolate_salty(geo, s, istats)
                  if ((istatt == 0).and.(istats == 0)) then !   assume same, no padding if /= 0 
                     v       = box%da(ixy)*dz(ilay)   ! volume associated with this sampling
                     t_acc   = t_acc + v*t
                     s_acc   = s_acc + v*s
                     vol_acc = vol_acc + v
                  endif
               enddo  ! ixy
            enddo     ! iz
c
            if (vol_acc > 1e-6) then
               temp(atlay,ibox) = t_acc/vol_acc
               salt(atlay,ibox) = s_acc/vol_acc
            endif
c
         enddo  !  ilay
c
c        ---- evaluate sediment layer; ignore slope contribution to jacobian of seabed 
c
         t_acc    = 0.0
         s_acc    = 0.0
         area_acc = 0.0  ! accepted area
         do ixy = 1, size(box%xy, 2)  ! horizontal mesh loop
            geo(1:2) = box%xy(:,ixy)
            call interpolate_wdepth(geo, cwd, istat) 
            if (istat /= 0) cycle       ! dry point or other exception, continue with next ixy
            geo(3) = min(cwd, cwd-1e-3) ! hover just over bottum
            call interpolate_temp(geo, t, istatt)
            call interpolate_salty(geo, s, istats)
            if ((istatt == 0).and.(istats == 0)) then !   assume same, no padding if /= 0 
               a        = box%da(ixy) ! area associated with this sampling, ignore slope
               t_acc    = t_acc + a*t
               s_acc    = s_acc + a*s
               area_acc = area_acc + a
            endif
         enddo                  ! ixy
c
         if (area_acc > 1e-6) then
            temp(nlayers,ibox)  = t_acc/area_acc    ! sediment goes into last layer position
            salt(nlayers,ibox)  = s_acc/area_acc    ! sediment goes into last layer position
         endif
      enddo ! ibox

      end subroutine get_box_phys_averages
      
      

c     ==========================  Flux analyzers  ========================================

      subroutine get_horizontal_water_fluxes(hflux)
c     -------------------------------------------------------
c     Evaluate horizontal water fluxes across loaded Faces at the 
c     current time. hflux is an allocated output buffer for results with
c     unit m3/s. Positive direction for fluxes across faces follows from the
c     face normal vector. Apply pad value 0.0 to faces to dry cells
c     Faulty interpolations (whatever reason) are padded with 0
c     -------------------------------------------------------
      real, intent(out)           :: hflux(:,0:) ! [m3/s] shape = (1:nlaýers, 0:nfaces-1)
      integer                     :: iface, ilay, nfaces, nlayers
      integer                     :: iz, ixy, nintgpts, istat
      type(AtlantisFace), pointer :: face
      type(AtlantisBox), pointer  :: left, right
      real                        :: geo(3), uvw(3), intg, area
c     -------------------------------------------------------
      hflux     = hydro_fill     ! default pad
      nfaces    = size(Faces)    ! 0:nfaces-1     
      do iface = 0, nfaces-1
         face  => Faces(iface)
         left  => Boxes(face%left_right(1))
         right => Boxes(face%left_right(2))
c        --- assess how many layers are wet-wet for face
         nlayers   = min(left%nwet, right%nwet)
         nintgpts  = size(face%xy, 2) ! number of horizontal integration points
         do ilay = 1, nlayers   ! layer loop, internal numbering
c           ---- perform surface integral for this layer at this face
            intg = 0.0
            do iz = 1, nzpts(ilay) ! intra layer loop
               geo(3) = zoffset(ilay) + (iz-1)*dz(ilay)    
               do ixy = 1, nintgpts  ! horizontal mesh loop
                  geo(1:2) = face%xy(:,ixy)
                  call interpolate_currents(geo, uvw, istat)
                  if (istat > 0) uvw = 0.0                 ! Faulty interpolations (whatever reason) are padded with 0
                  area = face%dl(ixy)*dz(ilay)            ! area element associated with this sampling
                  intg = intg + sum(face%normal*uvw)*area ! project current vector
               enddo
            enddo
            hflux(ilay, iface) = intg
         enddo  ! ilay
      enddo     ! iface

      end subroutine get_horizontal_water_fluxes



      subroutine get_vertical_water_fluxes(vflux)
c     -------------------------------------------------------
c     Evaluate vertical water fluxes for all Boxes between layers
c     at the current time. vflux is an allocated output buffer for results with
c     unit m3/s. Positive direction for fluxes across faces is upward. 
c     vflux(i,j) is the flux from layer i+1 to layer i (offset 1) in box j (offset 0)
c     Faulty interpolations (whatever reason) are padded with 0
c     Fluxes through bottom of a column: if box%nwet < nlaýers = size(zgrid_CC) 
c     then vflux(box, box%nwet:) = 0 else (if box%nwet = nlaýers) 
c     vflux(box, nlaýers) is calculated
c     -------------------------------------------------------
      real, intent(out)           :: vflux(:,0:) ! [m3/s] shape = (1:nlaýers, 0:nboxes-1)
      integer                     :: ibox, ilay
      integer                     :: ixy, nscan, istat
      integer                     :: nboxes, maxlayers, nintgpts
      type(AtlantisBox), pointer  :: box
      real                        :: geo(3), uvw(3), intg
c     ----------------------------------------------------
      vflux     = hydro_fill     ! default pad
      nboxes    = size(Boxes)    ! 0:nboxes-1
      maxlayers = size(zgrid_CC) 
      do ibox = 0, nboxes-1
         box  => Boxes(ibox)
c        ---- resolve number of layers to scan
         if (box%nwet == maxlayers) then
            nscan = maxlayers   ! allow of water flux through bottom a layer column
         else 
            nscan = box%nwet-1  ! water flux through bottom layer implicitly zero
         endif
         nintgpts  = size(box%xy, 2) ! number of horizontal integration points
         do ilay = 1, nscan          ! layer loop (possibly void)
            intg   = 0.0
            geo(3) = zgrid_bounds(ilay+1) ! bottom of layer ilay
            do ixy = 1, nintgpts  ! horizontal mesh loop
                  geo(1:2) = box%xy(:,ixy)
                  call interpolate_currents(geo, uvw, istat)
                  if (istat > 0) uvw = 0.0                 ! Faulty interpolations (whatever reason) are padded with 0
                  intg = intg - uvw(3)*box%da(ixy)        ! minus: IBMlib fluxes are positive down, Atlantis(?) up 
            enddo
            vflux(ilay, ibox) = intg
         enddo ! ilay
      enddo    ! ibox
 
      end subroutine get_vertical_water_fluxes
      

c     ==========================  basic I/O  ============================================== 

      subroutine print_box(abox,iunit)
c     --------------------------------------------
c     Log AtlantisBox instance abox to unit iunit
c     --------------------------------------------
      type(AtlantisBox), intent(in) :: abox
      integer, intent(in)           :: iunit
      integer                       :: i
c     --------------------------------------------
      write(iunit,'(50("-"))') 
      write(iunit,*) "AtlantisBox instance content"
      write(iunit,*) "shape = "
      call print_lonlat_polygon(abox%shape,iunit)
      write(iunit,*) "faces      = ", abox%faces
      write(iunit,*) "neighbors  = ", abox%neighbors
      write(iunit,*) "wet layers = ", abox%nwet
      write(iunit,*) "sampling points  xy  area[m2]"
      do i=1, size(abox%xy, 2)
         write(iunit,*) i, abox%xy(:,i), abox%da(i)
      enddo
      end subroutine print_box



      subroutine print_face(aface,iunit)
c     --------------------------------------------
c     Log AtlantisFace instance aface to unit iunit
c     --------------------------------------------
      type(AtlantisFace), intent(in) :: aface
      integer, intent(in)            :: iunit
      integer                        :: i
c     --------------------------------------------
      write(iunit,'(50("-"))') 
      write(iunit,*) "AtlantisFace instance content"
      write(iunit,*) "first point    = ", aface%nodes(:,1)
      write(iunit,*) "last  point    = ", aface%nodes(:,2)
      write(iunit,*) "face normal    = ", aface%normal
      write(iunit,*) "left/right box = ", aface%left_right
      write(iunit,*) "sampling points  xy  length[m]"
      do i=1, size(aface%xy, 2)
         write(iunit,*) i, aface%xy(:,i), aface%dl(i)
      enddo
      end subroutine print_face

c     ==========================  netCDF I/O  ============================================== 

      subroutine nf90_check_event(istat, message)
c     ----------------------------------------------------------------------------
c     respond to return status > 0 of netCDF calls
c     ----------------------------------------------------------------------------
      integer, intent(in)       :: istat
      character*(*), intent(in) :: message
c     ----------------------------------------------------------------------------
      if (istat /= NF90_NOERR) then
         write(*,*) "nf90_check_event: netCDF call failed"
         write(*,*) "netCDF error message = ", 
     +              trim(adjustl(NF90_STRERROR(istat)))
         write(*,*) "call context message = ", message
         stop 394
      endif
      end subroutine nf90_check_event


      subroutine initialize_input_files(ncfname_hydro, ncfname_temp, 
     +                                  ncfname_salt, toffset)
c     ----------------------------------------------------------------------------
c     Initialize netCDF files for Atlantis input data and setup static
c     auxillary arrays
c
c     Minimal example: ~/ModellingPackages/Atlantis/ForcingFiles/hydrodummy.cdf
c     Dimension order in arrays: in fortran reverse order compared to CDL def
c     ----------------------------------------------------------------------------
      character*(*), intent(in) :: ncfname_hydro
      character*(*), intent(in) :: ncfname_temp
      character*(*), intent(in) :: ncfname_salt
      type(clock), intent(in)   :: toffset     ! offset time for time vector
      integer                   :: ibox, nboxes, nlayers, nfaces
      integer                   :: istat, ifc, nscan, ilay, idest
      integer                   :: from_box, from_layer, max_dest
      integer                   :: to_box, to_layer
      type(AtlantisFace), pointer :: face
      type(AtlantisBox), pointer  :: left, right, box
      character(len=256)          :: cbuf
      double precision            :: xxx
c     ----------------------------------------------------------------------------
      nboxes  = size(Boxes)
      nlayers = size(zgrid_CC) 
      nfaces  = size(Faces)
c
c     First estimate upper limit for dest (= max_dest):
c           horizontal limit: max(nconn) + one vertical flux + size(Boxes): (bottom layer vertical fluxes connected to box,layer = 0,0, if present)
c     later adjust dest to actual needed size by scaning buffers
c
c     looping over all faces and boxes later determine the actual number of destinations (dest) needed
c     to compress the netCDF set as much as possible. max_dest >= dest
c     
      max_dest = 1
      do ibox = 0, nboxes-1
         max_dest = max(max_dest, size(Boxes(ibox)%neighbors))
      enddo
      max_dest = max_dest + nboxes + 1   
c
c     ======== setup aux buffers for writing a single time frame ========
c
      allocate( dest_b(max_dest, nlayers, nboxes)   )
      allocate( dest_k(max_dest, nlayers, nboxes)   )
      allocate( h_destmap(3, nlayers, 0:nfaces-1)  )  ! index map to easily fill exchange
      allocate( v_destmap(3, nlayers, 0:nboxes-1 ) )  ! index map to easily fill exchange
      dest_b = dest_fill  ! set _FillValue
      dest_k = dest_fill  ! set _FillValue
      dest   = 1   ! dest < max_dest
c
c     -------- scan over all connections --------
c
c     map to Atlantis layer indexing which is bottom-up (and offset 0), opposite the internal which is top down
c     positive flux direction if right-to-left
c
c     First scan horizontal connections and fill destination maps:
c
      do ifc = 0, nfaces-1
         face  => Faces(ifc)
         left  => Boxes(face%left_right(1))
         right => Boxes(face%left_right(2))
         nscan  = min(left%nwet, right%nwet)
         do ilay = 1, nscan  ! top-down
            from_box   = face%left_right(2)   ! right box, flux positive from right to left
            from_layer = right%nwet - ilay    ! Atlantis indexing
            to_box     = face%left_right(1)   ! left box
            to_layer   = left%nwet - ilay     ! Atlantis indexing
c           if ((ifc==14).and.(ilay==1)) then
c           ---- locate next vacant destination (idest) for this source ----
            do idest = 1, max_dest
               if (dest_b(idest, from_layer+1, from_box+1) == dest_fill)
     +             exit  ! 
            enddo
            if (idest > max_dest) then ! loop exhausted 
               write(*,*) "initialize_input_files: "//
     +                    "exhausted dest (face loop)"
               write(*,*) "max_dest=",max_dest
               write(*,*) "from(box,layer)=",from_box,from_layer
               write(*,*) "dest_b=", dest_b(:, from_layer+1, from_box+1)
               write(*,*) "dest_k=", dest_k(:, from_layer+1, from_box+1)
               stop 8449
            else    ! accept connection
               dest = max(idest, dest)                              ! update upper limit for idest
               dest_b(idest, from_layer+1, from_box+1) = to_box     ! +1: array declared 1:
               dest_k(idest, from_layer+1, from_box+1) = to_layer   ! +1: array declared 1:
c              .... hflux(ilay, ifc) -> exchange( h_destmap(:,ilay, ifc) )
               h_destmap(1, ilay, ifc) = idest         ! auxillary map for packing array exchange 
               h_destmap(2, ilay, ifc) = 1+from_layer  ! auxillary map for packing array exchange 
               h_destmap(3, ilay, ifc) = 1+from_box    ! auxillary map for packing array exchange 
            endif
         enddo
      enddo
c
c     Then vertical connections:
c           
      do ibox = 0, nboxes-1
         box => Boxes(ibox)
         if (box%nwet == nlayers) then ! assume bottom layer fluxes, possibly zero
            nscan = nlayers            ! allow of water flux through bottom a layer column
         else 
            nscan = box%nwet-1         ! water flux through bottom layer implicitly zero
         endif
         do ilay = 1, nscan            ! top-down scan
            if (ilay < nlayers) then
               from_box   = ibox
               from_layer = box%nwet - ilay -1  ! layer below this, Atlantis indexing
            else
               from_box   = 0           ! required source box/layer for flux through bottom
               from_layer = 0           ! required source box/layer for flux through bottom  
            endif
            to_box     = ibox              ! this box
            to_layer   = box%nwet - ilay   ! this layer, Atlantis indexing
c
c           ---- locate next vacant destination (idest) ----
c
            do idest = 1, max_dest
               if (dest_b(idest, from_layer+1, from_box+1) == dest_fill) 
     +            exit
            enddo
            if (idest > max_dest) then 
               write(*,*) "initialize_input_files: "//
     +                    "exhausted dest (box loop)"
               write(*,*) "max_dest=",max_dest
               write(*,*) "from(box,layer)=",from_box,from_layer
               write(*,*) "dest_b=", dest_b(:, from_layer+1, from_box+1)
               write(*,*) "dest_k=", dest_k(:, from_layer+1, from_box+1)
               write(*,*)
               stop 8448
            else   ! accept connection
               dest = max(idest, dest)                            ! update upper limit for idest
               dest_b(idest, from_layer+1, from_box+1) = to_box
               dest_k(idest, from_layer+1, from_box+1) = to_layer
c              .... vflux(ilay, ibox) -> exchange( v_destmap(:,ilay, ibox) )
c                   vflux(i,j) is the flux from layer i+1 to layer i (offset 1) in box j (offset 0)
               v_destmap(1, ilay, ibox) = idest         ! auxillary map for packing array exchange 
               v_destmap(2, ilay, ibox) = 1+from_layer  ! auxillary map for packing array exchange 
               v_destmap(3, ilay, ibox) = 1+from_box    ! auxillary map for packing array exchange
            endif
         enddo
      enddo
      if (verbose > 0) write(*,334) dest
 334  format("initialize_input_files: maximum connectivity = ",i6)
c
c     ================= define netCDF sets =================
c
c     1) ------------------ hydro data set ------------------
c
      istat = nf90_create(ncfname_hydro, NF90_CLOBBER, ncid_hydro)   ! overwrite, if exists
      call nf90_check_event(istat, 
     +               "h:creating file:"//trim(adjustl(ncfname_hydro)))
      
      if (verbose > 0) 
     +    write(*,336) "hydro", trim(adjustl(ncfname_hydro))
 336  format("initialize_input_files: creating ",a, 
     +       " input file: ",a)

c     --- create netCDF dimensions ---
      istat =  nf90_def_dim(ncid_hydro, "t", NF90_UNLIMITED,
     +                      t_dimid_hydro)
      call nf90_check_event(istat, "h:creating dimension t")
      istat =  nf90_def_dim(ncid_hydro, "b", nboxes, b_dimid_hydro)
      call nf90_check_event(istat, "h:creating dimension b")
      istat =  nf90_def_dim(ncid_hydro, "z", nlayers, z_dimid_hydro)
      call nf90_check_event(istat, "h:creating dimension z")
      istat =  nf90_def_dim(ncid_hydro, "dest", dest, dest_dimid_hydro)  ! use needed size dest, rather than max_dest
      call nf90_check_event(istat, "h:creating dimension dest")
      istat =  nf90_put_att(ncid_hydro, NF90_GLOBAL, "title", 
     +                  "hydrodynamics from IBMlib")
      call nf90_check_event(istat, "setting title attribute")
c     ----- define variables -----   

      istat =  nf90_def_var(ncid_hydro, "t", NF90_DOUBLE, 
     +                 (/t_dimid_hydro/), 
     +                  t_id_hydro)
      call nf90_check_event(istat, "h:creating variable t")
      cbuf = get_datetime(toffset)
      cbuf = "seconds since "//trim(adjustl(cbuf))
      istat =  nf90_put_att(ncid_hydro, t_id_hydro, "units", 
     +                  trim(adjustl(cbuf)) )
      call nf90_check_event(istat, "h:setting attribute t:units")
      istat =  nf90_put_att(ncid_hydro, t_id_hydro, 
     +                  "_FillValue", 0.0d0)
      call nf90_check_event(istat, "h:setting attribute t:_FillValue")
c
      istat =  nf90_def_var(ncid_hydro, "exchange", NF90_DOUBLE,
     +                 (/dest_dimid_hydro, z_dimid_hydro, 
     +                  b_dimid_hydro, t_dimid_hydro/),  
     +                  exchange_id_hydro)
      call nf90_check_event(istat, "h:creating variable exchange")
      istat =  nf90_put_att(ncid_hydro, exchange_id_hydro, 
     +                  "_FillValue", hydro_fill)
      call nf90_check_event(istat, "h:setting att exchange:_FillValue")
c
      istat =  nf90_def_var(ncid_hydro, "dest_b",   NF90_INT,
     +                 (/dest_dimid_hydro, z_dimid_hydro, 
     +                  b_dimid_hydro, t_dimid_hydro/), 
     +                  dest_b_id_hydro)
      call nf90_check_event(istat, "h:creating variable dest_b")
      istat =  nf90_put_att(ncid_hydro, dest_b_id_hydro, 
     +                  "_FillValue", dest_fill)
      call nf90_check_event(istat, "h:setting att dest_b:_FillValue")
c
      istat =  nf90_def_var(ncid_hydro, "dest_k",   NF90_INT,
     +                 (/dest_dimid_hydro, z_dimid_hydro, 
     +                  b_dimid_hydro, t_dimid_hydro/), 
     +                  dest_k_id_hydro)
      call nf90_check_event(istat, "h:creating variable dest_k")
      istat =  nf90_put_att(ncid_hydro, dest_k_id_hydro, 
     +                  "_FillValue", dest_fill)
      call nf90_check_event(istat, "h:setting att dest_k:_FillValue")
c
      istat =  nf90_enddef(ncid_hydro)  ! leave netCDF define mode
      call nf90_check_event(istat, "h:leaving def mode")

c
c     2) ------------------ temperature data set ------------------
c

      istat = nf90_create(ncfname_temp, NF90_CLOBBER, ncid_temp)   ! overwrite, if exists
      call nf90_check_event(istat, 
     +               "t:creating file:"//trim(adjustl(ncfname_temp)))
      
      if (verbose > 0) 
     +    write(*,336) "temp", trim(adjustl(ncfname_temp))
 
c     --- create netCDF dimensions ---
      istat =  nf90_def_dim(ncid_temp, "t", NF90_UNLIMITED,
     +                      t_dimid_temp)
      call nf90_check_event(istat, "t:creating dimension t")
      istat =  nf90_def_dim(ncid_temp, "b", nboxes, b_dimid_temp)
      call nf90_check_event(istat, "t:creating dimension b")
      istat =  nf90_def_dim(ncid_temp, "z", nlayers+1, z_dimid_temp) ! NB: wet layers + 1 substrate layer
      call nf90_check_event(istat, "t:creating dimension z")
      istat =  nf90_put_att(ncid_temp, NF90_GLOBAL, "title", 
     +                  "temperature from IBMlib")
      call nf90_check_event(istat, "t:setting title attribute")
c     ----- define variables -----   
      istat =  nf90_def_var(ncid_temp, "t", NF90_DOUBLE, 
     +                 (/t_dimid_temp/), 
     +                  t_id_temp)
      call nf90_check_event(istat, "t:creating variable t")
      cbuf = get_datetime(toffset)
      cbuf = "seconds since "//trim(adjustl(cbuf))
      istat =  nf90_put_att(ncid_temp, t_id_temp, "units", 
     +                  trim(adjustl(cbuf)) )
      call nf90_check_event(istat, "t:setting attribute t:units")
      istat =  nf90_put_att(ncid_temp, t_id_temp, 
     +                  "_FillValue", 0.d0)
      call nf90_check_event(istat, "t:setting attribute t:_FillValue")
c
      istat =  nf90_def_var(ncid_temp, "temperature", NF90_DOUBLE,
     +                 (/z_dimid_temp, 
     +                  b_dimid_temp, t_dimid_temp/),  
     +                  temperature_id_temp)
      call nf90_check_event(istat, "t:creating variable temperature")
      istat =  nf90_put_att(ncid_temp, temperature_id_temp, 
     +                  "_FillValue", temp_fill)
      call nf90_check_event(istat, 
     +     "t:setting att temperature:_FillValue")

      istat =  nf90_enddef(ncid_temp)  ! leave netCDF define mode
      call nf90_check_event(istat, "t:leaving def mode")
c
c     3) ------------------ salinity data set ------------------
c
      istat = nf90_create(ncfname_salt, NF90_CLOBBER, ncid_salt)   ! overwrite, if exists
      call nf90_check_event(istat, 
     +               "s:creating file:"//trim(adjustl(ncfname_salt)))
      
      if (verbose > 0) 
     +    write(*,336) "salt", trim(adjustl(ncfname_salt))
 
c     --- create netCDF dimensions ---
      istat =  nf90_def_dim(ncid_salt, "t", NF90_UNLIMITED,
     +                      t_dimid_salt)
      call nf90_check_event(istat, "s:creating dimension t")
      istat =  nf90_def_dim(ncid_salt, "b", nboxes, b_dimid_salt)
      call nf90_check_event(istat, "s:creating dimension b")
      istat =  nf90_def_dim(ncid_salt, "z", nlayers+1, z_dimid_salt) ! NB: wet layers + 1 substrate layer
      call nf90_check_event(istat, "s:creating dimension z")
      istat =  nf90_put_att(ncid_salt, NF90_GLOBAL, "title", 
     +                  "salinity from IBMlib")
      call nf90_check_event(istat, "s:setting title attribute")
c     ----- define variables -----   
      istat =  nf90_def_var(ncid_salt, "t", NF90_DOUBLE, 
     +                 (/t_dimid_salt/), 
     +                  t_id_salt)
      call nf90_check_event(istat, "s:creating variable t")
      cbuf = get_datetime(toffset)
      cbuf = "seconds since "//trim(adjustl(cbuf))
      istat =  nf90_put_att(ncid_salt, t_id_salt, "units", 
     +                  trim(adjustl(cbuf)) )
      call nf90_check_event(istat, "s:setting attribute t:units")
      istat =  nf90_put_att(ncid_salt, t_id_salt, 
     +                  "_FillValue", 0.d0)
      call nf90_check_event(istat, "s:setting attribute t:_FillValue")
c
      istat =  nf90_def_var(ncid_salt, "salinity", NF90_DOUBLE,
     +                 (/z_dimid_salt, 
     +                  b_dimid_salt, t_dimid_salt/),  
     +                  salinity_id_salt)
      call nf90_check_event(istat, "s:creating variable salinity")
      istat =  nf90_put_att(ncid_salt, salinity_id_salt, 
     +                  "_FillValue", salt_fill)
      call nf90_check_event(istat, 
     +     "s:setting att salinity:_FillValue")

      istat =  nf90_enddef(ncid_salt)  ! leave netCDF define mode
      call nf90_check_event(istat, "s:leaving def mode")


      end subroutine initialize_input_files



      subroutine close_input_files()
c     ----------------------------------------------------------------------------
c     Close netCDF files + deallocate tmp buffers
c     ----------------------------------------------------------------------------
      integer :: istat
c     ----------------------------------------------------------------------------
      istat =  nf90_close(ncid_hydro)
      call nf90_check_event(istat, "close_input_files: hydro")

      istat =  nf90_close(ncid_temp)
      call nf90_check_event(istat, "close_input_files: temp")

      istat =  nf90_close(ncid_salt)
      call nf90_check_event(istat, "close_input_files: salt")

      if (allocated(dest_b))    deallocate( dest_b    )
      if (allocated(dest_k))    deallocate( dest_k    )
      if (allocated(h_destmap)) deallocate( h_destmap )
      if (allocated(v_destmap)) deallocate( v_destmap )

      if (verbose > 0) write(*,*) "close_input_files: done"

      end subroutine close_input_files



      subroutine write_input_file_frames(time, vtransp, htransp,
     +                                   temp_avg, salt_avg)
c     -------------------------------------------------------
c     Write a frame to input files
c     -------------------------------------------------------
      integer, intent(in) :: time          ! [sec]
      real, intent(in)    :: htransp(:,0:)   ! accumulated hor transport [m3]          shape = (1:nlaýers, 0:nfaces-1)
      real, intent(in)    :: vtransp(:,0:)   ! accumulated ver transport [m3]          shape = (1:nlaýers, 0:nboxes-1)
      real, intent(in)    :: temp_avg(0:,0:) ! accumulated temperature average [deg C] shape = (0:nlaýers, 0:nboxes-1) = Atlantis layout
      real, intent(in)    :: salt_avg(0:,0:) ! accumulated salinity average [PSU]      shape = (0:nlaýers, 0:nboxes-1) = Atlantis layout
      integer             :: nboxes, nlayers, nfaces
      integer             :: ibox, ifc, ilay, id,il,ib, nt, istat
      real, allocatable   :: exchange(:,:,:)
c     ----------------------------------------------------------------------------
      nboxes  = size(Boxes)
      nlayers = size(zgrid_CC) 
      nfaces  = size(Faces)

c     ---- set exchange buffer ----
      allocate( exchange(dest, nlayers, nboxes) ) 
      exchange = hydro_fill            ! _FillValue
      do ifc = 0, nfaces-1
         do ilay = 1, nlayers     ! scan all, regardless zeros
            id = h_destmap(1, ilay, ifc)
            il = h_destmap(2, ilay, ifc)
            ib = h_destmap(3, ilay, ifc)
            exchange(id, il, ib) =  htransp(ilay, ifc)
          enddo
      enddo
      do ibox = 0, nboxes-1
         do ilay = 1, nlayers     ! scan all, regardless zeros
            id = v_destmap(1, ilay, ibox)
            il = v_destmap(2, ilay, ibox)
            ib = v_destmap(3, ilay, ibox)
            exchange(id, il, ib) =  vtransp(ilay, ibox)
          enddo
      enddo
c
c     ---- append arrays along unlimited dimension ----
c          By default, count(:numDims) = shape(values) and count(numDims + 1:) = 1, where numDims = size(shape(values)). 
c
c     1) hydro data set    
c
      istat = nf90_inquire_dimension(ncid_hydro, t_dimid_hydro, len=nt) ! inquire current number of records
      call nf90_check_event(istat, "h:inquiring t_dimid_hydro")


      istat = nf90_put_var(ncid_hydro, t_id_hydro, 1.d0*time,        ! cast time to double
     +                  start = (/nt+1/))                            ! implicitly count = (/1/)                              
      call nf90_check_event(istat, "h:put_var time")

      istat = nf90_put_var(ncid_hydro, exchange_id_hydro, 
     +                  exchange( :,     :,       :), 
     +                  start = (/1,     1,       1,      nt+1/))       
      call nf90_check_event(istat, "h:put_var exchange")

      istat = nf90_put_var(ncid_hydro, dest_b_id_hydro, 
     +                  dest_b(  1:dest, :,       :), 
     +                  start = (/1,     1,       1,      nt+1/))
     
      call nf90_check_event(istat, "h:put_var dest_b")

      istat = nf90_put_var(ncid_hydro, dest_k_id_hydro,  
     +                  dest_k(  1:dest, :,       :),
     +                  start = (/1,     1,       1,      nt+1/))
      call nf90_check_event(istat, "h:put_var dest_k")
c
c     2) temperature data set    
c
      istat = nf90_inquire_dimension(ncid_temp, t_dimid_temp, len=nt) ! inquire current number of records
      call nf90_check_event(istat, "t:inquiring t_dimid_temp")

      istat = nf90_put_var(ncid_temp, t_id_temp, 1.d0*time,          ! cast time to double
     +                  start = (/nt+1/))                            ! implicitly count = (/1/)
      call nf90_check_event(istat, "t:put_var time")

      istat = nf90_put_var(ncid_temp, temperature_id_temp, 
     +                  temp_avg,                                   
     +                  start = (/1,     1,       1,      nt+1/))       
      call nf90_check_event(istat, "t:put_var temperature")
c
c     3) salinity data set    
c
      istat = nf90_inquire_dimension(ncid_salt, t_dimid_salt, len=nt) ! inquire current number of records
      call nf90_check_event(istat, "s:inquiring t_dimid_salt")

      istat = nf90_put_var(ncid_salt, t_id_salt, 1.d0*time,          ! cast time to double
     +                  start = (/nt+1/))                            ! implicitly count = (/1/)
      call nf90_check_event(istat, "s:put_var time")

      istat = nf90_put_var(ncid_salt, salinity_id_salt, 
     +                  salt_avg,                                   
     +                  start = (/1,     1,       1,      nt+1/))       
      call nf90_check_event(istat, "s:put_var salinity")
c
c     --- clean-up buffers ---
c
      deallocate(exchange)

      end subroutine write_input_file_frames





      subroutine make_input_files(tstart, tend, dt_frames, dt_sampling,
     +                            toffset, hydrofname, tempfname, 
     +                            saltfname)
c     -------------------------------------------------------
c     Generate physics input files (water flow, salinity and temperature) for Atlantis
c     in netCDF format in the period [tstart, tend] with time frame interval dt_frame
c
c     tstart: start of simulation period
c     tend  : end of simulation period
c     dt_frames: interval between time frames dumped in simulation period 
c     dt_sampling: sampling time interval for generating frames; dt_sampling < dt_frames
c     toffset: time offset for time variable in physics input files
c     hydrofname: file name for hydrodynamic data
c     tempfname:  file name for temperature data
c     saltfname:  file name for salinity data
c     
c     Round number of sampling steps truncated, if dt_frames/dt_sampling is not an integer
c     Apply simple Eulerian flux integration to evaluate net transport
c     -------------------------------------------------------
      type(clock), intent(in)   :: tstart        ! start time of time frame generation
      type(clock), intent(in)   :: tend          ! end time of time frame generation
      integer, intent(in)       :: dt_frames     ! [sec] interval between time frames dumped in simulation period 
      integer, intent(in)       :: dt_sampling   ! [sec] sampling time interval for generating frames; dt_sampling < dt_frames    
      type(clock), intent(in)   :: toffset       ! time offset - used to reference time in physics input as seconds since offset
      character*(*), intent(in) :: hydrofname    ! file name for hydrodynamic data
      character*(*), intent(in) :: tempfname     ! file name for temperature data
      character*(*), intent(in) :: saltfname     ! file name for salinity data
c
      real, allocatable          :: hflux(:,:)    ! horizontal water flux     [m3/sec] shape = (1:nlaýers, 0:nfaces-1) 
      real, allocatable          :: vflux(:,:)    ! vertical water flux       [m3/sec] shape = (1:nlaýers, 0:nboxes-1) 
      real, allocatable          :: htransp(:,:)  ! accumulated hor transport [m3]     shape = (1:nlaýers, 0:nfaces-1)
      real, allocatable          :: vtransp(:,:)  ! accumulated ver transport [m3]     shape = (1:nlaýers, 0:nboxes-1) 
      real, allocatable          :: temp_now(:,:)   ! current temperature [deg C] 
      real, allocatable          :: salt_now(:,:)   ! current salinity [deg C] 
      real, allocatable          :: temp_avg(:,:)   ! accumulated temperature average [deg C]
      real, allocatable          :: salt_avg(:,:)   ! accumulated salinity average [PSU]
      integer                    :: nfaces, nboxes, nlayers, nframes
      integer                    :: ntinner, iframe, it, t, tdiff, istat
      type(clock), pointer       :: tnow
c     -------------------------------------------------------
      nlayers = size(zgrid_CC)
      nfaces  = size(Faces)   ! range = 0:nfaces-1
      nboxes  = size(Boxes)   ! range = 0:nboxes-1

      allocate( hflux(  nlayers, 0:nfaces-1) ) 
      allocate( vflux(  nlayers, 0:nboxes-1) )
      allocate( htransp(nlayers, 0:nfaces-1) ) 
      allocate( vtransp(nlayers, 0:nboxes-1) )

      allocate( temp_now(0:nlayers, 0:nboxes-1) )
      allocate( salt_now(0:nlayers, 0:nboxes-1) )
      allocate( temp_avg(0:nlayers, 0:nboxes-1) ) 
      allocate( salt_avg(0:nlayers, 0:nboxes-1) ) 

      ntinner = dt_frames/dt_sampling  ! inner time loop, truncate remainder, if any
      if (mod(dt_frames, dt_sampling) /= 0) then 
         write(*,*) "make_input_files: mod(dt_frames/dt_sampling) != 0" 
         write(*,*) "dt_frames   = ", dt_frames
         write(*,*) "dt_sampling = ", dt_sampling
         stop 491
      endif
      call get_period_length_sec(tstart, tend, tdiff)  !
      nframes = tdiff/dt_frames        ! outer time loop, truncate remainder, if any
      call get_period_length_sec(toffset, tstart, t)  ! t=tstart referenced to toffset
      tnow    => get_master_clock()
      call set_clock(tnow, tstart)
      if (verbose > 0) then
         write(*, 832) "start", get_datetime(tstart)
         write(*, 832) "end  ", get_datetime(tend)
         write(*, 833) "period between frames   ", dt_frames
         write(*, 833) "requested frame sampling",  dt_sampling
         write(*, 834) ntinner
      endif
      call initialize_input_files(hydrofname, tempfname, saltfname,
     +                            toffset)
      do iframe = 1, nframes     ! outer loop begin
         if (verbose > 0) write(*,831) iframe, nframes
         htransp = 0.0           ! reset accumulants
         vtransp = 0.0           ! reset accumulants
         temp_avg = 0.0          ! reset accumulants
         salt_avg = 0.0          ! reset accumulants
         do it = 1, ntinner      ! inner loop begin
            call update_physical_fields()
            call get_horizontal_water_fluxes(hflux)
            call get_vertical_water_fluxes(vflux)
            htransp = htransp + hflux*dt_sampling ! simple Eulerian flux integration
            vtransp = vtransp + vflux*dt_sampling ! simple Eulerian flux integration
            call get_box_phys_averages(temp_now, salt_now) 
            temp_avg = temp_avg + temp_now  
            salt_avg = salt_avg + salt_now
            call add_seconds_to_clock(tnow, dt_sampling) 
            t = t + dt_sampling   ! manually sync t and tnow ..
         enddo  ! it (inner time loop)
         temp_avg = temp_avg/ntinner  ! plain average over inner time steps
         salt_avg = salt_avg/ntinner  ! plain average over inner time steps
         call write_input_file_frames(t, vtransp, htransp, 
     +                                temp_avg, salt_avg) 
      enddo     ! iframe (outer time loop)
c
c     --- wind up ---
c
      call close_input_files()
      deallocate( hflux )
      deallocate( vflux )
      deallocate( htransp )     
      deallocate( vtransp )
      deallocate( temp_now )
      deallocate( salt_now )
      deallocate( temp_avg )
      deallocate( salt_avg )
c     
 831  format("make_input_files: writing frame",i6,"/",i6)
 832  format("make_input_files:",1x,a,1x,"time =",1x,a) 
 833  format("make_input_files:",1x,a,1x,i8,1x,"sec")    
 834  format("make_input_files: inner time steps =",i6)      
c      
      end subroutine make_input_files



c     ==========================   destructors   ==============================================    
                     
      
      subroutine deallocate_AtlantisBox(abox)
c     -------------------------------------------------------------------- 
      type(AtlantisBox), intent(inout) :: abox
c     -------------------------------------------------------------------- 
      call delete_lonlat_polygon(abox%shape)
      if (associated( abox%faces ))     deallocate( abox%faces ) 
      if (associated( abox%neighbors )) deallocate( abox%neighbors ) 
      if (associated( abox%xy ))        deallocate( abox%xy ) 
      if (associated( abox%da ))        deallocate( abox%da ) 
      end subroutine deallocate_AtlantisBox


      subroutine deallocate_AtlantisFace(aface)
c     --------------------------------------------------------------------
      type(AtlantisFace), intent(inout) :: aface
c     -------------------------------------------------------------------- 
      if (associated( aface%xy ))        deallocate( aface%xy ) 
      if (associated( aface%dl ))        deallocate( aface%dl ) 
      end subroutine deallocate_AtlantisFace



      subroutine close_atlantis_grid()
c     --------------------------------------------------------------------
c     Deallocate module data
c     --------------------------------------------------------------------
      integer :: i
c     --------------------------------------------------------------------
      if (allocated( zgrid_bounds )) deallocate( zgrid_bounds )
      if (allocated( zgrid_CC ))     deallocate( zgrid_CC )
      if (allocated( dz ))           deallocate( dz )         
      if (allocated( nzpts ))        deallocate( nzpts )         
      if (allocated( zoffset ))      deallocate( zoffset )
      
      if (allocated( Boxes )) then
         do i = 0, size(Boxes)-1   ! offset = 0
            call deallocate_AtlantisBox(Boxes(i))
         enddo
         deallocate( Boxes )
      endif
      
      if (allocated( Faces )) then
         do i = 0, size(Faces)-1   ! offset = 0
            call deallocate_AtlantisFace(Faces(i))
         enddo
         deallocate( Faces )
      endif

      end subroutine close_atlantis_grid    

  
      end module
