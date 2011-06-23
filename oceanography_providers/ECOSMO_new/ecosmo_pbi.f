ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     ECOSMO pbi 
c     ---------------------------------------------------
c     $Rev: 310 $
c     $LastChangedDate: 2011-02-22 15:16:23 +0100 (Tue, 22 Feb 2011) $
c     $LastChangedBy: asch $ 
c
c     PBI for ECOSMO data output (in 2011 binary format)
c
c     NOTE: this file must be compiled with UNFORMATTED_LITTLE_ENDIAN
c           specify how your compiler sets this in compiler_defaults.mk
c
c     TODO: fix ibio=6 load
c           non CC interpolate_turbulence       (remember to uncomment mesh_grid load)
c           non CC interpolate_turbulence_deriv (remember to uncomment mesh_grid load)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module physical_fields

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
      
      character*999             :: hydroDBpath ! hydrographic data sets
      character*3, parameter    :: prefix = "tl5" ! later move to input file
      integer, parameter        :: tag_lenght = 7 ! SSSYYMM like tl58001
      
      character(len=tag_lenght) :: cur_tag      ! tag SSSYYMM of current set in buffer
      integer                   :: cur_frame    ! current time frame in buffer
      integer                   :: phys_file = -1 ! physics file handler, assigned below
      integer                   :: bio_file  = -1 ! NPZD file handler, assigned below
      character*(*), parameter  :: not_set_c = "not set"
      integer, parameter        :: not_set_i = -9999  

c     --- 3D grids ---
                
      real,allocatable,target :: ccdepth0(:,:,:)   ! reference cell center depth water below surface [m] (dslm=0)
      real,allocatable,target :: acc_width0(:,:,:) ! accumulated water above this undisturbed layer [m] dim=nz+1   
      
c     --- 2D grids ---
      
      real,allocatable :: wdepth0(:,:)      ! reference depth at cell-center, positive  [m] - dry points negative
 

c     --- 1D grids ---

      real, allocatable :: levels(:)        ! center depth of layers for DSLM=0. [m] dim=nz 
      real, allocatable :: layer_width(:)   ! width for undisturbed layers [m] dim=nz 
      real, allocatable :: layer_faces(:)   ! accumulated width for undisturbed layers above this [m] dim=nz+1 
      

c     ------ ECOSMO auxillary dimensions and grids
      integer               :: khor,ndrei,kasor
      integer               :: nbio ! NPZD boxes
      integer, allocatable  :: iwet(:)  
      integer, allocatable  :: lazc(:)  
      integer, allocatable  :: indend(:) 
      integer, allocatable  :: iland(:,:,:)
      integer, allocatable  :: iizwei(:,:)
      integer, allocatable  :: iidrei(:,:,:)

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
      call reset_frame_handler()

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
      get_pbi_version =  "SUNFISH pbi version: $Rev: 310 $"
      end function

   


      subroutine read_grid_desc()
c     ---------------------------------------------------
c     Read grid definition which is contained in file 
c     given as tag grid_desc in the simulation file simulation_file
c     Read topography which is contained in file 
c     given as tag topography
c
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
      character*999        :: gd_fname, topo_fname
      type(control_file)   :: grid_ctrlfile
      integer              :: topo_file
      integer              :: ix,iy,iz,ihit,idx,nwet,nacc,nbott
      integer              :: iostat,iilo,nhoru,ntotu,izcheck
      integer              :: lw,lwa,lwe,i,j,k,lump
      real                 :: wd0
      real                 :: lambda1, dlambda, phi1, dphi ! LOCAL DUMMIES
      real,allocatable     :: dz(:)
      integer,allocatable  :: ldep(:)
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
    
c.....read internal ECOSMO dimensions      
      call read_control_data(grid_ctrlfile,"khor",  khor)  ! module dimension
      call read_control_data(grid_ctrlfile,"ndrei", ndrei) ! module dimension    
      call read_control_data(grid_ctrlfile,"kasor", kasor) ! module dimension
      call read_control_data(grid_ctrlfile,"nbio",  nbio)  ! module dimension
      write(*,234) "khor",  khor
      write(*,234) "ndrei", ndrei
      write(*,234) "kasor", kasor
      write(*,234) "nbio",  nbio

 229  format(a,a)    
 231  format("read_grid_desc: 3d grid dim (nx,ny,nz) = ", i4,i4,i4)
 232  format("read_grid_desc: lambda1 = ",f12.7,
     +                      " dlambda = ",f12.7," degrees")  
 233  format("read_grid_desc: phi1    = ",f12.7,
     +                      " dphi    = ",f12.7," degrees")  
 234  format("read_grid_desc: ECOSMO dimension ",a," = "i8)


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
     
c     --- 1D grids ---
      
      allocate( levels(nz)        ) 
      allocate( layer_width(nz)   )  ! width for undisturbed layers in meters
      allocate( layer_faces(nz+1) )  ! meters of undisturbed water above layer
      
      write(*,*) "read_grid_desc: allocate grid arrays: OK"   

      call close_control_file(grid_ctrlfile)

c
c.....read topography -> wdepth0 (only wety points provided, other assumed dry)
c
c     1) set bottom_layer
c     2) adjust ccdepth0, so that actual width of bottom layer is consistent 
c        with wdepth0 by stretching the last wet layer to the topographic bottom
c         
c
      call read_control_data(simulation_file, "topography", topo_fname)
      write(*,229) "read_grid_desc: reading topography from: ", 
     +           trim(adjustl(topo_fname))
      
      call find_free_IO_unit(topo_file)

c     -------- parsing code from ECOSMO: begin --------
c     dz= lower layer boundary vector [m]
c
c     common /ind/      iwet(khor1),ldep(khor),lazc(khor),indend(n)
c     common /field3p1/ dz(ilo)
      open(topo_file,file=topo_fname,status='old',iostat=iostat)
      if (iostat /= 0) then
         write(*,*) "Error opening topography file "//
     *               trim(adjustl(topo_fname))
         stop
      endif   
      allocate( dz(nz)          )  ! local array, needed for topography at start
      allocate( ldep(khor)      )  ! local array, needed for topography at start
      allocate( iwet(khor + 1)  )  ! module array, needed for time frames
      allocate( lazc(khor)      )  ! module array, needed for time frames
      allocate( indend(nx)      )  ! module array, needed for time frames
      allocate( iland(ny,nx,nz) )  ! module array, needed for time frames, NB: permuted indices
      allocate( iizwei(ny,nx)   )  ! module array, maybe needed later, NB: permuted indices
      allocate( iidrei(ny,nx,nz))  ! module array, maybe needed later, NB: permuted indices

      read (topo_file,'(i8)')   iilo
      read (topo_file,'(10f10.2)')(dz(j),j=1,nz)  ! -> levels
      read (topo_file,'(2i8)')  nhoru,ntotu ! not used, implicit integers
      read (topo_file,'(10i10)')iwet 
      read (topo_file,'(10i10)')ldep 
      read (topo_file,'(10i10)')lazc 
      read (topo_file,'(10i10)')indend 
c     ... there is more in ECOSMO topography file, but this
c         is not used so skip it      
      close (topo_file)


c---------------------------------------------------------------------
c     computation of help arrays on native ECOSMO grid:
c
c             ECOSMO dim  NS setup
c       i=1,ny    m         177     north to south scan
c       j=1,nx    n         207     west to east scan
c       k=1,nz   ilo        20      surface to bottum scan
c
c
c     iizwei(i,j)  = 2d vector with grid index of compressed 2-d fields 
c     iidrei(i,j,k)= 3d vector with grid cell index for compressed 3-d fields
c     iland(i,j,k) = land/water distribution dry grid cells are marked by 0
c     ibotlay(i,j) = number of max. vertical layers for each grid point
c                    ibotlay=0 for all dry grid columns
c     topo(i,j)    = topography array
c
c
c     Index transformation ECOSMO to IBMlib regular lonlat grid:
c
c        ECOSMO(i,j,k) = lonlat_data(j, ny+1-i, k)
c        ix = j
c        iy = ny+1-i
c        iz = k 
c---------------------------------------------------------------------
      
c     ... initialize wdepth0 + bottom_layer 
c         and ECOSMO aux arrays iizwei,iidrei,iland

      wdepth0      = -1.0  ! default dry
      bottom_layer =  0    ! default dry
      iizwei  =  0   
      iidrei  =  0   
      iland   =  0   

      lwe=0
      nwet=0
      do j=1,nx  ! replaced n-> nx
         lwa=lwe+1
         lwe=indend(j)
         do lw=lwa,lwe
            lump=lazc(lw)
            i=iwet(lw)
            iizwei(i,j)=lw
c            ibotlay(i,j)=lump                      ! ECOSMO assignment
c            topo(i,j)=dz(lump-1) + real(ldep(lw))  ! ECOSMO assignment
            
            ix = j               ! corresponding regular lonlat indices
            iy = ny+1-i          ! corresponding regular lonlat indices
            bottom_layer(ix, iy) = lump
            wdepth0(ix, iy)      = dz(lump-1) + real(ldep(lw))

            do k=1,lump
               nwet=nwet+1
               iidrei(i,j,k)=nwet
               iland(i,j,k)=1  
            enddo
         enddo
      enddo
                              
c     -------- parsing code from ECOSMO: end --------        
c
c.....set vertical layer structure: levels -> layer_width0 + acc_width0
c
      layer_faces(1)      = 0.
      layer_faces(2:nz+1) = dz
      layer_width         = layer_faces(2:nz+1) - layer_faces(1:nz)
      levels              = layer_faces(1:nz) + 0.5*layer_width
     
c     set acc_width0, layer_width0 from levels
c     last element acc_width0(nz+1) is total (dslm=0) max depth of grid
      
      acc_width0(:,:,1) = 0. ! sea surface (DSLM=0)
      do iz=1,nz 
         ccdepth0(:,:,iz)     = levels(iz)
         acc_width0(:,:,iz+1) = layer_faces(iz+1)
      enddo
      

      write(*,241) levels
      write(*,242) layer_width
      write(*,243) layer_faces

 241  format("read_grid_desc: levels = ",1000f6.1)     
 242  format("read_grid_desc: calc ref layer_width = ",1000f6.1)     
 243  format("read_grid_desc: calc ref layer_faces = ",1000f6.1)

c
c     readjust vertical width of grid cells encompassing the sea bed 
c     so that the topography is realistic based on arrays 
c     wdepth0 and bottom_layer read from ECOSMO topography file
c     assume read data is consistent with level structure 
c     readjust the arrays ccdepth0, acc_width0
c     bottom_layer(ix,iy) = last wet layer (0 for dry points) nx,ny
c
      do ix=1,nx
         do iy=1,ny
            iz = bottom_layer(ix,iy)   
            if (iz==0) cycle ! skip if this is a dry point
c...........check input consistency:                       
c           layer_faces(izcheck) < wdepth0(ix,iy) < layer_faces(izcheck+1)
c           capture boundary case wdepth0(ix,iy) ~ layer_faces(izcheck+1)
c           rely on left-to-right evaluation of if clause
            call search_sorted_list(wdepth0(ix,iy),layer_faces,izcheck) 
            if ((iz /= izcheck).and.
     +          (abs(wdepth0(ix,iy)-acc_width0(ix,iy,iz+1))>1.e-3)) then 
               write(*,*)"read_grid_desc: inconsistent input "
               write(*,*)"read_grid_desc: ix,iy=",ix,iy
               write(*,*)"read_grid_desc: iz,izcheck=",iz,izcheck
               write(*,*)"read_grid_desc: wdepth0=",wdepth0(ix,iy)
               stop
            endif 
c...........re-center mid point of bottom cell
            acc_width0(ix, iy, iz+1) = wdepth0(ix, iy)
            ccdepth0(ix, iy, iz)     = 0.5*(acc_width0(ix, iy,iz) + 
     +                                         acc_width0(ix, iy,iz+1)) 
c...........correct cell below (if any) for completeness             
            if (iz<nz) then 
               ccdepth0(ix, iy, iz+1) = 0.5*(acc_width0(ix, iy,iz+1) + 
     +                                       acc_width0(ix, iy,iz+2)) 
            endif
         enddo
      enddo
     
      write(*,*) "read_grid_desc: "//
     +           "readjusted bottom cells reflecting topogrphy"

      acc_width = acc_width0 ! allow certain grid functions to work on land
      ccdepth   = ccdepth0   ! allow certain grid functions to work on land
c
c     --- define the wetmask statically (no flooding/drying) --- 
c   
      where(bottom_layer > 0)
         wetmask = 1 ! wet
      elsewhere
         wetmask = 0 ! dry
      end where 

      write(*,*) "read_grid_desc: grid arrays initialized"
c
c     --- final clean up --- 
c   
      deallocate( dz   )
      deallocate( ldep )

c     ------------------------------------------       
      end subroutine read_grid_desc



      subroutine close_physical_fields()
c     ------------------------------------------ 
c     ------------------------------------------             
      if (allocated(ccdepth0))    deallocate( ccdepth0 )
      if (allocated(acc_width0))  deallocate( acc_width0 )
      if (allocated(acc_width))   deallocate( acc_width )
      if (allocated(wdepth0))     deallocate( wdepth0 )
 
      if (allocated(levels))      deallocate( levels )
      if (allocated(layer_width)) deallocate( layer_width )
      if (allocated(layer_faces)) deallocate( layer_faces )
 
c     ..... ECOSMO auxillary arrays .....
      if (allocated(iwet))        deallocate( iwet   )
      if (allocated(lazc))        deallocate( lazc   )
      if (allocated(indend))      deallocate( indend )
      if (allocated(iland))       deallocate( iland  )
      if (allocated(iizwei))      deallocate( iizwei )
      if (allocated(iidrei))      deallocate( iidrei )
    
      call close_mesh_grid()
      call reset_frame_handler()
      call close_horiz_grid_transf()
c     ------------------------------------------
      end subroutine 


      subroutine update_physical_fields(time, dt)
c     ------------------------------------------  
      type(clock), intent(in),optional :: time
      integer, intent(in),optional     :: dt
      logical                          :: update
      character(len=tag_lenght)        :: tag
      integer                          :: frame
      type(clock), pointer             :: aclock
c     ------------------------------------------
      aclock => get_master_clock()
      if (present(time)) then
         call set_master_clock(time)
      elseif (present(dt)) then
         call add_seconds_to_clock(aclock, dt)
      endif     
c
      call resolve_corresp_dataset(aclock, tag, frame)
      call update_dataset(tag, frame)
c     ------------------------------------------       
      end subroutine update_physical_fields




      subroutine resolve_corresp_dataset(aclock, tag, frame)
c     -------------------------------------------------------
c     Resolve tag and and data frame in file corresponding to aclock
c
c     In the ECOSMO/MEECE set, data are bundled by months with file names
c       SSSYYMM   (physics)
c       SSSYYMMb  (biogeochemistry)
c     where YY is the last two year digits and MM is the month-in-year
c     and SSS is a scenarie prefix. The tag corresponding to these
c     fiels is SSSYYMM. There is a Y2K problem here ...
c     E.g. physics for Jan 1980 for scenarie tl5 is in file: tl58001
c     corresponding to the tag "tl58001"
c     -------------------------------------------------------
      type(clock), intent(in)                :: aclock
      character(len=tag_lenght), intent(out) :: tag
      integer, intent(out)                   :: frame
      integer          :: year, month, day
      real             :: frame_r
c     -------------------------------------------------------
      call get_date_from_clock(aclock, year, month, day)
      write(tag,455) prefix, mod(year,100), month  
      frame = day   
 455  format(a3,i2.2,i2.2)
      end subroutine resolve_corresp_dataset




      subroutine update_dataset(tag, frame)
c     ------------------------------------------
c     Only open new data set file if necessary
c     Only load new frames if necessary
c     Only this subroutine updates buffer state descriptors
c     ------------------------------------------
      character(len=tag_lenght), intent(inout) :: tag
      integer, intent(inout)                   :: frame
      logical  :: needed
c     ------------------------------------------
      needed = ((cur_tag==not_set_c).or.(cur_frame==not_set_i))
c
      if (needed.or.(tag /= cur_tag)) then
         call reset_frame_handler() ! resets cur_frame+close old
         call open_data_files(tag)
         cur_tag = tag
         needed = .true. ! force update if accidentally frame == cur_frame
      endif
c      
      if (needed.or.(frame /= cur_frame)) then
         call load_data_frames(frame)
         cur_frame = frame
      endif  
c 
      end subroutine update_dataset





      subroutine reset_frame_handler()
c     ------------------------------------------
      call close_data_files()
      cur_tag         = not_set_c
      cur_frame       = not_set_i   
      end subroutine reset_frame_handler



      subroutine close_data_files()
c     ------------------------------------------
c     Only close, if data files are open
c     At close release unit to avoid another write opens it
c     ------------------------------------------
      logical :: isopen
c     ------------------------------------------
      inquire(unit=phys_file,opened=isopen)
      if (isopen) then
         close( phys_file )
         phys_file = not_set_i  ! release unit 
         write(*,*) "close_data_files: closed phys_file"
      endif

      inquire(unit=bio_file,opened=isopen)
      if (isopen) then
         close( bio_file )
         bio_file = not_set_i  ! release unit 
         write(*,*) "close_data_files: closed bio_file"
      endif
         
      end subroutine close_data_files


      subroutine open_data_files(tag)
c     ------------------------------------------
c     Open binary file set corresponding to 
c     tag == SSSYYMM 
c     phys_file -> file hydroDBpath/SSSYYMM
c     bio_file  -> file hydroDBpath/SSSYYMMb
c     ------------------------------------------
      character(len=tag_lenght),intent(in) :: tag
      character*999                        :: fname
      integer :: iostat
c     ------------------------------------------
      call close_data_files() ! in case they are open ...
c      
      call find_free_IO_unit(phys_file)
      write(fname,333) trim(adjustl(hydroDBpath)), tag
      open(phys_file, file=fname, 
     +                status="old", form="unformatted", 
     +                action="read", iostat=iostat)
      if (iostat /= 0) then
         write(*,*) "error ",iostat,"opening phys_file: ", 
     +               trim(adjustl(fname))
         stop
      endif
c
      call find_free_IO_unit(bio_file)
      write(fname,334) trim(adjustl(hydroDBpath)), tag
      open(bio_file, file=fname, 
     +               status="old", form="unformatted", 
     +               action="read", iostat=iostat)  
      if (iostat /= 0) then
         write(*,*) "error ",iostat,"opening bio_file: ", 
     +               trim(adjustl(fname))
         stop
      endif
       
 333  format(a,"/",a)
 334  format(a,"/",a,"b")

      write(*,*) "open_data_files: opened files for tag ", tag

      end subroutine open_data_files

 
      
      subroutine load_data_frames(frame)
c     ------------------------------------------
c     Verbose reading of time frame = frame
c
c     The current implementation is a little inefficient 
c     for backtracking - because data files are written
c     with access == "sequential", it is not possible robustly
c     to read them with access == "direct", picking 
c     a specific record. Therefore rewind + skip must be performed
c     Options to be tested later: 
c        1) try ifort option access == "stream" which
c           allows file positioning. This is not standard F90
c        2) Cache all time frames in a file (approx 50-260 Mb,
c           depending on how much is used). Only uncompress
c           frame(s) that are needed
c     Reading is somewhat fragile - do not try to coalesence 
c     (or split) read lines
c     
c     ibio=3 - flaggelates  [mgC/m**3]
c     ibio=4 - diatomes     [mgC/m**3]
c     ibio=5 - small zoo    [mgC/m**3]
c     ibio=6 - large zoo    [mgC/m**3]
c     ibio=7 - Detritus     [mgC/m**3]
c     ibio=8 - NH4          [mmolN/m**3]
c     ibio=9 - NO2          [mmolN/m**3]
c     ibio=10- NO3          [mmolN/m**3]
c     ibio=11- PO4          [mmolP/m**3]
c     ibio=12- SiO2         [mmolSi/m**3]
c     ibio=13- O2           [mmolO/m**3]
c     ibio=15- cyanopbacteria [mgC/m**3]
c     ibio=16- sediment 1 detritus [mgC/m**3]
c     ibio=17- sediment 2 silicate [mgC/m**3]

c     ------------------------------------------
      integer, intent(in) :: frame ! pick frame
      integer :: ir, nread, i1, ix,iy
      integer :: lw,lwa,lwe,i,j,k,lump,ibio
      real    :: dz,h
      integer :: mjarcc,lmoncc,ndaycc,iviercc,ipcc
c     ---- automatic arrays as reading buffers
      real :: zmit(khor)      ! sea surface elevation 
      real :: umit(ndrei)     ! u current component
      real :: vmit(ndrei)     ! v current component
      real :: wcmit(ndrei)    ! w current component
      real :: acmit(ndrei)    ! viscosity
      real :: szmit(ndrei)    ! schmid number
      real :: tcmit(ndrei)    ! temperarature
      real :: scmit(ndrei)    ! salinity
      real :: dump_nynx(ny,nx)  
      real :: dump_khor(khor)
      real :: Tc(ndrei,3:nbio) ! NPZD

c     ------------------------------------------ 
c
c     ======= figure out where file pointer is now =======
c
      if ((cur_frame == not_set_i).or.(frame < cur_frame)) then
         rewind(phys_file)
         nread = frame
      else
         nread = frame - cur_frame
      endif
c
c     ======= advance to the needed frame =======
c
      do ir=1,nread
         read(phys_file) mjarcc,lmoncc,ndaycc,iviercc,ipcc
         read(phys_file) zmit   
         read(phys_file) umit
         read(phys_file) vmit
         read(phys_file) wcmit
         read(phys_file) acmit
         read(phys_file) szmit
         read(phys_file) tcmit
         read(phys_file) scmit
         read(phys_file) dump_nynx  ! frimit(m,n)
         read(phys_file) dump_nynx  ! hismit(m,n)
         read(phys_file) dump_nynx  ! hisrmit(m,n)
         read(phys_file) dump_nynx  ! tismit(m,n)
         read(phys_file) dump_nynx  ! uimit(m,n)
         read(phys_file) dump_nynx  ! vimit(m,n)
         read(phys_file) dump_nynx  ! fqgmit(m,n)
         read(phys_file) dump_khor  ! fqrmit(khor
         read(phys_file) dump_khor  ! fqsmit(khor)
         read(phys_file) dump_khor  ! fqlmit(khor)
         read(phys_file) dump_nynx  ! qois(m,n)
         read(phys_file) dump_nynx  ! qiis(m,n)
         read(bio_file)  Tc
      enddo
c
c     ======= uncompress and parse last frame to 3D arrays =======
c      
c     Index transformation ECOSMO to IBMlib regular lonlat grid:
c
c        ECOSMO(i,j,k) = lonlat_data(j, ny+1-i, k)
c
c
c     ---- update  dslm(:,:) current sea surface elevation over reference [m] <- zmit   
c
c     zmit: assumed positive up 
      
      
      lwe=0
      do j=1,nx   ! replaced n-> nx
         lwa=lwe+1
         lwe=indend(j)
         do lw=lwa,lwe
            lump=lazc(lw)
            i=iwet(lw)  
c            ddz(i,j,1)=ddz1(i,j,1)+zac(lw) ! ECOSMO assignment
            ix = j               ! corresponding regular lonlat indices
            iy = ny+1-i          ! corresponding regular lonlat indices
            dslm(ix, iy) = zmit(lw)
         enddo
      enddo

c     ------ add DSLM to auxillary grid descriptors for wet points  
c            acc_width(ix,iy,1) = 0 (sea surface) 
c            update ccdepth
c            update acc_width
c            update wdepth
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
    
c     ---- update u(:,:,:) current [m/s] (positive east)  <- umit
c                 h = layer height at east side of cell
c                 u is on east side of cell

      u  = 0. ! ignore padval_u
      i1 = 1
      do j=1,nx ! replaced n-> nx
        do i=1,ny
           do k=1,nz 
              if (iland(i,j,k)>0) then                
c                 ug(i,j,k)=uc(i1)/((ddz(i,j,k)+ddz(i,j-1,k))/2.) ! ECOSMO assignment
                 ix = j         ! corresponding regular lonlat indices
                 iy = ny+1-i    ! corresponding regular lonlat indices
                 h  = acc_width(ix,iy,k+1) - acc_width(ix,iy,k)
                 h  = h + acc_width(ix+1,iy,k+1) - acc_width(ix+1,iy,k)
                 h  = h/2.  
                 u(ix,iy,k) = umit(i1)/h 
                 i1 = i1+1
              endif
           enddo
        enddo
      enddo

c     ---- update v(:,:,:)  current [m/s] (positive north)
c                 h = layer height at south side of cell
c                 v is on south side of cell

      v  = 0. ! ignore padval_v 
      i1 = 1
      do j=1,nx ! replaced n-> nx
        do i=1,ny
           do k=1,nz 
              if (iland(i,j,k)>0) then                
c                vg(i,j,k)=vc(b)/((ddz(i,j,k)+ddz(i-1,j,k))/2.) ! ECOSMO assignment
                 ix = j         ! corresponding regular lonlat indices
                 iy = ny+1-i    ! corresponding regular lonlat indices
                 h  = acc_width(ix,iy,k+1) - acc_width(ix,iy,k)
                 h  = h + acc_width(ix,iy-1,k+1) - acc_width(ix,iy-1,k)
                 h  = h/2.  
                 v(ix,iy,k) = vmit(i1)/h 
                 i1 = i1+1
              endif
           enddo
        enddo
      enddo

       
c     ---- update w(:,:,:) current [m/s] (positive down) <- wmit
c                 w is on upper cell face (horizontal mid cell)
c                 wmit is positive toward surface (opposite z orientation)
c                 flip sign of w: physical_fields sign convention is positive down, 
c                 ECOSMO convention is positive up

      w  = 0.  ! ignore padval_w - also imposes w(bottum) = 0
      i1 = 1
      do j=1,nx ! replaced n-> nx
        do i=1,ny
           do k=1,nz 
              if (iland(i,j,k)>0) then                
c                wg(i,j,k)=wc(c)/1000.   ! ECOSMO assignment
                 ix = j         ! corresponding regular lonlat indices
                 iy = ny+1-i    ! corresponding regular lonlat indices
                 w(ix,iy,k) = wcmit(i1)
                 i1 = i1+1
              endif
           enddo
        enddo
      enddo
     


      w = -w ! flip sign, physical_fields sign convention is positive down

c     ---- update  CC scalars:
c     ----     temp(:,:,:) Water temp. [Celcius]
c     ----     salinity(:,:,:)  Salinity. [psu]

c     temp/salinity init value: rely on mesh_grid padval_
      i1 = 1
      do j=1,nx ! replaced n-> nx
        do i=1,ny
           do k=1,nz 
              if (iland(i,j,k)>0) then                
                 ix = j         ! corresponding regular lonlat indices
                 iy = ny+1-i    ! corresponding regular lonlat indices
                 temp(ix,iy,k)     = tcmit(i1)
                 salinity(ix,iy,k) = scmit(i1)
                 i1 = i1+1
              endif
           enddo
        enddo
      enddo

c     ---- update vdiffus(:,:,:) vertical   diffusivity [m**2/s]     
c          vdiffus is on upper cell face (horizontal mid cell)
c          diffusion coefficient = schmid number (szmit)* viscosity (acmit) 
c           
      vdiffus = 0. ! ignore padval_vdiffus
      i1 = 1
      do j=1,nx ! replaced n-> nx
        do i=1,ny
           do k=1,nz 
              if (iland(i,j,k)>0) then                
                 ix = j         ! corresponding regular lonlat indices
                 iy = ny+1-i    ! corresponding regular lonlat indices
                 vdiffus(ix,iy,k) = scmit(i1)*acmit(i1)
                 i1 = i1+1
              endif
           enddo
        enddo
      enddo

      vdiffus(:,:,1) = vdiffus(:,:,2)  ! apply surface BC
                
      
c     ---- update hdiffus(:,:,:)    ! horizontal diffusivity [m**2/s]
c     keep init value == molecular_diffusivity 
c
c     ---- update zoo(:,:,:)        ! Zooplankton [10^-3 mol N/liter]
      zoo  = 0.0
      ibio = 6
      i1   = 1
      do j=1,nx ! replaced n-> nx
        do i=1,ny
           do k=1,nz 
              if (iland(i,j,k)>0) then                
                 ix = j         ! corresponding regular lonlat indices
                 iy = ny+1-i    ! corresponding regular lonlat indices
                 zoo(ix,iy,k) = Tc(i1,ibio)
                 i1 = i1+1
              endif
           enddo
        enddo
      enddo
      write(*,*) "loaded zoo group ", ibio, "to zoo"
 
      write(*,233) frame, mjarcc,lmoncc,ndaycc
 233  format("load_data_frames: loaded frame ",i3," (labels=",3i3,")")
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

      
