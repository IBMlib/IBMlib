ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     POLCOMS+ERSEM pbi for online coupling
c     ---------------------------------------------------
c     $Rev: 228 $
c     $LastChangedDate: 2011-01-26 02:52:59 +0100 (Wed, 26 Jan 2011) $
c     $LastChangedBy: asch $ 
c    
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module physical_fields

      use mesh_grid
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
      public :: get_master_clock ! from time_services
      public :: set_master_clock ! from time_services
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

c     --- online extensions ---
      public :: physical_data_bucket
      public :: transfer_data_to_ibmlib



c     ---- Exchange data type between IBMlib and POLCOMS+ERSEM
c          keep internal structure public
c          3D data pointers is expected to point to
c          buffers with the POLCOMS layout: (zdim,londim,latdim)
c          No constructor is provided - it is up to the external hydrodynamic model 
c          either allocate/point pointers in physical_data_bucket to data
c          buffers in a physical_data_bucket are only read (not modified)
c
      type physical_data_bucket
        real(kind=8),pointer :: ucur(:,:,:)   ! POLCOMS layout: (zdim,londim,latdim)
        real(kind=8),pointer :: vcur(:,:,:)   ! POLCOMS layout: (zdim,londim,latdim)
        real(kind=8),pointer :: wcur(:,:,:)   ! POLCOMS layout: (zdim,londim,latdim)
        real(kind=8),pointer :: ETW(:,:,:)    ! POLCOMS layout: (zdim,londim,latdim)
        real(kind=8),pointer :: nuv(:,:,:)    ! POLCOMS layout: (zdim,londim,latdim)
        real(kind=8),pointer :: depth(:,:,:)  ! POLCOMS layout: (zdim,londim,latdim)
        real(kind=8),pointer :: zbnd(:,:,:,:) ! POLCOMS layout: (2,zdim,londim,latdim) !MOMME: please check
      end type

c     -------------------- module data --------------------  
      
c
c     ------ data frame handler ------
     
      logical     :: data_in_buffers = .false.  ! first time setting
      logical     :: data_in_sync    = .false.  ! first time setting
      real        :: fill_value_real = 0.0      ! assumed setting
      
c     --- 3D grids ---
c     --- 2D grids ---    
c     --- 1D grids ---
    
c     ===================================================
                            contains
c     ===================================================
  
      subroutine init_physical_fields(time)
c     ------------------------------------------
c     Do not trigger data load
c     ------------------------------------------
      type(clock), intent(in),optional :: time
      real,parameter :: very_deep = 1000.0
c     ------------------------------------------
      if (present(time)) call set_master_clock(time)
      write(*,*) trim(get_pbi_version()) 

      call read_grid_desc()  ! invokes init_horiz_grid_transf + alloc core arrays  
      
c     ---- specials for this data set / version ----
      
      dslm             = 0.0        ! NA for sigma grids
      bottom_layer     = 0          ! default dry
      
      ccdepth(:,:,1)   = very_deep 
      ccdepth(:,:,2:)  = 2*very_deep
      acc_width(:,:,1) = 0.0        ! including surface layer (iz=1)
      acc_width(:,:,2) = 2*very_deep
      acc_width(:,:,3:)= 2*very_deep
      
      zoo  = 0.0 
      
      end subroutine init_physical_fields

 

      character*100 function get_pbi_version()  
      get_pbi_version =  "POLCOMS+ERSEM online pbi version: $Rev: 228 $"
      end function



      subroutine close_physical_fields()
c     ------------------------------------------ 
c     ------------------------------------------    
      call close_mesh_grid()
      call close_horiz_grid_transf()

      data_in_buffers = .false.
      data_in_sync    = .false.
c     ------------------------------------------
      end subroutine 



      subroutine update_physical_fields(time, dt)
c     ------------------------------------------
c     Data must have been pushed by transfer_data (data_in_buffers == TRUE)
c     ------------------------------------------  
      type(clock), intent(in),optional :: time
      integer, intent(in),optional     :: dt
      logical                          :: update
      integer                          :: frame ! number in that set
      type(clock), pointer             :: aclock

      integer :: i,j
c     ------------------------------------------  
      aclock => get_master_clock()
      if (present(time)) then
         call set_master_clock(time)
      elseif (present(dt)) then
         call add_seconds_to_clock(aclock, dt)
      endif    
      if (.not.data_in_buffers) then
         write(*,*) "update_physical_fields: data_in_buffers is FALSE"
         stop
      endif
      if (.not.data_in_sync) call syncronize_data()
c     ------------------------------------------       
      end subroutine update_physical_fields





      subroutine transfer_data_to_ibmlib(pol_data)
c     ------------------------------------------  
c     Push data to raw buffers (data_in_buffers -> .true.) 
c
c     Reshape and cast to single precision at this point,
c     but do not invoke data syncronization (data_in_sync -> .false.)
c
c     fortran: fastests index left most
c     ------------------------------------------   
      type(physical_data_bucket) :: pol_data
      integer                    :: mlon,mlat,mz
c     ------------------------------------------
c
c     --- probe array dimensions and test grid match ---
c
      mlon = size(pol_data%ucur, dim=2)
      mlat = size(pol_data%ucur, dim=3)
      mz   = size(pol_data%ucur, dim=1)
      if (mlon /= nx) then
         write(*,*) "transfer_data_to_ibmlib: grid mismatch (lon)"
         write(*,*) "mlon =", mlon, "/= nx =", nx
         stop
      endif
      if (mlat /= ny) then
         write(*,*) "transfer_data_to_ibmlib: grid mismatch (lat)"
         write(*,*) "mlat =", mlat, "/= ny =", ny
         stop
      endif
      if (mz   /= nz) then
         write(*,*) "transfer_data_to_ibmlib: grid mismatch (vertical)"
         write(*,*) "mz =", mz, "/= nz =", nz
         stop
      endif
c
c     --- reshape + recast ---
c
      u = real(reshape(pol_data%ucur,   (/mlon,mlat,mz/),
     +                 order=(/2,3,1/)), kind=4)
      v = real(reshape(pol_data%vcur,   (/mlon,mlat,mz/),
     +                 order=(/2,3,1/)), kind=4)
      w = real(reshape(pol_data%wcur,   (/mlon,mlat,mz/),
     +                 order=(/2,3,1/)), kind=4)
      temp = real(reshape(pol_data%ETW, (/mlon,mlat,mz/),
     +                    order=(/2,3,1/)), kind=4)
      vdiffus(:,:,1:nz) = real(reshape(pol_data%nuv, 
     +                         (/mlon,mlat,mz/),
     +                          order=(/2,3,1/)), kind=4)
      ccdepth = real(reshape(pol_data%depth, (/mlon,mlat,mz/),
     +                 order=(/2,3,1/)), kind=4)
c     ... no reshape necessary for zbnd -> wdepth, just type cast
      wdepth  = real(pol_data%zbnd(2,mz,:,:), kind=4) ! pick lower face, bottom layer
c
c     --- update data handlers --- 
c
      data_in_buffers = .true. 
      data_in_sync    = .false. ! allow asyncroneous transfer
c
      end subroutine transfer_data_to_ibmlib
      




      subroutine read_grid_desc()
c     ---------------------------------------------------
c     Read grid definition which is contained in file 
c     given as tag grid_desc in the simulation file simulation_file
c     ---------------------------------------------------
      character*999        :: gd_fname
      type(control_file)   :: grid_ctrlfile
      real                 :: lambda1, dlambda, phi1, dphi ! LOCAL DUMMIES
c     ---------------------------------------------------
      if (data_in_buffers) stop "read_grid_desc: unexpected"

      call read_control_data(simulation_file, "grid_desc", gd_fname)
      write(*,229) "read_grid_desc: reading grid descriptor from: ", 
     +           trim(adjustl(gd_fname))

      call open_control_file(gd_fname, grid_ctrlfile)
    
c.....grid scale/dimensions (held globally in module regular_lonlat_grid)
      
      call read_control_data(grid_ctrlfile,"lambda1", lambda1)
      call read_control_data(grid_ctrlfile,"dlambda", dlambda) 
      call read_control_data(grid_ctrlfile,"phi1",    phi1)
      call read_control_data(grid_ctrlfile,"dphi",    dphi) 
      call read_control_data(grid_ctrlfile,"nx",      nx) ! nx hosted in horizontal_representation
      call read_control_data(grid_ctrlfile,"ny",      ny) ! ny hosted in horizontal_representation
      call read_control_data(grid_ctrlfile,"nz",      nz) ! nz hosted in mesh_grid
      call close_control_file(grid_ctrlfile)

      write(*,231) nx, ny, nz
      write(*,232) lambda1, dlambda
      write(*,233) phi1,    dphi 
    
 229  format(a,a)    
 231  format("read_grid_desc: 3d grid dim (nx,ny,nz) = ", i4,i4,i4)
 232  format("read_grid_desc: lambda1 = ",f12.7,
     +                      " dlambda = ",f12.7," degrees")  
 233  format("read_grid_desc: phi1    = ",f12.7,
     +                      " dphi    = ",f12.7," degrees")  
     

      write(*,*) "read_grid_desc: horizontal initialization"
      call init_horiz_grid_transf(lambda1, phi1, dlambda, dphi)


      write(*,*) "read_grid_desc: allocate grid arrays: begin"
      call init_mesh_grid()  ! incl allocation of 3D arrays
  
      
      end subroutine read_grid_desc
      



      subroutine syncronize_data()
c     ------------------------------------------
c     Syncronize auxillary grid fields after data load. 
c     It is assumed that grids do not
c     change during a simulation (assertion not checked)
c     ------------------------------------------    
      integer :: ix,iy
      real,pointer  :: cc(:), acc(:)
c     ------------------------------------------ 
      write(*,*) "syncronize_data: begin"

c     ---- render all elements defined, only 1:nz loaded
      
      vdiffus(:,:,nz+1) = 1.e-9 
      where (vdiffus<1.e-9)
         vdiffus = 1.e-9
      end where
c
c          NB: horizontal derivatives NOT implemented
c          when spatially non constant hdiffus are applied,
c          interpolate_turbulence_deriv must be updated
c          currently no horizontal turbulent diffusivity: set to
c          molecular diffusivity lower limit ~ 1.e-9 m2/s   --

      hdiffus = 1.e-9   
  

c     ---- cell center depths are positive below water, so flip sign 
c          of input data, which has opposite convention
 
      ccdepth = -ccdepth  ! flip sign so wdepth is positive in water

c     ---- grid is sigma type so all vertical cells are wet, if any
c          depths are positive below water, so flip sign of input data, which
c          has opposite convention
c     
      do ix=1,nx
         do iy=1,ny
            if (abs(wdepth(ix,iy))<1.e-6) wdepth(ix,iy) = 1.0 ! we flip sign below
            if (is_fill_value(wdepth(ix,iy)))
     +          wdepth(ix,iy) = 1.0                           ! we flip sign below
         enddo
      enddo
      wdepth = -wdepth ! flip sign of whole array so wdepth is positive in water

c     --- define the wetmask auxillary  ---    
      where(wdepth>0)
         wetmask      = 1 ! wet
         bottom_layer = nz
      elsewhere
         wetmask      = 0 ! dry
         bottom_layer = 0
      end where 

c     ------ generate auxillary grid descriptors for wet points  
c            acc_width(ix,iy,1) = 0 (sea surface)  
      do ix=1,nx
         do iy=1,ny
            if (wetmask(ix,iy)>0) then
               cc  => ccdepth(ix,iy,:)
               acc => acc_width(ix,iy,:)
               call levels2accwidth(cc, acc)
            endif      
         enddo
      enddo
c
c     ---- enforce exact consistency between acc_width(:,:,nz+1) and wdepth 
c          (there are numerical diffs of order 10^-5)
c
      acc_width(:,:,nz+1) = wdepth 
           

      write(*,*) "syncronize_data: end"
      data_in_sync = .true.  ! update data handler

c$$$      do ix=1,nx
c$$$         do iy=1,ny
c$$$            if (wetmask(ix,iy)==0) 
c$$$     +          write(77,*) (ix-1)*dlambda + lambda1,
c$$$     +                     (iy-1)*dphi    + phi1
c$$$            if (wetmask(ix,iy)==0) write(77,*) ix,iy
c$$$         enddo
c$$$      enddo
c$$$      stop 84465

      end subroutine syncronize_data
      



      subroutine levels2accwidth(levels, faces)
c     ------------------------------------------
c     Compute position of layer faces from center levels
c     ------------------------------------------
      real, intent(in)  :: levels(:)  ! 1:nz
      real, intent(out) :: faces(:)   ! 1:nz+1
      integer :: iz
      real    :: layer_width
c     ------------------------------------------ 
      faces(1) = 0
      do iz = 1,size(levels)
         layer_width   = 2.0*(levels(iz) - faces(iz))
         faces(iz+1)   = faces(iz) + layer_width
      enddo
      end subroutine levels2accwidth




      logical function is_fill_value(rval)
c     ------------------------------------------
c     query function can be overloaded later on
c     Assume fill_value /= 0
c     ------------------------------------------
      real, intent(in) :: rval
      real             :: reldeff
c     ------------------------------------------ 
      reldeff = abs((rval - fill_value_real)/fill_value_real)
      is_fill_value = (abs(reldeff) < 1.e-6)
      end function is_fill_value

      end module

      
