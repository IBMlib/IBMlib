      module read_cmod_ergom
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Hacked version of read_cmod_ergom.f (svn rev 453) which loads silicate
c     as salinity from a parallel data set. Template for easy hacking.
c     Changes in 
c         read_frame(request_time)
c         get_s
c     Data Reading module for HBM/ERGOM time series provided in OPEC 
c     derived from read_cmod.f, based on readtmp.f90 provided summer 2011 by 
c     Mikhail Dobrynin at DMI and bio2cdf.f90 + cmod2cdf.f90 provided 2012 + 2013 by 
c     Zhenwen Wan at DMI. HBM was earlier called cmod.
c     Preserve units applied by HBM/ERGOM. 
c     
c     Main public arrays: area refers to subgrid number 
c       m1(area)%p(:,:,:):    3D to compressed 1D map of wet grid points (cmod indexing)
c       hz(area)%p(:)    :    cell center depth 
c       cellw(area)%p(:) :    cell width
c
c       integer data buffers: (loaded by read_frame)
c         u_i4(area)%p(:)     :    u current
c         v_i4(area)%p(:)     :    v current
c         w_i4(area)%p(:)     :    w current
c         z_i4(area)%p(:)     :    sea surf elev
c         s_i4(area)%p(:)     :    salinity
c         t_i4(area)%p(:)     :    temperature
c         wu_i4(area)%p(:)    :    wind u stress
c         wv_i4(area)%p(:)    :    wind v stress
c         eco_i4(area)%p(:,:) :    biogeochemistry
c
c         ERGOM State variables:
c          nh4(im,jm,km)  nh4                        unit == mmol/m3 (more precisely umol/kg)
c          no3(im,jm,km)  no3                        unit == mmol/m3 (more precisely umol/kg)
c          po4(im,jm,km)  po4                        unit == mmol/m3 (more precisely umol/kg)
c          dia(im,jm,km)  diatoms                    unit == mmol/m3 (more precisely umol/kg) 
c          fla(im,jm,km)  flagellate                 unit == mmol/m3 (more precisely umol/kg)
c          cya(im,jm,km)  cyanobacteria              unit == mmol/m3 (more precisely umol/kg)
c          zo1(im,jm,km)  micro-zooplankton          unit == mmol/m3 (more precisely umol/kg)
c          zo2(im,jm,km)  meso-zooplankton           unit == mmol/m3 (more precisely umol/kg)
c          odt(im,jm,km)  organic detritus           unit == mmol/m3 (more precisely umol/kg)
c          oxy(im,jm,km)  dissolved oxygen           unit == mmol/m3 (more precisely umol/kg)
c          pom(im,jm,km)  particular organic matter  unit == mmol/m3 (more precisely umol/kg)
c          dic(im,jm,km)  dissolved inorganic carbon unit == mmol/m3 (more precisely umol/kg)
c          alk(im,jm,km)  alkalinity                 unit == mmol/m3 (more precisely umol/kg)
c 
c         ERGOM Derived variables:
c          din(im,jm,km)  dissolved inorganic nitrogen unit == mmol/m3 (more precisely umol/kg)
c          chl(im,jm,km)  chlorophyl                   unit == mmol/m3 (more precisely umol/kg) 
c
c         Unpacking of ERGOM buffers: 
c            nh4(i,j,k)=eco(1,id3d(jm+1-j,i,k))/scale                
c            no3(i,j,k)=eco(2,id3d(jm+1-j,i,k))/scale
c            po4(i,j,k)=eco(3,id3d(jm+1-j,i,k))/scale                
c            dia(i,j,k)=eco(4,id3d(jm+1-j,i,k))/scale
c            fla(i,j,k)=eco(5,id3d(jm+1-j,i,k))/scale
c            cya(i,j,k)=eco(6,id3d(jm+1-j,i,k))/scale
c            zo1(i,j,k)=eco(7,id3d(jm+1-j,i,k))/scale
c            zo2(i,j,k)=eco(8,id3d(jm+1-j,i,k))/scale
c            odt(i,j,k)=eco(9,id3d(jm+1-j,i,k))/scale
c            oxy(i,j,k)=eco(10,id3d(jm+1-j,i,k))/scale
c            pom(i,j,k)=eco(11,id3d(jm+1-j,i,k))/scale
c            dic(i,j,k)=eco(12,id3d(jm+1-j,i,k))/scale
c            alk(i,j,k)=eco(13,id3d(jm+1-j,i,k))/scale
cc           --- set derived variables 
c            din(i,j,k)=no3(i,j,k) + nh4(i,j,k)  
c            chl(i,j,k)=2.0*(dia(i,j,k)+fla(i,j,k)+cya(i,j,k))   
c
c     changes: 
c        May 2013: implemented format change in ERGOM data frames, adding more state variables; read test on data_sample2 OK
c
c        free format -> fixed format
c        removed station+interpolation+nesting logic
c        removed exitme
c        changed unit cm/s -> m/s for currents
c        added array cellw giving cell widths in meters  
c
c     Compile:
c        ifort -c time_tools.f
c        ifort -convert big_endian -c read_cmod_ergom.f
c        ifort -convert big_endian time_tools.o read_cmod_ergom.f libtime/libtime77.a
c
c     Notes on conventions/units/arrays etc.
c        m1(iarea)%p(jy,jx,jz) is 3D to 1D map of cell (jx,jy,jz) -> idx (notice switch of x,y)
c        x scans west    -> east
c        y scans north   -> south 
c        z scans surface -> bottum
c        so that cmod grids start in upper left corner
c        wet pixels: m1(iarea)%p(jy,jx,jz) > 0
c        IBMlib mesh grid reindexing:  (jy,jx,jz)  ->  (ny+1-iy,ix,iz)
c        surface layer: m1(iarea)%p(jy,jx,1)
c        current staggering: 
c             u on east  side of cell
c             v on south side of cell
c        sign of z validated from tide table Tide07_DK.pdf close to Esbjerg
c
c     TODO: set undefs at dry points, define undef test condition
c 
c     LOG: Feb 2014: enabled physics-only mode (ecosystem variables not loaded/allocated)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use time_tools
      use netcdf  ! for HACK

c     directives ----------------------------------------------------------------
      implicit none
      private  ! default scope
c
c     local type defs -----------------------------------------------------------------
c
      type cmi1
      integer(4), pointer :: p(:)
      end type
      type cmi2
      integer(4), pointer :: p(:,:)
      end type
      type cmi3
      integer(4), pointer :: p(:,:,:)
      end type
      type cmr1
      real(8),    pointer :: p(:)
      end type
c
c     limits and other parameters -----------------------------------------------
c
      integer(4)            :: narea  = 2        ! number of ocean model domains  
      integer(4), parameter :: neco   = 13       ! number of biogeochem 3D stored variables
      real,       parameter :: dryval      = 0.0        ! assign this to dry cells in get_X
      real(8),    parameter :: undef       = 9999.9     ! assign this to variable
      real(8),    parameter :: test_undef  = 0.99*undef ! for test: if u>test_undef ...
      logical,    parameter :: ldebug = .false.  ! 


      character(len=3), parameter :: eco_labels(neco) =    ! there is only a slot here to raw variables
     +     (/"nh4", "no3", "po4", "dia", "fla",            ! further labels are resolved by get_eco_3D
     +       "cya", "zo1", "zo2", "odt", "oxy",            ! for derived variables
     +       "pom", "dic", "alk"/)

c
c     allocatable arrays --------------------------------------------------------
c
      integer(4),     allocatable, dimension(:)   :: mmx,nmx,kmx,iw2,iw3 ! dimensions 
      type (cmi3),    pointer,     dimension(:)   :: m1
      type (cmr1),    pointer,     dimension(:)   :: hz,cellw
      type (cmi1),    allocatable, dimension(:)   :: u_i4,v_i4,w_i4,z_i4  ! integer rep buffers for physics
      type (cmi1),    allocatable, dimension(:)   :: s_i4,t_i4            ! integer rep buffers for physics
      type (cmi1),    allocatable, dimension(:)   :: wu_i4,wv_i4          ! integer rep buffers for physics
      type (cmi2),    allocatable, dimension(:)   :: eco_i4               ! integer rep buffers for biogeochemisty
      type (cmi1),    allocatable, dimension(:)   :: benthos_i4           ! integer rep buffers for biogeochemisty

      real(8),        allocatable, dimension(:)   :: lambda1, phi1
      real(8),        allocatable, dimension(:)   :: dlam, dphi  
      character(len=8), allocatable               :: setname(:)

      real, allocatable :: databuf(:,:,:)   ! HACK
c
c     -----------------------  uncompression scale factors -----------------------
c
      real, parameter :: wuvscale = 1.0e7    ! wind scalefactor
      real, parameter :: uvscale  = 1.0e8    ! current scalefac.
      real, parameter :: stscale  = 1.0e7    ! salt/temp scfac.
      real, parameter :: wscale   = 1.0e11   ! vert.vel. scfac.
      real, parameter :: zscale   = 1.0e8    ! elevation scfac.
      real, parameter :: bscale   = 1.0e4    ! for biology conversion
c  
c     -----------------------  data handlers  ------------------------------------
c
      character(len=11)  :: current_tag      = " "  ! date tag of file set containing current frame    (same for phys+bio)
      integer            :: current_frame(4) = -1   ! current frames, as year, month, day, hour in day (same for phys+bio)
      integer            :: iunit_phy        = -1   ! assign at first read ()
      integer            :: iunit_bio        = -1   ! assign at first read ()
      character(len=999) :: hydroDBpath      = " "  ! including trailing /
      logical            :: include_bio      = .true. ! handle to control whether biodata is read (physics is always read)

      character(len=7), parameter :: physprefix = "physics"  ! new prefix May 2013
      character(len=6), parameter :: bioprefix  = "biodat"     
      character*(*), parameter    :: griddef    = 'grid_def.dat'
c      
c     -----------------------  set module public scope -----------------------
c
      public :: cmi1, cmi3, cmr1, cmi2
      public :: init_read_cmod_ergom
      public :: close_read_cmod_ergom
      public :: get_grid_descriptors
      public :: update_buffers
      public :: m1,hz,cellw
      public :: test_undef          
      public :: get_z, get_wu, get_wv              ! 2D physics 
      public :: get_u, get_v, get_w, get_s, get_t  ! 3D physics 
      public :: get_eco_3D                         ! 3D biogeochemistry
      public :: check_buffer_bounds
      public :: create_hydro_file_tag              ! open for debugging
      
c     ========================================================================
                                  contains 
c     ========================================================================


      subroutine init_read_cmod_ergom(path, readbio)
c----------------------------------------------------------------------------+
c     path is to grid descriptor files, without trailing slash
c     readbio controls wheter ecosystem varaible are read
c----------------------------------------------------------------------------+ 
      character*(*), intent(in) :: path ! applied as prefix to file names (without trailing slash)
      logical, intent(in)       :: readbio
      integer(4)                :: ios
      integer(4)                :: lun,ia,ie
      character(LEN=999)        :: cfgfile
      real(8)                   :: xmin, ymin, xmax, ymax
c----------------------------------------------------------------------------+ 
      hydroDBpath = trim(adjustl(path)) // "/" ! prepend linux delimiter
      write(*,*) "init_read_cmod_ergom: path = ", trim(hydroDBpath)
      write(*,*) "WARNING: HACKED VERSION ! (silicate shadow load)"

      include_bio = readbio
      if (include_bio) then
         write(*,*) "init_read_cmod_ergom: ecosystem data will be read"
      else
         write(*,*) "init_read_cmod_ergom: ecosystem data excluded"
      endif
         
c     allocate descriptor arrays:
      allocate( mmx(narea),nmx(narea),kmx(narea),iw2(narea),iw3(narea),          
     +     m1(narea),hz(narea),cellw(narea),
     +     lambda1(narea), phi1(narea), dlam(narea), dphi(narea),
     +     setname(narea))

c     -------------- Get the ocean model configurations ----------------
      lun = io_new_unit()
      open (unit=lun, file=trim(hydroDBpath)//griddef, 
     +      status='old', iostat=ios)
      if (ios /= 0) then
         write(*,*) "init_read_cmod_ergom: grid def file not found"
         write(*,*) "file name = ",trim(hydroDBpath)//griddef
         stop  ! removed default setting, require explicit input
      endif
      

      do ia=1,narea
         read(lun,*) nmx(ia), mmx(ia), kmx(ia)  ! reverserd order to x,y,z
         read(lun,*) lambda1(ia), dlam(ia)
         read(lun,*) phi1(ia), dphi(ia)
         read(lun,'(a8,a)') setname(ia), cfgfile
         setname(ia) = adjustl(setname(ia))
         read(lun,*)   ! read an empty line
c        ------  m1(ia)%p, hz(ia)%p, cellw(ia)%p are allocated in read_cmod_grid ------ 
         cfgfile = trim(hydroDBpath)//trim(adjustl(cfgfile))
         call read_cmod_grid(cfgfile,
     +           nmx(ia),mmx(ia),kmx(ia), ! NB: flipped order to xyz
     +           iw2(ia),iw3(ia),
     +           m1(ia)%p, hz(ia)%p, cellw(ia)%p) 
c 
c           -------- print grid info --------
c     
         xmin = lambda1(ia) 
         ymin = phi1(ia) 
         xmax = lambda1(ia) + nmx(ia)*dlam(ia)
         ymax = phi1(ia)    + mmx(ia)*dphi(ia)
         write(*,'(50("-"))') 
         write(*,*) "init_read_cmod_ergom: data set   : ", setname(ia)
         write(*,*) "init_read_cmod_ergom: grid number: ", ia
         write(*,*) "long range = [", xmin, xmax,"] ", "deg E"
         write(*,*) "lati range = [", ymin, ymax,"] ", "deg N"
         write(*,*) "dlon,dlat  = ", dlam(ia), dphi(ia), "deg"
         write(*,*) "grid dims  = ", nmx(ia), mmx(ia), kmx(ia)
         write(*,*) "init_read_cmod_ergom: grid config: ", 
     +                  trim(cfgfile)
         write(*,*) "surf/bulk wet points = ", iw2(ia),iw3(ia)  

      enddo
      close(unit=lun)      
      write(*,'(50("-"))') 

c     ---------- allocate hydrographic buffers ----------
      allocate(u_i4(narea), v_i4(narea), 
     +     w_i4(narea), z_i4(narea), s_i4(narea), t_i4(narea), 
     +     wu_i4(narea), wv_i4(narea))

      if (include_bio) allocate(eco_i4(narea), benthos_i4(narea)) 

      do ia=1,narea
         allocate( u_i4(ia)%p(0:iw3(ia)), 
     +        v_i4(ia)%p(0:iw3(ia)),w_i4(ia)%p(0:iw3(ia)), 
     +        z_i4(ia)%p(0:iw2(ia)),s_i4(ia)%p(0:iw3(ia)), 
     +        t_i4(ia)%p(0:iw3(ia)), wu_i4(ia)%p(0:iw2(ia)), 
     +        wv_i4(ia)%p(0:iw2(ia)) )
         
         if (include_bio) allocate( benthos_i4(ia)%p(0:iw2(ia))      )
         if (include_bio) allocate(     eco_i4(ia)%p(neco,0:iw3(ia)) )
      enddo

      allocate( databuf(nmx(1), mmx(1), kmx(1)) ) ! HACK - only valiod for area 1

      end subroutine init_read_cmod_ergom




      subroutine close_read_cmod_ergom()
c     -------------------------------------------------------------------
c     deallocate module data arrays/pointers
c     -------------------------------------------------------------------  
      integer(4)           :: ia
c     -------------------------------------------------------------------  
c     ---- physics ----
      deallocate( mmx,nmx,kmx,iw2,iw3 )
      do ia=1,narea          
         if (associated( m1(ia)%p ))    deallocate( m1(ia)%p )
         if (associated( hz(ia)%p ))    deallocate( hz(ia)%p )
         if (associated( cellw(ia)%p )) deallocate( cellw(ia)%p )
         if (associated( u_i4(ia)%p ))  deallocate( u_i4(ia)%p )
         if (associated( v_i4(ia)%p ))  deallocate( v_i4(ia)%p )
         if (associated( w_i4(ia)%p ))  deallocate( w_i4(ia)%p )
         if (associated( z_i4(ia)%p ))  deallocate( z_i4(ia)%p )
         if (associated( s_i4(ia)%p ))  deallocate( s_i4(ia)%p )
         if (associated( t_i4(ia)%p ))  deallocate( t_i4(ia)%p )
         if (associated( wu_i4(ia)%p )) deallocate( wu_i4(ia)%p )
         if (associated( wv_i4(ia)%p )) deallocate( wv_i4(ia)%p )
      enddo
      deallocate( m1 )
      deallocate( hz,cellw )
      deallocate( u_i4,v_i4,w_i4,z_i4 )  
      deallocate( s_i4,t_i4 )            
      deallocate( wu_i4,wv_i4 )

c     ---- biology ----
      if (include_bio) then
         do ia=1,narea 
            if (associated(benthos_i4(ia)%p)) 
     +         deallocate(benthos_i4(ia)%p)
            if (associated( eco_i4(ia)%p )) 
     +         deallocate( eco_i4(ia)%p)
         enddo   
         deallocate( benthos_i4, eco_i4)
      endif

      deallocate (databuf) ! HACK
c
c     --- reset data handlers ---
c
      current_tag   = " "  
      current_frame = -1   
      iunit_phy     = -1   
      iunit_bio     = -1   
      end subroutine close_read_cmod_ergom




      subroutine get_grid_descriptors(sname_in, data_set_handler,
     +                          lam1,dl, ph1,dp, nx,ny,nz) 
c     ----------------------------------------------------------------------------+
c     Interface query function to retrieve grid characteristics:
c
c     data_set_handler           : data set index <= narea
c     lambda1,dlambda, phi1,dphi : lower left grid point
c     nx,ny,nz                   : grid dimensions
c
c     ----------------------------------------------------------------------------+
      character*(*),intent(in) :: sname_in
      integer,intent(out)      :: data_set_handler ! area number for this set
      real,intent(out)         :: lam1,dl,ph1,dp
      integer,intent(out)      :: nx,ny,nz
      integer                  :: ncpm
      character(len=999)       :: sname
c     ----------------------------------------------------------------------------+  
      sname = trim(adjustl(sname_in)) ! avoid prefixed blanks gives trouble
      ncpm = min(len_trim(sname), len_trim(setname(1)))  
      do data_set_handler = 1, narea
         if (sname(1:ncpm) == setname(data_set_handler)(1:ncpm)) exit
      enddo
c     ---- standard exit condition for exhausted do loop: 
      if (data_set_handler == (narea+1)) then  
         write(*,*) "get_grid_descriptors : invalid data set name = ", 
     +              trim(sname)
         stop
      endif
      lam1 = lambda1(data_set_handler)
      dl   = dlam(data_set_handler)
      ph1  = phi1(data_set_handler)
      dp   = dphi(data_set_handler)
      nx   = nmx(data_set_handler)
      ny   = mmx(data_set_handler)
      nz   = kmx(data_set_handler)
      end subroutine get_grid_descriptors
      


      logical function update_buffers(now) 
c     ----------------------------------------------------------------------------+
c     Update buffers corresponding to requested time (year,month,day,hour)
c     Update is intelligent in the sense that file reading is only done, if needed
c     update_buffers returns TRUE, if data has changed since last call
c     so the return value can be used for starting update cascade and minimizing update overhead
c     Assume needed, if fname has changed (crude assumption)
c     Assume biogeochemical + physical data are syncronized file wise
c     ----------------------------------------------------------------------------+
      type(clock), intent(in)    :: now 
      integer                    :: year,month,day,hour
      integer                    :: request_time(4), tcur, treq, ios
      logical                    :: newfile, newtime
      character(len=999)         :: fname
      character(len=11)          :: tag     ! file name postfix
c     ----------------------------------------------------------------------------+
      call create_hydro_file_tag(now, tag, year, month, day, hour)
      request_time(1) = year
      request_time(2) = month
      request_time(3) = day
      request_time(4) = hour
      tcur = time_index( current_frame )
      treq = time_index( request_time  )
c
c     ---- if we need open a new file, prepare it for reading and close old
c          syncronize phys + bio data handling
c          
      newfile   = (trim(adjustl(tag)) /= trim(adjustl(current_tag)))
      if (newfile) then
         current_tag = tag
c        ---------- first physics ----------
         if (iunit_phy>0) close(iunit_phy)                               ! we opened it previously
         if (iunit_phy<0) iunit_phy = io_new_unit()                      ! first time, pick a free unit 
         fname = trim(hydroDBpath) // physprefix // tag          ! physprefix + tag are cropped
         open(iunit_phy,file=fname, form='unformatted',status='old',
     +               iostat=ios )                                ! keep same logical unit
         if (ios /= 0) then
            write(*,'(a,a)') "error opening data file ", 
     +                       trim(adjustl(fname))
            stop
         else
            write(*,'(a,a)') "update_buffers: opened new data file", 
     +                       trim(adjustl(fname))
         endif
         call reset_tempfile(iunit_phy) ! wind to first frame
c        ---------- then biochemistry ----------
         if (include_bio) then
            if (iunit_bio>0) close(iunit_bio) ! we opened it previously
            if (iunit_bio<0) iunit_bio = io_new_unit() ! first time, pick a free unit 
            fname = trim(hydroDBpath) // bioprefix // tag ! bioprefix + tag are cropped
            open(iunit_bio,file=fname, form='unformatted',status='old',
     +           iostat=ios )                                ! keep same logical unit
            if (ios /= 0) then
               write(*,'(a,a)') "error opening data file ", 
     +              trim(adjustl(fname))
               stop
            else
               write(*,'(a,a)') "update_buffers: opened new data file", 
     +              trim(adjustl(fname))
            endif
            call reset_tempfile(iunit_bio) ! wind to first frame
         endif   ! include_bio
c         
      endif
c      
      tcur = time_index( current_frame )
      treq = time_index( request_time  )
      newtime = (tcur /= treq)
c
c     decide to read or not and set update flag
c
      if (newtime.or.newfile) then  
         if (treq <= tcur) then
              call reset_tempfile(iunit_phy)                  ! will rewind if same frame is requested  
              if (include_bio) call reset_tempfile(iunit_bio) ! will rewind if same frame is requested  
         endif     
         call read_frame(request_time)       ! updates current_frame when loaded
         update_buffers = .true.
      else
         update_buffers = .false.
      endif

      end function update_buffers



      subroutine create_hydro_file_tag(clck,tag,year,month,day,hour)
c     ----------------------------------------------------------------------------+
c     Map a time point represented as onto a file tag
c     along with rounded time stamp (year,month,day,hour) with hour in [0;23]
c     Daily data are split in 4 chunks; data mapping: X = phydat/biodat
c        X00.YYYYMMDD: hour =  1, 2, 3, 4, 5, 6 @ YYYYMMDD
c        X06.YYYYMMDD: hour =  7, 8, 9,10,11,12 @ YYYYMMDD
c        X12.YYYYMMDD: hour = 13,14,15,16,17,18 @ YYYYMMDD
c        X18.YYYYMMDD: hour = 19,20,21,22,23    @ YYYYMMDD + 0 @ (YYYYMMDD + 1 day)
c     Notice that last frame in X18.YYYYMMDD (hour 0) carries date tag corresponding 
c     to YYYYMMDD + 1 day
c     ----------------------------------------------------------------------------+
      type(clock), intent(in)         :: clck
      character(len=11), intent(out)  :: tag
      integer, intent(out)            :: year,month,day,hour
      integer, parameter   :: qmap(0:23) = (/18, 
     +                                       0, 0, 0, 0, 0, 0, 
     +                                       6, 6, 6, 6, 6, 6,
     +                                       12,12,12,12,12,12,
     +                                       18,18,18,18,18   /)
      type(clock)          :: dummy_clock 
      integer              :: isec, yyyy,mm,dd
c     ----------------------------------------------------------------------------+      
      call set_clock(dummy_clock, clck)
      call add_seconds_to_clock(dummy_clock, -1800) ! allow fetching date part for tag
      call get_date_from_clock(dummy_clock, yyyy, mm, dd) ! for tag
      call set_clock(dummy_clock, clck)
      call add_seconds_to_clock(dummy_clock, 1800)  ! allow rounding by truncation
      call get_date_from_clock(dummy_clock, year, month, day)
      call get_second_in_day(dummy_clock, isec)
      hour = int(isec/3600) ! cast by truncation
      write(tag,433) qmap(hour),yyyy, mm, dd
 433  format(i2.2,".",i4.4,i2.2,i2.2)
c     ----------------------------------------------------------------------------+
      end subroutine create_hydro_file_tag


      subroutine reset_tempfile(lun)
c     ----------------------------------------------------------------------------+
c     Internal method to rewind and advance file pointer to first frame 
c     by reading "aufdat" header (just waste it)
c     Robust toward invoking it multiple times 
c     ----------------------------------------------------------------------------+
      integer,intent(in) :: lun
      character(len=11)  :: aufdat
c     ----------------------------------------------------------------------------+
      rewind(lun)
      read(lun) aufdat ! advance file pointer to first frame 
      write(*,*) "reset_tempfile: ready to read first frame"
      end subroutine reset_tempfile


      

      subroutine read_frame(request_time)
c     ----------------------------------------------------------------------------+
c     --------------------------------- READ DATA ---------------------------------
c     ----------------------------------------------------------------------------+
c
c     Verbose load of time frame request_time from opened tempfile 
c     associated with file handlers (iunit_phy,iunit_bio). 
c     read_frame will continue reading from current reading position
c     (optimal for forward time simulations). Requested time should be after
c     current reading position. It is the task of the calling unit 
c     attach lun to the appropriate file in proper mode, e.g 
c
c        open(lun,file=filename, form='unformatted',status='old')
c
c     When a new file is attached, reset_tempfile should always be invoked.
c
c     frame header example:
c
c      aufdat = V200103.154
c      ctim   = 2001.03.15 13:00:00
c      tim    = 36963.5416666668    = days since 1900.01.01 00:00:00
c
c     physics and biochemistry files have same layout; file content:
c
c         aufdat             
c         ctim1,tim1
c         data_segment1 (N=2 areas)
c         ctim2,tim2
c         data_segment2 (N=2 areas)
c         ...
c     it is assumed that data is time ordered: tim1 < tim2 < ...
c     ----------------------------------------------------------------------------+ 
      integer,intent(in)       :: request_time(4) ! year, month, day, hour in day
 
      integer(4)               :: ia,it,ios,istat
      integer(4)               :: ncid,varid
      character(len=11)        :: aufdat
      character(len=19)        :: ctim
      real(8)                  :: tim
      integer                  :: treq, tcur
      integer                  :: time_stamp(4)
      character(len=256)       :: fname
c     ----------------------------------------------------------------------------+      
c
c     ======== load physics: read time stamp from tempdat file: ========
c
      read(iunit_phy, iostat=ios) ctim,tim
      it = ios
      do while (it == 0)
c
         read(ctim,532) time_stamp
         if (ldebug) write(*,*) "reading phys: ", ctim
c        read vars from tempdat file:
         read(iunit_phy) (wu_i4(ia)%p, wv_i4(ia)%p, ia=1,narea)
         do ia=1,narea
            read(iunit_phy) u_i4(ia)%p, v_i4(ia)%p,
     +                      w_i4(ia)%p, z_i4(ia)%p
         enddo
         do ia=1,narea
            read(iunit_phy) s_i4(ia)%p, t_i4(ia)%p
         enddo

         if ( is_equal(time_stamp, request_time) ) then
            write(*,*) "read_frame: loaded physics frame ", ctim
            current_frame = time_stamp ! update data handler
            goto 744                   ! we got what we looked for - proceed to biogeochem
         endif                
         
c     in readtmp, here was the stations loop - removed all this (asc.09Sep2011)
c     read next time stamp from tempdat file:

         read(iunit_phy,iostat=ios) ctim,tim
         it = ios
      enddo                     ! it-loop (time frames in tempdat file)

c     capture unexpected case where requested data
c     is not in file, or data is otherwise corrupted
c     At successful read, we will exit in the do-while loop above     

      write(*,*) "read_frame: did not find requested physics time frame"
      write(*,*) "read_frame: requested time = ", request_time
      write(*,*) "read_frame: last record    = ", ctim
      stop ! fatal condition

 744  continue  ! entry point after reading physics

c
c     HACK: load silicate frame from daily means, bundled into monthly files
c
      fname = "biodat.H4_si_2005020100.daily_mean.0.nc"      ! HACK
      write(fname(14:17),'(i4)')   request_time(1) ! year    ! HACK
      write(fname(18:19),'(i2.2)') request_time(2) ! month   ! HACK
      istat = nf90_open(fname, NF90_NOWRITE, ncid)           ! HACK
      if (istat /= NF90_NOERR) stop "nf90_open: error"       ! HACK
      istat = nf90_inq_varid(ncid, 'si', varid)              ! HACK
      if (istat /= NF90_NOERR) stop "nf90_inq_varid: error"  ! HACK
      istat = nf90_get_var(ncid, varid, databuf,             ! HACK
     +                     start=(/1,1,1,request_time(3)/))  ! HACK
      if (istat /= NF90_NOERR) stop "nf90_get_var: error"    ! HACK
      istat = nf90_close(ncid)                               ! HACK
      if (istat /= NF90_NOERR) stop "nf90_close: error"      ! HACK

c
c     ======== load biogeochemistry: read time stamp from tempdat file: ========
c
      if (.not. include_bio) then
         write(*,*) "read_frame: omitted ecosystem vars ", ctim
         return    
      endif

      read(iunit_bio, iostat=ios) ctim,tim
      it = ios
      do while (it == 0)
c
         read(ctim,532) time_stamp
         if (ldebug) write(*,*) "reading bio: ", ctim
c        read vars from tempdat file:
         do ia=1,narea
            read(iunit_bio) eco_i4(ia)%p
            read(iunit_bio) benthos_i4(ia)%p
         enddo

         if ( is_equal(time_stamp, request_time) ) then
            write(*,*) "read_frame: loaded biogeo frame ", ctim
            current_frame = time_stamp ! update data handler
            return                   ! we got all (phys+bio) what we looked for - return
         endif                
         
c     in readtmp, here was the stations loop - removed all this (asc.09Sep2011)
c     read next time stamp from tempdat file:

         read(iunit_bio,iostat=ios) ctim,tim
         it = ios
      enddo                     ! it-loop (time frames in tempdat file)

c     capture unexpected case where requested data
c     is not in file, or data is otherwise corrupted
c     At successful read, we will exit in the do-while loop above     

      write(*,*) "read_frame: did not find requested biogeo time frame"
      write(*,*) "read_frame: requested time = ", request_time
      write(*,*) "read_frame: last record    = ", ctim
      stop ! fatal condition

 532  format(i4,1x,i2,1x,i2,1x,i2)

c     ----- local functions/subroutines -----
      contains
c     ----- local functions/subroutines -----


      logical function is_equal(iv1,iv2)
c     -------------------------------------------
c     Compare (iv1,iv2) element-by-element
c     F90 standard does not invoke this when seeing
c     the expression "iv1 == iv2"
c     -------------------------------------------
      integer,intent(in) :: iv1(:), iv2(:)
      integer            :: i
      if (size(iv1) /= size(iv2)) then
         is_equal = .false.
         return
      else
         is_equal =  .true. ! capture a null vector, if not captured as runtime exception
         do i=1,size(iv1)
            if (iv1(i) /= iv2(i)) then
               is_equal = .false.
               return ! do not process remaining elements
            endif 
         enddo
      endif
      end function is_equal ! local to read_frame

      end subroutine read_frame




      integer(8) function time_index(i4)
c     -------------------------------------------
c     Map a i4 vector (year, month, day, hour in day)
c     to a time index so that time_index(time1) > time_index(time2)
c     if and only if time1 > time2.
c     This is a easy way to circumvent leap-year/varying month length complication
c     The price is that time_index have jumps at year/month shifts
c     -------------------------------------------
      integer, intent(in) :: i4(4)
      integer, parameter  :: day2h   = 24
      integer, parameter  :: month2h = 31*day2h   ! multiplication factor must be >= 31  
      integer, parameter  :: year2h  = 12*month2h ! multiplication factor must be >= 366
c     -----------------------------------------      
      time_index = (i4(1)-2000.0)*year2h + (i4(2)-1)*month2h + 
     +              (i4(3)-1)*day2h + i4(4)
      end function time_index



      subroutine read_cmod_grid(fname,nx,ny,nz,iw2,iw3,m,ccd,cwd)    ! NB: flipped order to xyz
c----------------------------------------------------------------------------+
c     Merging of hzmod and read_wet from readtmp.f90
c----------------------------------------------------------------------------+  
      integer,      intent(in)  :: nx,ny,nz
      integer,      intent(out) :: iw2,iw3
      character(*), intent(in)  :: fname
      
      integer, pointer :: m(:,:,:)   ! 3D to 1D map        - allocated here
      real(8), pointer :: ccd(:)     ! cell center depths  - allocated here
      real(8), pointer :: cwd(:)     ! cell widths         - allocated here
   
      integer :: lun, ios
      integer :: i,j,l,l1,m1
      real(8) :: hs
c----------------------------------------------------------------------------+ 
c     alloc index array:
      allocate( m(ny,nx,nz) ) ! notice order

c     read index array from file:
      lun = io_new_unit()
      open(lun,file=trim(fname),iostat=ios)
      if (ios /= 0) then
         write(*,*) 'Error opening '//trim(fname)
         stop
      endif

      read(lun,*,iostat=ios) m
      if (ios /= 0) then
         write(*,*) "read_cmod_grid: error reading 3D to 1D map"
         write(*,*) "file name = ",trim(fname)
         stop  ! removed default setting, require explicit input
      endif



c     obtain No of 3D wet points and surface wet points:
      iw3 = maxval(m)
      iw2 = maxval(m(:,:,1))     ! scanned layer wise, starting at surface layer

c     allocate and read depth array:
      allocate( cwd(0:iw3) )
      read(lun,'(15f8.2)',i ostat=ios) cwd
      if (ios /= 0) then
         write(*,*) "read_cmod_grid: error cell width array"
         write(*,*) "file name = ",trim(fname)
         stop  ! removed default setting, require explicit input
      endif

      close(lun)
      
      allocate( ccd(0:iw3) )
c     
c     transform cell width array cwd(:) to cell center depths ccd(:)
c     
      allocate( ccd(0:iw3) )
      ccd(0) = cwd(0)    
      
c     flipped i,j from original loops

      do j=1,ny
         do i=1,nx
            do l=1,nz
               m1 = m(j,i,l)
               if (m1 /= 0) then
                  hs = 0.0
                  do l1=1,l
                     hs = hs + cwd(m(j,i,l1))
                  enddo
                  ccd(m1) = hs - cwd(m1)*0.5
               endif
            enddo
         enddo
      enddo

      end subroutine read_cmod_grid
c
c     ================ data fetch operators: get_X ================
c

      subroutine get_z(x, iset)    ! sea surface elevation [meters] - positive up (sign validated from tide table Tide07_DK.pdf)
c     -------------------------------------------------------
      real, intent(out)   :: x(:,:)
      integer, intent(in) :: iset   
c     -------------------------------------------------------
      call get_X_2D(x, iset, z_i4, zscale)
      end subroutine get_z

      subroutine get_wu(x, iset)    ! zonal wind stress [ N/m2 ] - positive toward East
c     -------------------------------------------------------
      real, intent(out)   :: x(:,:)
      integer, intent(in) :: iset   
c     -------------------------------------------------------
      call get_X_2D(x, iset, wu_i4, wuvscale)
      end subroutine get_wu


      subroutine get_wv(x, iset)    ! meridonal wind stress [ N/m2 ] - positive toward North
c     -------------------------------------------------------
      real, intent(out)   :: x(:,:)
      integer, intent(in) :: iset   
c     -------------------------------------------------------
      call get_X_2D(x, iset, wv_i4, wuvscale)
      end subroutine get_wv

      subroutine get_u(x, iset)    ! zonal current [meters/sec] - positive toward East
c     -------------------------------------------------------
      real, intent(out)   :: x(:,:,:)
      integer, intent(in) :: iset   
c     -------------------------------------------------------
      call get_X_3D(x, iset, u_i4, uvscale)
      end subroutine get_u


      subroutine get_v(x, iset)    ! meridonal current [meters/sec] - positive toward North
c     -------------------------------------------------------
      real, intent(out)   :: x(:,:,:)
      integer, intent(in) :: iset   
c     -------------------------------------------------------
      call get_X_3D(x, iset, v_i4, uvscale)
      end subroutine get_v


      subroutine get_w(x, iset)    ! vertical current [meters/sec] - positive up
c     -------------------------------------------------------
      real, intent(out)   :: x(:,:,:)
      integer, intent(in) :: iset   
c     -------------------------------------------------------
      call get_X_3D(x, iset, w_i4, wscale)
      end subroutine get_w

   

      subroutine get_s(x, iset)    !  (previously salinity [PSU])
c     -------------------------------------------------------      
c     HACKED VERSION - shadow load silicate frame as salinity
c     -------------------------------------------------------
      real, intent(out)   :: x(:,:,:)
      integer, intent(in) :: iset   
c     -------------------------------------------------------
      x = databuf ! HACK same layout
      end subroutine get_s


      subroutine get_t(x, iset)    ! temperature [*C]
c     -------------------------------------------------------
      real, intent(out)   :: x(:,:,:)
      integer, intent(in) :: iset   
c     -------------------------------------------------------
      call get_X_3D(x, iset, t_i4, stscale)
      end subroutine get_t

c 
c     --------------- generic 2D uncompression  --------------------
c

      subroutine get_X_2D(x, iset, i4buf, scale)
c     -------------------------------------------------------
c     Generic transformation from integer dry rep i4buf 
c     in set iset to 2D full rep x using rescaling by scale
c     -------------------------------------------------------   
      real, intent(out)       :: x(:,:)
      integer, intent(in)     :: iset   ! which area (data set) to pick
      type (cmi1), intent(in) :: i4buf(:)
      real, intent(in)        :: scale  ! int -> real transformation scale
      integer                 :: i,j,m
c     -------------------------------------------------------  
      if ((iset < 1).or.(iset > narea)) then
         write(*,*) "get_X_2D: invalid set", iset
      endif
      x = dryval
      do j=1,mmx(iset)
         do i=1,nmx(iset)
            m = m1(iset)%p(mmx(iset)+1-j,i,1)
            if(m > 0) x(i,j)=i4buf(iset)%p(m)/scale
         enddo
      enddo
      end subroutine get_X_2D

c 
c     --------------- generic 3D uncompression  --------------------
c

      subroutine get_X_3D(x, iset, i4buf, scale)
c     -------------------------------------------------------
c     Generic transformation from integer dry rep i4buf 
c     in set iset to 3D full rep x using rescaling by scale
c     -------------------------------------------------------   
      real, intent(out)       :: x(:,:,:)
      integer, intent(in)     :: iset     ! which area (data set) to pick
      type (cmi1), intent(in) :: i4buf(:)
      real, intent(in)        :: scale    ! int -> real transformation scale
      integer                 :: i,j,k,m
c     -------------------------------------------------------  
      if ((iset < 1).or.(iset > narea)) then
         write(*,*) "get_X_3D: invalid set", iset
      endif
      x = dryval
      do k=1,kmx(iset)
         do j=1,mmx(iset)
            do i=1,nmx(iset)
               m = m1(iset)%p(mmx(iset)+1-j,i,k)
               if(m > 0) x(i,j,k) = i4buf(iset)%p(m)/scale
            enddo
         enddo
      enddo
      end subroutine get_X_3D


      integer function get_eco_slot(property)
c     -------------------------------------------------------
c     Locate slot number of property in eco_labels
c     If property is absent return -1
c     -------------------------------------------------------
      character(len=3), intent(in) :: property
c     -------------------------------------------------------
      do get_eco_slot = 1, neco
         if (eco_labels(get_eco_slot) == property) return
      enddo
      get_eco_slot = -1 ! signal property is absent
      end function get_eco_slot


      recursive subroutine get_eco_3D(x, iset, property)
c     -------------------------------------------------------
c     Generic transformation from integer dry rep eco_i4 variable
c     corresponding to property in set iset to 3D full rep x 
c     return with x == dryval if include_bio is false
c     -------------------------------------------------------   
      real, intent(out)            :: x(:,:,:)
      integer, intent(in)          :: iset    ! which area (data set) to pick
      character(len=3), intent(in) :: property
      integer                      :: i,j,k,m, slot
      real, allocatable            :: x1(:,:,:), x2(:,:,:), x3(:,:,:)
c     -------------------------------------------------------  
      if ((iset < 1).or.(iset > narea)) then
         write(*,*) "get_eco_3D: invalid set", iset
         stop
      endif
      
      if (.not. include_bio) then
         if (ldebug) write(*,*) "get_eco_3D: include_bio == false"
         x = dryval 
         return  ! ecosystem arrays below are not allocated
      endif

c     --- first handle derived properties      
      if      (property == "din") then      ! dissolved inorganic nitrogen
         allocate( x1(size(x,1), size(x,2), size(x,3)) )
         allocate( x2(size(x,1), size(x,2), size(x,3)) )
         call get_eco_3D(x1, iset, "no3")
         call get_eco_3D(x2, iset, "nh4")
         x = x1 + x2
         deallocate( x1 )
         deallocate( x2 )
         return
      elseif  (property == "chl") then      ! chlorophyl
         allocate( x1(size(x,1), size(x,2), size(x,3)) )
         allocate( x2(size(x,1), size(x,2), size(x,3)) )
         allocate( x3(size(x,1), size(x,2), size(x,3)) )
         call get_eco_3D(x1, iset, "dia")
         call get_eco_3D(x2, iset, "fla")
         call get_eco_3D(x3, iset, "cya")
         x = 2.0*(x1 + x2 + x3)
         deallocate( x1 )
         deallocate( x2 )
         deallocate( x3 )
         return
      elseif  (property == "zoo") then      ! sum of micro and meso zooplankton
         allocate( x1(size(x,1), size(x,2), size(x,3)) )
         allocate( x2(size(x,1), size(x,2), size(x,3)) )  
         call get_eco_3D(x1, iset, "zo1")
         call get_eco_3D(x2, iset, "zo2")
         x = x1 + x2
         deallocate( x1 )
         deallocate( x2 )
         return
      endif
c     --- retrieve an ERGOM State variable
      slot = get_eco_slot(property) 
      if (slot < 1) then
         write(*,*) "get_eco_3D: unknown property", property
         stop
      endif
      do k=1,kmx(iset)
         do j=1,mmx(iset)
            do i=1,nmx(iset)
               m = m1(iset)%p(mmx(iset)+1-j,i,k)
               if(m > 0) x(i,j,k) = eco_i4(iset)%p(slot,m)/bscale   ! bscale applies to all eco vars
            enddo
         enddo
      enddo
      end subroutine get_eco_3D



      integer function io_new_unit ()
c----------------------------------------------------------------------------+
c     We have this one in run_context as well, but keep io_new_unit here
c     as local method now for the sake of code separation 
c     (try keeping read_cmod self contained and original)
c----------------------------------------------------------------------------+     
      integer            :: i
      integer, parameter :: io_min=7, io_max=130
      logical            :: file_is_open
c----------------------------------------------------------------------------+       
      file_is_open = .true.
      i = io_min-1
      do while (file_is_open)
         i = i+1
         if (i > io_max) then
            stop 'ERROR: Could not find new I/O unit number'
         endif
         inquire(unit=i,opened=file_is_open)
      enddo

      io_new_unit = i
      end function io_new_unit


c     ============= diagnostic tools =============

      subroutine check_buffer_bounds(iset)   
c     ------------------------------------------------------------------------
c     Assess min/max bounds of data in set iset in current integer rep buffers
c     by applying get_X(). Buffers must have been loaded
c     ------------------------------------------------------------------------
      integer, intent(in) :: iset
      real, allocatable   :: buf3d(:,:,:), buf2D(:,:)
      integer             :: isl
      character(len=3)    :: label
c     --------------------------------------------------------------------------
      allocate( buf2D(nmx(iset), mmx(iset))            )
      allocate( buf3D(nmx(iset), mmx(iset), kmx(iset)) )
c      
      write(*,*) "checking min/max vakues in data buffers:"
      call get_z(buf2D, iset)
      write(*,*) minval(buf2D), "< z < ", maxval(buf2D)
      call get_u(buf3D, iset)
      write(*,*) minval(buf3D), "< u < ", maxval(buf3D)
      call get_v(buf3D, iset)
      write(*,*) minval(buf3D), "< v < ", maxval(buf3D)
      call get_w(buf3D, iset)
      write(*,*) minval(buf3D), "< w < ", maxval(buf3D)
      call get_s(buf3D, iset)
      write(*,*) minval(buf3D), "< s < ", maxval(buf3D)
      call get_t(buf3D, iset)
      write(*,*) minval(buf3D), "< t < ", maxval(buf3D)
      call get_wu(buf2D, iset)
      write(*,*) minval(buf2D), "< wu < ", maxval(buf2D)
      call get_wv(buf2D, iset)
      write(*,*) minval(buf2D), "< wv < ", maxval(buf2D)
      do isl=1,neco
         label = eco_labels(isl)
         call get_eco_3D(buf3D, iset, label)
         write(*,*) minval(buf3D), " < ",label, " < ", maxval(buf3D)
      enddo
c
      deallocate (buf2D)
      deallocate (buf3D)
      end subroutine check_buffer_bounds
c
      end   ! module

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c$$$      program test_read_cmod 
c$$$      use read_cmod_ergom
c$$$      use time_tools
c$$$      implicit none
c$$$      integer :: gh,nx,ny,nz,lun=11, ix,iy,iz,idx1d,j
c$$$      integer :: request_time(4), isec
c$$$      integer :: year, month, day, hour
c$$$      real:: lambda1, dlambda, phi1, dphi,x,y
c$$$      character(len=256) :: fname
c$$$      type(clock)  :: now 
c$$$      character(len=11) :: tag
c$$$      logical       :: newdata
c$$$      real, allocatable :: u(:,:,:)
c$$$      
c$$$
c$$$      call init_read_cmod_ergom("../data_sample")
c$$$     
c$$$      call get_grid_descriptors("coarse",gh,lambda1,dlambda, ! fine/coarse
c$$$     +                          phi1,dphi, nx,ny,nz)
c$$$      write(*,*) "using:",gh,lambda1,dlambda, phi1,dphi,nx,ny,nz
c$$$
c$$$
c$$$c     ---- check file name mapping (need to make read_cmod_ergom: create_hydro_file_tag public ----    
c$$$
c$$$      do isec = 1, 86400, 100
c$$$         call set_clock(now, 2003, 4, 5, isec)
c$$$         call create_hydro_file_tag(now, tag, year, month, day, hour)
c$$$         write(*,*) isec, hour, tag
c$$$      enddo

c$$$          
c$$$c     ---- dump wet pixels: m1(iarea)%p(jy,jx,jz)>0 ----
c$$$
c$$$      do ix=1,nx
c$$$         do iy=1,ny
c$$$            idx1d = m1(gh)%p(ny+1-iy, ix, 1)
c$$$            x = lambda1 + (ix-1)*dlambda
c$$$            y = phi1    + (iy-1)*dphi
c$$$            if (idx1d>0) write(77,*) x,y ! wet map
c$$$         enddo
c$$$      enddo
c$$$      

c$$$      call set_clock(now, 2007, 1, 5, 22700)
c$$$      newdata = update_buffers(now)
c$$$      write(*,*) "newdata = ", newdata
c$$$      newdata = update_buffers(now)    ! test caching by iterative load
c$$$      write(*,*) "newdata = ", newdata
c$$$      allocate( u(nx,ny,nz) )
c$$$      call get_u(u, gh)
c$$$
c$$$c      write(*,*) u(5:7,6:8,1)


c$$$
c$$$      do ix=1,nx
c$$$         do iy=1,ny
c$$$            x = lambda1 + (ix-1)*dlambda
c$$$            y = phi1    + (iy-1)*dphi
c$$$            write(89,*) u(ix,iy,1)
c$$$         enddo
c$$$      enddo
c$$$

c$$$c     ---- vertical scan ----
c$$$      do ix=1,nx
c$$$         do iy=1,ny
c$$$            idx1d = m1(gh)%p(ny+1-iy, ix, 1)
c$$$            if (idx1d>0) then
c$$$               do iz=1,nz
c$$$                  j = m1(gh)%p(ny+1-iy, ix, iz)
c$$$                  if (j>0) write(*,*) hz(gh)%p(j), cellw(gh)%p(j)
c$$$               enddo
c$$$               stop 65
c$$$            endif        
c$$$         enddo
c$$$      enddo
c$$$
c$$$      call close_read_cmod_ergom()
c$$$      end program test_read_cmod 
c$$$     
c$$$      program test_file_mapping
c$$$      use read_cmod_ergom
c$$$      use time_tools
c$$$      implicit none
c$$$      type(clock)       :: now 
c$$$      character(len=11) :: tag
c$$$      integer           :: isec,year, month, day, hour
c$$$      do isec = 1, 86400, 100
c$$$         call set_clock(now, 2003, 4, 5, isec)
c$$$         call create_hydro_file_tag(now, tag, year, month, day, hour)
c$$$         write(*,566) isec/3600., tag, year, month, day, hour
c$$$      enddo
c$$$ 566  format(f8.4,1x,a11,3x,4i5)
c$$$      end program
c$$$
c$$$      call set_clock(now, 1995, 3, 20, int(16.1*3600))
c$$$      was_updated = update_buffers(now)
c$$$      write(*,*) "buffer update", was_updated
c$$$      do iset = 1, 2
c$$$         write(*,*) "checking buffer bounds for iset = ", iset
c$$$         call check_buffer_bounds(iset)
c$$$      enddo
c$$$      end program 
c$$$
c$$$      program test_read_cmod_minimal 
c$$$      use read_cmod_ergom
c$$$      use time_tools
c$$$      implicit none
c$$$      type(clock)  :: now 
c$$$      logical      :: was_updated
c$$$      integer      :: iset
c$$$      call init_read_cmod_ergom("~/DTU/OPEC/DMIdata/data_sample2")
c$$$      call set_clock(now, 1995, 3, 20, int(16.1*3600))
c$$$      was_updated = update_buffers(now)
c$$$      write(*,*) "buffer update", was_updated
c$$$      do iset = 1, 2
c$$$         write(*,*) "checking buffer bounds for iset = ", iset
c$$$         call check_buffer_bounds(iset)
c$$$      enddo
c$$$      end program 
