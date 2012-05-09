      module read_cmod
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Data Reading module for SUNFISH long physics time series based
c     on readtmp.f90 provided summer 2011 by Mikhail Dobrynin at DMI
c
c     Main public arrays: area refers to subgrid number 
c       m1(area)%p(:,:,:):    3D to compressed 1D map of wet grid points (cmod indexing)
c       hz(area)%p(:)    :    cell center depth 
c       cellw(area)%p(:) :    cell width
c
c       data buffers: (loaded by read_frame)
c         u(area)%p(:)     :    u current
c         v(area)%p(:)     :    v current
c         w(area)%p(:)     :    w current
c         z(area)%p(:)     :    sea surf elev
c         s(area)%p(:)     :    salinity
c         t(area)%p(:)     :    temperature
c         wu(area)%p(:)    :    wind u current
c         wv(area)%p(:)    :    wind v current
c
c     changes: 
c        free format -> fixed format
c        removed station+interpolation logic
c        removed exitme
c        changed unit cm/s -> m/s for currents
c        added array cellw giving cell widths in meters
c
c     TODO: 
c             
c        iso_points.tex ??
c
c     Compile:
c        ifort -convert big_endian -c read_cmod.f
c        ifort -convert big_endian -o test_read_cmod  read_cmod.f 
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
c        hz(iarea)%p(m1(iarea)%p(jy,jx,jz)) is cell center depth (AFTER hzmod)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     directives ----------------------------------------------------------------
      implicit none
      private  ! default scope

c     local type defs -----------------------------------------------------------------
      type cmi1
      integer(4), pointer :: p(:)
      end type
      type cmi3
      integer(4), pointer :: p(:,:,:)
      end type
      type cmr1
      real(8),    pointer :: p(:)
      end type

c     limits and other parameters -----------------------------------------------
      integer(4), parameter :: maxareas = 99 !max number of ocean model domains
      integer(4), parameter :: mainarea = 1 !number of the main area
      integer(4)            :: narea    = 2 !number of ocean model domains
      integer(4)            :: nfjord   = 0 !number of fjord model domains   
      real(8),    parameter :: undef    = 999.9 ! undefined
      real(8),    parameter :: nundef   = 999.0 ! nearly undefined (hmmm!)
      logical,    parameter :: ldebug   = .false.  ! 

c     allocatable arrays --------------------------------------------------------
      integer(4),     allocatable, dimension(:)   :: mmx,nmx,kmx,iw2,iw3 ! dimensions
      type (cmi3),    pointer,     dimension(:)   :: m1
      type (cmr1),    pointer,     dimension(:)   :: hz,cellw
      type (cmr1),    allocatable, dimension(:)   :: u,v,w,z,s,t,wu,wv    ! real rep buffers
      type (cmi1),    allocatable, dimension(:)   :: u_i4,v_i4,w_i4,z_i4  ! integer rep buffers
      type (cmi1),    allocatable, dimension(:)   :: s_i4,t_i4            ! integer rep buffers
      type (cmi1),    allocatable, dimension(:)   :: wu_i4,wv_i4          ! integer rep buffers
      real(8),        allocatable, dimension(:)   :: xs,ys,y0,x0,dxa,dya
      real(8),        allocatable, dimension(:)   :: y9,x9    

c     -----------------------  data handlers  ------------------------------------

      character(len=999) :: current_file    = " "   ! name of file containing current frame
      integer            :: current_frame(4) = -1   ! current frame, as year, month, day, hour in day
      integer            :: iunit            = -1   ! assign at first read 

c     Config and bathy file names -----------------------------------------------
      
      integer(4)                      :: timelevel_in, maxtimelevels
      integer(4)                      :: nestinglevels_in 
      integer(4)                      :: nestingto_in(maxareas)
      integer(4), allocatable         :: nestinglevels(:), timelevel(:)
      type (cmi1), pointer            :: nestingto(:), nestingfrom(:)
      logical, allocatable            :: onthislevel(:,:)
      character(LEN=256), allocatable :: cfgfile(:), bathyfile(:)

c     -----------------------  set module public scope -----------------------

      public :: cmi1, cmi3, cmr1
      public :: init_read_cmod
      public :: close_read_cmod
      public :: get_grid_descriptors  
      public :: update_buffers
      public :: m1,hz,cellw
      public :: u,v,w,z,s,t,wu,wv     ! do not make *_i4 buffers public
      public :: undef
c     ========================================================================
                                  contains 
c     ========================================================================


      subroutine init_read_cmod(path)
c----------------------------------------------------------------------------+
c     path is to grid descriptor files, without trailing slash
c----------------------------------------------------------------------------+ 
      character*(*), intent(in) :: path ! applied as prefix to file names (without trailing slash)

      character(LEN=256)   :: cfgfile_in, bathyfile_in
      integer(4)           :: kmax,nlevs,ios
      integer(4)           :: ii,iia,iii,iv, areafound
      integer(4)           :: lun,ia
      character(LEN=999)   :: prefix           
c     namelists 
      namelist /cfglist/ narea, nfjord
      namelist /cfgfn/   cfgfile_in, bathyfile_in, timelevel_in,    
     +                   nestinglevels_in, nestingto_in
c----------------------------------------------------------------------------+ 

      prefix = trim(adjustl(path)) // "/" ! apply linux delimiter
      write(*,*) "init_read_cmod: path = ", trim(prefix)
c     Get the ocean model configuration --------------------------------------
      lun = io_new_unit()
      open (unit=lun, file=trim(prefix)//'cfg.nml', 
     +      status='old', iostat=ios)
      if (ios /= 0) then
         write(*,*) "init_read_cmod: config file not found"
         write(*,*) "file name = ",trim(prefix)//'cfg.nml'
         stop  ! removed default setting, require explicit input
      endif


c     apply settings from cfg file:
      read (unit=lun, nml=cfglist, iostat=ios)
      if (ios /= 0) then
         stop 'Error: Cannot read basic numbers from cfg file'
      endif

      if (narea <= 0 .or. narea > maxareas) then
         stop 'Error: Invalid number of areas'
      endif
      
      if (nfjord > 0) then
         stop 'Error: Cannot handle fjords'
      endif

c     allocate and get file names, etc:
      allocate( cfgfile(narea), bathyfile(narea), timelevel(narea),         
     +     nestinglevels(narea), nestingto(narea),                         
     +     mmx(narea),nmx(narea),kmx(narea),iw2(narea),iw3(narea),          
     +     m1(narea),hz(narea),cellw(narea))

      do ia=1,narea
         cfgfile_in   = ' '
         bathyfile_in = ' '
         timelevel_in = 1
         if (ia == narea) then
            nestinglevels_in = 0
            nestingto_in(:)  = 0
         else
            nestinglevels_in = 1
            nestingto_in(:)  = min( ia+1, narea )
         endif
         read (unit=lun, nml=cfgfn, iostat=ios)
         if ( ios /= 0 ) then
            stop 'Error: Cannot read file names from cfg file'
         else
            cfgfile(ia)   = trim(cfgfile_in)
            bathyfile(ia) = trim(bathyfile_in)
            timelevel(ia) = timelevel_in
            nlevs         = nestinglevels_in
            if (nlevs < 0 .or. nlevs > maxareas) then
               write(*,*) 'Error: Invalid No. of levels for area ',ia
               stop 
            endif
            nestinglevels(ia) = nlevs

c
c           ifort compiler notice:associated(pt) == T even if pt has not been allocated    
c        
            if (nlevs > 0) then
               allocate( nestingto(ia)%p(nlevs) )
               nestingto(ia)%p(1:nlevs) = nestingto_in(1:nlevs)
            else
               nullify ( nestingto(ia)%p ) ! render nestingto(ia)%p in well defined state
            endif


            if (ia == mainarea .and. timelevel(ia) /= 1) then
               stop 'Error: Time level of main area must be 1'
            elseif (timelevel(ia) < 1) then
               stop 'Error: Time level must be a positive integer'
            endif

            if (ia == mainarea .and. nestinglevels(ia) < 1 
     +                         .and. narea > 1) then
               write(*,*) 'Warning: No nesting into main area'//
     +                    '... continue without'
            endif
         endif
      enddo
      close(unit=lun)
      
c     validate configuration:
      do ia=1,narea
         do nlevs=1,nestinglevels(ia)
            if (nestingto(ia)%p(nlevs) == ia    .or.   
     +          nestingto(ia)%p(nlevs) <= 0     .or.          
     +          nestingto(ia)%p(nlevs) > narea ) then
               write(*,*) 'Error: Invalid nesting for area ',ia
               stop 
            endif

            if (timelevel(ia) > timelevel(nestingto(ia)%p(nlevs))) then
              stop 'Error: Cannot nest to domains on lower time levels' 
            endif
         enddo
      enddo


c     Do the remaining nesting set up -------------------------------------------
      maxtimelevels = maxval(timelevel(:))
      allocate( onthislevel(narea,maxtimelevels), nestingfrom(narea) )
      allocate( nestingfrom(mainarea)%p(0:0) )
      nestingfrom(mainarea)%p(0) = 0
      do ia=2,narea
         ii = 0
         do iia=1,narea
            if (iia /= ia .and. nestinglevels(iia) > 0) then
               do iii=1,nestinglevels(iia)
                  if (nestingto(iia)%p(iii) == ia) ii = ii + 1
               enddo
            endif
         enddo
         allocate( nestingfrom(ia)%p(0:ii) )
         nestingfrom(ia)%p(0) = ii
         if (ii > 0) then
            iv = 0
            do iia=1,narea
               if (iia /= ia .and. nestinglevels(iia) > 0) then
                  do iii=1,nestinglevels(iia)
                     if (nestingto(iia)%p(iii) == ia) then
                        iv = iv + 1
                        nestingfrom(ia)%p(iv) = iia
                     endif
                  enddo
               endif
            enddo
            if (iv /= ii) then
               stop 'Error: Cannot figure out how to do the nesting'
            endif
         endif
      enddo
      onthislevel(:,:) = .false.
      do ia=1,narea
         onthislevel(ia,timelevel(ia)) = .true.
      enddo
c  
c     GET MODEL SIZE, WET POINT INDEX, DEPTHS AND ALLOCATE ARRAYS 
c
      do ia=1,narea
         call read_dima(trim(prefix)//cfgfile(ia),
     +                  mmx(ia),nmx(ia),kmx(ia))
         call read_wet(trim(prefix)//bathyfile(ia),
     +        iw2(ia),iw3(ia),mmx(ia),nmx(ia), 
     +        kmx(ia),m1(ia)%p,hz(ia)%p)
         allocate( cellw(ia)%p(0:iw3(ia)) ) ! hz is allocated in read_wet
         cellw(ia)%p = hz(ia)%p ! copy data before hzmod transforms hz
         call hzmod(m1(ia)%p,hz(ia)%p,iw3(ia),mmx(ia),nmx(ia),kmx(ia))
      enddo

      kmax = maxval(kmx)
      
      allocate( u(narea), v(narea), w(narea), z(narea), s(narea), 
     +     t(narea), wu(narea), wv(narea), u_i4(narea), v_i4(narea), 
     +     w_i4(narea), z_i4(narea), s_i4(narea), t_i4(narea), 
     +     wu_i4(narea), wv_i4(narea), x0(narea), x9(narea), 
     +     y0(narea), y9(narea), dxa(narea), dya(narea))

      do ia=1,narea
         allocate( u(ia)%p(0:iw3(ia)), v(ia)%p(0:iw3(ia)), 
     +        z(ia)%p(0:iw2(ia)),s(ia)%p(0:iw3(ia)), t(ia)%p(0:iw3(ia)),
     +        w(ia)%p(0:iw3(ia)), wu(ia)%p(0:iw3(ia)), 
     +        wv(ia)%p(0:iw3(ia)),u_i4(ia)%p(0:iw3(ia)), 
     +        v_i4(ia)%p(0:iw3(ia)),w_i4(ia)%p(0:iw3(ia)), 
     +        z_i4(ia)%p(0:iw2(ia)),s_i4(ia)%p(0:iw3(ia)), 
     +        t_i4(ia)%p(0:iw3(ia)), wu_i4(ia)%p(0:iw2(ia)), 
     +        wv_i4(ia)%p(0:iw2(ia)) )
      enddo

c     ---------------------------- INITIALIZE GRID -------------------------------
      do ia=1,narea
         call readgridmin(trim(prefix)//cfgfile(ia),
     +                    x0(ia),y0(ia),dxa(ia),dya(ia)) ! returns degrees,
         y9(ia) = y0(ia) - (mmx(ia)-1)*dya(ia)
         x9(ia) = x0(ia) + (nmx(ia)-1)*dxa(ia)
         if (ldebug) then
            write(*,*) 'grid number:',ia
            write(*,*) '  upper left =',x0(ia),y0(ia)
            write(*,*) '  lower right=',x9(ia),y9(ia)            
         endif
c 
c        -------- print grid info --------
c        
         write(*,*) "init_read_cmod: grid number: ", ia
         write(*,*) "init_read_cmod: grid ID    : ", trim(cfgfile(ia))
         write(*,*) "long range = [", x0(ia),x9(ia),"] ", "deg E"
         write(*,*) "lati range = [", y0(ia),y9(ia),"] ", "deg N"
         write(*,*) "dlon,dlat  = ", dxa(ia), dya(ia),    "deg"
         write(*,*) "grid dims  = ", nmx(ia), mmx(ia), kmx(ia)
         write(*,*) "surf/bulk wet points = ", iw2(ia),iw3(ia)
         
      enddo

      end subroutine init_read_cmod




      subroutine close_read_cmod()
c     -------------------------------------------------------------------
c     deallocate module data arrays/pointers
c     -------------------------------------------------------------------  
      integer(4)           :: ia
      deallocate( mmx,nmx,kmx,iw2,iw3 )
      deallocate( nestinglevels, timelevel )
      deallocate( onthislevel )
      deallocate( cfgfile, bathyfile )
c     -------------------------------------------------------------------  
      do ia=1,narea          
         if (associated( m1(ia)%p )) deallocate( m1(ia)%p )
         if (associated( hz(ia)%p )) deallocate( hz(ia)%p )
         if (associated( cellw(ia)%p )) deallocate( cellw(ia)%p )
         if (associated( u(ia)%p )) deallocate( u(ia)%p )
         if (associated( v(ia)%p )) deallocate( v(ia)%p )
         if (associated( w(ia)%p )) deallocate( w(ia)%p )
         if (associated( z(ia)%p )) deallocate( z(ia)%p )
         if (associated( s(ia)%p )) deallocate( s(ia)%p )
         if (associated( t(ia)%p )) deallocate( t(ia)%p )
         if (associated( wu(ia)%p )) deallocate( wu(ia)%p )
         if (associated( wv(ia)%p )) deallocate( wv(ia)%p )
         if (associated( u_i4(ia)%p )) deallocate( u_i4(ia)%p )
         if (associated( v_i4(ia)%p )) deallocate( v_i4(ia)%p )
         if (associated( w_i4(ia)%p )) deallocate( w_i4(ia)%p )
         if (associated( z_i4(ia)%p )) deallocate( z_i4(ia)%p )
         if (associated( s_i4(ia)%p )) deallocate( s_i4(ia)%p )
         if (associated( t_i4(ia)%p )) deallocate( t_i4(ia)%p )
         if (associated( wu_i4(ia)%p )) deallocate( wu_i4(ia)%p )
         if (associated( wv_i4(ia)%p )) deallocate( wv_i4(ia)%p )
         if (associated( nestingto(ia)%p )) deallocate( nestingto(ia)%p)
         if (associated( nestingfrom(ia)%p )) 
     +                   deallocate( nestingfrom(ia)%p )
      enddo
      
      deallocate( m1 )
      deallocate( hz,cellw )
      deallocate( u,v,w,z,s,t,wu,wv )   
      deallocate( u_i4,v_i4,w_i4,z_i4 )  
      deallocate( s_i4,t_i4 )            
      deallocate( wu_i4,wv_i4 )          
      deallocate( nestingto, nestingfrom )

c     --- reset data handlers ---
      current_file  = " "  
      current_frame = -1   
      iunit         = -1   

      end subroutine close_read_cmod




      subroutine get_grid_descriptors(data_set_id, data_set_handler,
     +                          lambda1,dlambda, phi1,dphi, nx,ny,nz) 
c     ----------------------------------------------------------------------------+
c     Interface query function to retrieve grid characteristics:
c
c       data_set_handler           : data set index, corrsponding to  data_set_id
c       lambda1,dlambda, phi1,dphi : lower left grid point
c       nx,ny,nz                   : grid dimensions
c
c     ----------------------------------------------------------------------------+
      character*(*),intent(in) :: data_set_id      ! map to cfgfile_in, currently "data_coarse"/"data_fine"
      integer,intent(out)      :: data_set_handler ! area number for this set
      real,intent(out)         :: lambda1, dlambda, phi1, dphi
      integer,intent(out)      :: nx,ny,nz
      integer                  :: ia
      character(len=256)       :: target_name
c     ----------------------------------------------------------------------------+     
      write(target_name,974) trim(adjustl(data_set_id))
 974  format("data_",a) 

      do ia = 1,narea
         if (trim(adjustl(target_name)) == 
     +       trim(adjustl(cfgfile(ia)))) then
            exit   ! do loop
         endif
      enddo

      if (ia>narea) then  ! standard condition when loop terminates without exit
         write(*,*) "get_grid_descriptors: no match for ",
     +               trim(adjustl(data_set_id))
         stop
      else
         data_set_handler = ia
         lambda1          = x0(ia)
         dlambda          = dxa(ia)
         phi1             = y9(ia)
         dphi             = dya(ia)
         nx               = nmx(ia)
         ny               = mmx(ia)
         nz               = kmx(ia)
      endif
      
      end subroutine get_grid_descriptors
                                  


      subroutine update_buffers(fname,year,month,day,hour)
c     ----------------------------------------------------------------------------+
c     Update buffers corresponding to requested time (year,month,day,hour)
c     with data from tempfile with name fname
c     Update is intelligent in the sense that file reading is only done, if needed
c     Assume needed, if fname has changed (crude assumption)
c     ----------------------------------------------------------------------------+
      character*(*), intent(in)  :: fname
      integer, intent(in)        :: year,month,day,hour
      integer                    :: request_time(4), tcur, treq, ios
      logical                    :: newfile, newtime
c     ----------------------------------------------------------------------------+
      request_time(1) = year
      request_time(2) = month
      request_time(3) = day
      request_time(4) = hour
      tcur = time_index( current_frame )
      treq = time_index( request_time  )
      newfile   = .false.
      if (trim(adjustl(fname)) /= trim(adjustl(current_file))) then
         newfile   = .true. 
         if (iunit>0) close(iunit)                               ! we opened it previously
         if (iunit<0) iunit = io_new_unit()                      ! first time, pick a free unit 
         open(iunit,file=fname, form='unformatted',status='old',
     +               iostat=ios )                                ! keep same logical unit
         if (ios /= 0) then
            write(*,'(a,a)') "error opening data file ", 
     +                       trim(adjustl(fname))
            stop
         else
            write(*,'(a,a)') "update_buffers: opened new data file", 
     +                       trim(adjustl(fname))
         endif
         call reset_tempfile(iunit)
         current_file = fname
      endif
c      
      tcur = time_index( current_frame )
      treq = time_index( request_time  )
      newtime = (tcur /= treq)
c
c     decide to read or not
c
      if (newtime.or.newfile) then  
         if (treq <= tcur) call reset_tempfile(iunit) ! will rewind if same frame is requested  
         call read_frame(iunit, request_time)
      endif

      end subroutine update_buffers



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


      

      subroutine read_frame(lun, request_time)
c     ----------------------------------------------------------------------------+
c     --------------------------------- READ DATA ---------------------------------
c     ----------------------------------------------------------------------------+
c
c     Verbose load of time frame request_time from opened tempfile 
c     with file handler lun. read_frame will continue reading from current reading position
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
c     file content:
c
c         aufdat             
c         ctim1,tim1
c         data_segment1 (N=2 areas)
c         ctim2,tim2
c         data_segment2 (N=2 areas)
c         ...
c     it is assumed that data is time ordered: tim1 < tim2 < ...
c     ----------------------------------------------------------------------------+ 
      integer,intent(in)       :: lun
      integer,intent(in)       :: request_time(4) ! year, month, day, hour in day
 
      integer(4)               :: ia,it,ios
      character(len=11)        :: aufdat
      character(len=19)        :: ctim
      real(8)                  :: tim
      integer                  :: treq, tcur
      integer                  :: time_stamp(4)
c     ----------------------------------------------------------------------------+      

c     read time stamp from tempdat file:
      read(lun,iostat=ios) ctim,tim
      it = ios
      do while (it == 0)
c
         read(ctim,532) time_stamp
 532     format(i4,1x,i2,1x,i2,1x,i2)
         if (ldebug) write(*,*) "reading : ", ctim
c        read vars from tempdat file:
         read(lun) (wu_i4(ia)%p, wv_i4(ia)%p, ia=1,narea)
         do ia=1,narea
            read(lun) u_i4(ia)%p, v_i4(ia)%p, w_i4(ia)%p, z_i4(ia)%p
         enddo
         do ia=1,narea
            read(lun) s_i4(ia)%p, t_i4(ia)%p
         enddo

         if ( is_equal(time_stamp, request_time) ) then
c        
c           invoke transformation including surface fields (asc.09Sep2011)
c
            do ia=1,narea
               call wcomptr( iw2(ia),iw3(ia),u(ia)%p,v(ia)%p,w(ia)%p,
     +              z(ia)%p,u_i4(ia)%p,v_i4(ia)%p,w_i4(ia)%p,z_i4(ia)%p, 
     +              s(ia)%p,t(ia)%p,s_i4(ia)%p,t_i4(ia)%p,wu(ia)%p,
     +              wv(ia)%p,wu_i4(ia)%p,wv_i4(ia)%p )
            enddo
            write(*,*) "read_frame: loaded frame ", ctim
            current_frame = time_stamp ! update data handler
            return                     ! we got what we looked for - exit
         endif                
         
c     in readtmp, here was the stations loop - removed all this (asc.09Sep2011)

c     read next time stamp from tempdat file:

         read(lun,iostat=ios) ctim,tim
         it = ios
      enddo                     ! it-loop (time frames in tempdat file)
      
c     capture unexpected case where requested data
c     is not in file, or data is otherwise corrupted
c     At successful read, we will exit in the do-while loop above     

      write(*,*) "read_frame: did not find requested time frame in file"
      write(*,*) "read_frame: requested time = ", request_time
      write(*,*) "read_frame: last record    = ", ctim
      stop ! fatal condition

      contains

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
      integer, parameter  :: year2h  = 366*40*24 ! multiplication factor must be >= 366
      integer, parameter  :: month2h = 31*24     ! multiplication factor must be >= 31
      integer, parameter  :: day2h   = 24
c     -----------------------------------------      
      time_index = (i4(1)-2000.0)*year2h + i4(2)*month2h + 
     +              i4(3)*day2h + i4(4)
      end function time_index


      subroutine hzmod(mm,hz,iw3,m,n,k)
c     ---------------------------------------------------------------
c     transform hz from cell width array to cell center depths
c     ---------------------------------------------------------------
      integer, intent(in)    :: iw3,m,n,k
      integer, intent(in)    :: mm(m,n,k)
      real(8), intent(inout) :: hz(0:iw3)
      
      integer :: i,j,l,l1,m1
      real(8) :: h1(1:iw3), hs
c     --------------------------------------------------------------   
      do i=1,m
         do j=1,n
            do l=1,k
               m1 = mm(i,j,l)
               if (m1 /= 0) then
                  hs = 0.0
                  do l1=1,l
                     hs = hs + hz(mm(i,j,l1))
                  enddo
                  h1(m1) = hs - hz(m1)*0.5
               endif
            enddo
         enddo
      enddo

      hz(1:) = h1(1:)
      
      end subroutine hzmod





      subroutine wcomptr(iwet2,iwet3,u,v,w,z,u_i4,v_i4,w_i4,z_i4,s,t,
     +                   s_i4,t_i4,wu,wv,wu_i4,wv_i4 )
c----------------------------------------------------------------------------+
c
c     compresses the calculated data
c
c     output units:
c      u,v,w   currents           = m/s 
c      t       temperature        = deg Celcius
c      s       salty              = PSU
c      z       sea sruf elevation = m    (not confirmed)
c      wu,wv   wind               = m/s
c      positive directions of currents,wind are east/north
c
c     NOTE: in readtmp.f90 current output unit is cm/s
c           divided uscal,wscal by factor 100 because IBMlib prefers SI units
c----------------------------------------------------------------------------+
      integer, intent(in)                 :: iwet2,iwet3
      real(8), intent(out), dimension(0:) :: u,v,w,s,t,z
      integer, intent(in),  dimension(0:) :: u_i4,v_i4,w_i4,s_i4,t_i4
      integer, intent(in),  dimension(0:) :: z_i4
      real(8), intent(out), dimension(0:), optional :: wu,wv
      integer, intent(in),  dimension(0:), optional :: wu_i4,wv_i4
      
      real(8), parameter :: uscal  = 1.d-8  ! previously 1.d-6  in readtmp.f90
      real(8), parameter :: wscal  = 1.d-13 ! previously 1.d-11 in readtmp.f90
      real(8), parameter :: tscal  = 1.d-7
      real(8), parameter :: zscal  = 1.d-8
      real(8), parameter :: wuscal = 1.d-7
      
      real(8), parameter :: scalarmax = 50.d0
      real(8), parameter :: levundef  = 99.d0

c     FIXME: Why is zero an undefined value for velocities???
c     Why is u,v,w initialized as undefined??? 
c     Pls, help comming up with a good answer.
      real(8), parameter :: velundef  =  0.d0 
      integer :: i
c----------------------------------------------------------------------------+ 
      
      u(0)  = 0.
      v(0)  = 0.
      w(0)  = 0.
      s(0)  = 0.
      t(0)  = 0.
c     FIXME: all four wu, wv, wu_i4, wv_i4 must be present; the sloppyness of
c     the present developer prevent him from validating this nicely.
      if (present(wu)) then 
         z(0)  = 0.
         wu(0) = 0.
         wv(0) = 0.
      else
         z(:) = 0.
      endif

      do i=0,iwet3
         u(i) = u_i4(i)*uscal
         v(i) = v_i4(i)*uscal
         w(i) = w_i4(i)*wscal
         s(i) = s_i4(i)*tscal
         t(i) = t_i4(i)*tscal
         if (abs(u(i)) == velundef)  u(i) = undef
         if (abs(v(i)) == velundef)  v(i) = undef
         if (abs(w(i)) == velundef)  w(i) = undef
         if (abs(s(i)) >  scalarmax) s(i) = undef
         if (abs(t(i)) >  scalarmax) t(i) = undef
      enddo

      if (present(wu)) then 
         do i=0,iwet2
            z(i)  = z_i4(i)*zscal
            wu(i) = wu_i4(i)*wuscal
            wv(i) = wv_i4(i)*wuscal
            if (abs(z(i))  >= levundef)  z(i) = undef
            if (abs(wu(i)) >= levundef) wu(i) = undef
            if (abs(wv(i)) >= levundef) wv(i) = undef
         enddo
      endif

      end subroutine wcomptr



      subroutine read_dima(filename,m,n,k)
c----------------------------------------------------------------------------+
c----------------------------------------------------------------------------+      
      character(len = *), intent(in)  :: filename
      integer,            intent(out) :: m,n,k   
      integer       :: ios, lun
      real(8)       :: t(300)
      character(3)  :: ch3
      character(80) :: line
c----------------------------------------------------------------------------+ 
      
      lun = io_new_unit()
      open(lun,file=trim(filename),status='old',iostat=ios)
      if (ios /= 0) then
         write(*,*) 'Error opening '//filename
         stop
      endif   

c     Read parameters from data file ------------------------------------------
      read(lun,'(2X,A3)') ch3
      do while (ch3 /= 'DIM')
         read(lun,'(2X,A3)') ch3
      enddo
      read(lun,'(14X,I4)') n
      read(lun,'(14X,I4)') m
      do while (ch3 /= 'SCH')
         read(lun,'(2X,A3)') ch3
      enddo
      k=0
      read(lun,'(A)') line
      ch3 = line(6:8)
      do while (ch3 == 'SCH')
         k=k+1
         read(line(17:26),'(F10.5)') t(k)
         read(lun,'(A)') line
         ch3 = line(6:8)
      enddo
      
      close(lun)
      
      end subroutine read_dima




      subroutine read_wet(fname,iw2,iw3,mmx,nmx,kmx,m,h)
c----------------------------------------------------------------------------+
c----------------------------------------------------------------------------+  
      integer,      intent(in)  :: mmx,nmx,kmx
      integer,      intent(out) :: iw2,iw3
      character(*), intent(in)  :: fname
      
      integer, pointer :: m(:,:,:)
      real(8), pointer :: h(:)    ! cell widths
   
      integer :: lun
      integer :: ios
c----------------------------------------------------------------------------+ 
c     alloc index array:
      allocate( m(mmx,nmx,kmx) )

c     read index array from file:
      lun = io_new_unit()
      open(lun,file=trim(fname),iostat=ios)
      if (ios /= 0) then
         write(*,*) 'Error opening '//trim(fname)
         stop
      endif

      read(lun,*) m

c     obtain No of 3D wet points and surface wet points:
      iw3 = maxval(m)
      iw2 = maxval(m(:,:,1))     ! scanned layer wise, starting at surface layer

c     allocate and read depth array:
      allocate( h(0:iw3) )
      read(lun,'(15f8.2)') h
      close(lun)

      end subroutine read_wet




      subroutine readgridmin(filename,xor,yor,dx,dy)
c----------------------------------------------------------------------------+       
c----------------------------------------------------------------------------+
      character(*), intent(in)  :: filename
      real(8),      intent(out) :: xor,yor,dy,dx
      
      integer       :: deg,minu,sec, ios
      character(80) :: string
      integer       :: lun
c----------------------------------------------------------------------------+ 

      lun = io_new_unit()
      open(lun,file=trim(filename),status='old',iostat=ios)
      if (ios /= 0) then
         write(*,*) 'Could not open data file'//trim(filename)
         stop
      endif

      read(lun,'(a80)',iostat=ios) string
      if (ios /= 0) then
         write(*,*) 'Could not read data file'//trim(filename)
         stop
      endif
      do while (string(3:5) /= 'POS')
         read(lun,'(a80)',iostat=ios) string
         if (ios /= 0) then
            write(*,*) 'Could not read data file'//trim(filename)
            stop
      endif
      enddo
      read(lun,'(a80)') string
      read(lun,'(a80)') string
      read(string(18:20),'(i3)') deg
      read(string(22:23),'(i2)') minu
      read(string(25:26),'(i2)') sec
      yor = (deg*60 + minu + sec/60.)/60.
      if (string(27:27) == 'S') yor = -yor
      read(lun,'(a80)') string
      read(string(18:20),'(i3)') deg
      read(string(22:23),'(i2)') minu
      read(string(25:26),'(i2)') sec
      xor = (deg*60 + minu + sec/60.)/60.
      if (string(27:27) == 'W') xor = -xor
      
      read(lun,'(a80)',iostat=ios) string
      if (ios /= 0) then
         write(*,*) 'Could not read data file'//trim(filename)
         stop
      endif
         
      do while (string(3:5) /= 'GIT')
         read(lun,'(a80)',iostat=ios) string
         if (ios /= 0) then
            write(*,*) 'Could not read data file'//trim(filename)
            stop
         endif
      enddo
      read(lun,'(a80)') string
      read(string(18:20),'(i3)') deg
      read(string(22:23),'(i2)') minu
      read(string(25:26),'(i2)') sec
      dy = (deg*60 + minu + sec/60.)/60.
      read(lun,'(a80)') string
      read(string(18:20),'(i3)') deg
      read(string(22:23),'(i2)') minu
      read(string(25:26),'(i2)') sec
      dx = (deg*60 + minu + sec/60.)/60.
      
      close(lun)
      
      end subroutine readgridmin


      
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

      end   ! module

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c$$$      program test_read_cmod 
c$$$      use read_cmod
c$$$      implicit none
c$$$      integer :: gh,nx,ny,nz,lun=11, ix,iy,iz,idx1d,j
c$$$      integer :: request_time(4)
c$$$      real:: lambda1, dlambda, phi1, dphi,x,y
c$$$      character(len=256) :: fname
c$$$      
c$$$      
c$$$      call init_read_cmod(".")
c$$$      call get_grid_descriptors("data_fine",gh,lambda1,dlambda, ! fine/coarse
c$$$     +                          phi1,dphi, nx,ny,nz)
c$$$      write(*,*) gh,lambda1,dlambda, phi1,dphi, nx,ny,nz
c$$$c
c$$$cc     ---- wet pixels: m1(iarea)%p(jy,jx,jz)>0
c$$$c      do ix=1,nx
c$$$c         do iy=1,ny
c$$$c            idx1d = m1(gh)%p(ny+1-iy, ix, 1)
c$$$c            x = lambda1 + (ix-1)*dlambda
c$$$c            y = phi1    + (iy-1)*dphi
c$$$c            if (idx1d>0) write(77,*) x,y ! wet map
c$$$c         enddo
c$$$c      enddo
c$$$      
c$$$      fname="~/DTU/SUNFISH/data_long_phys/test_data/archive.2001031600"
c$$$      call update_buffers(fname,2001, 3, 15, 20)
c$$$      call update_buffers(fname,2001, 3, 15, 21)
c$$$
c$$$c$$$      do ix=1,nx
c$$$c$$$         do iy=1,ny
c$$$c$$$            idx1d = m1(gh)%p(ny+1-iy, ix, 1)
c$$$c$$$            if (idx1d>0) then
c$$$c$$$               x = lambda1 + (ix-1)*dlambda
c$$$c$$$               y = phi1    + (iy-1)*dphi
c$$$c$$$               if (u(gh)%p(idx1d)<900) write(88,*) x,y
c$$$c$$$               if (v(gh)%p(idx1d)<900) write(89,*) x,y
c$$$c$$$c               write(88,*) idx1d,u(gh)%p(idx1d),v(gh)%p(idx1d)
c$$$c$$$c               write(89,*) x,y,w(gh)%p(idx1d),z(gh)%p(idx1d)
c$$$c$$$c               write(90,*) x,y,s(gh)%p(idx1d),t(gh)%p(idx1d)
c$$$c$$$c               write(91,*) x,y,wu(gh)%p(idx1d),wv(gh)%p(idx1d)
c$$$c$$$            endif
c$$$c$$$         enddo
c$$$c$$$      enddo
c$$$c$$$      stop
c$$$
c$$$      do ix=1,nx
c$$$         do iy=1,ny
c$$$            idx1d = m1(gh)%p(ny+1-iy, ix, 1)
c$$$            if (idx1d>0) then
c$$$               do iz=1,nz
c$$$                  j = m1(gh)%p(ny+1-iy, ix, iz)
c$$$                  if (j>0) write(*,*) hz(gh)%p(j)
c$$$               enddo
c$$$               stop 65
c$$$            endif        
c$$$         enddo
c$$$      enddo
c$$$c$$$               x = lambda1 + (ix-1)*dlambda
c$$$c$$$               y = phi1    + (iy-1)*dphi
c$$$c$$$               if (u(gh)%p(idx1d)<900) write(88,*) x,y
c$$$c$$$               if (v(gh)%p(idx1d)<900) write(89,*) x,y
c$$$c$$$c               write(88,*) idx1d,u(gh)%p(idx1d),v(gh)%p(idx1d)
c$$$c$$$c               write(89,*) x,y,w(gh)%p(idx1d),z(gh)%p(idx1d)
c$$$c$$$c               write(90,*) x,y,s(gh)%p(idx1d),t(gh)%p(idx1d)
c$$$c$$$c               write(91,*) x,y,wu(gh)%p(idx1d),wv(gh)%p(idx1d)
c$$$c$$$            endif
c$$$c$$$         enddo
c$$$c$$$      enddo
c$$$c$$$      stop
c$$$
c$$$      call close_read_cmod()
c$$$
c$$$      end program test_read_cmod 
      
