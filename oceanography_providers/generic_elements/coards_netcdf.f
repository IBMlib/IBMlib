ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     -----------------------------------------------------------------
c     Decorate nf90 netcdf with some the conventions of COARDS 
c     -----------------------------------------------------------------
c
c     NB: default scope = public (decause it augments netcdf)
c
c     * target netcdf variables stored as (x,y,z,t)
c     * maskerade if a float variable is stored as float or scaled integer      
c     * internalize netcdf handlers
c     * adapt COARDS layout variations to IBMlib layout for easy access
c     * extract wetmask(nx,ny) and bottom_layer(nx,ny)
c     * follow https://ferret.pmel.noaa.gov/Ferret/documentation/coards-netcdf-conventions
c     ifort -e90 -c -I/usr/local_intel/include coards_netcdf.f
c
c     Jun 2020: ifort executable crashes if an integer buffer is provided for NF90_SHORT data reading
c     currently apply this mapping (based on byte-based declaration):
c     Beware that although it is common for the KIND parameter to be the same as the number of bytes
c     stored in a variable of that KIND, it is not required by the Fortran standard.
c     -------------------------------------------------------------------      
c     netcdftype        netcdf handle   bits        map to intrinsic f90 type
c     -------------------------------------------------------------------          
c     byte      	NF90_BYTE 	8
c     char   		NF90_CHAR 	8
c     short 		NF90_SHORT 	16          INTEGER(kind=2)
c     int 		NF90_INT 	32          INTEGER(kind=4)
c     float 		NF90_FLOAT 	32          REAL (currently implicit cast)
c     double 		NF90_DOUBLE 	64          REAL (currently implicit cast)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module coards_netcdf

      use netcdf   ! decorate and reexport (scope = public)

      implicit none     

      type coards_xyzt_variable
      private                      
      integer            :: ncid          ! attached netcdf set ID 
      integer            :: varid         ! attached netcdf variable ID within set ncid
      character(len=999) :: varname       ! store reference for easier tracking
      integer            :: nx,ny,nz,nt   ! local copy of dimensions
      integer            :: storage_type  ! int_type / float_type
      real               :: scale, offset ! for integer storage, x_real = scale*x_int + offset
      end type

      integer, parameter :: i2type  = 1   ! NF90_SHORT mapping 
      integer, parameter :: i4type  = 2   ! NF90_INT mapping
      integer, parameter :: f4type  = 3   ! NF90_FLOAT mapping
      integer, parameter :: f8type  = 4   ! NF90_DOUBLE mapping
      character(len=6), parameter :: storage_type_names(4) =
     +    (/"i2type", "i4type", "f4type", "f8type"/)  ! for debugging
      
      interface load_xyz_frame
         module procedure load_xyz_frame_asreal
      end interface

c     =================================================================
                          contains
c     =================================================================

      
      subroutine resolve_time_parameters(ncid, time_vname, year, month,
     +                                   day, hour, minute, sec, tunit)
c     -----------------------------------------------------------------
c     Resolve time offset as and time unit in data set from time:unit attribute
c
c     1) Resolve time offset as (year,month,day,hour,minute,sec)
c     2) Resolve time unit tunit as number of seconds per time unit, i.e.
c        for time unit = sec   tunit = 1
c        for time unit = hours tunit = 3600      
c        for time unit = days  tunit = 86400      
c     -----------------------------------------------------------------
      integer, intent(in)         :: ncid
      character*(*), intent(in)   :: time_vname
      integer, intent(out)        :: year,month,day,hour,minute,sec
      integer, intent(out)        :: tunit ! number of seconds per time unit
      integer                     :: varid,ihyph,it
      character*999               :: toffstr, text, date, time
      integer                     :: start(99), nwords
c     -----------------------------------------------------------------
c
c     resolve the time offset for the time variable in the data set
c     by parsing the attribute "time:units"
c
      call NetCDFcheck( nf90_inq_varid(ncid, time_vname, varid) )       ! fixed name
      call NetCDFcheck( nf90_get_att(ncid, varid, "units", toffstr) )   ! retrieve time:units attribute
      
c     --------- known variants ---------
c     1) "seconds since 1970-01-01T00:00:00Z"
c     2) "seconds since 1970-01-01 00:00:00"
c     3) "hours since 2014-04-01 01:00:00"
c     4) "seconds since 2014-4-1 00:00:00"

      call tokenize(toffstr, start, nwords) ! in string_tools.f
      if (nwords == 3) then     ! assume variant 1) above
         it = index(toffstr, "T")
         if (it<1) then
            write(*,*) "resolve_time_parameters: nwords==3, but no T"
            write(*,*) "unit attr = >",trim(adjustl(toffstr)),"<"
         endif
         toffstr(it:it) = " "      ! remove "T" for easier parsing below
         call tokenize(toffstr, start, nwords) ! reparse, robust toward spaces
      endif
      
      if (nwords /= 4) then
         write(*,*) "resolve_time_parameters: unable to parse unit attr"
         write(*,*) "unit attr = >",trim(adjustl(toffstr)),"<"
         stop
      endif
c     --- capture unit name ---
      text = ""
      read(toffstr(start(1):start(2)-1), *) text
      if (trim(adjustl(text))     == "seconds") then
         tunit = 1
      elseif (trim(adjustl(text)) == "minutes") then
         tunit = 60
      elseif (trim(adjustl(text)) == "hours") then
         tunit = 3600
      elseif (trim(adjustl(text)) == "days") then
         tunit = 86400
      else
         write(*,*) "resolve_time_parameters: unable to parse time unit"
         write(*,*) "time unit = >",trim(adjustl(text)),"<"
         stop
      endif
c     --- check "since" is next word ---
      text = ""
      read(toffstr(start(2):start(3)-1), *) text
      if (trim(adjustl(text)) /= "since") then
         write(*,*) "resolve_time_parameters: parse error of unit attr"
         write(*,*) "unit attr = >",trim(adjustl(toffstr)),"<"
         stop
      endif
c     --- capture offset date ---
      date = toffstr(start(3):start(4)-1)   ! like 2014-04-01 or 2014-4-1 
      read(date, 708, err=832) year         !
      ihyph = index(date, "-", back=.true.) ! last "-" in date token
      read(date(6:ihyph-1), *, err=832) month 
      read(date(ihyph+1:),  *, err=832) day
c     --- capture offset time-in-day ---      
      time = toffstr(start(4):)
      read(time, 711, err=833) hour,minute,sec
 708  format(i4)             ! 4 digit year expected
 711  format(i2,1x,i2,1x,i2) ! match known time variants
      
      return
      
 832  write(*,*) "resolve_time_parameters: error parsing offset date"
      write(*,*) "unit attr = >",trim(adjustl(toffstr)),"<"
      stop     
 833  write(*,*) "resolve_time_parameters: error parsing offset time"
      write(*,*) "unit attr = >",trim(adjustl(toffstr)),"<"
      stop
      
      end subroutine resolve_time_parameters


      
      
      subroutine attach_variable(xyzt, varname, ncid)
c     -----------------------------------------------------------------
c     Initialize container xyzt by attaching to the variable 
c     varname in necdf set ncid
c     -----------------------------------------------------------------
      type(coards_xyzt_variable), intent(out)  :: xyzt
      character*(*), intent(in)                :: varname
      integer, intent(in)                      :: ncid
c      
      character(len=999)                       :: cdummy
      integer                                  :: xtype, ndims, i
      integer                                  :: dimids(99) ! allow graceful handing, if dim mismatch
      integer                                  :: natts,istat
c     -----------------------------------------------------------------
      xyzt%ncid    = ncid
      xyzt%varname = varname
      call NetCDFcheck( nf90_inq_varid(ncid, varname, xyzt%varid))

      call NetCDFcheck( nf90_inquire_variable(ncid, xyzt%varid, cdummy,
     +                  xtype, ndims, dimids, natts) )
      if (ndims /= 4) then
         write(*,*) "attach_variable: for name=", trim(adjustl(varname))
         write(*,*) "attach_variable: unexpected ndims=", ndims
         stop
      endif
c     
      call NetCDFcheck(nf90_inquire_dimension(ncid, dimids(1),
     +                 len=xyzt%nx))
      call NetCDFcheck(nf90_inquire_dimension(ncid, dimids(2),
     +                 len=xyzt%ny))
      call NetCDFcheck(nf90_inquire_dimension(ncid, dimids(3),
     +                 len=xyzt%nz))
      call NetCDFcheck(nf90_inquire_dimension(ncid, dimids(4),
     +                 len=xyzt%nt))
c
c     --- set storage_type and pick up offset/scale for integer compressed storage ---  
c       
      if   ((xtype == NF90_SHORT).or.(xtype == NF90_INT)) then
         if (xtype==NF90_SHORT) then
            xyzt%storage_type = i2type  ! NF90_SHORT
         else
            xyzt%storage_type = i4type  ! NF90_INT
         endif
c         
c        common for integer compressed storage 
c         
         istat = nf90_get_att(ncid, xyzt%varid,"scale_factor",
     +        xyzt%scale)       ! does possibly not exist
         if (istat /= nf90_noerr) then
             write(*,*) "attach_variable: scale_factor mandatory for"//
     +                  "integer-compressed floats"
             stop
         endif
c 
         istat = nf90_get_att(ncid, xyzt%varid,"add_offset",
     +        xyzt%offset)       ! does possibly not exist
         if (istat /= nf90_noerr) then
             write(*,*) "attach_variable: add_offset mandatory for"//
     +                  "integer-compressed floats"
             stop
         endif          
c         
      elseif (xtype == NF90_FLOAT) then
         xyzt%storage_type = f4type
         
      elseif (xtype == NF90_DOUBLE) then
         xyzt%storage_type = f8type
         
      else
         write(*,*) "attach_variable: for name=", trim(adjustl(varname))
         write(*,*) "attach_variable: unhandled storage type", xtype
         stop
      endif
     
      end subroutine attach_variable

      
      
      subroutine load_xyz_frame_asreal(xyzt, rbuffer, iframe)
c     -----------------------------------------------------------------
c     Load frame iframe from xyzt into rbuffer
c     Float types relies on build-in cast
c     integer types are transformed as x_real = scale*x_int + offset
c     -----------------------------------------------------------------
      type(coards_xyzt_variable), intent(in)  :: xyzt
      real, intent(out)                       :: rbuffer(:,:,:)
      integer, intent(in)                     :: iframe
c
      integer(kind=2), allocatable            :: i2buffer(:,:,:) 
      integer(kind=4), allocatable            :: i4buffer(:,:,:)
      integer                                 :: start4D(4), count4D(4)
c     -----------------------------------------------------------------
c     --- check 1 <= iframe <= nt ---
      if ((iframe < 1).or.(iframe > xyzt%nt)) then
         write(*,*) "load_xyz_frame_asreal: for name=",
     +               trim(adjustl(xyzt%varname))
         write(*,*) "illegal iframe   = ", iframe
         write(*,*) "available frames = ", xyzt%nt
         stop
      endif
      if ((size(rbuffer,1) /= xyzt%nx).or.
     +    (size(rbuffer,2) /= xyzt%ny).or.
     +    (size(rbuffer,3) /= xyzt%nz)) then
         write(*,*) "load_xyz_frame_asreal: for name=",
     +               trim(adjustl(xyzt%varname))
         write(*,*) "buffer shape mismatch"
         write(*,*) "shape(rbuffer) = ", shape(rbuffer)
         write(*,*) "expects ", xyzt%nx,xyzt%ny,xyzt%nz
         stop
      endif
      
c     ---- load data, depending on storage type
      
      start4D = (/1,       1,       1,       iframe/) 
      count4D = (/xyzt%nx, xyzt%ny, xyzt%nz, 1/)
c
      if ((xyzt%storage_type == f4type).or.
     +    (xyzt%storage_type == f8type)) then 
         
         call NetCDFcheck( nf90_get_var(xyzt%ncid, xyzt%varid, rbuffer,  ! just load, rely on build-in recast for different floats
     +        count=count4D, start=start4D) )
         
      elseif (xyzt%storage_type == i2type) then     ! recast integer as x_real = scale*x_int + offset
         
         allocate( i2buffer(xyzt%nx, xyzt%ny, xyzt%nz) )
         call NetCDFcheck( nf90_get_var(xyzt%ncid, xyzt%varid, i2buffer,
     +        count=count4D, start=start4D) )
         rbuffer = xyzt%scale*i2buffer + xyzt%offset  ! whole array transformation
         deallocate( i2buffer )
         
      elseif (xyzt%storage_type == i4type) then     ! recast integer as x_real = scale*x_int + offset
         
         allocate( i4buffer(xyzt%nx, xyzt%ny, xyzt%nz) )
         call NetCDFcheck( nf90_get_var(xyzt%ncid, xyzt%varid, i4buffer,
     +        count=count4D, start=start4D) )
         rbuffer = xyzt%scale*i4buffer + xyzt%offset  ! whole array transformation
         deallocate( i4buffer )
         
      endif
      
      end subroutine load_xyz_frame_asreal

      
      
      subroutine detect_topography(xyzt, wmask, botlayer)
c     -----------------------------------------------------------------
c     Detect topography from xyzt (first frame)
c     by looking for _Fill ...
c     Topography assummed static inbetween time frames
c     Match (by comparison == ) for values corresponding to
c     _FillValue (precedence) or missing_value      
c     Either attribute _FillValue or missing_value must be set
c     Otherwise all data is assumed wet
c     For integer storage, make test in storage representation, not casted
c     Apply same test for z/sigma type data
c     -----------------------------------------------------------------
      type(coards_xyzt_variable), intent(in)  :: xyzt
      integer, intent(out)                    :: wmask(:,:)    ! nx,ny
      integer, intent(out)                    :: botlayer(:,:) ! nx,ny
      real                                    :: dry_real
      integer                                 :: dry_int
      real, allocatable                       :: rbuf(:,:,:)
      integer(kind=2), allocatable            :: i2buf(:,:,:)
      integer(kind=4), allocatable            :: i4buf(:,:,:)      
      integer                                 :: istat, ix, iy, iz
      integer                                 :: start4D(4), count4D(4)
c     -----------------------------------------------------------------
      start4D = (/1,       1,       1,       1/)  ! probe topography for first time frame
      count4D = (/xyzt%nx, xyzt%ny, xyzt%nz, 1/)
c
c     --- set botlayer: split up test, depending on storage type ---
c
      if ((xyzt%storage_type == i2type) .or.
     +    (xyzt%storage_type == i4type)) then
         
c        ======== integer type storage ========
         istat = nf90_get_att(xyzt%ncid, xyzt%varid,
     +                        "_FillValue", dry_int) ! precedence
         if (istat /= nf90_noerr) then
            istat = nf90_get_att(xyzt%ncid, xyzt%varid,"missing_value",
     +                           dry_int) ! secondary choice
            if (istat /= nf90_noerr) then
               write(*,*) "detect_topography: warning: either "//
     +                    "_FillValue or missing_value must be set."
               write(*,*) "detect_topography: warning: Assume all wet"
               botlayer = xyzt%nz ! grid all wet
               goto 553           ! omit dry value scan of data
               
            endif  ! istat test 2
         endif     ! istat test 1
         
c     --- now the dry value has been identified; apply test on each data point

         if (xyzt%storage_type == i2type) then
            
            allocate( i2buf(xyzt%nx, xyzt%ny, xyzt%nz) )
            call NetCDFcheck( nf90_get_var(xyzt%ncid, xyzt%varid, i2buf,
     +                     count=count4D, start=start4D) )   
         
            do ix = 1, xyzt%nx
               do iy = 1, xyzt%ny
                  do iz = 1, xyzt%nz
                     if (i2buf(ix,iy,iz) == dry_int) exit
                  enddo
                  botlayer(ix,iy) = iz-1 ! last wet layer in this column
               enddo
            enddo
            deallocate( i2buf )
            
         else   ! storage_type == i4type 

            allocate( i4buf(xyzt%nx, xyzt%ny, xyzt%nz) )
            call NetCDFcheck( nf90_get_var(xyzt%ncid, xyzt%varid, i4buf,
     +                     count=count4D, start=start4D) )   
         
            do ix = 1, xyzt%nx
               do iy = 1, xyzt%ny
                  do iz = 1, xyzt%nz
                     if (i4buf(ix,iy,iz) == dry_int) exit
                  enddo
                  botlayer(ix,iy) = iz-1 ! last wet layer in this column
               enddo
            enddo
            deallocate( i4buf )
            
         endif                  ! storage_type == integer type

         
      else
c        ======== real type storage ========
         
         istat = nf90_get_att(xyzt%ncid, xyzt%varid,
     +                        "_FillValue", dry_real) ! precedence
         if (istat /= nf90_noerr) then
            istat = nf90_get_att(xyzt%ncid, xyzt%varid,"missing_value",
     +                           dry_real) ! secondary choice
            if (istat /= nf90_noerr) then
               write(*,*) "detect_topography: warning: either "//
     +                    "_FillValue or missing_value must be set."
               write(*,*) "detect_topography: warning: Assume all wet"
               botlayer = xyzt%nz ! grid all wet
               goto 553           ! omit dry value scan of data
               
            endif   ! istat test 2
         endif      ! istat test 1
         
c        --- now the dry value has been identified; apply test on each data point
         allocate( rbuf(xyzt%nx, xyzt%ny, xyzt%nz) )
         call NetCDFcheck( nf90_get_var(xyzt%ncid, xyzt%varid, rbuf,
     +                     count=count4D, start=start4D) )   
     
         do ix = 1, xyzt%nx
            do iy = 1, xyzt%ny
               do iz = 1, xyzt%nz
                  if (rbuf(ix,iy,iz) == dry_real) exit
               enddo
               botlayer(ix,iy) = iz-1 ! last wet layer in this column
            enddo
         enddo
         deallocate( rbuf )
         
      endif                     ! xyzt%storage_type
c
c     --- set wmask consistently with botlayer ---
c
      
 553  continue                ! reentry point, if dry value scan was omitted
      
      where(botlayer > 0)
         wmask      = 1       ! wet
      elsewhere
         wmask      = 0       ! dry
      end where

      
      end subroutine detect_topography                   


      
      subroutine NetCDFcheck(status)
c     -----------------------------------------------------------------
c     This subroutine supports the recommended reading 
c     style of NetCDF4:
c       call NetCDFcheck( nf90_get_var(ncid, varid, data_in) )
c     -----------------------------------------------------------------
      integer, intent (in) :: status
c     -----------------------------------------------------------------    
      if(status /= nf90_noerr) then 
         print *, trim(nf90_strerror(status))
         stop "NetCDFcheck:Stopped"
      end if
      end subroutine NetCDFcheck


      subroutine print_variable_state(xyzt, iunit)
c     -----------------------------------------------------------------
c     -----------------------------------------------------------------
      type(coards_xyzt_variable), intent(in) :: xyzt
      integer, intent(in)                    :: iunit
      write(iunit, *)    "print_variable_state:"
      write(iunit, *)  "ncid    = ", xyzt%ncid
      write(iunit, *)  "varid   = ", xyzt%varid
      write(iunit, *)  "varname = ", trim(xyzt%varname)
      write(iunit, *) "nx,ny,nz,nt = ",
     +     xyzt%nx,  xyzt%ny,  xyzt%nz,  xyzt%nt
      write(iunit, *)  "storage_type = ",
     +     storage_type_names(xyzt%storage_type)
      
      if ((xyzt%storage_type == i2type) .or.
     +    (xyzt%storage_type == i4type)) then
         write(iunit, *)  "scale  = ",  xyzt%scale
         write(iunit, *)  "offset = ",  xyzt%offset
      endif

      end subroutine print_variable_state
      
      end module 
