ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     NetCDF output writer module
c     ---------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
c     Allows configuration and writing of data to NetCDF files
C     As an initial framework, this module works on the 
c     assumption that the structure of the file is 2D, with
c     an unlimited time dimension and a limited dimension
c     corresponding to the particle ID in some form. All variables
c     stored are assumed to be 2D in nature ie functions of both
c     particles and time. Future developments maygeneralise this 
c     scheme, but it is hard to imagine many situations where
c     this is not sufficient.
c
c     Initial version by Mark Payne, DTU-Aqua, 1 October 2010
c
c     TODO: String support in NetCDF - good luck with that
c           Check outputs for legality
c           Can we setup the time_var as an array, rather than single variable
c------------------------------------------------------------  

      module netcdf_writer
      use output
      use netcdf
      use particles
      use input_parser
            
      implicit none
      private
      
      type netcdf_output_file
        private
        character (len=999) filename
        integer :: ncid
        integer, pointer :: varids(:)
        integer :: time_varid
        integer :: par_varid
        type(variable),pointer  :: par_var
        type(variable),pointer  :: time_var
        type(variable),pointer  :: vars(:)
        integer :: npars
        integer :: nvars
      end type netcdf_output_file
      public :: netcdf_output_file
      
      interface init_output
        module procedure init_output_netcdf
      end interface
      public :: init_output
      public :: set_var_range

      interface write_frame
         module procedure write_frame_netcdf   
      end interface
      public :: write_frame

      interface setup_output_from_file
         module procedure setup_output_from_file_netcdf
      end interface
      public :: setup_output_from_file

      integer, parameter :: NF90_SHORT_highlim = (2**15)-1
      integer, parameter :: NF90_SHORT_lowlim  = -NF90_SHORT_highlim+1
      integer, parameter :: NF90_SHORT_missval = NF90_FILL_SHORT

      integer, parameter :: NF90_INT_highlim = (2**31)-1
      integer, parameter :: NF90_INT_lowlim  = -NF90_INT_highlim+1
      integer, parameter :: NF90_INT_missval = NF90_FILL_INT
 
      
      contains
      
      
      subroutine init_output_netcdf(of,filename,vars,npars)
c------------------------------------------------------------  
      type(netcdf_output_file),intent(out) :: of
      character(*), intent(in) :: filename
      type(variable), intent(in) :: vars(:)
      integer, intent(in)    :: npars
      !----locals----
      type(particle) :: par
      integer i, ncid, ok, time_dimid, par_dimid,dimids(2)
      type(variable) :: time_var, par_var
c------------------------------------------------------------        
      !Initialise file object
      write(*,*) "Initialising netcdf outputfile : ",trim(filename)
      of%filename =filename
      of%nvars = size(vars)
      of%npars = npars
      allocate(of%vars(of%nvars))
      allocate(of%par_var)
      allocate(of%time_var)
      allocate(of%varids(of%nvars))
      of%vars = vars
      !Setup time and particle variables
      call get_metadata("POSIX",time_var)
      call get_metadata("tracerID",par_var)
      of%par_var = par_var
      of%time_var = time_var
      !Write a description of the file
      do i=1,of%nvars
          write(*,380) i, trim(get_name(of%vars(i))),
     +       trim(get_desc(of%vars(i))) 
      enddo
380   format("  Variable ",i2, ": ",a," (",a,")")      
      !Test to see what type of variable the particle counter is - must be integer
      if(get_type(of%par_var)/="int") then
        call abort_run("init_output_netcdf","Particle variable, '" 
     +       //trim(get_name(of%par_var)) //"' must be integer but is '"
     +       //get_type(of%par_var)//"'")
      end if

      ! Create the file. 
      call check( nf90_create(trim(of%filename), nf90_clobber, ncid))
      ! Define the dimensions - a particle dimension and an unlimted datetime dimension
      call check( nf90_def_dim(ncid, get_name(of%par_var),
     +        of%npars,par_dimid))
      call check( nf90_def_dim(ncid, get_name(of%time_var), 
     +                          NF90_UNLIMITED, time_dimid) )
      ! Define the coordinate variables    
      call check( nf90_def_var(ncid,get_name(of%par_var), 
     +     get_nf90_type(of%par_var),par_dimid, of%par_varid))
      call check( nf90_def_var(ncid,get_name(of%time_var),
     +     get_nf90_type(of%time_var),time_dimid,of%time_varid))

c      ! Now define the other variables
      dimids = (/ par_dimid,time_dimid /)
      do i=1,of%nvars
       call check(nf90_def_var(ncid,get_name(of%vars(i)),
     +        get_nf90_type(of%vars(i)),dimids,of%varids(i)))
      enddo
      ! Now define the variable attributes (metadata)
      do i=1,of%nvars
        call netcdf_define_attributes(ncid,of%vars(i),of%varids(i))
      enddo
      call netcdf_define_attributes(ncid,
     +          of%par_var,of%par_varid)
      call netcdf_define_attributes(ncid,
     +          of%time_var,of%time_varid)

      ! End define mode.
      call check( nf90_enddef(ncid) )
      ! Write the particle coordinate variable data. 
      call check( nf90_put_var(ncid, of%par_varid,(/(i,i=1,of%npars)/)))
      ! Close the file. 
      call check( nf90_close(ncid))

      end subroutine
      

      subroutine write_frame_netcdf(of,par_ens)
c     ------------------------------------------------------------  
      type(netcdf_output_file),intent(inout) :: of
      type(particle_ensemble), intent(in) :: par_ens
c     ------locals ------
      type(particle),pointer :: par
      type(polytype) :: bucket
      integer last, ipar, time_dimid, time_recs, i
      integer par_ID,ok
c     ------------------------------------          
      !Find the last particle in the ensemble
      call get_last_particle_number(par_ens,last)

      !Open NetCDF file for writing
      call check( nf90_open(of%filename, NF90_WRITE, of%ncid) )
      !Figure out how many time records there are already
      call check( nf90_inq_dimid(of%ncid,
     +          get_name(of%time_var),time_dimid))
      call check( nf90_inquire_dimension(of%ncid,
     +          time_dimid,len=time_recs))
      !Get the variable ids for each variable
      !This is a precautionary approach, incase the varids have changed
      do i=1,size(of%vars)
        call check( nf90_inq_varid(of%ncid,
     +       get_name(of%vars(i)),of%varids(i)))
      enddo 
      call check( nf90_inq_varid(of%ncid,
     +       get_name(of%time_var),of%time_varid))
      !Scan by particles
      do ipar=1,last
        !Get the particle object
        par=> get_particle(par_ens,ipar)
        !Write the particle in the new record, according to particle number
        call write_par_netcdf(of,par,time_recs+1)
      enddo
c
      !Finally, write the unlimited (time) variable data
      call get_property(of%time_var,bucket,ok)
      call prepare_variable_netcdf(of%time_var,bucket)
      call write_point_netcdf(of%ncid,of%time_varid,bucket,
     +          t_start=time_recs+1)

      !Close netcdf file    
      call check( nf90_close(of%ncid) )

      end subroutine
      

      subroutine write_par_netcdf(of,par,time_idx)
c------------------------------------------------------------  
c     Writes the parameters in a netcdf file that correspond
c     to a single particle. Parameters to write are specified 
c     in the output file object
      type(netcdf_output_file),intent(inout) :: of
      integer, intent(in) :: time_idx
      type(particle),intent(in) :: par
c     ------locals ------
      integer i, ok,p_idx
      type(polytype) :: bucket
c------------------------------------------------------------        
      !Figure out where to write the data
      call get_property(par,of%par_var,bucket,ok)
      p_idx = get_data_int(bucket)
      !Loop over supplied variables
      do i=1,size(of%vars)
         !Get the data
         call get_property(par,of%vars(i),bucket,ok)
         !Prepare the variable for output
         call prepare_variable_netcdf(of%vars(i),bucket) 
         !Now write the result
         call write_point_netcdf(of%ncid,of%varids(i),bucket,
     +           t_start=time_idx,p_start=p_idx)
      enddo
      end subroutine


      subroutine write_point_netcdf(ncid,varid,bucket,t_start,p_start)
c------------------------------------------------------------  
c     Writes the data contained in bucket to a single point in a netcdf file
      integer, intent(in) ::  varid, ncid
      type(polytype), intent(in) :: bucket
      integer, intent(in),optional :: t_start,p_start
c     ------locals ------
      integer, allocatable :: start(:),cnt(:)
      integer :: status
c------------------------------------------------------------        
      !Decide whether we are writing 1D or 2D data
      if(present(t_start) .and. present(p_start)) then
         allocate(start(2),cnt(2))
         start = (/p_start,t_start/)
         cnt=(/1,1/)
      else if (present(t_start)) then
         allocate(start(1),cnt(1))
         start = (/t_start/)
         cnt=(/1/)
      else if (present(p_start)) then
         allocate(start(1),cnt(1))
         start = (/p_start/)
         cnt=(/1/)
      else
         call abort_run("write_point_netcdf()","At least one of "
     +    // "par_start and time_start need to be specified.")
      endif
      !Writing depends in the first instance on the data type
      select case(get_type(bucket))
      case ("real")
         status= nf90_put_var(ncid, varid,
     +     (/get_data_real(bucket)/), start=start,count=cnt)
      case ("int")
         status = nf90_put_var(ncid, varid,
     +     (/get_data_int(bucket)/),  start=start,count=cnt)
      case default
        call abort_run("write_point_netcdf()","Data of type '" //
     +   trim(get_type(bucket))// "' is not currently handled.")
      end select
      !Check write status. If no error, continue. If data is not 
      !representable (nf90_erange), then do nothing and continue - 
      !this is equivalent to writing an NA
      if(status/=nf90_noerr.and.status/=nf90_erange ) then
         call check(status)
         stop 
      endif
      !Tidy up to finish with
      deallocate(start)
      deallocate(cnt)
      end subroutine write_point_netcdf


      subroutine prepare_variable_netcdf(var,bucket)
c------------------------------------------------------------  
c     Prepares an output variable for writing to NETCDF by taking
c     care of type conversions and the like
      type(variable),intent(in) :: var
      type(polytype),intent(inout) :: bucket
c     ------locals ------
      integer i, write_int,ok,missval
      real write_real, scaled,ulim,llim, offset,scalar
      logical write_byte
      character(len=999) bucket_name
c------------------------------------------------------------  
      offset = get_offset(var)
      scalar = get_scalar(var)
      bucket_name = get_name(bucket)
c     !We have to split the processing up according to the 
c     !type that will be stored in the netcdf file
      select case( get_nf90_type(var))
      case (NF90_SHORT, NF90_INT) !-------------------------------------
         !Avoid writing ints that are too large 
         if (get_nf90_type(var)==NF90_SHORT) then
           ulim = NF90_SHORT_highlim
           llim = NF90_SHORT_lowlim 
           missval = NF90_SHORT_missval
         elseif (get_nf90_type(var)==NF90_INT) then
           ulim = NF90_INT_highlim
           llim = NF90_INT_lowlim
           missval = NF90_INT_missval
         endif
        !Cast into (scaled) integer format
        select case(get_type(bucket))
        case ("real")
          scaled=(get_data_real(bucket)-offset)/scalar
        case ("int")
          scaled =(real(get_data_int(bucket))-offset)/scalar
        case default
          call abort_run("prepare_variable_netcdf()","Variable '"// 
     +    trim(get_name(var)) // "' is of"
     +     //" type '" //trim(get_type(bucket))//"' and cannot "
     +     //"be cast as NF90_SHORT or NF90_INT")
         end select
         !Check for out of bounds - set to missing value if so
         if(scaled >ulim .OR. scaled < llim) then
            write_int = missval 
         else
            write_int=nint(scaled)
         endif
         !Now reconstruct the bucket
         call construct(bucket,bucket_name,write_int)
      case (NF90_FLOAT, NF90_DOUBLE) !-------------------------------------
        !Cast into (scaled) real format
         select case(get_type(bucket))
         case ("real")
           write_real = (get_data_real(bucket)-offset)/scalar
         case ("int")
           write_real = (real(get_data_int(bucket))-offset)/scalar
         case default
           call abort_run("prepare_variable_netcdf()"," Variable '"// 
     +      trim(get_name(var)) // "' is of"
     +      //" type '" //trim(get_type(bucket))//"' and cannot "
     +      //"be cast as NF90_FLOAT or NF90_DOUBLE")
         end select
         !Now construct the bucket
         call construct(bucket,bucket_name,write_real)
      case default!----------------------------------------------------------
         call abort_run("prepare_variable_netcdf()","Variable '"
     +       //trim(get_name(var)) // "' is of unknown type ")
      end select
      end subroutine prepare_variable_netcdf


      subroutine netcdf_define_attributes(ncid,var,varid)
c------------------------------------------------------------  
      integer, intent(in) :: ncid,varid
      type(variable),intent(in) :: var
c------------------------------------------------------------        
      call check(nf90_put_att(ncid,varid,
     +           "long_name",get_desc(var)))
      call check(nf90_put_att(ncid,varid,
     +           "units",get_units(var)))
      call check(nf90_put_att(ncid,varid,
     +           "scale_factor",get_scalar(var)))
      call check(nf90_put_att(ncid,varid,
     +           "add_offset",get_offset(var)))
      !Define missing value, if appropriate
      select case(get_nf90_type(var))
      case (NF90_SHORT)
         call check(nf90_put_att(ncid,varid,
     +           "_FillValue",NF90_FILL_SHORT))
         call check(nf90_put_att(ncid,varid,
     +           "missing_value",NF90_SHORT_missval))
      case (NF90_INT)
         call check(nf90_put_att(ncid,varid,
     +           "_FillValue",NF90_FILL_INT))
         call check(nf90_put_att(ncid,varid,
     +           "missing_value",NF90_INT_missval))
      case (NF90_FLOAT)
         call check(nf90_put_att(ncid,varid,
     +           "_FillValue",NF90_FILL_FLOAT))
         call check(nf90_put_att(ncid,varid,
     +           "missing_value",NF90_FILL_FLOAT))
      case (NF90_DOUBLE)
         call check(nf90_put_att(ncid,varid,
     +           "_FillValue",NF90_FILL_DOUBLE))
         call check(nf90_put_att(ncid,varid,
     +           "missing_value",NF90_FILL_DOUBLE))
      end select

      end subroutine


      subroutine check(status)
c     -------------------------------------------------------
c     Checks that NetCDF operations worked
c     -------------------------------------------------------
      integer, intent ( in) :: status
c     --------------------------  
      if(status /= nf90_noerr) then 
         write(*,*) "NF90 error code : ",status
         write(*,*) trim(nf90_strerror(status))
         stop 
      end if
      end subroutine check 


      integer function get_nf90_type(var)
c     -------------------------------------------------------
c     Converts the type, text string of an output variable to
c     an integer constant definted by the NetCDF module. Some
c     default values for other types (real, int) are also 
c     accepted as a default mapping
c     -------------------------------------------------------
      type(variable), intent(in) :: var
      select case (get_type(var))
      case("NF90_FLOAT","real") 
         get_nf90_type=NF90_FLOAT
      case("NF90_INT","int") 
         get_nf90_type=NF90_INT
      case("NF90_DOUBLE") 
         get_nf90_type=NF90_DOUBLE
      case("NF90_SHORT") 
         get_nf90_type=NF90_SHORT
      case default
         call abort_run("get_nf90_type","Type '" //
     +    trim(get_type(var)) // "' ('" //
     +    trim(get_name(var))//"' variable) is not recognised by NF90")
      end select

      end function


      subroutine set_var_range(this,rng)
c     -------------------------------------------------------
c     Sets up output variable parameters that are specific 
c     to NetCDF operations worked
c     Specifying a range and a var_type overrules any offset and
c     scalar specification
c     -------------------------------------------------------
      type(variable), intent(inout) :: this
      real, intent(in) :: rng(2)
      real :: highx, lowx,scalar,offset
c     -------------------------------------------------------
      select case(get_nf90_type(this))
      case (NF90_SHORT)
         highx=NF90_SHORT_highlim
         lowx= NF90_SHORT_lowlim
      case (NF90_INT)
         highx=NF90_INT_highlim
         lowx= NF90_INT_lowlim
      case default
         call abort_run("set_var_range","Setting a range "
     +      // "does not have a meaning for variables that are "
     +      // "not NF90_SHORT or NF90_INT. Specify a scalar and "
     +      // "offset directly instead.")
      end select
      scalar = (maxval(rng)-minval(rng))/(highx-lowx)
      offset = minval(rng)-scalar*lowx
      call set_scalar(this,scalar)
      call set_offset(this,offset)
      end subroutine

  
      subroutine setup_output_from_file_netcdf(of,par_ens,
     +            ctrl,fname_tag,var_tag)
      !-------------------------------------------------------  
      !Creates an ascii_output_file by reading configuration
      !tags from a configuration file
      type(netcdf_output_file),intent(inout) :: of
      type(control_file), intent(in) :: ctrl
      character(*), intent(in) :: fname_tag, var_tag
      type(particle_ensemble), intent(in) :: par_ens
      !------locals ------
      integer :: nvars,i,ihit,nwords, start(256),npars
      real    :: var_range(2)
      character*999 :: strbuf,fname,var_type,var_name
      type(variable), allocatable :: vars(:)
      type(variable) :: time_var, par_var
      !-------------------------------------------------------  
      !Count the number of varables in the ctrlfile and allocate
      nvars = count_tags(ctrl,var_tag)                
      allocate(vars(nvars))

      !Read tags sequentially and create variables
      ihit=1
      do i=1,nvars
        call read_control_data(ctrl,var_tag,strbuf,ihit)
        ihit = ihit +1 
        !Split strbuf into variables as required
        call tokenize(strbuf, start, nwords)
        read(strbuf(start(1):),*) var_name 
        !Now create output variable
        call get_metadata(adjustl(var_name),vars(i))
        !Set variable type if specified
        if(nwords>=2) then
          read(strbuf(start(2):),*) var_type
          call set_type(vars(i),var_type)
        endif
        !Set variable range if specified
        if(nwords>=4) then
          read(strbuf(start(3):),*) var_range(1)
          read(strbuf(start(4):),*) var_range(2)
          call set_var_range(vars(i),var_range)
        endif

      enddo

      !Get filename 
      call read_control_data(ctrl,fname_tag,fname)
      fname = adjustl(trim(fname))

      !Setup file
      npars = get_ensemble_size(par_ens)
      call init_output(of,fname,vars,npars)

      deallocate(vars)
      end subroutine

      end module
