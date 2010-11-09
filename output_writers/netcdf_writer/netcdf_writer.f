ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     NetCDF output writer module
c     ---------------------------------------------------
c     $Rev: 42 $
c     $LastChangedDate: 2010-09-28 14:18:09 +0200 (ti, 28 sep 2010) $
c     $LastChangedBy: mpay $ 
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
            
      implicit none
      private
      
      type netcdf_output_file
        character (len=999) filename
        type(output_var),pointer  :: time_var
        integer :: n_particles
        type(output_var),pointer  :: par_var
        type(output_var),pointer  :: vars(:)
      end type netcdf_output_file
      public :: netcdf_output_file
      
      interface init_output
        module procedure init_output_netcdf
      end interface
      public :: init_output

      interface write_frame
         module procedure write_frame_netcdf   
      end interface
      public :: write_frame

      interface write_particles
         module procedure write_par_netcdf   
      end interface
      public :: write_particles

      interface construct
        module procedure construct_outvar_ncdf
      end interface
      public :: construct
      
      contains
      
      
      subroutine init_output_netcdf(of)
c------------------------------------------------------------  
      type(netcdf_output_file),intent(inout) :: of
      type(particle) :: par
      integer i, ncid, ok, time_dimid, par_dimid,varid,dimids(2)
c------------------------------------------------------------        
      !Get metadata and insert into variable slot
      write(*,*) "Initialising netcdf outputfile : ",trim(of%filename)
      do i=1,size(of%vars)
          !Get the metadata
          call get_metadata(par,of%vars(i),ok)
          if(ok/=0) then
            call abort("init_output_netcdf","Variable '"
     +       //trim(of%vars(i)%name)// "' is not available.")
          endif
          !Write a description of the file
          write(*,380) i, trim(of%vars(i)%metadata%name),
     +       trim(of%vars(i)%metadata%desc)
          !If no format information supplied, then default to that from the metadata
          if(len_trim(of%vars(i)%fmt)==0) then
            of%vars(i)%fmt=of%vars(i)%metadata%default_fmt
          endif
      enddo
380   format("  Variable ",i2, ": ",a," (",a,")")      

      !Test to see what type of variable the particle counter is - must be integer
      call get_metadata(par,of%par_var,ok)
      if(ok/=0) then 
        call abort("init_output_netcdf","Cannot find variable "// 
     +           of%par_var%name)
      end if
      if(of%par_var%metadata%type/="int") then
        call abort("init_output_netcdf","Particle variable, '" 
     +       //trim(of%par_var%name) //"' must be integer but is '"
     +       //of%par_var%metadata%type//"'")
      end if

      ! Create the file. 
      call check( nf90_create(trim(of%filename), nf90_clobber, ncid))
      ! Define the dimensions - a particle dimension and an unlimted datetime dimension
      call check( nf90_def_dim(ncid, of%time_var%name, 
     +                          NF90_UNLIMITED, time_dimid) )
      call check( nf90_def_dim(ncid, of%par_var%name,
     +        of%n_particles,par_dimid))
      ! Define the coordinate variables    
      call check( nf90_def_var(ncid,of%par_var%name, 
     +           of%par_var%var_type,par_dimid, of%par_var%varid))
      call check( nf90_def_var(ncid,of%time_var%name,
     +           of%time_var%var_type, time_dimid, of%time_var%varid))

      ! Now define the other variables
      dimids = (/ par_dimid,time_dimid /)
      do i=1,size(of%vars)
       call check(nf90_def_var(ncid,of%vars(i)%name,of%vars(i)%var_type,
     +        dimids,of%vars(i)%varid))
      enddo
      ! Now define the variable attributes (metadata)
      do i=1,size(of%vars)
        call netcdf_define_attributes(ncid,of%vars(i))
      enddo
      call netcdf_define_attributes(ncid,of%par_var)
      call netcdf_define_attributes(ncid,of%time_var)

      ! End define mode.
      call check( nf90_enddef(ncid) )
      ! Write the particle coordinate variable data. 
      call check( nf90_put_var(ncid, of%par_var%varid, 
     +     (/(i,i=1,of%n_particles)/)))
      ! Close the file. 
      call check( nf90_close(ncid) )

      end subroutine
      

      subroutine write_frame_netcdf(of,par_ens)
c------------------------------------------------------------  
      type(netcdf_output_file),intent(inout) :: of
      type(particle_ensemble), intent(in) :: par_ens
c     ------locals ------
      type(particle),pointer :: par
      integer last, ipar, time_dimid, time_recs, ncid,i
      integer par_ID,ok
c------------------------------------          
      !Find the last particle in the ensemble
      call get_last_particle_number(par_ens,last)

      !Open NetCDF file for writing
      call check( nf90_open(of%filename, NF90_WRITE, ncid) )
      !Figure out how many time records there are already
      call check( nf90_inq_dimid(ncid, of%time_var%name, time_dimid) )
      call check( nf90_inquire_dimension(ncid,time_dimid,len=time_recs))
      !Get the variable ids for each variable
      !I'm not sure if this is strictly necessary, but it seems wise
      do i=1,size(of%vars)
        call check( nf90_inq_varid(ncid,
     +     of%vars(i)%name, of%vars(i)%varid)) 
      enddo 
      call check( nf90_inq_varid(ncid,of%time_var%name,
     +      of%time_var%varid)) 
      !Scan by particles
      do ipar=1,last
        !Get the particle object
        par=> get_particle(par_ens,ipar)
        !Write the particle in the new record, according to particle number
        call write_par_netcdf(ncid,par,of,time_recs+1)
      enddo

      !Finally, write the unlimited (time) variable data
      call get_property(of%time_var,ok)
      call prepare_variable_netcdf(of%time_var)
      call write_point_netcdf(ncid,of%time_var,time_start=time_recs+1)

      !Close netcdf file    
      call check( nf90_close(ncid) )

      end subroutine
      

      subroutine write_par_netcdf(ncid,par,of,time_idx)
c------------------------------------------------------------  
c     Writes the parameters in a netcdf file that correspond
c     to a single particle. Parameters to write are specified 
c     in the output file object
      integer, intent(in) :: ncid,time_idx
      type(particle),intent(in) :: par
      type(netcdf_output_file),intent(inout) :: of
c     ------locals ------
      integer i, ok
c------------------------------------------------------------        
      !Figure out where to write the data
      call get_property(par,of%par_var,ok)
      !Loop over supplied variables
      do i=1,size(of%vars)
         !Get the data
         call get_property(par,of%vars(i),ok)
         !Prepare the variable for output
         call prepare_variable_netcdf(of%vars(i)) 
         !Now write the result
         call write_point_netcdf(ncid,of%vars(i),
     +           time_start=time_idx,par_start=of%par_var%dat%i)
      enddo
      end subroutine


      subroutine write_point_netcdf(ncid,var,time_start,par_start)
c------------------------------------------------------------  
c     Writes the formatted data contained in a single variable to a single
c     variable object to a single point in a netcdf file
      integer, intent(in) :: ncid
      type(output_var), intent(in) :: var
      integer, intent(in),optional :: time_start,par_start
c     ------locals ------
      integer, allocatable :: start(:),cnt(:)
c------------------------------------------------------------        
      !Decide whether we are writing 1D or 2D data
      if(present(time_start) .and. present(par_start)) then
         allocate(start(2),cnt(2))
         start = (/par_start,time_start/)
         cnt=(/1,1/)
      else if (present(time_start)) then
         allocate(start(1),cnt(1))
         start = (/time_start/)
         cnt=(/1/)
      else if (present(par_start)) then
         allocate(start(1),cnt(1))
         start = (/par_start/)
         cnt=(/1/)
      else
         call abort("write_point_netcdf()","At least one of "
     +    // "par_start and time_start need to be specified.")
      endif
      !Writing depends in the first instance on the variable type
      select case(var%dat%type)
      case ("real")
         call check( nf90_put_var(ncid, var%varid,(/var%dat%r/),
     +     start=start,count=cnt))
      case ("int")
         call check( nf90_put_var(ncid, var%varid,(/var%dat%i/),
     +     start=start,count=(/1,1/)))
      case default
        call abort("write_point_netcdf()"," Variable '"//
     +   trim(var%name) // "' is of type '" //trim(var%dat%type)//
     +   "' and is not currently handled.")
      end select
      !Tidy up to finish with
      deallocate(start)
      deallocate(cnt)
      end subroutine write_point_netcdf


      subroutine prepare_variable_netcdf(var)
c------------------------------------------------------------  
c     Prepares an output variable for writing to NETCDF by taking
c     care of type conversions and the like
      type(output_var),intent(inout) :: var
c     ------locals ------
      integer i, write_int,ok
      real write_real, scaled,ulim,llim
      logical write_byte
      type(polytype) :: bucket
c------------------------------------------------------------  
c     !We have to split the processing up according to the 
c     !type that will be stored in the netcdf file
      select case( var%var_type)
      case (NF90_SHORT, NF90_INT) !-------------------------------------
         !Avoid writing ints that are too large 
         if (var%var_type==NF90_SHORT) then
           ulim=2**15-1
           llim =-ulim+1  !-32767 is a missing value for NF90_SHORT 
         elseif (var%var_type==NF90_INT) then
           ulim=2**31-1
           llim =-ulim+1 
         endif
        !Cast into (scaled) integer format
        select case(var%dat%type)
        case ("real")
          scaled =(var%dat%r-var%offset)/var%scalar
          write_int=nint(min(ulim,max(llim,scaled)))
        case ("int")
          scaled =(real(var%dat%i)-var%offset)/var%scalar
          write_int=nint(min(ulim,max(llim,scaled)))
        case default
          call abort("prepare_variable_netcdf()","Variable '"// 
     +    trim(var%name) // "' is of"
     +     //" type '" //trim(var%dat%type)//"' and cannot "
     +     //"be cast as NF90_SHORT or NF90_INT")
         end select
         !Now construct the bucket
         call construct(bucket,var%dat%name,write_int)
      case (NF90_FLOAT, NF90_DOUBLE) !-------------------------------------
        !Cast into real format
         select case(var%dat%type)
         case ("real")
           write_real = real((var%dat%r-var%offset)/var%scalar)
         case ("int")
           write_real = real((var%dat%i-var%offset)/var%scalar)
         case default
           call abort("prepare_variable_netcdf()"," Variable '"// 
     +      trim(var%name) // "' is of"
     +      //" type '" //trim(var%dat%type)//"' and cannot "
     +      //"be cast as NF90_FLOAT or NF90_DOUBLE")
         end select
         !Now construct the bucket
         call construct(bucket,var%dat%name,write_real)
      case default!----------------------------------------------------------
         call abort("prepare_variable_netcdf()","Variable '"
     +       //var%name // "' is of unknown type ")
      end select
      !Replace the data in the variable with the corresponding bucket
      var%dat = bucket
      end subroutine prepare_variable_netcdf


      subroutine netcdf_define_attributes(ncid,var)
c------------------------------------------------------------  
      integer, intent(in) :: ncid
      type(output_var),intent(in) :: var
c------------------------------------------------------------        
      call check(nf90_put_att(ncid,var%varid,
     +           "long_name",var%metadata%desc))
      call check(nf90_put_att(ncid,var%varid,
     +           "units",var%metadata%units))
      call check(nf90_put_att(ncid,var%varid,
     +           "scale_factor",var%scalar))
      call check(nf90_put_att(ncid,var%varid,
     +           "add_offset",var%offset))
      end subroutine



      subroutine check(status)
c     -------------------------------------------------------
c     Checks that NetCDF operations worked
c     -------------------------------------------------------
      integer, intent ( in) :: status
c     --------------------------  
      if(status /= nf90_noerr) then 
         write(*,*) trim(nf90_strerror(status))
         stop 
      end if
      end subroutine check 


      subroutine construct_outvar_ncdf(this,name,fmt,var_type,
     +    scalar,offset,range)
c     -------------------------------------------------------
c     Sets up output variable parameters that are specific 
c     to NetCDF operations worked
c     Specifying a range and a var_type overrules any offset and
c     scalar specification
c     -------------------------------------------------------
      type(output_var), intent(inout) :: this
      character(len=*), intent(in) :: name
      character(len=*), intent(in),optional :: fmt
      integer, intent(in) :: var_type
      real, optional, intent(in) :: scalar, offset, range(2)
      real lowx,highx
c     -------------------------------------------------------
      this%name=name
      this%var_type=var_type
      this%scalar =1.0
      this%offset=0.0
      this%fmt=""
      if(present(fmt)) this%fmt=fmt
      if(present(scalar)) this%scalar= scalar
      if(present(offset)) this%offset= offset
      !If a range is specified override the scalar and offset
      if(present(range)) then
        select case(var_type)
        case (NF90_SHORT)
          highx=(2**15)-1
          lowx=-highx 
        case (NF90_INT)
          highx=(2**31)-1
          lowx=-highx 
        case default
          call abort("construct_outvar_ncdf","Supplying a range "
     +      // "does not have a meaning for variables that are "
     +      // "not NF90_SHORT or NF90_INT. Specify a scalar and "
     +      // "offset instead.")
        end select
        this%scalar=(maxval(range)-minval(range))/(highx-lowx)
        this%offset=minval(range)-this%scalar*lowx
      endif
      end subroutine construct_outvar_ncdf


  
      end module
