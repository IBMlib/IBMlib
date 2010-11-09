      module ascii_writer
      use output
      use particles
            
      implicit none
      private
      
      type ascii_output_file
        character (len=999) filename
        type(output_var),pointer  :: vars(:)
        character (len=1) :: separator 
      end type ascii_output_file
      public :: ascii_output_file
      
      interface init_output
        module procedure init_output_ascii
      end interface
      public :: init_output

      interface write_frame
         module procedure write_frame_ascii   
      end interface
      public :: write_frame

      interface write_particles
         module procedure write_pars_ascii   
      end interface
      public :: write_particles

      
      contains
      
      
      subroutine init_output_ascii(of)
c------------------------------------------------------------  
      type(ascii_output_file),intent(inout) :: of
      type(particle) :: par
      integer i, iunit, ok
c------------------------------------------------------------        
      !Get metadata and insert into variable slot
      write(*,*) "Initialising ASCII outputfile : ",trim(of%filename)
      do i=1,size(of%vars)
          !Get the metadata
          call get_metadata(par,of%vars(i),ok)
          if(ok/=0) then
            write(*,*) "ERROR: Variable "//trim(of%vars(i)%name)// 
     +       " is not available."
            stop
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

      !Setup file headers
      call find_free_io_unit(iunit)
      open(unit=iunit,file=trim(of%filename))
      do i=1,size(of%vars)
        write(iunit,"(a)",advance="no") trim(of%vars(i)%name)
        if (i/=size(of%vars)) then   !Write separator
          write(iunit,"(a)",advance="no")  trim(of%separator)
        endif
      enddo
      close(iunit)
      end subroutine
      

      subroutine write_frame_ascii(of,par_ens)
c------------------------------------------------------------  
      type(ascii_output_file),intent(inout) :: of
      type(particle_ensemble), intent(in) :: par_ens
c     ------locals ------
      type(particle),pointer :: par
      integer ipar,iunit,last
c------------------------------------          
      !Find the last particle in the ensemble
      call get_last_particle_number(par_ens,last)

      !Open file
      call find_free_io_unit(iunit)
      open(unit=iunit,file=trim(of%filename),position="append")
      
      !Scan by particles
      do ipar=1,last
        !Get the particle object
        par=> get_particle(par_ens,ipar)
        !Write the particle
        call write_line_ascii(iunit,par,of)
      enddo
      close(iunit)
      end subroutine
      

      subroutine write_pars_ascii(of,par_ens,par_vec)
c------------------------------------------------------------  
      type(ascii_output_file),intent(inout) :: of
      type(particle_ensemble), intent(in) :: par_ens
      integer, intent(in) :: par_vec(:)      !Write these particle numbers
c     ------locals ------
      type(particle),pointer :: par
      integer ipar,iunit,npars
c------------------------------------          
      !Open file
      call find_free_io_unit(iunit)
      open(unit=iunit,file=trim(of%filename),position="append")
      
      !Scan by particles
      do ipar=1,size(par_vec)
        !Get the particle object
        par=> get_particle(par_ens,par_vec(ipar))
        !Write the particle
        call write_line_ascii(iunit,par,of)
      enddo
      close(iunit)
      end subroutine

 
      subroutine write_line_ascii(iunit,par,of)
c------------------------------------------------------------  
c     Constructs and writes a line in an ascii file corresponding to a single particle
      integer, intent(in) :: iunit
      type(particle),intent(in) :: par
      type(ascii_output_file),intent(inout) :: of
c     ------locals ------
      character (len=9999) :: outstr,tmpstr,fmtstr  
      integer i, ok
      type(output_var) :: var
c------------------------------------------------------------        
      !Extract variables and place in output file object in polytypes
      do i=1,size(of%vars)
        call get_property(par,of%vars(i),ok)
      enddo
      
      !Now concatenate variables into a single string
      outstr=""
      do i=1,size(of%vars)
        var=of%vars(i)
        !Format the string according to the data type
        select case(var%dat%type)
        case ("real")
            write(tmpstr,fmt=trim(var%fmt)) var%dat%r
        case ("int")
            write(tmpstr,fmt=trim(var%fmt)) var%dat%i
        case ("log")
            write(tmpstr,fmt=trim(var%fmt)) var%dat%l
        case ("str")
            write(tmpstr,fmt=trim(var%fmt)) var%dat%s
        case default
          write(*,*) "ERROR: Variable '"// trim(var%name) // "' has an "
     +       //"unknown type: "//trim(var%dat%type)
          stop
        end select
        !Concatentate string, add separator
        outstr = trim(outstr) // trim(tmpstr)
        if(i/=size(of%vars)) outstr = trim(outstr) // of%separator
      enddo
      !Now write the line to the unit (without the last separator)
      write(iunit,*) trim(outstr)

      end subroutine
      
      end module
