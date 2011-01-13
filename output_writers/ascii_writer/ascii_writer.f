ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ---------------------------------------------------
c     ASCII output writer module
c     ---------------------------------------------------
c     $Rev$
c     $LastChangedDate$
c     $LastChangedBy$ 
c
c     Allows configuration and writing of data to ASCII files
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module ascii_writer
      use output
      use particles
      use input_parser
            
      implicit none
      private
      
      type ascii_output_file
        private
        character (len=999) filename
        integer :: iunit
        type(variable),pointer  :: vars(:)
        integer :: nvars
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

      interface setup_output_from_file
         module procedure setup_output_from_file_ascii
      end interface
      public :: setup_output_from_file
      
      contains
      
      
      subroutine init_output_ascii(of,filename,separator,vars)
c------------------------------------------------------------  
      type(ascii_output_file),intent(out) :: of
      character(*), intent(in) :: filename
      character(*), intent(in) :: separator
      type(variable), intent(in) :: vars(:)
      !-----locals-----
      integer i,  ok
c------------------------------------------------------------        
      !Initialise ASCII file
      write(*,*) "Initialising ASCII outputfile : ",trim(filename)
      of%filename = filename 
      of%separator= separator
      of%nvars = size(vars)
      allocate(of%vars(of%nvars))
      of%vars = vars
      !Display variable information
      do i=1,of%nvars
          !Write a description of the file
          write(*,380) i, trim(get_name(of%vars(i))),
     +       trim(get_desc(of%vars(i))) 
      enddo
380   format("  Variable ",i2, ": ",a," (",a,")")      

      !Setup file headers
      call find_free_io_unit(of%iunit)
      open(unit=of%iunit,file=trim(of%filename))
      do i=1,of%nvars
        write(of%iunit,"(a)",advance="no") trim(get_name(of%vars(i)))
        if (i/=of%nvars) then   !Write separator
          write(of%iunit,"(a)",advance="no")  trim(of%separator)
        endif
      enddo
      close(of%iunit)
      end subroutine
      

      subroutine write_frame_ascii(of,par_ens)
c------------------------------------------------------------  
      type(ascii_output_file),intent(inout) :: of
      type(particle_ensemble), intent(in) :: par_ens
c     ------locals ------
      type(particle),pointer :: par
      integer ipar,last
c------------------------------------          
      !Find the last particle in the ensemble
      call get_last_particle_number(par_ens,last)

      !Open file
      call find_free_io_unit(of%iunit)
      open(unit=of%iunit,file=trim(of%filename),position="append")
      
      !Scan by particles
      do ipar=1,last
        !Get the particle object
        par=> get_particle(par_ens,ipar)
        !Write the particle
        call write_par_ascii(of,par)
      enddo
      close(of%iunit)
      end subroutine
      

      subroutine write_pars_ascii(of,par_ens,par_vec)
c------------------------------------------------------------  
      type(ascii_output_file),intent(inout) :: of
      type(particle_ensemble), intent(in) :: par_ens
      integer, intent(in) :: par_vec(:)      !Write these particle numbers
c     ------locals ------
      type(particle),pointer :: par
      integer ipar,npars
c------------------------------------          
      !Open file
      call find_free_io_unit(of%iunit)
      open(unit=of%iunit,file=trim(of%filename),position="append")
      
      !Scan through particles
      do ipar=1,size(par_vec)
        !Get the particle object
        par=> get_particle(par_ens,par_vec(ipar))
        !Write the particle
        call write_par_ascii(of,par)
      enddo
      close(of%iunit)
      end subroutine

 
      subroutine write_par_ascii(of,par)
c------------------------------------------------------------  
c     Constructs and writes a line in an ascii file corresponding to a single particle
      type(ascii_output_file),intent(inout) :: of
      type(particle),intent(in) :: par
c     ------locals ------
      character (len=9999) :: outstr,tmpstr,fmtstr  
      integer i, ok
      type(variable) :: var
      type(polytype) :: bucket
c------------------------------------------------------------        
      !Loop through output variables to build the output line
      outstr=""
      !Extract variables and place in output file object in polytypes
      do i=1,of%nvars
        call get_property(par,of%vars(i),bucket,ok)
        !Concatentate string, add separator
        call get_data_string(bucket,get_fmt(of%vars(i)),tmpstr)
        outstr = trim(outstr) // trim(tmpstr)
        if(i/=size(of%vars)) outstr = trim(outstr) // of%separator
      enddo
      !Now write the line to the unit (without the last separator)
      write(of%iunit,*) trim(outstr)
      end subroutine

      
      subroutine setup_output_from_file_ascii(of,
     +                   ctrl,fname_tag,var_tag)
      !-------------------------------------------------------  
      !Creates an ascii_output_file by reading configuration
      !tags from a configuration file
      type(ascii_output_file),intent(inout) :: of
      type(control_file), intent(in) :: ctrl
      character(*), intent(in) :: fname_tag, var_tag
      !------locals ------
      integer :: nvars,i,ihit,nwords,start(256)
      character*999 :: strbuf,fname,fmt,var_name
      type(variable), allocatable :: vars(:)
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
        !Set variable fmt if specified
        if(nwords>=2) then
          read(strbuf(start(2):),*) fmt
          call set_fmt(vars(i),fmt)
        endif
      enddo

      !Get filename 
      call read_control_data(ctrl,fname_tag,fname)
      fname = adjustl(trim(fname))

      !Setup file
      call init_output(of,fname,char(9),vars)

      deallocate(vars)
      end subroutine

      end module
