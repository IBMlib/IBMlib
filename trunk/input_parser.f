
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ---------------------------------------------------
c                          Input parser
c     ---------------------------------------------------
c     This module provide utilities for input of limited control data with simple mark-up
c
c     Version stamp  :  23 Jan 2014 15:34:37
c     Template file  :  /home/asbjorn/ModellingPackages/Atlantis/hydrodynamics_from_IBMlib/trunk/build_tools/input_parser_template.py
c
c     Implements a simple "name = value" input format like:
c    
c     ! --- Simulation parameters for ... 
c     y = 3445      ! an optional trailing comment
c     x = 45.
c     zvector = 1 56 99
c     aname   = whatever/dev/null
c
c     separator/comment tokens are optional single-character symbols, defaulting to "=" and "!"
c     Line order, white spaces does not matter.
c     Reading is faciliated by these 3 subroutines
c     
c     subroutine open_control_file(filename, ctrlfile)
c     subroutine read_control_data(ctrlfile, name, variable)
c     subroutine close_control_file(ctrlfile)
c
c     read_control_data is overloaded in "variable", to use same
c     subroutine call, independent of datatype of "variable"
c
c     comments:        everything after (and including) "!" (or whatever commnet-token) is skipped
c     empty lines:     just skipped
c     malformed lines: just skipped - later pass a return code
c     entry not found: currently fatal - later pass a return code
c
c     usage example:
c         use input_parser
c         real               :: rvalue
c         integer            :: ivalue, ivector(3)
c         character*45       :: string
c         type(control_file) :: ctrlfile
c         call open_control_file("inputexample", ctrlfile)
c         call read_control_data(ctrlfile, "x",rvalue)
c         call read_control_data(ctrlfile, "y",ivalue)
c         call read_control_data(ctrlfile, "zvector",ivector)
c         call read_control_data(ctrlfile, "aname", string)
c         call close_control_file(ctrlfile)
c
c     auxillary processing tools:
c         integer function count_tags(ctrlfile, variablename)
c
c     log:
c          * 18jan2005    1) added count_tags function 
c                         2) added optional read offset (for sequential read)
c          * 23mar2006    drop trailing "!" when reading string (and other) values
c          * 07jan2014    optional separator/comment single-character symbols
c          * 07jan2014    fix parsing irregularity ("= xxx" was previously considered valid markup)
c                         and segregate out line parsing to separate subroutine
c        
c     SVN tags:
c           $Rev$
c           $LastChangedDate$
c           $LastChangedBy$
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module input_parser
      implicit none

      character, parameter        :: default_separator_token = "="
      character, parameter        :: default_comment_token   = "!"
      integer, parameter, private :: linelen                 = 256

      type control_file
         character(len=linelen)          :: filename
         character(len=linelen), pointer :: rawdata(:)  ! lines
         character                       :: separator_token
         character                       :: comment_token
      end type

      private :: find_free_IO_unit     

      interface read_control_data
        module procedure read_control_data_char
        module procedure read_control_data_int
        module procedure read_control_data_real
        module procedure read_control_data_real8
        module procedure read_control_data_log
        module procedure read_control_data_intvec
        module procedure read_control_data_realvec
        module procedure read_control_data_real8vec
      end interface



      contains 

      subroutine open_control_file(filename, filehandler,
     +                             separator_token, comment_token)
c=========================================================================================
c     This subroutine initializes a control_file
c  
c     Allocate space to rawdata (stored line wise) and read it
c     data content in file without parsing it
c     close the file.
c     separator_token and comment_token are optional arguments allowing to parse dialects
c     If not present, theyt default to default_separator_token and default_comment_token
c     respectively
c=========================================================================================
      character*(*), intent(in)       :: filename
      type(control_file), intent(out) :: filehandler 
      character, optional, intent(in) :: separator_token
      character, optional, intent(in) :: comment_token 
      integer                         :: iunit, nlines, OK
      character(len=linelen)          :: linebuf
c-------------------------------------------
      filehandler%filename = filename
      call find_free_IO_unit(iunit)
      open(iunit, FILE=filename, ACTION='READ', IOSTAT=OK, PAD='YES')
      if (OK /= 0) then
         write(*,*) "open_control_file: error opening file "//filename
         stop
      endif
      
c.....inquire the length of the input file (nchars)
      nlines = 0
      do 
         read(iunit,FMT='(A)',end=345) linebuf
         nlines = nlines + 1 
      enddo 
 345  continue
      rewind(iunit)
      
c.....allocate and read character buffer for file content   
c.....read format is repeated as needed
      allocate(filehandler%rawdata(nlines))
      read(iunit, FMT='(A)') filehandler%rawdata
      close(iunit)

c.....Set dialect
      if (present(separator_token)) then
         filehandler%separator_token = separator_token
      else
         filehandler%separator_token = default_separator_token 
      endif
      if (present(comment_token)) then
         filehandler%comment_token = comment_token
      else
         filehandler%comment_token = default_comment_token 
      endif

      
      end subroutine open_control_file



      subroutine close_control_file(filehandler)
c========================================================
c     Shut down control_file
c     Have this subroutine to allow future developments
c========================================================
      type(control_file), intent(inout) :: filehandler
c-------------------------------------------
      deallocate(filehandler%rawdata)
      end subroutine close_control_file



      subroutine find_free_IO_unit(isfree) 
c========================================================
c     Grap first unused file unit
c========================================================
      integer, intent(out) :: isfree
      logical              :: is_opened
c--------------------------------------------------------
      do isfree=7, 99, 1
         inquire(UNIT=isfree, OPENED=is_opened)
         if (.not.is_opened) return
      enddo
      end subroutine find_free_IO_unit


      subroutine parse_line(filehandler,lineno,valid,isep,iend,txt)
c========================================================
c     This subroutine will parse filehandler%rawdata(lineno)
c     The line is returned in txt (with leading blanks from
c     filehandler%rawdata(lineno) removed). 
c     The result of parsing is returned in valid:
c        valid == 0 : a syntactical valid line is contained in txt
c                     txt(1      : isep-1) contains keyword
c                     txt(isep+1 : iend  ) contains value
c        valid == 1 : lineno > size(filehandler%rawdata), line not parsed
c        valid == 2 : filehandler%rawdata(lineno) is not syntactical valid
c        For valid > 0 isep, iend is arbitrary
c========================================================
      type(control_file), intent(in)      :: filehandler
      integer, intent(in)                 :: lineno
      integer, intent(out)                :: valid, isep, iend
      character(len=linelen), intent(out) :: txt
      integer                             :: icomm
c--------------------------------------------------------
      if (lineno > size(filehandler%rawdata)) then
         isep  = 0
         iend  = 0
         valid = 1
         return
      endif
      
      txt   = adjustl(filehandler%rawdata(lineno))
      isep  = index(txt, filehandler%separator_token)
      icomm = index(txt, filehandler%comment_token)
      if (icomm == 0) then
         iend = linelen    ! no comment, no need to cut a comment
      else
         iend = icomm - 1  ! a coment is present; skip trailing ampersand
      endif
      
      if     (isep .le. 1) then         ! skip comments/blank/malformed lines
          valid = 2
          return
      elseif (iend .le. isep) then
          valid = 2
          return
      endif
      
c.....txt(1:isep-1) txt(isep+1:iend) are valid slices 
      valid = 0       ! signal syntactical valid line 
       
      end subroutine parse_line





      integer function count_tags(filehandler, name)       
c========================================================
c     Count the number of times tag appears (possibly zero)
c     in the data set hold by filehandler 
c========================================================
      type(control_file),intent(in) :: filehandler         
      character*(*),intent(in)      :: name                                    
      integer                       :: iline, isep, iend, valid
      character(len=linelen)        :: thisline, thisname
c-------------------------------------------
      count_tags = 0
      do iline = 1, size(filehandler%rawdata)

         call parse_line(filehandler,iline,valid,isep,iend,thisline)
         if (valid>0) cycle
               
         read(thisline(1:isep-1),'(A)') thisname

c........drop leading/trailing blanks
         if (trim(thisname) == trim(adjustl(name))) then
             count_tags = count_tags + 1
         endif
      enddo

      end function count_tags






      subroutine read_control_data_char(filehandler,
     +                                     name, value, next)       ! here goes tag
c========================================================
c     Read a specific item name from filehandler
c     into value
c
c     Optional argument next controls whether to read starting from a
c     specified line (default: start from the beginning). Allows reading
c     sequential values of a tag by storing returned value of next, because
c     the matched line is returned in next
c
c     Implementation notes:
c         string function index(a,b) returns 0 if b is not in a    
c========================================================
      type(control_file)    :: filehandler         
      character*(*)         :: name
      integer, optional     :: next
      character(len=*)    :: value    ! here goes decl
      integer               :: iline, istart, isep, iend, valid
      character(len=linelen) :: thisline, thisname, cbuf
c-------------------------------------------

c.....resolve which line to start from

      if (present(next)) then
          istart = next
      else
          istart = 1       ! default
      endif


c.....loop over lines in control file

      do iline = istart, size(filehandler%rawdata)

         call parse_line(filehandler,iline,valid,isep,iend,thisline)
         if (valid>0) cycle
                
         read(thisline(1:isep-1),'(A)') thisname

c........drop leading/trailing blanks
         if (trim(thisname) == trim(adjustl(name))) then
             cbuf = adjustl(thisline(isep+1:iend))    ! we need to remove leading blanks for string reading
             read(cbuf, '(A)') value   ! here goes readfmt
             if (present(next)) next = iline                   ! store matched line number
             return
         endif  
      enddo
c.....not being returned by here means name not found
c.....maybe a success flag should be passed 
      write(*,*) "read_control_data: item "//trim(name)// 
     &           " not in "//filehandler%filename
      stop
      end subroutine read_control_data_char



      subroutine read_control_data_int(filehandler,
     +                                     name, value, next)       ! here goes tag
c========================================================
c     Read a specific item name from filehandler
c     into value
c
c     Optional argument next controls whether to read starting from a
c     specified line (default: start from the beginning). Allows reading
c     sequential values of a tag by storing returned value of next, because
c     the matched line is returned in next
c
c     Implementation notes:
c         string function index(a,b) returns 0 if b is not in a    
c========================================================
      type(control_file)    :: filehandler         
      character*(*)         :: name
      integer, optional     :: next
      integer             :: value    ! here goes decl
      integer               :: iline, istart, isep, iend, valid
      character(len=linelen) :: thisline, thisname, cbuf
c-------------------------------------------

c.....resolve which line to start from

      if (present(next)) then
          istart = next
      else
          istart = 1       ! default
      endif


c.....loop over lines in control file

      do iline = istart, size(filehandler%rawdata)

         call parse_line(filehandler,iline,valid,isep,iend,thisline)
         if (valid>0) cycle
                
         read(thisline(1:isep-1),'(A)') thisname

c........drop leading/trailing blanks
         if (trim(thisname) == trim(adjustl(name))) then
             cbuf = adjustl(thisline(isep+1:iend))    ! we need to remove leading blanks for string reading
             read(cbuf, *) value   ! here goes readfmt
             if (present(next)) next = iline                   ! store matched line number
             return
         endif  
      enddo
c.....not being returned by here means name not found
c.....maybe a success flag should be passed 
      write(*,*) "read_control_data: item "//trim(name)// 
     &           " not in "//filehandler%filename
      stop
      end subroutine read_control_data_int



      subroutine read_control_data_real(filehandler,
     +                                     name, value, next)       ! here goes tag
c========================================================
c     Read a specific item name from filehandler
c     into value
c
c     Optional argument next controls whether to read starting from a
c     specified line (default: start from the beginning). Allows reading
c     sequential values of a tag by storing returned value of next, because
c     the matched line is returned in next
c
c     Implementation notes:
c         string function index(a,b) returns 0 if b is not in a    
c========================================================
      type(control_file)    :: filehandler         
      character*(*)         :: name
      integer, optional     :: next
      real                :: value    ! here goes decl
      integer               :: iline, istart, isep, iend, valid
      character(len=linelen) :: thisline, thisname, cbuf
c-------------------------------------------

c.....resolve which line to start from

      if (present(next)) then
          istart = next
      else
          istart = 1       ! default
      endif


c.....loop over lines in control file

      do iline = istart, size(filehandler%rawdata)

         call parse_line(filehandler,iline,valid,isep,iend,thisline)
         if (valid>0) cycle
                
         read(thisline(1:isep-1),'(A)') thisname

c........drop leading/trailing blanks
         if (trim(thisname) == trim(adjustl(name))) then
             cbuf = adjustl(thisline(isep+1:iend))    ! we need to remove leading blanks for string reading
             read(cbuf, *) value   ! here goes readfmt
             if (present(next)) next = iline                   ! store matched line number
             return
         endif  
      enddo
c.....not being returned by here means name not found
c.....maybe a success flag should be passed 
      write(*,*) "read_control_data: item "//trim(name)// 
     &           " not in "//filehandler%filename
      stop
      end subroutine read_control_data_real



      subroutine read_control_data_real8(filehandler,
     +                                     name, value, next)       ! here goes tag
c========================================================
c     Read a specific item name from filehandler
c     into value
c
c     Optional argument next controls whether to read starting from a
c     specified line (default: start from the beginning). Allows reading
c     sequential values of a tag by storing returned value of next, because
c     the matched line is returned in next
c
c     Implementation notes:
c         string function index(a,b) returns 0 if b is not in a    
c========================================================
      type(control_file)    :: filehandler         
      character*(*)         :: name
      integer, optional     :: next
      real(kind=8)        :: value    ! here goes decl
      integer               :: iline, istart, isep, iend, valid
      character(len=linelen) :: thisline, thisname, cbuf
c-------------------------------------------

c.....resolve which line to start from

      if (present(next)) then
          istart = next
      else
          istart = 1       ! default
      endif


c.....loop over lines in control file

      do iline = istart, size(filehandler%rawdata)

         call parse_line(filehandler,iline,valid,isep,iend,thisline)
         if (valid>0) cycle
                
         read(thisline(1:isep-1),'(A)') thisname

c........drop leading/trailing blanks
         if (trim(thisname) == trim(adjustl(name))) then
             cbuf = adjustl(thisline(isep+1:iend))    ! we need to remove leading blanks for string reading
             read(cbuf, *) value   ! here goes readfmt
             if (present(next)) next = iline                   ! store matched line number
             return
         endif  
      enddo
c.....not being returned by here means name not found
c.....maybe a success flag should be passed 
      write(*,*) "read_control_data: item "//trim(name)// 
     &           " not in "//filehandler%filename
      stop
      end subroutine read_control_data_real8



      subroutine read_control_data_log(filehandler,
     +                                     name, value, next)       ! here goes tag
c========================================================
c     Read a specific item name from filehandler
c     into value
c
c     Optional argument next controls whether to read starting from a
c     specified line (default: start from the beginning). Allows reading
c     sequential values of a tag by storing returned value of next, because
c     the matched line is returned in next
c
c     Implementation notes:
c         string function index(a,b) returns 0 if b is not in a    
c========================================================
      type(control_file)    :: filehandler         
      character*(*)         :: name
      integer, optional     :: next
      logical             :: value    ! here goes decl
      integer               :: iline, istart, isep, iend, valid
      character(len=linelen) :: thisline, thisname, cbuf
c-------------------------------------------

c.....resolve which line to start from

      if (present(next)) then
          istart = next
      else
          istart = 1       ! default
      endif


c.....loop over lines in control file

      do iline = istart, size(filehandler%rawdata)

         call parse_line(filehandler,iline,valid,isep,iend,thisline)
         if (valid>0) cycle
                
         read(thisline(1:isep-1),'(A)') thisname

c........drop leading/trailing blanks
         if (trim(thisname) == trim(adjustl(name))) then
             cbuf = adjustl(thisline(isep+1:iend))    ! we need to remove leading blanks for string reading
             read(cbuf, *) value   ! here goes readfmt
             if (present(next)) next = iline                   ! store matched line number
             return
         endif  
      enddo
c.....not being returned by here means name not found
c.....maybe a success flag should be passed 
      write(*,*) "read_control_data: item "//trim(name)// 
     &           " not in "//filehandler%filename
      stop
      end subroutine read_control_data_log



      subroutine read_control_data_intvec(filehandler,
     +                                     name, value, next)       ! here goes tag
c========================================================
c     Read a specific item name from filehandler
c     into value
c
c     Optional argument next controls whether to read starting from a
c     specified line (default: start from the beginning). Allows reading
c     sequential values of a tag by storing returned value of next, because
c     the matched line is returned in next
c
c     Implementation notes:
c         string function index(a,b) returns 0 if b is not in a    
c========================================================
      type(control_file)    :: filehandler         
      character*(*)         :: name
      integer, optional     :: next
      integer             :: value(:)    ! here goes decl
      integer               :: iline, istart, isep, iend, valid
      character(len=linelen) :: thisline, thisname, cbuf
c-------------------------------------------

c.....resolve which line to start from

      if (present(next)) then
          istart = next
      else
          istart = 1       ! default
      endif


c.....loop over lines in control file

      do iline = istart, size(filehandler%rawdata)

         call parse_line(filehandler,iline,valid,isep,iend,thisline)
         if (valid>0) cycle
                
         read(thisline(1:isep-1),'(A)') thisname

c........drop leading/trailing blanks
         if (trim(thisname) == trim(adjustl(name))) then
             cbuf = adjustl(thisline(isep+1:iend))    ! we need to remove leading blanks for string reading
             read(cbuf, *) value   ! here goes readfmt
             if (present(next)) next = iline                   ! store matched line number
             return
         endif  
      enddo
c.....not being returned by here means name not found
c.....maybe a success flag should be passed 
      write(*,*) "read_control_data: item "//trim(name)// 
     &           " not in "//filehandler%filename
      stop
      end subroutine read_control_data_intvec



      subroutine read_control_data_realvec(filehandler,
     +                                     name, value, next)       ! here goes tag
c========================================================
c     Read a specific item name from filehandler
c     into value
c
c     Optional argument next controls whether to read starting from a
c     specified line (default: start from the beginning). Allows reading
c     sequential values of a tag by storing returned value of next, because
c     the matched line is returned in next
c
c     Implementation notes:
c         string function index(a,b) returns 0 if b is not in a    
c========================================================
      type(control_file)    :: filehandler         
      character*(*)         :: name
      integer, optional     :: next
      real                :: value(:)    ! here goes decl
      integer               :: iline, istart, isep, iend, valid
      character(len=linelen) :: thisline, thisname, cbuf
c-------------------------------------------

c.....resolve which line to start from

      if (present(next)) then
          istart = next
      else
          istart = 1       ! default
      endif


c.....loop over lines in control file

      do iline = istart, size(filehandler%rawdata)

         call parse_line(filehandler,iline,valid,isep,iend,thisline)
         if (valid>0) cycle
                
         read(thisline(1:isep-1),'(A)') thisname

c........drop leading/trailing blanks
         if (trim(thisname) == trim(adjustl(name))) then
             cbuf = adjustl(thisline(isep+1:iend))    ! we need to remove leading blanks for string reading
             read(cbuf, *) value   ! here goes readfmt
             if (present(next)) next = iline                   ! store matched line number
             return
         endif  
      enddo
c.....not being returned by here means name not found
c.....maybe a success flag should be passed 
      write(*,*) "read_control_data: item "//trim(name)// 
     &           " not in "//filehandler%filename
      stop
      end subroutine read_control_data_realvec



      subroutine read_control_data_real8vec(filehandler,
     +                                     name, value, next)       ! here goes tag
c========================================================
c     Read a specific item name from filehandler
c     into value
c
c     Optional argument next controls whether to read starting from a
c     specified line (default: start from the beginning). Allows reading
c     sequential values of a tag by storing returned value of next, because
c     the matched line is returned in next
c
c     Implementation notes:
c         string function index(a,b) returns 0 if b is not in a    
c========================================================
      type(control_file)    :: filehandler         
      character*(*)         :: name
      integer, optional     :: next
      real(kind=8)        :: value(:)    ! here goes decl
      integer               :: iline, istart, isep, iend, valid
      character(len=linelen) :: thisline, thisname, cbuf
c-------------------------------------------

c.....resolve which line to start from

      if (present(next)) then
          istart = next
      else
          istart = 1       ! default
      endif


c.....loop over lines in control file

      do iline = istart, size(filehandler%rawdata)

         call parse_line(filehandler,iline,valid,isep,iend,thisline)
         if (valid>0) cycle
                
         read(thisline(1:isep-1),'(A)') thisname

c........drop leading/trailing blanks
         if (trim(thisname) == trim(adjustl(name))) then
             cbuf = adjustl(thisline(isep+1:iend))    ! we need to remove leading blanks for string reading
             read(cbuf, *) value   ! here goes readfmt
             if (present(next)) next = iline                   ! store matched line number
             return
         endif  
      enddo
c.....not being returned by here means name not found
c.....maybe a success flag should be passed 
      write(*,*) "read_control_data: item "//trim(name)// 
     &           " not in "//filehandler%filename
      stop
      end subroutine read_control_data_real8vec


      end module
