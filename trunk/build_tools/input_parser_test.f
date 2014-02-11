      program test
c --------------------------------------------------------
c ifort -e90 input_parser.f test.f ; a.out
c
c option e90: issue errors for non-standard  Fortran  90
c --------------------------------------------------------
      use input_parser
      implicit none

      real               :: rvalue,  rvector(2)
      real(kind=8)       :: r8value, r8vector(2)
      integer            :: ivalue, ivector(3),next
      character*45       :: string
      type(control_file) :: ctrlfile
c     --------------------------------
      call open_control_file("input_parser_inputexample", ctrlfile)
      call read_control_data(ctrlfile, "x", rvalue, next)
      write(*,*) "variable x appears ", 
     +           count_tags(ctrlfile, "x"), "times"
      write(*,*) "x        =", rvalue
      write(*,*) "x appeared in line ", next
      call read_control_data(ctrlfile, "x", rvalue, next+1)
      write(*,*) "x (next appearence) =", rvalue
      call read_control_data(ctrlfile, "y", ivalue)
      write(*,*) "y                 =", ivalue
      call read_control_data(ctrlfile, "x", ivalue)
      write(*,*) "x as int          =", ivalue
      call read_control_data(ctrlfile, "reals", rvector)
      write(*,*) "reals             =", rvector 
      call read_control_data(ctrlfile, "reals", r8vector)
      write(*,*) "reals (as kind=8) =", r8vector
      call read_control_data(ctrlfile, "ints", ivector)
      write(*,*) "ints              =", ivector
      call read_control_data(ctrlfile, "aname", string)
      write(*,*) "aname             =", string
      call close_control_file(ctrlfile)
 
      end program

    
    
