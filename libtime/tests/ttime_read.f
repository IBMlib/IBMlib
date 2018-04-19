c this is <ttime_read.f>
c------------------------------------------------------------------------------
c
c 02/02/99 by Thomas Forbriger (IfG Stuttgart)
c
c test libtime read-function
c
c REVISIONS and CHANGES
c    02/02/99   V1.0   Thomas Forbriger
c
c==============================================================================
c
      program ttime_read
      character*79 string
      integer date(7)
      string=' '
      print *,'enter ''end'' to exit'
      do while (string(1:4).ne.'end ')
        print *,' '
        print *,'enter timestring: '
        read(5, '(a)') string
        if (string(1:4).ne.'end ') then
          print *,string
          call time_read(string, date)
          call time_sprint(date,string)
          print *,string
        endif
      enddo
      stop
      end
c
c ----- END OF ttime_read.f -----
