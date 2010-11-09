ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Collection of misc run time utilities
c     organized as standard subroutine/function aggregation
c     (no module is defined)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine find_free_IO_unit(isfree) 
c========================================================
c     Grap first unused file unit
c========================================================
      integer :: isfree
      logical :: is_opened
      do isfree=7, 99, 1
         inquire(UNIT=isfree, OPENED=is_opened)
         if (.not.is_opened) return
      enddo
      end subroutine find_free_IO_unit

      subroutine abort_run(proc,message)
        character(len=*) proc, message
        write(*,* ) "HALT: "//trim(proc)//": "//trim(message)
        stop 
      end subroutine
