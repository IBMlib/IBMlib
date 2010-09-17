      program test
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ifort -c -I../.. vortex1.f  test.f; ifort test.o vortex1.o ../../constants.o ../../time_tools.o ../../libtime/libtime77.a
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use physical_fields
      implicit none
      real    :: xyz1(3),xyz2(3)
      logical :: anycross
      real    :: xyzref(3), xyzhit(3)
      
      xyz1 = (/3., 55., 1./)
      xyz2 = (/66., 653., 3./)
      
      call coast_line_intersection(xyz1, xyz2, anycross, xyzref, 
     +                                   xyzhit) 
      write(*,*) "anycross=",  anycross 
      write(*,*) "xyzref  =",  xyzref
      write(*,*) "xyzhit  =",  xyzhit
      end program
