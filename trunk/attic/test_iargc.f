      program test_iargc
c     ------------------------------------------------------------
c     Access portability of iargc() usage
c     ------------------------------------------------------------
c     ifort    -std90                    test_iargc.f    [OK]
c     gfortran -fall-intrinsics -std=f95 test_iargc.f    [OK]
c     gfortran -fall-intrinsics -std=gnu test_iargc.f    [OK]
c     pgf90    -Mstandard                test_iargc.f    [OK]
c     a.out yy ee rr
c     ------------------------------------------------------------
      implicit none
c      integer, external :: iargc
      integer :: iargc
      write(*,*) "#args=", iargc()
      end program
