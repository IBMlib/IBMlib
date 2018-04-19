c this is <time_util_warning.f>
c------------------------------------------------------------------------------
c $Id: time_util_warning.f,v 2.1 2000/08/06 14:08:47 thof Exp $
c
c 05/08/2000 by Thomas Forbriger (IfG Stuttgart)
c
c Print a warning message using FORTRAN i/o routines
c
c REVISIONS and CHANGES
c    05/08/2000   V1.0   Thomas Forbriger
c
cS
c==============================================================================
c
      subroutine time_util_warning(caller,text)
c
c declare parameters
      character*(*) caller,text
c
cE
c declare local variables
      integer index
c
c------------------------------------------------------------------------------
c go
      print 50,caller(1:index(caller,' ')-1),text
c
      return
   50 format('WARNING (',a,'): ',a)
      end
c
cS
c==============================================================================
c
      subroutine time_util_warning_n(caller,text,n)
c
c declare parameters
      character*(*) caller,text
      integer n
c
cE
c declare local variables
      integer index
c
c------------------------------------------------------------------------------
c go
      print 50,caller(1:index(caller,' ')-1),text,n
c
      return
   50 format('WARNING (',a,'): ',a,I6)
      end
c
c ----- END OF time_util_warning.f -----
