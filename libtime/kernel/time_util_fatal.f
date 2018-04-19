c this is <time_util_fatal.f>
c------------------------------------------------------------------------------
c $Id: time_util_fatal.f,v 2.1 2000/08/06 12:47:50 thof Exp $
c
c 05/08/2000 by Thomas Forbriger (IfG Stuttgart)
c
c Print a fatal error message using FORTRAN i/o routines and abort
c
c REVISIONS and CHANGES
c    05/08/2000   V1.0   Thomas Forbriger
c
cS
c==============================================================================
c
      subroutine time_util_fatal(caller,text)
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
      call abort()
c
      return
   50 format('ERROR (',a,'): ',a)
      end
c
c ----- END OF time_util_fatal.f -----
