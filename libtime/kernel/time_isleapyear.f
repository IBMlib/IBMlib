c this is <time_isleapyear.f> (extracted from ../libtime.f)
c automatically generated by "SPLITF.PL   V1.0   SPLIT Fortran source code"
c----------------------------------------------------------------------
c $Id: time_isleapyear.f,v 2.0 2000/08/05 21:58:30 thof Stab $
c
c 05/08/2000 by Thomas Forbriger (IfG Stuttgart)
c
c check for leap year
c
c REVISIONS and CHANGES
c    05/08/2000   V2.0   Thomas Forbriger
c
c ============================================================================
cS
      logical function time_isleapyear(year)
c
c is true if year is a leap-year (else false ;-))
c
c input:
c   year:   full qualified year value to be checked for being a leap-year
c
c last change: V2.00 (05/08/2000)
c
      integer year
cE
      integer iyear
      logical result
c 
      iyear=year
      call time_fullyear(iyear)
      result=(((mod(iyear,4).eq.0).and.(mod(iyear,100).ne.0)).or.
     &     (mod(iyear,400).eq.0)) 
      time_isleapyear=result
      return
      end
c
c ----- END OF <time_isleapyear.f> -----