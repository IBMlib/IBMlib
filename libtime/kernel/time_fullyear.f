c this is <time_fullyear.f> (extracted from ../libtime.f)
c automatically generated by "SPLITF.PL   V1.0   SPLIT Fortran source code"
c----------------------------------------------------------------------
c $Id: time_fullyear.f,v 2.1 2000/08/06 14:08:47 thof Exp $
c
c 05/08/100 by Thomas Forbriger (IfG Stuttgart)
c
c create a full qualified year value from a two-digit abbreviation
c
c REVISIONS and CHANGES
c    05/08/100   V1.0   Thomas Forbriger
c                V2.0   use language specific warning wrapper
c
c ============================================================================
cS
      subroutine time_fullyear(year)
c
c Makes year to be a full 4 digit year value. This is used to set
c the year to a meanigfull value. Notice that this routine does not
c make sense in combination with relative times.
c
c year < 70         ---> year := year+2000
c 69 < year < 100   ---> year := year+1900
c
c last change: V2.00 (05/08/2000)
c
      integer year
cE
      if (year.lt.70) year=year+2000
      if (year.lt.100) year=year+1900
c      if (year.lt.1970) 
c     &  call time_util_warning_n('time_fullyear',
c     &                           'spurious year value: ',year)
      return
      end
c
c ----- END OF <time_fullyear.f> -----
