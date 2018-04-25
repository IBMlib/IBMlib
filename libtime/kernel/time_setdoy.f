c this is <time_setdoy.f> (extracted from ../libtime.f)
c automatically generated by "SPLITF.PL   V1.0   SPLIT Fortran source code"
c----------------------------------------------------------------------
c $Id: time_setdoy.f,v 2.1 2000/08/06 14:08:47 thof Exp $
c
c 05/08/2000 by Thomas Forbriger (IfG Stuttgart)
c
c set doy index from given day and month index
c
c REVISIONS and CHANGES
c    05/08/2000   V1.0   Thomas Forbriger
c                 V2.0   call language specific warning handler
c
c ============================================================================
cS
      subroutine time_setdoy(day, month, date)
c 
c Set doy in date from day and month (january first is doy 1).
c This routine will call time_fullyear first!
c 
c This routine is senseless in combination with relative times!
c
c input:
c   day:    day index within month
c   month:  month index within year (contained in date)
c input/output:
c   date:   doy index in time record is changed according to given day and
c           month index and year value in date
c
c last change: V2.00 (05/08/2000)
c
      integer day, month, date(7)
cE
      integer result
      integer sel, mon
      logical time_isleapyear
      integer days(12,2)
      data days/31,28,31,30,31,30,31,31,30,31,30,31,
     &          31,29,31,30,31,30,31,31,30,31,30,31/
c 
      if (date(1).eq.0) then
        call time_util_warning('time_setdoy',
     &                         'do not use this with relative times')
        call time_util_warning('time_setdoy',
     &                         'routine skipped...')
      else
        call time_fullyear(date(1))
        if (time_isleapyear(date(1))) then
          sel=2
        else
          sel=1
        endif
        result=0
        mon=1
    1   if (mon.ge.month) goto 2
          result=result+days(mon,sel)
          mon=mon+1
          if (mon.gt.13) call time_util_fatal('time_setdoy',
     &                           'month value out of range')
          goto 1
    2   continue
        result=result+day
        date(2)=result
      endif
      return
      end
c
c ----- END OF <time_setdoy.f> -----