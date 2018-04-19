c this is <time_norm.f> (extracted from ../libtime.f)
c automatically generated by "SPLITF.PL   V1.0   SPLIT Fortran source code"
c----------------------------------------------------------------------
c $Id: time_norm.f,v 2.0 2000/08/05 21:58:30 thof Stab $
c
c 05/08/2000 by Thomas Forbriger (IfG Stuttgart)
c
c regularize any time record (relative or absolute)
c
c REVISIONS and CHANGES
c    05/08/2000   V2.0   Thomas Forbriger
c
c ============================================================================
cS
      subroutine time_norm(date)
c
c Regularize any time record (relative or absolute)
c (i.e. set all fields to correct value range)
c
c input/output:
c   date:   date record to be regularized
c
c last change: V2.00 (05/08/2000)
c
      integer date(7)
cE
      integer carry, i, diy
      integer limit(7)
      logical time_isleapyear
      data limit/-1,-1,24,60,60,1000,1000/
c
c set linear value ranges (from hours on)
      do i=7,3,-1
        carry=int(date(i)/limit(i))
        if (date(i).lt.0) carry=carry-1
        date(i)=date(i)-carry*limit(i)
        date(i-1)=date(i-1)+carry
      enddo
c work on date if not relative
      if (date(1).gt.0) then
c set full 4 digit year
        call time_fullyear(date(1))
        if (date(2).gt.0) then
          if (time_isleapyear(date(1))) then
            diy=366
          else
            diy=365
          endif
    1     if (date(2).le.diy) goto 2
            date(2)=date(2)-diy
            date(1)=date(1)+1
            if (time_isleapyear(date(1))) then
              diy=366
            else
              diy=365
            endif
            goto 1
    2     continue
        else
    3     if (date(2).gt.0) goto 4
            date(1)=date(1)-1
            if (time_isleapyear(date(1))) then
              diy=366
            else
              diy=365
            endif
            date(2)=date(2)+diy
            goto 3
    4     continue
        endif
      endif
      return
      end
c
c ----- END OF <time_norm.f> -----
