/* this is <ctime_sprint.c>
 * ----------------------------------------------------------------------------
 *
 * $Id: ctime_sprint.c,v 2.0 2000/08/06 12:49:10 thof Stab $
 *
 * 06/08/2000 by Thomas Forbriger (IfG Stuttgart)
 *
 * libtime C kernel function
 * 
 * NOTICE: This routine returns a pointer to a static character array. The
 * next call to time_sprint will destroy the contents of this string.
 *
 * REVISIONS and CHANGES
 *    06/08/2000   V1.0   Thomas Forbriger
 *
 * ============================================================================
 */

#include <libtime.h>
#include <stdio.h>

char *time_sprint(time_Ts Date)
{
  long int day;
  long int month;
  static char sdate[TIME_SLEN];

  time_norm(&Date);

  if (Date.year == 0L) {
    sprintf(sdate, "%03dd %02dh %02dm %02d.%03d%03ds", Date.doy,
      Date.hour, Date.minute, Date.second, Date.milsec, Date.micsec);
  } else {
    time_getdate(&day, &month, Date);
    sprintf(sdate, "%03d %02d/%02d/%04d %02d:%02d:%02d.%03d%03d",
      Date.doy, day, month, Date.year, 
      Date.hour, Date.minute, Date.second, Date.milsec, Date.micsec);
  }
  return(sdate);
} /* time_setdoy */
 
/* ----- END OF ctime_sprint.c ----- */
