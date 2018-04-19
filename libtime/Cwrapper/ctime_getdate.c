/* this is <ctime_getdate.c>
 * ----------------------------------------------------------------------------
 *
 * $Id: ctime_getdate.c,v 2.0 2000/08/06 12:49:10 thof Stab $
 *
 * 06/08/2000 by Thomas Forbriger (IfG Stuttgart)
 *
 * C-wrapper for FORTRAN time_getdate libtime kernel function
 *
 * REVISIONS and CHANGES
 *    06/08/2000   V1.0   Thomas Forbriger
 *
 * ============================================================================
 */

#include <libtime.h>

void time_getdate(long int *day, long int *month, time_Ts Date)
{
  time_Tu *u;
  extern int time_getdate__();
  u=(time_Tu *)&Date;
  time_getdate__(day, month, &u->farray);
} /* time_getdate */
 
/* ----- END OF ctime_getdate.c ----- */
