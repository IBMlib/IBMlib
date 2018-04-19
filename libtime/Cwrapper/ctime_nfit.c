/* this is <ctime_nfit.c>
 * ----------------------------------------------------------------------------
 *
 * $Id: ctime_nfit.c,v 2.0 2000/08/06 12:49:10 thof Stab $
 *
 * 06/08/2000 by Thomas Forbriger (IfG Stuttgart)
 *
 * C-wrapper for FORTRAN time_nfit libtime kernel function
 *
 * REVISIONS and CHANGES
 *    06/08/2000   V1.0   Thomas Forbriger
 *
 * ============================================================================
 */

#include <libtime.h>

void time_nfit(time_Ts Date1, time_Ts Date2, long int *n, time_Ts *Pfull)
{
  time_Tu *u1;
  time_Tu *u2;
  time_Tu *full;
  extern int time_nfit__();
  u1=(time_Tu *)&Date1;
  u2=(time_Tu *)&Date2;
  full=(time_Tu *)Pfull;
  time_nfit__(&u1->farray, &u2->farray, n, &full->farray);
} /* time_nfit */
 
/* ----- END OF ctime_nfit.c ----- */
