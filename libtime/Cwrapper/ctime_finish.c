/* this is <ctime_finish.c>
 * ----------------------------------------------------------------------------
 *
 * $Id: ctime_finish.c,v 2.0 2000/08/06 12:49:10 thof Stab $
 *
 * 06/08/2000 by Thomas Forbriger (IfG Stuttgart)
 *
 * C-wrapper for FORTRAN time_finish libtime kernel function
 *
 * REVISIONS and CHANGES
 *    06/08/2000   V1.0   Thomas Forbriger
 *
 * ============================================================================
 */

#include <libtime.h>

void time_finish(time_Ts *Pdate)
{
  time_Tu *u;
  extern int time_finish__();
  u=(time_Tu *)Pdate;
  time_finish__(&u->farray);
} /* time_finish */
 
/* ----- END OF ctime_finish.c ----- */
