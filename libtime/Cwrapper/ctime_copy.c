/* this is <ctime_copy.c>
 * ----------------------------------------------------------------------------
 *
 * $Id: ctime_copy.c,v 2.0 2000/08/06 12:49:10 thof Stab $
 *
 * 06/08/2000 by Thomas Forbriger (IfG Stuttgart)
 *
 * C-wrapper for FORTRAN time_copy libtime kernel function
 *
 * REVISIONS and CHANGES
 *    06/08/2000   V1.0   Thomas Forbriger
 *
 * ============================================================================
 */

#include <libtime.h>

void time_copy(time_Ts Date1, time_Ts *Pdate2)
{
  time_Tu *u1;
  time_Tu *u2;
  extern int time_copy__();
  u1=(time_Tu *)&Date1;
  u2=(time_Tu *)Pdate2;
  time_copy__(&u1->farray, &u2->farray);
} /* time_copy */
 
/* ----- END OF ctime_copy.c ----- */
