/* this is <ctime_sub.c>
 * ----------------------------------------------------------------------------
 *
 * $Id: ctime_sub.c,v 2.0 2000/08/06 12:49:10 thof Stab $
 *
 * 06/08/2000 by Thomas Forbriger (IfG Stuttgart)
 *
 * C-wrapper for FORTRAN time_sub libtime kernel function
 *
 * REVISIONS and CHANGES
 *    06/08/2000   V1.0   Thomas Forbriger
 *
 * ============================================================================
 */

#include <libtime.h>

void time_sub(time_Ts Date1, time_Ts Date2, time_Ts *Pdate3)
{
  time_Tu *u1;
  time_Tu *u2;
  time_Tu *u3;
  extern int time_sub__();
  u1=(time_Tu *)&Date1;
  u2=(time_Tu *)&Date2;
  u3=(time_Tu *)Pdate3;
  time_sub__(&u1->farray, &u2->farray, &u3->farray);
} /* time_sub */
 
/* ----- END OF ctime_sub.c ----- */
