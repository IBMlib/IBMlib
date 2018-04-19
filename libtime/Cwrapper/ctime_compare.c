/* this is <ctime_compare.c>
 * ----------------------------------------------------------------------------
 *
 * $Id: ctime_compare.c,v 2.0 2000/08/06 12:49:10 thof Stab $
 *
 * 06/08/2000 by Thomas Forbriger (IfG Stuttgart)
 *
 * C-wrapper for FORTRAN time_compare libtime kernel function
 *
 * REVISIONS and CHANGES
 *    06/08/2000   V1.0   Thomas Forbriger
 *
 * ============================================================================
 */

#include <libtime.h>

long int time_compare(time_Ts Date1, time_Ts Date2)
{
  time_Tu *u1;
  time_Tu *u2;
  extern integer time_compare__();
  u1=(time_Tu *)&Date1;
  u2=(time_Tu *)&Date2;
  return((long int)time_compare__(&u1->farray, &u2->farray));
} /* time_compare */
 
/* ----- END OF ctime_compare.c ----- */
