/* this is <ctime_norm.c>
 * ----------------------------------------------------------------------------
 *
 * $Id: ctime_norm.c,v 2.0 2000/08/06 12:49:10 thof Stab $
 *
 * 06/08/2000 by Thomas Forbriger (IfG Stuttgart)
 *
 * C-wrapper for FORTRAN time_norm libtime kernel function
 *
 * REVISIONS and CHANGES
 *    06/08/2000   V1.0   Thomas Forbriger
 *
 * ============================================================================
 */

#include <libtime.h>

void time_norm(time_Ts *Pdate)
{
  time_Tu *u;
  extern int time_norm__();
  u=(time_Tu *)Pdate;
  time_norm__(&u->farray);
} /* time_norm */
 
/* ----- END OF ctime_norm.c ----- */
