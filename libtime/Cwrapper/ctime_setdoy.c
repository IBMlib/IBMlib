/* this is <ctime_setdoy.c>
 * ----------------------------------------------------------------------------
 *
 * $Id: ctime_setdoy.c,v 2.0 2000/08/06 12:49:10 thof Stab $
 *
 * 06/08/2000 by Thomas Forbriger (IfG Stuttgart)
 *
 * C-wrapper for FORTRAN time_setdoy libtime kernel function
 *
 * REVISIONS and CHANGES
 *    06/08/2000   V1.0   Thomas Forbriger
 *
 * ============================================================================
 */

#include <libtime.h>

void time_setdoy(long int day, long int month, time_Ts *Pdate)
{
  time_Tu *u;
  extern int time_setdoy__();
  u=(time_Tu *)Pdate;
  time_setdoy__(&day, &month, &u->farray);
} /* time_setdoy */
 
/* ----- END OF ctime_setdoy.c ----- */
