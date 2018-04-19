/* this is <ctime_clear.c>
 * ----------------------------------------------------------------------------
 *
 * $Id: ctime_clear.c,v 2.0 2000/08/06 12:49:10 thof Stab $
 *
 * 06/08/2000 by Thomas Forbriger (IfG Stuttgart)
 *
 * C-wrapper for FORTRAN time_clear libtime kernel function
 *
 * REVISIONS and CHANGES
 *    06/08/2000   V1.0   Thomas Forbriger
 *
 * ============================================================================
 */

#include <libtime.h>

void time_clear(time_Ts* Pdate)
{
  time_Tu *u;
  extern int time_clear__();
  u=(time_Tu *)Pdate;
  time_clear__(&u->farray);
} /* time_clear */
 
/* ----- END OF ctime_clear.c ----- */
