/* this is <ctime_fullyear.c>
 * ----------------------------------------------------------------------------
 *
 * $Id: ctime_fullyear.c,v 2.0 2000/08/06 12:49:10 thof Stab $
 *
 * 06/08/2000 by Thomas Forbriger (IfG Stuttgart)
 *
 * C-wrapper for FORTRAN time_fullyear libtime kernel function
 *
 * REVISIONS and CHANGES
 *    06/08/2000   V1.0   Thomas Forbriger
 *
 * ============================================================================
 */

#include <libtime.h>

void time_fullyear(long int *year)
{
  extern int time_fullyear__();
  time_fullyear__(year);
} /* time_fullyear */
 
/* ----- END OF ctime_fullyear.c ----- */
