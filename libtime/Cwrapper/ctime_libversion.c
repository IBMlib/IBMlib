/* this is <ctime_libversion.c>
 * ----------------------------------------------------------------------------
 *
 * $Id: ctime_libversion.c,v 2.0 2000/08/06 12:49:10 thof Stab $
 *
 * 06/08/2000 by Thomas Forbriger (IfG Stuttgart)
 *
 * C-wrapper for FORTRAN time_libversion libtime kernel function
 *
 * REVISIONS and CHANGES
 *    06/08/2000   V1.0   Thomas Forbriger
 *
 * ============================================================================
 */

#include <libtime.h>

double time_libversion()
{
  extern doublereal time_libversion__();
  return((double)time_libversion__());
} /* time_libversion */
 
/* ----- END OF ctime_libversion.c ----- */
