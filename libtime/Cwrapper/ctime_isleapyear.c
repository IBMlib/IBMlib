/* this is <ctime_isleapyear.c>
 * ----------------------------------------------------------------------------
 *
 * $Id: ctime_isleapyear.c,v 2.0 2000/08/06 12:49:10 thof Stab $
 *
 * 06/08/2000 by Thomas Forbriger (IfG Stuttgart)
 *
 * C-wrapper for FORTRAN time_isleapyear libtime kernel function
 *
 * REVISIONS and CHANGES
 *    06/08/2000   V1.0   Thomas Forbriger
 *
 * ============================================================================
 */

#include <libtime.h>

long int time_isleapyear(long int year)
{
  extern logical time_isleapyear__(); 
  return ((long int)time_isleapyear__((integer *)&year));
} /* time_isleapyear */
 
/* ----- END OF ctime_isleapyear.c ----- */
