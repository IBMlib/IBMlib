/* this is <ctime_util_fatal.c>
 * ----------------------------------------------------------------------------
 *
 * $Id: ctime_util_fatal.c,v 2.1 2000/08/06 14:08:47 thof Exp $
 *
 * 06/08/2000 by Thomas Forbriger (IfG Stuttgart)
 *
 * C language specific fatal error handler
 *
 * REVISIONS and CHANGES
 *    06/08/2000   V1.0   Thomas Forbriger
 *
 * ============================================================================
 */

#include <libtime.h>
#include <stdio.h>

/* 
 * Fortran calling convention:
 */
int time_util_fatal__(char *caller, char *text, 
                      ftnlen caller_len, ftnlen text_len)
{
    int i;
    fputs("ERROR (", stderr);
    for (i=0; i<caller_len; i++) { fputc((int)*(caller++), stderr); }
    fputs("): ",stderr);
    for (i=0; i<text_len; i++) { fputc((int)*(text++), stderr); }
    fputc('\n', stderr);
    abort();
} /* time_util_fatal__ */

 
/* ----- END OF ctime_util_fatal.c ----- */
