/* this is <ctime_util_warning.c>
 * ----------------------------------------------------------------------------
 *
 * $Id: ctime_util_warning.c,v 2.3 2004/02/07 17:38:17 tforb Exp $
 *
 * 06/08/2000 by Thomas Forbriger (IfG Stuttgart)
 *
 * C language specific warning handler
 *
 * REVISIONS and CHANGES
 *    06/08/2000   V1.0   Thomas Forbriger
 *    09/08/2000   V1.1   both are WARNINGs not ERRORs
 *
 * ============================================================================
 */

#include <libtime.h>
#include <stdio.h>

/* 
 * Fortran calling convention:
 */
int time_util_warning__(char *caller, char *text, 
                        ftnlen caller_len, ftnlen text_len)
{
    int i;
    fputs("WARNING (", stderr);
    for (i=0; i<caller_len; i++) { fputc((int)*(caller++), stderr); }
    fputs("): ",stderr);
    for (i=0; i<text_len; i++) { fputc((int)*(text++), stderr); }
    fputc('\n', stderr);
    return 0;
} /* time_util_warning__ */

int time_util_warning_n__(char *caller, char *text, integer *n, 
                          ftnlen caller_len, ftnlen text_len)
{
    long int i;
    fputs("WARNING (", stderr);
    for (i=0; i<caller_len; i++) { fputc((int)*(caller++), stderr); }
    fputs("): ",stderr);
    for (i=0; i<text_len; i++) { fputc((int)*(text++), stderr); }
    i=(int)*n;
    fprintf(stderr," %d\n",i);
    return 0;
} /* time_util_warning_n__ */
 
/* ----- END OF ctime_util_warning.c ----- */
