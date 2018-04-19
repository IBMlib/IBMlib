/* this is <ctime_read.c>
 * ----------------------------------------------------------------------------
 *
 * $Id: ctime_read.c,v 2.2 2000/08/09 21:02:51 thof Exp $
 *
 * 06/08/2000 by Thomas Forbriger (IfG Stuttgart)
 *
 * fill a time_Ts structure from a character string
 *
 * REVISIONS and CHANGES
 *    06/08/2000   V1.0   Thomas Forbriger
 *    09/08/2000   V1.1   correct handling of relative dates
 *
 * ============================================================================
 */

#include <libtime.h>
#include <string.h>
#include <stdio.h>

int time_read(time_Ts *Date, char *String)
{
  double sec;
  char *ptr;
  char cmilmicsec[7];
  long int day, month, milmicsec;
  int i;

  time_clear(Date);
  ptr=String;

  ptr=strpbrk(ptr, "0123456789");
  RETURNERROR((ptr==NULL), "time_read", \
    "year is missing", EXIT_FAILURE)
  Date->year=strtol(ptr, &ptr, 10);

  ptr=strpbrk(ptr, "0123456789");
  RETURNERROR((ptr==NULL), "time_read", \
    "month is missing", EXIT_FAILURE)
  month=strtol(ptr, &ptr, 10);

  ptr=strpbrk(ptr, "0123456789");
  RETURNERROR((ptr==NULL), "time_read", \
    "day is missing", EXIT_FAILURE)
  day=strtol(ptr, &ptr, 10);

  ptr=strpbrk(ptr, "0123456789");
  if (ptr!=NULL) {
    Date->hour=strtol(ptr, &ptr, 10);
    ptr=strpbrk(ptr, "0123456789");
  }
  if (ptr!=NULL) {
    Date->minute=strtol(ptr, &ptr, 10);
    ptr=strpbrk(ptr, "0123456789");
  }
  if (ptr!=NULL) {
    Date->second=strtol(ptr, &ptr, 10);
    ptr=strpbrk(ptr, "0123456789");
  }
  if (ptr!=NULL) {
    for(i=0; i<6; i++) {
      if (ptr!=NULL) {
        if (isdigit(ptr[i])) {
          cmilmicsec[i]=ptr[i];
        } else {
          ptr=NULL;
          cmilmicsec[i]='0';
        }
      } else {
        cmilmicsec[i]='0';
      }
      cmilmicsec[6]='\0';
    }
    milmicsec=strtol(cmilmicsec, &ptr, 10);
    Date->milsec=(long int)(milmicsec/1000);
    Date->micsec=milmicsec-(Date->milsec*1000);
  }

  if ((month>0 && Date->year>0) || month>0) 
  { 
    time_fullyear(&Date->year); 
    time_setdoy(day, month, Date);
  }
  else
  {
    Date->doy=day;
  }

  time_norm(Date);
  return(EXIT_SUCCESS);
} /* time_read */
 
/* ----- END OF ctime_read.c ----- */
