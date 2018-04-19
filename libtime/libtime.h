/* this is <libtime.h>
 * ----------------------------------------------------------------------------
 * 
 * $Id: libtime.h,v 2.4 2000/08/09 21:02:51 thof Exp $
 *
 * 12/08/97 by Thomas Forbriger (IfG Stuttgart)
 *
 * some definitions and prototypes for the C libtime interface
 *
 * REVISIONS and CHANGES
 *    12/08/97   V1.0   Thomas Forbriger
 *    06/08/00   V2.0   do not include f2c - rather copy relevant definitions
 *                      added C linkage convention for C++ 
 *                      use namespace time_kernel
 *
 * ============================================================================
 */

#ifndef _TF_LIBTIME_H
#define _TF_LIBTIME_H

/* #include <f2c.h> */
/* #include <stdio.h> */

#ifdef __cplusplus
namespace time_kernel {
extern "C" {
#endif

/*
 * all f2c stuff that is needed here
 * =================================
 */

#ifndef F2C_INCLUDE
/* FORTRAN (f2c) types needed by the wrapper functions */
typedef long int integer;
typedef double doublereal;
typedef long int logical;
typedef long int ftnlen;
#endif

/*
 * a few tf-macros needed by the C specific functions
 * ==================================================
 */

/* return value of time_read on success */
#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif

/* return value of time_read on failure */
#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif

/* error handling macro used within time_read */
#ifndef RETURNERROR
#define RETURNERROR( EXPR , SUB, STR, CODE )\
  if ( EXPR ) { fprintf(stderr, "ERROR (%s):\n   %s\n", SUB, STR );\
      return(CODE); }
#endif

/*S*/
/*
 * some macro constants
 * ====================
 */

/* value returned by time_isleapyear in case year IS a leap-year */
#define TIME_ISLEAP (1)
/* value returned by time_isleapyear in case year IS NOT a leap-year */
#define TIME_ISNOTLEAP (0)

/* length of static string buffer in time_sprint */
#define TIME_SLEN (35)

/* 
 * time data structure
 * ===================
 */

/* standard structure to hold date record */
typedef struct {
  long int year;    /* year  (=0 for relative times)                */
  long int doy;     /* day within yaer (may be 0 for relative times */
  long int hour;    /* hour within day                              */
  long int minute;  /* minute within hour                           */
  long int second;  /* second within minute                         */
  long int milsec;  /* millisecond within second                    */
  long int micsec;  /* microsecond within millisecond               */
} time_Ts;

/*E*/

/* this unios is used to convert date records between different
 * representations and language specific functions */
typedef union {
  long int array[7];  /* linear access in C */
  integer farray[7];  /* linear access in FORTRAN */
  time_Ts s;          /* the structure that is fed in */
} time_Tu;

/*S*/
/*
 * wrapper function prototypes
 * ===========================
 */

/* void time_add(time_Ts date1, time_Ts date2, time_Ts *date3)
 *
 * date1:   input: any date record
 * date2:   input: date record (may be absolute if date1 is relative)
 * date3:   output: sum of date1 and date2
 */
void time_add(time_Ts, time_Ts, time_Ts *);

/* void time_clear(time_Ts *date)
 *
 * date:    input: any date record
 *          output: zero relative time
 */
void time_clear(time_Ts *);

/* long int time_compare(time_Ts date1, time_Ts date2)
 *
 * date1:   input: any date record
 * date2:   input: any date record (must be relative is date1 is relative,
 *                 must be absolute if date1 is absolute)
 * returns: 1  if date1 >  date2
 *          0  if date1 == date2
 *          -1 if date1 <  date2
 *          -2 when mixing relative and absolute date records
 */
long int time_compare(time_Ts, time_Ts);

/* void time_copy(time_Ts date1, time_Ts *date2)
 *
 * date1:   input: any date record
 * date2:   output: copy of date1
 */
void time_copy(time_Ts, time_Ts *);

/* void time_div(time_Ts date1, time_Ts *date2, long int n, long int *rest)
 *
 * date1:   input: any relative time
 * date2:   output: n-th fraction of date1
 * n:       input: divisor for date1
 * rest:    output: rest of division in mircoseconds
 *                  always: date1 >= (n*date2)
 */
void time_div(time_Ts, time_Ts *, long int, long int *);

/* void time_finish(time_Ts *date)
 *
 * date:    input: any date record
 *          output: fully qualified and regularized date record
 */
void time_finish(time_Ts *);

/* void time_fullyear(long int *year)
 *
 * year:    ainput: ny year value (may be a 2-digit abbreviation)
 *          output: a full qualified year value 
 */
void time_fullyear(long int *);

/* void time_getdate(long int *day, long int *month, time_Ts date)
 *
 * day:     output: day within month index of date
 * month:   output: month wihtin year index of date
 * date:    input: any absolute date record
 */
void time_getdate(long int*, long int*, time_Ts);

/* long int time_isleapyear(long int year)
 *
 * year:    input: full qualified year value to be checked
 * returns: TIME_ISLEAP       if argument is a leap-year
 *          TIME_ISNOLEAP     if argument is not a leap-year
 */
long int time_isleapyear(long int);

/* double time_libversion
 *
 * returns: version number of library kernel
 */
double time_libversion();

/* void time_mul(time_Ts date1, time_Ts *date2, long int n)
 *
 * date1:   input: any relative date record
 * date2:   output: n times date1
 * n:       input: factor to multiply date1 with
 */
void time_mul(time_Ts, time_Ts *, long int);

/* void time_nfit(time_Ts date1, time_Ts date2, long int *n, time_Ts *full)
 *
 * date1:   input: any relative time record
 * date2:   input: any relative time record
 * n:       output: number os date2 intervals the fit best into date1
 *                  so that abs((n*date2)-date1) <= date2/2
 * full:    output: full time span defined by n and date2 (full=n*date2)
 */
void time_nfit(time_Ts, time_Ts, long int *, time_Ts *);

/* void time_norm(time_Ts *date)
 *
 * date:    input: any date record
 *          output: regularized date record
 */
void time_norm(time_Ts *);

/* void time_setdoy(long int day, long int month, time_Ts *date)
 *
 * day:     input: day index within month
 * month:   input: month index within year
 * date:    input: any date record with year set
 *          output: has correct doy set from day and month 
 */
void time_setdoy(long int, long int, time_Ts *);

/* void time_sub(time_Ts date1, time_Ts date2, time_Ts date3)
 *
 * date1:   input: any date record
 * date2:   input: any date record
 * date3:   output: absolute (positive) difference between date1 and date2
 */
void time_sub(time_Ts, time_Ts, time_Ts *);

/*
 * prototypes of pure C functions
 * ==============================
 */

/* int time_read(time_Ts *date, char *string)
 *
 * string:  input: character representation of a time with the fields in the
 *                 following order:
 *                    year month day hour minute seconds
 *                 - the fields must be separated by non-numeric characters
 *                 - all fields except the field 'seconds' must be integer
 *                 - you may omit any number of trailing fields
 *                 - year AND month must be zero to specify a relative time
 * date:    output: full qualified and regularized date record specified by
 *                  string
 * returns: EXIT_SUCCESS on success
 *          EXIT_FAILURE on FAILURE
 *
 * NOTICE: time_read is not the direct inverse operation of time_sprint
 */
int time_read(time_Ts *, char *);

/* char *time_sprint(time_Ts date)
 *
 * date:    input: any date record
 * returns: a pointer to a static character string of at most TIME_SLEN
 *          characters length containing the ASCII text representation
 *          of date (NOTICE: the next call to time_sprint will overwrite
 *          the static character string)
 */
char *time_sprint(time_Ts);

/*E*/

#ifdef __cplusplus
}}
#endif

#endif /* _TF_LIBTIME_H */

/* ----- END OF libtime.h ----- */
