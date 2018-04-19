/* this is <libtime++.h>
 * ----------------------------------------------------------------------------
 *
 * $Id: libtime++.h,v 2.16 2005/07/15 15:36:54 tforb Exp $
 *
 * Copyright (c) by Thomas Forbriger (IfG Stuttgart)
 *
 * interface of class libtime++
 *
 * REVISIONS and CHANGES
 *    09/08/2000   V1.0   Thomas Forbriger
 *    10/08/2000   V1.1   started with operators
 *    06/09/2000   V1.2   finished operators and rearranged base class
 *                        constructors
 *    12/09/2000   V1.3   added absolute - absolute --> relative
 *                        implemented non-member multiplication
 *    22/12/2000   V1.4   changed namespace time to libtime (resolved conflict
 *                        with system time library)
 *    22/12/2003   V1.5   added function now()
 *    22/12/2003   V1.6   whooo - was a bit sloppy - corrected some errors
 *    27/12/2003   V1.7   added member function of TAbsoluteTime to
 *                        set doy explicitely
 *    13/01/2004   V1.8   added constructor from double for TRelativeTime
 *    02/02/2004   V1.9   constructor from double for TRelativeTime
 *                        was highly ambiguous
 *                        move conversion to doubel2time() function
 *    15/07/2005   V1.10  print usage information
 *
 * ============================================================================
 */
 
#ifndef TF_LIBTIME_H_
#define TF_LIBTIME_H_ \
  "$Id: libtime++.h,v 2.16 2005/07/15 15:36:54 tforb Exp $"
 
#include <string>
#include "libtime.h"

using std::string;

namespace libtime {

/*!
 * TBaseClassTime
 * ==============
 *
 * This is a class definition holding common parts od relative times and
 * absolute times. It's constructor is protected as it should not be
 * instantiated. It's an abstract class.
 */
class TBaseClassTime {
/*S*/
// member functions common to both flavours of time
  public:
    std::string timestring() const;
    long int hour() const;
    long int minute() const;
    long int second() const;
    long int milsec() const;
    long int micsec() const;
    double float_second() const;

    operator time_kernel::time_Ts() const;
    operator std::string() const;
/*E*/
// make constructors and member protected
  protected:
    TBaseClassTime(const std::string &Itimestring);
    TBaseClassTime(char *Itimestring);
    TBaseClassTime(const time_kernel::time_Ts &Itime_Ts);
    TBaseClassTime(const long int &year, const long int &doy,
                   const long int &hour=0,
                   const long int &minute=0, const long int &second=0,
                   const long int &milsec=0, const long int &micsec=0);
    void string_read(const std::string &timestring);
    void char_read(char *timestring);
    // well, there is no need to make norm accessible to public, as we will
    // always norm the structure ourselfs
    void norm();
    // here we hold the current values
    time_kernel::time_Ts Mtime_Ts;
};

/*S*/
  
//! class to contain absolute times
class TAbsoluteTime: public TBaseClassTime {
  friend class TRelativeTime;
  public:
    TAbsoluteTime(const std::string &Itimestring);
    TAbsoluteTime(char *Itimestring);
    TAbsoluteTime(const time_kernel::time_Ts &Itime_Ts);
    TAbsoluteTime(const long int &year=2000, const long int &month=1,
                  const long int &day=1, const long int &hour=0,
                  const long int &minute=0, const long int &second=0,
                  const long int &milsec=0, const long int &micsec=0);

    void setdoy(const long int &doy);
    void setdoy(const long int &day, const long int &month);
    void setdate(const long int &day, const long int &month);
    void getdate(long int &day, long int &month) const;
    bool isleapyear() const;

    long int year() const;
    long int doy() const;
    long int month() const;
    long int day() const;

    bool operator==(const TAbsoluteTime &A) const;
    bool operator!=(const TAbsoluteTime &A) const;
    bool operator<=(const TAbsoluteTime &A) const;
    bool operator>=(const TAbsoluteTime &A) const;
    bool operator< (const TAbsoluteTime &A) const;
    bool operator> (const TAbsoluteTime &A) const;

    TAbsoluteTime &operator+=(const TRelativeTime &A);
    TAbsoluteTime &operator-=(const TRelativeTime &A);
    TAbsoluteTime  operator+ (const TRelativeTime &A) const;
    TAbsoluteTime  operator- (const TRelativeTime &A) const;
    TRelativeTime  operator- (const TAbsoluteTime &A) const;

    TAbsoluteTime &operator= (const time_kernel::time_Ts &A);
    TAbsoluteTime &operator= (const std::string &timestring);
    TAbsoluteTime &operator= (char *timestring);
};

//! class to contain relative times
class TRelativeTime: public TBaseClassTime {
  friend class TAbsoluteTime;
  public:
    TRelativeTime(const std::string &Itimestring);
    TRelativeTime(char *Itimestring);
    TRelativeTime(const time_kernel::time_Ts &Itime_Ts);
    TRelativeTime(const int &days=0, const int &hour=0,
                  const int &minute=0, const int &second=0,
                  const int &milsec=0, const int &micsec=0);

    long int days() const;
    void nfit(const TRelativeTime &delta, long int &n, TRelativeTime &full)
      const;
    void div(const long int &n, TRelativeTime &frac, long int &rest)
      const;

    bool operator==(const TRelativeTime &A) const;
    bool operator!=(const TRelativeTime &A) const;
    bool operator<=(const TRelativeTime &A) const;
    bool operator>=(const TRelativeTime &A) const;
    bool operator< (const TRelativeTime &A) const;
    bool operator> (const TRelativeTime &A) const;

    TRelativeTime &operator+=(const TRelativeTime &A);
    TRelativeTime &operator-=(const TRelativeTime &A);
    TRelativeTime  operator+ (const TRelativeTime &A) const;
    TRelativeTime  operator- (const TRelativeTime &A) const;

    TAbsoluteTime  operator+ (TAbsoluteTime A) const;
    TAbsoluteTime  operator- (TAbsoluteTime A) const;

    TRelativeTime &operator*=(const long int &n);
    TRelativeTime &operator/=(const long int &n);
    TRelativeTime &operator%=(const long int &n);

    TRelativeTime  operator* (const long int &n) const;
    TRelativeTime  operator/ (const long int &n) const;
    TRelativeTime  operator% (const long int &n) const;

    long int       operator/ (const TRelativeTime &A) const;
    TRelativeTime  operator% (const TRelativeTime &A) const;

    TRelativeTime &operator= (const time_kernel::time_Ts &A);
    TRelativeTime &operator= (const std::string &timestring);
    TRelativeTime &operator= (char *timestring);
};

TRelativeTime operator* (const long int &n, const TRelativeTime &A);

/*E*/

/*
 * Here we go for inline functions
 * ===============================
 */

/*
 * TBaseClassTime
 * --------------
 */
inline TBaseClassTime::TBaseClassTime(const std::string &Itimestring)
  { string_read(Itimestring); }
inline TBaseClassTime::TBaseClassTime(char *Itimestring)
  { char_read(Itimestring); }

inline std::string TBaseClassTime::timestring() const
  { return(std::string(time_kernel::time_sprint(Mtime_Ts))); }

inline void TBaseClassTime::norm()
  { time_kernel::time_norm(&Mtime_Ts); }

inline TBaseClassTime::TBaseClassTime(const time_kernel::time_Ts &Itime_Ts):
  Mtime_Ts(Itime_Ts) { }

inline TBaseClassTime::TBaseClassTime(const long int &year,
          const long int &doy, const long int &hour,
          const long int &minute, const long int &second,
          const long int &milsec, const long int &micsec)
{
  Mtime_Ts.year   =year;
  Mtime_Ts.doy    =doy;
  Mtime_Ts.hour   =hour;
  Mtime_Ts.minute =minute;
  Mtime_Ts.second =second;
  Mtime_Ts.milsec =milsec;
  Mtime_Ts.micsec =micsec;
}

inline long int TBaseClassTime::hour()   const { return(Mtime_Ts.hour); }
inline long int TBaseClassTime::minute() const { return(Mtime_Ts.minute); }
inline long int TBaseClassTime::second() const { return(Mtime_Ts.second); }
inline long int TBaseClassTime::milsec() const { return(Mtime_Ts.milsec); }
inline long int TBaseClassTime::micsec() const { return(Mtime_Ts.micsec); }
inline double TBaseClassTime::float_second() const
  { return(Mtime_Ts.second+1.e-3*Mtime_Ts.milsec+1.e-6*Mtime_Ts.micsec); }

inline TBaseClassTime::operator time_kernel::time_Ts() const 
  { return(Mtime_Ts); }
inline TBaseClassTime::operator std::string() const { return(timestring()); }

/*
 * TAbsoluteTime
 * -------------
 */
inline TAbsoluteTime::TAbsoluteTime(const std::string &Itimestring):
  TBaseClassTime(Itimestring) { time_kernel::time_finish(&Mtime_Ts); }
inline TAbsoluteTime::TAbsoluteTime(char *Itimestring):
  TBaseClassTime(Itimestring) { time_kernel::time_finish(&Mtime_Ts); }

inline TAbsoluteTime::TAbsoluteTime(const time_kernel::time_Ts &Itime_Ts):
  TBaseClassTime(Itime_Ts) 
  { time_kernel::time_finish(&Mtime_Ts); }

inline TAbsoluteTime::TAbsoluteTime(const long int &year, 
           const long int &month, const long int &day, const long int &hour, 
           const long int &minute, const long int &second,
           const long int &milsec, const long int &micsec):
  TBaseClassTime(year, 1, hour, minute, second, milsec, micsec)
  { time_kernel::time_fullyear(&Mtime_Ts.year); setdate(day, month); norm(); }

inline void TAbsoluteTime::setdoy(const long int &doy)
  { Mtime_Ts.doy=doy; this->norm(); }

inline void TAbsoluteTime::setdoy(const long int &day, const long int &month)
  { time_kernel::time_setdoy(day, month, &Mtime_Ts); }

inline void TAbsoluteTime::setdate(const long int &day, const long int &month)
  { time_kernel::time_setdoy(day, month, &Mtime_Ts); }

inline void TAbsoluteTime::getdate(long int &day, long int &month) const
  { time_kernel::time_getdate(&day, &month, Mtime_Ts); }

inline bool TAbsoluteTime::isleapyear() const
  { return(time_kernel::time_isleapyear(Mtime_Ts.year) == TIME_ISLEAP); }

inline long int TAbsoluteTime::year()  const { return(Mtime_Ts.year); }
inline long int TAbsoluteTime::doy()   const { return(Mtime_Ts.doy); }
inline long int TAbsoluteTime::month() const
  { long int day, month; getdate(day, month); return(month); } 
inline long int TAbsoluteTime::day() const 
  { long int day, month; getdate(day, month); return(day); } 

inline bool TAbsoluteTime::operator==(const TAbsoluteTime &A) const
  { return(time_kernel::time_compare(Mtime_Ts,A.Mtime_Ts)==0); }
inline bool TAbsoluteTime::operator!=(const TAbsoluteTime &A) const
  { return(time_kernel::time_compare(Mtime_Ts,A.Mtime_Ts)!=0); }
inline bool TAbsoluteTime::operator<=(const TAbsoluteTime &A) const
  { return(time_kernel::time_compare(Mtime_Ts,A.Mtime_Ts)<=0); }
inline bool TAbsoluteTime::operator>=(const TAbsoluteTime &A) const
  { return(time_kernel::time_compare(Mtime_Ts,A.Mtime_Ts)>=0); }
inline bool TAbsoluteTime::operator< (const TAbsoluteTime &A) const
  { return(time_kernel::time_compare(Mtime_Ts,A.Mtime_Ts)< 0); }
inline bool TAbsoluteTime::operator> (const TAbsoluteTime &A) const
  { return(time_kernel::time_compare(Mtime_Ts,A.Mtime_Ts)> 0); }

inline TAbsoluteTime &TAbsoluteTime::operator+=(const TRelativeTime &A)
{ 
  TAbsoluteTime B(Mtime_Ts); 
  time_kernel::time_add(B.Mtime_Ts, A.Mtime_Ts, &Mtime_Ts);
  return(*this);
}
inline TAbsoluteTime &TAbsoluteTime::operator-=(const TRelativeTime &A)
{ 
  TAbsoluteTime B(Mtime_Ts); 
  time_kernel::time_sub(B.Mtime_Ts, A.Mtime_Ts, &Mtime_Ts);
  return(*this);
}

inline TAbsoluteTime TAbsoluteTime::operator+ (const TRelativeTime &A) const
  { TAbsoluteTime B(Mtime_Ts); return(B+=A); }
inline TAbsoluteTime TAbsoluteTime::operator- (const TRelativeTime &A) const
  { TAbsoluteTime B(Mtime_Ts); return(B-=A); }
inline TRelativeTime TAbsoluteTime::operator- (const TAbsoluteTime &A) const
{
    TRelativeTime B;
    time_kernel::time_sub(Mtime_Ts, A.Mtime_Ts, &B.Mtime_Ts);
    return(B); 
}

inline TAbsoluteTime &TAbsoluteTime::operator= (const time_kernel::time_Ts &A)
  { Mtime_Ts=A; time_kernel::time_finish(&Mtime_Ts); return(*this); }
inline TAbsoluteTime &TAbsoluteTime::operator= (const std::string &timestring)
  { string_read(timestring); return(*this); }
inline TAbsoluteTime &TAbsoluteTime::operator= (char *timestring)
  { char_read(timestring); return(*this); }
  
/*
 * TRelativeTime
 * -------------
 */
inline TRelativeTime::TRelativeTime(const string &Itimestring):
  TBaseClassTime("0/0/" + Itimestring) { }
inline TRelativeTime::TRelativeTime(char *Itimestring):
  TBaseClassTime("0/0/" + std::string(Itimestring)) { }

inline TRelativeTime::TRelativeTime(const time_kernel::time_Ts &Itime_Ts):
  TBaseClassTime(Itime_Ts) 
  { Mtime_Ts.year=0; norm(); }

inline TRelativeTime::TRelativeTime(const int &days, 
           const int &hour, 
           const int &minute, const int &second,
           const int &milsec, const int &micsec):
  TBaseClassTime(0, days, hour, minute, second, milsec, micsec)
  { norm(); }

inline long int TRelativeTime::days()  const { return(Mtime_Ts.doy); }

inline void TRelativeTime::nfit(const TRelativeTime &delta, long int &n, 
                 TRelativeTime &full) const
  { time_kernel::time_nfit(Mtime_Ts, delta.Mtime_Ts, &n, &full.Mtime_Ts); }

inline void TRelativeTime::div(const long int &n, TRelativeTime &frac, 
                long int &rest) const
  { time_kernel::time_div(Mtime_Ts, &frac.Mtime_Ts, n, &rest); }


inline bool TRelativeTime::operator==(const TRelativeTime &A) const
  { return(time_kernel::time_compare(Mtime_Ts,A.Mtime_Ts)==0); }
inline bool TRelativeTime::operator!=(const TRelativeTime &A) const
  { return(time_kernel::time_compare(Mtime_Ts,A.Mtime_Ts)!=0); }
inline bool TRelativeTime::operator<=(const TRelativeTime &A) const
  { return(time_kernel::time_compare(Mtime_Ts,A.Mtime_Ts)<=0); }
inline bool TRelativeTime::operator>=(const TRelativeTime &A) const
  { return(time_kernel::time_compare(Mtime_Ts,A.Mtime_Ts)>=0); }
inline bool TRelativeTime::operator< (const TRelativeTime &A) const
  { return(time_kernel::time_compare(Mtime_Ts,A.Mtime_Ts)< 0); }
inline bool TRelativeTime::operator> (const TRelativeTime &A) const
  { return(time_kernel::time_compare(Mtime_Ts,A.Mtime_Ts)> 0); }

inline TRelativeTime &TRelativeTime::operator+=(const TRelativeTime &A)
{ 
  TRelativeTime B(Mtime_Ts); 
  time_kernel::time_add(B.Mtime_Ts, A.Mtime_Ts, &Mtime_Ts);
  return(*this);
}
inline TRelativeTime &TRelativeTime::operator-=(const TRelativeTime &A)
{ 
  TRelativeTime B(Mtime_Ts); 
  time_kernel::time_sub(B.Mtime_Ts, A.Mtime_Ts, &Mtime_Ts);
  return(*this);
}

inline TRelativeTime TRelativeTime::operator+ (const TRelativeTime &A) const
  { TRelativeTime B(Mtime_Ts); return(B+=A); }
inline TRelativeTime TRelativeTime::operator- (const TRelativeTime &A) const
  { TRelativeTime B(Mtime_Ts); return(B-=A); }

inline TAbsoluteTime TRelativeTime::operator+ (TAbsoluteTime A) const
  { TRelativeTime B(Mtime_Ts); return(A+=B); }
inline TAbsoluteTime TRelativeTime::operator- (TAbsoluteTime A) const
  { TRelativeTime B(Mtime_Ts); return(A-=B); }
  
inline TRelativeTime &TRelativeTime::operator*=(const long int &n)
{ 
  time_kernel::time_Ts A(Mtime_Ts);
  time_kernel::time_mul(A, &Mtime_Ts, n); 
  return(*this);
}
inline TRelativeTime &TRelativeTime::operator/=(const long int &n)
{ 
  time_kernel::time_Ts B(Mtime_Ts);
  long int rest;
  time_kernel::time_div(B, &Mtime_Ts, n, &rest); 
  return(*this);
}
inline TRelativeTime &TRelativeTime::operator%=(const long int &n)
{ 
  time_kernel::time_Ts B(Mtime_Ts);
  long int rest;
  time_kernel::time_div(B, &Mtime_Ts, n, &rest); 
  time_kernel::time_clear(&Mtime_Ts);
  Mtime_Ts.micsec=rest;
  time_kernel::time_norm(&Mtime_Ts);
  return(*this);
}

inline TRelativeTime TRelativeTime::operator* (const long int &n) const
  { TRelativeTime B(Mtime_Ts); return(B*=n); }
inline TRelativeTime TRelativeTime::operator/ (const long int &n) const
  { TRelativeTime B(Mtime_Ts); return(B/=n); }
inline TRelativeTime TRelativeTime::operator% (const long int &n) const
  { TRelativeTime B(Mtime_Ts); return(B%=n); }

inline TRelativeTime operator* (const long int &n, const TRelativeTime &A)
  { TRelativeTime B(A); return(B*=n); }

inline long int TRelativeTime::operator/ (const TRelativeTime &A) const
{
  long int n;
  time_kernel::time_Ts full;
  time_kernel::time_nfit(Mtime_Ts, A.Mtime_Ts, &n, &full);
  return(n);
}
inline TRelativeTime TRelativeTime::operator% (const TRelativeTime &A) const
{
  long int n;
  time_kernel::time_Ts full, rest;
  time_kernel::time_nfit(Mtime_Ts, A.Mtime_Ts, &n, &full);
  time_kernel::time_sub(Mtime_Ts, full, &rest);
  return(TRelativeTime(rest));
}

inline TRelativeTime &TRelativeTime::operator= (const time_kernel::time_Ts &A)
  { Mtime_Ts=A; Mtime_Ts.year=0; norm(); return(*this); }
inline TRelativeTime &TRelativeTime::operator= (const std::string &timestring)
  { string_read("0/0/" + timestring); return(*this); }
inline TRelativeTime &TRelativeTime::operator= (char *timestring)
  { string_read("0/0/" + std::string(timestring)); return(*this); }

/*======================================================================*/
/*S*/

  //! time range
  class TRange {
    public:
      TRange(const TAbsoluteTime& begin,
             const TAbsoluteTime& end);
      TAbsoluteTime begin() const { return(Mbegin); }
      TAbsoluteTime end() const { return(Mend); }
      TRelativeTime size() const { return(Mend-Mbegin); }
      bool includes(const TAbsoluteTime&) const;
      bool includes(const TRange&) const;
      bool overlaps(const TRange&) const;
      TRange largestcommon(const TRange&) const;
      TRange smallestcommon(const TRange&) const;
      TRange& delay(const TRelativeTime&);
      TRange& advance(const TRelativeTime&);
      TRange delayedby(const TRelativeTime&) const;
      TRange advancedby(const TRelativeTime&) const;
    private:
      TAbsoluteTime Mbegin;
      TAbsoluteTime Mend;
  }; // class TRange

/*E*/ 

/*======================================================================*/
// error handling and exception class
// ----------------------------------

/*S*/
  struct Warning {
    public:
      static bool suppress_normal;
      static bool suppress_year;
      static bool suppress_any;
  }; // class Warning

  class Exception 
  {
    public:
      //! Creates exception with no explaining comments
      Exception();
      //! Creates an exception with an explanation message
      Exception(const char* message);
      //! Creates an exception with message and failed assertion
      Exception(const char* message, 
                const char* condition);
      //! Create with message, failed assertion, and code position
      Exception(const char* message, 
                const char* file,
                const int& line,
                const char* condition);
      //! Create with message and code position
      Exception(const char* message, 
                const char* file,
                const int& line);
      //! Screen report
      virtual void report() const;
      //! Issue a screen report on construction of exception
      static void report_on_construct();
      //! Issue NO screen report on construction of exception
      static void dont_report_on_construct();
    protected:
      //! Screen report
      void base_report() const;
    private:
      //! Shall we print to cerr at construction time?
      static bool Mreport_on_construct;
      //! pointer to message string
      const char* Mmessage;
      //! pointer to file name string
      const char* Mfile;
      //! pointer to line number in source file
      const int& Mline;
      //! pointer to assertion condition text string
      const char* Mcondition;
  }; // class Exception

/*======================================================================*/
//
// preprocessor macros
// ===================

/*! \brief Check an assertion and report by throwing an exception
 *
 * \ingroup group_error
 * \param C assert condition
 * \param M message of type char*
 * \param E exception class to throw
 */
#define libtime_Xassert(C,M,E) \
  if (!(C)) { throw( E ( M , __FILE__, __LINE__, #C )); }

/*! \brief Check an assertion and report by throwing an exception
 *
 * \ingroup group_error
 * \param C assert condition
 * \param M message of type char*
 */
#define libtime_assert(C,M) libtime_Xassert( C , M , libtime::Exception )

/*! \brief Abort and give a message
 *
 * \ingroup group_error
 * \param M message of type char*
 * \param E exception class to throw
 */
#define libtime_abort(M) \
  throw( libtime::Exception ( M , __FILE__, __LINE__ )) 
/*E*/ 
  
/*======================================================================*/
/*S*/

/*
 * some functions
 * --------------
 */
//! return system time
TAbsoluteTime now();
//! convert relative time to seconds
double time2double(const TRelativeTime&);
//! convert seconds to relative time
TRelativeTime double2time(const double& seconds);
/*E*/ 
  
/*======================================================================*/
/*S*/

/*
 * some string constants
 * --------------
 */
//! print information on time format
extern const char usage_time_format_string[];
/*E*/ 

} // namespace libtime
#endif // TF_LIBTIME_H_
 
/* ----- END OF libtime++.h ----- */
