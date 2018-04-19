/* this is <testtime++.cc>
 * ----------------------------------------------------------------------------
 *
 * $Id: testtime++.cc,v 1.7 2004/02/02 09:31:40 tforb Exp $
 *
 * 09/08/2000 by Thomas Forbriger (IfG Stuttgart)
 *
 * comprehensive test code for libtime C++ version
 *
 * REVISIONS and CHANGES
 *    09/08/2000   V1.0   Thomas Forbriger
 *    22/12/2000   V1.1   changed namespace time to libtime
 *    13/01/2004   V1.2   test constructor for TRelativeTime from double
 *
 * ============================================================================
 */

#include"libtime++.h"
#include<iostream>

using std::cout;
using std::endl;
using std::string;

/*
 * first of all declare some useful functions
 */
template<class X>
void init_from_string(const std::string &Initializer)
{
  X Instance(Initializer);
  std::cout << Instance.timestring() 
            << " initialized from string \"" 
            << Initializer << "\"" << endl;
}

void Ainit_from_value(const long int &year, 
           const long int &month, const long int &day, 
           const long int &hour=0, const long int &minute=0,
           const long int &second=0, const long int &milsec=0,
           const long int &micsec=0)
{
  libtime::TAbsoluteTime Instance(year,month,day,hour,minute,
                               second,milsec,micsec);
  std::cout << Instance.timestring() 
            << " initialized from values: " 
            << year << " " << month << " " << day << " " 
            << hour << " " << minute << " " << second << " " 
            << milsec << " " << micsec << endl;
}

void Rinit_from_value(const long int &days, 
           const long int &hour=0, const long int &minute=0,
           const long int &second=0, const long int &milsec=0,
           const long int &micsec=0)
{
  libtime::TRelativeTime Instance(days,hour,minute,
                               second,milsec,micsec);
  std::cout << Instance.timestring() 
            << " initialized from values: " 
            << days << " " 
            << hour << " " << minute << " " << second << " " 
            << milsec << " " << micsec << endl;
}

void init_from_seconds(const double& seconds)
{
  libtime::TRelativeTime value(libtime::double2time(seconds));
  std::cout << value.timestring() 
    << " initialized from "
    << seconds 
    << " seconds"
    << endl;
}

/*
 * main testing code
 */

int main() 
{
  std::cout << "Hello world!\n";

  std::cout << "\nTesting constructors"
            << "\n====================" << endl;

  std::cout << "\nTesting TAbsoluteTime"
            << "\n---------------------" << endl;
  init_from_string<libtime::TAbsoluteTime>("0/1/1");
  init_from_string<libtime::TAbsoluteTime>("70/1/1");
  init_from_string<libtime::TAbsoluteTime>("60/1/1");
  init_from_string<libtime::TAbsoluteTime>("160/1/1");
  init_from_string<libtime::TAbsoluteTime>("1978/1/1");
  init_from_string<libtime::TAbsoluteTime>("1978/13/1");
  init_from_string<libtime::TAbsoluteTime>("1978/1/2/3/4/5/6/7/8/9");

  std::cout << endl;
  Ainit_from_value(2000,1,1);
  Ainit_from_value(000,1,50);
  Ainit_from_value(000,1,50,12,23,34,45,56);

  std::cout << "\nTesting TRelativeTime"
            << "\n---------------------" << endl;
  init_from_string<libtime::TRelativeTime>("0");
  init_from_string<libtime::TRelativeTime>("8");
  init_from_string<libtime::TRelativeTime>("123.23.34.45.56789");
  init_from_string<libtime::TRelativeTime>("123.2345.34567.45678.56789");
  init_from_seconds(1.e-6);
  init_from_seconds(1.e-3);
  init_from_seconds(1.);
  init_from_seconds(60.);
  init_from_seconds(3600.);
  init_from_seconds(86400.);
  init_from_seconds(123.345763);
  init_from_seconds(871323.345763);

  std::cout << endl;
  Rinit_from_value(0);
  Rinit_from_value(2000,1,1);
  Rinit_from_value(0,24,60,60,1000,1000);
  Rinit_from_value(00,0,0,0,0,99999999L);

// junk:
  libtime::TRelativeTime Q(600,20,10,30,40,10);
  libtime::TRelativeTime OneSecond(0,0,0,1);
  libtime::TAbsoluteTime Now(2000,9,12,16,12);
  libtime::TAbsoluteTime Then(2000,10,2,10,0);
  cout << string(Q) << endl;
  cout << string(OneSecond) << endl;
  cout << string(Q%OneSecond) << endl;
  long int numbo=Q/OneSecond;
  cout << Q/OneSecond << endl;
  cout << string(OneSecond*numbo) << endl;
  numbo=(Now-Then)/OneSecond;
  cout << numbo << " seconds from "
       << string(Now) << endl << " to " << string(Then) << endl;
  Q=numbo*OneSecond;
  cout << "these are " << string(Q) << endl << " and lead to " 
       << string(Now+Q) << endl << " from " << string(Now) << endl;
}

/* ----- END OF testtime++.cc ----- */
