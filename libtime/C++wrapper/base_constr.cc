/* this is <base_constr.cc>
 * ----------------------------------------------------------------------------
 *
 * $Id: base_constr.cc,v 1.12 2004/02/07 17:38:15 tforb Exp $
 *
 * 09/08/2000 by Thomas Forbriger (IfG Stuttgart)
 *
 * TBaseClassTime constructor
 *
 * REVISIONS and CHANGES
 *    09/08/2000   V1.0   Thomas Forbriger
 *    06/09/2000   V1.1   added character array initializer
 *    22/12/2000   V1.2   changed namespace time to libtime (resolved conflict
 *                        with system time library)
 *    22/12/2003   V1.3   added function now()
 *                        correction: month in struct tm is in (0,11)
 *    02/02/2004   V1.4   moved conversion code double2time()
 *
 * ============================================================================
 */

#include "libtime++.h"
#include <cstring>
#include <ctime>
#include <iostream>

namespace libtime {

void TBaseClassTime::char_read(char *timestring)
{
  if (time_kernel::time_read(&Mtime_Ts, timestring)!=EXIT_SUCCESS)
  {
    std::cerr << "TBaseClassTime could not initialize time structure "
              << "from string:\n" << std::string(timestring) << "\n";
    std::abort();
  }
}

void TBaseClassTime::string_read(const std::string &timestring)
{
  char charstring[TIME_SLEN+2];
  std::strncpy(charstring, timestring.c_str(), TIME_SLEN+1);
  char_read(charstring);
}

}; // namespace libtime
 
/* ----- END OF base_constr.cc ----- */
