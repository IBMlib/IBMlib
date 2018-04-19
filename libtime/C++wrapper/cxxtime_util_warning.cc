/*! \file cxxtime_util_warning.cc
 * \brief issue warning (implementation)
 * 
 * ----------------------------------------------------------------------------
 * 
 * $Id: cxxtime_util_warning.cc,v 1.1 2004/02/07 17:38:16 tforb Exp $
 * \author Thomas Forbriger
 * \date 07/02/2004
 * 
 * issue warning (implementation)
 * 
 * Copyright (c) 2004 by Thomas Forbriger (BFO Schiltach) 
 * 
 * REVISIONS and CHANGES 
 *  - 07/02/2004   V1.0   Thomas Forbriger
 * 
 * ============================================================================
 */
#define TF_CXXTIME_UTIL_WARNING_CC_VERSION \
  "TF_CXXTIME_UTIL_WARNING_CC   V1.0   "
#define TF_CXXTIME_UTIL_WARNING_CC_CVSID \
  "$Id: cxxtime_util_warning.cc,v 1.1 2004/02/07 17:38:16 tforb Exp $"

#include <libtime++.h>
#include <string>
#include <iostream>

namespace libtime {
  bool Warning::suppress_normal=false;
  bool Warning::suppress_year=false;
  bool Warning::suppress_any=false;
} // namespace libtime

extern "C" {
/* 
 * Fortran calling convention:
 */
int time_util_warning__(char *caller, char *text, 
                        time_kernel::ftnlen caller_len, 
                        time_kernel::ftnlen text_len)
{
  std::string callerstring, textstring;
  int i;
  if (!(libtime::Warning::suppress_normal|| libtime::Warning::suppress_any))
  {
    for (i=0; i<caller_len; i++) { callerstring += *(caller++); }
    for (i=0; i<text_len; i++) { textstring += *(text++); }
    std::string message="ERROR ("+callerstring+"): "+textstring;
    std::cerr << message << std::endl;
  }
  return(0);
} /* time_util_warning__ */

int time_util_warning_n__(char *caller, char *text, 
                          time_kernel::integer *n, 
                          time_kernel::ftnlen caller_len, 
                          time_kernel::ftnlen text_len)
{
  std::string callerstring, textstring;
  int i;
  if (!(libtime::Warning::suppress_year|| libtime::Warning::suppress_any))
  {
    for (i=0; i<caller_len; i++) { callerstring += *(caller++); }
    for (i=0; i<text_len; i++) { textstring += *(text++); }
    std::string message="ERROR ("+callerstring+"): "+textstring;
    i=(int)*n;
    std::cerr << message << " " << i << std::endl;
  }
  return(0);
} /* time_util_warning_n__ */

} // extern "C"

/* ----- END OF cxxtime_util_warning.cc ----- */
