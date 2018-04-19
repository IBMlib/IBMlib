/*! \file cxxtime_util_fatal.cc
 * \brief fatal error handling (implementation)
 * 
 * ----------------------------------------------------------------------------
 * 
 * $Id: cxxtime_util_fatal.cc,v 1.1 2004/02/07 17:38:15 tforb Exp $
 * \author Thomas Forbriger
 * \date 07/02/2004
 * 
 * fatal error handling (implementation)
 * 
 * Copyright (c) 2004 by Thomas Forbriger (BFO Schiltach) 
 * 
 * REVISIONS and CHANGES 
 *  - 07/02/2004   V1.0   Thomas Forbriger
 * 
 * ============================================================================
 */
#define TF_CXXTIME_UTIL_FATAL_CC_VERSION \
  "TF_CXXTIME_UTIL_FATAL_CC   V1.0   "
#define TF_CXXTIME_UTIL_FATAL_CC_CVSID \
  "$Id: cxxtime_util_fatal.cc,v 1.1 2004/02/07 17:38:15 tforb Exp $"

#include <libtime++.h>
#include <string>

extern "C" {
/* 
 * Fortran calling convention:
 */
int time_util_fatal__(char *caller, char *text, 
                      time_kernel::ftnlen caller_len, 
                      time_kernel::ftnlen text_len)
{
  std::string callerstring, textstring;
  int i;
  for (i=0; i<caller_len; i++) { callerstring += *(caller++); }
  for (i=0; i<text_len; i++) { textstring += *(text++); }
  std::string message="ERROR ("+callerstring+"): "+textstring;
  throw(libtime::Exception(message.c_str()));
  return(0);
} /* time_util_fatal__ */

} // extern "C"

/* ----- END OF cxxtime_util_fatal.cc ----- */
