/*! \file now.cc
 * \brief return system time (implementation)
 * 
 * ----------------------------------------------------------------------------
 * 
 * $Id: now.cc,v 1.1 2004/02/07 17:38:17 tforb Exp $
 * \author Thomas Forbriger
 * \date 06/02/2004
 * 
 * return system time (implementation)
 * 
 * Copyright (c) 2004 by Thomas Forbriger (BFO Schiltach) 
 * 
 * REVISIONS and CHANGES 
 *  - 06/02/2004   V1.0   Thomas Forbriger
 * 
 * ============================================================================
 */
#define TF_NOW_CC_VERSION \
  "TF_NOW_CC   V1.0   "
#define TF_NOW_CC_CVSID \
  "$Id: now.cc,v 1.1 2004/02/07 17:38:17 tforb Exp $"

#include <libtime++.h>

namespace libtime {

TAbsoluteTime now()
{
  std::time_t nowtime=std::time(NULL);
  std::tm *nowtm=std::localtime(&nowtime);
  int year=nowtm->tm_year+1900;
  int month=nowtm->tm_mon+1;
  int day=nowtm->tm_mday;
  int hour=nowtm->tm_hour;
  int minute=nowtm->tm_min; 
  int second=nowtm->tm_sec;
  TAbsoluteTime thetime(year, month, day, hour, minute, second); 
  return(thetime);
} // now()

} // namespace libtime

/* ----- END OF now.cc ----- */
