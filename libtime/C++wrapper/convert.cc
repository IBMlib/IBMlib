/*! \file convert.cc
 * \brief convert time representation (implementation)
 * 
 * ----------------------------------------------------------------------------
 * 
 * $Id: convert.cc,v 1.1 2004/02/07 17:38:15 tforb Exp $
 * \author Thomas Forbriger
 * \date 06/02/2004
 * 
 * convert time representation (implementation)
 * 
 * Copyright (c) 2004 by Thomas Forbriger (BFO Schiltach) 
 * 
 * REVISIONS and CHANGES 
 *  - 06/02/2004   V1.0   Thomas Forbriger
 * 
 * ============================================================================
 */
#define TF_CONVERT_CC_VERSION \
  "TF_CONVERT_CC   V1.0   "
#define TF_CONVERT_CC_CVSID \
  "$Id: convert.cc,v 1.1 2004/02/07 17:38:15 tforb Exp $"

#include <libtime++.h>

namespace libtime {

TRelativeTime double2time(const double& seconds)
{
  time_kernel::time_Ts thetime_Ts(TRelativeTime(0));
  double remain=seconds;
  typedef long int li;
  thetime_Ts.second=li(seconds);
  remain-=double(thetime_Ts.second);
  remain*=1.e3;
  thetime_Ts.milsec=li(remain);
  remain-=double(thetime_Ts.milsec);
  remain*=1.e3;
  thetime_Ts.micsec=li(remain);
  return(TRelativeTime(thetime_Ts));
} // TRelativeTime double2time(const double& seconds)

double time2double(const TRelativeTime& rtime)
{
  double retval(0.);
  retval=rtime.float_second();
  retval += 60.*rtime.minute();
  retval += 3600.*rtime.hour();
  retval += 24.*3600.*rtime.days();
  return(retval);
} // double time2float(const TRelativeTime&)

} // namespace libtime

/* ----- END OF convert.cc ----- */
