/*! \file usage.cc
 * \brief print usage information (implementation)
 * 
 * ----------------------------------------------------------------------------
 * 
 * $Id: usage.cc,v 1.1 2005/07/15 15:36:54 tforb Exp $
 * \author Thomas Forbriger
 * \date 15/07/2005
 * 
 * print usage information (implementation)
 * 
 * Copyright (c) 2005 by Thomas Forbriger (BFO Schiltach) 
 * 
 * REVISIONS and CHANGES 
 *  - 15/07/2005   V1.0   Thomas Forbriger
 * 
 * ============================================================================
 */
#define TF_USAGE_CC_VERSION \
  "TF_USAGE_CC   V1.0"
#define TF_USAGE_CC_CVSID \
  "$Id: usage.cc,v 1.1 2005/07/15 15:36:54 tforb Exp $"

#include <libtime++.h>

namespace libtime {

  const char usage_time_format_string[]=
  {
    TF_LIBTIME_H_ "\n"
    "Times and dates can be specified as follows:" "\n"
    "1. Absolute time" "\n"
    "   defines a specific date, like 2005/7/15_12:30 specifies the" "\n"
    "   15th of July in 2005 at half past twelve. The general format" "\n"
    "   for absolute time is: yyyy/mm/dd/HH/MM/SS.SSSSSS" "\n"
    "2. Relative time" "\n"
    "   specifies a time range like 0/8/15, which means 8 hours and" "\n"
    "   fifteen minutes. The general format for relative time is" "\n"
    "   dd/HH/MM/SS.SSSSSS" "\n"
    "In the format strings given above, '/' serves as a field separator." "\n"
    "This character may be replaced by any non-digit character. Leading" "\n"
    "zeroes may be omitted in all fields, where yyyy means year, mm" "\n"
    "means month, dd means day, HH means hour, mm means minute and" "\n"
    "SS.SSSSSS means seconds. Seconds are specified as a floating point" "\n"
    "value with precision down to microseconds. Two-digit year values" "\n"
    "will be interpreted as years in the 20th century (for values larger" "\n"
    "than 69) or in the 21st century (for values smaller than 70)." "\n"
    "Fields from HH on may be omitted." "\n"
  }; // usage_time_format_string

}

/* ----- END OF usage.cc ----- */
