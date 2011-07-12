//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL: file:///home/fochler/svnrepos/cascade_full/main/trunk/src/FPT_compare.cpp $
//$LastChangedDate: 2008-02-22 14:53:42 +0100 (Fri, 22 Feb 2008) $
//$LastChangedRevision: 33 $
//$LastChangedBy: fochler $
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------


/*
     Comparion functions for floating point types (FPT)
     
     Oliver Fochler

     last modified: 24.10.2006
*/



#include <math.h>
#include "FPT_compare.h"

//------------------------------------------------------------------------
// check for equality of floating point numbers 
// up to a relative precision of FPT_COMP_PRECISION
// 
bool FPT_COMP_E(const double a, const double b)
{
  double c1 = fabs(a);
  double c2 = fabs(b);
  double c;
  c = (c1 > c2) ? c1 : c2;

  //----------------------------------------------------------
  // special cases
  if (c == 0.0)
    return true;        // both numbers are equal to 0.0

  // either a OR b equals 0.0 (a=b=0 is excluded via the preceding condition)
  // for this case an absolute comparison is used,
  // i.e. any number x with abs(x) < epsilon^2 is considered to be equal to 0.0
  // the relative treatment would always yield a != b for any b (however small) when a = 0 
  // note that the exponent 3 is somewhat arbitrary,
  // a more consistent treatment would be desirable 
  if (a == 0.0 || b == 0.0)  
    return ( fabs(a - b) < pow(FPT_COMP_PRECISION,3.0) );
  //----------------------------------------------------------
  
  return (fabs(a - b) < c*FPT_COMP_PRECISION);  //i.e. a = b when |a-b|/max(|a|,|b|) < epsilon
}

bool FPT_COMP_E(const float a, const float b)
{
  return FPT_COMP_E(static_cast<double>(a), static_cast<double>(b));
}
//------------------------------------------------------------------------


//------------------------------------------------------------------------
// check for a < b, where FPT_COMP_E(a,b) is false
//
bool FPT_COMP_L(const double a, const double b)
{
  return( (a < b) && !FPT_COMP_E(a,b) );
}

bool FPT_COMP_L(const float a, const float b)
{
  return FPT_COMP_L(static_cast<double>(a), static_cast<double>(b));
}
//------------------------------------------------------------------------


//------------------------------------------------------------------------
// check for a > b, where FPT_COMP_E(a,b) is false
//
bool FPT_COMP_G(const double a, const double b)
{
  return( (a > b) && !FPT_COMP_E(a,b) );
}

bool FPT_COMP_G(const float a, const float b)
{
  return FPT_COMP_G(static_cast<double>(a), static_cast<double>(b));
}
//------------------------------------------------------------------------


//------------------------------------------------------------------------
// check for a <= b
//
bool FPT_COMP_LE(const double a, const double b)
{
  return( (a < b) || FPT_COMP_E(a,b) );
}

bool FPT_COMP_LE(const float a, const float b)
{
  return FPT_COMP_LE(static_cast<double>(a), static_cast<double>(b));
}
//------------------------------------------------------------------------


//------------------------------------------------------------------------
// check for a >= b
//
bool FPT_COMP_GE(const double a, const double b)
{
  return( (a > b) || FPT_COMP_E(a,b) );
}

bool FPT_COMP_GE(const float a, const float b)
{
  return FPT_COMP_GE(static_cast<double>(a), static_cast<double>(b));
}
//------------------------------------------------------------------------
