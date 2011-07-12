//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL: file:///home/fochler/svnrepos/cascade_full/main/branches/restructure_branch/src/woodsaxon.cpp $
//$LastChangedDate: 2010-07-13 00:25:52 +0200 (Tue, 13 Jul 2010) $
//$LastChangedRevision: 126 $
//$LastChangedBy: fochler $
//---------------------------------------------
//---------------------------------------------

#include <math.h>
#include "woodsaxon.h"

double WoodSaxon::densityA( double b, double z, const WoodSaxon& _w )
{
  double gamma = _w.gamma;
  double n0A = _w.n0A;
  double RA = _w.RA;
  double dA = _w.dA;

  return gamma*n0A / ( exp(( sqrt( b*b + pow( double( gamma*z ), 2.0 ) ) - RA ) / dA ) + 1.0 );
}