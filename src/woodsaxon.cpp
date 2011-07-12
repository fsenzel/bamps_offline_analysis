//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
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