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
#include <string>
#include <iostream>
#include <fstream>
#include <list>


#include "integrand_time.h"
#include "configuration.h"
#include "random.h"
#include "woodsaxon.h"


void integrand_time::operator()( const int *ndim, const double xx[], const int *ncomp, double ff[] ) const
{
  double wgt;
  ff[0] = this->operator()( xx, wgt );
}


double integrand_time::operator()( const double x[], double wgt ) const  //wgt is a dummy variable needed by VEGAS routine, x[] is {x,z,y}
{
  double gama, velocity, RA0;
  double vt, gvt, c1, c2;
  double V, xmax, xmin, xx;
  double zAl, zAr, zBl, zBr, zmax, zmin, zz;
  double yA2, yB2, ymax, yy;
  double b1, b2;

  gama = woodSaxonParameter.gamma;
  velocity = woodSaxonParameter.velocity;
  RA0 = woodSaxonParameter.RA0;

  vt = velocity * time;
  gvt = gama * vt;

  c1 = 4.0 * pow( RA0, 2.0 ) / ( pow( bImp, 2.0 ) + 4.0 * pow( gvt, 2.0 ) );

  // choose a x[fm] value
  if ( 2.0*gvt > sqrt( bImp*( 2.0*RA0 - bImp ) ) )
    xmax = bImp / 2.0 + gvt * sqrt( c1 - 1.0 );
  else
    xmax = RA0;
  xmin = bImp - xmax;

  if ( fabs( xmax - xmin ) < 1.0e-8 )
    return 0.0;

  V = xmax - xmin;
  xx = V * double( x[1] ) + xmin;

  // choose a z[fm] value
  c1 = sqrt(( RA0 + xx ) * ( RA0 - xx ) ) / gama;
  c2 = sqrt(( RA0 + xx - bImp ) * ( RA0 - xx + bImp ) ) / gama;
  zAl = vt - c1;
  zAr = vt + c1;
  zBl = -vt - c2;
  zBr = -vt + c2;
  if ( zAl > zBl )
    zmin = zAl;
  else
    zmin = zBl;
  if ( zAr < zBr )
    zmax = zAr;
  else
    zmax = zBr;

  if ( fabs( zmax - zmin ) < 1.0e-8 )
    return 0.0;

  V = V * ( zmax - zmin );
  zz = ( zmax - zmin ) * double( x[2] ) + zmin;

  // choose a y[fm] value
  yA2 = ( RA0 + gama * ( zz - vt ) ) * ( RA0 - gama * ( zz - vt ) ) - xx * xx;
  yB2 = ( RA0 + gama * ( zz + vt ) ) * ( RA0 - gama * ( zz + vt ) ) - ( bImp - xx ) * ( bImp - xx );

  if (( fabs( yA2 ) < 1.0e-8 ) || ( fabs( yB2 ) < 1.0e-8 ) )
    return 0.0;

  if ( yA2 <= yB2 )
    ymax = sqrt( yA2 );
  else
    ymax = sqrt( yB2 );

  V = V * 2.0 * ymax;
  yy = 2.0 * ymax * double( x[3] ) - ymax;

  b1 = sqrt( xx * xx + yy * yy );
  b2 = sqrt(( xx - bImp ) * ( xx - bImp ) + yy * yy );

  return double( V*WoodSaxon::densityA( b1, zz - vt, woodSaxonParameter )*WoodSaxon::densityA( b2, zz + vt, woodSaxonParameter ) );
}
