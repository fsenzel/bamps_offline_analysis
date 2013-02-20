
#include <math.h>
#include <string>
#include <iostream>
#include "initialmodel.h"
#include "random.h"

using std::cout;
using std::endl;


void initialModel::setUniqueID( std::vector<Particle>& _particles )
{
  Particle::unique_id_counter = 0;
  for ( int i = 0; i < _particles.size(); i++ )
  {
    _particles[i].unique_id = Particle::unique_id_counter;
    ++Particle::unique_id_counter;
  }
}


initialModelWS::initialModelWS( const config& _config ) :
  A( _config.getA() ),
  Aatomic( _config.getAatomic() ),
  B( _config.getB() ),
  Batomic( _config.getBatomic() )
{
}

initialModelWS::~initialModelWS( )
{
}

void initialModelWS::sample_XYZ_WoodSaxon(double T, double &X, double &Y, double &Z) const
{
  double vt, gv, gvt, tc1, tc2, tc3, xx, bA, bB, zA, zB, fds, fd;
  double zmax, zmin, xmax, xmin, ymax, ymin, max;
  double c1, c2;

  double densityA_max;
  
  double m = 2.0;
  //   double m = 2.718;
  double L_z;// L_z is z component of radius where densityA is only 1/m-th of maximum value

  const double RA0 = WoodSaxonParameter.RA0; // just a shortcut
  double RR;

  vt = WoodSaxonParameter.velocity * fabs( T );
  gv = WoodSaxonParameter.gamma * WoodSaxonParameter.velocity;
  gvt = gv * fabs( T );

  c1 = sqrt(( RA0 + impactParameter ) * fabs( RA0 - impactParameter ) );
  tc1 = ( RA0 - c1 ) / 2.0 / gv;
  tc2 = sqrt( impactParameter * ( 2.0 * RA0 - impactParameter ) ) / 2.0 / gv;
  tc3 = ( RA0 + c1 ) / 2.0 / gv;

  xx = impactParameter / 2.0 + gvt * sqrt( RA0*RA0 / ( impactParameter*impactParameter / 4.0 + gvt*gvt ) - 1.0 );

  // z-axis in overlap region
  if (( RA0 > impactParameter ) && ( fabs( T ) >= tc1 ) && ( fabs( T ) <= tc3 ) ) 
  { 
    RR = RA0;
  }
  else
  {
    RR = sqrt(( RA0 + impactParameter - xx ) * ( RA0 - impactParameter + xx ) );
  }
  zmin = vt - RR / WoodSaxonParameter.gamma;
  zmax = -zmin;

  // x-axis in overlap region
  xmax = ( fabs( T ) <= tc2 ) ? RA0 : xx;
  xmin = impactParameter - xmax;

  // y-axis in overlap region
  ymax = sqrt( RA0*RA0 - impactParameter*impactParameter / 4.0 - gvt*gvt );
  ymin = -ymax;

  // the maximum of nA(s,z-vt)*nB(s-b,z+vt) in the overlap region
  max = densityA( impactParameter / 2.0, -WoodSaxonParameter.velocity * T ) * densityA( impactParameter / 2.0, WoodSaxonParameter.velocity * T );

  // sampling of position
  do
  {
    X = ( xmax - xmin ) * ran2() + xmin;
    Y = ( ymax - ymin ) * ran2() + ymin;
    Z = ( zmax - zmin ) * ran2() + zmin;

    bA = sqrt( X*X + Y*Y );
    bB = sqrt(( X - impactParameter )*( X - impactParameter ) + Y*Y );

    if (( bA > RA0 ) || ( bB > RA0 ) ) 
    {
      fds = -1.0;
    }
    else
    {
      zA = sqrt(( RA0 + bA ) * ( RA0 - bA ) ) / WoodSaxonParameter.gamma;
      zB = sqrt(( RA0 + bB ) * ( RA0 - bB ) ) / WoodSaxonParameter.gamma;
      c1 = fabs( Z - WoodSaxonParameter.velocity * T );
      c2 = fabs( Z + WoodSaxonParameter.velocity * T );
      fds = (( c1 > zA ) || ( c2 > zB ) ) ? -1.0 : densityA( bA, c1 ) * densityA( bB, c2 );
    }

    if ( fds > max ) cout << "problem in initialModel::sample_XYZ_WoodSaxon" << endl;

    fd = max * ran2();
  }
  while ( fds < fd );
}



void initialModelWS::sample_T_WoodSaxon(double &T) const
{
  T = distrTime->GenerateNumber();
  if (ran2() > 0.5) 
  {
    T = -T;
  }
}


void initialModelWS::sample_TXYZ_singleParticle( Particle& _particle )
{
  double T, X, Y, Z;
  
  sample_T_WoodSaxon(T);
  _particle.T = T;
  
  sample_XYZ_WoodSaxon(T, X,Y,Z);
  _particle.X = X - impactParameter / 2.0;    // shift into the correct coordinate system 
  _particle.Y = Y;
  _particle.Z = Z;
}




/**
 * This routine generates the probability distribution and the cdf for
 * the generation of random numbers for the time T.
 *
 * We assume full symmetry around zero and thus only calculate one
 * half of the problem.
 **/
void initialModelWS::generateTimeDistributionWS(double &T_AB)
{
  if ( !distrTime ) // only do this stuff if the distribution has not been computed before (i.e. the shared_ptr is empty)
  {
    int nn, ncut;
    double dt, tmp;
    double tgral, sd, chi2a; // for VEGAS
    double *dist, *tt, *distin;
    
    integrand_time ftime;
    ftime.setB( impactParameter );
    ftime.setWoodSaxonParameter( WoodSaxonParameter );
    
    double tmax = sqrt( 4 * pow( WoodSaxonParameter.RA0, 2 ) - pow( impactParameter, 2 ) ) / ( 2 * WoodSaxonParameter.velocity * WoodSaxonParameter.gamma );
    
    cout << "---- generateTimeDistributionWS:" << endl;
    //  cout << "   tmax = " << tmax << endl;
    
    nn = 100;
    dt = tmax / nn;
    tt = new double[nn+1];     // the time axis
    dist = new double[nn+1];   // the function values
    distin = new double[nn+1]; // the integrated function values
    
    // reset the function values, populate the time array:
    
    for ( int i = 0; i < nn+1; i++ ) { dist[i] = 0.0; }
    tt[nn] = 0.0;
    for ( int i = nn - 1; i >= 0; i-- )
    {
      tt[i] = tt[i+1] - dt;
    }
    
    // calculate the function values:
    
    ncut = 0;
    for ( int i = nn; i >= 0; i-- )
    {
      ftime.setTime( -tt[i] );
      vegas( 3, ftime, &tgral, &sd, &chi2a );
      dist[i] = double( tgral );
      
      if ( tgral < 1.0e-7 )
      {
        ncut = i;
        i = -1; // quit the loop
      }
    }
    if (ncut > 0) ncut--;
    
    // calculate the integral of the function values:
    
    distin[0] = 0.0;
    for ( int i = 1; i < nn+1; i++ )
    {
      distin[i] = distin[i-1] + 0.5 * ( dist[i-1] + dist[i] ) * dt;
    }
    
    T_AB = (2.0*distin[nn]) * 2.0 * WoodSaxonParameter.velocity / 10.;// 1/10 due to the rescale 1/fm²->1/10mb
    
    // normalize the integral:
    for ( int i = 1; i < nn+1; i++ ) { distin[i] /= distin[nn]; }
    
    // now we generate the random generator:
    distrTime.reset( new ranGen_Distr(&tt[ncut], &distin[ncut], nn-ncut+1, interp_cspline) );
    
    cout << "---- generateTimeDistributionWS: finished" << endl;
  }
}

////////////////////////////////////////////////////////////

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

  return double( V*woodSaxonParameter.densityA( b1, zz - vt)*woodSaxonParameter.densityA( b2, zz + vt ) );
}
