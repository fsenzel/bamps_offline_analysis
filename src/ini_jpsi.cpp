#include <math.h>
#include <iostream>
#include "ini_jpsi.h"
#include "random.h"
#include "FPT_compare.h"
#include "interpolation_iniJpsi.h"
#include "binning.h"
#include "particle.h"
#include "configuration.h"

using namespace std;

extern coord *addedPartcl;
extern int numberAdded;
extern int Ntest;

binning ptbins("output/ptbins.dat", 0.0, 5.0, 100);

ini_jpsi::ini_jpsi( const double sqrtS_arg, const double Bimp_arg, const double sigmaAbs_arg, const double agN_arg, const shadowModelJpsi shadowing_model_arg )
: sqrtS(sqrtS_arg), impact_parameter(Bimp_arg), sigmaAbs(sigmaAbs_arg), agN(agN_arg), shadowing_model(shadowing_model_arg)
{
  theInterpolation_dndptdy.init( sqrtS, impact_parameter, sigmaAbs, agN, shadowing_model );
}

void ini_jpsi::sample_jpsis()
{
  double total_number_jpsi_one_Au_collision = 0.0; // integration over dndydpt, gives the number of jpsis in one collision
  
  if( sqrtS == 200.0)
  {
    if( sigmaAbs == 2.8 && agN == 0.1 && shadowing_model == eps08 )
    {
      if( impact_parameter == 0.0)
        total_number_jpsi_one_Au_collision = 0.04716;
      else if( impact_parameter == 3.3)
        total_number_jpsi_one_Au_collision = 0.03716;
      else if( impact_parameter == 4.6)
        total_number_jpsi_one_Au_collision = 0.03030;
      else if( impact_parameter == 5.8)
        total_number_jpsi_one_Au_collision = 0.02374;
      else if( impact_parameter == 8.2)
        total_number_jpsi_one_Au_collision = 0.01191;
      else if( impact_parameter == 10.3)
        total_number_jpsi_one_Au_collision = 0.004515;
    }
    else if( sigmaAbs == 1.5 && agN == 0.1 && shadowing_model == eps08 )
    {
      if( impact_parameter == 0.0)
        total_number_jpsi_one_Au_collision = 0.057675;
      else if( impact_parameter == 8.2)
        total_number_jpsi_one_Au_collision = 0.0139615;
    }
    else if( sigmaAbs == 0.0 && agN == 0.1 && shadowing_model == eps08 )
    {
      if( impact_parameter == 0.0)
        total_number_jpsi_one_Au_collision = 0.073531;
      else if( impact_parameter == 8.2)
        total_number_jpsi_one_Au_collision = 0.0169033;
    }
    else if( sigmaAbs == 0.0 && agN == 0.0 && shadowing_model == none )
    {
      if( impact_parameter == 0.0)
        total_number_jpsi_one_Au_collision = 0.0862351;    
      else if( impact_parameter == 8.2)
        total_number_jpsi_one_Au_collision = 0.018879;
    }
  }
  
  
   
  double KInicharm = theConfig.getKIniCharm();
  
  int number_jpsi = int(total_number_jpsi_one_Au_collision * KInicharm * Ntest);
  double remainder = total_number_jpsi_one_Au_collision * KInicharm * Ntest - number_jpsi; // int() above rounds down, but we want to take the remainder also into account
  if(ran2() < remainder)
    number_jpsi++;

//   cout << "Jpsi: " << number_jpsi << endl;
  
  for(int i = 0; i < number_jpsi; i++)
  {
    numberAdded++;
    sample_one_jpsi(numberAdded);
  }
  
  
  ptbins.print();
}


// void ini_jpsi::sample_one_jpsi( const int partclNmb )
// {
//   double pt, y, phi;
//   
//   double MASS,PX,PY,PZ,E;
// 
//   MASS = Mjpsi;
// 
//   // get y and pt
//   sample_metropolis_dndptdy( pt, y );
//   
//   // angle of jpsi
//   phi = ran2() * 2.0 * M_PI;
// 
//   // momenta and energy
//   PX = pt * sin( phi );
//   PY = pt * cos( phi );
//   PZ = ( pow( pt, 2.0 ) + pow( MASS, 2.0 ) ) * pow( exp( y ) - exp( -y ) , 2.0 ) / 4.0;
//   // consider also negativ pz. y is only sampled for positiv y since tables are only for positiv y and y is symmetric around 0. Here, substitute randomly pz by -pz:
//   if ( ran2() < 0.5 )
//     PZ = -PZ;
// 
//   E = sqrt( pow( PX, 2.0 ) + pow( PY, 2.0 ) + pow( PZ, 2.0 ) + pow( MASS, 2.0 ) );
//   
//   
// //   // Kai has another definition of the J/psi's mass:
// //   const double M_jpsi = 3.6; // GeV
// //   // make Jpsi heavier by breaking energy conservation
// //   MASS = M_jpsi;
// //   E = sqrt( pow( PX, 2.0 ) + pow( PY, 2.0 ) + pow( PZ, 2.0 ) + pow( MASS, 2.0 ) );
// 
// 
//   double y_test = 0.5 * log( (E + PZ) / (E - PZ) );
//   double pt_test = sqrt( pow( PX, 2.0 ) + pow( PY, 2.0 ) );
//   
//   if( fabs(y_test) < 0.5 )
//     ptbins.add(pt_test);
// 
// }



void ini_jpsi::sample_one_jpsi( const int partclNmb )
{
  double pt, y, phi;

  addedPartcl[partclNmb].MASS = Mjpsi;
  addedPartcl[partclNmb].N_EVENT = partclNmb;
  addedPartcl[partclNmb].FLAVOR = 50;

  // get y and pt
  sample_metropolis_dndptdy( pt, y );
  
  // angle of jpsi
  phi = ran2() * 2.0 * M_PI;

  // momenta and energy
  addedPartcl[partclNmb].PX = pt * sin( phi );
  addedPartcl[partclNmb].PY = pt * cos( phi );
  addedPartcl[partclNmb].PZ = sqrt( ( pow( pt, 2.0 ) + pow( addedPartcl[partclNmb].MASS, 2.0 ) ) * pow( exp( y ) - exp( -y ) , 2.0 ) / 4.0 );
  // consider also negativ pz. y is only sampled for positiv y since tables are only for positiv y and y is symmetric around 0. Here, substitute randomly pz by -pz:
  if ( ran2() < 0.5 )
    addedPartcl[partclNmb].PZ = -addedPartcl[partclNmb].PZ;

  addedPartcl[partclNmb].E = sqrt( pow( addedPartcl[partclNmb].PX, 2.0 ) + pow( addedPartcl[partclNmb].PY, 2.0 ) + pow( addedPartcl[partclNmb].PZ, 2.0 ) + pow( addedPartcl[partclNmb].MASS, 2.0 ) );
  
  // this is now done earlier, by changing the Jpsi mass to 3.6 in configuration.h
//   // Kai has another definition of the J/psi's mass:
//   const double M_jpsi = 3.6; // GeV
//   // make Jpsi heavier by breaking energy conservation
//   addedPartcl[partclNmb].MASS = M_jpsi;
//   addedPartcl[partclNmb].E = sqrt( pow( addedPartcl[partclNmb].PX, 2.0 ) + pow( addedPartcl[partclNmb].PY, 2.0 ) + pow( addedPartcl[partclNmb].PZ, 2.0 ) + pow( addedPartcl[partclNmb].MASS, 2.0 ) );

 

}


void ini_jpsi::sample_metropolis_dndptdy( double& pt_arg, double& y_arg )
{
  double pt, y, r;
  double pt_new, y_new;
  double g = 0, g_new = 0;
  double ratio;

  // the ranges in which the variables pt and y need to be sampled
  const double pt_min = 0.0;
  const double pt_max = 5.0;
  // y is sampled in the positiv range only, since the distribution is symmetric around 0. Negativ y are taken into account at the conversin to pz
  const double y_min = 0.0;
  const double y_max = 5.0;

  // randomly select initial values of pt and y, such that
  do
  {
    pt = pt_min + ran2()*(pt_max - pt_min);
//     r = ran2();
//     pt = - log( 1.0 + r * ( exp( -pt_max ) - 1.0 ) ); // inverse function of integrated comparison function exp(-x). Holds only if pt_min = 0.0

    y = y_min + ran2() * ( y_max - y_min );

    g = theInterpolation_dndptdy.getdN( y, pt );
  }
  while ( FPT_COMP_E( g, 0.0 ) );


  // number of steps in the Markov chain
  const int n_steps = 50;

  // do n_steps steps
  // the current location is (pt,y)
  // one steps consists of
  //  - proposing a new point (pt',y') according to a certain proposal distribution
  //  - calculate the matrix element g(pt',y') at this new point
  //  - if g(pt',y') > g(pt,y) accept the new point (pt',y')
  //  - else accept the new point (pt',y') with a probability g(pt',y') / g(pt, y)
  //
  for ( int i = 0; i < n_steps; i++ )
  {
    do
    {
      pt_new = pt_min + ran2()*(pt_max - pt_min);          // propose new pt using a uniform distribution over the entire range
//       r = ran2();
//       pt_new = - log( 1.0 + r * ( exp( -pt_max ) - 1.0 ) ); // inverse function of integrated comparison function exp(-x). Holds only if pt_min = 0.0
    }
    while ( pt_new < pt_min || pt_new > pt_max );

    do
    {
      y_new = y_min + ran2() * ( y_max - y_min );   // propose new y using a uniform distribution over the entire range
    }
    while ( y_new < y_min || y_new > y_max );     // check that the new values are in range

    g_new = theInterpolation_dndptdy.getdN( y_new, pt_new );             // calculate the matrix element at the proposed point

    ratio = g_new / g;                                    // ratio of g(pt',y') / g(pt,y)

//     ratio = ratio *  exp( -pt )  /  exp( -pt_new ) ; // necessary if one does not use a symmetric propose fct. for pt

    if ( FPT_COMP_GE( ratio, 1.0 ) || ran2() < ratio )     // accept if g(pt',y') > g(pt,y) or with probability g(pt',y') / g(pt,y)
    {
      pt = pt_new;
      y = y_new;
      g = g_new;
    }

  }

  pt_arg = pt;
  y_arg = y;

}
