//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/offlineAnalysis/branches/heavy_quark_offline/src/configuration.cpp $
//$LastChangedDate: 2012-02-29 11:20:46 +0100 (Wed, 29 Feb 2012) $
//$LastChangedRevision: 423 $
//$LastChangedBy: uphoff $
//---------------------------------------------
//---------------------------------------------

#include <math.h>
#include <iostream>
#include "initialmodel_jpsi.h"
#include "random.h"
#include "FPT_compare.h"
#include "interpolation_iniJpsi.h"
#include "binning.h"
#include "particle.h"
#include "configuration.h"

using namespace std;
using namespace ns_casc;

initialModel_Jpsi::initialModel_Jpsi( const config& _config, WoodSaxon& _WoodSaxonParameter )
: initialModelWS(_config), sigmaAbs( _config.getSigmaAbs() ), agN( _config.getJpsiagN() ), shadowing_model( _config.getShadowingModel() ), nEventsAA( _config.getNaddedEvents() ), testparticles( _config.getTestparticles() * _config.getJpsiTestparticles() )
{
  double Tab;
  
  impactParameter = _config.getImpactParameter();
  sqrtS_perNN = _config.getSqrtS();
  if (!WoodSaxonParameter.Calculate( A, impactParameter, sqrtS_perNN))
  {
    std::string errMsg = "Impact parameter b too large. b > 2 R_A0";
    throw eJpsi_error( errMsg );
  }
  _WoodSaxonParameter = WoodSaxonParameter;
  
  cout << "======= Generating data sets for sampling of initial state =======" << endl;
  generateTimeDistributionWS(Tab);
  cout << "++++  Tab = " << Tab << "1/mb" << endl;
  cout << "==================================================================" << endl;
  
  theInterpolation_dndptdy.configure( sqrtS_perNN, impactParameter, sigmaAbs, agN, shadowing_model );
}


void initialModel_Jpsi::populateParticleVector( std::vector< Particle >& _particles )
{
  double total_number_jpsi_one_Au_collision = 0.0; // integration over dndydpt, gives the number of jpsis in one collision
  
  if( sqrtS_perNN == 200.0)
  {
    if( sigmaAbs == 2.8 && agN == 0.1 && shadowing_model == eps08 )
    {
      if( impactParameter == 0.0)
        total_number_jpsi_one_Au_collision = 0.04716;
      else if( impactParameter == 3.3)
        total_number_jpsi_one_Au_collision = 0.03716;
      else if( impactParameter == 4.6)
        total_number_jpsi_one_Au_collision = 0.03030;
      else if( impactParameter == 5.8)
        total_number_jpsi_one_Au_collision = 0.02374;
      else if( impactParameter == 8.2)
        total_number_jpsi_one_Au_collision = 0.01191;
      else if( impactParameter == 10.3)
        total_number_jpsi_one_Au_collision = 0.004515;
    }
    else if( sigmaAbs == 1.5 && agN == 0.1 && shadowing_model == eps08 )
    {
      if( impactParameter == 0.0)
        total_number_jpsi_one_Au_collision = 0.057675;
      else if( impactParameter == 8.2)
        total_number_jpsi_one_Au_collision = 0.0139615;
    }
    else if( sigmaAbs == 0.0 && agN == 0.1 && shadowing_model == eps08 )
    {
      if( impactParameter == 0.0)
        total_number_jpsi_one_Au_collision = 0.073531;
      else if( impactParameter == 8.2)
        total_number_jpsi_one_Au_collision = 0.0169033;
    }
    else if( sigmaAbs == 0.0 && agN == 0.0 && shadowing_model == none )
    {
      if( impactParameter == 0.0)
        total_number_jpsi_one_Au_collision = 0.0862351;
      else if( impactParameter == 8.2)
        total_number_jpsi_one_Au_collision = 0.018879;
    }
  }
  else if( sqrtS_perNN == 2760.0)
  {
    if( sigmaAbs == 2.8 && agN == 0.1 && shadowing_model == eps08 )
    {

    }
    else if( sigmaAbs == 1.5 && agN == 0.1 && shadowing_model == eps08 )
    {

    }
    else if( sigmaAbs == 0.0 && agN == 0.1 && shadowing_model == eps08 )
    {

    }
    else if( sigmaAbs == 0.0 && agN == 0.1 && shadowing_model == eps09 )
    {
      if( impactParameter == 3.6)
        total_number_jpsi_one_Au_collision = 0.46568;
      else if( impactParameter == 9.7)
        total_number_jpsi_one_Au_collision = 0.0815199;
    }
    else if( sigmaAbs == 0.0 && agN == 0.0 && shadowing_model == eps09 )
    {
      if( impactParameter == 3.6)
        total_number_jpsi_one_Au_collision = 0.46187;
    }
    else if( sigmaAbs == 0.0 && agN == 0.0 && shadowing_model == none )
    {
      if( impactParameter == 3.6)
        total_number_jpsi_one_Au_collision = 0.656866;
      else if( impactParameter == 9.7)
        total_number_jpsi_one_Au_collision = 0.101697;
    }
  }
  
//   total_number_jpsi_one_Au_collision = 0;
  if ( total_number_jpsi_one_Au_collision == 0.0 )
  {
    std::string errMsg = "Number of initial Jpsi not set.";
    throw eJpsi_error( errMsg );
  }

  int number_jpsi = int(total_number_jpsi_one_Au_collision * nEventsAA * testparticles);
  double remainder = total_number_jpsi_one_Au_collision * nEventsAA * testparticles - number_jpsi; // int() above rounds down, but we want to take the remainder also into account
  if(ran2() < remainder)
    number_jpsi++;
  
  for(int i = 0; i < number_jpsi; i++)
  {
    Particle tempParticle;
    
    sample_PXYZE_FLAV_singleParticle( tempParticle, _particles.size() );
    sample_TXYZ_singleParticle( tempParticle );
    
    _particles.push_back( tempParticle );
  }
}



void initialModel_Jpsi::sample_PXYZE_FLAV_singleParticle( Particle& _tempParticle, const int ParticleNumber )
{
  double pt, y, phi;

  _tempParticle.m = Particle::getMass( jpsi );
  _tempParticle.N_EVENT_pp = ParticleNumber + 1;
  _tempParticle.N_EVENT_AA = int( ran2() * nEventsAA ) + 1; // event_AA_tmp is integer in [1;nEventsAA]
  _tempParticle.FLAVOR = jpsi;

  // get y and pt
  sample_metropolis_dndptdy( pt, y );
  
  // angle of jpsi
  phi = ran2() * 2.0 * M_PI;

  // momenta and energy
  _tempParticle.Mom.Px() = pt * sin( phi );
  _tempParticle.Mom.Py() = pt * cos( phi );
  _tempParticle.Mom.Pz() = sqrt( ( pow( pt, 2.0 ) + pow( _tempParticle.m, 2.0 ) ) * pow( exp( y ) - exp( -y ) , 2.0 ) / 4.0 );
  // consider also negativ pz. y is only sampled for positiv y since tables are only for positiv y and y is symmetric around 0. Here, substitute randomly pz by -pz:
  if ( ran2() < 0.5 )
    _tempParticle.Mom.Pz() = -_tempParticle.Mom.Pz();

  _tempParticle.Mom.E() = sqrt( _tempParticle.Mom.vec2() + pow( _tempParticle.m, 2.0 ) );
}


void initialModel_Jpsi::sample_metropolis_dndptdy( double& pt_arg, double& y_arg )
{
  double pt, y, r;
  double pt_new, y_new;
  double g = 0, g_new = 0;
  double ratio;
  
    // the ranges in which the variables pt and y need to be sampled
  double pt_min, pt_max, y_min, y_max;
  if( sqrtS_perNN == 200 )
  {
    // the ranges in which the variables pt and y need to be sampled
    pt_min = 0.0;
    pt_max = 5.0;
    // y is sampled in the positiv range only, since the distribution is symmetric around 0. Negativ y are taken into account at the conversin to pz
    y_min = 0.0;
    y_max = 5.0;
  }
  else if( sqrtS_perNN == 2760 )
  {
    // the ranges in which the variables pt and y need to be sampled
    pt_min = 0.0;
    pt_max = 10.0;
    // y is sampled in the positiv range only, since the distribution is symmetric around 0. Negativ y are taken into account at the conversin to pz
    y_min = 0.0;
    y_max = 7.0;
  }
  else
  {
    std::string errMsg = "Error in  initialModel_Jpsi::sample_metropolis_dndptdy().";
    throw eJpsi_error( errMsg );
  }

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
