#include <iostream>
#include <math.h>
#include "hadronization_hq.h"
#include "particle.h"
#include "random.h"
#include "configuration.h"

using namespace ns_casc;
using namespace std;


void hadronization_hq::heavyQuarkFragmentation()
{
  double z;
  int flav;

  addedParticlesCopy = addedParticles; // copy all added particles, perform hadronization here (important to have a second instant for initial hadronization, simulating pp collisions)

  for ( int i = 0; i < addedParticles.size(); i++ )
  {
    flav = addedParticlesCopy[i].FLAVOR;
    if ( flav >= 7 && flav <= 10 )
    {
      z = getFragmentationZ ( flav );
      addedParticlesCopy[i].PX = addedParticlesCopy[i].PX * z;
      addedParticlesCopy[i].PY = addedParticlesCopy[i].PY * z;
      addedParticlesCopy[i].PZ = addedParticlesCopy[i].PZ * z;

      // D+ = c dbar, D0 = c ubar, D0bar = cbar u, D- = cbar d
      // B+ = bbar u, B0 = bbar d, B0bar = b dbar, B- = b ubar
      // for flavor use PDG particle codes:
      // 411  D+ ,   -411  D- ,   421  D0 ,  -421  D0bar
      // 521  B+ ,   -521  B- ,   511  B0 ,  -511  B0bar
      const double fraction = 0.31; // fraction of c->D+ / (c->D0  +  c->D+), from PHENIX paper 1005.1627
      if ( flav == 7 )
      {
        if ( ran2() < fraction )
          addedParticlesCopy[i].FLAVOR = dmeson_plus;
        else
          addedParticlesCopy[i].FLAVOR = dmeson_zero;
      }
      else
        if ( flav == 8 )
        {
          if ( ran2() < fraction )
            addedParticlesCopy[i].FLAVOR = dmeson_minus;
          else
            addedParticlesCopy[i].FLAVOR = dmeson_zero_bar;
        }
        else
          if ( flav == 9 )
          {
            if ( ran2() < fraction )
              addedParticlesCopy[i].FLAVOR = bmeson_zero_bar;
            else
              addedParticlesCopy[i].FLAVOR = bmeson_minus;
          }
          else
            if ( flav == 10 )
            {
              if ( ran2() < fraction )
                addedParticlesCopy[i].FLAVOR = bmeson_zero;
              else
                addedParticlesCopy[i].FLAVOR = bmeson_plus;
            }

      // Just leave the mass of the meson as the mass of the heavy quark to avoid any energy gain due to the larger mass. For the final eta distribution this has no effect anyhow since only the direction counts. However, the rapidity y distribution depends on the energy and therefore also on the heavy meson mass. If this is artificially increased here, it leads to an strange strong increase at y=0. However, previous results on RAA (and also v2) seems to not depend on this effect since it is done both initially and finally.
//       // give mesons there actual mass. However, by breaking energy conservation
//       addedParticlesCopy[i].m = ParticlePrototype::getMass( addedParticlesCopy[i].FLAVOR );

      addedParticlesCopy[i].E = sqrt ( pow ( addedParticlesCopy[i].PX, 2.0 ) + pow ( addedParticlesCopy[i].PY, 2.0 ) + pow ( addedParticlesCopy[i].PZ, 2.0 ) + pow ( addedParticlesCopy[i].m, 2.0 ) );
    }
  }
}


// samples z with the metropolis algoritm according to fragmentation function getFragmentationFunction(z)
double hadronization_hq::getFragmentationZ ( const int flav )
{
  double z;
  double z_new;
  double g = 0, g_new = 0;
  double ratio;
  double r;

  // the range in which the variable z needs to be sampled
  const double z_min = 0.0;
  const double z_max = 1.0;

  // select initial values of z
  do
  {
//     z = z_min + ran2() * ( z_max - z_min );
    z = pow ( ran2() , 1.0 / 17.0 ); // inverse function of integrated comparison function x^16. Thus, samples z according to z^16
    g = getFragmentationFunction ( z, flav );
  }
  while ( g == 0.0 );


  // number of steps in the Markov chain
  const int n_steps = 50;//50;

  // do n_steps steps
  // the current location is (t)
  // one steps consists of
  //  - proposing a new point (t') according to a certain proposal distribution
  //  - calculate the matrix element g(t') at this new point
  //  - if g(t') > g(t) accept the new point (t')
  //  - else accept the new point (t') with a probability g(t') / g(t) (last point is modified since function for sampling t is not symmetric
  //
  for ( int i = 0; i < n_steps; i++ )
  {
    do
    {
//       z_new = z_min + ran2() * ( z_max - z_min );      // propose new z using a uniform distribution over the entire range
      z_new = pow ( ran2() , 1.0 / 17.0 );
    }
    while ( z_new < z_min || z_new > z_max );

    g_new = getFragmentationFunction ( z_new, flav );           // calculate the matrix element at the proposed point

    ratio = g_new / g;                                    // ratio of g(z') / g(z)

    ratio = ratio * ( pow ( z, 16.0 ) ) / ( pow ( z_new, 16.0 ) ); // necessary if one does not use a symmetric propose fct. for t_new

    if ( ratio >= 1.0 || ran2() < ratio )     // accept if g(z') > g(z) or with probability "ratio"
    {
      z = z_new;
      g = g_new;
    }
  }

  return z;
}


// use Peterson fragmentation, cf. Peterson et al. Phys. Rev. D 27, 105-111 (1983)
double hadronization_hq::getFragmentationFunction ( const double z, const int flav )
{
  double D, epsilon;
//   const double m_q = 0.02;
  const double epsilon_c = 0.05; //pow( m_q / Mcharm , 2.0 ); // standard values: epsilon_charm = 0.05, epsilon_bottom = 0.005
  const double epsilon_b = 0.005;
  if ( flav == 7 || flav == 8 )
    epsilon = epsilon_c;
  else
    if ( flav == 9 || flav == 10 )
      epsilon = epsilon_b;
    else
      cout << "error in hadronization_hq::getFragmentationFunction( const double z, const int flav )" << endl;

  D = 1.0 / ( z * pow ( 1.0 - 1.0 / z - epsilon / ( 1.0 - z ) , 2.0 ) );

  return D;
}

