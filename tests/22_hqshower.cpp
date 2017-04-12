//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/trunk/tests/21_checkinit.cpp $
//$LastChangedDate: 2014-03-21 13:29:33 +0100 (Fri, 21 Mar 2014) $
//$LastChangedRevision: 1641 $
//$LastChangedBy: gallmei $
//---------------------------------------------
//---------------------------------------------

#include <fstream>
#include <iostream>
#include <vector>
#include <boost/regex.hpp>
#include <initialmodel_pythiashower.h>
#include <analysis.h>

using namespace std;

const double NaddedEvents = 1000000;
const SHOWER_TYPE shower_type = inclusive_shower_spectra;
const FLAVOR_TYPE initial_parton_flavor = bottom;
const string filename_prefix = "output/inclusiveBottomQuarkShower";
const double A = 207;
const double Aatomic = 82;
const double B = 207;
const double Batomic = 82;
const double sqrtS_perNN = 2760;
const double impactParameter = 3.4;
const double P0 = 80.0;
const double pt_min_v2RAA = P0/2.0;
const double pt_max_v2RAA = 3*P0;
const double nbins_v2RAA = 30;
uint32_t seed = 0;

/**
 * @brief The main routine of BAMPS
 */
int main(int argc, char *argv[])
{
  //--------------------------------------------------------------
  // initialize the random number generator
  // defined globally in random.h and random.cpp

  if (seed == 0) seed = ran2.findSeed();
  ran2.setSeed( seed );
  cout << "seed: " << seed << endl;
  //--------------------------------------------------------------

  Particle::setCharmMass(1.3);
  Particle::setBottomMass(4.6);

  analysisRapidityRange yRange;
  std::vector<analysisRapidityRange> rapidityRanges;
  yRange.reset( 0, 0.8 );
  rapidityRanges.push_back(yRange);
  yRange.reset( 0, 0.5 );
  rapidityRanges.push_back(yRange);
  yRange.reset( 0, 1.0 );
  rapidityRanges.push_back(yRange);
  yRange.reset( 0, 2.0 );
  rapidityRanges.push_back(yRange);

  WoodSaxon WoodSaxonParameter;
  initialModel_PYTHIAShower thePythiaShower( WoodSaxonParameter, shower_type, NaddedEvents, seed, filename_prefix,
                                             A, Aatomic, B, Batomic, sqrtS_perNN, impactParameter, P0, initial_parton_flavor );

  vector<ParticleOffline> particles, particles_HFPC, particles_HFEX, particles_GLSP;
  int N_HFPC = 0, N_HFEX = 0, N_GLSP = 0;

  for( int i = 0; i < NaddedEvents; i++ )
  {
    vector<Particle> showerParticles, initialPartonPair;
    thePythiaShower.getShowerEvent( showerParticles, initialPartonPair );

    for ( int i = 0; i < showerParticles.size(); i++ )
    {
      ParticleOffline tempParticle( showerParticles[i] );
      particles.push_back( tempParticle );
    }

    if( ( initialPartonPair[0].FLAVOR == bottom && initialPartonPair[1].FLAVOR == anti_bottom ) || ( initialPartonPair[1].FLAVOR == bottom && initialPartonPair[0].FLAVOR == anti_bottom ) )
    {
      N_HFPC++;
      for ( int i = 0; i < showerParticles.size(); i++ )
      {
        ParticleOffline tempParticle( showerParticles[i] );
        particles_HFPC.push_back( tempParticle );
      }
    }
    else if( initialPartonPair[0].FLAVOR == bottom || initialPartonPair[1].FLAVOR == anti_bottom ||  initialPartonPair[1].FLAVOR == bottom || initialPartonPair[0].FLAVOR == anti_bottom )
    {
      N_HFEX++;
      for ( int i = 0; i < showerParticles.size(); i++ )
      {
        ParticleOffline tempParticle( showerParticles[i] );
        particles_HFEX.push_back( tempParticle );
      }
    }
    else if( ( initialPartonPair[0].FLAVOR == gluon && ( initialPartonPair[1].FLAVOR == gluon || Particle::mapToGenericFlavorType( initialPartonPair[1].FLAVOR ) == light_quark || Particle::mapToGenericFlavorType( initialPartonPair[1].FLAVOR ) == anti_light_quark ) )
            || ( initialPartonPair[1].FLAVOR == gluon && ( initialPartonPair[0].FLAVOR == gluon || Particle::mapToGenericFlavorType( initialPartonPair[0].FLAVOR ) == light_quark || Particle::mapToGenericFlavorType( initialPartonPair[0].FLAVOR ) == anti_light_quark ) ) )
    {
      N_GLSP++;
      for ( int i = 0; i < showerParticles.size(); i++ )
      {
        ParticleOffline tempParticle( showerParticles[i] );
        particles_GLSP.push_back( tempParticle );
      }
    }
    else
    {
      cout << "Splitting from string between quarks..." << endl;
    }
  }

  v2RAA theV2RAA( "shower", filename_prefix, rapidityRanges );
  theV2RAA.setPtBinProperties( pt_min_v2RAA, pt_max_v2RAA, nbins_v2RAA );

  theV2RAA.computeFor( bottom, particles, particles.size(), "total", 0.0, v2jets );
  theV2RAA.computeFor( bottom, particles_HFPC, particles_HFPC.size(), "HFPC", 0.0, v2jets );
  theV2RAA.computeFor( bottom, particles_HFEX, particles_HFEX.size(), "HFEX", 0.0, v2jets );
  theV2RAA.computeFor( bottom, particles_GLSP, particles_GLSP.size(), "GLSP", 0.0, v2jets );

  cout << "N_HFPC" << sep << "N_HFEX" << sep << "N_GLSP" << endl;
  cout << N_HFPC << sep << N_HFEX << sep << N_GLSP << endl;
  cout << "N_sum = " << sep << N_HFPC + N_HFEX + N_GLSP << endl;
  return EXIT_SUCCESS;
}