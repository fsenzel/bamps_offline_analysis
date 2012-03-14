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
#include <algorithm> // for std::count

#include "initialmodel_pythia.h"
#include "particle.h"

using std::cout;
using std::endl;
using namespace ns_casc;


void initialModel_Pythia::sampleMomenta( std::vector< Particle >& _particles )
{
  int event_tmp = 0;
  int event_counter = 0;
  
  std::ifstream readParticles( filename_particleFile.c_str() );
  cout << "Read particle momenta from PYTHIA data file." << endl;
  for ( int i = 0; i < numberOfParticlesToGenerate; i++ )
  {
    Particle tempParticle; 
    int flavTemp;
    int n_event_pp_input;
    // structure of file
    // number of pythia event  is hard?   flavour  energy  px  py  pz  m
    readParticles >> n_event_pp_input >> tempParticle.HARD >> flavTemp >> tempParticle.E >> tempParticle.PX >> tempParticle.PY >> tempParticle.PZ >> tempParticle.m;
    tempParticle.FLAVOR = static_cast<FLAVOR_TYPE>( flavTemp );
    
    if( n_event_pp_input != event_tmp )
    {
      event_tmp = n_event_pp_input;
      event_counter++;
      tempParticle.N_EVENT_pp = event_counter;
    }
    else
    {
      tempParticle.N_EVENT_pp = event_counter;
    }
    
    if ( flavTemp <= 2 * Particle::max_N_light_flavor )
    {
      tempParticle.m = 0;
    }
    
    // to avoid rounding errors compute energy from momenta and mass E^2=p^2+m^2, for light partons this can be different from Pythia's value if in Pythia this parton had a mass
    tempParticle.E = sqrt( pow( tempParticle.PX, 2.0 ) + pow( tempParticle.PY, 2.0 ) + pow( tempParticle.PZ, 2.0 ) + pow( tempParticle.m, 2.0 ) );
    
    double pt = sqrt( pow( tempParticle.PX, 2.0 ) + pow( tempParticle.PY, 2.0 ) );
    
    if( pt >= minimumPT )
    {
      _particles.push_back( tempParticle );
    }
  }
  
  readParticles.close();
  
  // charm quarks in Pythia have a mass of 1.5 GeV
  // make charm quarks from pythia lighter, if Mcharm is not 1.5 GeV
  const double M_old_charm = 1.5;
  changeCharmMass( _particles, M_old_charm );
  
  // make bottom quarks from pythia lighter, if Mbottom is not 4.8 GeV
  const double M_old_bottom = 4.8;
  changeBottomMass( _particles, M_old_bottom );

}
