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

#include "initialmodel_mcatnlo.h"
#include "particle.h"

using std::cout;
using std::endl;
using namespace ns_casc;


void initialModel_Mcatnlo::sampleMomenta ( std::vector< Particle >& _particles )
{
  int event_tmp = 0;
  int event_counter = 0;

  std::ifstream readParticles ( filename_particleFile.c_str() );
  cout << "Read particle momenta from MC@NLO data file." << endl;
  for ( int i = 0; i < numberOfParticlesToGenerate; i++ )
  {
    Particle tempParticle;
    int flavTemp;
    int n_event_pp_input;
    // structure of file
    // flavour  energy  px  py  pz  m
    readParticles >> flavTemp >> tempParticle.E >> tempParticle.PX >> tempParticle.PY >> tempParticle.PZ >> tempParticle.m;
    tempParticle.FLAVOR = static_cast<FLAVOR_TYPE> ( flavTemp );

    tempParticle.HARD = 1;
    tempParticle.N_EVENT_pp = int ( ( i + 2 ) / 2 ); // pairs get same Event number_mcatnloFile

    if ( flavTemp <= 2 * Particle::max_N_light_flavor )
    {
      tempParticle.m = 0;
    }

    // to avoid rounding errors compute energy from momenta and mass E^2=p^2+m^2, for light partons this can be different from Pythia's value if in Pythia this parton had a mass
    tempParticle.E = sqrt ( pow ( tempParticle.PX, 2.0 ) + pow ( tempParticle.PY, 2.0 ) + pow ( tempParticle.PZ, 2.0 ) + pow ( tempParticle.m, 2.0 ) );

    double pt = sqrt ( pow ( tempParticle.PX, 2.0 ) + pow ( tempParticle.PY, 2.0 ) );

    if ( pt >= minimumPT )
    {
      _particles.push_back ( tempParticle );
    }
  }

  readParticles.close();

  // charm quarks in MC@NLO have a mass of 1.3 GeV
  // make charm quarks from pythia lighter, if Mcharm is not 1.3 GeV
  if ( Particle::N_heavy_flavor > 0 )
  {
    const double M_old_charm = 1.3;
    changeCharmMass ( _particles, M_old_charm );

    // make bottom quarks from pythia lighter, if Mbottom is not 4.6 GeV
    if ( Particle::N_heavy_flavor > 1 )
    {
      const double M_old_bottom = 4.6;
      changeBottomMass ( _particles, M_old_bottom );
    }
  }

}

