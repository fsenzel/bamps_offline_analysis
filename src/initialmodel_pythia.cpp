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
#include <algorithm> // for std::count

#include <sstream>
#include <list>

#include "configuration.h"
#include "initialmodel_pythia.h"
#include "random.h"
#include "particle.h"
#include "FPT_compare.h"

using std::cout;
using std::endl;
using namespace ns_casc;


initialModel_Pythia::initialModel_Pythia( const config& _config, WoodSaxon& _WoodSaxonParameter) :
  initialModelWS(_config) ,
  numberOfParticlesToGenerate ( 0 ),
  numberOfTestparticles ( _config.getTestparticles() )
{
  double Tab;

  impactParameter = _config.getImpactParameter();
  sqrtS_perNN = _config.getSqrtS();
  filename_pythiaParticleFile = _config.getPythiaParticleFile();
  if (!WoodSaxonParameter.Calculate( A, impactParameter, sqrtS_perNN))
  {
    std::string errMsg = "Impact parameter b too large. b > 2 R_A0";
    throw ePythia_error( errMsg );
  }
  WoodSaxonParameter.Calculate( A, impactParameter, sqrtS_perNN);
  _WoodSaxonParameter = WoodSaxonParameter;


  cout << "======= Generating data sets for sampling of initial state =======" << endl;
  generateTimeDistributionWS(Tab);
  cout << "++++  Tab = " << Tab << "1/mb" << endl;
  cout << "==================================================================" << endl;

  cout << "PYTHIA particle data file: " << filename_pythiaParticleFile << endl;

  std::ifstream countPythiaParticles( filename_pythiaParticleFile.c_str() );
  if ( countPythiaParticles.good() )
  {
    numberOfParticlesToGenerate = std::count( std::istreambuf_iterator<char>( countPythiaParticles ), std::istreambuf_iterator<char>(), '\n' ); // counts number of lines in PYTHIA particle data file
  }
  else
  {
    string errMsg = "Error at opening PYTHIA particle data file.";
    throw ePythia_error( errMsg );
  }

  countPythiaParticles.close();
}



void initialModel_Pythia::populateParticleVector( std::vector< ParticleOffline >& _particles )
{
  /**
  * Reserve memory for the Particle vector. Due to possible particle creation this size should be chosen somewhat larger
  * than the initial number of particles. Doing this, push_back operations to add new particles won't lead to internal
  * memory re-allocation in the vector, which could possibly be time consuming.
  * However push_back operations are always possible, even when the reserved size is used up (as long as there is physical
  * memory available of course). And modern day compilers probably optimize to the extent where it doesn't really matter -
  * but still, it's good practice.
  */
  _particles.reserve( static_cast<int>( numberOfParticlesToGenerate * 1.2 ) );
  /**
  * Now the particle vector is re-sized to hold the initial number of particles (this is NOT done when reserving memory!).
  * The particles are initialised with the standard constructor, actual attributes such as momentum etc. MUST be set later!
  */
  _particles.resize( numberOfParticlesToGenerate );

  sampleMomenta( _particles );
  samplePositions( _particles );
}




void initialModel_Pythia::sampleMomenta( std::vector< ParticleOffline >& _particles )
{

  std::ifstream readPythiaParticles( filename_pythiaParticleFile.c_str() );
  cout << "Read particle momentum from PYTHIA data file." << endl;
  for ( int i = 0; i < numberOfParticlesToGenerate; i++ )
  {
    int flavTemp;
    // structure of file
    // number of pythia event  is hard?   flavour  energy  px  py  pz  m
    readPythiaParticles >> _particles[i].N_EVENT_pp >> _particles[i].HARD >> flavTemp >> _particles[i].E >> _particles[i].PX >> _particles[i].PY >> _particles[i].PZ >> _particles[i].m;
    _particles[i].FLAVOR = static_cast<FLAVOR_TYPE>( flavTemp );
    
    if ( flavTemp <= 2 * Particle::max_N_light_flavor )
    {
      _particles[i].m = 0;
    }
    
    // to avoid rounding errors compute energy from momenta and mass E^2=p^2+m^2, for light partons this can be different from Pythia's value if in Pythia this parton had a mass
    _particles[i].E = sqrt( pow( _particles[i].PX, 2.0 ) + pow( _particles[i].PY, 2.0 ) + pow( _particles[i].PZ, 2.0 ) + pow( _particles[i].m, 2.0 ) );
  }
  
  // charm quarks in Pythia have a mass of 1.5 GeV
  // make charm quarks from pythia lighter, if Mcharm is not 1.5 GeV
  if( !FPT_COMP_E( Particle::Mcharm, 1.5 ) )
  {
    const double M_old = 1.5; // charm mass in PYTHIA
    const double M_new = Particle::Mcharm;
    if(M_new > M_old)
    {
      cout << "problem in rhic::init(), charm mass to high. Leads to nan GeV energy." << endl;
    }
    cout << "Make charm quarks lighter. Old mass=" << M_old << "  new mass=" << M_new << endl;
    for( int j = 0; j < _particles.size(); j++ )
    {
      if( _particles[j].FLAVOR == charm || _particles[j].FLAVOR ==  anti_charm )
      {
        const double pp_old = sqrt ( pow ( _particles[j].PX, 2.0 ) + pow ( _particles[j].PY, 2.0 ) + pow ( _particles[j].PZ, 2.0 ) );
        const double pp_new = sqrt ( pow ( pp_old, 2.0 ) + pow ( M_old, 2.0 ) - pow ( M_new, 2.0 ) );

        // scaling in order to conserve the energy when making quarks massles
        _particles[j].PX = _particles[j].PX * pp_new / pp_old;
        _particles[j].PY = _particles[j].PY * pp_new / pp_old;
        _particles[j].PZ = _particles[j].PZ * pp_new / pp_old;
        _particles[j].E = sqrt ( pow ( _particles[j].PX, 2.0 ) + pow ( _particles[j].PY, 2.0 ) + pow ( _particles[j].PZ, 2.0 ) + pow ( M_new, 2.0 ) );
        _particles[j].m = M_new;
      }
    }
  }
  
  // make bottom quarks from pythia lighter, if Mbottom is not 4.8 GeV
  if( !FPT_COMP_E( Particle::Mbottom, 4.8 ) )
  {
    const double M_old = 4.8;  // bottom mass from pythia
    const double M_new = Particle::Mbottom; // our bottom mass
    if(M_new > M_old)
    {
      cout << "problem in rhic::init(), bottom mass to high. Leads to nan GeV energy." << endl;
    }
    cout << "Make bottom quarks lighter. Old mass=" << M_old << "  new mass=" << M_new << endl;
    for( int j = 0; j < _particles.size(); j++ )
    {
      if( _particles[j].FLAVOR == bottom || _particles[j].FLAVOR == anti_bottom )
      {
        const double pp_old = sqrt ( pow ( _particles[j].PX, 2.0 ) + pow ( _particles[j].PY, 2.0 ) + pow ( _particles[j].PZ, 2.0 ) );
        const double pp_new = sqrt ( pow ( pp_old, 2.0 ) + pow ( M_old, 2.0 ) - pow ( M_new, 2.0 ) );

        // scaling in order to conserve the energy when making quarks massles
        _particles[j].PX = _particles[j].PX * pp_new / pp_old;
        _particles[j].PY = _particles[j].PY * pp_new / pp_old;
        _particles[j].PZ = _particles[j].PZ * pp_new / pp_old;
        _particles[j].E = sqrt ( pow ( _particles[j].PX, 2.0 ) + pow ( _particles[j].PY, 2.0 ) + pow ( _particles[j].PZ, 2.0 ) + pow ( M_new, 2.0 ) );
        _particles[j].m = M_new;
      }
    }
  }
}




void initialModel_Pythia::samplePositions( std::vector< ParticleOffline >& _particles )
{
  cout << "Start sampling of particle positions for particles from PYTHIA." << endl;
  // particles are already read from file in momentum()
  int nmb_of_events = _particles.back().N_EVENT_pp;
  int event_tmp = 0;
  double x_tmp, y_tmp, z_tmp, t_tmp;

  bool soft_event[nmb_of_events + 1]; // array to store which events are soft and which not (true if event is soft)
  for ( int i = 1; i <= nmb_of_events; i++ )
  {
    soft_event[i] = false;
  }
  int partcl_soft_event[nmb_of_events + 1]; // array to store one particle number of each event to obtain the position of the event afterwards for soft particles

  // sample positions for binary collisions/events where hard partons are produced. Sample also if there are also soft particles at this event (reason: scaling behavior of soft and hard partons differ)
  for ( int j = 0; j < numberOfParticlesToGenerate; j++ )
  {
    if ( _particles[j].HARD )
    {
      if ( _particles[j].N_EVENT_pp != event_tmp ) // first hard particle of an event
      {
        sample_TXYZ_one_partcl( _particles[j], soft_event[_particles[j].N_EVENT_pp] );

        x_tmp = _particles[j].X;
        y_tmp = _particles[j].Y;
        z_tmp = _particles[j].Z;
        t_tmp = _particles[j].T;

        partcl_soft_event[_particles[j].N_EVENT_pp] = j;

        event_tmp = _particles[j].N_EVENT_pp;
      }
      else // all other hard particles of same event get same positions
      {
        _particles[j].X = x_tmp;
        _particles[j].Y = y_tmp;
        _particles[j].Z = z_tmp;
        _particles[j].T = t_tmp;
      }
    }
  }

  // start insert for non vanishing b=Bimp
  if ( impactParameter != 0.0 ) // soft particles are deleted for b!=0, because treatment of these is written for b=0 and would cause an error
  {
    std::vector<ParticleOffline>::iterator iIt;
    //!! delete particles
    for ( iIt = _particles.begin(); iIt != _particles.end(); )
    {
      if ( ( *iIt ).HARD == false )
      {
        iIt = _particles.erase( iIt );
      }
      else
      {
        ++iIt;
      }
    }
    cout << "Soft particles have been deleted since b!=0 and this would cause an error!!!" << endl;
  }
  // end of insert
  else
  {
    int sum = 0;
    int event = 1;
    // get positions for soft partons from events which have besides hard also soft partons
    int soft_break = 0;
    for ( int j = 0; j < numberOfParticlesToGenerate; j++ )
    {
      if ( !_particles[j].HARD )
      {
        sum++;
      }
    }
    cout <<  "Soft particles:" << sum << endl;

    sum = 0;
    for ( int i = 1; i <= nmb_of_events; i++ )
    {
      if ( soft_event[i] )
      {
        sum++;
      }
    }
    if ( sum == 0 )
    {
      string errMsg = "ERROR: No soft events! This will result in an infinite loop!";
      throw ePythia_error( errMsg );
    }
    cout << "Soft events: " << sum << endl;

    for ( int s = 0; s <  numberOfParticlesToGenerate; s++ )
    {
      if ( !_particles[s].HARD )
      {
        while ( !soft_event[event] )
        {
          if ( event < nmb_of_events )
          {
            event++;
          }
          else
          {
            soft_break = s;   // there are more soft particles than events for it
            event = 1;
          }
        }

        // event is soft
        _particles[s].X = _particles[partcl_soft_event[event]].X;
        _particles[s].Y = _particles[partcl_soft_event[event]].Y;
        _particles[s].Z = _particles[partcl_soft_event[event]].Z;
        _particles[s].T = _particles[partcl_soft_event[event]].T;
        if ( event < nmb_of_events )
        {
          event++;
        }
        else
        {
          soft_break = s;   // there are more soft particles than events for it
          event = 1;
        }
      }
    }
  }
  cout << "Finished sampling of particle positions for particles from PYTHIA." << endl;
}




// sample position for only one particle with id=number
void initialModel_Pythia::sample_TXYZ_one_partcl( ParticleOffline& _particle, bool& soft )
{
  double T, X, Y, Z;
  double densityA_max;
  double L_z;
  double p_soft;

  sample_TXYZ_singleParticle( _particle );

  // sample if there are also soft partons at this position:
  const double sigma = 40.0 * 0.1; // p+p cross section in mb, converted to fm^2  
  densityA_max = densityA( impactParameter / 2.0, WoodSaxonParameter.velocity * T );

  L_z = 2.0 / WoodSaxonParameter.gamma * sqrt( pow( WoodSaxonParameter.RA, 2.0 ) - pow( X, 2.0 ) - pow( Y, 2.0 ) ); // only valid for central collision b(=impactParameter)=0
  p_soft = 1.0 / (sigma * densityA_max * L_z);

//   if(p_soft > 1.0)
//   {
//     cout << "error, p_soft>1 in init_pos() (but doen't matter; so it's clearly a soft event), p_soft=" << p_soft << endl;
//   }

//   if(impactParameter != 0.0) // b != 0 would cause an error since consideration above is only for b=0
//     soft = false;
//   else

  soft = (ran2() < p_soft);

}







// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;
