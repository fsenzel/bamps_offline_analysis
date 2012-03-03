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

#include "configuration.h"
#include "initialmodel_particlesFromFile.h"
#include "particle.h"
#include "FPT_compare.h"

using std::cout;
using std::endl;
using namespace ns_casc;

/**
* @brief Constructor for initial conditions for particles from file. Filename must be speciefied in derived class.
* @param[in] _config configuration class
* @param[out] _WoodSaxonParameter returns Wood Saxon Parameters
* @param[in] _minimumPT Particles with transverse momenta below this value are discared and not considered for simulation.
* @param[in] _nToGenerate Specifies number of particles to generate. A negative value is standard and indicates that all particles in the particle file are used.
*/
initialModel_ParticlesFromFile::initialModel_ParticlesFromFile( const string _filename_particleFile, const config& _config, WoodSaxon& _WoodSaxonParameter, const double _minimumPT, const int _nToGenerate ) :
  filename_particleFile(_filename_particleFile),
  initialModelWS(_config) ,
  numberOfParticlesToGenerate ( 0 ),
  numberOfTestparticles ( _config.getTestparticles() ),
  minimumPT( _minimumPT ),
  nEventsAA( _config.getNaddedEvents() )
{
  double Tab;

  impactParameter = _config.getImpactParameter();
  sqrtS_perNN = _config.getSqrtS();
  if (!WoodSaxonParameter.Calculate( A, impactParameter, sqrtS_perNN))
  {
    std::string errMsg = "Impact parameter b too large. b > 2 R_A0";
    throw eParticlesFromFile_error( errMsg );
  }
  _WoodSaxonParameter = WoodSaxonParameter;


  cout << "======= Generating data sets for sampling of initial state =======" << endl;
  generateTimeDistributionWS(Tab);
  cout << "++++  Tab = " << Tab << "1/mb" << endl;
  cout << "==================================================================" << endl;

  cout << "Particle data file: " << filename_particleFile << endl;

  std::ifstream countParticles( filename_particleFile.c_str() );
  if ( countParticles.good() )
  {
    int numberOfParticlesInDataFile = std::count( std::istreambuf_iterator<char>( countParticles ), std::istreambuf_iterator<char>(), '\n' ); // counts number of lines in particle data file
    
    if( _nToGenerate < 0 )
    {  
      numberOfParticlesToGenerate = numberOfParticlesInDataFile;
    }
    else
    {
      if( numberOfParticlesInDataFile < _nToGenerate )
      {
        cout << "Not enough particles in data file to sample all particles. Sample only " << numberOfParticlesInDataFile << " instead of " << _nToGenerate << endl;
        numberOfParticlesToGenerate = numberOfParticlesInDataFile;
      }
      else
      {
        numberOfParticlesToGenerate = _nToGenerate;
      }
    }
  }
  else
  {
    string errMsg = "Error at opening particle data file.";
    throw eParticlesFromFile_error( errMsg );
  }

  countParticles.close();
}



void initialModel_ParticlesFromFile::populateParticleVector( std::vector< Particle >& _particles )
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

  sampleMomenta( _particles );
  samplePositions( _particles );
}

void initialModel_ParticlesFromFile::samplePositions( std::vector< Particle >& _particles )
{
  cout << "Start sampling of particle positions for particles from file." << endl;
  // particles are already read from file in momentum()
  int nmb_of_events;
  if( _particles.size() > 0 )
    nmb_of_events = _particles.back().N_EVENT_pp;
  else
    nmb_of_events = 0;
  int event_tmp = 0;
  double x_tmp, y_tmp, z_tmp, t_tmp, event_AA_tmp;

  bool soft_event[nmb_of_events + 1]; // array to store which events are soft and which not (true if event is soft)
  for ( int i = 1; i <= nmb_of_events; i++ )
  {
    soft_event[i] = false;
  }
  int partcl_soft_event[nmb_of_events + 1]; // array to store one particle number of each event to obtain the position of the event afterwards for soft particles

  // sample positions for binary collisions/events where hard partons are produced. Sample also if there are also soft particles at this event (reason: scaling behavior of soft and hard partons differ)
  for ( int j = 0; j < _particles.size(); j++ )
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

        // pp event to which the particle belongs
        event_tmp = _particles[j].N_EVENT_pp;
        
        // sample heavy ion event to which the particle belongs
        event_AA_tmp = int( ran2() * nEventsAA ) + 1; // event_AA_tmp is integer in [1;nEventsAA]
        if(event_AA_tmp > nEventsAA)
          cout << " error in initialModel_ParticlesFromFile::samplePositions, event_AA_tmp=" << event_AA_tmp << endl;
        _particles[j].N_EVENT_AA = event_AA_tmp; // give particles the flag of the heavy ion event
      }
      else // all other hard particles of same event get same positions
      {
        _particles[j].X = x_tmp;
        _particles[j].Y = y_tmp;
        _particles[j].Z = z_tmp;
        _particles[j].T = t_tmp;
        
        _particles[j].N_EVENT_AA = event_AA_tmp; // give particles the flag of the heavy ion event
      }
    }
  }

  // start insert for non vanishing b=Bimp
  if ( impactParameter != 0.0 ) // soft particles are deleted for b!=0, because treatment of these is written for b=0 and would cause an error
  {
    // delete soft particles
    for(int j = 0; j < _particles.size(); j++ )
    {
      // delete last particle if soft
      if( !_particles[j].HARD )
      {
        // delete last particle if also not active otherwise switch position with particle to be deleted
        while( ( !_particles.back().HARD ) && ( j != _particles.size() - 1  ) ) // if particle j is the last particle in the particle list it is deleted here and the then last in the list below as well, which would be wrong.
        {
          _particles.pop_back();
        }
        _particles[j] = _particles.back();
        _particles.pop_back();
      }
    }
    cout << "Soft particles have been deleted since b!=0 and this would cause an error!" << endl;
  }
  // end of insert
  else
  {
    int sum = 0;
    int event = 1;
    // get positions for soft partons from events which have besides hard also soft partons
    int soft_break = 0;
    for ( int j = 0; j < _particles.size(); j++ )
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
      throw eParticlesFromFile_error( errMsg );
    }
    cout << "Soft events: " << sum << endl;

    for ( int s = 0; s < _particles.size(); s++ )
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
}



// sample position for only one particle with id=number
void initialModel_ParticlesFromFile::sample_TXYZ_one_partcl( Particle& _particle, bool& soft )
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


void initialModel_ParticlesFromFile::changeCharmMass( std::vector< Particle >& _particles, const double M_old )
{
  const double M_new = Particle::Mcharm;
  
  if( !FPT_COMP_E( M_new, M_old ) )
  {
    if(M_new > M_old)
    {
      std::string errMsg = "problem in initialModel_ParticlesFromFile::changeCharmMass, charm mass to high. Leads to nan GeV energy.";
      throw eParticlesFromFile_error( errMsg );
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
}

void initialModel_ParticlesFromFile::changeBottomMass( std::vector< Particle >& _particles, const double M_old )
{
  const double M_new = Particle::Mbottom; // our bottom mass
  
  if( !FPT_COMP_E( M_new, M_old ) )
  {
    if(M_new > M_old)
    {
      cout << M_new << "  " << M_old << endl;
      std::string errMsg = "problem in initialModel_ParticlesFromFile::changeBottomMass, bottom mass to high. Leads to nan GeV energy.";
      throw eParticlesFromFile_error( errMsg );
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


