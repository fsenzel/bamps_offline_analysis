//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------
//---------------------------------------------


#include <string>
#include <iostream>

#include "minijets.h"
#include "pythiaInitialDistribution.h"
#include "cgcInitialDistribution.h"
#include "initialdistribution.h"
#include "woodsaxon.h"
#include "configuration.h"

using namespace ns_casc;


initialDistribution::initialDistribution( const config* const _config, const INITIAL_STATE_TYPE _initialStateType ) :
    initialStateType( _initialStateType ),
    configObject( _config )
{
  A = _config->getA();
  Aatomic = _config->getAatomic();
  B = _config->getB();
  Batomic = _config->getBatomic();
  sqrtS = _config->getSqrtS();
  impactParameter = _config->getImpactParameter();
  numberOfTestparticles = _config->getTestparticles();
}



void initialDistribution::populateParticleVector( std::vector< ParticleOffline >& _particles, WoodSaxon& _wsParameter )
{
  switch ( initialStateType )
  {
    case miniJetsInitialState:
    {
      miniJets _miniJets( *configObject, _wsParameter, useStoredTables );
      WoodSaxonParameter = _wsParameter;
      _miniJets.populateParticleVector( _particles );
      break;
    }
    case pythiaInitialState:
    {
      pythiaInitialDistribution _pythiaInitialDistribution( *configObject, _wsParameter, useStoredTables );
      WoodSaxonParameter = _wsParameter;
      _pythiaInitialDistribution.populateParticleVector( _particles );
      break;
    }
    case cgcInitialState:
    {
      cgcInitialDistribution _cgcInitialDistribution( *configObject);
      _cgcInitialDistribution.populateParticleVector( _particles );
      break;
    }
    default:
      std::string errMsg = "Model for sampling the initial state not implemented yet!";
      throw eInitialState_error( errMsg );
      break;
  }
  
  ParticleOffline::unique_id_counter = 0;
  for ( int i = 0; i < _particles.size(); i++ )
  {
    _particles[i].unique_id = ParticleOffline::unique_id_counter;
    ++ParticleOffline::unique_id_counter;
  }
}