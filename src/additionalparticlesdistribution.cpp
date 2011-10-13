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
#include <math.h>
#include <vector>

#include "minijets.h"
#include "additionalparticlesdistribution.h"
#include "pythiaInitialDistribution.h"
#include "woodsaxon.h"
#include "configuration.h"

using namespace ns_casc;
using namespace std;


additionalParticlesDistribution::additionalParticlesDistribution( const config* const _config, const INITIAL_STATE_TYPE _initialStateType ) :
    initialStateType( _initialStateType ),
    configObject( _config ),
    numberOfParticlesToAdd( 0 )
{
  numberOfParticlesToAdd = configObject->getNumberOfParticlesToAdd();
  minimumPT = configObject->getMinimumPT();
  
  A = _config->getA();
  Aatomic = _config->getAatomic();
  B = _config->getB();
  Batomic = _config->getBatomic();
  sqrtS = _config->getSqrtS();
  impactParameter = _config->getImpactParameter();
  numberOfTestparticles = _config->getTestparticles();
}



void additionalParticlesDistribution::populateParticleVector( std::vector< ParticleOffline >& _particles, WoodSaxon& _wsParameter )
{
  switch ( initialStateType )
  {
  case miniJetsInitialState:
  {
    miniJets _miniJets( *configObject, _wsParameter, useStoredTables );
    WoodSaxonParameter = _wsParameter;
    _miniJets.populateParticleVector( _particles, numberOfParticlesToAdd, minimumPT );
    break;
  }
  case pythiaInitialState:
  {
    pythiaInitialDistribution _pythiaInitialDistribution( *configObject, _wsParameter, useStoredTables );
    WoodSaxonParameter = _wsParameter;
    _pythiaInitialDistribution.populateParticleVector( _particles, numberOfParticlesToAdd, minimumPT ); 
    break;
  }
  default:
    std::string errMsg = "Model for sampling the initial state not implemented yet!";
    throw eInitialState_error( errMsg );
    break;
  }

  for ( int i = 0; i < _particles.size(); i++ )
  {
    _particles[i].unique_id = ParticleOffline::unique_id_counter_added;
    --ParticleOffline::unique_id_counter_added;
  }
}



void additionalParticlesDistribution::prepareParticles( std::vector< ParticleOffline >& _particles )
{
  double max = 0;
  double dtt = 0;
  double eta_max = 5.0;

  double PT, y, cc;
  double shift;
  double sum = 0;
  for ( int j = 0; j < _particles.size(); j++ )
  {   
    // last interaction spacetime point is point of creation
    _particles[j].X_lastInt = _particles[j].X;
    _particles[j].Y_lastInt = _particles[j].Y;
    _particles[j].Z_lastInt = _particles[j].Z;
    _particles[j].T_lastInt = _particles[j].T;
    _particles[j].T_creation = _particles[j].T;
    
    
    dtt = fabs( _particles[j].Z ) / tanh( eta_max ) - _particles[j].T;  //tanh(eta)=z/t
    if ( dtt < configObject->getTimeshift() )
    {
      shift = configObject->getTimeshift();
    }
    else
    {
      // all particles should be shifted by dtt, but we cannot do that for gluons from cascade. Therefore we just shift that particle and hope that dtt-timeshiftGluons is negligible small...
      shift = dtt;
      cout << "timeshift larger than in original BAMPS run: " << dtt << "  difference: " <<  dtt - configObject->getTimeshift() << endl;
      cout << "shift only that particle by that timeshift." << endl;
    }

    _particles[j].T += shift;
    
    //formation time 1/p_T
    y = 0.5 * log(( _particles[j].E+_particles[j].PZ ) / ( _particles[j].E-_particles[j].PZ ) );
    PT = sqrt( pow( _particles[j].PX, 2 ) + pow( _particles[j].PY, 2 ) );
    dtt = 1 / PT * cosh( y ) * 0.197327;  //fm/c   //cosh(y) = gamma  (of that particle wrt motion in z-direction)

    cc = dtt / _particles[j].E;
    _particles[j].T = _particles[j].T + dtt;
    _particles[j].X = _particles[j].X + _particles[j].PX * cc;
    _particles[j].Y = _particles[j].Y + _particles[j].PY * cc;
    _particles[j].Z = _particles[j].Z + _particles[j].PZ * cc;

    _particles[j].init = true;
    
    sum += _particles[j].E;
  }
 
  cout << "** " << _particles.size() << " particles added with total energy: " << sum << endl;
  
  for(int i = 0; i < _particles.size(); i++)
  {
    _particles[i].X_init = _particles[i].X;
    _particles[i].Y_init = _particles[i].Y;
    _particles[i].Z_init = _particles[i].Z;
    
    _particles[i].E_init = _particles[i].E;
    _particles[i].PX_init = _particles[i].PX;
    _particles[i].PY_init = _particles[i].PY;
    _particles[i].PZ_init = _particles[i].PZ;
    
    _particles[i].X_traveled = 0.0;
  }
}
// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on; 
