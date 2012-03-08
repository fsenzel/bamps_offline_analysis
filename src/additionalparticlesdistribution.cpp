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

#include "additionalparticlesdistribution.h"
#include "configuration.h"
#include "initialmodel.h"
#include "initialmodel_mcatnlo.h"
#include "initialmodel_minijets.h"
#include "initialmodel_pythia.h"
#include "initialmodel_cgc.h"
#include "initialmodel_jpsi.h"

using namespace ns_casc;
using namespace std;


additionalParticlesDistribution::additionalParticlesDistribution( const config* const _config, const INITIAL_STATE_TYPE _initialStateType ) :
    initialStateType( _initialStateType ),
    configObject( _config ),
    numberOfParticlesToAdd( _config->getNumberOfParticlesToAdd() ),
    minimumPT( _config->getMinimumPT() ),
    minijet_P0( _config->getPtCutoff() ),
    impactParameter( _config->getImpactParameter() ),
    numberOfTestparticles( _config->getTestparticles() )
{
}



void additionalParticlesDistribution::populateParticleVector( std::vector< ParticleOffline >& _particles, WoodSaxon& _wsParameter )
{
  initialModel *initialmodel;
  
  switch ( initialStateType )
  {
    case miniJetsInitialState:
      double usedMinimumPT;
      // if minimum PT of added particles is larger than the minijet cut off it is not necessary to sample particles below the value of the former because they are not active in the simulation anyhow.
      if( minimumPT > minijet_P0 )
        usedMinimumPT = minimumPT;
      else
        usedMinimumPT = minijet_P0;
      initialmodel = new initialModel_minijets( *configObject, _wsParameter, usedMinimumPT, numberOfParticlesToAdd );
      break;
    case pythiaInitialState:
      initialmodel = new initialModel_Pythia( *configObject, _wsParameter, minimumPT, numberOfParticlesToAdd );
      break;
    case mcatnloInitialState:
      initialmodel = new initialModel_Mcatnlo( *configObject, _wsParameter, minimumPT, numberOfParticlesToAdd );
      break;
    default:
      std::string errMsg = "Model for sampling the initial state not implemented yet!";
      throw eInitialModel_error( errMsg );
      break;
  }
  
  std::vector<Particle> tempParticleVector;
  initialmodel->populateParticleVector( tempParticleVector );
  
  if( Particle::N_psi_states > 0 )
  {
    initialModel_Jpsi theIni_Jpsi( *configObject, _wsParameter );
    theIni_Jpsi.populateParticleVector( tempParticleVector );
  }
  
  _particles.reserve( tempParticleVector.size() );
  for ( int i = 0; i < tempParticleVector.size(); i++ )
  {
    ParticleOffline tempParticle( tempParticleVector[i] );
    _particles.push_back( tempParticle );
  }
  
  if( configObject->isStudyNonPromptJpsiInsteadOfElectrons() )
    deleteAllParticlesExceptBottom( _particles );

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

  double MT, y, cc;
  double shift;
  for ( int j = 0; j < _particles.size(); j++ )
  {   
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
    
    // last interaction spacetime point is point of creation
    _particles[j].X_lastInt = _particles[j].X;
    _particles[j].Y_lastInt = _particles[j].Y;
    _particles[j].Z_lastInt = _particles[j].Z;
    _particles[j].T_lastInt = _particles[j].T;
    
    //formation time 1/sqrt( p_T^2 + m^2) = 1/m_T
    y = 0.5 * log(( _particles[j].E+_particles[j].PZ ) / ( _particles[j].E-_particles[j].PZ ) );
    MT = sqrt( pow( _particles[j].PX, 2 ) + pow( _particles[j].PY, 2 ) + pow( _particles[j].m, 2 ) );   
    dtt = 1 / MT * cosh( y ) * 0.197327;  //fm/c   //cosh(y) = gamma  (of that particle wrt motion in z-direction)
    
    // additional formation time for Jpsi
    if( _particles[j].FLAVOR == jpsi )
      dtt += configObject->getJpsiFormationTime() * cosh( y ); //fm/c   //cosh(y) = gamma  (of that particle wrt motion in z-direction)

    cc = dtt / _particles[j].E;
    _particles[j].T = _particles[j].T + dtt;
    _particles[j].X = _particles[j].X + _particles[j].PX * cc;
    _particles[j].Y = _particles[j].Y + _particles[j].PY * cc;
    _particles[j].Z = _particles[j].Z + _particles[j].PZ * cc;

    _particles[j].init = true;
    
    // creation time is time at which the particle is allowed to scatter
    _particles[j].T_creation = _particles[j].T;
  }
  
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


void additionalParticlesDistribution::deleteAllParticlesExceptBottom( std::vector< ParticleOffline >& _particles )
{
  for(int j = 0; j < addedParticles.size(); j++ )
  {
    if( !( addedParticles[j].FLAVOR == bottom || addedParticles[j].FLAVOR == anti_bottom ) )
    {
      // delete last particle if also not active otherwise switch position with particle to be deleted
      while( !( addedParticles.back().FLAVOR == bottom || addedParticles.back().FLAVOR == anti_bottom ) && 
            ( j != addedParticles.size() - 1 ) ) // if particle j is the last particle in the particle list it is deleted here and the then last in the list below as well, which is not correct.
      {
        addedParticles.pop_back();
      }
      addedParticles[j] = addedParticles.back();
      addedParticles.pop_back();
    }
  }
}
