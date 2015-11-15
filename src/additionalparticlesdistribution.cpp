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
#include "initialmodel_pythiashower.h"
#include "FPT_compare.h"

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
    {
      double usedMinimumPT;
      // if minimum PT of added particles is larger than the minijet cut off it is not necessary to sample particles below the value of the former because they are not active in the simulation anyhow.
      if( minimumPT > minijet_P0 )
        usedMinimumPT = minimumPT;
      else
        usedMinimumPT = minijet_P0;
      initialmodel = new initialModel_minijets( *configObject, _wsParameter, usedMinimumPT, numberOfParticlesToAdd );
      break;
    }
    case pythiaInitialState:
      initialmodel = new initialModel_Pythia( *configObject, _wsParameter, minimumPT, numberOfParticlesToAdd );
      break;
    case mcatnloInitialState:
      initialmodel = new initialModel_Mcatnlo( *configObject, _wsParameter, minimumPT, numberOfParticlesToAdd );
      break;
    case onlyJpsiInitialState:
      initialmodel = new initialModel_Jpsi( *configObject, _wsParameter );
      break;
    case fixedShowerInitialState:
      initialmodel = new initialModel_PYTHIAShower( *configObject, _wsParameter, fixedShower, configObject->getInitialPartonPt(), static_cast<FLAVOR_TYPE>( configObject->getInitialPartonFlavor() ) );
      break;
    case fixedPartonInitialState:
      initialmodel = new initialModel_PYTHIAShower( *configObject, _wsParameter, fixedParton, configObject->getInitialPartonPt(), static_cast<FLAVOR_TYPE>( configObject->getInitialPartonFlavor() ) );
      break;
    case pythiaShowerInitialState:
      initialmodel = new initialModel_PYTHIAShower( *configObject, _wsParameter, pythiaShower, minijet_P0 );
      break;
    case photonShowerInitialState:
      initialmodel = new initialModel_PYTHIAShower( *configObject, _wsParameter, photonShower, minijet_P0 );
      break;
    case charmShowerInititalState:
      initialmodel = new initialModel_PYTHIAShower( *configObject, _wsParameter, heavyQuarkShower, minijet_P0, charm );
      break;
    case bottomShowerInitialState:
      initialmodel = new initialModel_PYTHIAShower( *configObject, _wsParameter, heavyQuarkShower, minijet_P0, bottom );
      break;
    default:
      std::string errMsg = "Model for sampling the initial state not implemented yet!";
      throw eInitialModel_error( errMsg );
      break;
  }
  
  std::vector<Particle> tempParticleVector;
  initialmodel->populateParticleVector( tempParticleVector );
  
  if( Particle::N_psi_states > 0 && initialStateType != onlyJpsiInitialState ) // for onlyJpsiInitialState Jpsi has been sampled above
  {
    initialModel_Jpsi theIni_Jpsi( *configObject, _wsParameter );
    theIni_Jpsi.populateParticleVector( tempParticleVector );
  }
  
  _particles.reserve( tempParticleVector.size() );
  for ( unsigned int i = 0; i < tempParticleVector.size(); i++ )
  {
    ParticleOffline tempParticle( tempParticleVector[i] );
    _particles.push_back( tempParticle );
  }
  
  if( configObject->isStudyNonPromptJpsiInsteadOfElectrons() )
    deleteAllParticlesExceptBottom( _particles );

  for ( unsigned int i = 0; i < _particles.size(); i++ )
  {
    _particles[i].unique_id = ParticleOffline::unique_id_counter_added;
    --ParticleOffline::unique_id_counter_added;
  }
}



void additionalParticlesDistribution::prepareParticles( std::vector< ParticleOffline >& _particles )
{
  double dtt = 0;
  double eta_max = 5.0;

  double MT, y;
  double shift;

  for ( unsigned int j = 0; j < _particles.size(); j++ )
  {   
    dtt = fabs( _particles[j].Pos.Z() ) / tanh( eta_max ) - _particles[j].Pos.T();  //tanh(eta)=z/t
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

    _particles[j].Pos.T() += shift;
    
    // last interaction spacetime point is point of creation
    _particles[j].lastInt = _particles[j].Pos;
    
    //formation time 1/sqrt( p_T^2 + m^2) = 1/m_T
    y = _particles[j].Mom.Rapidity();
    MT = _particles[j].Mom.Mt( _particles[j].m );
    dtt = 1 / MT * cosh( y ) * 0.197327;  //fm/c   //cosh(y) = gamma  (of that particle wrt motion in z-direction)
    
    // additional formation time for Jpsi
    if( _particles[j].FLAVOR == jpsi )
      dtt += configObject->getJpsiFormationTime() * cosh( y ); //fm/c   //cosh(y) = gamma  (of that particle wrt motion in z-direction)

    _particles[j].Propagate( _particles[j].Pos.T() + dtt );

    _particles[j].init = true;
    
    // creation time is time at which the particle is allowed to scatter
    _particles[j].T_creation = _particles[j].Pos.T();
  }
  
  for(unsigned int i = 0; i < _particles.size(); i++)
  {
    _particles[i].PosInit = _particles[i].Pos;
    _particles[i].MomInit = _particles[i].Mom;
    _particles[i].X_traveled = 0.0;
  }
}


void additionalParticlesDistribution::deleteAllParticlesExceptBottom( std::vector< ParticleOffline >& _particles )
{
  for(unsigned int j = 0; j < addedParticles.size(); j++ )
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