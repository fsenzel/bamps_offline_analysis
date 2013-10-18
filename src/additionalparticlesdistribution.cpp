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
#include "FPT_compare.h"

using namespace ns_casc;
using namespace std;

extern "C" {
  void shower_( const double* px, const double* py, const double* pz1, const double* pz2, int* pythiaFlavor1, int* pythiaFlavor2, double* tauf, uint32_t* seed );
  struct
  {
    double pa[100][6];
  } bamps_;
}

additionalParticlesDistribution::additionalParticlesDistribution( const config* const _config, const INITIAL_STATE_TYPE _initialStateType ) :
    initialStateType( _initialStateType ),
    configObject( _config ),
    numberOfParticlesToAdd( _config->getNumberOfParticlesToAdd() ),
    minimumPT( _config->getMinimumPT() ),
    minijet_P0( _config->getPtCutoff() ),
    impactParameter( _config->getImpactParameter() ),
    numberOfTestparticles( _config->getTestparticles() ),
    insertionTime( _config->getInsertionTime() ),
    seed( _config->getSeed() ),
    filename_prefix( _config->getStandardOutputDirectoryName() + "/" + _config->getJobName() ),
    initialPartonPt( _config->getInitialPartonPt() )
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
    case showerInitialState:
    {
      double usedMinimumPT;
      // if minimum PT of added particles is larger than the minijet cut off it is not necessary to sample particles below the value of the former because they are not active in the simulation anyhow.
      if( minimumPT > minijet_P0 )
        usedMinimumPT = minimumPT;
      else
        usedMinimumPT = minijet_P0;
      if( numberOfParticlesToAdd % 2 != 0 )
      {
        string errMsg = "For shower initial conditions an even number of initial particles is needed!";
        throw eInitialState_error( errMsg );
      }

      initialmodel = new initialModel_minijets( *configObject, _wsParameter, usedMinimumPT, numberOfParticlesToAdd );
      break;
    }
    case fixedPartonState:
    {
      double usedMinimumPT;
      // if minimum PT of added particles is larger than the minijet cut off it is not necessary to sample particles below the value of the former because they are not active in the simulation anyhow.
      if( minimumPT > minijet_P0 )
        usedMinimumPT = minimumPT;
      else
        usedMinimumPT = minijet_P0;
      if( numberOfParticlesToAdd % 2 != 0 )
      {
        string errMsg = "For shower initial conditions an even number of initial particles is needed!";
        throw eInitialState_error( errMsg );
      }

      initialmodel = new initialModel_minijets( *configObject, _wsParameter, usedMinimumPT, numberOfParticlesToAdd );
      break;
    }
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

  //     Showering of particles, if initialStateType = 5, unsure if before preparing or after...
  if( initialStateType == showerInitialState )
  {
    setEventID( _particles );
    initialShowerInitOutput( _particles );
    showerParticles( _particles );
  }
  else if( initialStateType == fixedPartonState )
  {
    for( unsigned int i = 0; i < _particles.size(); i += 2 )
    {
      _particles[i].Mom.Px() = initialPartonPt;
      _particles[i + 1].Mom.Px() = -initialPartonPt;
      _particles[i].Mom.Py() = _particles[i + 1].Mom.Py() = 0.0;
      _particles[i].Mom.Pz() = _particles[i + 1].Mom.Pz() = 0.0;
      _particles[i].Mom.E()  = _particles[i + 1].Mom.E()  = initialPartonPt;
      _particles[i].FLAVOR   = _particles[i + 1].FLAVOR   = gluon;
    }
    setEventID( _particles );
    initialShowerInitOutput( _particles );
  }
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


void additionalParticlesDistribution::showerParticles( vector< ParticleOffline >& _particles )
{
  vector<ParticleOffline> tempParticles;

  if( insertionTime != 0.0 )
  {
    for( int i = 0; i < _particles.size(); i += 2 )
    {
      vector<ParticleOffline> particleShower;
      double px = _particles[i].Mom.Px();
      double py = _particles[i].Mom.Py();
      double pz1 = _particles[i].Mom.Pz();
      double pz2 = _particles[i + 1].Mom.Pz();
      FLAVOR_TYPE flavor1 = _particles[i].FLAVOR;
      FLAVOR_TYPE flavor2 = _particles[i + 1].FLAVOR;
      particleShower = createShowerEvent( px, py, pz1, pz2, flavor1, flavor2 );

      for( int j = 0; j < particleShower.size(); j++ )
      {
        particleShower[j].N_EVENT_pp = _particles[i].N_EVENT_pp;
        particleShower[j].N_EVENT_AA = _particles[i].N_EVENT_AA;
        particleShower[j].Pos     = _particles[i].Pos;
        particleShower[j].PosInit = particleShower[j].Pos;
        particleShower[j].MomInit = particleShower[j].Mom;

        particleShower[j].init = true;

        tempParticles.push_back( particleShower[j] );
      }
    }
    _particles.clear();
    _particles = tempParticles;
    cout << "#### " << "After showering for " << insertionTime * 0.198 << " fm, " << _particles.size()
         << " showered particles out of " << numberOfParticlesToAdd << " particles added..." << endl;
  }
}


vector<ParticleOffline> additionalParticlesDistribution::createShowerEvent( const double _px, const double _py, const double _pz1, const double _pz2, const FLAVOR_TYPE _flavor1, const FLAVOR_TYPE _flavor2 )
{

// PYTHIA seed has max value 9E8
  while( seed > 900000000 )
    seed -= 900000000;

  vector<ParticleOffline> particlesToAdd;
  int flavor1, flavor2;

  switch( _flavor1 )
  {
  case down:
    flavor1 = 1;
    break;
  case anti_down:
    flavor1 = -1;
    break;
  case up:
    flavor1 = 2;
    break;
  case anti_up:
    flavor1 = -2;
    break;
  case strange:
    flavor1 = 3;
    break;
  case anti_strange:
    flavor1 = -3;
    break;
  case gluon:
    flavor1 = 21;
    break;
  default:
    cout << "Unknown flavor type...." << endl;
  }

  switch( _flavor2 )
  {
  case down:
    flavor2 = 1;
    break;
  case anti_down:
    flavor2 = -1;
    break;
  case up:
    flavor2 = 2;
    break;
  case anti_up:
    flavor2 = -2;
    break;
  case strange:
    flavor2 = 3;
    break;
  case anti_strange:
    flavor2 = -3;
    break;
  case gluon:
    flavor2 = 21;
    break;
  default:
    cout << "Unknown flavor type...." << endl;
  }

  int attempt = 0;
  do
  {
    particlesToAdd.clear();

    shower_( &_px, &_py, &_pz1, &_pz2, &flavor1, &flavor2, &insertionTime, &seed );

    int index = 0;
    while( bamps_.pa[index][0] != 0 )
    {
      ParticleOffline tempParticle;
      switch( static_cast<int>( bamps_.pa[index][0] ) )
      {
      case 21:
        tempParticle.FLAVOR = gluon;
        break;
      case 2:
        tempParticle.FLAVOR = up;
        break;
      case 1:
        tempParticle.FLAVOR = down;
        break;
      case -2:
        tempParticle.FLAVOR = anti_up;
        break;
      case -1:
        tempParticle.FLAVOR = anti_down;
        break;
      case 3:
        tempParticle.FLAVOR = strange;
        break;
      case -3:
        tempParticle.FLAVOR = anti_strange;
        break;
      default:
        cout << "Unknown flavor type:\t" << bamps_.pa[index][0] << endl;
      }

      tempParticle.Mom.Px() = bamps_.pa[index][1];
      tempParticle.Mom.Py() = bamps_.pa[index][2];
      tempParticle.Mom.Pz() = bamps_.pa[index][3];
      tempParticle.Mom.E() = sqrt( tempParticle.Mom.vec2() );
      tempParticle.m = 0.0;
      particlesToAdd.push_back( tempParticle );
      index++;
    }

    attempt++;
  }
  while( particlesToAdd.size() == 0 );

  if( attempt > 10 )
    cout << attempt << " Attempts needed to get an allowed shower." << endl;

  double sumE = 0.0;
  double E1 = sqrt( _px * _px + _py * _py + _pz1 * _pz1 );
  double E2 = sqrt( _px * _px + _py * _py + _pz2 * _pz2 );
  for( int i = 0; i < particlesToAdd.size(); i++ )
  {
    sumE += particlesToAdd[i].Mom.E();
  }
  if( FPT_COMP_GE( abs( sumE - ( E1 + E2 ) ) / sumE, 0.05 ) )
  {
    stringstream errMsg;
    errMsg << "Total energy of shower particles differs from energy of shower-initiating partons more than 5%:\t"
           << sumE << "\t" << E1 + E2 << "Unrecoverable error!" << endl;
    throw eInitialState_error( errMsg.str() );
  }

  return particlesToAdd;
}

void additionalParticlesDistribution::setEventID( vector< ParticleOffline >& _particles )
{
  for( int index = 0; index < _particles.size(); index += 2 )
  {
    _particles[index].N_EVENT_pp = _particles[index + 1].N_EVENT_pp = static_cast<int>( index / 2 );
    _particles[index].N_EVENT_AA = _particles[index + 1].N_EVENT_AA = static_cast<int>( index / 2 );
  }
}

void additionalParticlesDistribution::initialShowerInitOutput( vector< ParticleOffline > _particles )
{
  time_t end;
  time( &end );

  string filename = filename_prefix + "_" + "unshoweredParticles";
  fstream file( filename.c_str(), ios::out | ios::trunc );

  //---- print header if file is empty ----
  /*  file.seekp( 0, ios::end );
    long size = file.tellp();
    file.seekp( 0, ios::beg );
    if ( size == 0 )
      printHeader( file, rapidityDistribution, end );*/
  //---------------------------------------
  string sep = "\t";

  file << "#event" << sep << "px" << sep << "py" << sep << "pz" << sep << "E" << sep << "x"
       << sep << "y" << sep << "z" << sep << "t" << endl;
  for( int index = 0; index < _particles.size(); index++ )
  {
    file << _particles[index].N_EVENT_pp << sep 
         << _particles[index].Mom.Px() << sep 
         << _particles[index].Mom.Py() << sep 
         << _particles[index].Mom.Pz() << sep
         << _particles[index].Mom.E() << sep 
         << _particles[index].Pos.X() << sep 
         << _particles[index].Pos.Y() << sep 
         << _particles[index].Pos.Z() << sep
         << _particles[index].Pos.T() << endl;
  }
  file.close();

}


