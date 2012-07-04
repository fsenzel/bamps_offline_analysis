//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------
//---------------------------------------------



#ifndef INITIALDISTRIBUTION_H
#define INITIALDISTRIBUTION_H

#include <stdexcept>
#include <vector>

#include "configuration.h"
#include "particle.h"
#include "woodsaxon.h"



class additionalParticlesDistribution
{
public:
  additionalParticlesDistribution( const config*const _config, const INITIAL_STATE_TYPE _initialStateType = miniJetsInitialState );
  ~additionalParticlesDistribution() {};

  void populateParticleVector( std::vector< ParticleOffline >& _particles, WoodSaxon& _wsParameter );
  void prepareParticles( std::vector< ParticleOffline >& _particles );
 

private:
  
  void deleteAllParticlesExceptBottom( std::vector< ParticleOffline >& _particles );
  
  int numberOfParticlesToAdd;
  /** @brief Minimum p_T [GeV] of the added particles */
  double minimumPT;
  /** @brief Lower PT-cutoff [GeV] used for minijet initial conditions */
  double minijet_P0;
  const config * const configObject;
  double impactParameter;       //impact parameter in fm
  int numberOfTestparticles;    //number of testparticles per real particle
  INITIAL_STATE_TYPE initialStateType;
  
  /** @brief Cut-off time for shower evolution [GeV^-1] */
  double insertionTime;
  
  /** @brief Random number generator seed for fixing PYTHIA seed */
  long int seed;
  
  /** @brief Filename prefix, needed for initial unshowered particle output */
  string filename_prefix; 
  
  /** @brief Routine for converting a particle vector consisting of back-to-back
   *         parton pairs into dijet showers */
  void showerParticles( std::vector< ParticleOffline >& _particles ); 
  
  /** @brief Routine for setting event flag for additional particles */
  void setEventID( std::vector< ParticleOffline >& _particles );
  
  /** @brief Output routine of shower initiating partons before showering */
  void initialShowerInitOutput( std::vector< ParticleOffline > _particles );
  
  /** @brief Routine for creating one dijet shower out of PYTHIA */
  vector<ParticleOffline> createShowerEvent( const double _px, const double _py, const double _pz1, const double _pz2, const FLAVOR_TYPE _flavor1, const FLAVOR_TYPE _flavor2 );
    
};


/** @brief exception class for handling unexpected critical behaviour within generation of mini-jet initial distributions  */
class eInitialState_error : public std::runtime_error
{
public:
  explicit eInitialState_error( const std::string& what ) : std::runtime_error( what ) {};

  virtual ~eInitialState_error() throw() {};
};

#endif // INITIALDISTRIBUTION_H
