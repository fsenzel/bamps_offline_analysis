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

    INITIAL_STATE_TYPE initialStateType;

    const config * const configObject;

    int numberOfParticlesToAdd;

    /** @brief Minimum p_T [GeV] of the added particles */
    double minimumPT;

    /** @brief Lower PT-cutoff [GeV] used for minijet initial conditions */
    double minijet_P0;

    /** @brief impact parameter in fm */
    double impactParameter;

    /** @brief number of testparticles per real particle */
    int numberOfTestparticles;


};


/** @brief exception class for handling unexpected critical behaviour within generation of mini-jet initial distributions  */
class eInitialState_error : public std::runtime_error
{
public:
    explicit eInitialState_error( const std::string& what ) : std::runtime_error( what ) {};

    virtual ~eInitialState_error() throw() {};
};

#endif // INITIALDISTRIBUTION_H
