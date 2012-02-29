//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------
//---------------------------------------------

#ifndef INITIALMODEL_PYTHIA_H
#define INITIALMODEL_PYTHIA_H

#include <stdexcept>
#include <vector>
#include <string>

#include "initialmodel.h"
#include "configuration.h"
#include "particle.h"
#include "vegas.h"
#include "woodsaxon.h"


/**
 * @brief Class to provide the Pythia initialization
 */
class initialModel_Pythia : public initialModelWS
{
  public:
    initialModel_Pythia( const config& _config, WoodSaxon& _WoodSaxonParameter, const double _minimumPT = 0.0, const int _nToGenerate = -1 );
    ~initialModel_Pythia() {};
    
    void populateParticleVector( std::vector<Particle>& _particles );
    
        
  protected:
    /**
     * @brief Set the positions of the particles
     */
    void samplePositions( std::vector<Particle>& _particles );

    /**
     * @brief Set the momenta of the particles
     */
    void sampleMomenta( std::vector<Particle>& _particles );
    
    /** 
     * @brief Sampling of time and positions of one parton 
     **/
    void sample_TXYZ_one_partcl( Particle& _particle, bool& soft );
    
    int numberOfTestparticles;
    int numberOfParticlesToGenerate;
    
    /** @brief Particles with transverse momenta below this value are not considered for simulation. */
    double minimumPT;
    
    /** @brief Number of independent heavy ion events. Particles from different events might not allowed to scatter with each other. */
    int nEventsAA;
 
    
  private:
    std::string filename_pythiaParticleFile;
};



/** @brief exception class for handling unexpected critical behaviour within generation of mini-jet initial distributions  */
class ePythia_error : public std::runtime_error
{
  public:
    explicit ePythia_error(const std::string& what) : std::runtime_error(what) {};
    
    virtual ~ePythia_error() throw() {};
};




#endif // INITIALMODEL_PYTHIA_H
