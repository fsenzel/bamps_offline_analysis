//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------
//---------------------------------------------

#ifndef HEAVY_QUARK_INITIAL_DISTRIBUTION_H
#define HEAVY_QUARK_INITIAL_DISTRIBUTION_H

#include <vector>
#include <string>

#include "initialmodel.h"
#include "configuration.h"
#include "particle.h"
#include "woodsaxon.h"

enum SEPERATE_HEAVY_QUARK_INITIAL_MODEL { none, pythia, mcAtNLO };


class initialModel_heavyQuarks : public initialModelWS
{
  public:
    initialModel_heavyQuarks ( const config& _config, WoodSaxon& _WoodSaxonParameter, SEPERATE_HEAVY_QUARK_INITIAL_MODEL _seperateHeavyQuarkInitialModel = none );
    ~initialModel_heavyQuarks() {};
 
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
    
  private:
    SEPERATE_HEAVY_QUARK_INITIAL_MODEL seperateHeavyQuarkInitialModel;
    std::string filename_heavyQuarkParticleFile;
    
    double numberOfTestparticles_heavyQuarkParticleFile;
    int numberOfTestparticles_BAMPS;
    
    int numberOfParticlesToGenerate;
};

#endif 
