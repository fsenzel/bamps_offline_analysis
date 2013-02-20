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

#include <vector>
#include <string>

#include "initialmodel_particlesFromFile.h"
#include "configuration.h"
#include "particle.h"
#include "woodsaxon.h"


/**
 * @brief Class to provide the Pythia initialization
 */
class initialModel_Pythia : public initialModel_ParticlesFromFile
{
  public:
    initialModel_Pythia( const config& _config, WoodSaxon& _WoodSaxonParameter, const double _minimumPT = 0.0, const int _nToGenerate = -1 ) :
    initialModel_ParticlesFromFile( _config.getPythiaParticleFile(), _config, _WoodSaxonParameter, _minimumPT, _nToGenerate) {};

  private:
    /**
     * @brief Set the momenta of the particles
     */
    void sampleMomenta( std::vector<Particle>& _particles );

};


#endif // INITIALMODEL_PYTHIA_H
