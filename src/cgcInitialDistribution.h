//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------
//---------------------------------------------

#ifndef CGC_INITIAL_DISTRIBUTION_H
#define CGC_INITIAL_DISTRIBUTION_H

#include <string>

#include "initialstatemodel.h"
#include "configuration.h"
#include "particle.h"


class cgcInitialDistribution : public initialStateModel
{
  public:
    cgcInitialDistribution( const config& _config);
    ~cgcInitialDistribution() {};
    
    void populateParticleVector( std::vector<ParticleOffline>& _particles );
  
        
  private:
    std::string filename_cgcParticleFile;
    int numberOfParticlesToGenerate;
};


#endif 
