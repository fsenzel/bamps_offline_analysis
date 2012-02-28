//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------
//---------------------------------------------

#ifndef INITIALMODEL_CGC_H
#define INITIALMODEL_CGC_H

#include <string>

#include "initialmodel.h"
#include "configuration.h"
#include "particle.h"


class initialModel_CGC : public initialModel
{
  public:
    initialModel_CGC( const config& _config);
    ~initialModel_CGC() {};
    
    void populateParticleVector( std::vector<Particle>& _particles );
  
        
  private:
    std::string filename_cgcParticleFile;
    int numberOfParticlesToGenerate;
};


#endif // INITIALMODEL_CGC_H
