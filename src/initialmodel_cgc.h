//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


// at revision 1004, this is identical to
// full/branches/vector4D/src/initialmodel_cgc.h


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
