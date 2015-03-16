//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------
//---------------------------------------------

#ifndef INITIALMODEL_PARTICLESFROMFILE_H
#define INITIALMODEL_PARTICLESFROMFILE_H

#include <stdexcept>
#include <vector>
#include <string>

#include "initialmodel.h"
#include "configuration.h"
#include "particle.h"
#include "woodsaxon.h"


/**
 * @brief Class to provide the initialization which reads in particles from particle data file (eg. from PYTHIA or MC@NLO)
 */
class initialModel_ParticlesFromFile : public initialModelWS
{
public:
    initialModel_ParticlesFromFile( const string _filename_particleFile, const config& _config, WoodSaxon& _WoodSaxonParameter, const double _minimumPT = 0.0, const int _nToGenerate = -1 );
    ~initialModel_ParticlesFromFile() {};

    void populateParticleVector( std::vector<Particle>& _particles );


protected:
    /**
     * @brief Set the positions of the particles
     */
    void samplePositions( std::vector<Particle>& _particles );

    /**
     * @brief Set the momenta of the particles
     */
    virtual void sampleMomenta( std::vector<Particle>& _particles ) = 0;

    /**
     * @brief Sampling of time and positions of one parton
     **/
    void sample_TXYZ_one_partcl( Particle& _particle, bool& soft );

    void changeCharmMass( std::vector< Particle >& _particles, const double M_old );
    void changeBottomMass( std::vector< Particle >& _particles, const double M_old );

    int numberOfTestparticles;
    int numberOfParticlesToGenerate;

    /** @brief Particles with transverse momenta below this value are not considered for simulation. */
    double minimumPT;

    /** @brief Number of independent heavy ion events. Particles from different events might not allowed to scatter with each other. */
    int nEventsAA;

    std::string filename_particleFile;
};



/** @brief exception class for handling unexpected critical behaviour within generation of mini-jet initial distributions  */
class eParticlesFromFile_error : public std::runtime_error
{
public:
    explicit eParticlesFromFile_error(const std::string& what) : std::runtime_error(what) {};

    virtual ~eParticlesFromFile_error() throw() {};
};




#endif // INITIALMODEL_PARTICLESFROMFILE_H
