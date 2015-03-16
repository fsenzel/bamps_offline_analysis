//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------
//---------------------------------------------

#ifndef INITIALMODEL_JPSI_H
#define INITIALMODEL_JPSI_H

#include "configuration.h"
#include "initialmodel.h"
#include "interpolation_iniJpsi.h"



class initialModel_Jpsi : public initialModelWS
{
public:
    initialModel_Jpsi( const config& _config, WoodSaxon& _WoodSaxonParameter );

    void populateParticleVector( std::vector<Particle>& _particles );

private:
    void sample_metropolis_dndptdy(double& pt_arg, double& y_arg);
    void sample_PXYZE_FLAV_singleParticle( Particle& _tempParticle, const int ParticleNumber );

    interpolation_iniJpsi_dndptdy theInterpolation_dndptdy;

    double sigmaAbs, agN;
    shadowModelJpsi shadowing_model;

    int testparticles;

    /** @brief Number of independent heavy ion events. Particles from different events might not allowed to scatter with each other. */
    int nEventsAA;

};


/**
 * @brief exception class for handling unexpected critical behaviour
 * within generation of mini-jet initial distributions
 */
class eJpsi_error : public std::runtime_error
{
public:
    explicit eJpsi_error(const std::string& what) : std::runtime_error(what) {};

    virtual ~eJpsi_error() throw() {};
};

#endif // INITIALMODEL_JPSI_H
