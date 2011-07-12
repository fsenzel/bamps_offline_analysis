//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------
//---------------------------------------------

#ifndef PYTHIA_INITIAL_DISTRIBUTION_H
#define PYTHIA_INITIAL_DISTRIBUTION_H

#include <stdexcept>
#include <vector>
#include <string>

#include "initialstatemodel.h"
#include "configuration.h"
#include "particle.h"
#include "vegas.h"
#include "woodsaxon.h"
#include "integrand_time.h"



class pythiaInitialDistribution : public initialStateModel
{
  public:
    pythiaInitialDistribution( const config& _config, WoodSaxon& _WoodSaxonParameter, STORED_TABLE_USAGE _storedTableUsage = computeNewTables );
    ~pythiaInitialDistribution() {};
    
    void populateParticleVector( std::vector<Particle>& _particles );
    void populateParticleVector( std::vector<Particle>& _particles, const int _numberOfParticlesToGenerate, const double _minimumPT );

    
        
  private:
    void samplePositions( std::vector<Particle>& _particles );
    void sampleMomenta( std::vector<Particle>& _particles );
    
    void setDataFilesProperties( STORED_TABLE_USAGE _storedTableUsage );
    
    /** @brief sampling of the collision time of one particle */
    void sample_T_one_partcl( std::vector<Particle>& _particles, const int ) const;
    /** @brief sampling of the positions of one parton at the given time */
    void sample_XYZ_one_partcl( std::vector<Particle>& _particles, const int number, bool& soft ) const;
    
    void computeWoodSaxonParameters( const config& _config, WoodSaxon& _WoodSaxonParameter );
    void generateSamplingDataSets();
    bool samplingDataSetsExist() const;
    
    double generateTimeDistribution( int& _count );
    
    /** @brief nuclear density (Woods-Saxon distribution) */
    double densityA(double b, double z) const;   
        
    /** @brief interpolation routine */
    void polint( const double xa[], const double ya[], const int n, const double x, double *y, double *dy ) const;
    
    
    /** @brief mass number of nucleus A */
    double A;          
    /** @brief atomic number, i.e. number of protons, of nucleus A */
    double Aatomic;    
    /** @brief mass number of nucleus B */
    double B;
    /** @brief atomic number of nucleus B */
    double Batomic;
    double impactParameter;
    double sqrtS_perNN;
    int numberOfTestparticles;
    
    /** @brief minimum PT to be sampled (GeV) */
    double minimumPT;
    
    double gamma;
    double velocity;
    double RA0;
    double RA;
    double dA;
    double n0A;
    
    WoodSaxon WoodSaxonParameter;
    
    int nEntries_collisionTimes;
    int numberOfParticlesToGenerate;
    
    std::string filename_samplingData_collisionTimes;
    std::string filename_pythiaParticleFile;
};



/** @brief exception class for handling unexpected critical behaviour within generation of mini-jet initial distributions  */
class ePythia_error : public std::runtime_error
{
  public:
    explicit ePythia_error(const std::string& what) : std::runtime_error(what) {};
    
    virtual ~ePythia_error() throw() {};
};




#endif 
