//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------
//---------------------------------------------

#ifndef MINIJETS_H
#define MINIJETS_H

#include <stdexcept>
#include <vector>
#include <string>

#include "configuration.h"
#include "particle.h"
#include "vegas.h"
#include "woodsaxon.h"
#include "initialstatemodel.h"
#include "integrand_time.h"



class miniJets : public initialStateModel
{
  public:
    miniJets( const config& _config, WoodSaxon& _WoodSaxonParameter, STORED_TABLE_USAGE _storedTableUsage = computeNewTables );
    ~miniJets() {};
    
    void populateParticleVector( std::vector<Particle>& _particles );
    void populateParticleVector( std::vector<Particle>& _particles, const int _numberOfParticlesToGenerate, const double _minimumPT );
    
    /** @brief nuclear density (Woods-Saxon distribution) static version */
    static double densityA(double b, double z, const WoodSaxon& _w);   
        
  private:
    void samplePositions( std::vector<Particle>& _particles );
    void sampleMomenta( std::vector<Particle>& _particles );
    /** @brief PT-sampling according to the calculated d(sigma)/d(PT) => PX, PY */
    void sample_PXY( std::vector<Particle>& _particles ) const;
    /** @brief Sampling of Y1 and Y2 at given PT => PZ, E */
    void sample_PZE( std::vector<Particle>& _particles ) const;
    /** @brief sampling of the flavors of the parton pair at given PT,Y1 and Y2 */
    void sample_FLAV( std::vector<Particle>& _particles ) const;
    
    /** @brief sampling of the collision times */
    void sample_T( std::vector<Particle>& _particles ) const;
    /** @brief sampling of the positions of the parton pair at the given time */
    void sample_XYZ( std::vector<Particle>& _particles ) const;
    
    void computeWoodSaxonParameters( const config& _config, WoodSaxon& _WoodSaxonParameter );
    void generateSamplingDataSets();
    bool samplingDataSetsExist() const;
    
    double generateTimeDistribution( int& _count );
    double generatePtDistribution( int& count, int& count_fine );
    double max( double PT, double& y );
    
    /** @brief nuclear density (Woods-Saxon distribution) */
    double densityA(double b, double z) const;   
        
    /** @brief interpolation routine */
    void polint( const double xa[], const double ya[], const int n, const double x, double *y, double *dy ) const;
    
    void setDataFilesProperties( const config& _config, STORED_TABLE_USAGE _storedTableUsage );
    
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
    /** @brief lower PT-cutoff, GeV */
    double P0;
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
    int nEntries_PT;
    int nEntries_PT_fine;
    int numberOfParticlesToGenerate;
    
    std::string filename_samplingData_collisionTimes;
    std::string filename_samplingData_PT;
    std::string filename_samplingData_PT_fine;
    std::string filename_samplingData_maximumPT;
};





class integrand_distPT : public integrand
{
  public:
    integrand_distPT() {};
    integrand_distPT(const double pt_arg, const double sqrtS_arg) : pt(pt_arg), sqrtS(sqrtS_arg) {};
    ~integrand_distPT() {};
    
    /** @brief Overloaded operator() that makes integrand23 a function object - for use with CUBA-vegas */
    void operator()(const int*, const double [], const int*, double []) const;
    /** @brief Overloaded operator() that makes integrand23 a function object - for use with NR-vegas */
    double operator()(const double [], double) const;
    
    void setPt(const double pt_arg) {pt = pt_arg;}
    void setSqrtS(const double sqrtS_arg) {sqrtS = sqrtS_arg;}
    
  private:
    double pt;
    double sqrtS;
};





/** @brief exception class for handling unexpected critical behaviour within generation of mini-jet initial distributions  */
class eMiniJet_error : public std::runtime_error
{
  public:
    explicit eMiniJet_error(const std::string& what) : std::runtime_error(what) {};
    
    virtual ~eMiniJet_error() throw() {};
};



//---------------------------------------------//
// the cross section d(sigma)/d(PT2)d(y1)d(y2) //
//---------------------------------------------//
double csjet(double sqrtS, double PT, double Y1, double Y2);




#endif // MINIJETS_H
