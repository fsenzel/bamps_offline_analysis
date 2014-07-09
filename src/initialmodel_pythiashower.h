//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL: svn+ssh://senzel@th.physik.uni-frankfurt.de/home/bamps/svn/full/offlineAnalysis/trunk/src/initialmodel_pythiashower.h $
//$LastChangedDate: $
//$LastChangedRevision: -1 $
//$LastChangedBy: $
//---------------------------------------------
//---------------------------------------------

#ifndef INITIALMODEL_PYTHIASHOWER_H
#define INITIALMODEL_PYTHIASHOWER_H

#include <string>

#include "initialmodel.h"
#include "configuration.h"
#include "particle.h"

enum SHOWER_TYPE { fixedShower, fixedParton, pythiaShower };

class initialModel_PYTHIAShower : public initialModelWS
{
  public:
    initialModel_PYTHIAShower( const config& _config, WoodSaxon& _WoodSaxonParameter, const double _minimumPT );
    initialModel_PYTHIAShower( const config& _config, WoodSaxon& _WoodSaxonParameter, const double _initialPartonPt, const int _initialPartonFlavor, const SHOWER_TYPE _shower_type );
    ~initialModel_PYTHIAShower() {};
    
    void populateParticleVector( std::vector<Particle>& _particles );
        
  private:

    /** @brief Random number generator seed for fixing PYTHIA seed */
    uint32_t seed;
  
    /** @brief Filename prefix, needed for initial unshowered particle output */
    string filename_prefix; 
  
    /** @brief Routine for creating one dijet shower out of PYSHOW */
    vector<Particle> getFixedShowerEvent( const double _px, const double _py, const double _pz, const FLAVOR_TYPE _flavorA, const FLAVOR_TYPE _flavorB );

    /** @brief Routine for creating one event out of PYTHIA */
    void getPythiaShowerEvent( vector<Particle> &_particles, vector<Particle> &_initialPartons );
    
    /** @brief lower PT-cutoff of PYTHIA spectrum (in GeV) */
    double P0;
  
    /** @brief Number of heavy ion collision events, set on top of the offline reconstruction. 
    * The number of particles in one such heavy ion collision event is equal to < number of produced particles 
    * in pp > * Ntest * Nbin. Consequently, 1 would mean that one adds as many particles as there are in a heavy 
    * ion collision times Ntest, making it analogously to a standard BAMPS simulation.
    * 
    * It is only used if the number of particles to sample is not explicitly given (in which case it is negative)
    */
    int nEventsToGenerate;
  
    /** @brief Type of shower: shower with fixed initial partons, fixed single parton, or full pythia di-jet event. */
    SHOWER_TYPE shower_type;
    
    /** @brief Transverse momentum of initial parton pair */
    double initialPartonPt;
    
    /** @brief Flavor of initial parton pair */
    int initialPartonFlavor;
    
    FLAVOR_TYPE mapToPYTHIAflavor ( int _pythiaFlavor )
    {
      switch( _pythiaFlavor )
      {
        case 21:
          return gluon;
        case 2:
          return up;
        case 1:
          return down;
        case -2:
          return anti_up;
        case -1:
          return anti_down;
        case 3:
          return strange;
        case -3:
          return anti_strange;
        default:
//           std::cout << "Unknown pythia flavor:\t" << _pythiaFlavor << std::endl;
          return allFlavors;
      }
  };

};


/** 
 * @brief exception class for handling unexpected critical behaviour
 * within generation of mini-jet initial distributions 
 */
class ePYTHIAShower_error : public std::runtime_error
{
  public:
    explicit ePYTHIAShower_error(const std::string& what) : std::runtime_error(what) {};
    
    virtual ~ePYTHIAShower_error() throw() {};
};


#endif // INITIALMODEL_PYTHIASHOWER_H
