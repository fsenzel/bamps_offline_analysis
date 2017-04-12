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

enum SHOWER_TYPE { fixedShower, fixedParton, inclusiveShower = 2, photonShower = 3, charmQuarkShower = 4, bottomQuarkShower = 5 };

class initialModel_PYTHIAShower : public initialModelWS
{
  public:
    initialModel_PYTHIAShower( const config& _config, WoodSaxon& _WoodSaxonParameter, const SHOWER_TYPE _shower_type, const double _ptCutOff );
    initialModel_PYTHIAShower( WoodSaxon& _WoodSaxonParameter, const SHOWER_TYPE _shower_type,
                               const int _Nevents, const uint32_t _seed, const string _filename_prefix,
                               const double _A, const double _Aatomic, const double _B, const double _Batomic,
                               const double _sqrtS_perNN, const double _impactParameter,
                               const double _ptCutOff, const FLAVOR_TYPE _initialPartonFlavor = light_parton );

    ~initialModel_PYTHIAShower() {};
    
    void populateParticleVector( std::vector<Particle>& _particles );
    void populateParticleVector( std::vector<Particle>& _particles, std::vector<Particle>& _initial_partons );

    /** @brief Routine for creating one parton shower event out of PYTHIA */
    void getShowerEvent( vector<Particle> &_particles, vector<Particle> &_initial_parton_pair );

private:

    /** @brief Random number generator seed for fixing PYTHIA seed */
    uint32_t seed;
  
    /** @brief Filename prefix, needed for initial unshowered particle output */
    string filename_prefix; 
  
    /** @brief lower PT-cutoff of PYTHIA spectrum (in GeV) */
    double ptCutOff;
  
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

    /** @brief Flavor of initial parton pair */
    FLAVOR_TYPE initial_parton_flavor;
    
    FLAVOR_TYPE mapToPYTHIAflavor( int _pythiaFlavor )
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
        case 4:
          return charm;
        case 5:
          return bottom;
        case -4:
          return anti_charm;
        case -5:
          return anti_bottom;
        case 22:
          return photon;
        default:
//           std::cout << "Unknown pythia flavor:\t" << _pythiaFlavor << std::endl;
          return allFlavors;
      }
    };

    int mapToBAMPSflavor( FLAVOR_TYPE _bampsFlavor )
    {
      switch( _bampsFlavor )
      {
        case gluon:
          return 21;
        case up:
          return 2;
        case down:
          return 1;
        case strange:
          return 3;
        case charm:
          return 4;
        case bottom:
          return 5;
        case anti_up:
          return -2;
        case anti_down:
          return -1;
        case anti_strange:
          return -3;
        case anti_charm:
          return -4;
        case anti_bottom:
          return -5;
        case photon:
          return 22;
        case light_parton:
          return 101;
        case allFlavors:
          return 111;
        default:
          string errMsg = "Unknown Pythia flavor for converting to BAMPS. Unrecoverable error!";
          throw eInitialModel_error( errMsg );
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
