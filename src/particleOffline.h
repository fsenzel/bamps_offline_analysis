//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file
 * @brief Declarations for the Offline Particle class
 */

#ifndef PARTICLEOFFLINE_H
#define PARTICLEOFFLINE_H

#include <vector>

#include "particle.h"


/**
 * @brief Provides properties of a particle needed in the offline reconstruction of BAMPS events.
 *
 * This class extends the Particle class that is used in the standard BAMPS simulations to include some
 * variables that are needed for the offline reconstruction.
 */
class ParticleOffline : public Particle
{
public:
  /** @brief Provide standard constructor (for completeness) */
  ParticleOffline() : 
    Particle(), 
    temperature(0), 
    initially_produced( true ), 
    jpsi_dissociation_number( -1 ),
    lambda_added( -1 ), lambda_added_old( -1 ), rate_added( -1 ),
    T_creation( 0 ), 
    X_traveled( 0 ),
    PosInit( ),
    MomInit( ),
    lastInt( ),    
    rate( 0 ), ratev( 0 ),
    isAlreadyInAddedParticles( 0 ),
    rate_added_32( 0.0 ),
    collision_tag(false),
    production_time(0.0)
  {};
    
  ParticleOffline( const Particle& _particle ) : 
    Particle( _particle ), 
    temperature(0), 
    initially_produced( true ), 
    jpsi_dissociation_number( -1 ),
    lambda_added( -1 ), lambda_added_old( -1 ), rate_added( -1 ),
    T_creation( 0 ), 
    X_traveled( 0 ),
    PosInit( ),
    MomInit( ),
    lastInt( ),    
    rate( 0 ), ratev( 0 ),
    isAlreadyInAddedParticles( 0 ),
    rate_added_32( 0.0 ),
    collision_tag(false),
    production_time(0.0)
  {};
    
  /** @brief counter for unique particle IDs of added particles (static) */
  static long int unique_id_counter_added;

  /** @brief Temperature of the surrounding medium, needed for J/psi melting */
  double temperature;
    
  /** @brief Whether the particle was initially produced or later in a secondary process */
  bool initially_produced;
    
  /** @brief Flag for discerning (background) particles that have already produced a photon */
  bool collision_tag;  
  
  /** @brief Producing time of the photon. This may prove useful in discriminating the building up of v2. */
  double production_time;   
  
  /** @brief Unique number of jpsi dissociation such that the same c+cbar do not reunite directly */
  int jpsi_dissociation_number;
    
  /** @brief If c+cbar form a Jpsi the variable N_EVENT_pp of the cbar is stored in this variable to be still accessible, In particular if the Jpsi dissociates again. */
  int N_EVENT_Cbar;

  /** @brief Mean free path of the added particle in this time step */
  double lambda_added; // fm
    
  /** @brief Mean free path of the added particle in the previous time step */
  double lambda_added_old; // fm
    
  /** @brief Rate of the added particle in this time step */
  double rate_added; // 1/fm

  /** @brief 3->2 rate of the added particle in this time step */
  double rate_added_32; // 1/fm

  /** stuff special to offline reconstruction */
  double T_creation;
  double X_traveled;//fm
  VectorTXYZ PosInit;
  VectorEPxPyPz MomInit;
  VectorTXYZ lastInt;
    
  double rate, ratev;

  /** @brief vector which holds information in which event this medium particle is already in added particles list */
  std::vector< bool > isAlreadyInAddedParticles;

    
  static int mapToPDGCodes( const FLAVOR_TYPE _flav )
  {
    switch ( _flav )
    {
      case gluon:
        return 21;
        break;
      case up:
        return 2;
        break;
      case anti_up:
        return -2;
        break;
      case down:
        return 1;
        break;
      case anti_down:
        return -1;
        break;
      case strange:
        return 3;
        break;
      case anti_strange:
        return -3;
        break;
      case charm:
        return 4;
        break;
      case anti_charm:
        return -4;
        break;
      default:
        return 0;
        break;
    }      
  }
    
private:
};

/**
 * @brief Provides basic properties of a particle, used for electrons from heavy flavor decays.
 * @author Oliver Fochler
 *
 * Use only basic properties of the ParticlePrototype class to minimize memory allocation. This is possible because heavy flavor electrons are not propagated through the medium. 
 */
// class ParticleHFelectron : public ParticlePrototype
class ParticleHFelectron : public ParticleOffline
{
};


#endif
