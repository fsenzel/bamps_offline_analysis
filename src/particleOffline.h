//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: svn+ssh://gallmei@th.physik.uni-frankfurt.de/home/bamps/svn/full/branches/vector4D/src/particle.h $
//$LastChangedDate: 2012-07-09 14:50:33 +0200 (Mo, 09. Jul 2012) $
//$LastChangedRevision: 714 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file
* @brief Declarations for the Offline Particle class
*/

#ifndef PARTICLEOFFLINE_H
#define PARTICLEOFFLINE_H

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
    ParticleOffline() : Particle(), T_creation( 0 ), 
			PosInit( 0,0,0,0 ), X_traveled( 0 ),
			MomInit( 0,0,0,0 ), 
			lastInt( 0,0,0,0 ),
    rate( 0 ), ratev( 0 ), temperature(0), initially_produced( true ), jpsi_dissociation_number( -1 ) {};
    
    ParticleOffline( const Particle& _particle ) : 
      Particle( _particle ), T_creation( 0 ), 
      PosInit( 0,0,0,0 ), X_traveled( 0 ),
      MomInit( 0,0,0,0 ),
      lastInt( 0,0,0,0 ),	   
      rate( 0 ), ratev( 0 ), temperature(0), initially_produced( true ), jpsi_dissociation_number( -1 ) {};
    
    /** @brief counter for unique particle IDs of added particles (static) */
    static long int unique_id_counter_added;

    /** @brief Temperature of the surrounding medium, needed for J/psi melting */
    double temperature;
    
    /** @brief Whether the particle was initially produced or later in a secondary process */
    bool initially_produced;
    
    /** @brief Unique number of jpsi dissociation such that the same c+cbar do not reunite directly */
    int jpsi_dissociation_number;
    
    /** @brief If c+cbar form a Jpsi the variable N_EVENT_pp of the cbar is stored in this variable to be still accessible, In particular if the Jpsi dissociates again. */
    int N_EVENT_Cbar;
    
    /** stuff special to offline reconstruction */
    double T_creation;
  double X_traveled;//fm
  VectorTXYZ PosInit;
  VectorEPxPyPz MomInit;
  VectorTXYZ lastInt;
    
    double rate, ratev;
    
    
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
