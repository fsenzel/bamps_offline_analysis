//---------------------------------------------------------------------------------------
//provided by subversion
//---------------------------------------------------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------

/** @file
* @brief Declarations for the Particle class
*/

#ifndef PARTICLE_H
#define PARTICLE_H

#include "particleprototype.h"
#include <vector>



/**
* @brief Provides properties of a particle needed in the simulations.
*
* This class encapsulates properties of particles that are needed for the simulation process, such as position and momentum variables.
* It is derived from ParticlePrototype (and extends it) that provides general properties of particles.
*/
class Particle : public ParticlePrototype
{
  public:
    /** @brief Provide standard constructor (for completeness) */
    Particle() : ParticlePrototype(), eta( 0 ), md2g( 0 ), md2q( 0 ), N_EVENT_pp( 0 ), HARD( true ), 
    N_EVENT_AA( 0 ), edge( -1 ), coll_id( -1 ), collisionTime( 0 ), collisionPartner( -1 ), PXold( 0 ), 
    PYold( 0 ), PZold( 0 ), as22( 0 ), as23( 0 ), rate23v( 0 ), rate32v( 0 ), rate22v( 0), cs22( 0 ), cs23( 0 ),
    lambda_scaled( 0 ), free( true ), init( true ), 
//     md2g_scaled_22( 0 ),md2q_scaled_22( 0 ), md2g_scaled_23( 0 ), md2q_scaled_23( 0 ), 
    step( 0 ), tstep( 0 ), taustep( 0 ) {};
    
    /** @brief space time rapidity \eta */
    double eta;
    
    /** @brief gluonic screening mass currently associated with the particle object, scaled by alpha_s as md2g / alpha_s [GeV^2] */
    double md2g;
    /** @brief quark screening mass currently associated with the particle object, scaled by alpha_s as md2g / alpha_s [GeV^2] */
    double md2q;
    
    /** @brief Pythia event number */
    int N_EVENT_pp;
    /** @brief Pythia hard or soft scattering */
    bool HARD; // true/1 if parton comes from hard scattering
    
    /** @brief Event number of heavy ion collision to which particle belongs. Necessary if the number of added particles is much larger than the number of particles which would be present in a event according to the test particles number of offline particles. This is only important if one considers scatterings among the added particles.  */
    int N_EVENT_AA;

    /** @brief index of edge cell the particle belongs to, edge = -1 corresponds to no edge cell */
    short int edge;
    
    /** @brief Unique ID of the last collision of the particle */
    long coll_id;
    
    /** @brief Flag for free streaming, true for particles in regions with energy density < config::freezeOutEnergyDensity or for particles outside the grid */
    bool free;
    /** @brief Flag for discerning particles that are still within their initial formation time */
    bool init;

    /** @brief collision ordering time [fm] (geometric collisions) */
    double collisionTime;
    /** @brief collision partner (geometric collisions) */
    int collisionPartner;
    
    /** @brief Momentum prior to geometric collision, needed for particles in "edge cell"*/
    double PXold;
    /** @brief Momentum prior to geometric collision, needed for particles in "edge cell"*/
    double PYold;
    /** @brief Momentum prior to geometric collision, needed for particles in "edge cell"*/
    double PZold;
    
    /** @brief Mean alpha_s for 2->2 interactions associated with this particle, averaged over cell in previous time step */ 
    double as22;
    /** @brief Mean alpha_s for 2->3 interactions associated with this particle, averaged over cell in previous time step */
    double as23;
    /** @brief Previous mean 2->3 rate (in GeV) associated with this particle, averaged over cell in the second to last time step */ 
    double rate23v;
    /** @brief Previous mean 3->2 rate (in GeV) associated with this particle, averaged over cell in the second to last time step */ 
    double rate32v;
    /** @brief Previous mean 2->2 rate (in GeV) associated with this particle, averaged over cell in the second to last time step */ 
    double rate22v;
    /** @brief Mean 2->2 cross section (1/GeV^2) associated with this particle, averaged over cell in previous time step */
    double cs22;
    /** @brief Mean 2->3 cross section (1/GeV^2) associated with this particle, averaged over cell in previous time step */
    double cs23;
    /** @brief Mean lambda_scaled associated with this particle, averaged over cell in previous time step */
    double lambda_scaled;
//     /** @brief Mean md2g (scaled with s) from 2->2 interactions, averaged over cell in previous time step */
//     double md2g_wo_as_scaled_22;
//     /** @brief Mean md2q (scaled with s) from 2->2 interactions, averaged over cell in previous time step */
//     double md2q_wo_as_scaled_22;
//     /** @brief Mean md2g (scaled with s) from 2->3 interactions, averaged over cell in previous time step */
//     double md2g_wo_as_scaled_23;
//     /** @brief Mean md2q (scaled with s) from 2->3 interactions, averaged over cell in previous time step */
//     double md2q_wo_as_scaled_23;
    
    int step,tstep,taustep;//fm
    
  private:
};



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
    ParticleOffline() : Particle(), T_creation( 0 ), X_init( 0 ), Y_init( 0 ), Z_init( 0 ), X_traveled( 0 ),
    PX_init( 0 ), PY_init( 0 ), PZ_init( 0 ), E_init( 0 ),
    X_lastInt( 0 ), Y_lastInt( 0 ), Z_lastInt( 0 ), T_lastInt( 0 ),
    lambda_added( -1 ), lambda_added_old( -1 ), rate_added( -1 ),
    Eold( 0 ), rate( 0 ), ratev( 0 ), temperature(0), initially_produced( true ), jpsi_dissociation_number( -1 ),
    isAlreadyInAddedParticles( 0 ) {};
    
    ParticleOffline( const Particle& _particle ) : Particle( _particle ), T_creation( 0 ), X_init( 0 ), Y_init( 0 ), Z_init( 0 ), X_traveled( 0 ),
    PX_init( 0 ), PY_init( 0 ), PZ_init( 0 ), E_init( 0 ),
    X_lastInt( 0 ), Y_lastInt( 0 ), Z_lastInt( 0 ), T_lastInt( 0 ),
    lambda_added( -1 ), lambda_added_old( -1 ), rate_added( -1 ),
    Eold( 0 ), rate( 0 ), ratev( 0 ), temperature(0), initially_produced( true ), jpsi_dissociation_number( -1 ),
    isAlreadyInAddedParticles( 0 ) {};
    
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
    
    /** @brief Mean free path of the added particle in this time step */
    double lambda_added; // fm
    
    /** @brief Mean free path of the added particle in the previous time step */
    double lambda_added_old; // fm
    
    /** @brief Rate of the added particle in this time step */
    double rate_added; // 1/fm
    
    /** stuff special to offline reconstruction */
    double T_creation;
    double X_init, Y_init,Z_init, X_traveled;//fm
    double PX_init, PY_init,PZ_init, E_init;//GeV
    double X_lastInt, Y_lastInt, Z_lastInt, T_lastInt;
    
    double Eold;
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
    
    /** @brief vector which holds information in which event this medium particle is already in added particles list */
    std::vector< bool > isAlreadyInAddedParticles;
    
  private:
};

/**
* @brief Provides basic properties of a particle, used for electrons from heavy flavor decays.
*
* Use only basic properties of the ParticlePrototype class to minimize memory allocation. This is possible because heavy flavor electrons are not propagated through the medium. 
*/
// class ParticleHFelectron : public ParticlePrototype
class ParticleHFelectron : public ParticleOffline
{
};

#endif
