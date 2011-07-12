//---------------------------------------------------------------------------------------
//provided by subversion
//---------------------------------------------------------------------------------------
//$HeadURL: file:///home/fochler/svnrepos/cascade_full/main/branches/restructure_branch/src/particle.cpp $
//$LastChangedDate: 2010-07-19 21:11:09 +0200 (Mon, 19 Jul 2010) $
//$LastChangedRevision: 144 $
//$LastChangedBy: fochler $
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------

/** @file
 * @brief definitions for the Particle class
 */

/** @author Oliver Fochler */

#include "configuration.h"
#include "particle.h"

/**
* Initialize static member for unique numbering of particles
*/
long int Particle::unique_id_counter = 0;
long int Particle::unique_id_counter_added = -1;


Particle::Particle() : T_creation( 0 ),
    X_init( 0 ), Y_init( 0 ), Z_init( 0 ), X_traveled( 0 ),
    PX_init( 0 ), PY_init( 0 ), PZ_init( 0 ), E_init( 0 ),
    X_lastInt( 0 ), Y_lastInt( 0 ), Z_lastInt( 0 ), T_lastInt( 0 ),
    T( 0 ), X( 0 ), Y( 0 ), Z( 0 ),
    E( 0 ), PX( 0 ), PY( 0 ), PZ( 0 ),
    eta( 0 ),
    m( 0 ),
    md2g( 0 ),
    md2q( 0 ),
    FLAVOR( gluon ),
    cell_id( -5 ),
    dead( false ),
    unique_id( -1 ),
    edge( false ),
    coll_id( -1 ),
    collisionTime( ns_casc::infinity ),
    collisionPartner( -1 ),
    isProjectile( false ),
    Eold( 0 ), PXold( 0 ), PYold( 0 ), PZold( 0 ),
//     as22( 0 ), as23( 0 ),
    rate( 0 ), ratev( 0 ),
//     rate23( 0 ), rate32( 0 ), rate22( 0 ),
//     rate23v( 0 ), rate32v( 0 ), rate22v( 0 ),
//     cs22( 0 ), cs23( 0 ), cs22t( 0 ), cs23t( 0 ),
//     lambda_scaled( 0 ),
//     md2g_scaled_22( 0 ),
//     md2q_scaled_22( 0 ),
//     md2g_scaled_23( 0 ),
    md2q_scaled_23( 0 ),
    free( false ),
    init( false )
//     step( -1 ), tstep( -1 ), taustep( -1 )
{

}




/**
 * Writes the 4-momentum vector to an array passed as argument arg. No checking as to correctness of array size is done.
 * Original values in arg are overwritten.
 *
 * @param[out] arg
 */
void Particle::getMomentumArray( double* const arg ) const
{
  arg[0] = E;
  arg[1] = PX;
  arg[2] = PY;
  arg[3] = PZ;
}


/**
 * Writes the space-time vector to an array passed as argument arg. No checking as to correctness of array size is done.
 * Original values in arg are overwritten.
 *
 * @param[out] arg
 */
void Particle::getCoordinateArray( double* const arg ) const
{
  arg[0] = T;
  arg[1] = X;
  arg[2] = Y;
  arg[3] = Z;
}

// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on; 
