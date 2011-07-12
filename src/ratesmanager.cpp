//---------------------------------------------------------------------------------------
//provided by subversion
//---------------------------------------------------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------


#include "ratesmanager.h"

/**
 * Interface routine to store (collect) a probability for a given 2->X interaction process
 *
 * @param[in] _type Interaction process type id (221, 3280, etc., see class interactionType)
 * @param[in] _P Probability to be stored
 * @return 1 for successful execution
 */
int ratesManager::add( const GENERIC_COLL_TYPE genericType, const FLAVOR_TYPE _F1, const FLAVOR_TYPE _F2, const double _P )
{
  if ( isNormalized )
  {
    std::string errMsg = "ratesManager::add called after having normalized the rates via ratesManager::normalizeRates. Unrecoverable error.";
    throw eRatesManager_error( errMsg );
  }

  int index = interactionType::getIndexFromProcessType( interactionType::getInclusiveProcessType( _F1, _F2, genericType ), genericType );


  if ( genericType == c22 )
  {
    rates22[index] += _P;
  }
  else if ( genericType == c23 )
  {
    rates23[index] += _P;
  }
  else if ( genericType == c32 )
  {
    rates32[index] += _P;
  }
  else
  {
    std::string errMsg = "Generic interaction type not found in ratesManager::add. Unrecoverable error.";
    throw eRatesManager_error( errMsg );
  }

  return 1;
}




/**
* Interface routine to store (collect) a probability for a given 3->X interaction process
*
* @param[in] _type Interaction process type id (221, 3280, etc., see class interactionType)
* @param[in] _P Probability to be stored
* @return 1 for successful execution
*/
int ratesManager::add( const GENERIC_COLL_TYPE genericType, const FLAVOR_TYPE _F1, const FLAVOR_TYPE _F2, const FLAVOR_TYPE _F3, const double _P )
{
  if ( isNormalized )
  {
    std::string errMsg = "ratesManager::add called after having normalized the rates via ratesManager::normalizeRates. Unrecoverable error.";
    throw eRatesManager_error( errMsg );
  }

  int index = interactionType::getIndexFromProcessType( interactionType::getInclusiveProcessType( _F1, _F2, _F3, genericType ), genericType );

  if ( genericType == c32 )
  {
    rates32[index] += _P;
  }
  else
  {
    std::string errMsg = "Generic interaction type not found in ratesManager::add (3->2). Unrecoverable error.";
    throw eRatesManager_error( errMsg );
  }

  return 1;
}






/**
 * Normalize the accumulated probabilties such that one gets the interaction rate per single particle and per time step dt.
 * For this the accumulated probabilities for a given process is weighted with the number of gluons (quarks, anti-quarks)
 * involved in the initial state of this interaction.
 * e.g. the rate for a single gluons in ggg -> gg processes would be: 3 * sum(probabilities) / (Ng * dt)
 * See notes dating 27.05.2010 for more details.
 *
 * @param[in] Ng Number of gluons for which the interaction probabilities have been accumulated
 * @param[in] Nq Number of quarks for which the interaction probabilities have been accumulated
 * @param[in] Nqbar Number of anti-quarks for which the interaction probabilities have been accumulated
 * @param[in] dt Time step dt in fm/c
 * @return 1 for successful execution
 */
int ratesManager::normalizeRates( const int Ng, const int Nq, const int Nqbar, const double dt )
{
  if ( isNormalized )
  {
    std::string errMsg = "ratesManager::normalizeRates called twice. Unrecoverable error.";
    throw eRatesManager_error( errMsg );
  }
  else
  {

    int type = 0;

    ratesGluons.assign( 3, 0 );
    ratesQuarks.assign( 3, 0 );
    ratesAntiQuarks.assign( 3, 0 );

    for ( int i = 0; i < rates22.size(); i++ )
    {
      type = interactionType::getInclusiveProcessTypeFromIndex( i, c22 );
      ratesGluons[0] += interactionType::getInvolvedInitialInclusive( gluon, type ) * rates22[i];
      ratesQuarks[0] += interactionType::getInvolvedInitialInclusive( quark, type ) * rates22[i];
      ratesAntiQuarks[0] += interactionType::getInvolvedInitialInclusive( anti_quark, type ) * rates22[i];
      rates22[i] /= dt;  // this yields the rate for the given process - not per particle as for the entries in ratesGluons etc.
    }

    for ( int i = 0; i < rates23.size(); i++ )
    {
      type = interactionType::getInclusiveProcessTypeFromIndex( i, c23 );
      ratesGluons[1] += interactionType::getInvolvedInitialInclusive( gluon, type ) * rates23[i];
      ratesQuarks[1] += interactionType::getInvolvedInitialInclusive( quark, type ) * rates23[i];
      ratesAntiQuarks[1] += interactionType::getInvolvedInitialInclusive( anti_quark, type ) * rates23[i];
      rates23[i] /= dt;  // this yields the rate for the given process - not per particle as for the entries in ratesGluons etc.
    }

    for ( int i = 0; i < rates32.size(); i++ )
    {
      type = interactionType::getInclusiveProcessTypeFromIndex( i, c32 );
      ratesGluons[2] += interactionType::getInvolvedInitialInclusive( gluon, type ) * rates32[i];
      ratesQuarks[2] += interactionType::getInvolvedInitialInclusive( quark, type ) * rates32[i];
      ratesAntiQuarks[2] += interactionType::getInvolvedInitialInclusive( anti_quark, type ) * rates32[i];
      rates32[i] /= dt;  // this yields the rate for the given process - not per particle as for the entries in ratesGluons etc.
    }

    for ( int i = 0 ; i < 3; i++ )
    {
      ratesGluons[i] /= ( dt * Ng );
      ratesQuarks[i] /= ( dt * Nq );
      ratesAntiQuarks[i] /= ( dt * Nqbar );
    }

    ratesGluonsTotal = ratesGluons[0] + ratesGluons[1] + ratesGluons[2];
    ratesQuarksTotal = ratesQuarks[0] + ratesQuarks[1] + ratesQuarks[2];
    ratesAntiQuarksTotal = ratesAntiQuarks[0] + ratesAntiQuarks[1] + ratesAntiQuarks[2];

    isNormalized = true;

    return 1;
  }
}



int ratesManager::normalizeRates()
{
  if ( isNormalized )
  {
    std::string errMsg = "ratesManager::normalizeRates called twice. Unrecoverable error.";
    throw eRatesManager_error( errMsg );
  }
  else
  {
    ratesGluons.assign( 3, 0 );
    ratesQuarks.assign( 3, 0 );
    ratesAntiQuarks.assign( 3, 0 );

    ratesGluonsTotal = 0;
    ratesQuarksTotal = 0;
    ratesAntiQuarksTotal = 0;

    isNormalized = true;

    return 1;
  }
}



/**
 * Get the rate for a certain particle type
 *
 * @param[in] particleType Particle type for which the rate is requested (gluon, quark, anti-quark)
 * @param[in] unit (Optional, default = fm) Whether the result should be returned in 1/fm or in GeV
 * @return rate
 */
double ratesManager::getRate( const FLAVOR_TYPE particleType, const UNIT_TYPE unit ) const
{
  if ( !isNormalized )
  {
    std::string errMsg = "ratesManager::getRate called without prior normalization of rates. Unrecoverable error.";
    throw eRatesManager_error( errMsg );
  }
  else
  {
    double conversionFactor = 1;
    if ( unit == GeV )
    {
      conversionFactor = 0.197;
    }

    FLAVOR_TYPE pType = particleType;
    if ( pType == up || pType == down || pType == strange )
    {
      pType = quark;
    }
    else if ( pType == anti_up || pType == anti_down || pType == anti_strange )
    {
      pType = anti_quark;
    }

    switch ( pType )
    {
    case gluon:
      return ratesGluonsTotal * conversionFactor;
      break;
    case quark :
      return ratesQuarksTotal * conversionFactor;
      break;
    case anti_quark:
      return ratesAntiQuarksTotal * conversionFactor;
      break;
    default:
      std::string errMsg = "Particle type not found in ratesManager::getRate. Unrecoverable error.";
      throw eRatesManager_error( errMsg );
      break;
    }
  }
}


/**
 * Get the rate for a certain particle type and a certain generic collision type (22, 23, 32)
 * e.g. to get the rate of quarks in 2->2 collisions, call ratesManager::getRate( quark, c22, fm )
 *
 * @param[in] particleType Particle type for which the rate is requested (gluon, quark, anti-quark)
 * @param[in] collType Generic collision type for which the rate is requested (c22, c23, c32)
 * @param[in] unit (Optional, default = fm) Whether the result should be returned in 1/fm or in GeV
 * @return rate
 */
double ratesManager::getRate( const FLAVOR_TYPE particleType, const GENERIC_COLL_TYPE collType, const UNIT_TYPE unit ) const
{
  if ( !isNormalized )
  {
    std::string errMsg = "ratesManager::getRate called without prior normalization of rates. Unrecoverable error.";
    throw eRatesManager_error( errMsg );
  }
  else
  {
    double conversionFactor = 1;
    if ( unit == GeV )
    {
      conversionFactor = 0.197;
    }

    const int cIndex = static_cast<int>( collType );

    FLAVOR_TYPE pType = particleType;
    if ( pType == up || pType == down || pType == strange )
    {
      pType = quark;
    }
    else if ( pType == anti_up || pType == anti_down || pType == anti_strange )
    {
      pType = anti_quark;
    }

    switch ( pType )
    {
    case gluon:
      return ratesGluons[cIndex] * conversionFactor;
      break;
    case quark :
      return ratesQuarks[cIndex] * conversionFactor;
      break;
    case anti_quark:
      return ratesAntiQuarks[cIndex] * conversionFactor;
      break;
    default:
      std::string errMsg = "Particle type not found in ratesManager::getRate. Unrecoverable error.";
      throw eRatesManager_error( errMsg );
      break;
    }
  }
}



/**
 * Get the rate for a certain inclusive process (9231 etc.).
 * NB: this rate is NOT per particle as all other rates but per process
 *
 * @param[in] inclCollType Inclusive collision type (9221, 9326 etc)
 * @param[in] genCollType Generic collision type for which the rate is requested (c22, c23, c32)
 * @param[in] unit (Optional, default = fm) Whether the result should be returned in 1/fm or in GeV
 * @return Rate for the inclusive process. NOT per particle but per gq -> X etc..
 */
double ratesManager::getRateInclusive( const int inclCollType, const GENERIC_COLL_TYPE genCollType, const UNIT_TYPE unit ) const
{
  if ( !isNormalized )
  {
    std::string errMsg = "ratesManager::getRateInclusive called without prior normalization of rates. Unrecoverable error.";
    throw eRatesManager_error( errMsg );
  }
  else
  {
    double conversionFactor = 1;
    if ( unit == GeV )
    {
      conversionFactor = 0.197;
    }

    int index = interactionType::getIndexFromProcessType( inclCollType, genCollType );

    switch ( genCollType )
    {
    case c22:
      return rates22[index] * conversionFactor;
      break;
    case c23:
      return rates23[index] * conversionFactor;
      break;
    case c32:
      return rates32[index] * conversionFactor;
      break;
    default:
      std::string errMsg = "Generic type not found in ratesManager::getRateInclusive. Unrecoverable error.";
      throw eRatesManager_error( errMsg );
      break;
    }
  }
}







/**
 * Get the mean free path for a certain particle type
 *
 * @param[in] particleType Particle type for which the rate is requested (gluon, quark, anti-quark)
 * @param[in] unit (Optional, default = fm) Whether the result should be returned in 1/fm or in GeV
 * @return rate
 */
double ratesManager::getLambda( const FLAVOR_TYPE particleType, const UNIT_TYPE unit ) const
{
  double R = getRate( particleType, unit );

  if ( R > 0 )
  {
    return ( 1 / R );
  }
  else
  {
    return 0;
  }
}


/**
 * Get the mean free path for a certain particle type and a certain generic collision type (22, 23, 32)
 * e.g. to get the mean free path of quarks in 2->2 collisions, call ratesManager::getRate( quark, c22, fm )
 *
 * @param[in] particleType Particle type for which the rate is requested (gluon, quark, anti-quark)
 * @param[in] collType Generic collision type for which the rate is requested (c22, c23, c32)
 * @param[in] unit (Optional, default = fm) Whether the result should be returned in 1/fm or in GeV
 * @return rate
 */
double ratesManager::getLambda( const FLAVOR_TYPE particleType, const GENERIC_COLL_TYPE collType, const UNIT_TYPE unit ) const
{
  double R = getRate( particleType, collType, unit );

  if ( R > 0 )
  {
    return ( 1 / R );
  }
  else
  {
    return 0;
  }
}



void ratesManager::clear()
{
  isNormalized = false;
  rates22.assign( interactionType::indexProcessesInclusive22.size(), 0 );
  rates23.assign( interactionType::indexProcessesInclusive23.size(), 0 );
  rates32.assign( interactionType::indexProcessesInclusive32.size(), 0 );

  ratesGluons.assign( 3, 0 );
  ratesQuarks.assign( 3, 0 );
  ratesAntiQuarks.assign( 3, 0 );
  ratesGluonsTotal = 0;
  ratesQuarksTotal = 0 ;
  ratesAntiQuarksTotal = 0;
  
  nCollectedParticleBasedRates_gluons = 0;
  nCollectedParticleBasedRates_quarks = 0;
  nCollectedParticleBasedRates_antiQuarks = 0;
}



void ratesManager::addParticleBasedRates( const Particle& _particle, const UNIT_TYPE unit )
{
    std::string errMsg = "ratesManager::addParticleBasedRates should not be called";
    throw eRatesManager_error( errMsg );
  
//    double conversionFactor = 1;
//    if ( unit == GeV )
//    {
//     conversionFactor = 1 / 0.197;
//   }
//   
//   if ( !isNormalized )
//   {
//     std::string errMsg = "ratesManager::addParticleBasedRates only works for previously finalized instances";
//     throw eRatesManager_error( errMsg );
//   }
//   
//   if ( _particle.FLAVOR == gluon )
//   {
//     ratesGluons[0] += _particle.rate22 * conversionFactor;
//     ratesGluons[1] += _particle.rate23 * conversionFactor;
//     ratesGluons[2] += _particle.rate32 * conversionFactor;
//     ratesGluonsTotal += ( _particle.rate22 + _particle.rate23 + _particle.rate32 ) * conversionFactor;
//     ++nCollectedParticleBasedRates_gluons;
//   }
//   else
//   {
//     int _F = static_cast<int>( _particle.FLAVOR );
//     if ( _F % 2 == 1 )
//     {
//       ratesQuarks[0] += _particle.rate22 * conversionFactor;
//       ratesQuarks[1] += _particle.rate23 * conversionFactor;
//       ratesQuarks[2] += _particle.rate32 * conversionFactor;
//       ratesQuarksTotal += ( _particle.rate22 + _particle.rate23 + _particle.rate32 ) * conversionFactor;
//       ++nCollectedParticleBasedRates_quarks;
//     }
//     else
//     {
//       ratesAntiQuarks[0] += _particle.rate22 * conversionFactor;
//       ratesAntiQuarks[1] += _particle.rate23 * conversionFactor;
//       ratesAntiQuarks[2] += _particle.rate32 * conversionFactor;
//       ratesAntiQuarksTotal += ( _particle.rate22 + _particle.rate23 + _particle.rate32 ) * conversionFactor;
//       ++nCollectedParticleBasedRates_antiQuarks;
//     }    
//   }
}


void ratesManager::prepareParticleBasedAverages()
{
  if ( nCollectedParticleBasedRates_gluons > 0 )
  {
    ratesGluons[0] /= static_cast<double>( nCollectedParticleBasedRates_gluons );
    ratesGluons[1] /= static_cast<double>( nCollectedParticleBasedRates_gluons );
    ratesGluons[2] /= static_cast<double>( nCollectedParticleBasedRates_gluons );
    ratesGluonsTotal /= static_cast<double>( nCollectedParticleBasedRates_gluons );
  }
  else
  {
    ratesGluons[0] = 0;
    ratesGluons[1] = 0;
    ratesGluons[2] = 0;
    ratesGluonsTotal = 0;
  }
  
  if ( nCollectedParticleBasedRates_quarks > 0 )
  {
    ratesQuarks[0] /= static_cast<double>( nCollectedParticleBasedRates_quarks );
    ratesQuarks[1] /= static_cast<double>( nCollectedParticleBasedRates_quarks );
    ratesQuarks[2] /= static_cast<double>( nCollectedParticleBasedRates_quarks );
    ratesQuarksTotal /= static_cast<double>( nCollectedParticleBasedRates_quarks );
  }
  else
  {
    ratesQuarks[0] = 0;
    ratesQuarks[1] = 0;
    ratesQuarks[2] = 0;
    ratesQuarksTotal = 0;
  }
  
  if ( nCollectedParticleBasedRates_antiQuarks > 0 )
  {
    ratesAntiQuarks[0] /= static_cast<double>( nCollectedParticleBasedRates_antiQuarks );
    ratesAntiQuarks[1] /= static_cast<double>( nCollectedParticleBasedRates_antiQuarks );
    ratesAntiQuarks[2] /= static_cast<double>( nCollectedParticleBasedRates_antiQuarks );
    ratesAntiQuarksTotal /= static_cast<double>( nCollectedParticleBasedRates_antiQuarks );
  }
  else
  {
    ratesAntiQuarks[0] = 0;
    ratesAntiQuarks[1] = 0;
    ratesAntiQuarks[2] = 0;
    ratesAntiQuarksTotal = 0;
  }
}




ratesManager& ratesManager::operator+=( const ratesManager & rhs )
{
  if ( !( rhs.isNormalized && ( *this ).isNormalized ) )
  {
    std::string errMsg = "ratesManager::operator+= only works for previously finalized instances";
    throw eRatesManager_error( errMsg );
  }

  ratesGluonsTotal += rhs.ratesGluonsTotal;
  ratesQuarksTotal += rhs.ratesQuarksTotal;
  ratesAntiQuarksTotal += rhs.ratesAntiQuarksTotal;

  for ( int i = 0; i < 3; i++ )
  {
    ratesGluons[i] += rhs.ratesGluons[i];
    ratesQuarks[i] += rhs.ratesQuarks[i];
    ratesAntiQuarks[i] += rhs.ratesAntiQuarks[i];
  }

  return ( *this );
}




ratesManager& ratesManager::operator/=( const double & arg )
{
  if ( !isNormalized )
  {
    std::string errMsg = "ratesManager::operator/= only works for previously finalized instances";
    throw eRatesManager_error( errMsg );
  }

  ratesGluonsTotal /= arg;
  ratesQuarksTotal /= arg;
  ratesAntiQuarksTotal /= arg;

  for ( int i = 0; i < 3; i++ )
  {
    ratesGluons[i] /= arg;
    ratesQuarks[i] /= arg;
    ratesAntiQuarks[i] /= arg;
  }

  return ( *this );
}



ratesManager& ratesManager::operator*=( const double & arg )
{
  if ( !isNormalized )
  {
    std::string errMsg = "ratesManager::operator*= only works for previously finalized instances";
    throw eRatesManager_error( errMsg );
  }

  ratesGluonsTotal *= arg;
  ratesQuarksTotal *= arg;
  ratesAntiQuarksTotal *= arg;

  for ( int i = 0; i < 3; i++ )
  {
    ratesGluons[i] *= arg;
    ratesQuarks[i] *= arg;
    ratesAntiQuarks[i] *= arg;
  }

  return ( *this );
}


// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;
