//---------------------------------------------------------------------------------------
//provided by subversion
//---------------------------------------------------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------


#ifndef RATESMANAGER_H
#define RATESMANAGER_H


#include <vector>
#include "particle.h"
#include "interactiontype.h"


/**
 * @brief enum type for requests in terms of fm or GeV
 *
 * Mapped according to:
 * 0 = fm
 * 1 = GeV
 */
enum UNIT_TYPE {fm, GeV};



/**
 * @brief Class encapsulating the management of scattering rates
 *
 */
class ratesManager
{
public:
  /** @brief Constructor */
  ratesManager() : isNormalized( false ), rates22( interactionType::indexProcessesInclusive22.size(), 0 ),
      rates23( interactionType::indexProcessesInclusive23.size(), 0 ), rates32( interactionType::indexProcessesInclusive32.size(), 0 ),
      ratesGluons( 3, 0 ), ratesQuarks( 3, 0 ), ratesAntiQuarks( 3, 0 ),
      ratesGluonsTotal( 0 ), ratesQuarksTotal( 0 ), ratesAntiQuarksTotal( 0 ),
      nCollectedParticleBasedRates_gluons( 0 ), nCollectedParticleBasedRates_quarks( 0 ), nCollectedParticleBasedRates_antiQuarks( 0 ) {};
  /** @brief Destructor */
  ~ratesManager() {};
  
  
  ratesManager& operator+=( const ratesManager& rhs );
  ratesManager& operator/=( const double& arg );
  ratesManager& operator*=( const double& arg );
  
  const ratesManager operator+(const ratesManager& rhs) const { return ratesManager(*this) += rhs; }
  const ratesManager operator/(const double& arg) const { return ratesManager(*this) /= arg; }
  const ratesManager operator*(const double& arg) const { return ratesManager(*this) *= arg; }
  
  

  /** @brief Clears all stored rates etc. */
  void clear();
  
  /** @brief Add sampled probability for 2->2 and 2->3 processes */
  int add( const GENERIC_COLL_TYPE genericType, const FLAVOR_TYPE _F1, const FLAVOR_TYPE _F2, const double _P );
  /** @brief Add sampled probability for 3->2 processes */
  int add( const GENERIC_COLL_TYPE genericType, const FLAVOR_TYPE _F1, const FLAVOR_TYPE _F2, const FLAVOR_TYPE _F3, double _P );
  
  void addParticleBasedRates( const Particle& _particle, const UNIT_TYPE unit = GeV );
  void prepareParticleBasedAverages();
  
  /** @brief Normalize the stored rates such that the interaction rate per particle and per time interval can be given */
  int normalizeRates( const int Ng, const int Nq, const int Nqbar, const double dt );
  /** @brief Set all rates to 0 and set normalized flag */
  int normalizeRates();
  
  
  /** @brief Get the rate for a given particle type */
  double getRate( const FLAVOR_TYPE particleType, const UNIT_TYPE unit = fm ) const;
  /** @brief Get the rate for a given particle type in given collision types (22, 23, 32) */
  double getRate( const FLAVOR_TYPE particleType, const GENERIC_COLL_TYPE collType, const UNIT_TYPE unit = fm ) const;
  /** @brief Get the rate per inclusie process */
  double getRateInclusive( const int inclCollType, const GENERIC_COLL_TYPE genCollType, const UNIT_TYPE unit = fm ) const;

  /** @brief Get the mean free path for a given particle type */
  double getLambda( const FLAVOR_TYPE particleType, const UNIT_TYPE unit = fm ) const;
  /** @brief Get the mean free path for a given particle type in given collision types (22, 23, 32) */
  double getLambda( const FLAVOR_TYPE particleType, const GENERIC_COLL_TYPE collType, const UNIT_TYPE unit = fm ) const;


  /** @brief Switch indicating whether the collection phase of probabilities has been completed and the rates have been normalized */
  bool isNormalized;

private:
  /** @brief Collect the probabilities for all 2->2 interaction types (each index corresponds to a interaction type as set in class interactionType) */
  std::vector<double> rates22;
  /** @brief Collect the probabilities for all 2->3 interaction types (each index corresponds to a interaction type as set in class interactionType) */
  std::vector<double> rates23;
  /** @brief Collect the probabilities for all 2->3 interaction types (each index corresponds to a interaction type as set in class interactionType) */
  std::vector<double> rates32;

  /** @brief Stores rates for gluons after ratesManager::normalizeRates has been called. Entries are: 0 = 2->2, 1 = 2->3, 2 = 3->2 */
  std::vector<double> ratesGluons;
  /** @brief Stores rates for quarks after ratesManager::normalizeRates has been called. Entries are: 0 = 2->2, 1 = 2->3, 2 = 3->2 */
  std::vector<double> ratesQuarks;
  /** @brief Stores rates for anti-quarks after ratesManager::normalizeRates has been called. Entries are: 0 = 2->2, 1 = 2->3, 2 = 3->2 */
  std::vector<double> ratesAntiQuarks;

  /** @brief Stores total (22 + 23 + 32) rate for gluons  after ratesManager::normalizeRates has been called. */
  double ratesGluonsTotal;
  /** @brief Stores total (22 + 23 + 32) rate for quarks after ratesManager::normalizeRates has been called. */
  double ratesQuarksTotal;
  /** @brief Stores total (22 + 23 + 32) rate for anti-quarks  after ratesManager::normalizeRates has been called. */
  double ratesAntiQuarksTotal;
  
  int nCollectedParticleBasedRates_gluons;
  int nCollectedParticleBasedRates_quarks;
  int nCollectedParticleBasedRates_antiQuarks;  
};



/** @brief exception class for handling unexpected behaviour when managing interaction rates */
class eRatesManager_error : public std::runtime_error
{
public:
  explicit eRatesManager_error( const std::string& what ) : std::runtime_error( what ) {};

  virtual ~eRatesManager_error() throw() {};
};

#endif // RATESMANAGER_H
// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;
