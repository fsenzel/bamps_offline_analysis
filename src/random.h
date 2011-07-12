//---------------------------------------------------------------------------------------
//provided by subversion
//---------------------------------------------------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------


#ifndef RANDOM_NG_H
#define RANDOM_NG_H

#include <stdint.h>
#include <boost/random.hpp>

/** 
 * @brief This class encapsulates the generation of random number
 *
 * This class provides a consistent interface to possibly different random number generators (RNG) and 
 * also handles the setting of seeds for these generators.
 */
class randomNumberGenerator
{
  public:
    /** 
     * @brief Standard constructor 
     * * This standard constructor initializes the boost::uniform_01 object with the RNG and sets a standard seed.
     */
    randomNumberGenerator() : ran_boost( rng_mt_boost ) { mySeed = setSeed(); }
    /**
    * @brief Constructor that already sets the seed
    * This constructor initializes the boost::uniform_01 object with the RNG and sets a user provided seed.
    */
    randomNumberGenerator( const uint32_t s ) : mySeed( s ), ran_boost( rng_mt_boost ) { rng_mt_boost.seed( mySeed ); }
    
    /** @brief Compute a nice seed and set it */
    uint32_t setSeed();
    /** @brief Set seed provided by user */
    void setSeed( const uint32_t s )  { rng_mt_boost.seed( s );  mySeed = s; }
    
    /** 
    * @brief Get random number
    * @return random floating point number uniformly distributed in [0,1)
    */
    double getRan() { return this->operator()(); }
    
    /**
    * @brief Get random number via overloaded ()-operator
    * @return random floating point number uniformly distributed in [0,1)
    */
    double operator()() { return ran_boost(); }
    
  private:
    /** @brief Stores the seed for internal use */
    uint32_t mySeed;
    
    /**
    * @brief A random number generator (RNG) from the Boost libraries 
    *
    * The Boost RNG is used to obtain random numbers (integer type). It must be seeded prior to use, c.f. setSeed().
    * boost::mt199937 uses a Mersenne twister algorithm.
    * Other reasonable choices for the RNG would be boost::lagged_fibonacci607 or boost::kreutzer1986.
    * See http://www.boost.org/doc/libs/release/libs/random
    */
    boost::mt19937 rng_mt_boost;
    
    
    /** 
    * @brief A distribution generator from the Boost libraries 
    *
    * This maps the random numbers from a Boost RNG (type boost:mt19937 in this case) to uniformly distributed floating point
    * numbers in the interval [0,1).
    * The pointer needs to be initialized in the constructors with the actual RNG (rng_mt_boost in our case).
    */
    boost::uniform_01<boost::mt19937&, double> ran_boost;
};



/** @brief Random number generator (RNG) object, GLOBAL 
*
* Declared as extern here. There MUST be ONE definition in a cpp-file that includes this header in the global 
* scope (it's in random.cpp).
* All files that include this header have then access to the random number generator.
*/
extern randomNumberGenerator ran2;

#endif
