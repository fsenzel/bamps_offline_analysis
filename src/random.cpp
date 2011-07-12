//---------------------------------------------------------------------------------------
//provided by subversion
//---------------------------------------------------------------------------------------
//$HeadURL: file:///home/fochler/svnrepos/cascade_full/main/branches/restructure_branch/src/random.cpp $
//$LastChangedDate: 2009-07-19 15:31:52 +0200 (Sun, 19 Jul 2009) $
//$LastChangedRevision: 75 $
//$LastChangedBy: fochler $
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <time.h>
#include <boost/random.hpp>

#include "random.h"


/** @brief definition of global RNG object, defined extern in random.h */
randomNumberGenerator ran2;


/**
* This routine computes and sets a (hopefully) unique seed from either /dev/urandom provided by the system or the system time as a fallback.
* Using /dev/urandom is prefered and should work on all *nix system. /dev/urandom is a special device that holds (pseudo) random bits
* generated from user input, network traffic and black voodoo. Use /dev/random for even more unpredictable voodoo.
*
* @return The computed seed.
*/
uint32_t randomNumberGenerator::setSeed()
{
  uint32_t s;
  
  // open /dev/urandom as an ifstream
  std::ifstream urand("/dev/urandom", std::ios::binary );
  
  // get the seed from /dev/urandom, fall back to a seed generated from the current time if /dev/urandom could not have been opened.
  if ( urand.is_open() && urand.good() )
  {
    // read random characters from the input stream assigned to /dev/urandom
    // for this the integer type s is reinterpreted as a character pointer and sizeof(unit32_t) characters are read into s
    urand.read( reinterpret_cast<char *>(&s), sizeof(uint32_t) );
  }
  else
  {
    s = static_cast<uint32_t>( time(NULL) );
  }
  
  // use s obtained from /dev/urandom to seed the Boost RNG
  rng_mt_boost.seed( s );
  
  // return the seed for the user to have fun with
  return s;
}