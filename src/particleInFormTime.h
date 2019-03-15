//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/offlineAnalysis/branches/stochasticLPM/src/particleOffline.h $
//$LastChangedDate: 2019-03-06 15:54:36 +0100 (Wed, 06. Mar 2019) $
//$LastChangedRevision: 2944 $
//$LastChangedBy: senzel $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file
 * @brief Declarations for the Offline Particle class
 */

#ifndef PARTICLEINFORMTIME_H
#define PARTICLEINFORMTIME_H

#include <vector>
#include <limits>

#include "particle.h"
#include "globalsettings.h"

class ParticleInFormTime : public Particle
{
public:
  // /** @brief Provide standard constructor (for completeness) */
  // ParticleInFormTime( const double ) : 
  //   Particle(), 
  //   tInitEmission( -1.0 ),
  //   deltaMomMother( VectorEPxPyPz( 0.0, 0.0, 0.0, 0.0 ) ),
  //   MomAtEmission( VectorEPxPyPz( 0.0, 0.0, 0.0, 0.0 ) )
  //   {};
    
  ParticleInFormTime( const Particle& _particle, const double _tInitEmission, const VectorEPxPyPz _deltaMomMother, const VectorEPxPyPz _momAtEmission, const VectorEPxPyPz _momMotherAtEmission ) : 
    Particle( _particle ), 
    tInitEmission( _tInitEmission ),
    deltaMomMother( _deltaMomMother ),
    MomAtEmission( _momAtEmission ),
    MomMotherAtEmission( _momMotherAtEmission ),
    NscattDuringTau( 1.0 )
  {};
    
  double NscattDuringTau;

  VectorEPxPyPz getDeltaMomMother() const { return deltaMomMother; };
  VectorEPxPyPz getMomAtEmission() const { return MomAtEmission; };
  VectorEPxPyPz getMomMotherAtEmission() const { return MomMotherAtEmission; };

  double getFormationTimeCoherentState( const Particle _mother ) const
  {
    const double kt2 = this->Mom.TransverseMomentumToVectorSquared( _mother.Mom );

    if( kt2 == 0.0 )
      return std::numeric_limits<double>::infinity();
    else
    {
      return this->tInitEmission + ( this->Mom.E() / kt2 * ns_casc::InvGevToFm ); // fm
    }
  };

private:

  double tInitEmission;
  VectorEPxPyPz deltaMomMother;
  VectorEPxPyPz MomAtEmission;
  VectorEPxPyPz MomMotherAtEmission;
  
};


#ifdef BAMPS_DECLARE_ALLOCATOR
BAMPS_DECLARE_ALLOCATOR(ParticleInFormTime);
#else
#warning "ALLOCATOR not defined"
#endif

#endif
