//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------
//---------------------------------------------

// at revision 902, this file is identical to full/branches/vector4D/src/...
//
// the method 'setLongitudinalGeometry' is new!!!
//
// The original version was written for 'ParticleOffline', but since
// this is derived from 'Particle', we can match it to the latter.

#include <vector>
#include <string>
#include <math.h>


#include "ringstructure.h"
#include "ringcontainer.h"


ringStructure::ringStructure( const int _nRings, const double _centralRadius, const double _deltaR ) : numberOfRings( _nRings ),
    centralRingRadius( _centralRadius ), deltaR( _deltaR )
{
  rings.resize( _nRings );

  rings[0].relocate( 0, _centralRadius );

  totalRadius = _centralRadius + ( _nRings - 1 ) * deltaR;
  for ( unsigned int i = 1; i < rings.size(); i++ )
  {
    rings[i].relocate( rings[i-1].maxRadius, rings[i-1].maxRadius + deltaR );
  }
}



void ringStructure::resize( const int _nRings, const double _centralRadius, const double _deltaR )
{
  numberOfRings = _nRings;
  centralRingRadius = _centralRadius;
  deltaR = _deltaR;
  
  rings.clear();
  rings.resize( _nRings );

  rings[0].relocate( 0, _centralRadius );

  totalRadius = _centralRadius + ( _nRings - 1 ) * deltaR;
  for ( unsigned int i = 1; i < rings.size(); i++ )
  {
    rings[i].relocate( rings[i-1].maxRadius, rings[i-1].maxRadius + deltaR );
  }
}

void ringStructure::setLongitudinalGeometry(const double _y_left, const double _y_right, const double _t)
{
  y_left = _y_left;
  y_right = _y_right;
  timenow = _t;
  delta_z = timenow * ( tanh( y_right ) - tanh( y_left ) );
}


int ringStructure::getIndex( const double _xt ) const
{
  int index = getIndexPure( _xt );
  if ( index >= rings.size() )
  {
    return ( static_cast<int>( rings.size() ) - 1 );
  }
  else
  {
    return index;
  }
}



int ringStructure::getIndexPure( const double _xt ) const
{
  if ( _xt < 0 )
  {
    std::string errMsg = "transverse position xt < 0";
    throw eRingStructure_error( errMsg );
  }
  
  if ( _xt < centralRingRadius )
  {
    return 0;
  }
  else
  {
//     if ( _xt > totalRadius )
//     {
//       std::string errMsg = "transverse position xt > R (R = total radius of ring structure)";
//       throw eRingStructure_error( errMsg );
//     }
    
    int index = static_cast<int>( ( _xt - centralRingRadius ) / deltaR ) + 1;  // +1 since index 0 is for the central ring
    return index;
  }
}





int ringStructure::getIndex( const Particle& _particle ) const
{
  double xt = _particle.Pos.Perp();
  return getIndex( xt );
}


ringContainer& ringStructure::operator[]( const int index )
{
  if ( index < 0 || index >= numberOfRings )
  {
    std::string errMsg = "index out of range in ringStructure";
    throw eRingStructure_error( errMsg );
  }
  
  return rings[ index ];
}



void ringStructure::addRates( const Particle& _particle )
{
  double xt = _particle.Pos.Perp();
  addRates( xt, _particle );
}


void ringStructure::addRates( const double _xt, const Particle& _particle )
{
  int index = getIndexPure( _xt );
  if ( index < static_cast<int>( rings.size() ) )
  {
    rings[ index ].addRates( _particle );
  }
}


void ringStructure::addParticle( const Particle& _particle )
{
  double xt = _particle.Pos.Perp();
  addParticle( xt, _particle );
}


void ringStructure::addParticle( const double _xt, const Particle& _particle )
{
  if( _particle.FLAVOR > 2*Particle::max_N_light_flavor )
  {
    std::string errMsg = "Heavy quark added to ringStructure which cannot deal with massive particles";
    throw eRingStructure_error( errMsg );
  }
  else
  {
    int index = getIndexPure( _xt );
    if ( index < static_cast<int>( rings.size() ) )
    {
      rings[ index ].addParticle( _particle );
    }
  }
}

void ringStructure::addParticleInFormGeom( const double _xt, const Particle& _particle, const double _time )
{
  if( _particle.FLAVOR > 2*Particle::max_N_light_flavor )
  {
    std::string errMsg = "Heavy quark added to ringStructure which cannot deal with massive particles";
    throw eRingStructure_error( errMsg );
  }
  else
    rings[ getIndex( _xt ) ].addParticleInFormGeom( _particle, _time );
}


void ringStructure::addParticleInFormGeom( const Particle& _particle, const double _time )
{
  if( _particle.FLAVOR > 2*Particle::max_N_light_flavor )
  {
    std::string errMsg = "Heavy quark added to ringStructure which cannot deal with massive particles";
    throw eRingStructure_error( errMsg );
  }
  else
    rings[ getIndex( _particle ) ].addParticleInFormGeom( _particle, _time );
}


void ringStructure::clear()
{
  for ( unsigned int i = 0; i < rings.size(); i++ )
  {
    rings[i].clear();
  }
}


void ringStructure::prepareAverages( const double _dz, const int _Ntest )
{
  for ( unsigned int i = 0; i < rings.size(); i++ )
  {
    rings[i].prepareAverages( _dz, _Ntest );
  }
}


// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;
