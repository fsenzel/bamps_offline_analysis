//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------
//---------------------------------------------

#include <vector>
#include <string>
#include <math.h>


#include "ringstructure.h"
#include "ringcontainer.h"


ringStructure::ringStructure( const int _nRings, const double _centralRadius, const double _deltaR ) : numberOfRings( _nRings ),
    centralRingRadius( _centralRadius ), deltaR( _deltaR ), y_left(0), y_right(0), delta_z(0), timenow(0)
{
  rings.resize( _nRings );

  rings[0].relocate( 0, _centralRadius );

  totalRadius = _centralRadius + ( _nRings - 1 ) * deltaR;
  for ( int i = 1; i < rings.size(); i++ )
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
  for ( int i = 1; i < rings.size(); i++ )
  {
    rings[i].relocate( rings[i-1].maxRadius, rings[i-1].maxRadius + deltaR );
  }
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


void ringStructure::setLongitudinalGeometry(const double _y_left, const double _y_right, const double _t)
{
  y_left = _y_left;
  y_right = _y_right;
  timenow = _t;
  delta_z = timenow * ( tanh( y_right ) - tanh( y_left ) );
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
  double xt = sqrt( pow( _particle.X, 2 ) + pow( _particle.Y, 2 ) );
  
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
  double xt = sqrt( pow( _particle.X, 2 ) + pow( _particle.Y, 2) );
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
  double xt = sqrt( pow( _particle.X, 2 ) + pow( _particle.Y, 2) );
  addParticle( xt, _particle );
}


void ringStructure::addParticle( const double _xt, const Particle& _particle )
{
  int index = getIndexPure( _xt );
  if ( index < static_cast<int>( rings.size() ) )
  {
    rings[ index ].addParticle( _particle );
  }
}


void ringStructure::clear()
{
  for ( int i = 0; i < rings.size(); i++ )
  {
    rings[i].clear();
  }
}


void ringStructure::prepareAverages( const double _dz, const int _Ntest )
{
  for ( int i = 0; i < rings.size(); i++ )
  {
    rings[i].prepareAverages( _dz, _Ntest );
  }
}


// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;
