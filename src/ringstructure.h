//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

// at revision 902, this file is identical to full/branches/vector4D/src/ringstructure.h
//
// the method 'setLongitudinalGeometry' is new !!!

#ifndef RINGSTRUCTURE_H
#define RINGSTRUCTURE_H

#include <vector>
#include <string>
#include <stdexcept>

#include "ringcontainer.h"
#include "particle.h"
#include "ratesmanager.h"



class ringStructure
{
public:
  ringStructure() : 
    delta_z(0), 
    y_left(0), y_right(0), 
    timenow(0),
    numberOfRings( 0 ), 
    centralRingRadius( 0 ), 
    totalRadius( 0 ), 
    deltaR( 0 )
  { 
    rings.resize(0); 
  }
  
  ringStructure( const int _nRings, const double _centralRadius, const double _deltaR );

  ~ringStructure() {};
  
  void resize( const int _nRings, const double _centralRadius, const double _deltaR );
  
  ringContainer& getRing( const double _xt ) { return rings[ getIndex( _xt ) ]; }
  ringContainer& operator[]( const int _index );

  void setLongitudinalGeometry( const double _y_left, const double _y_right, const double _t );

  
  int getIndex( const double _xt ) const;
  int getIndexPure( const double _xt ) const;
  int getIndex( const Particle& _particle ) const;
  int size() const { return numberOfRings; }
  double getCentralRadius() const { return centralRingRadius; }
  double getDeltaR() const { return deltaR; }
  
  void clear();  
  
  void addParticle( const double _xt, const Particle& _particle );
  void addParticle( const Particle& _particle );
  void addParticleInFormGeom( const double _xt, const Particle& _particle, const double _time );
  void addParticleInFormGeom( const Particle& _particle, const double _time );
  void addRates( const double _xt, const Particle& _particle );
  void addRates( const Particle& _particle );
  
  void prepareAverages( const double _dz, const int _Ntest );
 
  
  double delta_z;
  double y_left, y_right;
  double timenow;
  

 
private:
  std::vector<ringContainer> rings;
  
  int numberOfRings;
  double centralRingRadius;
  double totalRadius;
  double deltaR;
};




/** @brief exception class */
class eRingStructure_error : public std::runtime_error
{
  public:
    explicit eRingStructure_error(const std::string& what) : std::runtime_error(what) {};
    
    virtual ~eRingStructure_error() throw() {};
};


#endif // RINGSTRUCTURE_H
// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on; 
