//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL: svn+ssh://fochler@th.physik.uni-frankfurt.de/home/fochler/svnrepos/cascade_full/main/branches/LHC_initial_conditions/src/ringcontainer.h $
//$LastChangedDate: 2010-12-27 18:48:53 +0100 (Mon, 27 Dec 2010) $
//$LastChangedRevision: 256 $
//$LastChangedBy: fochler $
//---------------------------------------------
//---------------------------------------------



#ifndef RINGCONTAINER_H
#define RINGCONTAINER_H

#include "ratesmanager.h"
#include "particle.h"


class ringContainer
{
public:
  ringContainer();
  ringContainer( const double _minR, const double _maxR );
  ~ringContainer() {};
  
  void clear();
  void relocate( const double _minR, const double _maxR );
  
  void addParticle( const Particle& _particle );
  void addParticleInFormGeom( const Particle& _particle, const double _time );
  void addRates( const Particle& _particle );
  
  double getAveraged_md2g() const; 
  double getAveraged_md2q() const;
  
  double getAveraged_v_x() const { if ( numberOfParticles > 0 ) return ( v_x / numberOfParticles ); else return 0; } 
  double getAveraged_v_y() const { if ( numberOfParticles > 0 ) return ( v_y / numberOfParticles ); else return 0; } 
  double getAveraged_v_z() const { if ( numberOfParticles > 0 ) return ( v_z / numberOfParticles ); else return 0; } 
  double getAveraged_v_r() const { if ( numberOfParticles > 0 ) return ( v_r / numberOfParticles ); else return 0; } 
  double getGamma() const;
  
  double getAccumulated_E() const { if ( numberOfParticles > 0 ) return ( E ); else return 0; } 
  double getAveraged_E() const { if ( numberOfParticles > 0 ) return ( E / numberOfParticles ); else return 0; } 
  double getAveraged_p_z() const { if ( numberOfParticles > 0 ) return ( p_z / numberOfParticles ); else return 0; } 
  double getAveraged_p_r() const { if ( numberOfParticles > 0 ) return ( p_r / numberOfParticles ); else return 0; } 
  double getAveraged_pr2_over_E() const { if ( numberOfParticles > 0 ) return ( pr2_over_E / numberOfParticles ); else return 0; } 
  double getAveraged_pz2_over_E() const { if ( numberOfParticles > 0 ) return ( pz2_over_E / numberOfParticles ); else return 0; } 
  double getAveraged_pr_pz_over_E() const { if ( numberOfParticles > 0 ) return ( pr_pz_over_E / numberOfParticles ); else return 0; } 
  
  double getVolume( const double _dz) const;
  
  void prepareAverages( const double _dz, const int _Ntest );
  
  double getParticleDensity() const;
  double getGluonDensity() const;
  double getQuarkDensity() const;
  
  double getEnergyDensity() const;
  double getEffectiveTemperature() const;
  
  double transformEnergyToComovingFrame( double _P[4] ) const;
  
  double minRadius;
  double maxRadius;
  double deltaR;
  
  ratesManager rates;
    
    
// private:
  bool averagesPrepared;
  
  int numberOfParticles;
  int numberOfGluons;
  int numberOfQuarks;
  int numberOfActiveParticles;

  double md2g;
  double md2q;
  
  double v_x;
  double v_y;
  double v_z;
  double v_r;
  
  double E;
  double inverseE_gluons;
  double inverseE_quarks;
  double p_z;
  double p_r;
  double pr2_over_E;
  double pz2_over_E;
  double pr_pz_over_E;
  
  double gamma;
  double energyDensity;
  double particleDensity;
  double gluonDensity;
  double quarkDensity;
  
  double volume;

  int numberOfCollectedRateObjects;
};



/** @brief exception class */
class eRingContainer_error : public std::runtime_error
{
  public:
    explicit eRingContainer_error(const std::string& what) : std::runtime_error(what) {};
    
    virtual ~eRingContainer_error() throw() {};
};


#endif // RINGCONTAINER_H
