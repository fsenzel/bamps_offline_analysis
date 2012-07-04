//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


#ifndef INITIALSHOWER_H
#define INITIALSHOWER_H

#include <stdexcept>
#include <vector>
#include <string>

#include "initialmodel.h"
#include "configuration.h"
#include "particle.h"
#include "vegas.h"
#include "woodsaxon.h"
#include "rangen_distr.h"
#include "pdfinterface.h"



/**
 * @brief Class to provide the shower initialization
 */
class initialShower
{
public:
  initialShower( const config& _config, WoodSaxon& _WoodSaxonParameter, const double _minimumPT, const int _nToGenerate = -1 );
  ~initialShower();
  
private:
  
  /** @brief Cut-off time for shower evolution */
  double insertionTime;

  void showerParticles( vector<Particle> _particles );
  vector<Particle> getShowerParticles( const double _px, const double _py, const double _pz1, const double _pz2, const FLAVOR_TYPE _flavor1, const FLAVOR_TYPE _flavor2 );
  vector<Particle> createShower(  int flavor1, int flavor2, double px, double py, double pz1, double pz2);
  
  int nEventsToGenerate;
  
  long int seed;
};

/** 
 * @brief exception class for handling unexpected critical behaviour
 * within generation of mini-jet initial distributions 
 */
class eShower_error : public std::runtime_error
{
  public:
    explicit eShower_error(const std::string& what) : std::runtime_error(what) {};
    
    virtual ~eShower_error() throw() {};
};


#endif // INITIALSHOWER_H
