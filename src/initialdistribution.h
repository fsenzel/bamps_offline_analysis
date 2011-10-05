//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/restructureOutput/src/initialdistribution.h $
//$LastChangedDate: 2011-07-12 17:55:28 +0200 (Tue, 12 Jul 2011) $
//$LastChangedRevision: 2 $
//$LastChangedBy: fochler $
//---------------------------------------------
//---------------------------------------------



#ifndef INITIALDISTRIBUTION_H
#define INITIALDISTRIBUTION_H

#include <stdexcept>
#include <vector>

#include "configuration.h"
#include "particle.h"
#include "woodsaxon.h"

// enum INITIAL_STATE_TYPE { miniJetsInitialState, pythiaInitialState, cgcInitialState };


class initialDistribution
{
public:
  initialDistribution( const config* const _config, const INITIAL_STATE_TYPE _initialStateType = miniJetsInitialState );
  ~initialDistribution() {};

  void populateParticleVector( std::vector< ParticleOffline >& _particles, WoodSaxon& _wsParameter );
  

private:
  const config * const configObject;
  WoodSaxon WoodSaxonParameter;
  double A;          //mass number of nucleus A
  double Aatomic;    //atomic number, i.e. number of protons, of nucleus A
  double B;          //mass number of nucleus B
  double Batomic;    //atomic number of nucleus B
  double sqrtS;      //c.m. energy per NN pair, GeV
  double impactParameter;       //impact parameter in fm
  int numberOfTestparticles;    //number of testparticles per real particle
  INITIAL_STATE_TYPE initialStateType;
};


/** @brief exception class for handling unexpected critical behaviour within generation of mini-jet initial distributions  */
class eInitialState_error : public std::runtime_error
{
public:
  explicit eInitialState_error( const std::string& what ) : std::runtime_error( what ) {};

  virtual ~eInitialState_error() throw() {};
};

#endif // INITIALDISTRIBUTION_H
