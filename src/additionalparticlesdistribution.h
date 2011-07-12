//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL: file:///home/fochler/svnrepos/cascade_full/main/branches/restructure_branch/src/initialdistribution.h $
//$LastChangedDate: 2010-07-13 00:25:52 +0200 (Tue, 13 Jul 2010) $
//$LastChangedRevision: 126 $
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

enum INITIAL_STATE_TYPE { miniJetsInitialState, pythiaInitialState, cgcInitialState };


class additionalParticlesDistribution
{
public:
  additionalParticlesDistribution( const config*const _config, const INITIAL_STATE_TYPE _initialStateType = miniJetsInitialState );
  ~additionalParticlesDistribution() {};

  void populateParticleVector( std::vector< Particle >& _particles, WoodSaxon& _wsParameter );
  void prepareParticles( std::vector< Particle >& _particles );
 

private:
  int numberOfParticlesToAdd;
  double minimumPT;
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
