//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------
//---------------------------------------------


#ifndef INITIALSTATEMODEL_H
#define INITIALSTATEMODEL_H

#include <vector>
#include "particle.h"


enum STORED_TABLE_USAGE { computeNewTables, useStoredTables };

/** @brief abstract base class for all models that generate initial particle distributions */
class initialStateModel
{
public:
  /** @brief Pure virtual function. Any derived class must at least specify this routine. */
  virtual void populateParticleVector( std::vector<Particle>& _particles ) = 0;
};

#endif // INITIALSTATEMODEL_H


// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on; 
