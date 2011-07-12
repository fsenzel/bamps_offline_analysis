//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL: file:///home/fochler/svnrepos/cascade_full/main/branches/LHC_initial_conditions/src/minijets.h $
//$LastChangedDate: 2010-12-14 11:09:56 +0100 (Tue, 14 Dec 2010) $
//$LastChangedRevision: 241 $
//$LastChangedBy: fochler $
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
