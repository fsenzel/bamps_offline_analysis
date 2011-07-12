//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL: file:///home/fochler/svnrepos/cascade_full/main/trunk/src/lorentz.h $
//$LastChangedDate: 2008-02-22 14:53:42 +0100 (Fri, 22 Feb 2008) $
//$LastChangedRevision: 33 $
//$LastChangedBy: fochler $
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------


#ifndef LORENTZ_H
#define LORENTZ_H

#include "particle.h"

void lorentz(const double beta[4], const double A[4], double LA[4]);

#endif
