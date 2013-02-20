//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL: file:///home/fochler/svnrepos/cascade_box/main/branches/quark_branch/src/interpolation_gc.h $
//$LastChangedDate: 2009-06-23 12:06:07 +0200 (Tue, 23 Jun 2009) $
//$LastChangedRevision: 90 $
//$LastChangedBy: fochler $
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------

/** @file
 * @brief Declarations for the class interpolation_iniQNLO.
 */


#ifndef INTERPOLATION_NJPSI_H
#define INTERPOLATION_NJPSI_H

#include "interpolation_n_dimensions.h"

using std::string;


class interpolation_nJpsi : interpolation1d // based on interpolation2d
{
  public:
    void configure();

    double getN ( const double T ) const;

  private:

};

#endif
