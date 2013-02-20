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


#ifndef INTERPOLATIONINIJPSI_H
#define INTERPOLATIONINIJPSI_H

/**
@author Jan Uphoff and Oliver Fochler
 */


#include "interpolation_n_dimensions.h"

using std::string;

enum shadowModelJpsi {none, eks98, eps08, eps09};


class interpolation_iniJpsi_dndptdy : interpolation2d // based on interpolation2d
{
  public:
    void configure ( const double sqrtS_arg, const double Bimp_arg, const double sigmaAbs_arg = 1.5, const double agN = 0.1, const shadowModelJpsi shadowing_model = eps08 );

    double getdN ( const double y, const double pt ) const;

  private:
    double impact_parameter, sigmaAbs, sqrtS;

};

#endif
