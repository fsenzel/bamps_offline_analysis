//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------

/** @file 
 * @brief Declarations for the class interpolation_iniQNLO.
 */


#ifndef INTERPOLATIONINIJPSI_H
#define INTERPOLATIONINIJPSI_H

#include "interpolation_n_dimensions.h"

using std::string;

enum shadowModelJpsi {none, eks98, eps08, eps09};


class interpolation_iniJpsi_dndptdy : interpolation2d // based on interpolation2d
{
  public:
    void configure( const double sqrtS_arg, const double Bimp_arg, const double sigmaAbs_arg = 1.5, const double agN = 0.1, const shadowModelJpsi shadowing_model = eps08 );

    double getdN(const double y, const double pt) const;

  private:
    double impact_parameter, sigmaAbs, sqrtS;

};

#endif
