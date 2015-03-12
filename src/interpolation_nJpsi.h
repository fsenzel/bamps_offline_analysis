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


#ifndef INTERPOLATION_NJPSI_H
#define INTERPOLATION_NJPSI_H

#include "interpolation_n_dimensions.h"

using std::string;


class interpolation_nJpsi : interpolation1d // based on interpolation2d
{
  public:
    void configure();

    double getN(const double T) const;

  private:

};

#endif
