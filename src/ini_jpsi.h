#ifndef INI_JPSI_H
#define INI_JPSI_H

#include "configuration.h"
#include "interpolation_iniJpsi.h"



class ini_jpsi
{
public:
    ini_jpsi(const double sqrtS_arg, const double Bimp_arg, const double sigmaAbs_arg, const double agN, const shadowModelJpsi shadowing_model, const double KInicharm_arg  );
    void sample_jpsis();

  private:
    void sample_metropolis_dndptdy(double& pt_arg, double& y_arg);
    void sample_one_jpsi(const int nmb);
    
    interpolation_iniJpsi_dndptdy theInterpolation_dndptdy;
    
    double impact_parameter, sigmaAbs, sqrtS, agN;
    shadowModelJpsi shadowing_model;
    
    double KInicharm;
    
};

#endif // INI_HQ_NLO_H
