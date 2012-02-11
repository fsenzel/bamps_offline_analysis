#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <math.h>

#include "interpolation_iniJpsi.h"

using namespace std;

void interpolation_iniJpsi_dndptdy::configure( const double sqrtS_arg, const double Bimp_arg, const double sigmaAbs_arg, const double agN, const shadowModelJpsi shadowing_model ) 
{
  string name;
  
  sqrtS = sqrtS_arg;
  impact_parameter = Bimp_arg;
  sigmaAbs = sigmaAbs_arg;
  
  // total numbers of tabulated points on the a, b "axes"
  n_i = 21; // y
  n_j = 21; // pt
  
  // spacings of tabulated values in a, b-direction
  delta_a = 0.25;
  delta_b = 0.25;

  // the smallest tabulated values of a, b, c and d
  a_start = 0.0;
  b_start = 0.0;
 
  
  name = "ini_jpsi";
  
  if( sqrtS == 200.0)
    name = name + "_rhic200";
  else if( sqrtS == 5500.0)
    name = name + "_lhc55";
  else if( sqrtS == 2760.0)
    name = name + "_lhc276";
  
  if( impact_parameter == 0.0)
    name = name + "_b00";
  else if( impact_parameter == 3.3)
    name = name + "_b33";
  else if( impact_parameter == 4.6)
    name = name + "_b46";
  else if( impact_parameter == 5.8)
    name = name + "_b58";
  else if( impact_parameter == 8.2)
    name = name + "_b82";
  else if( impact_parameter == 10.3)
    name = name + "_b103";
  
  if( sigmaAbs == 2.8)
    name = name + "_sigmaAbs28";
  else if( sigmaAbs == 1.5)
    name = name + "_sigmaAbs15";
  else if( sigmaAbs == 0.0)
    name = name + "_sigmaAbs0";
  
  if( agN == 0.1 )
    name = name + "_agN01";
  else if( agN == 0.0 )
    name = name + "_agN0";
  
  if( shadowing_model == eps08 )
    name = name + "_eps08";
  else if( shadowing_model == none )
    name = name + "_noshad";
  
  filename = "data/"+name+".dat";
  
  loadData();
}

// not properly normalized
double interpolation_iniJpsi_dndptdy::getdN(const double y, const double pt) const
{
  double cs = getInterpolatedData( y, pt);
  
  if(cs >= 0.0)
    return cs;
  else
    return 0.0; 
}

