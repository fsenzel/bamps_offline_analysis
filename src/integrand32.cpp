//---------------------------------------------------------------------------------------
//provided by subversion
//---------------------------------------------------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------


#include <math.h>

#include "integrand32.h"
#include "FPT_compare.h"



void integrand32::set_beta_vec(const double x_[])
{
  beta_vec[0] = 0;
  for(int i=1;i<=3;i++)
  {
    beta_vec[i] = x_[i];
    //beta_vec[i] = 0;  // for testing consistency with old implementation (also see scattering32.cpp)
  }
  
  beta_abs = sqrt ( pow(beta_vec[1],2.0) + pow(beta_vec[2],2.0) + pow(beta_vec[3],2.0) );
}


// overloaded operator() that makes integrand32 a function object - for use with CUBA-vegas
void integrand32::operator()(const int *ndim, const double xx[], const int *ncomp, double ff[]) const
{
  double wgt;
  ff[0] = this->operator()(xx, wgt);  
}


// overloaded operator() that makes integrand32 a function object - for use with NR-vegas
double integrand32::operator()(const double xx[], double wgt) const
{
  double sin_gamma = sqrt( 1.0 - pow(cos_gamma_int,2.0) );   //sin(gamma)
  
  double u = xx[1];                       // cos(theta)                                   
  double us = sqrt( 1.0 - pow(u,2.0) );  // sin(theta)
  
  double V_phi = M_PI;  // volume for phi integration, as phi is sampled from 0 to M_PI rather than from 0 to 1
  double phi = M_PI*xx[2];
  double v = cos(phi);
  
  double cos_delta = sin_gamma * us * v + cos_gamma_int * u;   // cos(delta) = sin(gamma)*sin(theta)*cos(phi) + cos(gamma)*cos(theta)
  double sin_delta_2 = 1.0 - pow( cos_delta,2.0 );             // sin(delta)^2 = 1 - cos(delta)^2
  
  double qt2 = pow(E1_int,2.0) * pow(us,2.0);                  // qt^2 = E1^2 * sin(theta)^2 
  double kt2 = pow(E3_int,2.0) * sin_delta_2;                  // kt^2 = E3^2 * sin(delta)^2 
  //vector product qt*kt = -E1*E3*sin(theta)*( cos(gamma)*sin(theta) - sin(gamma)*cos(theta)*cos(phi) )
  double qtkt = -E1_int * E3_int * us * ( cos_gamma_int*us - sin_gamma*u*v );  
  
  // Theta (capital!) is the angle between the boost-vector beta' (here: beta_vec) and the vector p1'
  double cos_Theta;
  if ( FPT_COMP_E(beta_abs,0.0) )
    cos_Theta = 0;
  else
    cos_Theta = 1/beta_abs * ( beta_vec[1]*us*v + beta_vec[2]*us*sqrt(1-pow(v,2.0)) + beta_vec[3]*u );

  
  // the constraint for kt^2 depending on E3, Theta, delta, beta', sqrt(s) and lambda
  double constraint = 1/lambda_scaled_int * E3_int/sqrt(1 - pow(beta_abs,2.0)) * (1 + (beta_abs * cos_delta * cos_Theta)); 
  
  if (kt2 > constraint)
    return ( V_phi * qt2 / pow((qt2+md2_int),2.0) * 1/(kt2 * (kt2 + qt2 + md2_int - 2.0*qtkt)) );  // the matrix element
  else
    return 0;
}