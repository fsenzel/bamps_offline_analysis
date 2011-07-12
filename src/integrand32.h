//---------------------------------------------------------------------------------------
//provided by subversion
//---------------------------------------------------------------------------------------
//$HeadURL: file:///home/fochler/svnrepos/cascade_full/main/branches/restructure_branch/src/integrand32.h $
//$LastChangedDate: 2009-07-20 11:50:59 +0200 (Mon, 20 Jul 2009) $
//$LastChangedRevision: 77 $
//$LastChangedBy: fochler $
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------


#ifndef INTEGRAND32_H
#define INTEGRAND32_H

/**
@author Oliver Fochler
*/

#include "vegas.h"

#define NDIM 2
#define NCOMP 1


// Function object derived from prototype integrand (defined in vegas.h).
// It calculates (via the overloaded ()-operator) the matrix element for 3->2 scattering.

class integrand32 : public integrand
{
public:
  integrand32() {};
  ~integrand32() {};
  
  // overloaded operator() that makes integrand32 a function object - for use with NR-vegas
  double operator()(const double [], double) const;
  
  // verloaded operator() that makes integrand32 a function object - for use with CUBA-vegas
  void operator()(const int*, const double [], const int*, double []) const;
  
  void set_E1(const double x_) {E1_int = x_;}
  void set_E3(const double x_) {E3_int = x_;}
  void set_cos_gamma(const double x_) {cos_gamma_int = x_;}
  void set_md2(const double x_) {md2_int = x_;}
  void set_sqrtS(const double x_) {sqrtS_int = x_;}
  void set_lambda(const double x_) {lambda_scaled_int = x_;}
  
  void set_beta_vec(const double[]);
                   
  
private:
  // auxiliary member variables needed for the calculation of the matrix element
  // set via the public routines integrand32::set_???(???)
  double beta_abs;
  double beta_vec[4];
  
  double E1_int;
  double E3_int;
  double cos_gamma_int;
  double md2_int;
  double sqrtS_int;
  double lambda_scaled_int;
};


#endif
