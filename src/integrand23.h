//---------------------------------------------------------------------------------------
//provided by subversion
//---------------------------------------------------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------

/** @file 
* @brief Declarations for class integrand23.
*/

#ifndef INTEGRAND23_H
#define INTEGRAND23_H

/**
@author Oliver Fochler
*/

#include "vegas.h"

#define NDIM 2
#define NCOMP 1


// Function object derived from prototype integrand (defined in vegas.h).
// It calculates (via the overloaded ()-operator) the matrix element for 3->2 scattering.
/**
 * @brief Used for calculating (and integrating over) the matrix element for 2->3 scatterings
 *
 * This class serves as a function object and is derived from the prototype integrand (see vegas.h).
 * Via the overloaded ()-operator an object of this class can be "called" just as it was a function.
 * Utility routines are provided to set certain parameters.
 */
class integrand23 : public integrand
{
public:
  /** @brief standard constructor */
  integrand23() {};
  /** @brief standard destructor */
  ~integrand23() {};
  
  /** @brief Overloaded operator() that makes integrand23 a function object - for use with NR-vegas */
  double operator()(const double [], double) const;
  
  /** @brief Overloaded operator() that makes integrand23 a function object - for use with CUBA-vegas */
  void operator()(const int*, const double [], const int*, double []) const;
  
  /** @brief Utility routine to set the Debye mass (squared) scaled by 1/s
   * @param[in] x_ Debye mass squared in GeV^2 divided by mandelstam s (GeV^2)
   **/
  void set_md2(const double x_) {md2_int = x_;}
  
  /** @brief Utility routine to set the mean free path
  * @param[in] x_ sqrt(s) in GeV * mean free path in 1/GeV (dimensionless)
  **/
  void set_lambda(const double x_) {lambda_scaled_int = x_;}
  
  /** @brief Utility routine to set the boost velocity from lab frame to CMS
  * @param[in] x_ boost velocity beta from lab frame to CMS
  **/
  void set_beta(const double x_) {beta_int = x_;}
  
  /** @brief Utility routine to set the angle between beta and CMS-axis
  * @param[in] x_ cos(theta) where theta is the angle between the boost velocity vector and the axis of the CMS system
  **/
  void set_cos_theta(const double x_) {cos_theta_int = x_;}
                   
  
private:
  /** @brief debye mass squared scaled by 1/s (dimensionless) */
  double md2_int;
  /** @brief mean free path scaled by sqrt(s) (dimensionless) */
  double lambda_scaled_int;
  /** @brief boost velocity from lab frame to CMS of 2->3 scattering */
  double beta_int;
  /** @brief cos(theta) where theta is the angle between the boost velocity vector and the axis of the CMS system */
  double cos_theta_int;
};

#endif
