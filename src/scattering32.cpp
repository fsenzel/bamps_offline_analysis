//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL: file:///home/fochler/svnrepos/cascade_full/main/branches/restructure_branch/src/scattering32.cpp $
//$LastChangedDate: 2010-07-06 16:24:24 +0200 (Tue, 06 Jul 2010) $
//$LastChangedRevision: 115 $
//$LastChangedBy: fochler $
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------

//-----------------------------------------------------------
// Select the method for numerical integration
// Only one of the choices below must be defined!
//
#define NR_VEGAS
// #define CUBA_VEGAS
// #define CUBA_SUAVE
// #define CUBA_DIVONNE
//-----------------------------------------------------------

//-----------------------------------------------------------
// parameters for the integration via the Cuba routines
#define EPSREL 1e-2
#define EPSABS 1e-3
#define VERBOSE 0
//bits 0&1: output level, from 0 (no output) to 3 (lots of information)
//bit 2: 0 = all sets of samples collected on a subregion contribute, 1 = only the last and largest sample contributes
//bit 3: 0 = Sobol quasi-random numbers are used, 1 = Mersenne-Twister pseudo-random numbers are used

#define MINEVAL 0
#define MAXEVAL 600

// Vegas
#define NSTART 50
#define NINCREASE 100

// Suave
#define NNEW 200
#define FLATNESS 50

// Divonne
#define KEY1 -1
#define KEY2 -1
#define KEY3 0
#define MAXPASS 2
#define BORDER 0
#define MAXCHISQ 10
#define MINDEVIATION 0.25
#define NEXTRA 0
#define PEAKFINDER 0
//-----------------------------------------------------------



#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>

#include "scattering32.h"

#include "integrate.h"
#include "lorentz.h"
#include "random.h"
#include "FPT_compare.h"
#include "prefactors23.h"
//#include "distributions.h"  //only needed when a gaussian proposal function is used in scattering32::getMomenta32_metropolis(..)


using namespace std;


/**
 * This standard constructor sets all pointers to NULL and sets the standard values:
 *- #direct_estimate = true
 *- #collision = true;
 *- #total_ratio = 1.;
 */
scattering32::scattering32()
{
  P1 = NULL;
  P2 = NULL;
  P3 = NULL;
  P1cell = NULL;
  P2cell = NULL;
  P3cell = NULL;
  P1cm = NULL;
  P2cm = NULL;
  P3cm = NULL;

  direct_estimate = true;
  collision = true;
  total_ratio = 1;
}



/**
 *
 * @param vx x-component of the collective velocity of the cell
 * @param vy y-component of the collective velocity of the cell
 * @param vz z-component of the collective velocity of the cell
 * @param P1arg[] momentum vector of particle 1
 * @param P2arg[] momentum vector of particle 2
 * @param P3arg[] momentum vector of particle 3
 * @param sqrtS_arg center of mass energy, sqrt(s)
 * @param md2g_scaled_arg squared debye mass, scaled by s
 * @param lambda_scaled_arg mean free path, scaled by sqrt(s)
 */
scattering32::scattering32( const double vx, const double vy, const double vz,
                            const double P1arg[], const double P2arg[], const double P3arg[],
                            const FLAVOR_TYPE F1_arg, const FLAVOR_TYPE F2_arg, const FLAVOR_TYPE F3_arg,
                            const double sqrtS_arg, const double md2g_scaled_arg, const double lambda_scaled_arg )
    : F1( F1_arg ), F2( F2_arg ), F3( F3_arg ), sqrtS( sqrtS_arg ), md2g_scaled( md2g_scaled_arg ), lambda_scaled( lambda_scaled_arg )
{
  P1 = new double[4];
  P2 = new double[4];
  P3 = new double[4];
  P1cell = new double[4];
  P2cell = new double[4];
  P3cell = new double[4];
  P1cm = new double[4];
  P2cm = new double[4];
  P3cm = new double[4];

  cell_beta_vec[0] = 0.0;
  cell_beta_vec[1] = vx;
  cell_beta_vec[2] = vy;
  cell_beta_vec[3] = vz;

  for ( int i = 0; i <= 3; i++ )
  {
    P1[i] = P1arg[i];
    P2[i] = P2arg[i];
    P3[i] = P3arg[i];
  }

  E1_selected = 0.0;
  E3_selected = 0.0;
  cos_gamma = 0.0;
  N = 0;

  // boost momenta to rest frame of the cell
  lorentz( cell_beta_vec, P1, P1cell );
  lorentz( cell_beta_vec, P2, P2cell );
  lorentz( cell_beta_vec, P3, P3cell );

  // boost momenta to the centre of mass system of the colliding particles
  double totE = P1cell[0] + P2cell[0] + P3cell[0];
  for ( int j = 1; j <= 3; j++ )
  {
    beta_vec[j] = ( P1cell[j] + P2cell[j] + P3cell[j] ) / totE;
  }
  beta_abs = sqrt( pow( beta_vec[1], 2.0 ) + pow( beta_vec[2], 2.0 ) + pow( beta_vec[3], 2.0 ) );

  lorentz( beta_vec, P1cell, P1cm );
  lorentz( beta_vec, P2cell, P2cm );
  lorentz( beta_vec, P3cell, P3cm );


//   for(int i=0; i<7; i++)
//     for(int j=0; j<4; j++)
//       rotated_beta[i][j] = 0.0;

  // rotated beta for the different combinations of p1 and p3
  rotate_beta( P1cm, P3cm, rotated_beta[1] );
  rotate_beta( P1cm, P2cm, rotated_beta[2] );
  rotate_beta( P2cm, P3cm, rotated_beta[3] );
  rotate_beta( P2cm, P1cm, rotated_beta[4] );
  rotate_beta( P3cm, P2cm, rotated_beta[5] );
  rotate_beta( P3cm, P1cm, rotated_beta[6] );

  direct_estimate = true;
  collision = true;
  total_ratio = 1;
}


scattering32::~scattering32()
{
  delete[] P1;
  P1 = NULL;
  delete[] P2;
  P2 = NULL;
  delete[] P3;
  P3 = NULL;
  delete[] P1cell;
  P1cell = NULL;
  delete[] P2cell;
  P2cell = NULL;
  delete[] P3cell;
  P3cell = NULL;
  delete[] P1cm;
  P1cm = NULL;
  delete[] P2cm;
  P2cm = NULL;
  delete[] P3cm;
  P3cm = NULL;
}



/**
 * This method sets all necessary parameters for a given particle triplet. Previous values are deleted or overwritten.
 * Using this method an scattering32 object can be re-used for multiple particle triplets, thus reducing the need to constantly
 * creating new objects.
 *
 * Either this method or the constructor taking the same arguments MUST be called prior to any other methods of the class!
 *
 * @param vx x-component of the collective velocity of the cell
 * @param vy y-component of the collective velocity of the cell
 * @param vz z-component of the collective velocity of the cell
 * @param P1arg[] momentum vector of particle 1
 * @param P2arg[] momentum vector of particle 2
 * @param P3arg[] momentum vector of particle 3
 * @param sqrtS_arg center of mass energy, sqrt(s)
 * @param md2g_scaled_arg squared debye mass, scaled by s
 * @param lambda_scaled_arg mean free path, scaled by sqrt(s)
 * @return The absolute value of the boost velocity between cell frame and CMS (#beta_vec).
 */
double scattering32::setParameter( const double vx, const double vy, const double vz,
                                   const double P1arg[], const double P2arg[], const double P3arg[],
                                   const FLAVOR_TYPE F1_arg, const FLAVOR_TYPE F2_arg, const FLAVOR_TYPE F3_arg,
                                   const double sqrtS_arg, const double md2g_scaled_arg, const double lambda_scaled_arg )
{
  delete[] P1;
  P1 = NULL;
  delete[] P2;
  P2 = NULL;
  delete[] P3;
  P3 = NULL;
  delete[] P1cell;
  P1cell = NULL;
  delete[] P2cell;
  P2cell = NULL;
  delete[] P3cell;
  P3cell = NULL;
  delete[] P1cm;
  P1cm = NULL;
  delete[] P2cm;
  P2cm = NULL;
  delete[] P3cm;
  P3cm = NULL;

  cell_beta_vec[0] = 0.0;
  cell_beta_vec[1] = vx;
  cell_beta_vec[2] = vy;
  cell_beta_vec[3] = vz;

  P1 = new double[4];
  P2 = new double[4];
  P3 = new double[4];
  P1cell = new double[4];
  P2cell = new double[4];
  P3cell = new double[4];
  P1cm = new double[4];
  P2cm = new double[4];
  P3cm = new double[4];

  for ( int i = 0; i <= 3; i++ )
  {
    P1[i] = P1arg[i];
    P2[i] = P2arg[i];
    P3[i] = P3arg[i];
  }

  F1 = F1_arg;
  F2 = F2_arg;
  F3 = F3_arg;
  sqrtS = sqrtS_arg;
  md2g_scaled = md2g_scaled_arg;
  lambda_scaled = lambda_scaled_arg;

  E1_selected = 0.0;
  E3_selected = 0.0;
  cos_gamma = 0.0;
  N = 0;

  // boost momenta to rest frame of the cell
  lorentz( cell_beta_vec, P1, P1cell );
  lorentz( cell_beta_vec, P2, P2cell );
  lorentz( cell_beta_vec, P3, P3cell );

  // boost momenta to the centre of mass system of the colliding particles
  double totE = P1cell[0] + P2cell[0] + P3cell[0];
  for ( int j = 1; j <= 3; j++ )
  {
    beta_vec[j] = ( P1cell[j] + P2cell[j] + P3cell[j] ) / totE;
  }
  beta_abs = sqrt( pow( beta_vec[1], 2.0 ) + pow( beta_vec[2], 2.0 ) + pow( beta_vec[3], 2.0 ) );

  lorentz( beta_vec, P1cell, P1cm );
  lorentz( beta_vec, P2cell, P2cm );
  lorentz( beta_vec, P3cell, P3cm );

//   for(int i=0; i<7; i++)
//     for(int j=0; j<4; j++)
//       rotated_beta[i][j] = 0.0;

  // rotated beta for the different combinations of p1 and p3
  rotate_beta( P1cm, P3cm, rotated_beta[1] );
  rotate_beta( P1cm, P2cm, rotated_beta[2] );
  rotate_beta( P2cm, P3cm, rotated_beta[3] );
  rotate_beta( P2cm, P1cm, rotated_beta[4] );
  rotate_beta( P3cm, P2cm, rotated_beta[5] );
  rotate_beta( P3cm, P1cm, rotated_beta[6] );

  direct_estimate = true;
  collision = true;
  total_ratio = 1;

  return beta_abs;
}



/**
 * scattering32::getIntegral32(..) provides the integration of the matrix element for 3->2 processes
 * It is:
 * I32 = \int_{0}^{1} d\cos(theta) \int_{0}^{\M_PI} d\phi
 *         \frac{ q_t^2 }{ (q_t^2 + m_D^2)^2  k_t^2 [k_t^2 + q_t^2 - 2 \vec{k}_t\vec{q}_t + m_D^2 ] }        (1)
 *
 *
 * The following angles are relevant for the calculation of I32. Considered for the case, where particle 3 is absorbed by particle 2.
 * Other possible combinations are considered in getIntegral32(..) by renaming the variables accessed by the Monte-Carlo integration.
 * See notes for more details!
 * - p1,p2,p3 are the initial momentum vectors, p1' and p2' are the final momentum vectors
 * - vector p1 parallel to z-axis
 * - vectors p1,p2,p3 are located in the x-z-plane
 * - theta is the angle between p1' and p1
 * - gamma is the angle between p3 and p1
 * - phi is the angle between e_{x} and p1'*e_{x}
 * - the vector p1' can be decomposed into p1' = sqrt(s)/2 * ( sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta) )
 * - p1'*p3 = sqrt(s)/2*E3*cos(delta)
 * - delta is the angle between p1' and p3, cos(delta) = sin(gamma)*sin(theta)*cos(phi) + cos(gamma)*cos(theta)
 *
 * The LPM cutoff requires that
 *
 * kt^2 > 1/lambda * E3 * ( 1 + beta' * sqrt(s)/2 * cos(delta) * cos(Theta) )    (2)
 *
 * where Theta (capital!) is the angle between the boost-vector beta' (here: beta) and the vector p1'. For comparison with 2->3 note that
 * cosh(y) = E3 / kt,  tanh(y) = p1'*p3/E3 = sqrt(s)/2*cos(delta)  and  cos(Theta) = beta'*p1'/abs(beta'*p1').
 * Furthermore E3 > 1/lambda must be fulfilled, as E3 > kt and the gamma-factor (not the angle gamma) is larger than 1.
 * See notes for more details!
 *
 * theIntegrand is a function object of type integrand32 derived from type integrand. It can be passed to the VEGAS integration routine
 * simply via vegas(2, theIntegrand, &tgral, &sd, &chi2a).
 * Auxiliary members of integrand32 (E1_int, E3_int, cos_gamma_int, beta_vec, beta_abs, md2_int, sqrtS_int) are set
 * via integrand32::set_XX(..) prior to integration via vegas.
 * Results are stored in I32[..].
 *
 * @return Integral over the 3->2 matrix element I32.
 */
double scattering32::getIntegral32_vegas( int& initialStateIndex )
{
  theIntegrand.set_md2( md2g_scaled );
  theIntegrand.set_sqrtS( sqrtS );

  // set lambda for use in the integrand
  theIntegrand.set_lambda( lambda_scaled );

  double E1s = P1cm[0] / sqrtS;          // CMS energy of particle 1 scaled by sqrt(s)
  double E2s = P2cm[0] / sqrtS;          // CMS energy of particle 2 scaled by sqrt(s)
  double E3s = P3cm[0] / sqrtS;          // CMS energy of particle 3 scaled by sqrt(s)

  // calculate cos(gamma) for the possible combinations
  double cos_gamma_12 = ( P1cm[1] * P2cm[1] + P1cm[2] * P2cm[2] + P1cm[3] * P2cm[3] ) / ( P1cm[0] * P2cm[0] );  // cos(gamma) = p1*p2 / (abs(p1)*abs(p2))
  double cos_gamma_13 = ( P1cm[1] * P3cm[1] + P1cm[2] * P3cm[2] + P1cm[3] * P3cm[3] ) / ( P1cm[0] * P3cm[0] );  // cos(gamma) = p1*p3 / (abs(p1)*abs(p3))
  double cos_gamma_23 = ( P2cm[1] * P3cm[1] + P2cm[2] * P3cm[2] + P2cm[3] * P3cm[3] ) / ( P2cm[0] * P3cm[0] );  // cos(gamma) = p2*p3 / (abs(p2)*abs(p3))


  //--------------------------------------------------------------
  // parameters needed for calls to the integration routines
  int neval;  // actual number of integrand evaluations needed
  int fail;   // 0 = desired accuracy was reached, 1 = accuracy not reached, -1 = dimension out of range
  double intResult[NCOMP]; // result of the integration over the unit hypercube, NCOMP = #components, 1 in our case
  double error[NCOMP];     // presumed absolute error of integral
  double prob[NCOMP];      // xi^2 probability that error is NOT a reliable estimate of the true integration error
  //--------------------------------------------------------------

  //--------------------------------------------------------------
  // create the functionoid that handles the integration
  // will be called later with: integrate( theIntegrand, neval, fail, intResult, error, prob );
  // the integration routine is chosen via precompiler defined switches
#ifdef CUBA_VEGAS
  integrate_vegas integrate( NDIM, NCOMP, EPSREL, EPSABS, VERBOSE, MINEVAL, MAXEVAL, NSTART, NINCREASE );
#endif

#ifdef CUBA_SUAVE
  integrate_suave integrate( NDIM, NCOMP, EPSREL, EPSABS, VERBOSE, MINEVAL, MAXEVAL, NNEW, FLATNESS );
#endif

#ifdef CUBA_DIVONNE
  const int ngiven = 0;
  const int ldxgiven = 0;
  double * xgiven;
  int nregions = -1;
  integrate_divonne integrate( NDIM, NCOMP, EPSREL, EPSABS, VERBOSE, MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
                               BORDER, MAXCHISQ, MINDEVIATION, ngiven, ldxgiven, xgiven, NEXTRA, PEAKFINDER, nregions );
#endif

#ifdef NR_VEGAS
  integrate_nr_vegas integrate( NDIM );
#endif
  //--------------------------------------------------------------

  vector<double> I32( 7, 0 );   // vector holding the integration results for the various possible combinations, initialized with zeros
  // I32[0] is not used

  // absorption of gluon 3
  if ( F3 == gluon && E3s > 1 / lambda_scaled ) // otherwise LPM constraint (2) could not be fulfilled at all
  {
    // gluon 3 absorped by particle 2
    theIntegrand.set_E1( E1s );
    theIntegrand.set_E3( E3s );
    theIntegrand.set_cos_gamma( cos_gamma_13 );
    theIntegrand.set_beta_vec( rotated_beta[1] );

    integrate( theIntegrand, neval, fail, intResult, error, prob );
    I32[1] = intResult[0];


    // gluon 3 absorped by particle 1
    theIntegrand.set_E1( E2s );
    theIntegrand.set_cos_gamma( cos_gamma_23 );
    theIntegrand.set_beta_vec( rotated_beta[3] );

    integrate( theIntegrand, neval, fail, intResult, error, prob );
    I32[3] = intResult[0];
  }

  // absorption of particle 2
  if ( F2 == gluon && E2s > 1 / lambda_scaled )  // otherwise LPM constraint (2) could not be fulfilled at all
  {
    // gluon 2 absorped by gluon 3
    theIntegrand.set_E1( E1s );
    theIntegrand.set_E3( E2s );
    theIntegrand.set_cos_gamma( cos_gamma_12 );
    theIntegrand.set_beta_vec( rotated_beta[2] );

    integrate( theIntegrand, neval, fail, intResult, error, prob );
    I32[2] = intResult[0];

    // gluon 2 absorped by particle 1
    theIntegrand.set_E1( E3s );
    theIntegrand.set_cos_gamma( cos_gamma_23 );
    theIntegrand.set_beta_vec( rotated_beta[5] );

    integrate( theIntegrand, neval, fail, intResult, error, prob );
    I32[5] = intResult[0];
  }

  // absorption of gluon 1
  if ( F1 == gluon && E1s > 1 / lambda_scaled )  // otherwise LPM constraint (2) could not be fulfilled at all
  {
    // gluon 1 absorped by particle 3
    theIntegrand.set_E1( E2s );
    theIntegrand.set_E3( E1s );
    theIntegrand.set_cos_gamma( cos_gamma_12 );
    theIntegrand.set_beta_vec( rotated_beta[4] );

    integrate( theIntegrand, neval, fail, intResult, error, prob );
    I32[4] = intResult[0];

    //gluon 1 absorped by particle 2
    theIntegrand.set_E1( E3s );
    theIntegrand.set_cos_gamma( cos_gamma_13 );
    theIntegrand.set_beta_vec( rotated_beta[6] );

    integrate( theIntegrand, neval, fail, intResult, error, prob );
    I32[6] = intResult[0];
  }

  double sumI32 = I32[1] + I32[2] + I32[3] + I32[4] + I32[5] + I32[6];  // total results for I32

  // Randomly pick one of the possible combinations according to their weights.
  // E1_selected, E3_selected, cos_gamma and N are priate members of scattering32, the approriate values are set below
  double who = sumI32 * ran2();
  vector<double> pp( 7, 0 );
  for ( int i = 1; i <= 6; i++ )
    pp[i] = pp[i-1] + I32[i];

  unsigned int _F1;
  unsigned int _F2;

  if ( who < pp[1] )
  {
    E1_selected = E1s;
    E3_selected = E3s;
    cos_gamma = cos_gamma_13;
    N = 1;

    _F1 = std::min( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );
    _F2 = std::max( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );
  }
  else if ( who < pp[2] )
  {
    E1_selected = E1s;
    E3_selected = E2s;
    cos_gamma = cos_gamma_12;
    N = 2;

    _F1 = std::min( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F3 ) );
    _F2 = std::max( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F3 ) );
  }
  else if ( who < pp[3] )
  {
    E1_selected = E2s;
    E3_selected = E3s;
    cos_gamma = cos_gamma_23;
    N = 3;

    _F1 = std::min( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );
    _F2 = std::max( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );
  }
  else if ( who < pp[4] )
  {
    E1_selected = E2s;
    E3_selected = E1s;
    cos_gamma = cos_gamma_12;
    N = 4;

    _F1 = std::min( static_cast<unsigned int>( F2 ), static_cast<unsigned int>( F3 ) );
    _F2 = std::max( static_cast<unsigned int>( F2 ), static_cast<unsigned int>( F3 ) );
  }
  else if ( who < pp[5] )
  {
    E1_selected = E3s;
    E3_selected = E2s;
    cos_gamma = cos_gamma_23;
    N = 5;

    _F1 = std::min( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F3 ) );
    _F2 = std::max( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F3 ) );
  }
  else
  {
    E1_selected = E3s;
    E3_selected = E1s;
    cos_gamma = cos_gamma_13;
    N = 6;

    _F1 = std::min( static_cast<unsigned int>( F2 ), static_cast<unsigned int>( F3 ) );
    _F2 = std::max( static_cast<unsigned int>( F2 ), static_cast<unsigned int>( F3 ) );
  }

  double total_prefactor = 0;

  initialStateIndex = -1;

  if (( _F1 + _F2 ) == 0 ) // g+g+g -> g+g, g+g+g -> q+qbar
  {
    prefactor23_gg_ggg pObj1;
    prefactor23_gg_qqbarg pObj2;
    
    initialStateIndex = 0;

    total_prefactor = pObj1.prefactor() * pObj1.symFactor32() + pObj2.prefactor() * pObj2.symFactor32();
  }
  else if ( _F1 == _F2 )  // q+q+g -> q+q, qbar+qbar+g -> qbar+qbar
  {
    prefactor23_qq_qqg pObj1;
    
    initialStateIndex = 1;

    total_prefactor = pObj1.prefactor() * pObj1.symFactor32();
  }
  else if (( _F1 * _F2 ) == 0 ) // g+q+g -> g+q, g+qbar+g -> g+qbar
  {
    prefactor23_qg_qgg pObj1;
    
    initialStateIndex = 2;

    total_prefactor = pObj1.prefactor() * pObj1.symFactor32();
  }
  else if (( _F2 - _F1 ) == 1 && ( _F2 % 2 ) == 0 )  // q+qbar+g -> q+qbar, q+qbar+g -> q'+qbar', q+qbar+g -> g+g
  {
    prefactor23_qqbar_qqbarg pObj1;
    prefactor23_qqbar_qqbarDashg pObj2;
    prefactor23_qqbar_ggg pObj3;
    
    initialStateIndex = 3;

    total_prefactor = pObj1.prefactor() * pObj1.symFactor32() + pObj2.prefactor() * pObj2.symFactor32() + pObj3.prefactor() * pObj3.symFactor32();
  }
  else // q+q'+g -> q+q', q+qbar'+g -> q+qbar'
  {
    prefactor23_qqdash_qqdashg pObj1;
    
    initialStateIndex = 4;

    total_prefactor = pObj1.prefactor() * pObj1.symFactor32();
  }

  return sumI32 * total_prefactor;
}




/**
 * The integral over the matrix element, I32, is estimated using an upper limit of the real matrix element, which can be
 * integrated analytically. The theta-function is dealt with by finding areas of integration such that the resulting integral
 * is larger than the integral with the theta-function. See notes.
 *
 * Using this method should be faster, though likely less accurate, than employing direct numerical integration
 * via #getIntegral32_vegas.
 *
 * Thus the estimated I32 is larger than the real I32. Using a number of sampled points the method #get_I32estimate_ratio returns
 * an estimate of the ratio of the real to the estimated I32. Either
 * - I32 is directly corrected by this ratio (when #direct_estimate = true)
 *
 * OR
 *
 * - the estimated I32 is used (when #direct_estimate = false).
 *
 * In the latter case the method #getMomenta32_estimate_rejection introduces a further decision whether a collision really
 * takes place.
 *
 * @return An estimate for I32, with one of the following characteristics
 * - corrected by the estimated ratio to the real result, when #direct_estimate = true
 * - not corrected by the estimated ratio to the real result, when #direct_estimate = false
 */
double scattering32::getIntegral32_fast( int& initialStateIndex )
{
  const double lower_limit = 1.0e-4;

  double E1s = P1cm[0] / sqrtS;          // CMS energy of particle 1 scaled by sqrt(s)
  double E2s = P2cm[0] / sqrtS;          // CMS energy of particle 2 scaled by sqrt(s)
  double E3s = P3cm[0] / sqrtS;          // CMS energy of particle 3 scaled by sqrt(s)

  // calculate cos(gamma) for the possible combinations
  double cos_gamma_12 = ( P1cm[1] * P2cm[1] + P1cm[2] * P2cm[2] + P1cm[3] * P2cm[3] ) / ( P1cm[0] * P2cm[0] );  // cos(gamma) = p1*p2 / (abs(p1)*abs(p2))
  double cos_gamma_13 = ( P1cm[1] * P3cm[1] + P1cm[2] * P3cm[2] + P1cm[3] * P3cm[3] ) / ( P1cm[0] * P3cm[0] );  // cos(gamma) = p1*p3 / (abs(p1)*abs(p3))
  double cos_gamma_23 = ( P2cm[1] * P3cm[1] + P2cm[2] * P3cm[2] + P2cm[3] * P3cm[3] ) / ( P2cm[0] * P3cm[0] );  // cos(gamma) = p2*p3 / (abs(p2)*abs(p3))

  vector<double> I32( 7, 0 );   // vector holding the integration results for the various possible combinations, initialized with zeros
  // I32[0] is not used
  vector<double> ratio( 7, 0 ); // vector holding the ratio of the estimated result to the real result, initialized with zeros
  // ratio[0] is not used

  // absorption of gluon 3
  if ( F3 == gluon && E3s > 1 / lambda_scaled ) // otherwise LPM constraint (2) could not be fulfilled at all
  {
    // gluon 3 absorped by particle 2
    I32[1] = get_I32estimate( E1s, E3s, cos_gamma_13 );
    if ( I32[1] > lower_limit )
      ratio[1] = get_I32estimate_ratio( E1s, E3s, cos_gamma_13, 1 );

    // gluon 3 absorped by particle 1
    I32[3] = get_I32estimate( E2s, E3s, cos_gamma_23 );
    if ( I32[3] > lower_limit )
      ratio[3] = get_I32estimate_ratio( E2s, E3s, cos_gamma_23, 3 );
  }

  // absorption of gluon 2
  if ( F2 == gluon && E2s > 1 / lambda_scaled )  // otherwise LPM constraint (2) could not be fulfilled at all
  {
    // gluon 2 absorped by particle 3
    I32[2] = get_I32estimate( E1s, E2s, cos_gamma_12 );
    if ( I32[2] > lower_limit )
      ratio[2] = get_I32estimate_ratio( E1s, E2s, cos_gamma_12, 2 );

    // gluon 2 absorped by particle 1
    I32[5] = get_I32estimate( E3s, E2s, cos_gamma_23 );
    if ( I32[5] > lower_limit )
      ratio[5] = get_I32estimate_ratio( E3s, E2s, cos_gamma_23, 5 );
  }

  // absorption of gluon 1
  if ( F1 == gluon && E1s > 1 / lambda_scaled )  // otherwise LPM constraint (2) could not be fulfilled at all
  {
    // gluon 1 absorped by particle 3
    I32[4] = get_I32estimate( E2s, E1s, cos_gamma_12 );
    if ( I32[4] > lower_limit )
      ratio[4] = get_I32estimate_ratio( E2s, E1s, cos_gamma_12, 4 );

    //gluon 1 absorped by particle 2
    I32[6] = get_I32estimate( E3s, E1s, cos_gamma_13 );
    if ( I32[6] > lower_limit )
      ratio[6] = get_I32estimate_ratio( E3s, E1s, cos_gamma_13, 6 );
  }


  double sumI32, total_ratio;
  vector<double> pp( 7, 0 );

  if ( direct_estimate )
  {
    sumI32 = I32[1] * ratio[1] + I32[2] * ratio[2] + I32[3] * ratio[3] + I32[4] * ratio[4] + I32[5] * ratio[5] + I32[6] * ratio[6];  // total result for I32

    for ( int i = 1; i <= 6; i++ )
      pp[i] = pp[i-1] + I32[i] * ratio[i];
  }
  else
  {
    sumI32 = I32[1] + I32[2] + I32[3] + I32[4] + I32[5] + I32[6];  // total result for I32
    total_ratio = ( I32[1] * ratio[1] + I32[2] * ratio[2] + I32[3] * ratio[3] + I32[4] * ratio[4] + I32[5] * ratio[5] + I32[6] * ratio[6] ) / sumI32;

    for ( int i = 1; i <= 6; i++ )
      pp[i] = pp[i-1] + I32[i];
  }


  double total_prefactor = 0;
  if ( sumI32 > 0 )
  {

    // Randomly pick one of the possible combinations according to their weights.
    // E1_selected, E3_selected, cos_gamma and N are priate members of scattering32, the approriate values are set below
    double who = sumI32 * ran2();


    unsigned int _F1;
    unsigned int _F2;

    if ( who < pp[1] )
    {
      E1_selected = E1s;
      E3_selected = E3s;
      cos_gamma = cos_gamma_13;
      N = 1;

      _F1 = std::min( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );
      _F2 = std::max( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );
    }
    else if ( who < pp[2] )
    {
      E1_selected = E1s;
      E3_selected = E2s;
      cos_gamma = cos_gamma_12;
      N = 2;

      _F1 = std::min( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F3 ) );
      _F2 = std::max( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F3 ) );
    }
    else if ( who < pp[3] )
    {
      E1_selected = E2s;
      E3_selected = E3s;
      cos_gamma = cos_gamma_23;
      N = 3;

      _F1 = std::min( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );
      _F2 = std::max( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );
    }
    else if ( who < pp[4] )
    {
      E1_selected = E2s;
      E3_selected = E1s;
      cos_gamma = cos_gamma_12;
      N = 4;

      _F1 = std::min( static_cast<unsigned int>( F2 ), static_cast<unsigned int>( F3 ) );
      _F2 = std::max( static_cast<unsigned int>( F2 ), static_cast<unsigned int>( F3 ) );
    }
    else if ( who < pp[5] )
    {
      E1_selected = E3s;
      E3_selected = E2s;
      cos_gamma = cos_gamma_23;
      N = 5;

      _F1 = std::min( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F3 ) );
      _F2 = std::max( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F3 ) );
    }
    else
    {
      E1_selected = E3s;
      E3_selected = E1s;
      cos_gamma = cos_gamma_13;
      N = 6;

      _F1 = std::min( static_cast<unsigned int>( F2 ), static_cast<unsigned int>( F3 ) );
      _F2 = std::max( static_cast<unsigned int>( F2 ), static_cast<unsigned int>( F3 ) );
    }


    if (( _F1 + _F2 ) == 0 ) // g+g+g -> g+g, g+g+g -> q+qbar
    {
      prefactor23_gg_ggg pObj1;
      prefactor23_gg_qqbarg pObj2;
      
      initialStateIndex = 0;

      total_prefactor = pObj1.prefactor() * pObj1.symFactor32() + pObj2.prefactor() * pObj2.symFactor32();
    }
    else if ( _F1 == _F2 )  // q+q+g -> q+q, qbar+qbar+g -> qbar+qbar
    {
      prefactor23_qq_qqg pObj1;
      
      initialStateIndex = 1;

      total_prefactor = pObj1.prefactor() * pObj1.symFactor32();
    }
    else if (( _F1 * _F2 ) == 0 ) // g+q+g -> g+q, g+qbar+g -> g+qbar
    {
      prefactor23_qg_qgg pObj1;
      
      initialStateIndex = 2;

      total_prefactor = pObj1.prefactor() * pObj1.symFactor32();
    }
    else if (( _F2 - _F1 ) == 1 && ( _F2 % 2 ) == 0 )  // q+qbar+g -> q+qbar, q+qbar+g -> q'+qbar', q+qbar+g -> g+g
    {
      prefactor23_qqbar_qqbarg pObj1;
      prefactor23_qqbar_qqbarDashg pObj2;
      prefactor23_qqbar_ggg pObj3;

      initialStateIndex = 3;
      
      total_prefactor = pObj1.prefactor() * pObj1.symFactor32() + pObj2.prefactor() * pObj2.symFactor32() + pObj3.prefactor() * pObj3.symFactor32();
    }
    else // q+q'+g -> q+q', q+qbar'+g -> q+qbar', qbar+qbar'+g-> qbar+qbar'
    {
      prefactor23_qqdash_qqdashg pObj1;
      
      initialStateIndex = 4;

      total_prefactor = pObj1.prefactor() * pObj1.symFactor32();
    }
  }
  else
  {
    initialStateIndex = -1;
    E1_selected = E1s;
    E3_selected = E3s;
    cos_gamma = cos_gamma_13;
    N = 1;
  }
  return sumI32 * total_prefactor;
}



/**
 * For a given choice (p1,p3) of particles 1 and 3 this routine calculates an upper estimate of the integral I32.
 *
 * @param E1 Energy of particle 1
 * @param E3 Energy of particle 3
 * @param cosgamma_arg cos(gamma) for the given choice (p1,p3)
 * @return Estimated I32 for given choice of particle 1 and 3
 */
double scattering32::get_I32estimate( const double E1, const double E3, const double cosgamma_arg ) const
{
  double F;

  double cosgamma = fabs( cosgamma_arg );
  double singamma = sqrt( 1 - pow( cosgamma, 2.0 ) );
  double C = 1 / ( pow( lambda_scaled, 2 ) * pow( E3, 2 ) * ( 1 + beta_abs ) );

  double G1 = sqrt( 1 - C ) / cosgamma;
  double G3 = singamma;
  double G4 = sqrt( C ) * singamma + sqrt( 1 - C ) * cosgamma;

  double a = sqrt( 1 + md2g_scaled / pow( E1, 2.0 ) );
  double Q = pow( lambda_scaled, 2.0 ) * ( 1 + beta_abs ) / ( pow( E1, 2.0 ) * md2g_scaled );

  F = H( a, G3 ) - H( a, 0 );

  if ( G1 > G3 )
    F += H( a, min( 1.0, G1 ) ) - H( a, G3 );

  if ( 1 > G1 && G4 > max( G1, G3 ) )
    F += H( a, G4 ) - H( a, max( G1, G3 ) );

  return M_PI * Q * F;
}


/**
 * Utility method used by # get_I32estimate. It calculates the antiderivative of the estimated matrix element at given u.
 *
 * @param a
 * @param u
 * @return Value of the antiderivative of the estimated matrix element at given u.
 */
double scattering32::H( const double a, const double u ) const
{
  double a2 = pow( a, 2.0 );
  double h = ( a2 + 1 ) / ( 4 * pow( a, 3.0 ) ) * log(( a + u ) / ( a - u ) ) - ( a2 - 1 ) * u / ( 2 * a2 * ( a2 - pow( u, 2.0 ) ) ) ;
  return h;
}


/**
 * Estimates the ratio of the real result for I32 to the estimated result for I32. This is done by sampling a certain number of points
 * according to the estimated matrix element (which is always larger than the real) and taking the ratio of the function values
 * at these points. An average over the so obtained values for the ratio should give a somewhat reliable result.
 *
 * @param[in] E1 Energy of particle 1
 * @param[in] E3 Energy of particle 3
 * @param[in] cosgamma_arg cos(gamma) for the given choice (p1,p3)
 * @param[in] N_arg Number indicating the given choice (p1,p3)
 * @return Estimated ratio between real I32 and estimated I32
 */
double scattering32::get_I32estimate_ratio( const double E1, const double E3, const double cosgamma_arg, const int N_arg )
{
  const int samples = 60;
  const double A = pow( lambda_scaled, 2.0 ) * ( 1 + beta_abs ) / md2g_scaled;

  double max;
  double qt2;
  double u, phi , us;
  double g_estimate, g_real, g_compare;
  double ratio = 0;

  // the ranges in which the variables u and phi need to be sampled
  const double u_min = 0.0;
  const double u_max = 1.0;
  const double phi_min = 0.0;
  const double phi_max = M_PI;

  const double cosgamma = fabs( cosgamma_arg );
  const double singamma = sqrt( 1 - pow( cosgamma, 2.0 ) );
  const double C = 1 / ( pow( lambda_scaled, 2 ) * pow( E3, 2 ) * ( 1 + beta_abs ) );

  double sin_delta_2_limit;

  if ( md2g_scaled < 0.25 )
    max = A / ( 4 * md2g_scaled );
  else
    max = A * 0.25 / pow( 0.25 + md2g_scaled, 2.0 );

  for ( int i = 0; i < samples; i++ )
  {
    do
    {
      u = u_min + ran2() * ( u_max - u_min );
      phi = phi_min + ran2() * ( phi_max - phi_min );
      us = sqrt( 1.0 - pow( u, 2.0 ) );

      if ( u < singamma )
        sin_delta_2_limit = 1;
      else
        sin_delta_2_limit = 1 - pow(( u * cosgamma - us * singamma ), 2.0 );

      if ( sin_delta_2_limit < C )
        g_estimate = 0;
      else
      {
        qt2 = pow( E1, 2.0 ) * pow( us, 2.0 );
        g_estimate = A * qt2 / pow( qt2 + md2g_scaled, 2.0 );
      }

      g_compare = max * ran2();

    }
    while ( g_estimate < g_compare );

    cos_gamma = cosgamma_arg;
    E1_selected = E1;
    E3_selected = E3;
    N = N_arg;
    g_real = getMatrixElement( u, phi );

    if ( g_estimate > 0.000001 )
    {    
      ratio += g_real / g_estimate;
    }
  }

  ratio = ratio / samples;

  cos_gamma = 0;
  E1_selected = 0;
  E3_selected = 0;
  N = 0;
  return ratio;
}




/**
 * Samples new momenta of the outgoing particles according to the matrix element.
 *
 * Wrapper, calls a certain implementation of the sampling in u and phi.
 *
 * @param[out] u u=cos(theta), theta = angle between p1' and p1
 * @param[out] phi azimuthal angle between e_x and (p1' e_x)
 * @param[in] picked_ratio
 * @return Count how many time the comparison function was below the actual function. Only needed for rejection method. Should be zero.
 */
int scattering32::getMomenta32( double& u, double& phi, const double picked_ratio, int& typ, FLAVOR_TYPE& F1arg, FLAVOR_TYPE& F2arg )
{
  int error = 0;  // see return value

  if ( direct_estimate )
  {
    error = getMomenta32_metropolis( u, phi ); // fast

    // slow (as the comparison function in use might be off by some orders of magnitude, leading to high rejection probabilities), use is discouraged
    // error = getMomenta32_rejection(u,phi);
  }
  else
  {
    error = getMomenta32_estimate_rejection( u, phi, picked_ratio );
  }


  switch ( N )
  {
  case 1:
    F1arg = F1;
    F2arg = F2;
    break;
  case 2:
    F1arg = F1;
    F2arg = F3;
    break;
  case 3:
    F1arg = F1;
    F2arg = F2;
    break;
  case 4:
    F1arg = F2;
    F2arg = F3;
    break;
  case 5:
    F1arg = F1;
    F2arg = F3;
    break;
  case 6:
    F1arg = F2;
    F2arg = F3;
    break;
  default:
    F1arg = F1;
    F2arg = F2;
    cout << "error in scattering32::getMomenta23(...), switch(N) did resolve to default" << endl;
    break;
  }

  // sort F1 and F2 such that comparisons below are easier
  // convert to unsigned int to allow for some math that speeds up the decisions
  unsigned int _F1 = std::min( static_cast<unsigned int>( F1arg ), static_cast<unsigned int>( F2arg ) );
  unsigned int _F2 = std::max( static_cast<unsigned int>( F1arg ), static_cast<unsigned int>( F2arg ) );

  if (( _F1 + _F2 ) == 0 ) // g+g -> g+g+g, g+g -> q+qbar+g
  {
    prefactor23_gg_ggg pObj1;
    prefactor23_gg_qqbarg pObj2;
    double pf1 = pObj1.prefactor() * pObj1.symFactor32();
    double pf2 = pObj2.prefactor() * pObj2.symFactor32();

    if ( ran2() < pf1 / ( pf1 + pf2 ) )
    {
      typ = 321; // g+g+g -> g+g
    }
    else
    {
      typ = 322; // g+g+g -> q+qbar
      sampleFlavor23( F1arg, F2arg );
    }
  }
  else if ( _F1 == _F2 )  // q+q+g -> q+q, qbar+qbar+g -> qbar+qbar
  {
    typ = 327;
  }
  else if (( _F1 * _F2 ) == 0 ) // g+q+g -> g+q, g+qbar+g -> g+qbar
  {
    typ = 323;
  }
  else if (( _F2 - _F1 ) == 1 && ( _F2 % 2 ) == 0 )  // q+qbar+g -> q+qbar, q+qbar+g -> q'+qbar', q+qbar+g -> g+g
  {
    prefactor23_qqbar_qqbarg pObj1;
    prefactor23_qqbar_qqbarDashg pObj2;
    prefactor23_qqbar_ggg pObj3;

    double pf1 = pObj1.prefactor() * pObj1.symFactor32();
    double pf2 = pObj2.prefactor() * pObj2.symFactor32();
    double pf3 = pObj3.prefactor() * pObj3.symFactor32();

    double select = ran2() * ( pf1 + pf2 + pf3 );

    if ( select < pf1 )  // q+qbar+g -> q+qbar
    {
      typ = 324;
    }
    else if ( select < ( pf1 + pf2 ) )  // q+qbar+g -> q'+qbar'
    {
      typ = 325;
      sampleFlavor23( F1arg, F2arg, F1 );  // sample flavor excluding F1 (or anti-F1 = F2)
    }
    else  // q+qbar+g -> g+g
    {
      typ = 326;
      F1arg = gluon;
      F2arg = gluon;
    }
  }
  else // q+q'+g -> q+q', q+qbar'+g -> q+qbar'
  {
    typ = 328;
  }

  return error;  // only necessary when using the rejection method, counts how many times the comparison function is below the actual function
}



/**
 * Samples of new momenta according to the matrix element using the rejection method. Slow, since comparison function might be significantly
 * larger than the sampled function at some points. Use #getMomenta32_metropolis instead.
 *
 * @param[out] u_arg u=cos(theta), theta = angle between p1' and p1
 * @param[out] phi_arg azimuthal angle between e_x and (p1' e_x)
 * @return Count how many time the comparison function was below the actual function. Only needed for rejection method. Should be zero.
 */
int scattering32::getMomenta32_rejection( double& u_arg, double& phi_arg ) const
{
  // E1_selected, E3_selected, cos_gamma and N are private members of scattering32 and have been set in scattering32::getIntegral32_vegas()

  int error = 0;  // see return value

  double u, phi;
  double g = 0, gr = 0;

  // the ranges in which the variables u and phi need to be sampled
  const double u_min = 0.0;
  const double u_max = 1.0;
  const double phi_min = 0.0;
  const double phi_max = M_PI;

  //double beta_abs = 0;   // for testing consistency with old implementation (beta_abs = 0 needs then to be set in integrand32.cpp also)
  double max = lambda_scaled / E3_selected * 1 / ( 4 * pow( md2g_scaled, 4.0 ) ) * sqrt( 1 - pow( beta_abs, 2.0 ) ) / ( 1 - beta_abs );

  do
  {
    u = u_min + ran2() * ( u_max - u_min );
    phi = phi_min + ran2() * ( phi_max - phi_min );

    g = getMatrixElement( u, phi );

    if ( g > max )
    {
      //cout << "error in get32    g = " << g << "   max = " << max <<  endl;
      error++;
    }

    gr = max * ran2();

  }
  while ( g < gr );

  u_arg = u;
  phi_arg = phi;

  return error; // only necessary when using the rejection method, counts how many times the comparison function is below the actual function
}



/**
 * Routine for sampling u = cos(theta) and phi according to the matrix element (see scattering32::getMatrixElement(..)).
 * An implementation of the metropolis algorithm is used for sampling the distribution. This algorithm does not depend on
 * the knowledge of an absolute normalisation of the distribution
 *
 * Furthermore one needs no comparison function which makes this routine a lot faster than the version implemented in
 * scattering32::getMomenta32(..) as the comparison function in use there can be off the actual distribution by several orders
 * of magnitude, leading to high rejection probabilites for the sampled points.
 *
 * The computational expense is fixed by number of steps in the Markov chain (#n_steps).
 *
 * @param[out] u_arg u=cos(theta), theta = angle between p1' and p1
 * @param[out] phi_arg azimuthal angle between e_x and (p1' e_x)
 * @return By definition 0, provided for compability with #getMomenta32_rejection
 */
int scattering32::getMomenta32_metropolis( double& u_arg, double& phi_arg ) const
{
  double u, phi;
  double u_new, phi_new;
  double g = 0, g_new = 0;
  double ratio;

  // the ranges in which the variables u and phi need to be sampled
  const double u_min = 0.0;
  const double u_max = 1.0;
  const double phi_min = 0.0;
  const double phi_max = M_PI;

  // only needed when using gaussian distributions for proposing new steps in the markov chain
  // width
  //double sigma_u = (u_max - u_min) / 3.0;
  //double sigma_phi = (phi_max - phi_min) / 3.0;

  // randomly select initial values of u and phi, such that
  do
  {
    u = u_min + ran2() * ( u_max - u_min );
    phi = phi_min + ran2() * ( phi_max - phi_min );

    g = getMatrixElement( u, phi );
  }
  while ( FPT_COMP_E( g, 0.0 ) );


  // number of steps in the Markov chain
  const int n_steps = 50;

  // do n_steps steps
  // the current location is (u,phi)
  // one steps consists of
  //  - proposing a new point (u',phi') according to a certain proposal distribution
  //  - calculate the matrix element g(u',phi') at this new point
  //  - if g(u',phi') > g(u,phi) accept the new point (u',phi')
  //  - else accept the new point (u',phi') with a probability g(u',phi') / g(u, phi)
  //
  for ( int i = 0; i < n_steps; i++ )
  {
    do
    {
      //u_new = distributions::gaussian(u,sigma_u);    // propose new u using a gaussian with width sigma_u around the current value of u
      u_new = u_min + ran2() * ( u_max - u_min );      // propose new u using a uniform distribution over the entire range

    }
    while ( u_new < u_min || u_new > u_max );

    do
    {
      //phi_new = distributions::gaussian(phi,sigma_phi);   // propose new phi using a gaussian with width sigma_phi around the current value of phi
      phi_new = phi_min + ran2() * ( phi_max - phi_min );   // propose new phi using a uniform distribution over the entire range
    }
    while ( phi_new < phi_min || phi_new > phi_max );     // check that the new values are in range

    g_new = getMatrixElement( u_new, phi_new );           // calculate the matrix element at the proposed point

    ratio = g_new / g;                                    // ratio of g(u',phi') / g(u,phi)

    if ( FPT_COMP_GE( ratio, 1.0 ) || ran2() < ratio )    // accept if g(u',phi') > g(u,phi) or with probability g(u',phi') / g(u,phi)
    {
      u = u_new;
      phi = phi_new;
      g = g_new;
    }

  }

  u_arg = u;
  phi_arg = phi;

  return 0;  // only necessary when using the rejection method, always 0 when using Metropolis
}



/**
 * Sampling new momenta according to the upper estimate of the matrix element. Only used when #direct_estimate is set to false!
 * An additional decision whether the collision really takes place needs to be made. For this the ratio of the estimated matrix element
 * to the real matrix element at the sampled point is compared to the ratio of the chosen random number to the probability of the collision
 * (the latter computed with the estimated I32).
 *
 * @param[out] u_arg u=cos(theta), theta = angle between p1' and p1
 * @param[out] phi_arg azimuthal angle between e_x and (p1' e_x)
 * @param[in] picked_ratio Ratio of the chosen random number to the probability p32 (obtained from the estimated I32!).
 * @return Count how many time the comparison function was below the actual function. Only needed for rejection method. Should be zero.
 */
int scattering32::getMomenta32_estimate_rejection( double& u_arg, double& phi_arg, const double picked_ratio )
{
  const double A = pow( lambda_scaled, 2.0 ) * ( 1 + beta_abs ) / md2g_scaled;

  double max;
  double qt2;
  double g_estimate = 0, g_compare = 0, g_real = 0;
  double u, phi , us;

  int error = 0;  // see return value

  // the ranges in which the variables u and phi need to be sampled
  const double u_min = 0.0;
  const double u_max = 1.0;
  const double phi_min = 0.0;
  const double phi_max = M_PI;

  double sin_gamma = sqrt( 1.0 - pow( cos_gamma, 2.0 ) );
  const double C = 1 / ( pow( lambda_scaled, 2 ) * pow( E3_selected, 2 ) * ( 1 + beta_abs ) );

  double sin_delta_2_limit;

  if ( md2g_scaled < 0.25 )
    max = A / ( 4 * md2g_scaled );
  else
    max = A * 0.25 / pow( 0.25 + md2g_scaled, 2.0 );

  do
  {
    u = u_min + ran2() * ( u_max - u_min );
    phi = phi_min + ran2() * ( phi_max - phi_min );
    us = sqrt( 1.0 - pow( u, 2.0 ) );

    if ( u < sin_gamma )
      sin_delta_2_limit = 1;
    else
      sin_delta_2_limit = 1 - pow(( u * cos_gamma - us * sin_gamma ), 2.0 );

    if ( sin_delta_2_limit < C )
      g_estimate = 0;
    else
    {
      qt2 = pow( E1_selected, 2.0 ) * pow( us, 2.0 );
      g_estimate = A * qt2 / pow( qt2 + md2g_scaled, 2.0 );
    }

    if ( g_estimate > max )
    {
      error++;
    }

    g_compare = max * ran2();

  }
  while ( g_estimate < g_compare );

  g_real = getMatrixElement( u, phi );

  if (( picked_ratio*g_estimate ) < g_real )
    collision = true;
  else
    collision = false;

  u_arg = u;
  phi_arg = phi;

  return error; // only necessary when using the rejection method, counts how many times the comparison function is below the actual function
}





/**
 * Calculates the Matrix elment for 3->2 scatterings at given u and phi. The theta function stemming from the LPM-cutoff is taken into
 * account.
 *
 * @param[in] u u=cos(theta), theta = angle between p1' and p1
 * @param[in] phi azimuthal angle between e_x and (p1' e_x)
 * @return Matrix element for 3->2
 */
double scattering32::getMatrixElement( const double u, const double phi ) const
{
  // E1_selected, E3_selected, cos_gamma and N are private members of scattering32 and have been set in scattering32::getIntegral32_vegas()
  double qt2, kt2, qtkt;
  double us, v, constraint;
  double cos_delta, sin_delta_2;
  double cos_Theta;
  double g;

  //double beta_abs = 0;   // for testing consistency with old implementation (beta_abs = 0 needs then to be set in integrand32.cpp as well)
  double sin_gamma = sqrt( 1.0 - pow( cos_gamma, 2.0 ) );

  us = sqrt( 1.0 - pow( u, 2.0 ) );
  v = cos( phi );

  cos_delta = sin_gamma * us * v + cos_gamma * u;     // cos(delta) = sin(gamma)*sin(theta)*cos(phi) + cos(gamma)*cos(theta)
  sin_delta_2 = 1.0 - pow( cos_delta, 2.0 );          // sin(delta)^2 = 1 - cos(delta)^2

  qt2 = pow( E1_selected, 2.0 ) * pow( us, 2.0 );            // qt^2 = E1^2 * sin(theta)^2
  kt2 = pow( E3_selected, 2.0 ) * sin_delta_2;               // kt^2 = E3^2 * sin(delta)^2
  //vector product qt*kt = -E1*E3*sin(theta)*( cos(gamma)*sin(theta) - sin(gamma)*cos(theta)*cos(phi) )
  qtkt = -E1_selected * E3_selected * us * ( cos_gamma * us - sin_gamma * u * v );

  // Theta (capital!) is the angle between the boost-vector beta' (here: beta_vec) and the vector p1'
  if ( FPT_COMP_E( beta_abs, 0.0 ) )
    cos_Theta = 0;
  else
    cos_Theta = 1 / beta_abs * ( rotated_beta[N][1] * us * v + rotated_beta[N][2] * us * sqrt( 1 - pow( v, 2.0 ) ) + rotated_beta[N][3] * u );

  // the constraint for kt^2 depending on E3, Theta, delta, beta', sqrt(s) and lambda
  constraint = 1 / lambda_scaled * E3_selected / sqrt( 1 - pow( beta_abs, 2.0 ) ) * ( 1 + ( beta_abs * cos_delta * cos_Theta ) );

  if ( kt2 > constraint )
    g = qt2 / pow(( qt2 + md2g_scaled ), 2.0 ) * 1 / ( kt2 * ( kt2 + qt2 + md2g_scaled - 2.0 * qtkt ) );  // the matrix element
  else
    g = 0;

  return g;
}



/**
 * Sets new momenta for the outgoing particle according to u and phi. The results are written to the vectors P1_arg and P2_arg.
 *
 * @param[out] P1_arg[] Momentum vector of outgoing particle 1
 * @param[out] P2_arg[] Momentum vector of outgoing particle 2
 * @param[in] u u=cos(theta), theta = angle between p1' and p1
 * @param[in] phi azimuthal angle between e_x and (p1' e_x)
 */
void scattering32::setNewMomenta32( double P1_arg[], double P2_arg[], const double u, const double phi ) const
{
  double *temp1, *temp3;
  double P1_selected[4], P3_selected[4];

  switch ( N )
  {
  case 1:
    temp1 = P1cm;
    temp3 = P3cm;
    break;
  case 2:
    temp1 = P1cm;
    temp3 = P2cm;
    break;
  case 3:
    temp1 = P2cm;
    temp3 = P3cm;
    break;
  case 4:
    temp1 = P2cm;
    temp3 = P1cm;
    break;
  case 5:
    temp1 = P3cm;
    temp3 = P2cm;
    break;
  case 6:
    temp1 = P3cm;
    temp3 = P1cm;
    break;
  default:
    temp1 = P1cm;
    temp3 = P3cm;
    cout << "error in scattering32::setNewMomenta32(...), switch(N) did resolve to default" << endl;
    break;
  }

  for ( int i = 0; i < 4; i++ )
  {
    P1_selected[i] = temp1[i];
    P3_selected[i] = temp3[i];
  }

  double ex[4], ey[4], ez[4];
  double p1p3 = 0.0;

  // unit vector ez[] in z-direction
  // z-direction is given by p1cm[]
  for ( int i = 1; i <= 3; i++ )
  {
    ez[i] = P1_selected[i] / P1_selected[0];
    p1p3 += P1_selected[i] * P3_selected[i];  // scalar product between p1cm[] and p3cm[], needed for construction of unit vector in x-direction
  }

  // unit vector ex[] in x-direction
  // x-direction roughly points in the direction of p3cm[] (see notes)
  // construction is via: ex[] = p3cm[] - ( (p1cm[]*p3cm[]) / (p1cm*p1cm) ) * p1cm[]
  //
  p1p3 = p1p3 / pow( P1_selected[0], 2 );

  double norm = 0;
  for ( int j = 1; j <= 3; j++ )
  {
    ex[j] = P3_selected[j] - p1p3 * P1_selected[j];
    norm += ex[j] * ex[j];
  }
  // normalize to unit vector
  for ( int k = 1;k <= 3;k++ )
    ex[k] = ex[k] / sqrt( norm );

  // unit vector ey[] in y-direction
  // ey[] is contructed via vector product ez[] x ex[]
  ey[1] = ez[2] * ex[3] - ez[3] * ex[2];
  ey[2] = ex[1] * ez[3] - ex[3] * ez[1];
  ey[3] = ez[1] * ex[2] - ez[2] * ex[1];

  // u = cos(theta)
  double us;
  if ( u > 1.0 )
    us = 0.0;
  else
    us = sqrt( 1.0 - pow( u, 2.0 ) );  // us = sin(theta)

  // set new momentum vectors according to:
  //         p1'cm[] = sqrt(s)/2 * ( sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)
  //         p3'cm[] = - p1'cm[]
  // (see notes).
  for ( int i = 1; i <= 3; i++ )
  {
    P1_selected[i] = sqrtS / 2.0 * ( us * cos( phi ) * ex[i] + us * sin( phi ) * ey[i] + u * ez[i] );
    P3_selected[i] = -P1_selected[i];
  }
  P1_selected[0] = sqrtS / 2.0;
  P3_selected[0] = sqrtS / 2.0;

  // boost new momentum vectors back to original frame
  double inv_beta_cell[4];
  for ( int j = 1;j <= 3;j++ )
    inv_beta_cell[j] = -cell_beta_vec[j];

  double inv_beta[4];
  for ( int j = 1;j <= 3;j++ )
    inv_beta[j] = -beta_vec[j];

  lorentz( inv_beta, P1_selected, P1cell );
  lorentz( inv_beta, P3_selected, P2cell );

  lorentz( inv_beta_cell, P1cell, P1_arg );
  lorentz( inv_beta_cell, P2cell, P2_arg );
}



/**
 * Rotates the boost vector #beta_vec to a reference frame given by P1_arg and P3_arg.
 *
 * @param[in] P1_arg[]
 * @param[in] P3_arg[]
 * @param[out] beta_new[] The rotated vector is returned as the argument beta_new[]
 */
void scattering32::rotate_beta( const double P1_arg[], const double P3_arg[], double beta_new[] )
{
  //------------------
  // unit vectors for the reference frame used in 3->2 routines are constructed
  //
  double ex[4], ey[4], ez[4];

  double p1p3 = 0.0;

  // unit vector ez[] in z-direction
  // z-direction is given by p1_arg[]
  for ( int i = 1; i <= 3; i++ )
  {
    ez[i] = P1_arg[i] / P1_arg[0];
    p1p3 += P1_arg[i] * P3_arg[i];    // scalar product between p1_arg[] and p3_arg[], needed for construction of unit vector in x-direction
  }

  // unit vector ex[] in x-direction
  // x-direction roughly points in the direction of p3_arg[] (see notes)
  // construction is via: ex[] = p3_arg[] - ( (p1_arg[]*p3_arg[]) / (p1_arg*p1_arg) ) * p1_arg[]
  p1p3 = p1p3 / pow( P1_arg[0], 2 );

  double norm = 0;
  for ( int j = 1; j <= 3; j++ )
  {
    ex[j] = P3_arg[j] - p1p3 * P1_arg[j];
    norm += ex[j] * ex[j];
  }
  // normalize to unit vector
  for ( int k = 1;k <= 3;k++ )
    ex[k] = ex[k] / sqrt( norm );

  // unit vector ey[] in y-direction
  // ey[] is contructed via vector product ez[] x ex[]
  ey[1] = ez[2] * ex[3] - ez[3] * ex[2];
  ey[2] = ex[1] * ez[3] - ex[3] * ez[1];
  ey[3] = ez[1] * ex[2] - ez[2] * ex[1];
  //
  //------------------


  // the boost vector (given as beta_vec, a private member of scattering32) is transformed (rotated) to the new reference frame
  // given by ex, ey and ez
  // construction is via: a'[i] = V[i,j]*a[j] (sum over j) with V[i,j] = scalar_product(e'_[i],e_[j])
  // e'_[1] for example corresponds to the vector ex[], whereas e_[1] would be the normal unit vector (1,0,0)
  //
  beta_new[1] = beta_vec[1] * ex[1] + beta_vec[2] * ex[2] + beta_vec[3] * ex[3];
  beta_new[2] = beta_vec[1] * ey[1] + beta_vec[2] * ey[2] + beta_vec[3] * ey[3];
  beta_new[3] = beta_vec[1] * ez[1] + beta_vec[2] * ez[2] + beta_vec[3] * ez[3];
}


// returns N, thus the number indicating the actual order of P1, P2, P3
int scattering32::getOrder() const
{
  return N;
}


double scattering32::getRatio() const
{
  return total_ratio;
}


bool scattering32::getCollisionStatus() const
{
  return collision;
}
// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;
