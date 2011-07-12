//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
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
#include <math.h>

#include "scattering23.h"
#include "lorentz.h"
#include "FPT_compare.h"
#include "random.h"
#include "integrand23.h"
#include "integrate.h"
#include "prefactors23.h"

using std::cout;
using std::endl;


/**
 * This basic constructor sets all pointers to NULL and only sets up the interpolation23 object.
 * It MUST be used with scattering23::setParameter!
 *
 * @param[in] theI23_arg Pointer to an interpolation23 object.
 */
scattering23::scattering23( const interpolation23 * const theI23_arg )
    : theI23( theI23_arg ), P1( NULL ), P2( NULL ), P1cell( NULL ), P2cell( NULL ), P1cm( NULL ), P2cm( NULL )
{
}


/**
*
* @param[in] theI23_arg Pointer to an interpolation23 object.
* @param[in] vx x-component of the collective velocity of the cell
* @param[in] vy y-component of the collective velocity of the cell
* @param[in] vz z-component of the collective velocity of the cell
* @param[in] P1arg[] momentum vector of particle 1
* @param[in] P2arg[] momentum vector of particle 2
* @param[in] sqrtS_arg center of mass energy, sqrt(s)
* @param[in] md2g_scaled_arg squared debye mass, scaled by 1/s
* @param[in] lambda_scaled_arg mean free path, scaled by sqrt(s)
*/
scattering23::scattering23( const interpolation23 * const theI23_arg, const double vx, const double vy, const double vz, const double P1arg[], const double P2arg[],
                            const FLAVOR_TYPE F1_arg, const FLAVOR_TYPE F2_arg,
                            const double sqrtS_arg, const double md2g_scaled_arg, const double lambda_scaled_arg )
    : theI23( theI23_arg ), F1( F1_arg ), F2( F2_arg ), sqrtS( sqrtS_arg ), md2g_scaled( md2g_scaled_arg ), lambda_scaled( lambda_scaled_arg )
{
  P1 = new double[4];
  P2 = new double[4];
  P1cell = new double[4];
  P2cell = new double[4];
  P1cm = new double[4];
  P2cm = new double[4];

  cell_beta_vec[0] = 0.0;
  cell_beta_vec[1] = vx;
  cell_beta_vec[2] = vy;
  cell_beta_vec[3] = vz;

  for ( int i = 0; i <= 3; i++ )
  {
    P1[i] = P1arg[i];
    P2[i] = P2arg[i];
  }

  // boost momenta to rest frame of the cell
  lorentz( cell_beta_vec, P1, P1cell );
  lorentz( cell_beta_vec, P2, P2cell );

  double totE = P1cell[0] + P2cell[0];
  for ( int i = 1; i <= 3; i++ )
    beta_vec[i] = ( P1cell[i] + P2cell[i] ) / totE;
  beta = sqrt( pow( beta_vec[1], 2.0 ) + pow( beta_vec[2], 2.0 ) + pow( beta_vec[3], 2.0 ) );

  lorentz( beta_vec, P1cell, P1cm );
  lorentz( beta_vec, P2cell, P2cm );

  double abs_P1cm = sqrt( pow( P1cm[1], 2.0 ) + pow( P1cm[2], 2.0 ) + pow( P1cm[3], 2.0 ) );
  cos_theta = ( P1cm[1] * beta_vec[1] + P1cm[2] * beta_vec[2] + P1cm[3] * beta_vec[3] ) / ( beta * abs_P1cm );

  theta = acos( cos_theta );
}


scattering23::~scattering23()
{
  delete[] P1;
  P1 = NULL;
  delete[] P2;
  P2 = NULL;
  delete[] P1cell;
  P1cell = NULL;
  delete[] P2cell;
  P2cell = NULL;
  delete[] P1cm;
  P1cm = NULL;
  delete[] P2cm;
  P2cm = NULL;
}


/**
* This method sets all necessary parameters for a given particle pair. Previous values are deleted or overwritten.
* Using this method an scattering32 object can be re-used for multiple particle pair, thus reducing the need to constantly
* creating new objects.
*
* Either this method or the constructor taking the same arguments MUST be called prior to any other methods of the class!
*
* @param[in] vx x-component of the collective velocity of the cell
* @param[in] vy y-component of the collective velocity of the cell
* @param[in] vz z-component of the collective velocity of the cell
* @param[in] P1arg[] momentum vector of particle 1
* @param[in] P2arg[] momentum vector of particle 2
* @param[in] sqrtS_arg center of mass energy, sqrt(s)
* @param[in] md2g_scaled_arg squared debye mass, scaled by s
* @param[in] lambda_scaled_arg mean free path, scaled by sqrt(s)
* @return The absolute value of the boost velocity between cell frame and CMS (#beta_vec).
*/
double scattering23::setParameter( const double vx, const double vy, const double vz, const double P1arg[], const double P2arg[],
                                   const FLAVOR_TYPE F1_arg, const FLAVOR_TYPE F2_arg,
                                   const double sqrtS_arg, const double md2g_scaled_arg, const double lambda_scaled_arg )
{
  delete[] P1;
  P1 = NULL;
  delete[] P2;
  P2 = NULL;
  delete[] P1cell;
  P1cell = NULL;
  delete[] P2cell;
  P2cell = NULL;
  delete[] P1cm;
  P1cm = NULL;
  delete[] P2cm;
  P2cm = NULL;

  F1 = F1_arg;
  F2 = F2_arg;
  P1 = new double[4];
  P2 = new double[4];
  P1cell = new double[4];
  P2cell = new double[4];
  P1cm = new double[4];
  P2cm = new double[4];

  cell_beta_vec[0] = 0.0;
  cell_beta_vec[1] = vx;
  cell_beta_vec[2] = vy;
  cell_beta_vec[3] = vz;

  for ( int i = 0; i <= 3; i++ )
  {
    P1[i] = P1arg[i];
    P2[i] = P2arg[i];
  }

  sqrtS = sqrtS_arg;
  md2g_scaled = md2g_scaled_arg;
  lambda_scaled = lambda_scaled_arg;

  // boost momenta to rest frame of the cell
  lorentz( cell_beta_vec, P1, P1cell );
  lorentz( cell_beta_vec, P2, P2cell );

  double totE = P1cell[0] + P2cell[0];
  for ( int i = 1; i <= 3; i++ )
    beta_vec[i] = ( P1cell[i] + P2cell[i] ) / totE;
  beta = sqrt( pow( beta_vec[1], 2.0 ) + pow( beta_vec[2], 2.0 ) + pow( beta_vec[3], 2.0 ) );

  lorentz( beta_vec, P1cell, P1cm );
  lorentz( beta_vec, P2cell, P2cm );

  double abs_P1cm = sqrt( pow( P1cm[1], 2.0 ) + pow( P1cm[2], 2.0 ) + pow( P1cm[3], 2.0 ) );
  cos_theta = ( P1cm[1] * beta_vec[1] + P1cm[2] * beta_vec[2] + P1cm[3] * beta_vec[3] ) / ( beta * abs_P1cm );

  theta = acos( cos_theta );

  return beta;
}


/**
 * @return Integral over the 2->3 matrix element, I23.
 */
double scattering23::getIntegral23( int& initialStateIndex ) const
{
  double I23_gg_ggg = 0;

  // sort F1 and F2 such that comparisons below are easier
  unsigned int _F1 = std::min( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );
  unsigned int _F2 = std::max( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );

  double total_prefactor = 0;
  initialStateIndex = -1;
  

  if (( _F1 + _F2 ) == 0 ) // g+g -> g+g+g, g+g -> q+qbar+g
  {
    prefactor23_gg_ggg pObj1;
    prefactor23_gg_qqbarg pObj2;
    
    initialStateIndex = 0;

    total_prefactor = pObj1.prefactor() * pObj1.symFactor23() + pObj2.prefactor() * pObj2.symFactor23();
  }
  else if ( _F1 == _F2 )  // q+q -> q+q+g, qbar+qbar -> qbar+qbar+g
  {
    prefactor23_qq_qqg pObj1;
    
    initialStateIndex = 1;

    total_prefactor = pObj1.prefactor() * pObj1.symFactor23();
  }
  else if (( _F1 * _F2 ) == 0 ) // g+q -> g+q+g, g+qbar -> g+qbar+g
  {
    prefactor23_qg_qgg pObj1;
    
    initialStateIndex = 2;

    total_prefactor = pObj1.prefactor() * pObj1.symFactor23();
  }
  else if (( _F2 - _F1 ) == 1 && ( _F2 % 2 ) == 0 )  // q+qbar -> q+qbar+g, q+qbar -> q'+qbar'+g, q+qbar -> g+g+g
  {
    prefactor23_qqbar_qqbarg pObj1;
    prefactor23_qqbar_qqbarDashg pObj2;
    prefactor23_qqbar_ggg pObj3;
    
    initialStateIndex = 3;

    total_prefactor = pObj1.prefactor() * pObj1.symFactor23() + pObj2.prefactor() * pObj2.symFactor23() + pObj3.prefactor() * pObj3.symFactor23();
  }
  else // q+q'+g -> q+q'+g, q+qbar' -> q+qbar'+g
  {
    prefactor23_qqdash_qqdashg pObj1;
    
    initialStateIndex = 4;

    total_prefactor = pObj1.prefactor() * pObj1.symFactor23();
  }


  if ( total_prefactor > 0 )
  {
    if ( true )
    {
      I23_gg_ggg = theI23->getI23( log( md2g_scaled ), log( lambda_scaled ), beta, fabs( cos_theta ) );
    }
    else
    {
      integrand23 theIntegrand;
      theIntegrand.set_md2( md2g_scaled );
      theIntegrand.set_lambda( lambda_scaled );
      theIntegrand.set_cos_theta( fabs( cos_theta ) );
      theIntegrand.set_beta( beta );

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

      integrate( theIntegrand, neval, fail, intResult, error, prob );

      I23_gg_ggg = ( 9.0 / ( 2.0 * M_PI ) * intResult[0] );
    }
  }

  return ( total_prefactor * I23_gg_ggg );
}


/**
* Samples new momenta of the outgoing particles according to the matrix element using the rejection method.
*
* @param[out] pt1 Sampled qt (in GeV)
* @param[out] pt3 Sampled kt (in GeV)
* @param[out] y Sampled rapidity of emmited gluon.
* @param[out] phi Sampled angle of the outgoing gluon.
* @param[out] pz1 Sampled longitudinal momentum transfer.
* @return Count how many time the comparison function was below the actual function. Only needed for rejection method. Should be zero.
*/
int scattering23::getMomenta23( double& pt1, double& pt3, double& y, double& phi, double& pz1, int& typ, FLAVOR_TYPE& F1arg, FLAVOR_TYPE& F2arg ) const
{
  int error = 0;   //0 is returned if no error occurs, 1 otherwise

  const double epsilon = 1.0e-10;
  double deltaTrans_estimate = 0.001;

  bool solvable;

  double max;
  double qt2, qt, kt2, kt, g, gr;
  double c, delta, delta1, delta2;
  double dTf1, dTf2, deltaTrans_factor = 0;   // this is the variable that will hold 1/|dF/dy'_1| where F=0
  //(the factor arising from the transformation of the delta function, hence the name)
  double pz11, pz12, E1, pz3, E3, qtkt;
  double A, B, C;
  double y_min, y_max, kin_y_constraint;
  double V_y;

  double lambda_scaled2 = pow( lambda_scaled, 2.0 );
  if ( 1 / lambda_scaled2 >= 0.25 )
  {
    cout << "lambda_scaled2 = " << lambda_scaled2 << endl;
    std::string errMsg = "1 / lambda_scaled2 >= 0.25";
    throw eScatt23_error( errMsg );
  }
  
  double ymin_kt_independent = 0, ymax_kt_independent = 0;
  maximal_Y_range_23( ymin_kt_independent, ymax_kt_independent );
  double maxV_y = ymax_kt_independent - ymin_kt_independent;

  double max0 = maxV_y / ( md2g_scaled * deltaTrans_estimate );

  double maxV_kt2 = log( lambda_scaled2 / 4.0 );
  double maxV_qt2 = log( 1 / ( 4 * md2g_scaled ) + 1 );



  F1arg = F1;
  F2arg = F2;
  // sort F1 and F2 such that comparisons below are easier
  // convert to unsigned int to allow for some math that speeds up the decisions
  unsigned int _F1 = std::min( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );
  unsigned int _F2 = std::max( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );



  if (( _F1 + _F2 ) == 0 ) // g+g -> g+g+g, g+g -> q+qbar+g
  {
    prefactor23_gg_ggg pObj1;
    prefactor23_gg_qqbarg pObj2;
    double pf1 = pObj1.prefactor() * pObj1.symFactor23();
    double pf2 = pObj2.prefactor() * pObj2.symFactor23();
    
    if ( ran2() < pf1 / (pf1 + pf2) )
    {
      typ = 231; // g+g -> g+g+g
    }
    else
    {
      typ = 232; // g+g -> q+qbar+g
      sampleFlavor23( F1arg, F2arg );
    }
  }
  else if ( _F1 == _F2 )  // q+q -> q+q+g, qbar+qbar -> qbar+qbar+g
  {
    typ = 237;
  }
  else if ( (_F1 * _F2) == 0 ) // g+q -> g+q+g, g+qbar -> g+qbar+g
  {
    typ = 233;
  }
  else if (( _F2 - _F1 ) == 1 && ( _F2 % 2 ) == 0 )  // q+qbar -> q+qbar+g, q+qbar -> q'+qbar'+g, q+qbar -> g+g+g
  {
    prefactor23_qqbar_qqbarg pObj1;
    prefactor23_qqbar_qqbarDashg pObj2;
    prefactor23_qqbar_ggg pObj3;

    double pf1 = pObj1.prefactor() * pObj1.symFactor23();
    double pf2 = pObj2.prefactor() * pObj2.symFactor23();
    double pf3 = pObj3.prefactor() * pObj3.symFactor23();
    
    double select = ran2() * ( pf1 + pf2 + pf3 );
    
    if ( select < pf1 )  // q+qbar -> q+qbar+g
    {
      typ = 234;
    }
    else if ( select < ( pf1 + pf2 ) )  // q+qbar -> q'+qbar'+g
    {
      typ = 235;
      sampleFlavor23( F1arg, F2arg, F1 );  // sample flavor excluding F1 (or anti-F1 = F2)
    }
    else  // q+qbar -> g+g+g
    {
      typ = 236;
      F1arg = gluon;
      F2arg = gluon;
    }
  }
  else // q+q'-> q+q'+g, q+qbar' -> q+qbar'+g
  {
    typ = 238; // qq' -> qq', qqbar' -> qqbar'
  }


  do
  {
    kt2 = 1 / lambda_scaled2 * exp( ran2() * maxV_kt2 );
    qt2 = md2g_scaled * ( exp( ran2() * maxV_qt2 ) - 1.0 );

    qt = sqrt( qt2 );
    kt = sqrt( kt2 );

    solvable = true;
    //<<------------------------------------------
    //
    // limits for y are either determined from the kinematic constraint cosh(y) < sqrt(s)/(2k_t)
    // or from the requirement that
    //                                  cosh(y) + A*sinh(y) < B
    // stemming from the modelling of the LPM effect (see notes), where
    //
    //                                  A = beta*cos(theta)
    // and
    //                                  B = k_t*lambda*sqrt(1-beta^2).
    //
    // It is: B > 0 and 0 <= A <= 1.
    //
    A = beta * cos_theta;
    B = kt * lambda_scaled * sqrt( 1 - pow( beta, 2.0 ) );
    kin_y_constraint = 1 / ( 2 * kt );

    // as B > 0 and 0 <= A <= 1 this should not evaluate to "true"
    if ( FPT_COMP_L( A, -1.0 ) || FPT_COMP_G( A, 1.0 ) || FPT_COMP_L( B, 0.0 ) )
    {
      cout << "encountered unexpected values when determining limits for y" << endl;
      cout << "A = " << A << endl;
      cout << "B = " << B << endl;

      pz1 = 0;
      pt1 = 0;
      pt3 = 0;
      return 1;
    }
    else
    {
      //<<---------------------------------------------
      // some special cases...
      if ( FPT_COMP_E( B, 0.0 ) )              // B = 0
      {
        solvable = false;
      }

      if ( FPT_COMP_E( A, 0.0 ) )              // A = 0
      {
        if ( FPT_COMP_LE( B, 1.0 ) )         // A = 0  and  B <= 1
        {
          solvable = false;
        }
        else                                 // A = 0  and  B > 1
        {
          y_min = -acosh( B );
          y_max = acosh( B );
        }
      }
      //----------------------------------------------->>

      //<<-----------------------------------------------
      // the general case...
      else
      {
        if ( FPT_COMP_LE( B, ( sqrt( 1 + A )*sqrt( 1 - A ) ) ) ) //B <= (sqrt(1+A)*sqrt(1-A)) -> no solution
        {
          solvable = false;
        }
        else
        {
          y_min = log(( B - sqrt( pow( A, 2.0 ) + pow( B, 2.0 ) - 1 ) ) / ( 1 + A ) );
          y_max = log(( B + sqrt( pow( A, 2.0 ) + pow( B, 2.0 ) - 1 ) ) / ( 1 + A ) );
        }
      }
      //----------------------------------------------->>
    }

    //<<-----------------------------------------------
    // check whether the constraints on y from kinematics
    // are stronger than those imposed by the LPM cut-off
    if ( FPT_COMP_LE( kin_y_constraint, 1.0 ) )  // no solution in this case
    {
      solvable = false;
    }
    else
    {
      if ( FPT_COMP_L( acosh( kin_y_constraint ), y_max ) )
        y_max = acosh( kin_y_constraint );
      if ( FPT_COMP_G( -acosh( kin_y_constraint ), y_min ) )
        y_min = -acosh( kin_y_constraint );
    }
    //----------------------------------------------->>

    if ( solvable )
    {
      V_y = y_max - y_min;
      y = y_min + ( y_max - y_min ) * ran2();
    }
    else
    {
      V_y = 0;
      y = 0;
    }
    //----------------------------------------------->>

    phi = M_PI * ran2();

    // Now that q_t, k_t, y and phi have been fixed, proceed by calculating
    // the argument of the sigma_23 integral at these values.
    // See notes (starting 24.10.2006) for more details.
    //
    E3 = kt * cosh( y );     //some short hand notations
    pz3 = kt * sinh( y );
    qtkt = qt * kt * cos( phi );

    /* To evaluate the value of the matrix element, firstly the term 1/|dF/dy'_1| has to be calculated where F/s = 0:
    // Determine p_{1z}' from the condition F/s = 0:
    // After squaring one obtains the equation
    //
    //                 A*p_{1z}'^2 + B*p_{1z}' + C = 0                                      (1)
    //
    // with
    //                 A = (1-E_{3}')^2 - p_{3z}' = 1 - 2*E_{3}' + k_t^2
    //
    //              B = P_{3z}' * ( 1 - 2*E_{3}' + 2*q_t*k_t*cos(phi) )
    //
    //                 C = (1-E_{3}')^2 *q_t^2 - 1/4 * ( 1 - 2*E_{3}' + 2*q_t*k_t*cos(phi) )^2
    //
    // The general solution then reads
    //
    //                 p_{1z}' = ( -B +- sqrt( B^2 - 4*A*C ) ) / 2*A
    //
    // wheras in the case of A = 0 the solution is just
    //
    //                 p_{1z}' = -C/B.
    //
    // Additionally the condition
    //
    //                 1 - 2*E_{3}' + 2*q_t*k_t*cos(theta) - 2*p_{3z}'*p_{1z}' >= 0          (2)
    //
    // (equivalent to   B/p_{3z}' - 2*p_{3z}'*p_{1z}' >= 0)
    //
    // has to be fulfilled (see notes).
    */

    A = 1.0 - 2.0 * E3 + kt2;
    B = pz3 * ( 1.0 - 2.0 * E3 + 2.0 * qtkt );
    C = pow(( 1.0 - E3 ), 2.0 ) * qt2 - 0.25 * pow(( 1.0 - 2.0 * E3 + 2.0 * qtkt ) , 2.0 );

    if ( FPT_COMP_E( A, 0.0 ) )     // check for special case A = 0
    {
      if ( FPT_COMP_E( B, 0.0 ) ) // if A = 0 and B = 0 there is no solution
      {
        deltaTrans_factor = 0;
      }
      else
      {
        pz11 = -C / B;           // the solution for p_{z1}' in the case of A = 0
        delta1 = B / pz3 - 2.0 * pz3 * pz11;
        dTf1 = 0;
        dTf2 = 0;
        if ( FPT_COMP_GE( delta1, 0.0 ) && pz11 > 0 )  // check whether condition (2) is fulfilled
        {
          E1 = sqrt( qt2 + pow( pz11, 2.0 ) );
          c = 2 * fabs( pz11 - E3 * pz11 + E1 * pz3 ); // 1/c is the argument of the sum stemming from the transformation of the delta function
          // i.e., fabs(dF/dy_1') = 2*c where F=0
          if ( c < epsilon )
            deltaTrans_factor = 1.0 / epsilon;
          else
            deltaTrans_factor = 1.0 / c;
          dTf1 = deltaTrans_factor;
        }
        else
        {
          deltaTrans_factor = 0;
        }
      }
    }
    else                                        // the general case, i.e. A != 0
    {
      delta = pow( B, 2.0 ) - 4.0 * A * C;    // the argument of sqrt(..) in the general solution for p_{z1}'
      if ( delta < 0.0 )
      {
        deltaTrans_factor = 0;
      }
      else
      {
        // solution 1
        pz11 = ( -B + sqrt( delta ) ) / ( 2.0 * A );
        delta1 = B / pz3 - 2.0 * pz3 * pz11;    // condition (2) requires delta1 >= 0 in order for pz11 to be a valid solution

        // solution 2
        pz12 = ( -B - sqrt( delta ) ) / ( 2.0 * A );
        delta2 = B / pz3 - 2.0 * pz3 * pz12;    // condition (2) requires delta2 >= 0 in order for pz12 to be a valid solution

        deltaTrans_factor = 0;
        dTf1 = 0;
        dTf2 = 0;
        if ( delta1 >= 0.0 && pz11 > 0 )
        {
          E1 = sqrt( qt2 + pow( pz11, 2.0 ) );
          c = 2 * fabs( pz11 - E3 * pz11 + E1 * pz3 );
          if ( c < epsilon )
            dTf1 = 1.0 / epsilon;
          else
            dTf1 = 1.0 / c;
          deltaTrans_factor += dTf1;

        }
        if ( delta2 >= 0.0 && pz12 > 0 )
        {
          E1 = sqrt( qt2 + pow( pz12, 2.0 ) );
          c = 2 * fabs( pz12 - E3 * pz12 + E1 * pz3 );
          if ( c < epsilon )
            dTf2 = 1.0 / epsilon;
          else
            dTf2 = 1.0 / c;
          deltaTrans_factor += dTf2;
        }
      }
    }


    g = V_y * deltaTrans_factor * qt2 / pow(( qt2 + md2g_scaled ), 2 ) / kt2 / ( kt2 + qt2 + md2g_scaled - 2.0 * qtkt );

    max = max0 / ( qt2 + md2g_scaled ) / kt2;
    if ( g > max )
    {
      //cout << "error in get23" << endl;
      ++error;
    }

    gr = max * ran2();

  }
  while ( g < gr );

  if ( ran2() < dTf1 / ( dTf1 + dTf2 ) )
  {
    pz1 = pz11 * sqrtS;
  }
  else
  {
    pz1 = pz12 * sqrtS;
  }

  pt1 = qt * sqrtS;
  pt3 = kt * sqrtS;

  return error;
}


/**
* Sets new momenta for the outgoing particle according to the values sampled in getNewMomenta23.
* The results are written to the vectors P1, P2 and P3.
*
* @param[out] P1[] Momentum vector of outgoing particle 1
* @param[out] P2[] Momentum vector of outgoing particle 2
* @param[out] P3[] Momentum vector of outgoing particle 2
* @param[in] R1[] Space-Time vector of ingoing particle 1.
* @param[in] R2[] Space-Time vector of ingoing particle 2.
* @param[in] PT1 Momentum transfer qt
* @param[in] PT3 Momentum transfer kt
* @param[in] y3 Rapidity of emitted gluon.
* @param[in] phi azimuthal angle
* @param[in] PZ1 Longitudinal momentum transfer.
*/
void scattering23::setNewMomenta23( double P1[4], double P2[4], double P3[4], const double R1[4], const double R2[4],
                                    const double PT1, const double PT3, const double y3, const double phi, const double PZ1 )
{
  double PZ3, sinus, cosinus;
  double c = 0;
  double R1cell[4], R2cell[4], R1cm[4], R2cm[4], P3cm[4], P3cell[4];
  double PP[4], TT[4], transv[4];

  lorentz( cell_beta_vec, R1, R1cell );
  lorentz( cell_beta_vec, R2, R2cell );
  lorentz( beta_vec, R1cell, R1cm );
  lorentz( beta_vec, R2cell, R2cm );

  rotation( P1cm, R1cm, R2cm, PP, TT );

  for ( int i = 1;i <= 3;i++ )
  {
    P1cm[i] = PT1 * TT[i] + PZ1 * PP[i];
    c += P1cm[i] * P1cm[i];
  }
  P1cm[0] = sqrt( c );

  for ( int j = 1;j <= 3;j++ )
    TT[j] = -TT[j];
  sinus = sin( phi );
  cosinus = cos( phi );
  transv[1] = TT[1] * cosinus + ( PP[2] * TT[3] - PP[3] * TT[2] ) * sinus;
  transv[2] = TT[2] * cosinus + ( PP[3] * TT[1] - PP[1] * TT[3] ) * sinus;
  transv[3] = TT[3] * cosinus + ( PP[1] * TT[2] - PP[2] * TT[1] ) * sinus;

  for ( int j = 1;j <= 3;j++ )
    TT[j] = transv[j];

  PZ3 = PT3 * sinh( y3 );

  c = 0;
  for ( int k = 1;k <= 3;k++ )
  {
    P3cm[k] = PT3 * TT[k] + PZ3 * PP[k];
    c += P3cm[k] * P3cm[k];
  }
  P3cm[0] = sqrt( c );

  c = 0;
  for ( int i = 1;i <= 3;i++ )
  {
    P2cm[i] = -( P1cm[i] + P3cm[i] );
    c += P2cm[i] * P2cm[i];
  }
  P2cm[0] = sqrt( c );

  // boost new momentum vectors back to original frame
  double inv_beta_cell[4];
  for ( int j = 1;j <= 3;j++ )
    inv_beta_cell[j] = -cell_beta_vec[j];

  double inv_beta[4];
  for ( int j = 1;j <= 3;j++ )
    inv_beta[j] = -beta_vec[j];

  lorentz( inv_beta, P1cm, P1cell );
  lorentz( inv_beta, P2cm, P2cell );
  lorentz( inv_beta, P3cm, P3cell );

  lorentz( inv_beta_cell, P1cell, P1 );
  lorentz( inv_beta_cell, P2cell, P2 );
  lorentz( inv_beta_cell, P3cell, P3 );
}



/**
* Returns the maximal and minimal values of the rapidity y possible for ANY k_t (within the appropriate ranges). These
* values are returned via the ymin and ymax passed to maximal_Y_range_23 as references.
*
* y is constrained by kinematics
*
*    cosh(y) < 1/(2*k_t)
*
* and by the LPM cutoff condition (see notes of 18.10.06)
*
*    cosh(y) + A*sinh(y) < B
*
* with the solution
*                    ln((B - sqrt(A^2 + B^2 - 1))/(1+A)) < y < ln((B + sqrt(A^2 + B^2 - 1))/(1+A))    (1)
*
* where
*                                  A = beta*cos(theta)
* and
*                                  B = k_t*lambda*sqrt(1-beta^2).
*
* The extremal values for y are therefore given by the intersections of the kinematic and the LPM constraints (see notes
* of 30.01.07 for more details). The location of the intersection points on the k_t (the B) axis are thus given by the
* solutions to the equations
*                                  acosh(g/(2*B)) == ln((B + sqrt(A^2 + B^2 - 1))/(1+A))     (2)
*                                  -acosh(g/(2*B)) == ln((B - sqrt(A^2 + B^2 - 1))/(1+A))    (3)
* with g = lambda*sqrt(1-beta^2) = lambda/gamma.
* The solution to any of these two equations can be written as (notes)
*
*                                  B = 1/sqrt(2) * sqrt( g - A2 +/- sqrt( A^2*(A^2+g*(g-2)) ) )
*
* The appropriate sign has to be fixed by substituting the values into the acosh(..) and ln(..) and comparing (see below).
*
* @param[out] ymin Returns the minmial y value.
* @param[out] ymax Returns the maximal y value.
* @return Whether a valid range could be successfully determined.
*/
bool scattering23::maximal_Y_range_23( double& ymin, double& ymax ) const
{
  double A = beta * cos_theta;             // definition of A and g as given above
  double A2 = pow( A, 2 );
  double g = lambda_scaled * sqrt( 1 - pow( beta, 2 ) );

  bool success = true;

  ymin = 0;
  ymax = 0;

  // special case A=0
  if ( FPT_COMP_E( A, 0 ) )
  {
    double B = sqrt( g / 2.0 );  // solution of acosh(B) = acosh(g/(2*B))
    ymax = acosh( B );
    ymin = -acosh( B );
  }

  // the general case
  else
  {
    double B1 = 1 / sqrt( 2.0 ) * sqrt( g - A2 + sqrt( A2 * ( A2 + g * ( g - 2 ) ) ) );  // solution to eqns. (2) and (3)
    double B2 = 1 / sqrt( 2.0 ) * sqrt( g - A2 - sqrt( A2 * ( A2 + g * ( g - 2 ) ) ) );  // solution to eqns. (2) and (3)

    //cout << "kt_1 = " << B1 / g << "   kt_2 = " << B2 / g << endl;

    double intersection1, intersection2;
    double yacosh1, yacosh2, y1, y2;
    //<<-----------------------------------------------
    // determine intersection point 1
    y1 = log(( B1 + sqrt( A2 + pow( B1, 2 ) - 1 ) ) / ( 1 + A ) );          // y corresponding to B1 and the boundary (1) with + sign
    y2 = log(( B1 - sqrt( A2 + pow( B1, 2 ) - 1 ) ) / ( 1 + A ) );          // y corresponding to B1 and the boundary (1) with - sign
    yacosh1 = acosh( g / ( 2 * B1 ) );                         // y corresponding to B1 and the kinematic upper boundary
    yacosh2 = -acosh( g / ( 2 * B1 ) );                         // y corresponding to B1 and the kinematic lower boundary

    // check which of {y1,y2} equals any of {yacosh1,yacosh2}
    if ( FPT_COMP_E( y1, yacosh1 ) )
      intersection1 = y1;
    else if ( FPT_COMP_E( y2, yacosh1 ) )
      intersection1 = y2;
    else if ( FPT_COMP_E( y1, yacosh2 ) )
      intersection1 = y1;
    else if ( FPT_COMP_E( y2, yacosh2 ) )
      intersection1 = y2;
    else
    {
      intersection1 = 0;
      success = false;
    }
    //----------------------------------------------->>

    //<<-----------------------------------------------
    // determine intersection point2
    y1 = log(( B2 + sqrt( A2 + pow( B2, 2 ) - 1 ) ) / ( 1 + A ) );          // y corresponding to B2 and the boundary (1) with + sign
    y2 = log(( B2 - sqrt( A2 + pow( B2, 2 ) - 1 ) ) / ( 1 + A ) );          // y corresponding to B2 and the boundary (1) with - sign
    yacosh1 = acosh( g / ( 2 * B2 ) );                        // y corresponding to B2 and the kinematic upper boundary
    yacosh2 = -acosh( g / ( 2 * B2 ) );                        // y corresponding to B2 and the kinematic lower boundary

    // check which of {y1,y2} equals any of {yacosh1,yacosh2}
    if ( FPT_COMP_E( y1, yacosh1 ) )
      intersection2 = y1;
    else if ( FPT_COMP_E( y2, yacosh1 ) )
      intersection2 = y2;
    else if ( FPT_COMP_E( y1, yacosh2 ) )
      intersection2 = y1;
    else if ( FPT_COMP_E( y2, yacosh2 ) )
      intersection2 = y2;
    else
    {
      intersection2 = 0;
      success = false;
    }
    //----------------------------------------------->>

    if ( intersection1 > intersection2 )
    {
      ymax = intersection1;
      ymin = intersection2;
    }
    else
    {
      ymax = intersection2;
      ymin = intersection1;
    }
  }

  return success;
}


/** @todo document this function */
void scattering23::rotation( const double P[4], const double R1[4], const double R2[4], double PP[4], double TT[4] ) const
{

  double min = 1.0, phi, sinus, cosinus, c1, c2, c3;
  double transv[4];
  int mini = -1;


  c1 = c2 = 0.0;
  for ( int i = 1;i <= 3;i++ )
  {
    c1 += P[i] * P[i];
    c2 += P[i] * ( R2[i] - R1[i] );
  }

  c3 = 0.0;
  for ( int j = 1;j <= 3;j++ )
  {
    PP[j] = P[j] / sqrt( c1 );
    TT[j] = c1 * ( R1[j] - R2[j] ) + c2 * P[j];
    c3 += TT[j] * TT[j];
  }
  c3 = sqrt( c3 );

  // r2-r1 // p1 => TT[] = 0
  if ( c3 < 1.0e-8 )
  {

    for ( int i = 1;i <= 3;i++ )
    {
      if ( fabs( PP[i] ) < min )
      {
        min = fabs( PP[i] );
        mini = i;
      }
    }
    switch ( mini )
    {
    case 1:
      TT[1] = 0.0;
      TT[2] = PP[3];
      TT[3] = -PP[2];
      break;
    case 2:
      TT[1] = -PP[3];
      TT[2] = 0.0;
      TT[3] = PP[1];
      break;
    case 3:
      TT[1] = PP[2];
      TT[2] = -PP[1];
      TT[3] = 0.0;
      break;
    default:
      cout << "Error in rotation()" << endl;
    }
  }
  //----------------------

  phi = 2.0 * M_PI * ran2();
  sinus = sin( phi );
  cosinus = cos( phi );
  transv[1] = TT[1] * cosinus + ( PP[2] * TT[3] - PP[3] * TT[2] ) * sinus;
  transv[2] = TT[2] * cosinus + ( PP[3] * TT[1] - PP[1] * TT[3] ) * sinus;
  transv[3] = TT[3] * cosinus + ( PP[1] * TT[2] - PP[2] * TT[1] ) * sinus;

  c3 = 0.0;
  for ( int j = 1;j <= 3;j++ )
    c3 += transv[j] * transv[j];
  c3 = sqrt( c3 );

  for ( int i = 1;i <= 3;i++ )
    TT[i] = transv[i] / c3;
}
