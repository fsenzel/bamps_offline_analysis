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
#include <iostream>

#include "integrand23.h"
#include "FPT_compare.h"

using std::cout;
using std::endl;



// overloaded operator() that makes integrand23 a function object - for use with CUBA-vegas
void integrand23::operator()(const int *ndim, const double xx[], const int *ncomp, double ff[]) const
{
  double wgt;
  ff[0] = this->operator()(xx, wgt);  
}


/**
 * This function operator returns the matrix element 2->3 for given qt2, kt2, y and phi (qt2 and kt2 scaled by s to be dimensionless).
 * Arguments are expected to be in the range [0,1], thus indicating the relative position e.g. between qt2_min and qt2_max.
 * This scaled treatment is needed for use with numerical integration routines.
 *
 * @param[in] xx[] unit-offset array with function arguments (xx[1] = qt2, xx[2] = kt2, xx[3] = y, xx[4] = phi)
 * @param[in] wgt there for compatibility reasons - not used at the moment
 */
double integrand23::operator()(const double xx[], double wgt) const
{
  const double epsilon=1.0e-10;
  double V;
  double delta,delta1,delta2;
  double qt2,qt,kt2,kt,y,phi;
  double E1,pz11,pz12,pz1,E3,pz3,qtkt;

  //--------------------------------------------
  // q_t^2 is integrated from 0 to 1/4
  double qt2_min = 0;
  double qt2_max = 0.25;
  if (qt2_min >= qt2_max)
  {
    return 0;
  }
      
  V = (qt2_max - qt2_min);
  qt2 = qt2_min + (qt2_max - qt2_min)*xx[1];
  qt = sqrt(qt2);
  //--------------------------------------------


  //--------------------------------------------
  // k_t^2 is integrated from 1/lambda^2 to 1/4
  double kt2_min = log(1.0/pow(lambda_scaled_int,2.0)); //temporary change
  double kt2_max = log(0.25);                 //temporary change
/*  double kt2_min = 1.0/pow(lambda_scaled_int,2.0);
  double kt2_max = 0.25;            */     
  if (kt2_min >= kt2_max)
  {
    return 0;
  }

  V = V*(kt2_max - kt2_min);
  kt2 = kt2_min + (kt2_max - kt2_min)*xx[2];
  kt2 = exp(kt2); //temporary change
  kt = sqrt(kt2);
  //--------------------------------------------


  //<<------------------------------------------
  //
  // Integration limits for y are either determined from the kinematic constraint cosh(y) < sqrt(s)/(2k_t)
  // or from the requirement that 
  //                                  cosh(y) + A*sinh(y) < B
  // stemming from the modelling of the LPM effect, where
  //
  //                                  A = beta'*cos(theta)
  // and
  //                                  B = k_t*lambda*sqrt(1-beta'^2).
  //
  // It is: B > 0 and 0 <= A <= 1. For the physical meaning of beta' and theta see above.
  //
  double A = beta_int * cos_theta_int;
  double B = kt*lambda_scaled_int * sqrt(1-pow(beta_int,2.0));
  double y_min, y_max;
 
  
  double kin_y_constraint = 1/(2*kt);


  // as B > 0 and 0 <= A <= 1 this should not evaluate to "true"
  if ( FPT_COMP_L(A,0.0) || FPT_COMP_G(A,1.0) || FPT_COMP_L(B,0.0) )
  {
    cout << "encountered unexpected values when determining limits for y integration" << endl;
    cout << "A = " << A << endl;
    cout << "B = " << B << endl;

    return 0;
  }
  else
  {  
      //<<---------------------------------------------
      // some special cases...
    //
    if ( FPT_COMP_E(B,0.0) )                 // B = 0
    {
      return 0;
    }
  
    if ( FPT_COMP_E(A,0.0) )                 // A = 0
    {
      if ( FPT_COMP_LE(B,1.0) )            // A = 0  and  B <= 1
      {
        return 0;
      }
      else                                 // A = 0  and  B > 1
      {
        y_min = -acosh(B);
        y_max = acosh(B);
      }
    }
    else
    {
      if ( FPT_COMP_LE(B,(sqrt(1+A)*sqrt(1-A))) ) //B <= (sqrt(1+A)*sqrt(1-A)) -> no solution 
      {
        return 0;
      }
      else
      {
        y_min = log( (B - sqrt(pow(A,2.0)+pow(B,2.0)-1))/(1+A) );
        y_max = log( (B + sqrt(pow(A,2.0)+pow(B,2.0)-1))/(1+A) );
      }
    }
      //----------------------------------------------->>
  }

  //<<-----------------------------------------------
  // check whether the constraints on y from kinematics 
  // are stronger than those imposed by the LPM cut-off
  if ( FPT_COMP_LE(kin_y_constraint,1.0) )  // no solution in this case
  {
    return 0;
  }
  else
  {
    if ( FPT_COMP_G(-acosh(kin_y_constraint),y_max) || FPT_COMP_L(acosh(kin_y_constraint),y_min) )
    {
      return 0;  
    }
    else
    {
      if ( FPT_COMP_L(acosh(kin_y_constraint),y_max) )
        y_max = acosh(kin_y_constraint);
      if ( FPT_COMP_G(-acosh(kin_y_constraint),y_min) )
        y_min = -acosh(kin_y_constraint);
    }
  }
  //----------------------------------------------->>
  
  V = V*(y_max - y_min);
  y = y_min + (y_max - y_min)*xx[3];
  //y = y_max*xx[3];  //achtung!!! tempor�re �nderung!!! wieder r�ckg�ngig machen!!!!
  //------------------------------------------>>


  //<<------------------------------------------
  // phi is integrated from 0 to Pi
  double phi_min = 0;
  double phi_max = M_PI;
  if (phi_min >= phi_max)
  {
    return 0;
  }

  V = V*(phi_max - phi_min);
  phi = phi_min + (phi_max - phi_min)*xx[4];
  //------------------------------------------>>


  // Now that q_t, k_t, y and phi have been fixed, proceed by calculating 
  // the argument of the sigma_23 integral at these values.
  // See notes (starting 24.10.2006) for more details. 
  //
  E3 = kt*cosh(y);         //some short hand notations
  pz3 = kt*sinh(y);
  qtkt = qt*kt*cos(phi);
  
  // To evaluate the integral, firstly the term 1/|dF/dy'_1| has to be calculated where F/s = 0:
  // Determine p_{1z}' from the condition F/s = 0:
  // After squaring one obtains the equation
  //
  //                 A*p_{1z}'^2 + B*p_{1z}' + C = 0                                      (1)
  // 
  // with
  //                 A = (1-E_{3}')^2 - p_{3z}' = 1 - 2*E_{3}' + k_t^2
  //
  //	             B = P_{3z}' * ( 1 - 2*E_{3}' + 2*q_t*k_t*cos(phi) )
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
  //

  double c;
  double deltaTrans_factor = 0;   // this is the variable that will hold 1/|dF/dy'_1| where F=0
                                  //(the factor arising from the transformation of the delta function, hence the name)
  
  A = 1.0 - 2.0*E3 + kt2;
  B = pz3 * ( 1.0 - 2.0*E3 + 2.0*qtkt );
  double C = pow((1.0 - E3),2.0)*qt2 - 0.25*pow( ( 1.0 - 2.0*E3 + 2.0*qtkt ) ,2.0);
  
  if ( FPT_COMP_E(A,0.0) )        // check for special case A = 0
  {
    if( FPT_COMP_E(B,0.0) )     // if A = 0 and B = 0 there is no solution 
    {
      return 0;
    }
    else 
    {
      pz1 = -C/B;             // the solution for p_{z1}' in the case of A = 0 
      delta1 = B/pz3 - 2.0*pz3*pz1;
      if( FPT_COMP_GE(delta1,0.0) )  // check whether condition (2) is fulfilled
      {
        E1 = sqrt(qt2 + pow(pz1,2.0));
        c = 2*fabs(pz1 - E3*pz1 + E1*pz3); // 1/c is the argument of the sum stemming from the transformation of the delta function
	                                         // i.e., fabs(dF/dy_1') = 2*c where F=0
        if(c < epsilon) 
          deltaTrans_factor = 1.0/epsilon;
        else 
          deltaTrans_factor = 1.0/c;
      }
      else 
      {
        return 0;
      }
    }
  }
  else                                        // the general case, i.e. A != 0 
  {
    delta = pow(B,2.0) - 4.0*A*C;           // the argument of sqrt(..) in the general solution for p_{z1}'
    if(delta < 0.0) 
    {
      return 0;
    }
    else
    {
	  // solution 1
      pz11 = (-B + sqrt(delta)) / (2.0*A);
      delta1 = B/pz3 - 2.0*pz3*pz11;    // condition (2) requires delta1 >= 0 in order for pz11 to be a valid solution
	  
	  // solution 2
      pz12 = (-B - sqrt(delta)) / (2.0*A);
      delta2 = B/pz3 - 2.0*pz3*pz12;    // condition (2) requires delta2 >= 0 in order for pz12 to be a valid solution
	  
      if(delta1 >= 0.0)
      {
        E1 = sqrt(qt2 + pow(pz11,2.0));
        c = 2*fabs(pz11 - E3*pz11 + E1*pz3);
        if(c < epsilon) 
          deltaTrans_factor += 1.0/epsilon;
        else 
          deltaTrans_factor += 1.0/c;
      }
      if(delta2 >= 0.0)
      {
        E1 = sqrt(qt2 + pow(pz12,2.0));
        c = 2*fabs(pz12 - E3*pz12 + E1*pz3);
        if(c < epsilon) 
          deltaTrans_factor += 1.0/epsilon;
        else 
          deltaTrans_factor += 1.0/c;
      }
    }
  }
  
  return ( V * deltaTrans_factor * qt2 / pow((qt2+md2_int),2.0) * 1/(kt2 + qt2 + md2_int - 2.0*qtkt) ); //temporary change
  //return ( V * deltaTrans_factor * qt2 / pow((qt2+md2_int),2.0) * 1/(kt2 * (kt2 + qt2 + md2_int - 2.0*qtkt)) );
}

