#include <iostream>
#include <math.h>
#include "interpolation.h"

using std::cout;
using std::endl;

interpolation::interpolation()
{
}

interpolation::interpolation(const vector<double>& x, const vector<double>& y)
{
  spline_construction(x, y);
}


interpolation::interpolation(const vector<double>& x, const vector<double>& y, const double S_initial, const double S_final)
{
  spline_construction(x, y, S_initial, S_final);
}


interpolation::~interpolation()
{
}


void interpolation::spline_construction( const vector<double>& x, const vector<double>& y, const double S_initial, const double S_final)
{
  if ( x.size() != y.size() )
  {
    cout << "FATAL ERROR: vector sizes do not match in interpolation::spline_construction" << endl;
    return;
  }
    
  int n = x.size();
  vector<double> u(n, 0.0);   // auxiliary vector
  S.resize(n, 0.0);

  double p, qn, sig;
  
  if ( S_initial > 1e19 )
  {
    S[0] = 0;
    u[0] = 0;
  }
  else
  {
    S[0] = -0.5;
    u[0] = ( 3 / ( x[1] - x[0] ) ) * (( y[1] - y[0] ) / ( x[1] - x[0] ) - S_initial );
  }

  // decomposition loop of the tridiagonal algorithm
  // S and u are used for temporary storage
  for ( int i = 1; i < (n - 1); i++ )
  {
    sig = ( x[i] - x[i-1] ) / ( x[i+1] - x[i-1] );
    p = sig * S[i-1] + 2;
    S[i] = ( sig - 1 ) / p;
    u[i] = ( y[i+1] - y[i] ) / ( x[i+1] - x[i] ) - ( y[i] - y[i-1] ) / ( x[i] - x[i-1] );
    u[i] = ( 6 * u[i] / ( x[i+1] - x[i-1] ) ) - sig * u[i-1] / p;
  }

  if ( S_final > 1e19 )
  {
    qn = 0;
    u[n-1] = 0;
  }
  else
  {
    qn = -0.5;
    u[n-1] = ( 3 / ( x[n-1] - x[n-2] ) ) * ( S_final - ( y[n-1] - y[n-2] ) / ( x[n-1] - x[n-2] ) );
  }

  //the last (n-1) point in S is calculated  
  S[n-1] = ( u[n-1] - qn * u[n-2] ) / ( qn * S[n-2] + 1 );
  
  // backsubstitution loop of the tridiagonal algorithm
  for ( int k = n - 2; k >= 0; k-- )
  {
    S[k] = S[k] * S[k+1] + u[k];
  }
  
}



double interpolation::spline_interpolation(const vector<double>& xa, const vector<double>& ya, const double x, const int index) const
{
  double a,b;
  unsigned int n = 0;
  
  if (index == -1 || index >= xa.size() )
  {
    n = bisection(xa, x);
  }
  else
  {
    n = index;
  }
  
  if ( S.size() == 0 )
  {
    cout << "FATAL ERROR in spline_interpolation: spline_construction needs to be called first" << endl;
  }
  
  double h = xa[n+1] - xa[n];
  if (h == 0.0)
  {
    cout << "ERROR in spline_interpolation: x-values must be distinct" << endl;
    return 1e-20;
  }
  
  a = ( xa[n+1] - x ) / h;
  b = ( x - xa[n] ) / h;
  
  double y = a * ya[n] + b * ya[n+1] + ( (pow(a,3.0) - a)*S[n] + (pow(b,3.0) - b)*S[n+1] ) * pow(h,2.0) / 6;
  
  return y;
}



int interpolation::bisection(const vector<double>& xa, const double x) const
{
  unsigned int n_low = 0;
  unsigned int n_high = xa.size() - 1;
  unsigned int n;
  
  while( (n_high - n_low) > 1 )
  {
    n = (n_high + n_low) / 2;
    if (xa[n] > x)
    {
      n_high = n;        
    }
    else
    {
      n_low = n;
    }
  }
  
  return n_low;
}
