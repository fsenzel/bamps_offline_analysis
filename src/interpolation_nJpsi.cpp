#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <math.h>

#include "interpolation_nJpsi.h"

using namespace std;

void interpolation_nJpsi::configure()
{
  // total numbers of tabulated points on the a, b "axes"
  n_i = 201; // T
  
  // spacings of tabulated values in a, b-direction
  delta_a = 0.002;

  // the smallest tabulated values of a, b, c and d
  a_start = 0.1;
  
  filename = "data/n_Jpsi_36_eq.dat";
  
  loadData();
}

// not properly normalized
double interpolation_nJpsi::getN(const double T) const
{
  double n, temp;
  
  temp = T;
  
  if( T < a_start)
  {
    temp = a_start;
    cout << "error in interpolation_nJpsi::getN, T=" << T << endl;
  }
  else if(T > a_start + delta_a*n_i )
  {
    temp = a_start + delta_a*n_i;
    cout << "error in interpolation_nJpsi::getN, T=" << T << endl;
  }
    
  n = getInterpolatedData( temp );
  
//   if( T >= a_start && T <= a_start + delta_a*n_i )  
//     n = getInterpolatedData( T );
//   else
//   {
//     cout << "error in interpolation_nJpsi::getN, T=" << T << endl;
//     n = 0;
//   }
  
  return n;
}

