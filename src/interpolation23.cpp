//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <math.h>

#include "interpolation23.h"
#include "interpolation.h"
#include "FPT_compare.h"

using namespace std;


interpolation23::interpolation23()
{
  minus_infinity = -9999;      // value with which "-inf" entries will be substituted
  
  // total numbers of tabulated points on the a, b,c and d-"axes"
  n_i = 250;
  n_j = 181;
  n_k_low = 18;
  n_k_high = 4;
  n_k = n_k_low + n_k_high;
  n_l = 21;
  
  // spacings of tabulated values in a, b, c and d-direction
  // in c-direction there are different spacings for small than for large c-values
  delta_a = 0.1;
  delta_b = 0.1;
  delta_c_low = 0.2;
  delta_c_high = 1.0;
  delta_d = 0.05;

  // the smallest tabulated values of a, b, c and d
  a_start = -15;
  b_start = -3.0;
  c_start = 0;
  d_start = 0;
  
  c_separator = c_start + (n_k_low-1) * delta_c_low;
  
  expected_entries = n_i * n_j * n_k * n_l;
  //allocate memory based on the known size (# of entries) of the table
  //a matter of efficieny - not necessary
  I23data.reserve( expected_entries + 1 );

  if ( readTables() == I23_SUCCESS )
  {
    cout << "I23 readout successful" << endl;
  }
}


interpolation23::~interpolation23()
{
  I23data.clear();
}



I23_ERROR_CODE interpolation23::readTables()
{
  string filename = "data/I23_table.dat";   //data file holding the tabulated values, see below for required structure
  fstream file( filename.c_str(), ios::in );

  if ( file )
  {
    string line;
    stringstream inStream;
    double content;

    do
    {
      getline( file, line );
      if ( line.find( "#", 0 ) == string::npos && line.length() != 0 )    //ignore empty lines and lines with leading "#"
      {
        if ( line.find( "-inf", 0 ) == string::npos )  //the string "-inf" can't be pushed into a vector<double>
        {
          inStream.str( "" );
          inStream.clear();

          inStream.str( line );
          inStream >> content;

          if ( !inStream.fail() )
          {
            I23data.push_back( content );
          }
          else
          {
            cout << "severe error in interpolation23::readTables() - ill-formated data file" << endl;
            return I23_READOUT_ERROR;
          }
        }
        else     //i.e. when string "-inf" is found
        {
          I23data.push_back( minus_infinity );  //push_back(minus_infinity) instead of "-inf"
        }

      }
    }
    while ( !file.eof() );
  }
  else
  {
    // throw an exception if the file could not have been read
    throw eI23_read_error( "error in interpolation23::readTables() - could not open " + filename );
    return I23_FILE_ERROR;
  }

  // The number of tabulated values read from file must be equal to the expected number, otherwise
  // the interpolation routine will fail to work correctly!
  if ( I23data.size() != expected_entries )
  {
    cout << I23data.size() << " instead of " << expected_entries << "entries" << endl;
    // throw an exception if the number of entries in the file is not right
    throw eI23_read_error( "error in interpolation23::interpolation23() - unexpected number of entries in I23 table" );
    return I23_SIZE_ERROR;
  }
  else
    return I23_SUCCESS;
}



/**
 * 4-dimensional interpolation / extrapolation of the tabulated ln(I23)-values at a given point (a,b,c,d),
 * with a=ln(md2), b=ln(lambda), c=beta and d=theta.
 * This routine serves as a wrapper to different methods of interpolation (linear, polynomial, cubic splines). 
 *
 * @param[in] a ln(md2_scaled)
 * @param[in] b ln(lambda_scaled)
 * @param[in] c beta
 * @param[in] d cos(theta)
 * @return I23 = exp(f), the value of I23 at (a,b,c,d) found via interpolation from the tabulated values of f
 */
double interpolation23::getI23( const double a, const double b, const double ca, const double d ) const
{
  double I23;
  const double c = log ( 1 / sqrt( 1 - pow(ca,2.0) ) );
  if ( !inRange( a, b, c, d ) )
  {
    //cout << "not in range" << endl;
    return 0;
  }
  else
  {
    //if ( c < 2.0 && d < 0.85 && a > -3.0 && a < 5.0 )
    if ( true )
    {
      I23 = getI23_linearInterpolation( a, b, c, d );    
//       if ( I23 < 1 ) 
//       {
//         I23 = getI23_splineInterpolation( a, b, c, d );
//       }
    }
    else
    {
      I23 = getI23_splineInterpolation( a, b, c, d );
    }
    
    
    if (I23 < 0.1)
      return 0;
    else
      return I23;
  }
}



/**
 * 4-dimensional linear interpolation of the tabulated ln(I23)-values at a given point (a,b,c,d),
 * with a=ln(md2), b=ln(lambda), c=ln(gamma)=ln(1/sqrt(1-beta^2)) and d=cos(theta).
 * Even though purely linear interpolation can also be achieved by setting the order N=2 in the routine for polynomial interpolation,
 * this one might be a bit faster due to less overhead.
 * This routine is closely tailored to the known structure of the tabulated data (stored in interpolation23::I23data).
 * See below for a short documentation of the required data structure.
 *
 * @param[in] a ln(md2_scaled)
 * @param[in] b ln(lambda_scaled)
 * @param[in] c ln(gamma) = ln( 1/sqrt(1-beta^2) )
 * @param[in] d cos(theta)
 * @return I23 = exp(f), the value of I23 at (a,b,c,d) found via interpolation from the tabulated values of f
 */
double interpolation23::getI23_linearInterpolation( const double a, const double b, const double c, const double d ) const
{
  double delta_a_local = delta_a;
  // find indices i,j,k and l such that
  // a[i] <= a < a[i+1], b[j] <= b < b[j+1], c[k] <= c < c[k+1] and d[l] <= d < d[l+1]
  const int i = int(( a - a_start ) / delta_a );
  const int j = int(( b - b_start ) / delta_b );
  int k_low, k_high;
  if ( c < c_separator )
  {
    k_low = int(( c - c_start ) / delta_c_low );
    k_high = 0; 
  } else
  {
    k_low = n_k_low - 1;
    k_high = int(( c - c_separator ) / delta_c_high );
  }
  int l = int(( d - d_start ) / delta_d );

  //acount for the possibility cos(theta) = 1
  if ( l == n_l-1 )
  {
    l -= 1;
  }

  if ( i < 0 || j < 0 || j >= n_j || k_low < 0 || k_high < 0 || l < 0 || l >= n_l )
  {
    cout << "unexpected - (i,j,k,l) not in range" << endl;
    cout << n_i << "  " << n_j << "  " << n_k << "  " << n_l << endl;
    cout << i << "  " << j << "  " << k_low + k_high << "  " << l << endl;
    cout << a << "  " << b << "  " << c << "  " << d << endl;

    return 0;
  }

  //cout << " i = " << i << "   j = " << j << "   k = " << k << "   l = " << l << endl;

  // find a[i], b[j], c[k] and d[l]
  double a_i = a_start + i * delta_a;
  const double b_j = b_start + j * delta_b;
  const double c_j = c_start + k_low * delta_c_low + k_high * delta_c_high;
  const double d_l = d_start + l * delta_d;

  //cout << " a_i = " << a_i << "   b_j = " << b_j << "   c_k = " << c_k << "   d_l = " << d_l << endl;

  // temporary vector objects used to initialize the multi-dimensional arrays based on STL vectors
  vector<double> t1( 2, 0 );
  vector< vector<double> > t2( 2, t1 );
  vector< vector< vector<double> > > t3( 2, t2 );

  vector< vector< vector< vector<double> > > > f( 2, t3 );

  // a = ln(md2_scaled)=ln(md2/s) might be out of range due to small s
  // in this case the tabulated values are extrapolated linearly in a direction
  if ( i > (n_i - 1) )
  {
    //cout << "extrapolation" << endl;
    int delta_ex = 5;

    for ( int p = 0; p <= 1; p++ )
      for ( int q = 0; q <= 1; q++ )
        for ( int r = 0; r <= 1; r++ )
          for ( int s = 0; s <= 1; s++ )
            f[p][q][r][s] = I23data[ get_index( (n_i - 1) + (( p-1 )*delta_ex ), j+q, (k_low+k_high)+r, l+s )];

    a_i = a_start + ( (n_i - 1) - delta_ex ) * delta_a;
    delta_a_local = delta_a * delta_ex;
  }
  else  // the normal case for all values in range
  {
    // The 4-dimensional array f[p][q][r][s], where p,q,r,s = {0,1}, holds the tabulated values
    // at the corners of the 4-dim "cube" surrounding the point (a,b,c,d).
    // f[0][0][0][0] corresponds to the tabulated value at (a[i],b[j],c[k],d[l]),
    // whereas f[1][1][1][1] corresponds to the value at (a[i+1],b[j+1],c[k+1],d[l+1]).
    for ( int p = 0; p <= 1; p++ )
      for ( int q = 0; q <= 1; q++ )
        for ( int r = 0; r <= 1; r++ )
          for ( int s = 0; s <= 1; s++ )
            f[p][q][r][s] = I23data[ get_index( i+p, j+q, (k_low+k_high)+r, l+s )];
  }

  // f_a[q][r][s] holds the values of f[p][q][r][s] linearly interpolated in a-direction.
  // Similar to the notation above, f_a[0][0][0] corresponds to the value at (a,b[j],c[k],d[l])
  // and f_a[1][1][1] corresponds to the value at (a,b[j+1],c[k+1],d[l+1]).
  // Note that these points are at a not at a[i+p]!
  vector< vector< vector<double> > > f_a( 2, t2 );
  for ( int q = 0; q <= 1; q++ )
    for ( int r = 0; r <= 1; r++ )
      for ( int s = 0; s <= 1; s++ )
        f_a[q][r][s] = f[0][q][r][s] + ( a - a_i ) / delta_a_local * ( f[1][q][r][s] - f[0][q][r][s] );

      // f_ab[r][s] holds the values of f_a[q][r][s] linearly interpolated in b-direction.
  vector< vector<double> > f_ab( 2, t1 );
  for ( int r = 0; r <= 1; r++ )
    for ( int s = 0; s <= 1; s++ )
      f_ab[r][s] = f_a[0][r][s] + ( b - b_j ) / delta_b * ( f_a[1][r][s] - f_a[0][r][s] );


  // f_abc[s] holds the values of f_ab[r][s] linearly interpolated in c-direction.
  vector<double> f_abc( 2, 0 );
  for ( int s = 0; s <= 1; s++ )
  {
    if ( c < c_separator )
    {
      f_abc[s] = f_ab[0][s] + ( c - c_j ) / delta_c_low * ( f_ab[1][s] - f_ab[0][s] );
    }
    else
    {
      f_abc[s] = f_ab[0][s] + ( c - c_j ) / delta_c_high * ( f_ab[1][s] - f_ab[0][s] );
    }      
  }

  // f_abcd is the value of f_abc[s] linearly interpolated in d-direction,
  // thus the final result!
  double f_abcd = f_abc[0] + ( d - d_l ) / delta_d * ( f_abc[1] - f_abc[0] );

  return exp( f_abcd );
}





/**
 * This routine uses a polynomial interpolation / extrapolation instead of pure linear interpolation. The degree in each direction is 
 * controlled by the variables ord_? set early in the routine
 *
 * @param[in] a ln(md2_scaled)
 * @param[in] b ln(lambda_scaled)
 * @param[in] c ln(gamma) = ln( 1/sqrt(1-beta^2) )
 * @param[in] d cos(theta)
 * @return I23 = exp(f), the value of I23 at (a,b,c,d) found via extrapolation from the tabulated values of f
 */
double interpolation23::getI23_polyInterpolation( const double a, const double b, const double c, const double d ) const
{
  // number of input points for extrapolation in ?-direction, the degree of the polynomial used is N-1
  const int ord_i = 4;
  const int ord_j = 4;
  const int ord_k = 4;  
  const int ord_l = 4;

  //---------- variables needed for the use of the polynomial extrapolation routines
  double y, dy;
  
  // find indices i,j,k and l such that
  // a[i] <= a < a[i+1], b[j] <= b < b[j+1], c[k] <= c < c[k+1] and d[l] <= d < d[l+1]
  int i = 0, j = 0, k = 0, l = 0;
     
  //--------- a-direction --------------
  if ( a < a_start )
  {
    i = -1;  //to catch error
  }
  else if ( a > (a_start + (n_i - 1) * delta_a) )
  {
    i = n_i - ord_i;
  }
  else
  {
    i = (( a - a_start ) / delta_a );
    i = i - (ord_i / 2) + 1;
  }

  if ( i < 0 )
  {
    i = 0;
  }
  if ( i + ord_i > n_i )
  {
    i = n_i - ord_i; 
  }
  //--------- a-direction --------------
  
  //--------- b-direction --------------
  if ( b < b_start )
  {
    j = -1;  //to catch error
  }
  else if ( b > (b_start + (n_j - 1) * delta_b) )
  {
    j = n_j - ord_j;
  }
  else
  {
    j = (( b - b_start ) / delta_b );
    j = j - (ord_j / 2) + 1;
  }

  if ( j < 0 )
  {
    j = 0;
  }
  if ( j + ord_j > n_j )
  {
    j = n_j - ord_j; 
  }
  //--------- b-direction --------------
  
  //--------- c-direction --------------
  if ( c < c_start )
  {
    k = -1;  //to catch error
  }
  else if ( c > (c_separator + (n_k_high - 1) * delta_c_high) )
  {
    k = n_k - ord_k;
  }
  else
  {
    if ( c < c_separator )
    {
      k = int(( c - c_start ) / delta_c_low );
    } 
    else
    {
      k = n_k_low - 1 + int(( c - c_separator ) / delta_c_high );
    }
    k = k - (ord_k / 2) + 1;
  }

  if ( k < 0 )
  {
    k = 0;
  }
  if ( k + ord_k > (n_k) )
  {
    k = n_k - ord_k; 
  }
  //--------- c-direction --------------
  
  //--------- d-direction --------------
  if ( d < d_start )
  {
    l = -1;  //to catch error
  }
  else if ( d > (d_start + (n_l - 1) * delta_d) )
  {
    l = n_l - ord_l;
  }
  else
  {
    l = (( d - d_start ) / delta_d );
    l = l - (ord_l / 2) + 1;
  }

  if ( l < 0 )
  {
    l = 0;
  }
  if ( l + ord_l > n_l )
  {
    l = n_l - ord_l; 
  }
  //--------- d-direction --------------
  
  if ( i < 0 || j < 0 || j > n_j || k < 0 || k > (n_k) || l < 0 || l > n_l )
  {
    cout << "unexpected - (i,j,k,l) not in range" << endl;
    cout << i << "  " << j << "  " << k << "  " << l << endl;
    cout << a << "  " << b << "  " << c << "  " << d << endl;

    return 0;
  }

  //cout << " i = " << i << "   j = " << j << "   k = " << k << "   l = " << l << endl;

  
  //vectors used for interpolation / unit-offset !!!
  vector<double> xa( ord_k + 1, 0.0);
  vector<double> ya( ord_k + 1, 0.0);
  
  // temporary vector objects used to initialize the multi-dimensional arrays based on STL vectors
  vector<double> t1( ord_l, 0 );
  vector< vector<double> > t2( ord_k, t1 );  // N points in c-direction needed
  vector< vector< vector<double> > > t3( ord_j, t2 );

  vector< vector< vector< vector<double> > > > f( ord_i, t3 );
  vector< vector< vector< vector<double> > > > c_values_matrix( ord_i, t3 );
  
  for ( int p = 0; p < ord_i; p++ )
  {
    for ( int q = 0; q < ord_j; q++ )
    {
      for ( int s = 0; s < ord_l; s++ )
      {
        
        int first_r_inf = ord_k;
        
        for ( int r = ord_k - 1; r >= 0; r-- )
        {
          f[p][q][r][s] = I23data[ get_index( i+p, j+q, k+r, l+s )];
          if ( (k + r) < n_k_low )
          {
            c_values_matrix[p][q][r][s] = c_start + (k + r) * delta_c_low;
          }
          else
          {
            c_values_matrix[p][q][r][s] = c_separator + (k + r - n_k_low + 1) * delta_c_high;
          }
          
          if ( FPT_COMP_E(f[p][q][r][s],-9999) )
          {
            first_r_inf = r;
          }
        }
        
        
        if ( first_r_inf <= (ord_k/2 - 1) )
        {
          for (int r = 0; r < ord_k; r++)
          {
            f[p][q][r][s] = -9999;
          }
        }        
        else if ( first_r_inf < ord_k)
        {
          for (int r = 1; r <= first_r_inf; r++)
          {
            f[p][q][ord_k - r][s] = f[p][q][first_r_inf - r][s];
            if ( (k + first_r_inf - r) < n_k_low )
            {
              c_values_matrix[p][q][ord_k - r][s] = c_start + (k + first_r_inf - r) * delta_c_low;
            }
            else
            {
              c_values_matrix[p][q][ord_k - r][s] = c_separator + (k + first_r_inf - r - n_k_low + 1) * delta_c_high;
            }
          }
          for (int r = 0; r < ord_k - first_r_inf; r++)
          {
            f[p][q][r][s] = I23data[ get_index( i+p, j+q, k - (ord_k - first_r_inf) + r, l+s )];
            if ( (k - (ord_k - first_r_inf) + r) < n_k_low )
            {
              c_values_matrix[p][q][r][s] = c_start + (k - (ord_k - first_r_inf) + r) * delta_c_low;
            }
            else
            {
              c_values_matrix[p][q][r][s] = c_separator + (k - (ord_k - first_r_inf) + r - n_k_low + 1) * delta_c_high;
            }
          }
        }
        
           
      }
    }
  }

  vector< vector<double> > t4(ord_j,t1);
  vector< vector< vector<double> > > f_c( ord_i, t4 );
    xa.resize( ord_k + 1, 0.0 );
  ya.resize( ord_k + 1, 0.0 );
  for ( int p = 0; p < ord_i; p++ )
  {
    for ( int q = 0; q < ord_j; q++ )
    {
      for ( int s = 0; s < ord_l; s++ )
      {
        
        for ( int kk = 1; kk < ord_k + 1; kk++ )
        {
          xa[kk] = c_values_matrix[p][q][kk-1][s];
          ya[kk] = f[p][q][kk-1][s];
        }
        
        polynomialInterpolation( xa, ya, ord_k, c, y, dy );
        f_c[p][q][s] = y;
      }
    }
  }

  
  // f_ab[r][s] holds the values of f_a[q][r][s] linearly interpolated in b-direction.
  vector< vector<double> > f_ca( ord_j, t1 );
  xa.resize( ord_i + 1, 0.0 );
  ya.resize( ord_i + 1, 0.0 );
  for ( int ii = 1; ii < ord_i + 1; ii++ )   //write the x-values for the polynomial extrapolation
  {
    xa[ii] = a_start + (i + ii - 1) * delta_a;
  }
  for ( int q = 0; q < ord_j; q++ )
  {
    for ( int s = 0; s < ord_l; s++ )
    {
      
      for ( int ii = 1; ii < ord_i + 1; ii++ )
      {
        ya[ii] = f_c[ii-1][q][s];
      }
      polynomialInterpolation( xa, ya, ord_i, a, y, dy );
      f_ca[q][s] = y;
    }
  }


  // f_abc[s] holds the values of f_ab[r][s] polynomially extrapolated in c-direction.
  vector<double> f_cab( ord_l, 0 );
  xa.resize( ord_j + 1, 0.0 );
  ya.resize( ord_j + 1, 0.0 );
  for ( int jj = 1; jj < ord_j + 1; jj++ )   //write the x-values for the polynomial extrapolation
  {
    xa[jj] = b_start + (j + jj - 1) * delta_b;
  }
  for ( int s = 0; s < ord_l; s++ )
  {
    
    for ( int jj = 1; jj < ord_j + 1; jj++ )
    {
      ya[jj] = f_ca[jj-1][s];
    }

    polynomialInterpolation( xa, ya, ord_j, b, y, dy );

    f_cab[s] = y;
  }

  // f_abcd is the value of f_abc[s] linearly interpolated in d-direction,
  // thus the final result!
  xa.resize( ord_l + 1, 0.0 );
  ya.resize( ord_l + 1, 0.0 );
  for ( int ll = 1; ll < ord_l + 1; ll++ )   //write the x-values (beta) for the polynomial extrapolation
  {
    xa[ll] = d_start + (l + ll - 1) * delta_d;
    ya[ll] = f_cab[ll-1];
  }
  polynomialInterpolation( xa, ya, ord_l, d, y, dy );
  
  double f_cabd = y;

  return exp( f_cabd );
}



/**
 * This routine uses cubic splines for interpolation / extrapolation of the value for I23 at given a,b,c and d.
 *
 * @param[in] a ln(md2_scaled)
 * @param[in] b ln(lambda_scaled)
 * @param[in] c ln(gamma) = ln( 1/sqrt(1-beta^2) )
 * @param[in] d cos(theta)
 * @return I23 = exp(f), the value of I23 at (a,b,c,d) found via inter-/ extrapolation from the tabulated values of f
 */
double interpolation23::getI23_splineInterpolation( const double a, const double b, const double c, const double d ) const
{
  // number of input points for extrapolation in ?-direction
  const int ord_i = 4;
  const int ord_j = 4;
  const int ord_k = 4;  
  const int ord_l = 4;
  
  // the indices i,j,k and l are set such that e.g. the range a[ i, i+ord_i-1 ] contains the value a etc.
  int i, j, k, l;  
  // the indices ?_location are set such that e.g. a[i_location] <= a < a[i_location+1], i_location must be within [i,i+ord_i-2]
  int i_location, j_location, k_location, l_location;

  //>>>------ find indices in a-direction --------------
  if ( a < a_start )   //to catch error (should never happen)
  {
    i = -1; 
  }
  else if ( a > (a_start + (n_i - 1) * delta_a) )  // when a is out of the tabulated range
  {
    i = n_i - ord_i;
    i_location = n_i - 1;
  }
  else  // the standard case
  {
    i_location = (( a - a_start ) / delta_a );  // find i_location such that a[i_location] <= a < a[i_location+1]
    i = i_location - (ord_i / 2) + 1;           // extend range to the left of i_location
  }
  
  if ( i < 0 ) // in case the range has been extended to the left too much
  {
    i = 0;
  }
  if ( i + ord_i > n_i )  // in case the range has been extended to the left too little
  {
    i = n_i - ord_i; 
  }
  //--------- find indices in a-direction -----------<<<
  
  //>>>------ find indices in b-direction --------------
  if ( b < b_start )   //to catch error (should never happen)
  {
    j = -1;
  }
  else if ( b > (b_start + (n_j - 1) * delta_b) )  // when b is out of the tabulated range
  {
    j = n_j - ord_j;
    j_location = n_j - 1;
  }
  else  // the standard case
  {
    j_location = (( b - b_start ) / delta_b );   // find j_location such that b[j_location] <= b < b[j_location+1]               
    j = j_location - (ord_j / 2) + 1;            // extend range to the left of j_location
  }

  if ( j < 0 )        // in case the range has been extended to the left too much and has become smaller than 0
  {
    j = 0;
  }
  if ( j + ord_j > n_j ) // in case the range has been extended to the left too little and has become larger than the final index
  {
    j = n_j - ord_j; 
  }
  //--------- find indices in b-direction -----------<<<
  
  //>>>------ find indices in c-direction --------------
  if ( c < c_start )   //to catch error (should never happen)
  {
    k = -1; 
  }
  else if ( c > (c_separator + (n_k_high - 1) * delta_c_high) )  // when c is out of the tabulated range
  {
    k = (n_k) - ord_k;
    k_location = n_k - 1;
  }
  else  // the standard case
  {
    // find k_location such that c[k_location] <= c < c[k_location+1]
    if ( c < c_separator )
    {
      k_location = int(( c - c_start ) / delta_c_low );
    } 
    else
    {
      k_location = n_k_low - 1 + int(( c - c_separator ) / delta_c_high );
    }
    k = k_location - (ord_k / 2) + 1;           // extend range to the left of k_location
  }

  if ( k < 0 )   // in case the range has been extended to the left too much and has become smaller than 0
  {
    k = 0;
  }
  if ( k + ord_k > n_k )  // in case the range has been extended to the left too little and has become larger than the final index
  {
    k = n_k - ord_k; 
  }
  //--------- find indices in c-direction -----------<<<
  
  //>>>------ find indices in d-direction --------------
  if ( d < d_start )  //to catch error (should never happen)
  {
    l = -1;
  }
  else if ( d > (d_start + (n_l - 1) * delta_d) )   // when d is out of the tabulated range
  {
    l = n_l - ord_l;
    l_location = n_l - 1;
  }
  else   // the standard case
  {
    l_location = (( d - d_start ) / delta_d );  // find k_location such that c[k_location] <= c < c[k_location+1]
    l = l_location - (ord_l / 2) + 1;           // extend range to the left of k_location
  }

  if ( l < 0 )    // in case the range has been extended to the left too much and has become smaller than 0
  {
    l = 0;
  }
  if ( l + ord_l > n_l )   // in case the range has been extended to the left too little and has become larger than the final index
  {
    l = n_l - ord_l; 
  }
  //--------- find indices in d-direction -----------<<<
  
  
  // temporary vector objects used to initialize the multi-dimensional arrays based on STL vectors
  vector<double> t1( ord_k, 0.0 );
  vector< vector<double> > t2( ord_j, t1 );
  vector< vector< vector<double> > > f_d( ord_i, t2 ); // 3-dim matrix that will contain the results interpolated in d-direction
  
  vector<double> t3( ord_j, 0.0 );
  vector< vector<double> > f_dc( ord_i, t3);           // 2-dim matrix that will contain the results interpolated in d- and c-direction
  
  vector<double> f_dcb( ord_i, 0.0 );                  // 1-dim vector that will contain the results interpolated in d-, c- and b-direction
  double f_dcba = 0.0;                                 // double that will contain the result interpolated in all directions, a, b, c and d
  
  vector<double> x_values(ord_l, 0.0);   // vector that will contain the x-values for the current interpolation steps
  vector<double> y_values(ord_l, 0.0);   // vector that will contain the y-values for the current interpolation steps
  
  // in case the selected subrange contains the value -9999 (used to represent -infinity, i.e. I23 = exp(-inf) = 0) the subrange needs
  // to be shifted adequately. i_shift etc. will contain the magnitude of such a shift.
  int i_shift, j_shift, k_shift, l_shift;
  // in case the shift denoted by i_shift etc. would bring the subrange (of length ord_i) out of the tabulated range, the subrange needs
  // to be resized (i.e. shrinked). i_values_length etc. denotes the current size of a subrange  
  int i_values_length, j_values_length, k_values_length, l_values_length;
  int k_repetition = 0, l_repetition = 0;
  
  
  /* The interpolation is implemented as follows:
  In principle the 4-dimensional interpolation is done by successive interpolation 1-dim interpolations in one direction. First
  the values are interpolated in d-direction at every point in a,b and c, the so obtained results are then interpolated in c-direction
  and so on.
  This is implemented via three nested loops, the innermost being responsible for the d-direction. Enclosed by c and b-loops. After the
  nested loops the remaining interpolation in a-direction is carried out.
  Each interpolation itself is performed by passing a subset of values to routines doing cubic spline interpolation.
  
  Some specific cases require special attention. 
  In general, when the selected subset in a given direction contains values of -inf (-9999) at its end or beginning, the 
  interpolation should be redone using a shifted subset. For the two innermost interpolations (d and c-direction) this is done
  by enclosing them in a do-while loop and redoing the whole interpolation in that direction with shifted subsets if necessary.
  For the two outer interpolations (b and a-direction) this effort is avoided, instead the subset is simply reduced in since such
  that it does not contain the -inf values any more.      
  */
  
  // loop for p, i.e. a-direction
  for ( int p = 0; p < ord_i; p++ )
  {

    // loop for q, i.e. b-direction 
    for ( int q = 0; q < ord_j; q++ )
    {
      k_shift = 0;
      k_repetition = 0;
      do   // do-while-loop to eventually re-do interpolation in c-direction when the given subset contains values of -inf
      {
        
        //loop for r, i.e. c-direction
        for ( int r = 0; r < ord_k; r++ )
        {
          //>>>--------------- interpolation in d-direction ------------------
          
          // check whether d[i] <= d < d[i+1] such that f(d[i]) = -inf and f(d[i+1]) = -inf, then f(d) = -inf is assumed
          if( I23data[ get_index( i+p, j+q, k+k_shift+r, l_location )] < -9000 
              && ( l_location == (n_l - 1) || I23data[ get_index( i+p, j+q, k+k_shift+r, l_location + 1)] < -9000 ) )
          {
            f_d[p][q][r] = -9999;
          }
          else  // the general case
          { 
            l_shift = 0;
            l_repetition = 0;
            do  // do-while-loop to eventually re-do interpolation in d-direction when the given subset contains values of -inf
            {                                
              l_values_length = ord_l;
              
              // in case the shift in l-direction is such that the new subrange would exceed the range of tabulated values, the length of
              // the subrange (l_values_length) is adjusted
              if ( l + l_shift + ord_l > n_l ) 
              {
                l_values_length = n_l - (l + l_shift);
                x_values.resize(l_values_length, 0.0);
                y_values.resize(l_values_length, 0.0);
              }
              // adjustment of l_values_length in case the shift brings the beginning of the subrange to negative indices
              else if ( l + l_shift < 0 )
              {
                l_values_length = ord_l + (l + l_shift);
                l_shift = -l;
                x_values.resize(l_values_length, 0.0);
                y_values.resize(l_values_length, 0.0); 
              }
              
              //write the x- and y-values (subrange) used in the interpolation / extrapolation
              for ( int s = 0; s < l_values_length; s++ ) 
              {
                x_values[s] = d_start + (l + l_shift + s) * delta_d;
                y_values[s] = I23data[ get_index( i+p, j+q, k+k_shift+r, l+l_shift+s )];
              }
              
              if (l_repetition > 2)
              {
                int temp_l_shift = getShift(y_values);
                
                if ( temp_l_shift != 0 )
                {
                  l_values_length = ord_l - abs(temp_l_shift);
                  x_values.resize(l_values_length, 0.0);
                  y_values.resize(l_values_length, 0.0);
                  
                  if(temp_l_shift < 0)
                    temp_l_shift = 0;
        
                  for ( int s = 0; s < l_values_length; s++ ) 
                  {
                    x_values[s] = d_start + (l + l_shift + temp_l_shift + s) * delta_d;
                    y_values[s] = I23data[ get_index( i+p, j+q, k+k_shift+r, l+l_shift+temp_l_shift+s )];
                  } 
                }
              }
            
              l_shift = getShift(y_values);  // determine whether a shift is necessary
              ++l_repetition; 
            } while(l_shift != 0);
          
            // in case there is only one entry left in the subrange, this value is taken as an approximation to f(d)
            if (y_values.size() == 1)
            {
              f_d[p][q][r] = y_values[0];
            }
            else  // in general a cubic spline interpolation is performed
            {
              interpolation c_spline_int(x_values, y_values);
              f_d[p][q][r] = c_spline_int.spline_interpolation(x_values, y_values, d);
            }
          }
          //--------------- interpolation in d-direction ------------------<<<        
        
        } //end of loop for r, i.e. c-direction
      
        //>>>--------------- interpolation in c-direction ------------------
        x_values.resize(ord_k, 0.0);
        y_values.resize(ord_k, 0.0);
      
        // check whether c[k] <= c < c[k+1] such that f(c[k]) = -inf and f(c[k+1]) = -inf, then f(c) = -inf is assumed
        if( f_d[p][q][k_location-(k+k_shift)] < -9000 && ( k_location == (n_k - 1) || f_d[p][q][k_location-(k+k_shift)+1] < -9000 ) )
        {
          f_dc[p][q] = -9999;
        }
        else
        {
          k_values_length = ord_k;
          
            // in case the shift in k-direction is such that the new subrange would exceed the range of tabulated values, the length of
            // the subrange (k_values_length) is adjusted
          if ( k + k_shift + ord_k > n_k )
          {
            k_values_length = n_k - (k + k_shift);
            x_values.resize(k_values_length, 0.0);
            y_values.resize(k_values_length, 0.0);
          }
            // adjustment of l_values_length in case the shift brings the beginning of the subrange to negative indices
          else if ( k + k_shift < 0 )
          {
            k_values_length = ord_k + (k + k_shift);
            k_shift = -k;
            x_values.resize(k_values_length, 0.0);
            y_values.resize(k_values_length, 0.0); 
          }
          
                                      
          //write the x- and y-values (subrange) used in the interpolation / extrapolation
          for ( int r = 0; r < k_values_length; r++ ) 
          {
            if ( (k + k_shift + r) < n_k_low )
            {
              x_values[r] = c_start + (k + k_shift + r) * delta_c_low;
            }
            else
            {
              x_values[r] = c_separator + (k + k_shift + r - n_k_low + 1) * delta_c_high;
            }
            
            y_values[r] = f_d[p][q][r];
          }
          
          if (k_repetition > 2)
          {
            int temp_k_shift = getShift(y_values);
                
            if ( temp_k_shift != 0 )
            {
              k_values_length = ord_k - abs(temp_k_shift);
              x_values.resize(k_values_length, 0.0);
              y_values.resize(k_values_length, 0.0);
                  
              if(temp_k_shift < 0)
                temp_k_shift = 0;
        
              for ( int r = 0; r < k_values_length; r++ ) 
              {
                if ( (k + k_shift + r) < n_k_low )
                {
                  x_values[r] = c_start + (k + k_shift + temp_k_shift + r) * delta_c_low;
                }
                else
                {
                  x_values[r] = c_separator + (k + k_shift + temp_k_shift + r - n_k_low + 1) * delta_c_high;
                }
                
                y_values[r] = f_d[p][q][r+temp_k_shift];
              } 
            }
          }
          
          k_shift = getShift(y_values);   // determine whether a shift is necessary
              
          // in case there is only one entry left in the subrange, this value is taken as an approximation to f(c)
          if (y_values.size() == 1)
          {
            f_dc[p][q] = y_values[0];
          }
          else  // in general a cubic spline interpolation is performed
          {
            interpolation c_spline_int(x_values, y_values);
            f_dc[p][q] = c_spline_int.spline_interpolation(x_values, y_values, c);
          }
        }
      //--------------- interpolation in c-direction ------------------<<<
        ++k_repetition;
      } while(k_shift != 0);  // this while-loop encloses the whole interpolation in d-direction!!!
      
    } //end of loop for q, i.e. b-direction 
    
    //>>>--------------- interpolation in b-direction ------------------
    x_values.resize(ord_j, 0.0);
    y_values.resize(ord_j, 0.0);
      
    // check whether b[j] <= b < b[j+1] such that f(b[j]) = -inf and f(b[j+1]) = -inf, then f(b) = -inf is assumed
    if( f_dc[p][j_location-j] < -9000 && ( j_location == (n_j - 1) || f_dc[p][j_location-j+1] < -9000 ) )
    {
      f_dcb[p] = -9999;
    }
    else
    {    
      j_shift = 0;
      j_values_length = ord_j;
      
      // This do-while-loop does NOT repeat the previous interpolation if -inf-values are spotted in the selected subrange.
      // Only the subrange is shortened.
      do 
      {
        if(j_shift < 0)
          j_shift = 0;
        
        // here the subrange is selected first...
        for ( int q = 0; q < j_values_length; q++ ) 
        {
          x_values[q] = b_start + (j + j_shift + q) * delta_b;
          y_values[q] = f_dc[p][q+j_shift];
        }
            
        //... then a possible shift is determined
        j_shift = getShift(y_values);
      
        // j_values_length is reduced by j_shift; abs() because getShift can return positive and negative values
        if ( j_shift != 0 )
        {
          j_values_length = ord_j - abs(j_shift);
          x_values.resize(j_values_length, 0.0);
          y_values.resize(j_values_length, 0.0);
        }

      } while (j_shift != 0);
      
      // in case there is only one entry left in the subrange, this value is taken as an approximation to f(b)
      if (y_values.size() == 1)
      {
        f_dcb[p] = y_values[0];  
      }
      else  // in general a cubic spline interpolation is performed
      {
        interpolation c_spline_int(x_values, y_values);
        f_dcb[p] = c_spline_int.spline_interpolation(x_values, y_values, b);
      }
    }
    //--------------- interpolation in b-direction ------------------<<<

  } //end of loop for p, i.e. a-direction
  
  //>>>--------------- interpolation in a-direction ------------------    
  // check whether b[j] <= b < b[j+1] such that f(b[j]) = -inf and f(b[j+1]) = -inf, then f(b) = -inf is assumed
  if( f_dcb[i_location-i] < -9000 && ( i_location == (n_i - 1) || f_dcb[i_location-i+1] < -9000 ) )
  {
    f_dcba = -9999;
  }
  else
  {
    x_values.resize(ord_i, 0.0);
    y_values.resize(ord_i, 0.0);
            
    i_values_length = ord_i;
    i_shift = 0;
    
    // This do-while-loop does NOT repeat the previous interpolation if -inf-values are spotted in the selected subrange.
    // Only the subrange is shortened.
    do
    {
      if(i_shift < 0)
        i_shift = 0;
      
      // here the subrange is selected first...
      for ( int p = 0; p < i_values_length; p++ ) 
      {
        x_values[p] = a_start + (i + i_shift + p) * delta_a;
        y_values[p] = f_dcb[i_shift + p];
      }
            
      //... then a possible shift is determined
      i_shift = getShift(y_values);
      
      // i_values_length is reduced by i_shift; abs() because getShift can return positive and negative values
      if ( i_shift != 0 )
      {
        i_values_length = ord_i - abs(i_shift);
        x_values.resize(i_values_length, 0.0);
        y_values.resize(i_values_length, 0.0);
      }

    } while(i_shift != 0);
    
    // in case there is only one entry left in the subrange, this value is taken as an approximation to f(a)
    if (y_values.size() == 1)
    {
      f_dcba = y_values[0];
    }
    else  // in general a cubic spline interpolation is performed
    {
      interpolation c_spline_int(x_values, y_values);
      f_dcba = c_spline_int.spline_interpolation(x_values, y_values, a);
    }
  }
  //--------------- interpolation in a-direction ------------------<<<

  
  return exp( f_dcba );
}






int interpolation23::get_index( const int i, const int j, const int k, const int l ) const
{
  /*
  Utility routine to get the index number (line number - 1) of the point P=(a[i],b[j],c[k],d[l]) such that
  I23data[index] is the tabulated value corresponding to the point P.
  */

  return (( i * n_j * n_k * n_l ) + ( j * n_k * n_l ) + ( k * n_l ) + l );
  //add 1 to get the line number of the desired value in the data file (I23data vector has offset 0, line numbers have offset 1)
}



bool interpolation23::inRange( const double ln_md2, const double ln_lambda, const double beta, const double cos_theta ) const
{
  if ( ln_lambda <= -3.0 )
  {
    return false;
  }
  else
  {
    return true;
  }
}


int interpolation23::getShift(const vector<double>& values) const
{
  int shift = 0;
  
  if( values[0] < -9000 )
  {
    if ( values[values.size()-1] < -9000 )
      return 0;
    
    do
    {
      shift++;
    } while( shift < values.size() && values[shift] < -9000 );
  }
  else if ( values[values.size()-1] < -9000 )
  {
    do
    {
      shift--;
    } while( shift > -values.size() && values[values.size() - 1 + shift] < -9000 );
  }
  
  return shift;
}


bool interpolation23::canInterpolate(const double a, const double b, const double ca, const double d) const
{
  const double c = log ( 1 / sqrt( 1 - pow(ca,2.0) ) );
  
  if ( a < -6.0 || a > (a_start + (n_i+1)*delta_a) || c > (c_start + n_k_low*delta_c_low + (n_k_high+1)*delta_c_high) )
  {
    return false;
  }
  else
  {
    return true;
  }
  
}



/*
    Structure of the data file

It is:
    a = ln(md2)
    b = ln(lambda)
    c = beta
    d = cos(theta)
    f = ln(I23)
where md2 and lambda are scaled by sqrt(s)!

Empty lines and lines with leading "#" are ignored. Each line should contain a single numerical value of type double
or the string "-inf" which is internally converted. Additional entries (separated by tabs) in a line would be ignored.

The entries must be the function values at fixed points (a[i],b[j],c[k],d[l]) such that the following structure is satisfied:

a[0]     b[0]     c[0]     d[0]     f(a[0],b[0],c[0],d[0])
..
a[0]     b[0]     c[0]     d[n_l]   f(a[0],b[0],c[0],d[n_l])
a[0]     b[0]     c[1]     d[0]     f(a[0],b[0],c[1],d[0])
..
a[0]     b[0]     c[1]     d[n_l]   f(a[0],b[0],c[1],d[n_l])
....
a[0]     b[0]     c[n_k]   d[0]     f(a[0],b[0],c[n_k],d[0])
..
a[0]     b[0]     c[n_k]   d[n_l]   f(a[0],b[0],c[n_k],d[n_l])
a[0]     b[1]     c[0]     d[0]     f(a[0],b[1],c[0],d[0])
.....
a[0]     b[n_j]   c[n_k]   d[n_l]   f(a[0],b[n_j],c[n_k],d[n_l])
............
a[n_i]   b[n_j]   c[n_k]   d[n_l]   f(a[n_i],b[n_j],c[n_k],d[n_l])

Here the n_i, n_j .. are the numbers of tabulated points in the corresponding directions (-1 due to offset 0 vectors), see
interpolation23::get_index(..). The data file must only contain the function values f(...)!

*/





/**
 * Polynomial interpolation and extrapolation as found in Numerical Recipes. Arrays are unit-offset!
 *
 * @param[in] xa[] Array that holds the x-values x[1] .. x[n]
 * @param[in] ya[] Array that holds the y-values y[1] .. y[n]
 * @param[in] n Number of points to inter- / extrapolate from. n-1 is the degree of the polynomial used.
 * @param[in] x Point at which the interpolation polynomial should be evaluated.
 * @param[out] y The inter- / extrapolated value y(x)
 * @param[out] dy An error estimate
 */
void interpolation23::polynomialInterpolation( const vector<double>& xa, const vector<double>& ya,
                                               const int n, const double x, double& y, double& dy ) const
{

  int ns = 1;
  double den, dif, dift, ho, hp, w;
  double *c, *d;

  dif = fabs( x - xa[1] );

  c = new double[n+1];
  d = new double[n+1];

  for ( int i = 1; i <= n; i++ )
  {
    if (( dift = fabs( x - xa[i] ) ) < dif )
    {
      ns = i;
      dif = dift;
    }
    c[i] = ya[i];
    d[i] = ya[i];
  }

  y = ya[ns--];

  for ( int m = 1; m < n; m++ )
  {
    for ( int i = 1; i <= n - m; i++ )
    {
      ho = xa[i] - x;
      hp = xa[i+m] - x;
      w = c[i+1] - d[i];
      if (( den = ho - hp ) == 0.0 )  //this error can occur only if two input xa's are (to within roundoff) identical
      {
        cout << "Error in routine polint" << endl;
      }
      den = w / den;
      d[i] = hp * den;           // c's and d's are updated
      c[i] = ho * den;
    }

    y += ( dy = ( 2 * ns < ( n - m ) ? c[ns+1] : d[ns--] ) );
  }

  delete[] c;
  delete[] d;
}


