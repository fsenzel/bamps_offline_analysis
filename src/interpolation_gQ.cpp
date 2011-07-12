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

#include "interpolation_gQ.h"
#include "FPT_compare.h"
#include "configuration.h"

using namespace std;


interpolation_gQ::interpolation_gQ()
{

}

void interpolation_gQ::configure(const partclType species_arg, const bool asRun_arg)
{
  asRun = asRun_arg;
  species = species_arg;
  
  // tables are generated for M = 1.5 GeV for charm and M = 4.8 GeV for bottom
  if(species == charm)
  {
    if(Mcharm == 1.5)
    {
      Mass = Mcharm;
      // total numbers of tabulated points on the a, b "axes"
      n_i = 113;
      n_j = 225;
      // spacings of tabulated values in a, b, c and d-direction
      // in c-direction there are different spacings for small than for large c-values
      delta_a = 0.1;
      delta_b = 0.02;
    }
    else
      cout << "IgcData table generated for M_c = 1.5 GeV and not " << Mcharm << endl;
  }
  else if(species == bottom)
  {
    if(Mbottom == 4.8)
    {
      Mass = Mbottom;
      // total numbers of tabulated points on the a, b "axes"
      n_i = 113;
      n_j = 217;
      // spacings of tabulated values in a, b, c and d-direction
      // in c-direction there are different spacings for small than for large c-values
      delta_a = 0.1;
      delta_b = 0.01;
    }
    else
      cout << "IgbData table generated for M_b = 4.8 GeV and not " << Mbottom << endl;
  }
  else
    cout << "error in interpolation_gQ::configure()" << endl;

  // the smallest tabulated values of a, b, c and d
  a_start = log( 0.0001 );
  b_start = log( pow(Mass,2.0)+0.001);
  
  
  expected_entries = n_i * n_j;
  //allocate memory based on the known size (# of entries) of the table
  //a matter of efficieny - not necessary
  IgQdata.reserve( expected_entries + 1 );

  if ( readTables() == IgQ_SUCCESS )
  {
    cout << "IgQ readout successful" << endl;
  }
}


interpolation_gQ::~interpolation_gQ()
{
  IgQdata.clear();
}



IgQ_ERROR_CODE interpolation_gQ::readTables()
{
  string filename;   //data file holding the tabulated values, see below for required structure
  if(species == charm)
  {
    if(asRun)
      filename = "data/Igc_15_asRun_table.dat";
    else
//       filename = "data/Igc_15_asConsttzero_table.dat";
      filename = "data/Igc_15_asConst_table.dat";
  }
  else if(species == bottom)
  {
    if(asRun)
      filename = "data/Igb_48_asRun_table.dat";
    else
//       filename = "data/Igb_48_asConsttzero_table.dat"; 
      filename = "data/Igb_48_asConst_table.dat";
  }
    
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
            IgQdata.push_back( content );
          }
          else
          {
            cout << "severe error in interpolation_gQ::readTables() - ill-formated data file" << endl;
            return IgQ_READOUT_ERROR;
          }
        }
        else     //i.e. when string "-inf" is found
        {
          cout << "error  in interpolation_gQ::readTables() - inf entries" << endl;
          return IgQ_READOUT_ERROR;
        }

      }
    }
    while ( !file.eof() );
  }
  else
  {
    // throw an exception if the file could not have been read
    throw eIgQ_read_error( "error in interpolation_gQ::readTables() - could not open " + filename );
    return IgQ_FILE_ERROR;
  }

  // The number of tabulated values read from file must be equal to the expected number, otherwise
  // the interpolation routine will fail to work correctly!
  if ( IgQdata.size() != expected_entries )
  {
    cout << IgQdata.size() << " instead of " << expected_entries << "entries" << endl;
    // throw an exception if the number of entries in the file is not right
    throw eIgQ_read_error( "error in interpolation_gQ::interpolation_gQ() - unexpected number of entries in IgQ table" );
    return IgQ_SIZE_ERROR;
  }
  else
    return IgQ_SUCCESS;
}



/**
 * 4-dimensional interpolation / extrapolation of the tabulated ln(IgQ)-values at a given point (a,b,c,d),
 * with a=ln(md2), b=ln(lambda), c=beta and d=theta.
 * This routine serves as a wrapper to different methods of interpolation (linear, polynomial, cubic splines). 
 *
 * @param[in] a ln(md2_scaled)
 * @param[in] b ln(lambda_scaled)
 * @param[in] c beta
 * @param[in] d cos(theta)
 * @return IgQ = exp(f), the value of IgQ at (a,b,c,d) found via interpolation from the tabulated values of f
 */
double interpolation_gQ::getIgQ( const double a, const double b) const
{
  
  return getIgQ_linearInterpolation( a, b);
  
}



/**
 * 2-dimensional linear interpolation of the tabulated ln(IgQ)-values at a given point (a,b),
 * with a=ln(md2), b=ln(s).
 * This routine is closely tailored to the known structure of the tabulated data (stored in interpolation_gQ::IgQdata).
 * See below for a short documentation of the required data structure.
 *
 * @param[in] a ln(md2)
 * @param[in] b ln(s)
 * @return IgQ = exp(f), the value of IgQ at (a,b) found via interpolation from the tabulated values of f
 */
double interpolation_gQ::getIgQ_linearInterpolation( const double a, const double b) const
{
  // find indices i,j such that
  // a[i] <= a < a[i+1], b[j] <= b < b[j+1]
  int i = int(( a - a_start ) / delta_a ); 
  int j = int(( b - b_start ) / delta_b );

  if ( i < 0 || j < 0)
  {
    cout << "unexpected - (i,j) not in range" << endl;
    cout << n_i << "  " << n_j << "  "<< endl;
    cout << i << "  " << j   << endl;
    cout << a << "  " << b  << endl;
    cout << a_start << "  " << b_start  << endl;

//     return 0;
    if(i < 0)
      i = 0;
    if(j < 0)
      j = 0;
  }

  // temporary vector objects used to initialize the multi-dimensional arrays based on STL vectors
  vector<double> t1( 2, 0 );
  vector< vector<double> > f( 2, t1 );

  int delta_ex_a,delta_ex_b;
  int index_i, index_j;
  
  if ( i < (n_i - 1) ) // extrapolated linearly in a direction
  {
    index_i = i+1;
    delta_ex_a = 1;
    
  }
  else
  {
    index_i = (n_i - 1);
    delta_ex_a = 2;
  }
  
  if ( j < (n_j - 1) ) // extrapolated linearly in b direction
  {
    index_j = j+1;
    delta_ex_b = 1;
  }
  else
  {
    index_j = (n_j - 1);
    delta_ex_b = 2;
  }

  for ( int p = 0; p <= 1; p++ )
    for ( int q = 0; q <= 1; q++ )
      f[p][q] = IgQdata[ get_index( index_i + (( p-1 )*delta_ex_a ), index_j + (( q-1 )*delta_ex_b ) )]; // without extrapolation this becomes simply IgQdata[ get_index( i+p, j+q )]

  const double a_i = a_start + ( index_i - delta_ex_a ) * delta_a; // without extrapolation this is simply a_i = a_start + i * delta_a
  const double delta_a_local = delta_a * delta_ex_a; // without extrapolation this is simply delta_a
  
  const double b_j = b_start + ( index_j - delta_ex_b ) * delta_b; // without extrapolation this is simply b_i = b_start + i * delta_b
  const double delta_b_local = delta_b * delta_ex_b; // without extrbpolbtion this is simply deltb_b

  // f_a[q] holds the values of f[p][q] linearly interpolated in a-direction.
  // Similar to the notation above, f_a[0] corresponds to the value at (a,b[j])
  // and f_a[1] corresponds to the value at (a,b[j+1]).
  // Note that these points are at a not at a[i+p]!
  vector<double> f_a( 2, 0 );
  for ( int q = 0; q <= 1; q++ )
    f_a[q] = f[0][q] + ( a - a_i ) / delta_a_local * ( f[1][q] - f[0][q] );

  // f_ab holds the values of f_a[q] linearly interpolated in b-direction.
  // thus the final result!
  double f_ab = f_a[0] + ( b - b_j ) / delta_b_local * ( f_a[1] - f_a[0] );

  return exp( f_ab );
}







int interpolation_gQ::get_index( const int i, const int j) const
{
  /*
  Utility routine to get the index number (line number - 1) of the point P=(a[i],b[j],c[k],d[l]) such that
  IgQdata[index] is the tabulated value corresponding to the point P.
  */

  return (( i * n_j ) + ( j ) );
  //add 1 to get the line number of the desired value in the data file (IgQdata vector has offset 0, line numbers have offset 1)
}




