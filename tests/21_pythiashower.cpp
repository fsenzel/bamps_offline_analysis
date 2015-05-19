//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/trunk/tests/21_checkinit.cpp $
//$LastChangedDate: 2014-03-21 13:29:33 +0100 (Fri, 21 Mar 2014) $
//$LastChangedRevision: 1641 $
//$LastChangedBy: gallmei $
//---------------------------------------------
//---------------------------------------------

#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <exception>
#include <stdint.h>
#include <execinfo.h>
#include <signal.h>
#include <vector>
#include <boost/regex.hpp>

#include "particleOffline.h"
#include "FPT_compare.h"
#include "globalsettings.h"

using namespace std;
//using std::cout;
//using std::endl;
//using std::vector;

const double Ntest = 15;
const double NaddedEvents = 5000;


class analysisRapidityRange
{
  public:
    analysisRapidityRange() : yleft(0), yright(0) {};
    analysisRapidityRange( const double _left, const double _right ) : yleft(_left), yright(_right) {};
  
    void reset( const double _left, const double _right ) { yleft = _left; yright = _right; }
    double ydelta() const { return (yright - yleft); }
    
    double yleft;
    double yright;
};


class v2RAA
{
public:
  v2RAA( string name_arg, string filename_prefix_arg, std::vector<analysisRapidityRange> rapidityRanges_arg, const double pt_min_arg = 5.0, const double pt_max_arg = 40.0, const int n_g_arg = 35 );
  
  void setPtBinProperties( const double pt_min_arg, const double pt_max_arg, const int n_g_arg )
  {
    pt_min = pt_min_arg;
    pt_max = pt_max_arg;
    n_g = n_g_arg;
  }
  
  void computeFor( const FLAVOR_TYPE _flavTypeToComputeFor, vector< ParticleOffline >& _particles, const int n_particles, string additionalNameTag );

private:
  
  double pt_min,pt_max;  
  int n_c,n_b,n_g,eta_bins; 
  
  string name,filename_prefix;
  
  std::vector<analysisRapidityRange> rapidityRanges;
};


v2RAA::v2RAA( string name_arg, string filename_prefix_arg, vector< analysisRapidityRange > rapidityRanges_arg, const double pt_min_arg, const double pt_max_arg, const int n_g_arg ):
    name( name_arg ), filename_prefix( filename_prefix_arg ), rapidityRanges( rapidityRanges_arg ), pt_min( pt_min_arg ), pt_max( pt_max_arg ), n_g( n_g_arg )
{
  eta_bins = rapidityRanges.size();
}


void v2RAA::computeFor( const FLAVOR_TYPE _flavTypeToComputeFor, vector<ParticleOffline>& _particles, const int n_particles, string additionalNameTag )
{
  double eta, pp, pt, v2, xt;
  double sinAlpha, alpha;
  int dummy, flavor, n_bins;
  int alphaIndex;
  
  double _pt_min, _pt_max;

  string filename_v2, filename_v2_summed, filename_v2_tot, filename_yield, filename_pt_angleDependence, type;
  
  type = Particle::getName( _flavTypeToComputeFor );
  
  n_bins = n_g;
  _pt_max = pt_max;
  _pt_min = pt_min;
  
  // avoid problem with binning pt logaritmitically: cannot deal with pt = 0
  if( _pt_min < 0.1 )
    _pt_min = 0.1;
  
  const double d_ln_pt = ( log( _pt_max ) - log( _pt_min ) ) / n_bins;

  double v2sum[eta_bins];
  int NmbInRange[eta_bins];
  int NmbInnerRegion = 0;
  for ( int j = 0;j < eta_bins;j++ )
  {
    v2sum[j] = 0.0;
    NmbInRange[j] = 0;
  }

  double ptBinsV2[eta_bins][n_bins+1];
  int ptBinsNmb[eta_bins][n_bins+1];
  int ptBinsInnerRegion[n_bins+1];
  for ( int j = 0;j < n_bins + 1;j++ )
  {
    ptBinsInnerRegion[j] = 0;
    for ( int i = 0;i < eta_bins;i++ )
    {
      ptBinsV2[i][j] = 0.0;
      ptBinsNmb[i][j] = 0.0;
    }
  }
  
  const double deltaAlpha = 15; // degrees
  const int nAlphaBins = 6;  // 90° / 15°
  double ptBinsAngleDep[eta_bins][nAlphaBins][n_bins+1];
  for ( int i = 0; i < eta_bins; i++ )
  {
    for ( int j = 0; j < nAlphaBins; j++ )
    {
      for ( int k = 0; k < n_bins + 1; k++ )
      {
        ptBinsAngleDep[i][j][k] = 0;
      }
    }
  }
  

  // compute v2 and bin it into pt bins
  for ( int i = 0; i < n_particles; i++ )
  {
    pt = _particles[i].Mom.Pt();
    xt = _particles[i].Pos.Pt();

    sinAlpha = _particles[i].Mom.Py() / pt;
    alpha = asin( fabs( sinAlpha ) );
    alpha = alpha * 180 / M_PI;
    
    alphaIndex = static_cast<int>( alpha / deltaAlpha );
    if ( alphaIndex >= nAlphaBins )
    {
      alphaIndex = nAlphaBins - 1;
    }
        
    // for most particle species we are interested in the pseudorapidity (for massless particle it does not matter anyhow)
    eta = _particles[i].Mom.Pseudorapidity(_particles[i].m);

    v2 = ( pow( _particles[i].Mom.Px(), 2.0 ) - pow( _particles[i].Mom.Py(), 2.0 ) ) / pow( pt, 2.0 );

    flavor = _particles[i].FLAVOR;
    FLAVOR_TYPE genFlavor = ParticleOffline::mapToGenericFlavorType( static_cast<FLAVOR_TYPE>( flavor ) );
    
    if( ( _flavTypeToComputeFor == flavor ) || 
          ( _flavTypeToComputeFor == light_quark && ( genFlavor == light_quark || genFlavor == anti_light_quark ) ) ||
          ( _flavTypeToComputeFor == genFlavor )
      ) 
    {
      // individually check for each rapidity range whether this particle needs to be binned
      for ( int yRangeIndex = 0; yRangeIndex < eta_bins; yRangeIndex++ )
      {
        if ( fabs( eta ) >= rapidityRanges[yRangeIndex].yleft && fabs( eta ) <= rapidityRanges[yRangeIndex].yright )
        {
          v2sum[yRangeIndex] += v2;
          NmbInRange[yRangeIndex]++;

          if ( pt <= _pt_max && pt > _pt_min )
          {
            dummy = int(( log( pt ) - log( _pt_min ) ) / d_ln_pt );
            ptBinsV2[yRangeIndex][dummy] += v2;
            ptBinsNmb[yRangeIndex][dummy]++;
            ptBinsAngleDep[yRangeIndex][alphaIndex][dummy]++;
          }  
        }
      }
    }
  }

  int binMax = 0;
  int binMin = n_particles;
  for ( int k = 0; k < n_bins + 1; k++ )
  {
    if ( ptBinsNmb[0][k] > binMax )
    {
      binMax = ptBinsNmb[0][k];
    }
    if ( ptBinsNmb[0][k] < binMin && ptBinsNmb[0][k] != 0 )
    {
      binMin = ptBinsNmb[0][k];
    }
  }

  // file output
  double pt_out;

  filename_v2_summed = filename_prefix + "_" + type + "_" + additionalNameTag + "_v2_pt_summed_" + name;
  filename_v2_tot = filename_prefix + "_" + type + "_" + additionalNameTag + "_v2_tot_" + name;
  filename_yield = filename_prefix + "_" + type + "_" + additionalNameTag + "_yield_pt_" + name;
  filename_pt_angleDependence = filename_prefix + "_" + type + "_" + additionalNameTag + "_pt_angular_dependence_" + name;

  fstream print_v2_summed( filename_v2_summed.c_str(), ios::out | ios::trunc );
  fstream print_v2_tot( filename_v2_tot.c_str(), ios::out | ios::trunc );
  fstream print_yield( filename_yield.c_str(), ios::out | ios::trunc );
  fstream print_pt_angleDependence( filename_pt_angleDependence.c_str(), ios::out | ios::trunc );

//   cout << "total v2 of " << type << " in Y=+-0.35 = " << v2sum[0]/NmbInRange[0] << endl;

  // print total v2
  print_v2_tot << "# total v2 of " << type << endl;
  print_v2_tot << "# bin statistics for 0.35 mid-rapidity:  Avg per bin=" << double( NmbInRange[0] ) / n_c << "   Min=" << binMin << "   Max=" << binMax << endl;
  print_v2_tot << "# total v2, v2_sum and number in range for different rapidity bins" << endl;

  print_v2_tot << _pt_min;
  print_v2_tot.width( 15 );
  for ( int i = 0;i < eta_bins;i++ )
  {
    if ( NmbInRange[i] > 0 )
    {
      print_v2_tot <<  v2sum[i] / NmbInRange[i];
    }
    else
    {
      print_v2_tot <<  0;
    }
    print_v2_tot.width( 15 );
  }
  for ( int i = 0;i < eta_bins;i++ )
  {
    print_v2_tot <<  v2sum[i];
    print_v2_tot.width( 15 );
    print_v2_tot <<  NmbInRange[i];
    print_v2_tot.width( 15 );
  }
  print_v2_tot << endl;

  print_v2_tot << _pt_max;
  print_v2_tot.width( 15 );
  for ( int i = 0;i < eta_bins;i++ )
  {
    if ( NmbInRange[i] > 0 )
    {
      print_v2_tot <<  v2sum[i] / NmbInRange[i];
    }
    else
    {
      print_v2_tot <<  0;
    }
    print_v2_tot.width( 15 );
  }
  for ( int i = 0;i < eta_bins;i++ )
  {
    print_v2_tot <<  v2sum[i];
    print_v2_tot.width( 15 );
    print_v2_tot <<  NmbInRange[i];
    print_v2_tot.width( 15 );
  }
  print_v2_tot << endl;


  // print summed output, v2 is not computed, but summed v2 and the number in one bin
  print_v2_summed << "# summed v2 of " << type << endl;
  print_v2_summed << "#";
  print_v2_summed.width( 14 );
  print_v2_summed << "pt       summed v_2 and number in bin for different rapidity bins" << endl;

  for ( int k = 0;k < n_bins + 1;k++ )
  {
    print_v2_summed.width( 15 );
    pt_out = exp( double( k ) * d_ln_pt + log( _pt_min ) + d_ln_pt / 2.0 ); // +dpt/2.0 to shift value in the middle of the intervall, important for small number of bins
    print_v2_summed << pt_out;
    for ( int i = 0;i < eta_bins;i++ )
    {
      print_v2_summed.width( 15 );
      print_v2_summed << ptBinsV2[i][k];
      print_v2_summed.width( 10 );
      print_v2_summed << ptBinsNmb[i][k];
    }
    print_v2_summed << endl;
  }

  // print yield for RAA
  print_yield << "# " << type << " yield distribution" << endl;
  print_yield << "#";
  print_yield.width( 14 );
  print_yield << "pt       yield for different rapidity bins" << endl;

  for ( int k = 0;k < n_bins + 1;k++ )
  {
    print_yield.width( 15 );
    pt_out = exp( double( k ) * d_ln_pt + log( _pt_min ) + d_ln_pt / 2.0 ); // +dpt/2.0 to shift value in the middle of the intervall, important for small number of bins
    print_yield << pt_out;
    const double dpt = exp( double( k ) * d_ln_pt + log( _pt_min ) + d_ln_pt ) - exp( double( k ) * d_ln_pt + log( _pt_min ) );
    for ( int i = 0;i < eta_bins;i++ )
    {
      const double delta_eta = 2.0 * ( rapidityRanges[i].yright - rapidityRanges[i].yleft );
      
      double nInBin = double( ptBinsNmb[i][k] ) / Ntest / dpt / delta_eta / NaddedEvents;
      
      print_yield.width( 15 );
      print_yield << nInBin;
    }
    print_yield << endl;
  }
  
  
  // print yield for RAA for different angles with respect to the reaction plane
  print_pt_angleDependence << "# " << type << " yield distribution for different angles (alpha) with respect to the reaction plane" << endl;
  print_pt_angleDependence << "#";
  print_pt_angleDependence.width( 14 );
  print_pt_angleDependence << "pt       yield for different rapidity bins" << endl;
  for ( int j = 0; j < nAlphaBins; j++ )
  {
    print_pt_angleDependence << "#alpha in [ " << j * deltaAlpha << "°, " << (j+1)*deltaAlpha << "° ] "<< endl;
    for ( int k = 0;k < n_bins + 1;k++ )
    {
      print_pt_angleDependence.width( 15 );
      pt_out = exp( double( k ) * d_ln_pt + log( _pt_min ) + d_ln_pt / 2.0 ); // +dpt/2.0 to shift value in the middle of the intervall, important for small number of bins
      print_pt_angleDependence << pt_out;
      const double dpt = exp( double( k ) * d_ln_pt + log( _pt_min ) + d_ln_pt ) - exp( double( k ) * d_ln_pt + log( _pt_min ) );
      for ( int i = 0;i < eta_bins;i++ )
      {
        const double delta_eta = 2.0 * ( rapidityRanges[i].yright - rapidityRanges[i].yleft );
        
        print_pt_angleDependence.width( 15 );
        print_pt_angleDependence << double( ptBinsAngleDep[i][j][k] ) / Ntest / dpt / delta_eta;
      }
      print_pt_angleDependence << endl;
    }    
    print_pt_angleDependence << endl;
    print_pt_angleDependence << endl;
  }
  
}

void readFile( string _filename, vector<ParticleOffline>& _particles )
{
  std::fstream file ( _filename.c_str(), ios::in );

  if ( file )
  {
    string line;
    stringstream inStream;
    double dummy;
    vector<double> lineContent;

    do
    {
      getline ( file, line );
      if ( line.find ( "#", 0 ) == string::npos && line.length() != 0 )   //ignore empty lines and lines with leading "#"
      {
        inStream.str ( "" );
        inStream.clear();
        lineContent.clear();

        inStream.str ( line );
        while ( !inStream.eof() )
        {
          inStream >> dummy;
          lineContent.push_back ( dummy );
        }

        if ( lineContent.size() == 10 )
        {
          ParticleOffline tempParticle;
          tempParticle.Mom = VectorEPxPyPz( lineContent[4], lineContent[1], lineContent[2], lineContent[3] );
          tempParticle.Pos = VectorTXYZ( lineContent[8], lineContent[5], lineContent[6], lineContent[7] );
          tempParticle.FLAVOR = static_cast<FLAVOR_TYPE>( lineContent[9] );

          _particles.push_back( tempParticle );
        }
        else
        {
          string errMsg = "severe error in readFile() - line too short or too long";
          cout << errMsg << endl;
        }

        if ( inStream.fail() )
        {
          string errMsg = "severe error in readFile() - ill-formated data file";
          cout << errMsg << endl;
        }
      }
    }
    while ( !file.eof() );
  }
  else
  {
    string errMsg = "severe error in readFile() - file not found";
    cout << errMsg << endl; 
  }
}

/**
 * @brief The main routine of BAMPS
 */
int main(int argc, char *argv[])
{
  if ( !( argc == 3 ) )
  {
    cout << "Usage: ./21_pythiashower <input filename> <minimum pt of spectrum>" << endl;
  }
  {
    stringstream filename( argv[1] );
    stringstream minPt_string( argv[2] );
    double minPt;
    minPt_string >> minPt;
//     std::cout << "min pt of binning output = " << minPt << std::endl;
  
    string filename_woSuffix, basename, outputFolder;
    boost::regex suffix_regex( "_unshoweredParticles" );
    boost::regex prefix_regex( "input/" );
    boost::smatch matches_prefix;
    boost::smatch matches_suffix;
    if ( boost::regex_search( filename.str(), matches_suffix, suffix_regex ) )
    {
      filename_woSuffix = boost::regex_replace( filename.str(), suffix_regex, "" );
    }
    if ( boost::regex_search( filename_woSuffix, matches_prefix, prefix_regex ) )
    {
      basename = boost::regex_replace( filename_woSuffix, prefix_regex, "" );
    }
    
    outputFolder = "output/ptMin" + minPt_string.str() + "GeV/";
    
    analysisRapidityRange yRange;
    std::vector<analysisRapidityRange> rapidityRanges;
    yRange.reset( 0, 0.8 );
    rapidityRanges.push_back(yRange);
    yRange.reset( 0, 0.5 );
    rapidityRanges.push_back(yRange);
    yRange.reset( 0, 1.0 );
    rapidityRanges.push_back(yRange);
    yRange.reset( 0, 2.0 );
    rapidityRanges.push_back(yRange);

    std::vector<ParticleOffline> unshowered_partons;
    readFile( filename.str(), unshowered_partons );

    v2RAA theV2RAA( "initial", outputFolder + basename, rapidityRanges );
  
    const double pt_min_v2RAA = minPt;
    double pt_max_v2RAA, nbins_v2RAA;
    if( pt_min_v2RAA < 35 )
    {
      pt_max_v2RAA = 55.0;
      nbins_v2RAA = 67; // good for pt_min = 0
      if( FPT_COMP_GE( pt_min_v2RAA, 6.0 ) )
        nbins_v2RAA = 25;
      else if( FPT_COMP_GE( pt_min_v2RAA, 3.0 ) )
        nbins_v2RAA = 34;
    }
    else
    {
      pt_max_v2RAA = 150.0;
      nbins_v2RAA = 60;
    
      if( FPT_COMP_GE( pt_min_v2RAA, 70.0 ) )
        nbins_v2RAA = 52;
      else if( FPT_COMP_GE( pt_min_v2RAA, 90.0 ) )
        nbins_v2RAA = 45;
    }
  
    theV2RAA.setPtBinProperties( pt_min_v2RAA, pt_max_v2RAA, nbins_v2RAA );
    theV2RAA.computeFor( gluon, unshowered_partons, unshowered_partons.size(), "unshowered_partons" );
    theV2RAA.computeFor( light_quark, unshowered_partons, unshowered_partons.size(), "unshowered_partons" );
//     theV2RAA.computeFor( up, unshowered_partons, unshowered_partons.size(), "unshowered_partons" );
//     theV2RAA.computeFor( down, unshowered_partons, unshowered_partons.size(), "unshowered_partons" );
//     theV2RAA.computeFor( strange, unshowered_partons, unshowered_partons.size(), "unshowered_partons" );
//     theV2RAA.computeFor( anti_up, unshowered_partons, unshowered_partons.size(), "unshowered_partons" );
//     theV2RAA.computeFor( anti_down, unshowered_partons, unshowered_partons.size(), "unshowered_partons" );
//     theV2RAA.computeFor( anti_strange, unshowered_partons, unshowered_partons.size(), "unshowered_partons" );

    return EXIT_SUCCESS;
  }
}
