//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------
//---------------------------------------------


#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <string>
#include <sstream>
#include <time.h>

#include "particle.h"
#include "analysis.h"
#include "binning.h"
#include "binning2.h"
#include "configuration.h"

#include <stdio.h> // for getenv()
#include <stdlib.h> // for getenv()
#include <time.h>


using namespace std;
using namespace ns_casc;


analysis::analysis( config* const c ):
    theConfig( c ),
    rings( c->getRingNumber(), c->getCentralRingRadius(), c->getDeltaR() ),
    centralRingsCopyFromCascade( c->getRingNumber(), c->getCentralRingRadius(), c->getDeltaR() )
{
  //---get time and date---
  time( &start );
  //-----------------------

//   tstep[0]=6.1;        //fm/c
//   tstep[1]=infinity;   //fm/c
//   nTimeSteps = 2;      //6 timesteps during the cascade, plus 1 before

  //---- times for output of data --------
  tstep[0] = 0.1;
  tstep[1] = 0.5;      //fm/c
  tstep[2] = 1.0;      //fm/c
  tstep[3] = 1.5;      //fm/c
  tstep[4] = 2.0;      //fm/c
  tstep[5] = 2.5;      //fm/c
  tstep[6] = 3.0;      //fm/c
  tstep[7] = 3.5;      //fm/c
  tstep[8] = 4.0;      //fm/c
  tstep[9] = 4.5;      //fm/c
  tstep[10] = 5.0;      //fm/c
  tstep[11] = 5.5;      //fm/c
  tstep[12] = 6.0;      //fm/c
  tstep[13] = 6.5;      //fm/c
  tstep[14] = 7.0;      //fm/c
  tstep[15] = 7.5;      //fm/c
  tstep[16] = 8.0;      //fm/c
  tstep[17] = 9.0;      //fm/c
  tstep[18] = 10.0;
  tstep[19] = infinity; //fm/c
  nTimeSteps = 20;
  //--------------------------------------


  //---- times for output of data --------
  int tempCount_tstep = 0;
  tstep_movie[tempCount_tstep] = 0.1;
  do {
    ++tempCount_tstep;
    tstep_movie[tempCount_tstep] = tstep_movie[tempCount_tstep - 1] + 0.1;
  } while (tstep_movie[tempCount_tstep] <= 10.0 );
  
  ++tempCount_tstep;
  tstep_movie[tempCount_tstep] = infinity;
  
  nTimeSteps_movie = tempCount_tstep + 1;
  //--------------------------------------

  //------- filename prefix ---------
  char * temp = getenv( "PBS_JOBID" );
  if ( temp != NULL && !theConfig->isLocalCluster() )
  {
    string jobID( temp );

    filename_prefix = "/local/" + jobID + "/" + theConfig->getJobName();
  }
  else
  {
    filename_prefix = "./output/" + theConfig->getJobName();
  }
  //--------------------------------------

  // oscar movie output
  string oscarName = "";
  oscarName = filename_prefix + "_background.oscar";
  oscarBackground.open( oscarName.c_str(), ios::out );
  
  oscarName = filename_prefix + "_jets.oscar";
  oscarJets.open( oscarName.c_str(), ios::out );
 //--------------------------------------
 
  // mfp output for jets
  string mfpName = filename_prefix + "_mfp_jets";
  mfpJetsOutputFile.open( mfpName.c_str(), ios::out );
  //--------------------------------------
 

  // output of central energy densities
  string centralDensitiesName = filename_prefix + "_central_density";
  centralDensitiesOutputFile.open( centralDensitiesName.c_str(), ios::out );
  //--------------------------------------


  jetTracking_PT = 10.0;

  //---- initialisation of PT-binning ----
  minPT = 1.4;
  maxPTSoft = 3.0;
  maxPT = 34.4;
  binWidthPT = 1.0;
  numberBinsPT = int(( maxPT - minPT + 0.001 ) / binWidthPT );
  ptBinsDY1_gluons = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
  ptBinsDY2_gluons = new vector<double>[nTimeSteps+2];
  ptBinsDY3_gluons = new vector<double>[nTimeSteps+2];
  ptBinsDY4_gluons = new vector<double>[nTimeSteps+2];
  ptBinsDY5_gluons = new vector<double>[nTimeSteps+2];
  ptBinsAll_gluons = new vector<double>[nTimeSteps+2];

  ptBinsDY1_quarks = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
  ptBinsDY2_quarks = new vector<double>[nTimeSteps+2];
  ptBinsDY3_quarks = new vector<double>[nTimeSteps+2];
  ptBinsDY4_quarks = new vector<double>[nTimeSteps+2];
  ptBinsDY5_quarks = new vector<double>[nTimeSteps+2];
  ptBinsAll_quarks = new vector<double>[nTimeSteps+2];
  
  ptBinsDY1_ups = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
  ptBinsDY2_ups = new vector<double>[nTimeSteps+2];
  ptBinsDY3_ups = new vector<double>[nTimeSteps+2];
  ptBinsDY4_ups = new vector<double>[nTimeSteps+2];
  ptBinsDY5_ups = new vector<double>[nTimeSteps+2];
  ptBinsAll_ups = new vector<double>[nTimeSteps+2];

  ptBinsDY1_downs = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
  ptBinsDY2_downs = new vector<double>[nTimeSteps+2];
  ptBinsDY3_downs = new vector<double>[nTimeSteps+2];
  ptBinsDY4_downs = new vector<double>[nTimeSteps+2];
  ptBinsDY5_downs = new vector<double>[nTimeSteps+2];
  ptBinsAll_downs = new vector<double>[nTimeSteps+2];
  
  ptBinsDY1_stranges = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
  ptBinsDY2_stranges = new vector<double>[nTimeSteps+2];
  ptBinsDY3_stranges = new vector<double>[nTimeSteps+2];
  ptBinsDY4_stranges = new vector<double>[nTimeSteps+2];
  ptBinsDY5_stranges = new vector<double>[nTimeSteps+2];
  ptBinsAll_stranges = new vector<double>[nTimeSteps+2];
  
  ptBinsDY1_anti_ups = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
  ptBinsDY2_anti_ups = new vector<double>[nTimeSteps+2];
  ptBinsDY3_anti_ups = new vector<double>[nTimeSteps+2];
  ptBinsDY4_anti_ups = new vector<double>[nTimeSteps+2];
  ptBinsDY5_anti_ups = new vector<double>[nTimeSteps+2];
  ptBinsAll_anti_ups = new vector<double>[nTimeSteps+2];
  
  ptBinsDY1_anti_downs = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
  ptBinsDY2_anti_downs = new vector<double>[nTimeSteps+2];
  ptBinsDY3_anti_downs = new vector<double>[nTimeSteps+2];
  ptBinsDY4_anti_downs = new vector<double>[nTimeSteps+2];
  ptBinsDY5_anti_downs = new vector<double>[nTimeSteps+2];
  ptBinsAll_anti_downs = new vector<double>[nTimeSteps+2];
  
  ptBinsDY1_anti_stranges = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
  ptBinsDY2_anti_stranges = new vector<double>[nTimeSteps+2];
  ptBinsDY3_anti_stranges = new vector<double>[nTimeSteps+2];
  ptBinsDY4_anti_stranges = new vector<double>[nTimeSteps+2];
  ptBinsDY5_anti_stranges = new vector<double>[nTimeSteps+2];
  ptBinsAll_anti_stranges = new vector<double>[nTimeSteps+2];
  
  ptBinsDY1_all = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
  ptBinsDY2_all = new vector<double>[nTimeSteps+2];
  ptBinsDY3_all = new vector<double>[nTimeSteps+2];
  ptBinsDY4_all = new vector<double>[nTimeSteps+2];
  ptBinsDY5_all = new vector<double>[nTimeSteps+2];
  ptBinsAll_all = new vector<double>[nTimeSteps+2];

  for ( int j = 0; j < nTimeSteps + 2; j++ )          //+2 because of initial and final timesteps
  {
    for ( int i = 0; i <= numberBinsPT; i++ )
    {
      ptBinsDY1_gluons[j].push_back( 0 );
      ptBinsDY2_gluons[j].push_back( 0 );
      ptBinsDY3_gluons[j].push_back( 0 );
      ptBinsDY4_gluons[j].push_back( 0 );
      ptBinsDY5_gluons[j].push_back( 0 );
      ptBinsAll_gluons[j].push_back( 0 );

      ptBinsDY1_quarks[j].push_back( 0 );
      ptBinsDY2_quarks[j].push_back( 0 );
      ptBinsDY3_quarks[j].push_back( 0 );
      ptBinsDY4_quarks[j].push_back( 0 );
      ptBinsDY5_quarks[j].push_back( 0 );
      ptBinsAll_quarks[j].push_back( 0 );
      
      ptBinsDY1_ups[j].push_back( 0 );
      ptBinsDY2_ups[j].push_back( 0 );
      ptBinsDY3_ups[j].push_back( 0 );
      ptBinsDY4_ups[j].push_back( 0 );
      ptBinsDY5_ups[j].push_back( 0 );
      ptBinsAll_ups[j].push_back( 0 );
      
      ptBinsDY1_downs[j].push_back( 0 );
      ptBinsDY2_downs[j].push_back( 0 );
      ptBinsDY3_downs[j].push_back( 0 );
      ptBinsDY4_downs[j].push_back( 0 );
      ptBinsDY5_downs[j].push_back( 0 );
      ptBinsAll_downs[j].push_back( 0 );
      
      ptBinsDY1_stranges[j].push_back( 0 );
      ptBinsDY2_stranges[j].push_back( 0 );
      ptBinsDY3_stranges[j].push_back( 0 );
      ptBinsDY4_stranges[j].push_back( 0 );
      ptBinsDY5_stranges[j].push_back( 0 );
      ptBinsAll_stranges[j].push_back( 0 );
      
      ptBinsDY1_anti_ups[j].push_back( 0 );
      ptBinsDY2_anti_ups[j].push_back( 0 );
      ptBinsDY3_anti_ups[j].push_back( 0 );
      ptBinsDY4_anti_ups[j].push_back( 0 );
      ptBinsDY5_anti_ups[j].push_back( 0 );
      ptBinsAll_anti_ups[j].push_back( 0 );
      
      ptBinsDY1_anti_downs[j].push_back( 0 );
      ptBinsDY2_anti_downs[j].push_back( 0 );
      ptBinsDY3_anti_downs[j].push_back( 0 );
      ptBinsDY4_anti_downs[j].push_back( 0 );
      ptBinsDY5_anti_downs[j].push_back( 0 );
      ptBinsAll_anti_downs[j].push_back( 0 );
      
      ptBinsDY1_anti_stranges[j].push_back( 0 );
      ptBinsDY2_anti_stranges[j].push_back( 0 );
      ptBinsDY3_anti_stranges[j].push_back( 0 );
      ptBinsDY4_anti_stranges[j].push_back( 0 );
      ptBinsDY5_anti_stranges[j].push_back( 0 );
      ptBinsAll_anti_stranges[j].push_back( 0 );
      
      ptBinsDY1_all[j].push_back( 0 );
      ptBinsDY2_all[j].push_back( 0 );
      ptBinsDY3_all[j].push_back( 0 );
      ptBinsDY4_all[j].push_back( 0 );
      ptBinsDY5_all[j].push_back( 0 );
      ptBinsAll_all[j].push_back( 0 );
    }
  }
  //write the bin labels
  for ( int i = 0; i < numberBinsPT; i++ )
  {
    ptBinLabels.push_back( minPT + ( i * binWidthPT ) + ( binWidthPT / 2 ) );
  }
  cout << "number of bins: " << numberBinsPT << "  binWidth: " << binWidthPT << endl;

  
  //------ initialisation of rapidity binning ------
  minY = -6.0;
  maxY = 6.0;
  binWidthY = 0.1;
  numberBinsY = int(( maxY - minY + 0.00001 ) / binWidthY );
    
  yBins_gluon = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
  yBins_up = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
  yBins_down = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
  yBins_strange = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
  yBins_anti_up = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
  yBins_anti_down = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
  yBins_anti_strange = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
  
  for ( int j = 0; j < nTimeSteps + 2; j++ )          //+2 because of initial and final timesteps
  {
    for ( int i = 0; i <= numberBinsY; i++ )
    {
      yBins_gluon[j].push_back( 0 );
      yBins_up[j].push_back( 0 );
      yBins_down[j].push_back( 0 );
      yBins_strange[j].push_back( 0 );
      yBins_anti_up[j].push_back( 0 );
      yBins_anti_down[j].push_back( 0 );
      yBins_anti_strange[j].push_back( 0 );
    }
  }
  
  //write the bin labels
  for ( int i = 0; i < numberBinsY; i++ )
  {
    yBinLabels.push_back( minY + ( i * binWidthY ) + ( binWidthY / 2 ) );
  }
  cout << "number of y-bins: " << numberBinsY << "  binWidth: " << binWidthY << endl;
  //--------------------------------------------------


  //------ initialisation of transverse energy binning ------
  //------ use the same bins as for rapidity binning ------
  transverseEnergyGluons = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
  transverseEnergyQuarks = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
  transverseEnergyAntiQuarks = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
  
  for ( int j = 0; j < nTimeSteps + 2; j++ )          //+2 because of initial and final timesteps
  {
    for ( int i = 0; i <= numberBinsY; i++ )
    {
      transverseEnergyGluons[j].push_back( 0 );
      transverseEnergyQuarks[j].push_back( 0 );
      transverseEnergyAntiQuarks[j].push_back( 0 );
    }
  }
  //--------------------------------------------------

  
  //---- initialisation of softPT-binning ----
  maxPTSoft;
  binWidthSoftPT = 0.1;
  numberBinsSoftPT = int(( maxPTSoft + 0.001 ) / binWidthSoftPT );
  ptBinsSoftAll_gluons = new vector<double>[nTimeSteps+1];
  ptBinsSoftAll_quarks = new vector<double>[nTimeSteps+1];
  ptBinsSoftAll_ups = new vector<double>[nTimeSteps+1];
  ptBinsSoftAll_downs = new vector<double>[nTimeSteps+1];
  ptBinsSoftAll_stranges = new vector<double>[nTimeSteps+1];
  ptBinsSoftAll_anti_ups = new vector<double>[nTimeSteps+1];
  ptBinsSoftAll_anti_downs = new vector<double>[nTimeSteps+1];
  ptBinsSoftAll_anti_stranges = new vector<double>[nTimeSteps+1];
  ptBinsSoftAll_all = new vector<double>[nTimeSteps+1];
  for ( int j = 0; j < nTimeSteps + 1; j++ )          //+1 because of initial timestep
    for ( int i = 0; i < numberBinsSoftPT; i++ )
    {
      ptBinsSoftAll_gluons[j].push_back( 0 );
      ptBinsSoftAll_quarks[j].push_back( 0 );
      ptBinsSoftAll_ups[j].push_back( 0 );
      ptBinsSoftAll_downs[j].push_back( 0 );
      ptBinsSoftAll_stranges[j].push_back( 0 );
      ptBinsSoftAll_anti_ups[j].push_back( 0 );
      ptBinsSoftAll_anti_downs[j].push_back( 0 );
      ptBinsSoftAll_anti_stranges[j].push_back( 0 );
      ptBinsSoftAll_all[j].push_back( 0 );
    }
  //write the bin labels
  for ( int i = 0; i < numberBinsSoftPT; i++ )
  {
    ptSoftBinLabels.push_back(( i * binWidthSoftPT ) + ( binWidthSoftPT / 2 ) );
  }
  //--------------------------------------
}


analysis::~analysis()
{
  //---- get time and calculate real time needed by simulation ----
  time_t end;
  time( &end );
  int secs = ( int )difftime( end, start );
  int hours = secs / 3600, minutes = ( secs % 3600 ) / 60, seconds = secs - 3600 * hours - 60 * minutes;
  cout << hours << "h" << minutes << "m" << seconds << "s" <<  endl;
  //---------------------------------------------------------------

  oscarBackground.close();
  oscarJets.close();
  jetTrackerOutput();

  delete[] ptBinsDY1_gluons;
  delete[] ptBinsDY2_gluons;
  delete[] ptBinsDY3_gluons;
  delete[] ptBinsDY4_gluons;
  delete[] ptBinsDY5_gluons;
  delete[] ptBinsAll_gluons;
  
  delete[] ptBinsDY1_quarks;
  delete[] ptBinsDY2_quarks;
  delete[] ptBinsDY3_quarks;
  delete[] ptBinsDY4_quarks;
  delete[] ptBinsDY5_quarks;
  delete[] ptBinsAll_quarks;
  
  delete[] ptBinsDY1_ups;
  delete[] ptBinsDY2_ups;
  delete[] ptBinsDY3_ups;
  delete[] ptBinsDY4_ups;
  delete[] ptBinsDY5_ups;
  delete[] ptBinsAll_ups;
  
  delete[] ptBinsDY1_downs;
  delete[] ptBinsDY2_downs;
  delete[] ptBinsDY3_downs;
  delete[] ptBinsDY4_downs;
  delete[] ptBinsDY5_downs;
  delete[] ptBinsAll_downs;
  
  delete[] ptBinsDY1_stranges;
  delete[] ptBinsDY2_stranges;
  delete[] ptBinsDY3_stranges;
  delete[] ptBinsDY4_stranges;
  delete[] ptBinsDY5_stranges;
  delete[] ptBinsAll_stranges;
  
  delete[] ptBinsDY1_anti_ups;
  delete[] ptBinsDY2_anti_ups;
  delete[] ptBinsDY3_anti_ups;
  delete[] ptBinsDY4_anti_ups;
  delete[] ptBinsDY5_anti_ups;
  delete[] ptBinsAll_anti_ups;
  
  delete[] ptBinsDY1_anti_downs;
  delete[] ptBinsDY2_anti_downs;
  delete[] ptBinsDY3_anti_downs;
  delete[] ptBinsDY4_anti_downs;
  delete[] ptBinsDY5_anti_downs;
  delete[] ptBinsAll_anti_downs;
  
  delete[] ptBinsDY1_anti_stranges;
  delete[] ptBinsDY2_anti_stranges;
  delete[] ptBinsDY3_anti_stranges;
  delete[] ptBinsDY4_anti_stranges;
  delete[] ptBinsDY5_anti_stranges;
  delete[] ptBinsAll_anti_stranges;

  delete[] ptBinsDY1_all;
  delete[] ptBinsDY2_all;
  delete[] ptBinsDY3_all;
  delete[] ptBinsDY4_all;
  delete[] ptBinsDY5_all;
  delete[] ptBinsAll_all;
  
  delete[] ptBinsSoftAll_gluons;
  delete[] ptBinsSoftAll_quarks;
  delete[] ptBinsSoftAll_ups;
  delete[] ptBinsSoftAll_downs;
  delete[] ptBinsSoftAll_stranges;
  delete[] ptBinsSoftAll_anti_ups;
  delete[] ptBinsSoftAll_anti_downs;
  delete[] ptBinsSoftAll_anti_stranges;
  delete[] ptBinsSoftAll_all;
}



void analysis::collectPtDataInitial()
{
  ptDistribution( gluon, addedParticles, addedParticles.size(), 0 );
  ptDistribution( light_quark, addedParticles, addedParticles.size(), 0 );
  ptDistribution( allFlavors, addedParticles, addedParticles.size(), 0 );
  ptSoftDistribution( gluon, particles_atTimeNow, particles_atTimeNow.size(), 0 );
  ptSoftDistribution( light_quark, particles_atTimeNow, particles_atTimeNow.size(), 0 );
  ptSoftDistribution( allFlavors, particles_atTimeNow, particles_atTimeNow.size(), 0 );
}



void analysis::collectPtData( const int step )
{
  //are issued starting with 0
  //0 is reserved for initial output, thus step is incremented
  ptDistribution( gluon, addedParticles, addedParticles.size(), step + 1 );
  ptDistribution( light_quark, addedParticles, addedParticles.size(), step + 1 );
  ptDistribution( allFlavors, addedParticles, addedParticles.size(), step + 1 );
  ptSoftDistribution( gluon, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
  ptSoftDistribution( light_quark, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
  ptSoftDistribution( allFlavors, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
}


void analysis::collectYDataInitial()
{
  collectYData( -1 );
}


void analysis::collectYData( const int step )
{
  //are issued starting with 0
  //0 is reserved for initial output, thus step is incremented
  yDistribution( gluon, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
  yDistribution( up, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
  yDistribution( down, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
  yDistribution( strange, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
  yDistribution( anti_up, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
  yDistribution( anti_down, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
  yDistribution( anti_strange, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
}


void analysis::collectEtDataInitial()
{
  collectEtData( -1 );
}


void analysis::collectEtData( const int step )
{
  //are issued starting with 0
  //0 is reserved for initial output, thus step is incremented
  transverseEnergyDistribution( gluon, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
  transverseEnergyDistribution( light_quark, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
  transverseEnergyDistribution( anti_light_quark, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
}





void analysis::initialOutput()
{
  computeV2RAA( "initial", 0 );
  particleOutput( 0 );
}



void analysis::intermediateOutput( const int nn )
{
  string name;
  stringstream ss;

  ss << tstep[nn];
  name = ss.str() + "fm";

  computeV2RAA( name, tstep[nn] );
}



void analysis::finalOutput( const double _stoptime )
{
  computeV2RAA( "final", _stoptime );
  onePartclCorrelations();
  twoPartclCorrelations();
  printPtSpectra( gluon );
  printPtSpectra( light_quark );
  printPtSpectra( allFlavors );
  printSoftPtSpectra( gluon );
  printSoftPtSpectra( light_quark );
  printSoftPtSpectra( allFlavors );
  
  particleOutput( nTimeSteps );
  
  printYDistribution();
}



void analysis::movieOutput( const int step, const int jumpSteps )
{
  writePartclMovie( addedParticles, addedParticles.size(), oscarJets, step, jumpSteps );
//   writePartclMovie( particles_atTimeNow, particles_atTimeNow.size(), oscarBackground, step );
}


void analysis::movieOutputMedium( const int step, const int jumpSteps )
{
  writePartclMovie( particles_atTimeNow, particles_atTimeNow.size(), oscarBackground, step, jumpSteps );
}



void analysis::printYDistribution()
{
  time_t end;
  time( &end );

  string filename = filename_prefix + "_" + "_yDistribution";
  fstream file( filename.c_str(), ios::out | ios::trunc );

  //---- print header if file is empty ----
  file.seekp( 0, ios::end );
  long size = file.tellp();
  file.seekp( 0, ios::beg );
  if ( size == 0 )
    printHeader( file, rapidityDistribution, end );
  //---------------------------------------
  
  
  for ( int j = 0; j <= nTimeSteps; j++ )
  {
    if ( tstep[j-1] <= theConfig->getRuntime() || j == 0 || j == nTimeSteps )
    {
      if ( j == 0 )
      {
        file << "#t = 0" << endl;
      }
      else
      {
        file << "#t = " << tstep[j-1] << endl;
      }
      for ( int i = 0; i < numberBinsY; i++ )
      {
        file << yBinLabels[i] << sep << yBins_gluon[j][i] << sep << yBins_up[j][i] << sep << yBins_down[j][i] 
        << sep << yBins_strange[j][i] << sep << yBins_anti_up[j][i] << sep << yBins_anti_down[j][i] << sep 
        << yBins_anti_strange[j][i] << endl; 
      }
      file << endl << endl;
    }
  }
  
  file.close();
  
  
  //-------- quark numbers summed over all rapidities
  
  filename =  filename_prefix + "_" + "_quarkNumbers";
  file.open( filename.c_str(), ios::out | ios::trunc );
  
  //---- print header if file is empty ----
  file.seekp( 0, ios::end );
  size = file.tellp();
  file.seekp( 0, ios::beg );
  if ( size == 0 )
    printHeader( file, quarkNumbers, end );
  //---------------------------------------
  
  double nGluons, nUps, nDowns, nStranges, nAntiUps, nAntiDowns, nAntiStranges;
  
  for ( int j = 0; j <= nTimeSteps; j++ )
  {
    nGluons = nUps = nDowns = nStranges = nAntiUps = nAntiDowns = nAntiStranges = 0;
    if ( tstep[j-1] <= theConfig->getRuntime() || j == 0 || j == nTimeSteps )
    {
      if ( j == 0 )
      {
        file << 0 << sep;
      }
      else
      {
        file << tstep[j-1] << sep;
      }
      
      for ( int i = 0; i < numberBinsY; i++ )
      {
        nGluons += yBins_gluon[j][i];
        nUps += yBins_up[j][i];
        nDowns += yBins_down[j][i];
        nStranges += yBins_strange[j][i];
        nAntiUps += yBins_anti_up[j][i];
        nAntiDowns += yBins_anti_down[j][i];
        nAntiStranges += yBins_anti_strange[j][i];
      }
      file << nGluons << sep << nUps << sep << nDowns << sep << nStranges << sep << nAntiUps << sep << nAntiDowns << sep << nAntiStranges << endl;    
    }
  }
  
  
  file.close();
}




void analysis::printPtSpectra( const FLAVOR_TYPE _flavTypeToComputeFor )
{
  time_t end;
  time( &end );

  vector<double> * _ptBinsAll;
  vector<double> * _ptBinsDY5;
  vector<double> * _ptBinsDY4;
  vector<double> * _ptBinsDY3;
  vector<double> * _ptBinsDY2;
  vector<double> * _ptBinsDY1;

  if ( _flavTypeToComputeFor == gluon )
  {
    _ptBinsAll = ptBinsAll_gluons;
    _ptBinsDY5 = ptBinsDY5_gluons;
    _ptBinsDY4 = ptBinsDY4_gluons;
    _ptBinsDY3 = ptBinsDY3_gluons;
    _ptBinsDY2 = ptBinsDY2_gluons;
    _ptBinsDY1 = ptBinsDY1_gluons;
  }
  else if ( _flavTypeToComputeFor == light_quark )
  {
    _ptBinsAll = ptBinsAll_quarks;
    _ptBinsDY5 = ptBinsDY5_quarks;
    _ptBinsDY4 = ptBinsDY4_quarks;
    _ptBinsDY3 = ptBinsDY3_quarks;
    _ptBinsDY2 = ptBinsDY2_quarks;
    _ptBinsDY1 = ptBinsDY1_quarks;
  }
  else if ( _flavTypeToComputeFor == up )
  {
    _ptBinsAll = ptBinsAll_ups;
    _ptBinsDY5 = ptBinsDY5_ups;
    _ptBinsDY4 = ptBinsDY4_ups;
    _ptBinsDY3 = ptBinsDY3_ups;
    _ptBinsDY2 = ptBinsDY2_ups;
    _ptBinsDY1 = ptBinsDY1_ups;
  }
  else if ( _flavTypeToComputeFor == down )
  {
    _ptBinsAll = ptBinsAll_downs;
    _ptBinsDY5 = ptBinsDY5_downs;
    _ptBinsDY4 = ptBinsDY4_downs;
    _ptBinsDY3 = ptBinsDY3_downs;
    _ptBinsDY2 = ptBinsDY2_downs;
    _ptBinsDY1 = ptBinsDY1_downs;
  }
  else if ( _flavTypeToComputeFor == strange )
  {
    _ptBinsAll = ptBinsAll_stranges;
    _ptBinsDY5 = ptBinsDY5_stranges;
    _ptBinsDY4 = ptBinsDY4_stranges;
    _ptBinsDY3 = ptBinsDY3_stranges;
    _ptBinsDY2 = ptBinsDY2_stranges;
    _ptBinsDY1 = ptBinsDY1_stranges;
  }
  else if ( _flavTypeToComputeFor == anti_up )
  {
    _ptBinsAll = ptBinsAll_anti_ups;
    _ptBinsDY5 = ptBinsDY5_anti_ups;
    _ptBinsDY4 = ptBinsDY4_anti_ups;
    _ptBinsDY3 = ptBinsDY3_anti_ups;
    _ptBinsDY2 = ptBinsDY2_anti_ups;
    _ptBinsDY1 = ptBinsDY1_anti_ups;
  }
  else if ( _flavTypeToComputeFor == anti_down )
  {
    _ptBinsAll = ptBinsAll_anti_downs;
    _ptBinsDY5 = ptBinsDY5_anti_downs;
    _ptBinsDY4 = ptBinsDY4_anti_downs;
    _ptBinsDY3 = ptBinsDY3_anti_downs;
    _ptBinsDY2 = ptBinsDY2_anti_downs;
    _ptBinsDY1 = ptBinsDY1_anti_downs;
  }
    else if ( _flavTypeToComputeFor == anti_strange )
  {
    _ptBinsAll = ptBinsAll_anti_stranges;
    _ptBinsDY5 = ptBinsDY5_anti_stranges;
    _ptBinsDY4 = ptBinsDY4_anti_stranges;
    _ptBinsDY3 = ptBinsDY3_anti_stranges;
    _ptBinsDY2 = ptBinsDY2_anti_stranges;
    _ptBinsDY1 = ptBinsDY1_anti_stranges;
  }
  else if ( _flavTypeToComputeFor == allFlavors )
  {
    _ptBinsAll = ptBinsAll_all;
    _ptBinsDY5 = ptBinsDY5_all;
    _ptBinsDY4 = ptBinsDY4_all;
    _ptBinsDY3 = ptBinsDY3_all;
    _ptBinsDY2 = ptBinsDY2_all;
    _ptBinsDY1 = ptBinsDY1_all;
  }
  else
  {
    string errMsg = "error in ptDistribution, flavor not specified";
    throw eAnalysis_error( errMsg );
  }

  string type;
  if ( _flavTypeToComputeFor == gluon )
  {
    type = "gluon";
  }
  else if ( _flavTypeToComputeFor == light_quark || _flavTypeToComputeFor == anti_light_quark )
  {
    type = "quark";
  }
  if ( _flavTypeToComputeFor == up )
  {
    type = "up";
  }
  if ( _flavTypeToComputeFor == down )
  {
    type = "down";
  }
  if ( _flavTypeToComputeFor == strange )
  {
    type = "strange";
  }
  if ( _flavTypeToComputeFor == anti_up )
  {
    type = "anti_up";
  }
  if ( _flavTypeToComputeFor == anti_down )
  {
    type = "anti_down";
  }
  if ( _flavTypeToComputeFor == anti_strange )
  {
    type = "anti_strange";
  }
  else if ( _flavTypeToComputeFor == allFlavors )
  {
    type = "allFlavors";
  }

  string filename = filename_prefix + "_" + type + "_spectra.f2";
  fstream file( filename.c_str(), ios::out | ios::trunc );

  //---- print header if file is empty ----
  file.seekp( 0, ios::end );
  long size = file.tellp();
  file.seekp( 0, ios::beg );
  if ( size == 0 )
    printHeader( file, ptSpectrum, end );
  //---------------------------------------

  //---------------------- y in [-0.5,0.5] ---------------------
  file << "#y in [-0.5,0.5]" << endl;
  for ( int i = 0; i < numberBinsPT; i++ )
  {
    file << ptBinLabels[i] << sep;
    for ( int j = 0; j <= nTimeSteps; j++ )
    {
      if ( j == 0 || j == nTimeSteps )
        file << _ptBinsDY1[j][i] << sep;
      else if ( tstep[j-1] <= theConfig->getRuntime() )
        file << _ptBinsDY1[j][i] << sep;
    }
    file << endl;
  }
  //------------------------------------------------------------
  file << endl << endl;

  //---------------------- y in [-1.0,1.0] ---------------------
  file << "#y in [-1.0,1.0]" << endl;
  for ( int i = 0;i < numberBinsPT;i++ )
  {
    file << ptBinLabels[i] << sep;
    for ( int j = 0; j <= nTimeSteps; j++ )
    {
      if ( j == 0 || j == nTimeSteps )
        file << _ptBinsDY2[j][i] << sep;
      else if ( tstep[j-1] <= theConfig->getRuntime() )
        file << _ptBinsDY2[j][i] << sep;
    }
    file << endl;
  }
  //------------------------------------------------------------
  file << endl << endl;


  //---------------------- y in [-1.5,1.5] ---------------------
  file << "#y in [-1.5,1.5]" << endl;
  for ( int i = 0; i < numberBinsPT; i++ )
  {
    file << ptBinLabels[i] << sep;
    for ( int j = 0; j <= nTimeSteps; j++ )
    {
      if ( j == 0 || j == nTimeSteps )
        file << _ptBinsDY3[j][i] << sep;
      else if ( tstep[j-1] <= theConfig->getRuntime() )
        file << _ptBinsDY3[j][i] << sep;
    }
    file << endl;
  }
  //------------------------------------------------------------
  file << endl << endl;

  //---------------------- y in [-2.0,2.0] ---------------------
  file << "#y in [-2.0,2.0]" << endl;
  for ( int i = 0; i < numberBinsPT; i++ )
  {
    file << ptBinLabels[i] << sep;
    for ( int j = 0; j <= nTimeSteps; j++ )
    {
      if ( j == 0 || j == nTimeSteps )
        file << _ptBinsDY4[j][i] << sep;
      else if ( tstep[j-1] <= theConfig->getRuntime() )
        file << _ptBinsDY4[j][i] << sep;
    }
    file << endl;
  }
  //------------------------------------------------------------
  file << endl << endl;

  //---------------------- y in [-2.5,2.5] ---------------------
  file << "#y in [-2.5,2.5]" << endl;
  for ( int i = 0;i < numberBinsPT;i++ )
  {
    file << ptBinLabels[i] << sep;
    for ( int j = 0; j <= nTimeSteps; j++ )
    {
      if ( j == 0 || j == nTimeSteps )
        file << _ptBinsDY5[j][i] << sep;
      else if ( tstep[j-1] <= theConfig->getRuntime() )
        file << _ptBinsDY5[j][i] << sep;
    }
    file << endl;
  }
  //------------------------------------------------------------
  file << endl << endl;

  //---------------------- y in [-inf,inf] ---------------------
  file << "#y in [-inf,inf]" << endl;
  for ( int i = 0;i < numberBinsPT;i++ )
  {
    file << ptBinLabels[i] << sep;
    for ( int j = 0; j <= nTimeSteps; j++ )
    {
      if ( j == 0 || j == nTimeSteps )
        file << _ptBinsAll[j][i] << sep;
      else if ( tstep[j-1] <= theConfig->getRuntime() )
        file << _ptBinsAll[j][i] << sep;
    }
    file << endl;
  }
  //------------------------------------------------------------

  file.close();
  //---------------------------------------------------------------

}



void analysis::printSoftPtSpectra( const FLAVOR_TYPE _flavTypeToComputeFor )
{
  time_t end;
  time( &end );

  vector<double> * _ptBinsSoftAll;

  if ( _flavTypeToComputeFor == gluon )
  {
    _ptBinsSoftAll = ptBinsSoftAll_gluons;
  }
  else if ( _flavTypeToComputeFor == light_quark )
  {
    _ptBinsSoftAll = ptBinsSoftAll_quarks;
  }
  else if ( _flavTypeToComputeFor == up )
  {
    _ptBinsSoftAll = ptBinsSoftAll_ups;
  }
  else if ( _flavTypeToComputeFor == down )
  {
    _ptBinsSoftAll = ptBinsSoftAll_downs;
  }
  else if ( _flavTypeToComputeFor == strange )
  {
    _ptBinsSoftAll = ptBinsSoftAll_stranges;
  }
  else if ( _flavTypeToComputeFor == anti_up )
  {
    _ptBinsSoftAll = ptBinsSoftAll_anti_ups;
  }
  else if ( _flavTypeToComputeFor == anti_down )
  {
    _ptBinsSoftAll = ptBinsSoftAll_anti_downs;
  }
  else if ( _flavTypeToComputeFor == anti_strange )
  {
    _ptBinsSoftAll = ptBinsSoftAll_anti_stranges;
  }
  else if ( _flavTypeToComputeFor == allFlavors )
  {
    _ptBinsSoftAll = ptBinsSoftAll_all;
  }
  else
  {
    string errMsg = "error in printSoftPtSpectra, flavor not specified";
    throw eAnalysis_error( errMsg );
  }

  string type;
  if ( _flavTypeToComputeFor == gluon )
  {
    type = "gluon";
  }
  else if ( _flavTypeToComputeFor == light_quark || _flavTypeToComputeFor == anti_light_quark )
  {
    type = "quark";
  }
  if ( _flavTypeToComputeFor == up )
  {
    type = "up";
  }
  if ( _flavTypeToComputeFor == down )
  {
    type = "down";
  }
  if ( _flavTypeToComputeFor == strange )
  {
    type = "strange";
  }
  if ( _flavTypeToComputeFor == anti_up )
  {
    type = "anti_up";
  }
  if ( _flavTypeToComputeFor == anti_down )
  {
    type = "anti_down";
  }
  if ( _flavTypeToComputeFor == anti_strange )
  {
    type = "anti_strange";
  }
  else if ( _flavTypeToComputeFor == allFlavors )
  {
    type = "allFlavors";
  }

  string filename = filename_prefix + "_" + type + "_soft_spectra.f2";
  fstream file( filename.c_str(), ios::out | ios::trunc );

  //---- print header if file is empty ----
  file.seekp( 0, ios::end );
  long size = file.tellp();
  file.seekp( 0, ios::beg );
  if ( size == 0 )
    printHeader( file, ptSpectrum, end );
  //---------------------------------------


  //---------------------- y in [-inf,inf] ---------------------
  file << "#y in [-inf,inf]" << endl;
  for ( int i = 0;i < numberBinsSoftPT;i++ )
  {
    file << ptSoftBinLabels[i] << sep;
    for ( int j = 0; j <= nTimeSteps; j++ )
    {
      if ( j == 0 || j == nTimeSteps )
        file << _ptBinsSoftAll[j][i] << sep;
      else if ( tstep[j-1] <= theConfig->getRuntime() )
        file << _ptBinsSoftAll[j][i] << sep;
    }
    file << endl;
  }
  //------------------------------------------------------------

  file.close();
  //---------------------------------------------------------------

}




void analysis::onePartclCorrelations()
{
  const double eta_max = 1.0;
//  const double eta_max = 10000.0;//

  string filename;
  double pp, pt_fin, pt_init, eta, dpt_vec, dpt_scal, rt_init, dE;


  filename = filename_prefix + "_ptIni_ptFin_1correl";
  binning2d ptIniptFin( filename, 0.0, 11.1, 50, 0.0, 11.1, 50 );

  filename = filename_prefix + "_ptIni_dpt_1correl";
  binning2d ptInidpt( filename, 0.0, 11.1, 50, 0.0, 11.1, 50 );

  filename = filename_prefix + "_posIni_dpt_1correl";
  binning2d posInidpt( filename, 0.0, 8.6, 50, 0.0, 10.1, 70 );

  filename = filename_prefix + "_dx_dpt_1correl";
  binning2d dxdpt( filename, 0.0, 6.1, 50, 0.0, 8.1, 70 );

  filename = filename_prefix + "_dx_dE_1correl";
  binning2d dxdE( filename, 0.0, 6.1, 50, -2.1, 8.1, 70 );

  filename = filename_prefix + "_dx_dE_precise_1correl";
  binning2d dxdE_precise( filename, 0.0, 6.1, 50, 0.0, 0.21, 50 );

  filename = filename_prefix + "_posIni_ptIni_ptFinAv_1correl";
  binningValues2d posIniptIniptFin( filename, 0.0, 8.1, 35, 0.0, 15.0, 50 );

  filename = filename_prefix + "_posIni_ptIni_dptAv_1correl";
  binningValues2d posIniptInidpt( filename, 0.0, 8.1, 35, 0.0, 15.0, 50 );

//   filename = filename_prefix + "_posIni_ptIni_dptptIAv_1correl";
  filename = "/dev/null";
  binningValues2d posIniptInidptptI( filename, 0.0, 8.1, 35, 0.0, 15.0, 50 );

  filename = filename_prefix + "_ptFin_posIni";
  binningValues ptFinposIni( filename, 0.0, 15.0, 30 );
  filename = filename_prefix + "_ptIni_posIni";
  binningValues ptIniposIni( filename, 0.0, 15.0, 30 );


//   filename = filename_prefix + "_ptFin_ptFin_test";
  filename = "/dev/null";
  binningValues2d ptFinptFin( filename, 0.0, 11.1, 50, 0.0, 11.1, 50 );



  double q_hat_scal_sum = 0.0, q_hat_vec_sum = 0.0;
  int count_q_hat_sum = 0;

  for ( int i = 1;i <= addedParticles.size();i++ )
  {
//     if(addedParticles[i].T <= time)
//     {
    pp = sqrt( pow( addedParticles[i].E, 2.0 ) - pow( addedParticles[i].m, 2.0 ) );

    // pseudorapidity
    eta = 0.5 * log(( pp + addedParticles[i].PZ ) / ( pp - addedParticles[i].PZ ) );

    if ( fabs( eta ) <= eta_max )
    {
      pt_init = sqrt( pow( addedParticles[i].PX_init, 2.0 ) + pow( addedParticles[i].PY_init, 2.0 ) );
      pt_fin = sqrt( pow( addedParticles[i].PX, 2.0 ) + pow( addedParticles[i].PY, 2.0 ) );
      dpt_vec = sqrt( pow( addedParticles[i].PX_init - addedParticles[i].PX, 2.0 ) + pow( addedParticles[i].PY_init - addedParticles[i].PY, 2.0 ) );
      dpt_scal = sqrt( pow( addedParticles[i].PX_init, 2.0 ) + pow( addedParticles[i].PY_init, 2.0 ) ) - sqrt( pow( addedParticles[i].PX, 2.0 ) + pow( addedParticles[i].PY, 2.0 ) );
      rt_init = sqrt( pow( addedParticles[i].X_init, 2.0 ) + pow( addedParticles[i].Y_init, 2.0 ) );
      dE = addedParticles[i].E_init - addedParticles[i].E;

      ptIniptFin.add( pt_init, pt_fin );
      ptInidpt.add( pt_init, dpt_vec );
      posInidpt.add( rt_init, dpt_vec );
      dxdpt.add( addedParticles[i].X_traveled, dpt_vec );
      dxdE.add( addedParticles[i].X_traveled, dE );
      dxdE_precise.add( addedParticles[i].X_traveled, dE );

      posIniptIniptFin.add( rt_init, pt_init, pt_fin );
      posIniptInidpt.add( rt_init, pt_init, dpt_scal );
      posIniptInidptptI.add( rt_init, pt_init, dpt_scal / pt_init );
      ptFinposIni.add( pt_fin, rt_init );
      ptIniposIni.add( pt_init, rt_init );

      ptFinptFin.add( pt_init, pt_fin, dpt_scal );

      if ( addedParticles[i].X_traveled > 0.0 )
      {
        q_hat_scal_sum += pow( dpt_scal, 2.0 ) / addedParticles[i].X_traveled;
        q_hat_vec_sum += pow( dpt_vec, 2.0 ) / addedParticles[i].X_traveled;
        count_q_hat_sum++;
      }
    }
//     }
  }

  ptIniptFin.print();
  ptInidpt.print();
  posInidpt.print();
  dxdpt.print();
  dxdE.print();
  dxdE_precise.print();
  posIniptIniptFin.print();
  posIniptInidpt.print();
  posIniptInidptptI.print();
  ptFinposIni.print();
  ptIniposIni.print();

  ptFinptFin.print();

  filename = filename_prefix + "_q_hat";
  fstream print_q_hat( filename.c_str(), ios::out | ios::trunc );

//   cout << "Final output: q_hat = " << q_hat_scal_sum/count_q_hat_sum << endl;
  print_q_hat << q_hat_scal_sum / count_q_hat_sum << endl << q_hat_vec_sum / count_q_hat_sum << endl;

}

void analysis::twoPartclCorrelations()
{
  const double eta_max = 1.0;
//  const double eta_max = 10000.0; //!

  const double pt_min_trig = 4.0; //GeV, cut on trigger particle, which is considered for angle correlations dphi
  const double pt_min_assoc = 2.0; //GeV, cut on associated particle, which is considered for angle correlations dphi
  const double pt_min_dNdr = 5.0;//GeV, cut on particle for dN/dr_t

  string filename;
  double pp, eta, rt_init;
  double dpt_ini, dpt_fin, dp_ini, dp_fin, dphi_fin, dphi_ini, scal_prod, length1, length2, pt_fin_i, pt_fin_j;

  filename = filename_prefix + "_dpt_partons";
  fstream partdpt( filename.c_str(), ios::out | ios::app );


  filename = filename_prefix + "_dptIni_dptFin_2correl";
  binning2d dptIdptF2( filename, 0.0, 8.1, 60, 0.0, 10.5, 70 );

  filename = filename_prefix + "_dptIni_dptFin_posIniAv_2correl";
  binningValues2d dptIdptFposI2( filename, 0.0, 8.1, 60, 0.0, 10.5, 70 );

//   filename = filename_prefix + "_dpIni_dpFin_2correl";
//   binning2d dpIdpF2(filename, 0.0, 10.5, 40, 0.0, 10.5, 40);

  filename = filename_prefix + "_posIni_dptFin_2correl";
  binning2d posIdptF2( filename, 0.0, 8.1, 35, 0.0, 10.5, 80 );

  filename = filename_prefix + "_posIni_dptIni_2correl";
  binning2d posIdptI2( filename, 0.0, 8.1, 35, 0.0, 10.5, 80 );



  filename = filename_prefix + "_dphiFin_4_2_2correl";
  binning dphiF42( filename, 0.0, M_PI, 80 );
  filename = filename_prefix + "_dphiIni_4_2_2correl";
  binning dphiI42( filename, 0.0, M_PI, 80 );
  filename = filename_prefix + "_dphiFin_2_2_2correl";
  binning dphiF22( filename, 0.0, M_PI, 80 );
  filename = filename_prefix + "_dphiIni_2_2_2correl";
  binning dphiI22( filename, 0.0, M_PI, 80 );
  filename = filename_prefix + "_dphiFin_all_2correl";
  binning dphiFall( filename, 0.0, M_PI, 80 );
  filename = filename_prefix + "_dphiIni_all_2correl";
  binning dphiIall( filename, 0.0, M_PI, 80 );
  filename = filename_prefix + "_dphiFin_rt3_2correl";
  binning dphiFrt3( filename, 0.0, M_PI, 80 );
  filename = filename_prefix + "_dphiIni_rt3_2correl";
  binning dphiIrt3( filename, 0.0, M_PI, 80 );

  filename = filename_prefix + "_dphiFin_posIni_all_2correl";
  binning2d dphiFposIniall( filename, 0.0, M_PI, 50, 0.0, 8.1, 35 );
  filename = filename_prefix + "_dphiIni_posIni_all_2correl";
  binning2d dphiIposIniall( filename, 0.0, M_PI, 50, 0.0, 8.1, 35 );

  filename = filename_prefix + "_dphi_iso";
  binning dphiiso( filename, 0.0, M_PI, 80 );

  filename = filename_prefix + "_dN_drInit_all";
  binning dNdrI_all( filename, 0.0, 15.0, 60 );
  filename = filename_prefix + "_dN_drInit_pt5";
  binning dNdrI_pt5( filename, 0.0, 15.0, 60 );
  filename = filename_prefix + "_dN_drInit_ptbar5";
  binning dNdrI_ptbar5( filename, 0.0, 15.0, 60 );
  filename = filename_prefix + "_dN_drInit_ptbar5_dpt02";
  binning dNdrI_ptbar5_dpt02( filename, 0.0, 15.0, 60 );

  // 2 particle correlations, charm and anti-charm which are produced in same reaction
  for ( int i = 1;i <= addedParticles.size();i++ )
  {
//     if(addedParticles[i].T <= time)
//     {
    pp = sqrt( pow( addedParticles[i].E, 2.0 ) - pow( addedParticles[i].m, 2.0 ) );

    // pseudorapidity
    eta = 0.5 * log(( pp + addedParticles[i].PZ ) / ( pp - addedParticles[i].PZ ) );

    if ( fabs( eta ) <= eta_max )
    {
      for ( int j = i + 1;j <= addedParticles.size();j++ )
      {
//         if ( addedParticles[i].N_EVENT == addedParticles[j].N_EVENT )
        if ( true )
        {

          //           if(addedParticles[i].T <= time)
          //           {
          pp = sqrt( pow( addedParticles[j].E, 2.0 ) - pow( addedParticles[j].m, 2.0 ) );

          // pseudorapidity
          eta = 0.5 * log(( pp + addedParticles[j].PZ ) / ( pp - addedParticles[j].PZ ) );

          if ( fabs( eta ) <= eta_max )
          {
//                 dpt_ini = sqrt( pow(addedParticles[i].PX_init-addedParticles[j].PX_init,2.0) + pow(addedParticles[i].PY_init-addedParticles[j].PY_init,2.0) );
//                 dpt_fin = sqrt( pow(addedParticles[i].PX-addedParticles[j].PX,2.0) + pow(addedParticles[i].PY-addedParticles[j].PY,2.0) );
            dpt_ini = fabs( sqrt( pow( addedParticles[i].PX_init, 2.0 ) + pow( addedParticles[i].PY_init, 2.0 ) ) - sqrt( pow( addedParticles[j].PX_init, 2.0 ) + pow( addedParticles[j].PY_init, 2.0 ) ) );
            dpt_fin = fabs( sqrt( pow( addedParticles[i].PX, 2.0 ) + pow( addedParticles[i].PY, 2.0 ) ) - sqrt( pow( addedParticles[j].PX, 2.0 ) + pow( addedParticles[j].PY, 2.0 ) ) );

            dp_ini = sqrt( pow( addedParticles[i].PX_init - addedParticles[j].PX_init, 2.0 ) + pow( addedParticles[i].PY_init - addedParticles[j].PY_init, 2.0 ) + pow( addedParticles[i].PZ_init - addedParticles[j].PZ_init, 2.0 ) );
            dp_fin = sqrt( pow( addedParticles[i].PX - addedParticles[j].PX, 2.0 ) + pow( addedParticles[i].PY - addedParticles[j].PY, 2.0 ) + pow( addedParticles[i].PZ - addedParticles[j].PZ, 2.0 ) );
            rt_init = sqrt( pow( addedParticles[i].X_init, 2.0 ) + pow( addedParticles[i].Y_init, 2.0 ) ); // same for both particles

            dptIdptF2.add( dpt_ini, dpt_fin );
            dptIdptFposI2.add( dpt_ini, dpt_fin, rt_init );
//                 dpIdpF2.add(dp_ini,dp_fin);
            posIdptF2.add( rt_init, dpt_fin );
            posIdptI2.add( rt_init, dpt_ini );



//                 if(dpt_fin > 7.0)
//                 {
//                   partdpt.width(15);
//                   partdpt << dpt_fin;
//                   partdpt.width(15);
//                   partdpt << fabs( sqrt( pow(addedParticles[i].PX,2.0) + pow(addedParticles[i].PY,2.0) ) - sqrt( pow(addedParticles[j].PX,2.0) + pow(addedParticles[j].PY,2.0) ) );
//
//                   partdpt.width(25);
//                   partdpt << addedParticles[i].PX;
//                   partdpt.width(15);
//                   partdpt << addedParticles[i].PY;
//                   partdpt.width(15);
//                   partdpt << addedParticles[i].PZ;
//
//                   partdpt.width(25);
//                   partdpt << addedParticles[j].PX;
//                   partdpt.width(15);
//                   partdpt << addedParticles[j].PY;
//                   partdpt.width(15);
//                   partdpt << addedParticles[j].PZ;
//
//                   partdpt.width(25);
//                   partdpt << addedParticles[i].PX_init;
//                   partdpt.width(15);
//                   partdpt << addedParticles[i].PY_init;
//                   partdpt.width(15);
//                   partdpt << addedParticles[i].PZ_init;
//
//                   partdpt.width(25);
//                   partdpt << addedParticles[j].PX_init;
//                   partdpt.width(15);
//                   partdpt << addedParticles[j].PY_init;
//                   partdpt.width(15);
//                   partdpt << addedParticles[j].PZ_init;
//
//                   partdpt.width(25);
//                   partdpt << addedParticles[j].X_init;
//                   partdpt.width(15);
//                   partdpt << addedParticles[j].Y_init;
//                   partdpt.width(15);
//                   partdpt << addedParticles[j].Z_init;
//
//                   partdpt.width(25);
//                   partdpt << rt_init << endl;
//
//                 }


            // azimuthal angle between two charm quarks
            scal_prod = addedParticles[i].PX_init * addedParticles[j].PX_init + addedParticles[i].PY_init * addedParticles[j].PY_init;
            length1 = sqrt( pow( addedParticles[i].PX_init, 2.0 ) + pow( addedParticles[i].PY_init, 2.0 ) );
            length2 = sqrt( pow( addedParticles[j].PX_init, 2.0 ) + pow( addedParticles[j].PY_init, 2.0 ) );
            dphi_ini = acos( scal_prod / length1 / length2 );

            scal_prod = addedParticles[i].PX * addedParticles[j].PX + addedParticles[i].PY * addedParticles[j].PY;
            length1 = sqrt( pow( addedParticles[i].PX, 2.0 ) + pow( addedParticles[i].PY, 2.0 ) );
            length2 = sqrt( pow( addedParticles[j].PX, 2.0 ) + pow( addedParticles[j].PY, 2.0 ) );
            dphi_fin = acos( scal_prod / length1 / length2 );

            // total angle between two charm quarks
//                 scal_prod = addedParticles[i].PX_init*addedParticles[j].PX_init + addedParticles[i].PY_init*addedParticles[j].PY_init + addedParticles[i].PZ_init*addedParticles[j].PZ_init;
//                 length1 = sqrt( pow(addedParticles[i].PX_init,2.0) + pow(addedParticles[i].PY_init,2.0) + pow(addedParticles[i].PZ_init,2.0)  );
//                 length2 = sqrt( pow(addedParticles[j].PX_init,2.0) + pow(addedParticles[j].PY_init,2.0) + pow(addedParticles[j].PZ_init,2.0)  );
//                 dphi_ini = acos( scal_prod / length1 / length2 );
//
//                 scal_prod = addedParticles[i].PX*addedParticles[j].PX + addedParticles[i].PY*addedParticles[j].PY + addedParticles[i].PZ*addedParticles[j].PZ;
//                 length1 = sqrt( pow(addedParticles[i].PX,2.0) + pow(addedParticles[i].PY,2.0) + pow(addedParticles[i].PZ,2.0)  );
//                 length2 = sqrt( pow(addedParticles[j].PX,2.0) + pow(addedParticles[j].PY,2.0) + pow(addedParticles[j].PZ,2.0)  );
//                 dphi_fin = acos( scal_prod / length1 / length2 );

//                 dphi_fin = cos(dphi_fin);
//                 dphi_ini = cos(dphi_ini);

            dphiFall.add( dphi_fin );
            dphiIall.add( dphi_ini );

            if ( sqrt( pow( addedParticles[i].PX, 2.0 ) + pow( addedParticles[i].PY, 2.0 ) ) >= pt_min_assoc && sqrt( pow( addedParticles[j].PX, 2.0 ) + pow( addedParticles[j].PY, 2.0 ) )  >= pt_min_assoc )
            {
              dphiF22.add( dphi_fin );
              dphiI22.add( dphi_ini );

              if (( sqrt( pow( addedParticles[i].PX, 2.0 ) + pow( addedParticles[i].PY, 2.0 ) ) >= pt_min_trig && sqrt( pow( addedParticles[j].PX, 2.0 ) + pow( addedParticles[j].PY, 2.0 ) )  >= pt_min_assoc ) ||
                  ( sqrt( pow( addedParticles[i].PX, 2.0 ) + pow( addedParticles[i].PY, 2.0 ) ) >= pt_min_assoc && sqrt( pow( addedParticles[j].PX, 2.0 ) + pow( addedParticles[j].PY, 2.0 ) )  >= pt_min_trig ) )
              {
                dphiF42.add( dphi_fin );
                dphiI42.add( dphi_ini );
              }
            }

            if ( rt_init )
            {
              dphiFrt3.add( dphi_fin );
              dphiIrt3.add( dphi_ini );
            }

            dphiFposIniall.add( dphi_fin, rt_init );
            dphiIposIniall.add( dphi_ini, rt_init );



            pt_fin_i = sqrt( pow( addedParticles[i].PX, 2.0 ) + pow( addedParticles[i].PY, 2.0 ) );
            pt_fin_j = sqrt( pow( addedParticles[j].PX, 2.0 ) + pow( addedParticles[j].PY, 2.0 ) );

            dNdrI_all.add( rt_init );
            if ( pt_fin_i > pt_min_dNdr )
              dNdrI_pt5.add( rt_init );
            if (( pt_fin_i + pt_fin_j ) / 2.0 > pt_min_dNdr )
              dNdrI_ptbar5.add( rt_init );
            if ((( pt_fin_i + pt_fin_j ) / 2.0 > pt_min_dNdr ) && ( dpt_fin / (( pt_fin_i + pt_fin_j ) / 2.0 ) < 0.2 ) )
              dNdrI_ptbar5_dpt02.add( rt_init );

          }
          //           }
          break;
        }
      }

    }
//     }
  }



  for ( int j = 1;j <= addedParticles.size();j++ )
  {

    pp = sqrt( pow( addedParticles[j].E, 2.0 ) - pow( addedParticles[j].m, 2.0 ) );

    // pseudorapidity
    eta = 0.5 * log(( pp + addedParticles[j].PZ ) / ( pp - addedParticles[j].PZ ) );

    double x1 = 0.0;
    double y1 = 1.0;

    if ( fabs( eta ) <= eta_max )
    {
      scal_prod = x1 * addedParticles[j].PX + y1 * addedParticles[j].PY;
      length1 = sqrt( pow( x1, 2.0 ) + pow( y1, 2.0 ) );
      length2 = sqrt( pow( addedParticles[j].PX, 2.0 ) + pow( addedParticles[j].PY, 2.0 ) );
      dphi_fin = acos( scal_prod / length1 / length2 );
      dphiiso.add( dphi_fin );
    }
  }

  dptIdptF2.print();
  dptIdptFposI2.print();
//   dpIdpF2.print();
  posIdptF2.print();
  posIdptI2.print();
  dphiF42.print();
  dphiI42.print();
  dphiF22.print();
  dphiI22.print();
  dphiFall.print();
  dphiIall.print();
  dphiiso.print();
  dphiFposIniall.print();
  dphiIposIniall.print();
  dphiFrt3.print();
  dphiIrt3.print();

  dNdrI_all.print();
  dNdrI_pt5.print();
  dNdrI_ptbar5.print();
  dNdrI_ptbar5_dpt02.print();

}




// compute v2 of gluons and charm quarks
void analysis::computeV2RAA( string name, const double _outputTime )
{
  v2RAA theV2RAA( theConfig, name, filename_prefix );

  theV2RAA.computeFor( gluon, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
  theV2RAA.computeFor( up, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
  theV2RAA.computeFor( down, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
  theV2RAA.computeFor( strange, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
  theV2RAA.computeFor( anti_up, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
  theV2RAA.computeFor( anti_down, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
  theV2RAA.computeFor( anti_strange, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
  theV2RAA.computeFor( light_quark, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
  
  theV2RAA.computeFor( gluon, addedParticles, addedParticles.size(), "jets", _outputTime, v2jets );
  theV2RAA.computeFor( light_quark, addedParticles, addedParticles.size(), "jets", _outputTime, v2jets );
  theV2RAA.computeFor( up, addedParticles, addedParticles.size(), "jets", _outputTime, v2jets );
  theV2RAA.computeFor( down, addedParticles, addedParticles.size(), "jets", _outputTime, v2jets );
  theV2RAA.computeFor( strange, addedParticles, addedParticles.size(), "jets", _outputTime, v2jets );
  theV2RAA.computeFor( anti_up, addedParticles, addedParticles.size(), "jets", _outputTime, v2jets );
  theV2RAA.computeFor( anti_down, addedParticles, addedParticles.size(), "jets", _outputTime, v2jets );
  theV2RAA.computeFor( anti_strange, addedParticles, addedParticles.size(), "jets", _outputTime, v2jets );
}



v2RAA::v2RAA( config * const c, string name_arg, string filename_prefix_arg ):
    theConfig( c ), name( name_arg ), filename_prefix( filename_prefix_arg )
{
  pt_min = 5.0;
  pt_max = 40.0;
  n_g = 35;
  
  pt_min_background = 0;
  pt_max_background = 5.0;
  n_g_background = 25;

  eta_bins = 6; // Number of eta bins: 0: 0.35, 1: 0.5, 2: 0.75, 3: 1.0, 4: 1.5, 5: 2.0
}




void v2RAA::computeFor( const FLAVOR_TYPE _flavTypeToComputeFor, vector<ParticleOffline>& _particles, const int n_particles, string additionalNameTag, const double _outputTime, const v2Type _v2type )
{
  double eta, pp, pt, v2, xt;
  double sinAlpha, alpha;
  int dummy, flavor, n_bins;
  int alphaIndex;
  
  double _pt_min, _pt_max;

  string filename_v2, filename_v2_summed, filename_v2_tot, filename_yield, filename_pt_angleDependence, type;

  // Number of bins
  if ( _flavTypeToComputeFor == gluon )
  {
    n_bins = n_g;
    type = "gluon";
  }
  else if ( _flavTypeToComputeFor == light_quark || _flavTypeToComputeFor == anti_light_quark )
  {
    n_bins = n_g;
    type = "quark";
  }
  else if ( _flavTypeToComputeFor == up )
  {
    n_bins = n_g;
    type = "up";    
  }
  else if ( _flavTypeToComputeFor == down )
  {
    n_bins = n_g;
    type = "down";    
  }
  else if ( _flavTypeToComputeFor == strange )
  {
    n_bins = n_g;
    type = "strange";    
  }
  else if ( _flavTypeToComputeFor == anti_up )
  {
    n_bins = n_g;
    type = "anti_up";    
  }
  else if ( _flavTypeToComputeFor == anti_down )
  {
    n_bins = n_g;
    type = "anti_down";    
  }
  else if ( _flavTypeToComputeFor == anti_strange )
  {
    n_bins = n_g;
    type = "anti_strange";    
  }
  
  if ( _v2type == v2background )
  {
    n_bins = n_g_background;
    _pt_max = pt_max_background;
    _pt_min = pt_min_background;
  }
  else
  {
    n_bins = n_g;
    _pt_max = pt_max;
    _pt_min = pt_min;
  }
  
  const double dpt = ( _pt_max - _pt_min ) / n_bins;

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
    pp = sqrt( pow( _particles[i].E, 2.0 ) - pow( _particles[i].m, 2.0 ) );
    pt = sqrt( pow( _particles[i].PX, 2.0 ) + pow( _particles[i].PY, 2.0 ) );
    xt = sqrt( pow( _particles[i].X, 2.0 ) + pow( _particles[i].Y, 2.0 ) );

    sinAlpha = _particles[i].PY / pt;
    alpha = asin( fabs( sinAlpha ) );
    alpha = alpha * 180 / M_PI;
    
    alphaIndex = static_cast<int>( alpha / deltaAlpha );
    if ( alphaIndex > nAlphaBins )
    {
      alphaIndex = nAlphaBins - 1;
    }
        
    // pseudorapidity
    eta = 0.5 * log(( pp + _particles[i].PZ ) / ( pp - _particles[i].PZ ) );

    v2 = ( pow( _particles[i].PX, 2.0 ) - pow( _particles[i].PY, 2.0 ) ) / pow( pt, 2.0 );

    flavor = _particles[i].FLAVOR;

    if (( pt <= _pt_max && pt >= _pt_min ) &&  
      ( _flavTypeToComputeFor == flavor || 
      ( _flavTypeToComputeFor == ParticleOffline::mapToGenericFlavorType( static_cast<FLAVOR_TYPE>( flavor ) ) ) ) )
    {
      if ( fabs( eta ) <= 2.0 )
      {
        v2sum[5] += v2;
        NmbInRange[5]++;

        dummy = int(( pt - _pt_min ) / dpt );
        ptBinsV2[5][dummy] += v2;
        ptBinsNmb[5][dummy]++;
        ptBinsAngleDep[5][alphaIndex][dummy]++;

        if ( fabs( eta ) <= 1.5 )
        {
          v2sum[4] += v2;
          NmbInRange[4]++;

          dummy = int(( pt - _pt_min ) / dpt );
          ptBinsV2[4][dummy] += v2;
          ptBinsNmb[4][dummy]++;
          ptBinsAngleDep[4][alphaIndex][dummy]++;

          if ( fabs( eta ) <= 1.0 )
          {
            v2sum[3] += v2;
            NmbInRange[3]++;

            dummy = int(( pt - _pt_min ) / dpt );
            ptBinsV2[3][dummy] += v2;
            ptBinsNmb[3][dummy]++;
            ptBinsAngleDep[3][alphaIndex][dummy]++;

            if ( fabs( eta ) <= 0.75 )
            {
              v2sum[2] += v2;
              NmbInRange[2]++;

              dummy = int(( pt - _pt_min ) / dpt );
              ptBinsV2[2][dummy] += v2;
              ptBinsNmb[2][dummy]++;
              ptBinsAngleDep[2][alphaIndex][dummy]++;

              if ( fabs( eta ) <= 0.5 )
              {
                if ( xt < 1.5 )
                {
                  NmbInnerRegion++;
                  ptBinsInnerRegion[dummy]++; 
                }
                
                v2sum[1] += v2;
                NmbInRange[1]++;

                dummy = int(( pt - _pt_min ) / dpt );
                ptBinsV2[1][dummy] += v2;
                ptBinsNmb[1][dummy]++;
                ptBinsAngleDep[1][alphaIndex][dummy]++;

                if ( fabs( eta ) <= 0.35 )
                {
                  v2sum[0] += v2;
                  NmbInRange[0]++;

                  dummy = int(( pt - _pt_min ) / dpt );
                  ptBinsV2[0][dummy] += v2;
                  ptBinsNmb[0][dummy]++;
                  ptBinsAngleDep[0][alphaIndex][dummy]++;
                }
              }
            }
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

  filename_v2 = filename_prefix + "_" + type + "_" + additionalNameTag + "_v2_pt_" + name;
  filename_v2_summed = filename_prefix + "_" + type + "_" + additionalNameTag + "_v2_pt_summed_" + name;
  filename_v2_tot = filename_prefix + "_" + type + "_" + additionalNameTag + "_v2_tot_" + name;
  filename_yield = filename_prefix + "_" + type + "_" + additionalNameTag + "_yield_pt_" + name;
  filename_pt_angleDependence = filename_prefix + "_" + type + "_" + additionalNameTag + "_pt_angular_dependence_" + name;

  fstream print_v2( filename_v2.c_str(), ios::out | ios::trunc );
  fstream print_v2_summed( filename_v2_summed.c_str(), ios::out | ios::trunc );
  fstream print_v2_tot( filename_v2_tot.c_str(), ios::out | ios::trunc );
  fstream print_yield( filename_yield.c_str(), ios::out | ios::trunc );
  fstream print_pt_angleDependence( filename_pt_angleDependence.c_str(), ios::out | ios::trunc );

//   cout << "total v2 of " << type << " in Y=+-0.35 = " << v2sum[0]/NmbInRange[0] << endl;

  // print total v2
  print_v2_tot << "# total v2 of " << type << endl;
  print_v2_tot << "# t = " << _outputTime << endl;
  print_v2_tot << "# bin statistics for 0.35 mid-rapidity:  Avg per bin=" << double( NmbInRange[0] ) / n_c << "   Min=" << binMin << "   Max=" << binMax << endl;
  print_v2_tot << "#mid.rap. betw. +-Y for  Y=0.35  Y=0.5  Y=0.75  Y=1  Y=1.5  Y=2   v2 sum and number in range  for mid.rap. betw. +-Y for  Y=0.35  Y=0.5  Y=0.75  Y=1  Y=1.5  Y=2" << endl;

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


  // print v2 distribution as function of pt
  print_v2 << "#v2 of " << type << endl;
  print_v2 << "# t = " << _outputTime << endl;
//   print_v2 << "#simulation parameter:" << endl;
//   print_v2 << "#testparticles= " << theConfig->getTestparticles() << endl;
//   print_v2 << "#runtime= " << theConfig->getRuntime() << endl;
//   print_v2 << "#sqrtS= " << theConfig->getSqrtS() << " GeV" << endl;
//   print_v2 << "#b= " << theConfig->getImpactParameter() << " fm" << endl;
//   print_v2 << "#cascade data input folder: " << theConfig->getPathdirCascadeData() << endl;
//   print_v2 << "#data file for charm quarks: " << theConfig->getPythiaParticleFileCharm() << endl;
//   print_v2 << "#initial seed for random generator ran2(): " << initialSeed << endl;
//   print_v2 << "#" << endl;
  print_v2 << "#";
  print_v2.width( 14 );
  print_v2 << "pt";
  print_v2.width( 15 );
  print_v2 << "v_2 for mid.rap. betw. +-Y for  Y=0.35  Y=0.5  Y=0.75  Y=1  Y=1.5  Y=2" << endl;

  for ( int k = 0;k < n_bins + 1;k++ )
  {
    print_v2.width( 15 );
    pt_out = double( k ) * dpt + _pt_min + dpt / 2.0; // +dpt/2.0 to shift value in the middle of the intervall, important for small number of bins
    print_v2 << pt_out;
    for ( int i = 0;i < eta_bins;i++ )
    {
      print_v2.width( 15 );
      if ( ptBinsNmb[i][k] > 0 )
      {
        print_v2 << ptBinsV2[i][k] / ptBinsNmb[i][k];
      }
      else
      {
        print_v2 << 0;
      }
    }
    print_v2 << endl;
  }


  // print summed output, v2 is not computed, but summed v2 and the number in one bin
  print_v2_summed << "# summed v2 of " << type << endl;
  print_v2_summed << "# t = " << _outputTime << endl;
  print_v2_summed << "#";
  print_v2_summed.width( 14 );
  print_v2_summed << "pt";
  print_v2_summed.width( 15 );
  print_v2_summed << "summed v_2 and number in bin for mid.rap. betw. +-Y for  Y=0.35  Y=0.5  Y=0.75  Y=1  Y=1.5  Y=2" << endl;

  for ( int k = 0;k < n_bins + 1;k++ )
  {
    print_v2_summed.width( 15 );
    pt_out = double( k ) * dpt + _pt_min + dpt / 2.0; // +dpt/2.0 to shift value in the middle of the intervall, important for small number of bins
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
  print_yield << "# t = " << _outputTime << endl;
  print_yield << "#";
  print_yield.width( 14 );
  print_yield << "pt";
  print_yield.width( 15 );
  print_yield << "yield for mid.rap. betw. +-Y for  Y=0.35  Y=0.5  Y=0.75  Y=1  Y=1.5  Y=2   central_region (|y|<0.5, r < 1.5)" << endl;

  for ( int k = 0;k < n_bins + 1;k++ )
  {
    print_yield.width( 15 );
    pt_out = double( k ) * dpt + _pt_min + dpt / 2.0; // +dpt/2.0 to shift value in the middle of the intervall, important for small number of bins
    print_yield << pt_out;
    for ( int i = 0;i < eta_bins;i++ )
    {
      print_yield.width( 15 );
//       if ( particleType == charm || particleType == bottom || particleType == charmBottom )
//         print_yield << double( ptBinsNmb[i][k] ) / theConfig->getTestparticles() / theConfig->getKIniCharm();
//       else
      print_yield << double( ptBinsNmb[i][k] ) / theConfig->getTestparticles();
    }
    print_yield.width( 15 );
    print_yield << double( ptBinsInnerRegion[k] ) / theConfig->getTestparticles();    
    print_yield << endl;
  }
  
  
  // print yield for RAA for different angles with respect to the reaction plane
  print_pt_angleDependence << "# " << type << " yield distribution for different angles (alpha) with respect to the reaction plane" << endl;
  print_pt_angleDependence << "# t = " << _outputTime << endl;
  print_pt_angleDependence << "#";
  print_pt_angleDependence.width( 14 );
  print_pt_angleDependence << "pt";
  print_pt_angleDependence.width( 15 );
  print_pt_angleDependence << "yield for mid.rap. betw. +-Y for  Y=0.35  Y=0.5  Y=0.75  Y=1  Y=1.5  Y=2" << endl;
  for ( int j = 0; j < nAlphaBins; j++ )
  {
    print_pt_angleDependence << "#alpha in [ " << j * deltaAlpha << "°, " << (j+1)*deltaAlpha << "° ] "<< endl;
    for ( int k = 0;k < n_bins + 1;k++ )
    {
      print_pt_angleDependence.width( 15 );
      pt_out = double( k ) * dpt + _pt_min + dpt / 2.0; // +dpt/2.0 to shift value in the middle of the intervall, important for small number of bins
      print_pt_angleDependence << pt_out;
      for ( int i = 0;i < eta_bins;i++ )
      {
        print_pt_angleDependence.width( 15 );
        print_pt_angleDependence << double( ptBinsAngleDep[i][j][k] ) / theConfig->getTestparticles();
      }
      print_pt_angleDependence << endl;
    }    
    print_pt_angleDependence << endl;
    print_pt_angleDependence << endl;
  }
  
  

//   cout << type << ": " << NmbInRange[5] << endl;
}





void analysis::volumeMidrap( const int step ) const
{
  string name;
  stringstream ss;

  ss << step;
  name = "step" + ss.str();


  string filename_z = filename_prefix + "_zMidrap_" + name;
  string filename_t = filename_prefix + "_tMidrap_" + name;

  double y, xt;
  const double y_max = 0.5;

  int numberOfBins = 80;

  binning zBins( filename_z, numberOfBins );
  binning tBins( filename_t, numberOfBins );

  for ( int i = 0;i < particles_atTimeNow.size();i++ )
  {
//     y = 0.5*log( (partclAtTimeNow[i].E + partclAtTimeNow[i].PZ) / (partclAtTimeNow[i].E - partclAtTimeNow[i].PZ) );
//     if(fabs(y) <= y_max)
//     {
//       zBins.add( fabs(partclAtTimeNow[i].Z) );
//       xt = sqrt(pow(partclAtTimeNow[i].X,2.0) + pow(partclAtTimeNow[i].Y,2.0));
//       tBins.add(xt);
//     }

    if ( particles_atTimeNow[i].T <= tstep[step] )
    {

      y = 0.5 * log(( particles_atTimeNow[i].E + particles_atTimeNow[i].PZ ) / ( particles_atTimeNow[i].E - particles_atTimeNow[i].PZ ) );
      if ( fabs( y ) <= y_max )
      {
        zBins.add( fabs( particles_atTimeNow[i].Z ) );
        xt = sqrt( pow( particles_atTimeNow[i].X, 2.0 ) + pow( particles_atTimeNow[i].Y, 2.0 ) );
        tBins.add( xt );
      }



//       zBins.add( fabs(particles_atTimeNow[i].Z) );
//       xt = sqrt(pow(particles_atTimeNow[i].X,2.0) + pow(particles_atTimeNow[i].Y,2.0));
//       tBins.add(xt);
    }

  }

  zBins.print();
  tBins.print();
}


void analysis::ptDistribution( const FLAVOR_TYPE _flavTypeToComputeFor, vector<ParticleOffline>& _particles, const int n_particles, const int step )
{
  double pt, y;
  vector<double> * _ptBinsAll;
  vector<double> * _ptBinsDY5;
  vector<double> * _ptBinsDY4;
  vector<double> * _ptBinsDY3;
  vector<double> * _ptBinsDY2;
  vector<double> * _ptBinsDY1;

  if ( _flavTypeToComputeFor == gluon )
  {
    _ptBinsAll = ptBinsAll_gluons;
    _ptBinsDY5 = ptBinsDY5_gluons;
    _ptBinsDY4 = ptBinsDY4_gluons;
    _ptBinsDY3 = ptBinsDY3_gluons;
    _ptBinsDY2 = ptBinsDY2_gluons;
    _ptBinsDY1 = ptBinsDY1_gluons;
  }
  else if ( _flavTypeToComputeFor == light_quark )
  {
    _ptBinsAll = ptBinsAll_quarks;
    _ptBinsDY5 = ptBinsDY5_quarks;
    _ptBinsDY4 = ptBinsDY4_quarks;
    _ptBinsDY3 = ptBinsDY3_quarks;
    _ptBinsDY2 = ptBinsDY2_quarks;
    _ptBinsDY1 = ptBinsDY1_quarks;
  }
  else if ( _flavTypeToComputeFor == up )
  {
    _ptBinsAll = ptBinsAll_ups;
    _ptBinsDY5 = ptBinsDY5_ups;
    _ptBinsDY4 = ptBinsDY4_ups;
    _ptBinsDY3 = ptBinsDY3_ups;
    _ptBinsDY2 = ptBinsDY2_ups;
    _ptBinsDY1 = ptBinsDY1_ups;
  }
  else if ( _flavTypeToComputeFor == down )
  {
    _ptBinsAll = ptBinsAll_downs;
    _ptBinsDY5 = ptBinsDY5_downs;
    _ptBinsDY4 = ptBinsDY4_downs;
    _ptBinsDY3 = ptBinsDY3_downs;
    _ptBinsDY2 = ptBinsDY2_downs;
    _ptBinsDY1 = ptBinsDY1_downs;
  }
  else if ( _flavTypeToComputeFor == strange )
  {
    _ptBinsAll = ptBinsAll_stranges;
    _ptBinsDY5 = ptBinsDY5_stranges;
    _ptBinsDY4 = ptBinsDY4_stranges;
    _ptBinsDY3 = ptBinsDY3_stranges;
    _ptBinsDY2 = ptBinsDY2_stranges;
    _ptBinsDY1 = ptBinsDY1_stranges;
  }
  else if ( _flavTypeToComputeFor == anti_up )
  {
    _ptBinsAll = ptBinsAll_anti_ups;
    _ptBinsDY5 = ptBinsDY5_anti_ups;
    _ptBinsDY4 = ptBinsDY4_anti_ups;
    _ptBinsDY3 = ptBinsDY3_anti_ups;
    _ptBinsDY2 = ptBinsDY2_anti_ups;
    _ptBinsDY1 = ptBinsDY1_anti_ups;
  }
  else if ( _flavTypeToComputeFor == anti_down )
  {
    _ptBinsAll = ptBinsAll_anti_downs;
    _ptBinsDY5 = ptBinsDY5_anti_downs;
    _ptBinsDY4 = ptBinsDY4_anti_downs;
    _ptBinsDY3 = ptBinsDY3_anti_downs;
    _ptBinsDY2 = ptBinsDY2_anti_downs;
    _ptBinsDY1 = ptBinsDY1_anti_downs;
  }
    else if ( _flavTypeToComputeFor == anti_strange )
  {
    _ptBinsAll = ptBinsAll_anti_stranges;
    _ptBinsDY5 = ptBinsDY5_anti_stranges;
    _ptBinsDY4 = ptBinsDY4_anti_stranges;
    _ptBinsDY3 = ptBinsDY3_anti_stranges;
    _ptBinsDY2 = ptBinsDY2_anti_stranges;
    _ptBinsDY1 = ptBinsDY1_anti_stranges;
  }
  else if ( _flavTypeToComputeFor == allFlavors )
  {
    _ptBinsAll = ptBinsAll_all;
    _ptBinsDY5 = ptBinsDY5_all;
    _ptBinsDY4 = ptBinsDY4_all;
    _ptBinsDY3 = ptBinsDY3_all;
    _ptBinsDY2 = ptBinsDY2_all;
    _ptBinsDY1 = ptBinsDY1_all;
  }
  else
  {
    string errMsg = "error in ptDistribution, flavor not specified";
    throw eAnalysis_error( errMsg );
  }


  FLAVOR_TYPE genFlavor;
  for ( int j = 0; j < n_particles; j++ )
  {
    pt = sqrt( pow( _particles[j].PX, 2 ) + pow( _particles[j].PY, 2 ) );
    y = 0.5 * log(( _particles[j].E + _particles[j].PZ ) / ( _particles[j].E - _particles[j].PZ ) );

    genFlavor = ParticleOffline::mapToGenericFlavorType( _particles[j].FLAVOR );
    if ( _flavTypeToComputeFor == allFlavors || _particles[j].FLAVOR == _flavTypeToComputeFor || 
      ( _flavTypeToComputeFor == light_quark && ( genFlavor == light_quark || genFlavor == anti_light_quark ) ) )
    {
      //------------------------ y in [-inf,inf] -----------------------
      if ( pt <= maxPT && pt >= minPT )
      {
        if ( pt == maxPT )
        {
          ++_ptBinsAll[step][numberBinsPT - 1];
        }
        else
        {
          if ( pt == minPT )
            ++_ptBinsAll[step][0];
          else
            ++_ptBinsAll[step][int(( pt - minPT )/binWidthPT )];
        }

        //----------------------------------------------------------------

        //------------------------ y in [-2.5,2.5] -----------------------
        if ( fabs( y ) <= 2.5 )
        {
          if ( pt == maxPT )
          {
            ++_ptBinsDY5[step][numberBinsPT - 1];
          }
          else
          {
            if ( pt == minPT )
              ++_ptBinsDY5[step][0];
            else
              ++_ptBinsDY5[step][int(( pt - minPT )/binWidthPT )];
          }
        }
        //----------------------------------------------------------------

        //------------------------ y in [-2.0,2.0] -----------------------
        if ( fabs( y ) <= 2.0 )
        {
          if ( pt == maxPT )
          {
            ++_ptBinsDY4[step][numberBinsPT - 1];
          }
          else
          {
            if ( pt == minPT )
              ++_ptBinsDY4[step][0];
            else
              ++_ptBinsDY4[step][int(( pt - minPT )/binWidthPT )];
          }
        }
        //----------------------------------------------------------------

        //------------------------ y in [-1.5,1.5] -----------------------
        if ( fabs( y ) <= 1.5 )
        {
          if ( pt == maxPT )
          {
            ++_ptBinsDY3[step][numberBinsPT - 1];
          }
          else
          {
            if ( pt == minPT )
              ++_ptBinsDY3[step][0];
            else
              ++_ptBinsDY3[step][int(( pt - minPT )/binWidthPT )];
          }
        }
        //----------------------------------------------------------------

        //------------------------ y in [-1.0,1.0] -----------------------
        if ( fabs( y ) <= 1.0 )
        {
          if ( pt == maxPT )
          {
            ++_ptBinsDY2[step][numberBinsPT - 1];
          }
          else
          {
            if ( pt == minPT )
              ++_ptBinsDY2[step][0];
            else
              ++_ptBinsDY2[step][int(( pt - minPT )/binWidthPT )];
          }
        }
        //----------------------------------------------------------------

        //------------------------ y in [-0.5,0.5] -----------------------
        if ( fabs( y ) <= 0.5 )
        {
          if ( pt == maxPT )
          {
            ++_ptBinsDY1[step][numberBinsPT - 1];
          }
          else
          {
            if ( pt == minPT )
              ++_ptBinsDY1[step][0];
            else
              ++_ptBinsDY1[step][int(( pt - minPT )/binWidthPT )];
          }
        }
      }
      //----------------------------------------------------------------
    }

  }
}



void analysis::ptSoftDistribution( const FLAVOR_TYPE _flavTypeToComputeFor, vector<ParticleOffline>& _particles, const int n_particles, const int step )
{
  double pt;

  vector<double>* _ptBinsSoftAll;

  if ( _flavTypeToComputeFor == gluon )
  {
    _ptBinsSoftAll = ptBinsSoftAll_gluons;
  }
  else if ( _flavTypeToComputeFor == light_quark )
  {
    _ptBinsSoftAll = ptBinsSoftAll_quarks;
  }
  else if ( _flavTypeToComputeFor == up )
  {
    _ptBinsSoftAll = ptBinsSoftAll_ups;
  }
  else if ( _flavTypeToComputeFor == down )
  {
    _ptBinsSoftAll = ptBinsSoftAll_downs;
  }
  else if ( _flavTypeToComputeFor == strange )
  {
    _ptBinsSoftAll = ptBinsSoftAll_stranges;
  }
  else if ( _flavTypeToComputeFor == anti_up )
  {
    _ptBinsSoftAll = ptBinsSoftAll_anti_ups;
  }
  else if ( _flavTypeToComputeFor == anti_down )
  {
    _ptBinsSoftAll = ptBinsSoftAll_anti_downs;
  }
  else if ( _flavTypeToComputeFor == anti_strange )
  {
    _ptBinsSoftAll = ptBinsSoftAll_anti_stranges;
  }
  else if ( _flavTypeToComputeFor == allFlavors )
  {
    _ptBinsSoftAll = ptBinsSoftAll_all;
  }
  else
  {
    string errMsg = "error in ptSoftDistribution, flavor not specified";
    throw eAnalysis_error( errMsg );
  }

  FLAVOR_TYPE genFlavor;
  for ( int j = 0; j < n_particles; j++ )
  {
    pt = sqrt( pow( _particles[j].PX, 2 ) + pow( _particles[j].PY, 2 ) );

    genFlavor = ParticleOffline::mapToGenericFlavorType( _particles[j].FLAVOR );
    if ( _flavTypeToComputeFor == allFlavors || _particles[j].FLAVOR == _flavTypeToComputeFor || 
      ( _flavTypeToComputeFor == light_quark && ( genFlavor == light_quark || genFlavor == anti_light_quark ) ) )
    {
      //------------------------ y in [-inf,inf] -----------------------
      if ( pt <= maxPTSoft && pt >= 0 )
      {
        if ( pt == maxPTSoft )
        {
          ++_ptBinsSoftAll[step][numberBinsSoftPT - 1];
        }
        else
        {
          if ( pt == 0 )
            ++_ptBinsSoftAll[step][0];
          else
            ++_ptBinsSoftAll[step][int( pt/binWidthSoftPT )];
        }
      }
      //----------------------------------------------------------------
    }
  }
}



void analysis::yDistribution( const FLAVOR_TYPE _flavTypeToComputeFor, vector<ParticleOffline>& _particles, const int n_particles, const int step )
{
  double y;
  vector<double> * _yBins;


  if ( _flavTypeToComputeFor == gluon )
  {
    _yBins = yBins_gluon;
  }
  else if ( _flavTypeToComputeFor == up )
  {
    _yBins = yBins_up;
  }
  else if ( _flavTypeToComputeFor == down )
  {
    _yBins = yBins_down;
  }
  else if ( _flavTypeToComputeFor == strange )
  {
    _yBins = yBins_strange;
  }
  else if ( _flavTypeToComputeFor == anti_up )
  {
    _yBins = yBins_anti_up;
  }
  else if ( _flavTypeToComputeFor == anti_down )
  {
    _yBins = yBins_anti_down;
  }
  else if ( _flavTypeToComputeFor == anti_strange )
  {
    _yBins = yBins_anti_strange;
  }
  else
  {
    string errMsg = "error in yDistribution, flavor not specified";
    throw eAnalysis_error( errMsg );
  }

  for ( int j = 0; j < n_particles; j++ )
  {
    y = 0.5 * log(( _particles[j].E + _particles[j].PZ ) / ( _particles[j].E - _particles[j].PZ ) );

    if ( _particles[j].FLAVOR == _flavTypeToComputeFor )
    {
      if ( y <= maxY && y >= minY)
      {
        if ( y == maxY )
        {
          ++_yBins[step][numberBinsPT - 1];
        }
        else
        {
          if ( y == minY )
            ++_yBins[step][0];
          else
            ++_yBins[step][int(( y - minY )/binWidthY )];
        }
      }
      //----------------------------------------------------------------
    }
  }
}




void analysis::transverseEnergyDistribution( const FLAVOR_TYPE _flavTypeToComputeFor, vector<ParticleOffline>& _particles, const int n_particles, const int step )
{
  double Et, y;
  vector<double> * _EtBins;
  
  if ( _flavTypeToComputeFor == gluon )
  {
    _EtBins = transverseEnergyGluons;
  }
  else if ( _flavTypeToComputeFor == light_quark )
  {
    _EtBins = transverseEnergyQuarks;
  }
  else if ( _flavTypeToComputeFor == anti_light_quark )
  {
    _EtBins = transverseEnergyAntiQuarks;
  }
  else
  {
    string errMsg = "error in yDistribution, flavor not specified";
    throw eAnalysis_error( errMsg );
  }
  
  for ( int j = 0; j < n_particles; j++ )
  {
    Et = sqrt( pow( _particles[j].PX, 2) + pow( _particles[j].PY, 2) );
    y = 0.5 * log(( _particles[j].E + _particles[j].PZ ) / ( _particles[j].E - _particles[j].PZ ) );
    
    if ( ParticleOffline::mapToGenericFlavorType( _particles[j].FLAVOR ) == _flavTypeToComputeFor )
    {
      if ( y <= maxY && y >= minY)
      {
        if ( y == maxY )
        {
          _EtBins[step][numberBinsPT - 1] += Et;
        }
        else
        {
          if ( y == minY )
            _EtBins[step][0] += Et;
          else
            _EtBins[step][int(( y - minY )/binWidthY )] += Et;
        }
      }
      //----------------------------------------------------------------
    }
  }
}




void analysis::removeJetEvent_in( const int entity_ID )
{
  jetTracker[entity_ID].pop_back();
}


void analysis::makeJetTrackerCopy()
{
  jetTracker_copy = jetTracker;
}


void analysis::restoreJetTracker()
{
  jetTracker = jetTracker_copy;
}


void analysis::exchangeJetID( const int oldID, const int newID )
{
  int entity_index = 0;
  while ( entity_index < jetTracker.size() && jetTracker[entity_index].back().jet_ID_out != -addedParticles[oldID].unique_id )
  {
    entity_index++;
  }

  if ( entity_index < jetTracker.size() )
  {
    jetTracker[entity_index].back().jet_ID_out = -addedParticles[newID].unique_id;
  }
  else
  {
    cout << "error in exchangeJetID()" << endl;
  }

}



void analysis::addJetEvent_initial( const int jetID )
{
  jetTrackerSingleEvent tempEvent;
  tempEvent.jet_ID_in = -1;
  tempEvent.jet_ID_out = -addedParticles[jetID].unique_id;
  tempEvent.coll_type = initial_jet;
  
  tempEvent.flavor_in = -1;
  tempEvent.flavor_out = static_cast<int>( addedParticles[jetID].FLAVOR );

  tempEvent.R_proj[0] = addedParticles[jetID].T;
  tempEvent.R_proj[1] = addedParticles[jetID].X;
  tempEvent.R_proj[2] = addedParticles[jetID].Y;
  tempEvent.R_proj[3] = addedParticles[jetID].Z;

  tempEvent.P_proj_out[0] = addedParticles[jetID].E;
  tempEvent.P_proj_out[1] = addedParticles[jetID].PX;
  tempEvent.P_proj_out[2] = addedParticles[jetID].PY;
  tempEvent.P_proj_out[3] = addedParticles[jetID].PZ;

  vector<jetTrackerSingleEvent> tempVec;
  tempVec.push_back( tempEvent );
  jetTracker.push_back( tempVec );
}


void analysis::addJetEvents_final()
{
  double pt;
  for ( int i = 0; i < addedParticles.size(); i++ )
  {
    pt = sqrt( pow( addedParticles[i].PX, 2.0 ) + pow( addedParticles[i].PY, 2.0 ) );
    if ( pt > jetTracking_PT )
    {
      jetTrackerSingleEvent tempEvent;
      tempEvent.jet_ID_in = -addedParticles[i].unique_id;
      tempEvent.jet_ID_out = -1;
      tempEvent.R_proj[0] = addedParticles[i].T;
      tempEvent.R_proj[1] = addedParticles[i].X;
      tempEvent.R_proj[2] = addedParticles[i].Y;
      tempEvent.R_proj[3] = addedParticles[i].Z;

      tempEvent.P_proj_in[0] = addedParticles[i].E;
      tempEvent.P_proj_in[1] = addedParticles[i].PX;
      tempEvent.P_proj_in[2] = addedParticles[i].PY;
      tempEvent.P_proj_in[3] = addedParticles[i].PZ;
      
      tempEvent.flavor_out = -1;
      tempEvent.flavor_in = static_cast<int>( addedParticles[i].FLAVOR );
      
      tempEvent.coll_type = final_jet;


      int entity_index = 0;
      while ( entity_index < jetTracker.size() && jetTracker[entity_index].back().jet_ID_out != -addedParticles[i].unique_id )
      {
        entity_index++;
      }

      if ( entity_index < jetTracker.size() )
      {
        jetTracker[entity_index].push_back( tempEvent );
      }
    }
  }

}



int analysis::addJetEvent_in( const int ID_1, const int ID_2, const int added_ID, const jetTrackerCollType coll_type,
                              const double cross_section, const int cell_ID, const double lambda )
{
//   double jet_pt = sqrt( pow( addedParticles[added_ID].PX, 2.0 ) + pow( addedParticles[added_ID].PY, 2.0 ) );
//   double pt1 = sqrt( pow( particles_atTimeNow[ID_1].PX, 2.0 ) + pow( particles_atTimeNow[ID_1].PY, 2.0 ) );
//   double pt2 = -1;
//   if ( ID_2 > 0 )
//   {
//     pt2 = sqrt( pow( particles_atTimeNow[ID_2].PX, 2.0 ) + pow( particles_atTimeNow[ID_2].PY, 2.0 ) );
//   }

  int jetID = added_ID;
  int partner1 = ID_1;
  int partner2 = ID_2;


  jetTrackerSingleEvent tempEvent;
  tempEvent.jet_ID_in = -addedParticles[jetID].unique_id;
  tempEvent.R_proj[0] = addedParticles[jetID].T;
  tempEvent.R_proj[1] = addedParticles[jetID].X;
  tempEvent.R_proj[2] = addedParticles[jetID].Y;
  tempEvent.R_proj[3] = addedParticles[jetID].Z;
  tempEvent.flavor_in = static_cast<int>( addedParticles[jetID].FLAVOR );

  tempEvent.P_proj_in[0] = addedParticles[jetID].E;
  tempEvent.P_proj_in[1] = addedParticles[jetID].PX;
  tempEvent.P_proj_in[2] = addedParticles[jetID].PY;
  tempEvent.P_proj_in[3] = addedParticles[jetID].PZ;

  tempEvent.P1_in[0] = particles_atTimeNow[partner1].E;
  tempEvent.P1_in[1] = particles_atTimeNow[partner1].PX;
  tempEvent.P1_in[2] = particles_atTimeNow[partner1].PY;
  tempEvent.P1_in[3] = particles_atTimeNow[partner1].PZ;

  if ( coll_type == c3to2 )
  {
    tempEvent.P2_in[0] = particles_atTimeNow[partner2].E;
    tempEvent.P2_in[1] = particles_atTimeNow[partner2].PX;
    tempEvent.P2_in[2] = particles_atTimeNow[partner2].PY;
    tempEvent.P2_in[3] = particles_atTimeNow[partner2].PZ;
  }

  tempEvent.lambda = lambda;
  tempEvent.xSection = cross_section;
  tempEvent.cell_ID = cell_ID;
  tempEvent.coll_type = coll_type;
  


  int entity_index = 0;
  while ( entity_index < jetTracker.size() && jetTracker[entity_index].back().jet_ID_out != -addedParticles[jetID].unique_id )
  {
    entity_index++;
  }

  if ( entity_index < jetTracker.size() )
  {
    jetTracker[entity_index].push_back( tempEvent );
  }
  else
  {
    vector< jetTrackerSingleEvent > tempEntity;
    tempEntity.push_back( tempEvent );
    jetTracker.push_back( tempEntity );
  }

  return entity_index;
}


void analysis::addJetEvent_out( const int entity_ID, const int added_ID, const int ID_2, const int ID_3, const jetTrackerCollType coll_type )
{
  double pt1 = sqrt( pow( addedParticles[added_ID].PX, 2.0 ) + pow( addedParticles[added_ID].PY, 2.0 ) );
  double pt2 = sqrt( pow( particles_atTimeNow[ID_2].PX, 2.0 ) + pow( particles_atTimeNow[ID_2].PY, 2.0 ) );

  double pt3 = -1;
  if ( coll_type == c2to3 )
  {
    pt3 = sqrt( pow( addedParticles[ID_3].PX, 2.0 ) + pow( addedParticles[ID_3].PY, 2.0 ) );
  }

  int jetID;
  int partner1 = -1, partner2 = -1;

  jetID = added_ID;
  partner1 = ID_2;
  partner2 = ID_3;
  
  if ( pt3 > pt1 )
  {
    jetID = ID_3;
    partner2 = added_ID;
  }

  if ( entity_ID != -1 )
  {
    jetTracker[entity_ID].back().jet_ID_out = -addedParticles[jetID].unique_id;

    jetTracker[entity_ID].back().P_proj_out[0] = addedParticles[jetID].E;
    jetTracker[entity_ID].back().P_proj_out[1] = addedParticles[jetID].PX;
    jetTracker[entity_ID].back().P_proj_out[2] = addedParticles[jetID].PY;
    jetTracker[entity_ID].back().P_proj_out[3] = addedParticles[jetID].PZ;
    
    jetTracker[entity_ID].back().flavor_out = static_cast<int>( addedParticles[jetID].FLAVOR );
    

    jetTracker[entity_ID].back().P1_out[0] = particles_atTimeNow[partner1].E;
    jetTracker[entity_ID].back().P1_out[1] = particles_atTimeNow[partner1].PX;
    jetTracker[entity_ID].back().P1_out[2] = particles_atTimeNow[partner1].PY;
    jetTracker[entity_ID].back().P1_out[3] = particles_atTimeNow[partner1].PZ;

    if ( coll_type == c2to3 )
    {
      jetTracker[entity_ID].back().P2_out[0] = addedParticles[partner2].E;
      jetTracker[entity_ID].back().P2_out[1] = addedParticles[partner2].PX;
      jetTracker[entity_ID].back().P2_out[2] = addedParticles[partner2].PY;
      jetTracker[entity_ID].back().P2_out[3] = addedParticles[partner2].PZ;
    }
  }
  else if ( sqrt( pow( addedParticles[jetID].PX, 2.0 ) + pow( addedParticles[jetID].PY, 2.0 ) ) > jetTracking_PT )
  {
    jetTrackerSingleEvent tempEvent;
    tempEvent.jet_ID_in = -1;
    tempEvent.jet_ID_out = -addedParticles[jetID].unique_id;
    tempEvent.coll_type = production;
    tempEvent.R_proj[0] = addedParticles[jetID].T;
    tempEvent.R_proj[1] = addedParticles[jetID].X;
    tempEvent.R_proj[2] = addedParticles[jetID].Y;
    tempEvent.R_proj[3] = addedParticles[jetID].Z;
    tempEvent.flavor_in = -1;
    tempEvent.flavor_out = static_cast<int>( addedParticles[jetID].FLAVOR );

    tempEvent.P_proj_out[0] = addedParticles[jetID].E;
    tempEvent.P_proj_out[1] = addedParticles[jetID].PX;
    tempEvent.P_proj_out[2] = addedParticles[jetID].PY;
    tempEvent.P_proj_out[3] = addedParticles[jetID].PZ;

    tempEvent.P1_out[0] = particles_atTimeNow[partner1].E;
    tempEvent.P1_out[1] = particles_atTimeNow[partner1].PX;
    tempEvent.P1_out[2] = particles_atTimeNow[partner1].PY;
    tempEvent.P1_out[3] = particles_atTimeNow[partner1].PZ;

    if ( coll_type == c2to3 )
    {
      tempEvent.P2_out[0] = addedParticles[partner2].E;
      tempEvent.P2_out[1] = addedParticles[partner2].PX;
      tempEvent.P2_out[2] = addedParticles[partner2].PY;
      tempEvent.P2_out[3] = addedParticles[partner2].PZ;
    }

    vector< jetTrackerSingleEvent > tempEntity;
    tempEntity.push_back( tempEvent );
    jetTracker.push_back( tempEntity );
  }
}



void analysis::particleOutput( const int step )
{
  string name;
  stringstream ss;

  if ( step == 0 )
    name = "initial";
  else if ( step == nTimeSteps )
    name = "final";
  else
  {
    ss << step;
    name = "step" + ss.str();
  }

  //creates filename, for example: "./output/run32_step1.f1", "./output/run32_initial.f1" or the likes
  string filename;
  char * temp = getenv( "PBS_JOBID" );
  if ( temp != NULL )
  {
    string jobID( temp );
    filename = "/local/" + jobID + "/" + theConfig->getJobName() + "_" + name + ".f1";
  }
  else
  {
    filename = "./output/" + theConfig->getJobName() + "_" + name + ".f1";
  }
  fstream file( filename.c_str(), ios::out | ios::trunc );

  //---- print header if file is empty ----
  time_t end;
  time( &end );

  file.seekp( 0, ios::end );
  long size = file.tellp();
  file.seekp( 0, ios::beg );
  if ( size == 0 )
    printHeader( file, all, end );
  //---------------------------------------

  for ( int i = 0; i < addedParticles.size(); i++ )
  {
    file << i << sep << addedParticles[i].unique_id << sep << addedParticles[i].cell_id << sep << addedParticles[i].FLAVOR << sep << addedParticles[i].T << sep << addedParticles[i].X << sep
    << addedParticles[i].Y << sep  << addedParticles[i].Z << sep << addedParticles[i].E << sep << addedParticles[i].PX << sep << addedParticles[i].PY << sep
    << addedParticles[i].PZ << sep << addedParticles[i].md2g << sep << addedParticles[i].md2q << endl;
  }
  file.close();
}



void analysis::writePartclMovie( vector< ParticleOffline >& _particles, const int n_particles, fstream& _oscar, const int step, const int jumpSteps )
{
  const string sep = "  ";
  const int width = 14;
  double cc, dt;
  double t, x, y, z;
  const double zero = 0.0;

  double time;

  const int nOutput = 1; // every nOutput-th particle is writen to particle file, eg nOutput=2 every second particle, nOutput=1 for every particle
  int nCount_selected = nOutput; // dummy variable for counting

  if ( step == 0 )
  {
    time = 0.0;
  }
  else if ( step == nTimeSteps_movie - 1 )
  {
    return;
  }
  else
  {
    time = tstep_movie[step];
  }

  // determine number of timesteps, nTimeSteps is too large if runtime is shorter than last timestep
  int numberOfTimeSteps = 0;
  while ( theConfig->getRuntime() >= tstep_movie[numberOfTimeSteps] && numberOfTimeSteps < nTimeSteps_movie - 1 )
  {
    numberOfTimeSteps++;
  }
  numberOfTimeSteps -= jumpSteps;

  int numberActiveParticles = 0;
  for ( int i = 0; i < n_particles; i++ )
  {
    if ( _particles[i].T_creation <= time )
    {
      if ( nCount_selected != nOutput ) // do not write particle to file, increase count by one
      {
        nCount_selected++;
      }
      else // write particle to file
      {
        numberActiveParticles++;
        nCount_selected = 1;
      }
    }
  }
  nCount_selected = nOutput; // reset dummy variable for counting



  // <<---------------------------------------------------
  // oscar output

  if ( step == 0 ) // header, no particles yet
  {
    // file header
    _oscar << "OSC1997A    " << endl;
    _oscar << "final_id_p_x" << endl;

    // code name etc.
    // code_name version (aproj, zproj)+(atarg, ztarg) refframe, ebeam, ntestpart
    _oscar << ::std::scientific << ::std::uppercase
    << "BAMPS        0.2    (" << ::std::setw( 3 ) << int( theConfig->getA() ) << "," << ::std::setw( 6 )
    <<  int( theConfig->getAatomic() ) << ")+(" << ::std::setw( 3 ) << int( theConfig->getB() ) << ","
    << ::std::setw( 6 ) << int( theConfig->getBatomic() ) << ")  CMS " << ::std::setprecision( 4 )
    << theConfig->getSqrtS() << "  " << ::std::setw( 8 ) << theConfig->getTestparticles() << endl;
  }
  else // all other timesteps
  {
    // event header for this timestep
    // event npart bimp phi #timesteps timestep
    _oscar << "          1  " << ::std::setw( 10 ) << numberActiveParticles << "  " << ::std::resetiosflags( ::std::ios::scientific )
    << ::std::setw( 8 ) << ::std::fixed << ::std::setprecision( 3 ) << theConfig->getImpactParameter() << "  "
    << ::std::setw( 8 ) << ::std::fixed << ::std::setprecision( 3 ) << 0.000 << "  "  << ::std::setw( 4 )
    << numberOfTimeSteps << "  " << ::std::setw( 4 ) << step << endl;


    // write particle data
    // ipart, id, px, py, pz, p0, mass, x, y, z, t, flag#1, x_f, y_f, z_f, t_f, iflag#1
    for ( int i = 0; i < n_particles; i++ )
    {
      if ( _particles[i].T_creation <= time )
      {
        // find out if particle should be written to file (in case nOutput is not 1)
        if ( nCount_selected != nOutput ) // do not write particle to file, increase count by one
        {
          nCount_selected++;
        }
        else // write particle to file
        {
          _oscar << ::std::setw( 11 ) << i+1;
          _oscar << ::std::setw( 12 ) << ParticleOffline::mapToPDGCodes( _particles[i].FLAVOR ); //PDG group codes
          nCount_selected = 1;

          _oscar << ::std::scientific << ::std::uppercase << ::std::setprecision( 6 )
          << setw( width ) << _particles[i].PX << setw( width ) << _particles[i].PY << setw( width ) << _particles[i].PZ << setw( width )
          << _particles[i].E << setw( width ) << _particles[i].m;

          if ( _particles[i].T == time )
          {
            t = _particles[i].T;
            x = _particles[i].X;
            y = _particles[i].Y;
            z = _particles[i].Z;
          }
          else if ( _particles[i].T > time )
          {
            dt = time - _particles[i].T;
            cc = dt / _particles[i].E;
            t = _particles[i].T + dt;
            x = _particles[i].X + _particles[i].PX * cc;
            y = _particles[i].Y + _particles[i].PY * cc;
            z = _particles[i].Z + _particles[i].PZ * cc;
          }
          else
          {
            cout << "error in write movie particle data, particles from the past" << endl;
            cout << "timestep=" << time << "   time _particles=" << _particles[i].T << endl;
            dt = time - _particles[i].T;
            cc = dt / _particles[i].E;
            t = _particles[i].T + dt;
            x = _particles[i].X + _particles[i].PX * cc;
            y = _particles[i].Y + _particles[i].PY * cc;
            z = _particles[i].Z + _particles[i].PZ * cc;
          }

          _oscar << setw( width ) << x << setw( width ) << y << setw( width ) << z << setw( width ) << t;

          //       _oscar << "0.0" << sep << "0.0" << sep << "0.0" << sep << "0.0" << sep << "0.0" << sep << "0" << endl;
          //         _oscar << zero << sep << zero << sep << zero << sep << zero << sep << zero << ::std::setw( 10 ) << int(zero) << endl;
          //         _oscar << zero << sep << _particles[i].X_lastInt << sep << _particles[i].Y_lastInt << sep << _particles[i].Z_lastInt << sep << _particles[i].T_lastInt << ::std::setw( 10 ) << int(zero) << endl;
          _oscar << setw( width ) << zero << setw( width ) << _particles[i].X_lastInt << setw( width ) << _particles[i].Y_lastInt << setw( width ) << _particles[i].Z_lastInt << setw( width ) << _particles[i].T_lastInt << ::std::setw( 10 ) << -(_particles[i].unique_id) << endl;
          // --------------------------------------------------->>
        }
      }
    }
  }
}






void analysis::printHeader( fstream & f, const anaType mode, const time_t end )
{
  switch ( mode )
  {
  case ptSpectrum:
    f << "#pt-spectra " << endl;
    break;
  case ptSpectrumSoft:
    f << "#soft pt-spectra " << endl;
    break;
  case jets:
    f << "#jet tracking informatino " << endl;
    break;
  case all:
    f << "#information on particle locations, momenta etc. at step indicated by filename" << endl;
    f << "#initial = 0 fm/c" << endl;
    for ( int i = 0;i < nTimeSteps - 1;i++ )
      f << "#step " << ( i + 1 ) << " = " << tstep[i] << " fm/c" << endl;
    f << "#final = " << theConfig->getRuntime() << " fm/c" << endl;
    break;
  case rapidityDistribution:
    f << "#rapidity distribution " << endl;
    f << "#one block per timestep " << endl;
    break;
  case quarkNumbers:  
    f << "#quark numbers summed over all rapidities" << endl;
    break;
  default:
    f << "#undefined ";
  }

  f << "#start: " << ctime( &start );
  f << "#end: " << ctime( &end );
  f << "#" << endl;
  f << "#simulation parameter:" << endl;
  f << "#testparticles= " << theConfig->getTestparticles() << endl;
  f << "#runtime= " << theConfig->getRuntime() << endl;
  f << "#sqrtS= " << theConfig->getSqrtS() << " GeV" << endl;
  f << "#P0= " << theConfig->getPtCutoff() << " GeV" << endl;
  f << "#b= " << theConfig->getImpactParameter() << " fm" << endl;
  f << "#(" << theConfig->getA() << "," << theConfig->getAatomic() << ") on ("
  << theConfig->getB() << "," << theConfig->getBatomic() << ")" << endl;
  f << "#seed for random generator ran2(): " << seed << endl;
  f << "#" << endl;

  stringstream ss;

  switch ( mode )
  {
  case ptSpectrum:
    f << "#numbers NOT yet corrected for width of bins and number of testparticles!" << endl;
    f << "#binWidth= " << binWidthPT << endl;
    f << "#testparticles= " << theConfig->getTestparticles() << endl;
    f << "#" << endl;
    f << "#pt [GeV]" << sep;
    for ( int j = 0; j <= nTimeSteps; j++ )
    {
      if ( j == 0 )
        f << "initial" << sep;
      else if ( j == nTimeSteps )
        f << "final = " << theConfig->getRuntime() << " fm/c";
      else if ( tstep[j-1] <= theConfig->getRuntime() )
      {
        ss << tstep[j-1];
        f << ss.str() << " fm/c" << sep;
        ss.str( "" );
      }
    }
    f << endl;
    break;
  case ptSpectrumSoft:
    f << "#numbers NOT yet corrected for width of bins and number of testparticles!" << endl;
    f << "#binWidth= " << binWidthSoftPT << endl;
    f << "#testparticles= " << theConfig->getTestparticles() << endl;
    f << "#" << endl;
    f << "#pt [GeV]" << sep;
    for ( int j = 0; j <= nTimeSteps; j++ )
    {
      if ( j == 0 )
        f << "initial" << sep;
      else if ( j == nTimeSteps )
        f << "final = " << theConfig->getRuntime() << " fm/c";
      else if ( tstep[j-1] <= theConfig->getRuntime() )
      {
        ss << tstep[j-1];
        f << ss.str() << " fm/c" << sep;
        ss.str( "" );
      }
    }
    f << endl;
    break;
  case jets:
    f << "#jetID_in" << sep << "jetID_out" << sep << "coll_type" << sep
    << "P_in[0]" << sep << "P_in[1]" << sep << "P_in[2]" << sep  << "P_in[3]" << sep
    << "P_out[0]" << sep << "P_out[1]" << sep << "P_out[2]" << sep  << "P_out[3]" << sep
    << "R[0]" << sep << "R[1]" << sep << "R[2]" << sep << "R[3]" << sep
    << "xSection" << sep << "lambda" << endl;
    break;
  case rapidityDistribution:
    f << "# y   gluon   up   down   strange   anti-up   anti-down   anti-strange" << endl;
    break;
  case quarkNumbers:
    f << "# t   gluon   up   down   strange   anti-up   anti-down   anti-strange" << endl;
    break;
  case all:
    f << "#ID" << sep << "uniqueID" << sep << "Cell ID" << sep << "Flavor" << sep << "t" << sep << "x" << sep << "y" << sep << "z" << sep
    << "E"  << sep << "Px"  << sep << "Py"  << sep << "Pz" << sep << "md2g/alpha_s [GeV^2]" << sep
    << "md2q/alpha_s [GeV^2]" << endl;
    break;
  default:
    f << endl;
  }
  f << "#" << endl;
}



void analysis::mfpJetsOutput( const int step, const int jumpSteps )
{
  const string sep = "\t ";
  double time;

  if ( step == 0 )
  {
    time = 0.0;
  }
  else if ( step == nTimeSteps_movie - 1 )
  {
    return;
  }
  else
  {
    time = tstep_movie[step];
  }
  
  mfpJetsOutputFile << time;
  for ( int i = 0; i < rings.size(); i++ )
  {
    if ( rings[i].collectedGluon != 0 )
    {
      mfpJetsOutputFile << sep << rings[i].lambdaGluon / rings[i].collectedGluon;
    }
    else
    {
      mfpJetsOutputFile << sep << 0; 
    }
    if ( rings[i].collectedQuark != 0 )
    {
      mfpJetsOutputFile << sep << rings[i].lambdaQuark / rings[i].collectedQuark;
    }
    else
    {
      mfpJetsOutputFile << sep << 0; 
    }
    mfpJetsOutputFile << sep << rings[i].collectedGluon << sep << rings[i].collectedQuark;  
  }
  mfpJetsOutputFile << endl;
  
  rings.clear();

}



void analysis::printCentralDensities(const double _time)
{
  const string sep = "\t ";
  
  centralDensitiesOutputFile << _time << sep << centralRingsCopyFromCascade.y_left << sep << centralRingsCopyFromCascade.y_right; 
  for ( int i = 0; i < centralRingsCopyFromCascade.size(); i++ )
  {
    centralDensitiesOutputFile << sep << centralRingsCopyFromCascade[i].getEnergyDensity() << sep 
    << centralRingsCopyFromCascade[i].getGluonDensity() << sep << centralRingsCopyFromCascade[i].getQuarkDensity(); 
  }
  centralDensitiesOutputFile << endl;

  centralRingsCopyFromCascade.clear();
}





jetTrackerSingleEvent::jetTrackerSingleEvent()
{
  for ( int i = 0; i < 4; i++ )
  {
    R_proj[i] = 0;
    P_proj_in[i] = P_proj_out[i] = 0;
    P1_in[i] = P2_in[i] = 0;
    P1_out[i] = P2_out[i] = 0;
  }
  xSection = -1;
  lambda = -1;
  cell_ID = -1;
}


jetTrackerSingleEvent::~jetTrackerSingleEvent()
{

}



void analysis::jetTrackerOutput()
{
  //---- get time and calculate real time needed by simulation ----
  time_t end;
  time( &end );
  int secs = ( int )difftime( end, start );
  int hours = secs / 3600, minutes = ( secs % 3600 ) / 60, seconds = secs - 3600 * hours - 60 * minutes;
  cout << hours << "h" << minutes << "m" << seconds << "s" <<  endl;
  //---------------------------------------------------------------
  
  
  //---- write output file for pt-spectra -------------------------
  char * temp = getenv( "PBS_JOBID" );
  string filename;
  if ( temp != NULL )
  {
    string jobID( temp );
    filename = "/local/" + jobID + "/" + theConfig->getJobName() + ".f4";
  }
  else
  {
    filename = "./output/" + theConfig->getJobName() + ".f4";
  }
  fstream file( filename.c_str(), ios::out | ios::trunc );
  
  //---- print header if file is empty ----
  file.seekp( 0, ios::end );
  long size = file.tellp();
  file.seekp( 0, ios::beg );
  if ( size == 0 )
    printHeader( file, jets, end );
  //---------------------------------------
  
  for ( int jet = 0; jet < jetTracker.size(); jet++ )
  {
    for ( int event = 0; event < jetTracker[jet].size(); event++ )
    {
      file << jetTracker[jet][event].jet_ID_in << sep << jetTracker[jet][event].jet_ID_out << sep
      << jetTracker[jet][event].coll_type << sep;
      for ( int i = 0; i < 4; i++ )
      {
        file << jetTracker[jet][event].P_proj_in[i] << sep;
      }
      for ( int i = 0; i < 4; i++ )
      {
        file << jetTracker[jet][event].P_proj_out[i] << sep;
      }
      for ( int i = 0; i < 4; i++ )
      {
        file << jetTracker[jet][event].R_proj[i] << sep;
      }
      file << jetTracker[jet][event].xSection << sep << jetTracker[jet][event].lambda 
      << sep << jetTracker[jet][event].flavor_in << sep << jetTracker[jet][event].flavor_out << endl;
    }
    file << endl << endl;
  }
  
}


analysisRingStructure::analysisRingStructure( const int _nRings, const double _centralRadius, const double _deltaR ) : numberOfRings( _nRings ),
    centralRingRadius( _centralRadius ), deltaR( _deltaR )
{
  rings.resize( _nRings );

  rings[0].relocate( 0, _centralRadius );

  totalRadius = _centralRadius + ( _nRings - 1 ) * deltaR;
  for ( int i = 1; i < rings.size(); i++ )
  {
    rings[i].relocate( rings[i-1].maxRadius, rings[i-1].maxRadius + deltaR );
  }
}



void analysisRingStructure::resize( const int _nRings, const double _centralRadius, const double _deltaR )
{
  numberOfRings = _nRings;
  centralRingRadius = _centralRadius;
  deltaR = _deltaR;
  
  rings.clear();
  rings.resize( _nRings );

  rings[0].relocate( 0, _centralRadius );

  totalRadius = _centralRadius + ( _nRings - 1 ) * deltaR;
  for ( int i = 1; i < rings.size(); i++ )
  {
    rings[i].relocate( rings[i-1].maxRadius, rings[i-1].maxRadius + deltaR );
  }
}



int analysisRingStructure::getIndex( const double _xt ) const
{
  int index = getIndexPure( _xt );
  if ( index >= rings.size() )
  {
    return ( static_cast<int>( rings.size() ) - 1 );
  }
  else
  {
    return index;
  }
}



int analysisRingStructure::getIndexPure( const double _xt ) const
{
  if ( _xt < 0 )
  {
    std::string errMsg = "transverse position xt < 0";
    throw eAnalysis_error( errMsg );
  }
  
  if ( _xt < centralRingRadius )
  {
    return 0;
  }
  else
  {
//     if ( _xt > totalRadius )
//     {
//       std::string errMsg = "transverse position xt > R (R = total radius of ring structure)";
//       throw eAnalysis_error( errMsg );
//     }
    
    int index = static_cast<int>( ( _xt - centralRingRadius ) / deltaR ) + 1;  // +1 since index 0 is for the central ring
    return index;
  }
}





int analysisRingStructure::getIndex( const ParticleOffline& _particle ) const
{
  double xt = sqrt( pow( _particle.X, 2 ) + pow( _particle.Y, 2 ) );
  
  return getIndex( xt );
}



analysisRingContainer& analysisRingStructure::operator[]( const int index )
{
  if ( index < 0 || index >= numberOfRings )
  {
    std::string errMsg = "index out of range in analysisRingStructure";
    throw eAnalysis_error( errMsg );
  }
  
  return rings[ index ];
}



analysisRingStructure& analysisRingStructure::operator+=( analysisRingStructure& rhs )
{
  if ( !( rhs.size() == ( *this ).size() ) )
  {
    std::string errMsg = "analysisRingStructure::operator+= only works for equally sized structures";
    throw eAnalysis_error( errMsg );
  }

  for ( int i = 0; i < numberOfRings; i++ )
  {
    rings[i] += rhs[i];
  }

  return ( *this );
}



analysisRingContainer& analysisRingContainer::operator+=( const analysisRingContainer & rhs )
{
  lambdaGluon += rhs.lambdaGluon;
  lambdaQuark += rhs.lambdaQuark;
  collectedGluon += rhs.collectedGluon;
  collectedQuark += rhs.collectedQuark;
  
  return ( *this );
}



// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;
