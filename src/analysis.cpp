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
#include "FPT_compare.h"

#include <stdio.h> // for getenv()
#include <stdlib.h> // for getenv()
#include <time.h>


using namespace std;
using namespace ns_casc;

extern int IX, IY;

namespace ns_heavy_quarks
{
  extern int jpsi_dissociation;
  extern int jpsi_dissociation_from_temperature;
  extern int jpsicreation;
  extern int charmAnnihil;
}


analysis::analysis( config* const c ):
    theConfig( c ),
    rings( c->getRingNumber(), c->getCentralRingRadius(), c->getDeltaR() ),
    centralRingsCopyFromCascade( c->getRingNumber(), c->getCentralRingRadius(), c->getDeltaR() ),
    filename_prefix( c->getStandardOutputDirectoryName() + "/" + c->getJobName() ),
    v2output( c->isV2RAAoutput() ),
    dndyOutput( c->isDndyOutput() ),
    particleCorrelationsOutput( false ),
    hadronization_hq( c->isHadronizationHQ() ),
    mesonDecay( c->isMesonDecay() ),
//     charmTestJet( theConfig->isCharmTestJet() ),
    v2outputIntermediateSteps( c->isV2RAAoutputIntermediateSteps() ),
    studyParticleOutput( c->doOutput_detailedParticleOutput() ),
    outputScheme( c->getOutputScheme() )
{
  //---get time and date---
  time( &start );
  //-----------------------
  
  handle_output_studies( outputScheme );
  
  // to create the tsteps use such a bash script:
  // for (( I=50; $I <= 80; I++ )); do   J=$(echo "scale=1; ($I+1)/10" | bc); echo "tstep[$I]=$J;" ; done
  if( studyJpsi ) // jpsi evolution: more timesteps
  {
    tstep[0]=.005;
    tstep[1]=.006; // choose to be small that no analysis is performed, just start for all files at 0.2 in next step
    tstep[2]=.20;
    tstep[3]=.25;
    tstep[4]=.30;
    tstep[5]=.35;
    tstep[6]=.40;
    tstep[7]=.45;
    tstep[8]=.50;
    tstep[9]=.55;
    tstep[10]=.60;
    tstep[11]=.65;
    tstep[12]=.70;
    tstep[13]=.75;
    tstep[14]=.80;
    tstep[15]=.85;
    tstep[16]=.90;
    tstep[17]=.95;
    tstep[18]=1.00;
    tstep[19]=1.1;
    tstep[20]=1.2;
    tstep[21]=1.3;
    tstep[22]=1.4;
    tstep[23]=1.5;
    tstep[24]=1.6;
    tstep[25]=1.7;
    tstep[26]=1.8;
    tstep[27]=1.9;
    tstep[28]=2.0;
    tstep[29]=2.1;
    tstep[30]=2.2;
    tstep[31]=2.3;
    tstep[32]=2.4;
    tstep[33]=2.5;
    tstep[34]=2.6;
    tstep[35]=2.7;
    tstep[36]=2.8;
    tstep[37]=2.9;
    tstep[38]=3.0;
    tstep[39]=3.3;
    tstep[40]=3.7;
    tstep[41]=4.0;
    tstep[42]=4.3;
    tstep[43]=4.7;
    tstep[44]=5.0;
    tstep[45]=5.3;
    tstep[46]=5.7;
    tstep[47]=6.0;
    tstep[48]=6.3;
    tstep[49]=6.7;
    tstep[50]=7.0;
    tstep[51]=7.3;
    tstep[52]=7.7;
    tstep[53]=8.0;
    tstep[54]=8.5;
    tstep[55]=9.0;
    tstep[56]=9.5;
    tstep[57]=10.0;
    tstep[58]=10.5;                                                                                                                                                                     
    tstep[59]=11.0;                                                                                                                                                                     
    tstep[60]=11.5;                                                                                                                                                                     
    tstep[61]=12.0;                                                                                                                                                                     
    tstep[62]=12.5;                                                                                                                                                                     
    tstep[63]=13.0;                                                                                                                                                                     
    tstep[64]=13.5;                                                                                                                                                                     
    tstep[65]=14.0;                                                                                                                                                                     
    tstep[66]=14.5;                                                                                                                                                                     
    tstep[67]=15.0;                                                                                                                                                                     
    tstep[68]=15.5;                                                                                                                                                                     
    tstep[69]=16.0;                                                                                                                                                                     
    tstep[70]=16.5;                                                                                                                                                                     
    tstep[71]=17.0;                                                                                                                                                                     
    tstep[72]=17.5;                                                                                                                                                                     
    tstep[73]=18.0;                                                                                                                                                                     
    tstep[74]=18.5;                                                                                                                                                                     
    tstep[75]=19.0;                                                                                                                                                                     
    tstep[76]=19.5;                                                                                                                                                                     
    tstep[77]=20.0;
    tstep[78] = infinity; //fm/c
    nTimeSteps = 79;
    
    if( theConfig->getRuntime() > 20.0 )
    {
      tstep[78]=21;
      tstep[79]=22;
      tstep[80]=23;
      tstep[81]=24;
      tstep[82]=25;
      tstep[83]=26;
      tstep[84]=27;
      tstep[85]=28;
      tstep[86]=29;
      tstep[87]=30;
      tstep[88]=31;
      tstep[89]=32;
      tstep[90]=33;
      tstep[91]=34;
      tstep[92]=35;
      tstep[93]=36;
      tstep[94]=37;
      tstep[95]=38;
      tstep[96]=39;
      tstep[97]=40;
      tstep[98]=42;
      tstep[99]=44;
      tstep[100]=46;
      tstep[101]=48;
      tstep[102]=50;
      tstep[103]=52;
      tstep[104]=54;
      tstep[105]=56;
      tstep[106]=58;
      tstep[107]=60;
      tstep[108]=62;
      tstep[109]=64;
      tstep[110]=66;
      tstep[111]=68;
      tstep[112]=70;
      tstep[113]=72;
      tstep[114]=74;
      tstep[115]=76;
      tstep[116]=78;
      tstep[117]=80;
      tstep[118] = infinity; //fm/c
      nTimeSteps = 119;
    }
  }
  else
  {
    tstep[0] = 0.1;      //fm/c
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
  }
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
 
  
  //--------------------------------------
  string name_fug, name_temp, mfpName, centralDensitiesName, oscarName_jet, oscarName_background;
  if ( studyJpsi )
    name_fug = filename_prefix + "_jpsi_fugacity";
  else
    name_fug = "/dev/null";
  
  if( studyTempInTube )
    name_temp = filename_prefix + "_temperature";
  else
    name_temp = "/dev/null";
  
  if( studyJets )
    mfpName = filename_prefix + "_mfp_jets";
  else
    mfpName = "/dev/null";
  
  if( studyCentralDensity )
    centralDensitiesName = filename_prefix + "_central_density";
  else
    centralDensitiesName = "/dev/null";
  
  // movie output
  if( theConfig->doOutput_movieOutputBackground() )
    oscarName_background = filename_prefix + "_background.oscar";
  else
    oscarName_background = "/dev/null";
  
  if( theConfig->doOutput_movieOutputJets() )
    oscarName_jet = filename_prefix + "_jets.oscar";
  else
    oscarName_jet = "/dev/null";

  printJpsiFugacity.open( name_fug.c_str(), ios::out );
  printTempInTube.open( name_temp.c_str(), ios::out );
  mfpJetsOutputFile.open( mfpName.c_str(), ios::out );
  centralDensitiesOutputFile.open( centralDensitiesName.c_str(), ios::out );
  oscarBackground.open( oscarName_background.c_str(), ios::out );
  oscarJets.open( oscarName_jet.c_str(), ios::out );
  //--------------------------------------

  jetTracking_PT = 10.0;
  
  

  if( studyPtSpectra )
  {
    //---- initialisation of PT-binning ----
    minPT = 1.4;
    maxPTSoft = 3.0;
    maxPT = 34.4;
    binWidthPT = 1.0;
    numberBinsPT = int(( maxPT - minPT + 0.001 ) / binWidthPT );
    
    tArrayOfDoubleVec tmpArray;
    for ( unsigned int i = 0; i < rapidityRanges.size(); i++ )
    {
      tmpArray.reset( new vector<double>[nTimeSteps+2] );     //+2 because of initial and final timesteps
      ptBins_gluons.push_back( tmpArray );
      tmpArray.reset( new vector<double>[nTimeSteps+2] );     //+2 because of initial and final timesteps
      ptBins_quarks.push_back( tmpArray );
      tmpArray.reset( new vector<double>[nTimeSteps+2] );     //+2 because of initial and final timesteps
      ptBins_ups.push_back( tmpArray );
      tmpArray.reset( new vector<double>[nTimeSteps+2] );     //+2 because of initial and final timesteps
      ptBins_downs.push_back( tmpArray );
      tmpArray.reset( new vector<double>[nTimeSteps+2] );     //+2 because of initial and final timesteps
      ptBins_stranges.push_back( tmpArray );
      tmpArray.reset( new vector<double>[nTimeSteps+2] );     //+2 because of initial and final timesteps
      ptBins_anti_ups.push_back( tmpArray );
      tmpArray.reset( new vector<double>[nTimeSteps+2] );     //+2 because of initial and final timesteps
      ptBins_anti_downs.push_back( tmpArray );
      tmpArray.reset( new vector<double>[nTimeSteps+2] );     //+2 because of initial and final timesteps
      ptBins_anti_stranges.push_back( tmpArray );
      tmpArray.reset( new vector<double>[nTimeSteps+2] );     //+2 because of initial and final timesteps
      ptBins_all.push_back( tmpArray );
    }

    for ( unsigned int ny = 0; ny < rapidityRanges.size(); ny++ )
    {
      for ( int nt = 0; nt < nTimeSteps + 2; nt++ )
      {
        ptBins_gluons[ny][nt].resize( numberBinsPT + 1, 0 );
        ptBins_quarks[ny][nt].resize( numberBinsPT + 1, 0 );
        ptBins_ups[ny][nt].resize( numberBinsPT + 1, 0 );
        ptBins_downs[ny][nt].resize( numberBinsPT + 1, 0 );
        ptBins_stranges[ny][nt].resize( numberBinsPT + 1, 0 );
        ptBins_anti_ups[ny][nt].resize( numberBinsPT + 1, 0 );
        ptBins_anti_downs[ny][nt].resize( numberBinsPT + 1, 0 );
        ptBins_anti_stranges[ny][nt].resize( numberBinsPT + 1, 0 );
        ptBins_all[ny][nt].resize( numberBinsPT + 1, 0 );
      }
    }
    
    //write the bin labels
    for ( int i = 0; i < numberBinsPT; i++ )
    {
      ptBinLabels.push_back( minPT + ( i * binWidthPT ) + ( binWidthPT / 2 ) );
    }
    cout << "number of bins: " << numberBinsPT << "  binWidth: " << binWidthPT << endl;
    //---- initialisation of PT-binning ----
    
    
    //---- initialisation of softPT-binning ----
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
  
  
  
  
  
  
  if( studyYDistribution || studyEtSpectra )
  {
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
  }

  if( studyEtSpectra )
  {
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
  }


  if( studyJpsi )
  {
    numberJpsi_all_time = new int[nTimeSteps+1];
    numberJpsi_ini_time = new int[nTimeSteps+1];
    numberJpsi_midPseudoRap_all_time = new int[nTimeSteps+1];
    numberJpsi_midPseudoRap_ini_time = new int[nTimeSteps+1];
    numberJpsi_forwardPseudoRap_all_time = new int[nTimeSteps+1];
    numberJpsi_forwardPseudoRap_ini_time = new int[nTimeSteps+1];
    numberJpsi_midNormRap_all_time = new int[nTimeSteps+1];
    numberJpsi_midNormRap_ini_time = new int[nTimeSteps+1];
    numberJpsi_forwardNormRap_all_time = new int[nTimeSteps+1];
    numberJpsi_forwardNormRap_ini_time = new int[nTimeSteps+1];
  //   numberJpsi_midSpaceTimeRap_time = new int[nTimeSteps+1];
    numberJpsiProd_time = new int[nTimeSteps+1];
    numberJpsiDiss_time = new int[nTimeSteps+1];
    numberJpsiDissTd_time = new int[nTimeSteps+1];
    numberCCbGG_time = new int[nTimeSteps+1];
    timestepAnalysed = new bool[nTimeSteps+1];
  //   charmJetEnergy = new double[nTimeSteps+1];

    for ( int i = 0; i < nTimeSteps + 1; i++ )
    {
      numberJpsi_all_time[i] = numberJpsi_ini_time[i] = numberJpsiDiss_time[i] = numberJpsiDissTd_time[i] = 
      numberJpsi_midPseudoRap_all_time[i] = numberJpsi_midPseudoRap_ini_time[i] = numberJpsi_forwardPseudoRap_all_time[i] = numberJpsi_forwardPseudoRap_ini_time[i] = 
      numberJpsi_midNormRap_all_time[i] = numberJpsi_midNormRap_ini_time[i] = numberJpsi_forwardNormRap_all_time[i] = numberJpsi_forwardNormRap_ini_time[i] = 
      numberJpsiProd_time[i] = numberCCbGG_time[i] = 0;
      timestepAnalysed[i] = false;
  //     charmJetEnergy[i] = 0.0;
    }
    
    theInterpolation_nJpsi.configure();
  }
  
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
}


void analysis::handle_output_studies( OUTPUT_SCHEME _outputScheme )
{
  // By default switches are set off (please do not edit here, but define output scheme and change below, for more information see https://th.physik.uni-frankfurt.de/~bamps/cgi-bin/trac/wiki/RestrictedWiki/AnalysisOutputScheme):
  studyHQ = false; // heavy quarks
  studyJpsi = false; // jpsi
  studyTempInTube = false; // temperature in a centred tube in heavy ion collision
  studyTempAndVelocity = false; // write out temperature and velocity for each cell
  studyPtSpectra = false; // transverse momentum pt spectra (partly the same as in v2RAA class)
  studyEtSpectra = false; // transverse energy spectra
  studyYDistribution = false; // rapidity distribution
  studyJets = false; // study jets
  studyCentralDensity = false; // density in central part of collision
  studyBackground = false; // print also properties like v2, RAA of background
  
  
  //---- defining standard rapidity ranges ----
  // only use positiv ranges since the investigated collision systems usually are symmetric in +-y and we therefore only compare the absolute value of y
  analysisRapidityRange yRange;
  yRange.reset( 0, 0.5 );
  rapidityRanges.push_back(yRange);
  yRange.reset( 0, 0.8 );
  rapidityRanges.push_back(yRange);
  yRange.reset( 0, 1.0 );
  rapidityRanges.push_back(yRange);
  yRange.reset( 0, 1.5 );
  rapidityRanges.push_back(yRange);
  yRange.reset( 0, 2.0 );
  rapidityRanges.push_back(yRange);
  yRange.reset( 0, 2.5 );
  rapidityRanges.push_back(yRange);  
  yRange.reset( 0, infinity );
  rapidityRanges.push_back(yRange);  
  //---- defining rapidity ranges ----

  
  // add a new case for your outpute scheme which you can create in configuration.h
  switch ( _outputScheme )
  {
    case phenix_hq_electrons:
      studyHQ = true;
      
      rapidityRanges.clear();
      yRange.reset( 0, 0.35 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.5 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 1.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 2.0 );
      rapidityRanges.push_back(yRange);
      break;
    case alice_hq_electrons:
      studyHQ = true;
      
      rapidityRanges.clear();
      yRange.reset( 0, 0.8 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 1.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 2.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.5 );
      rapidityRanges.push_back(yRange);
      break;
    case alice_hq_muons:
      studyHQ = true;

      rapidityRanges.clear();
      yRange.reset( 2.5, 4.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.35 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.5 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.8 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 1.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 2.0 );
      rapidityRanges.push_back(yRange);
      break;
    case alice_hq_dmesons:
      studyHQ = true;
      
      rapidityRanges.clear();
      yRange.reset( 0, 0.5 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.8 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 1.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 2.0 );
      rapidityRanges.push_back(yRange);
      break;
    case cms_hq_nonPromptJpsi:
      studyHQ = true;
      
      rapidityRanges.clear();
      yRange.reset( 0, 2.4 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.35 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.5 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 1.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 2.0 );
      rapidityRanges.push_back(yRange);
      break;
    case phenix_jpsi:
      studyJpsi = true;
      studyHQ = true;
      studyTempInTube = true;

      rapidityRanges.clear();
      yRange.reset( 0, 0.35 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 1.2, 2.2 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.5 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 1.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 2.0 );
      rapidityRanges.push_back(yRange);
      break;
    case alice_jpsi:
      studyJpsi = true;
      studyHQ = true;
      studyTempInTube = true;

      rapidityRanges.clear();
      yRange.reset( 0, 0.9 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 2.5, 4.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.35 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.5 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 1.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 2.0 );
      rapidityRanges.push_back(yRange);
      break;
    case cms_jpsi:
      studyJpsi = true;
      studyHQ = true;
      studyTempInTube = true;

      rapidityRanges.clear();
      yRange.reset( 0, 2.4 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.9 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 2.5, 4.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.35 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.5 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 1.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 2.0 );
      rapidityRanges.push_back(yRange);
      break;
    default:
      break;
  }

}



void analysis::collectPtDataInitial()
{
  if( studyPtSpectra )
  {
    ptDistribution( gluon, addedParticles, addedParticles.size(), 0 );
    ptDistribution( light_quark, addedParticles, addedParticles.size(), 0 );
    ptDistribution( allFlavors, addedParticles, addedParticles.size(), 0 );
    ptSoftDistribution( gluon, particles_atTimeNow, particles_atTimeNow.size(), 0 );
    ptSoftDistribution( light_quark, particles_atTimeNow, particles_atTimeNow.size(), 0 );
    ptSoftDistribution( allFlavors, particles_atTimeNow, particles_atTimeNow.size(), 0 );
  }
}



void analysis::collectPtData( const int step )
{
  //are issued starting with 0
  //0 is reserved for initial output, thus step is incremented
  if( studyPtSpectra )
  {
    ptDistribution( gluon, addedParticles, addedParticles.size(), step + 1 );
    ptDistribution( light_quark, addedParticles, addedParticles.size(), step + 1 );
    ptDistribution( allFlavors, addedParticles, addedParticles.size(), step + 1 );
    ptSoftDistribution( gluon, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
    ptSoftDistribution( light_quark, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
    ptSoftDistribution( allFlavors, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
  }
}


void analysis::collectYDataInitial()
{
  if( studyYDistribution )
    collectYData( -1 );
}


void analysis::collectYData( const int step )
{
  //are issued starting with 0
  //0 is reserved for initial output, thus step is incremented
  if( studyYDistribution )
  {
    yDistribution( gluon, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
    yDistribution( up, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
    yDistribution( down, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
    yDistribution( strange, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
    yDistribution( anti_up, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
    yDistribution( anti_down, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
    yDistribution( anti_strange, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
  }
}


void analysis::collectEtDataInitial()
{
  if( studyEtSpectra )
    collectEtData( -1 );
}


void analysis::collectEtData( const int step )
{
  //are issued starting with 0
  //0 is reserved for initial output, thus step is incremented
  if( studyEtSpectra )
  {
    transverseEnergyDistribution( gluon, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
    transverseEnergyDistribution( light_quark, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
    transverseEnergyDistribution( anti_light_quark, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
  }
}





void analysis::initialOutput()
{
  if ( v2output )
    computeV2RAA( "initial", 0  );

  if ( studyJpsi ) // to consider charm annihaltion is just useful if added particles can scatter
  { 
    jpsiEvolution( 0 );
    ini_charm_correlations();
    writeJpsiFugacityOutput( 0 );
  }
  
  if( studyTempInTube)
     writeTempInTube( 0 );
  
  if ( dndyOutput )
    print_dndy( "initial" );
  
  if ( studyParticleOutput )
    particleOutput( 0 );
  
}



void analysis::intermediateOutput( const int nn )
{
  string name;
  stringstream ss;

  ss << tstep[nn];
  name = ss.str() + "fm";

  if ( v2output && v2outputIntermediateSteps )
    computeV2RAA( name, tstep[nn] );

  if ( studyJpsi ) // to consider charm annihaltion is just useful if added particles can scatter
  {
    jpsiEvolution( nn + 1 );
    writeJpsiFugacityOutput( nn + 1 );
  }
  
  if( studyTempInTube)
     writeTempInTube( nn + 1 ); 
  
  if ( studyTempAndVelocity ) // hydro
    writeTempAndVel( nn + 1  );
  
//   if ( charmTestJet )
//    analyseCharmTestJetEvolution( nn + 1 );
  
  if ( dndyOutput )
    print_dndy( name );
}



void analysis::finalOutput( const double _stoptime )
{
  
  if( studyPtSpectra )
  {  
    printPtSpectra( gluon );
    printPtSpectra( light_quark );
    printPtSpectra( allFlavors );
    printSoftPtSpectra( gluon );
    printSoftPtSpectra( light_quark );
    printSoftPtSpectra( allFlavors );
  }
  
  if ( studyParticleOutput )
    particleOutput( nTimeSteps );
  
  if( studyYDistribution )
    printYDistribution();
  
  if ( v2output )
    computeV2RAA( "final", _stoptime );

  if ( particleCorrelationsOutput )
  {
    onePartclCorrelations();
    twoPartclCorrelations();
  }

  if ( studyJpsi ) // to consider charm annihaltion is just useful if added particles can scatter
  {
    printJpsiEvolution();
//     jpsi_correlations();
  }
  
//   if ( charmTestJet )
//    analyseCharmTestJet();
  
  if( hadronization_hq && mesonDecay )
    analyseAngleDe();
  
  if ( dndyOutput )
    print_dndy( "final" );
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

  // populate temporary _ptBins with temporary "references" to the actual data structures according to the requested flavor
  tVecOfArrayOfDoubleVec _ptBins;
  for ( unsigned int yRange = 0; yRange < rapidityRanges.size(); yRange++ )
  {
    switch ( _flavTypeToComputeFor )
    {
      case gluon:
        _ptBins.push_back( ptBins_gluons[yRange] );
        break;
      case light_quark:
        _ptBins.push_back( ptBins_quarks[yRange] );
        break;
      case up:
        _ptBins.push_back( ptBins_ups[yRange] );
        break;
      case down:
        _ptBins.push_back( ptBins_downs[yRange] );
        break;
      case strange:
        _ptBins.push_back( ptBins_stranges[yRange] );
        break;
      case anti_up:
        _ptBins.push_back( ptBins_anti_ups[yRange] );
        break;
      case anti_down:
        _ptBins.push_back( ptBins_anti_downs[yRange] );
        break;  
      case anti_strange:
        _ptBins.push_back( ptBins_anti_stranges[yRange] );
        break;
      case allFlavors:
        _ptBins.push_back(ptBins_all[yRange]);
        break;
      default:
        string errMsg = "error in ptDistribution, flavor not specified";
        throw eAnalysis_error( errMsg );
        break;
    }
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

  //------------- the actual output ---------------
  for ( unsigned int yRangeIndex = 0; yRangeIndex < rapidityRanges.size(); yRangeIndex++ )
  {
    file << "#y in +- [" << rapidityRanges[yRangeIndex].yleft << ", " << rapidityRanges[yRangeIndex].yright << "]" << endl;
    for ( int i = 0; i < numberBinsPT; i++ )
    {
      file << ptBinLabels[i] << sep;
      for ( int j = 0; j <= nTimeSteps; j++ )
      {
        if ( j == 0 || j == nTimeSteps )
          file << _ptBins[yRangeIndex][j][i] << sep;
        else if ( tstep[j-1] <= theConfig->getRuntime() )
          file << _ptBins[yRangeIndex][j][i] << sep;
      }
      file << endl;
    }
    file << endl << endl;
  }
  //-------------------------------------------------

  file.close();
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
  double pt_fin, pt_init, eta, dpt_vec, dpt_scal, rt_init, dE;


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

  for ( unsigned int i = 0;i < addedParticles.size();i++ )
  {
//     if(addedParticles[i].T <= time)
//     {
    eta = addedParticles[i].Mom.Pseudorapidity( addedParticles[i].m );

    if ( fabs( eta ) <= eta_max )
    {
      pt_init = addedParticles[i].MomInit.Pt();
      pt_fin  = addedParticles[i].Mom.Pt();

      dpt_vec = (addedParticles[i].MomInit - addedParticles[i].Mom).Pt();

      dpt_scal = pt_init - pt_fin;
      rt_init  = addedParticles[i].PosInit.Pt();

      dE = addedParticles[i].MomInit.E() - addedParticles[i].Mom.E();

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
  double eta, rt_init;
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
  for ( unsigned int i = 0;i < addedParticles.size();i++ )
  {
//     if(addedParticles[i].T <= time)
//     {
    eta = addedParticles[i].Mom.Pseudorapidity( addedParticles[i].m );

    if ( fabs( eta ) <= eta_max )
    {
      for ( unsigned int j = i + 1;j < addedParticles.size();j++ )
      {
//         if ( addedParticles[i].N_EVENT_pp == addedParticles[j].N_EVENT_pp )
        if ( true )
        {

          //           if(addedParticles[i].T <= time)
          //           {
          eta = addedParticles[j].Mom.Pseudorapidity( addedParticles[j].m );

          if ( fabs( eta ) <= eta_max )
          {
//                 dpt_ini = sqrt( pow(addedParticles[i].MomInit.Px()-addedParticles[j].MomInit.Px(),2.0) + pow(addedParticles[i].MomInit.Py()-addedParticles[j].MomInit.Py(),2.0) );
//                 dpt_fin = sqrt( pow(addedParticles[i].Mom.Px()-addedParticles[j].Mom.Px(),2.0) + pow(addedParticles[i].Mom.Py()-addedParticles[j].Mom.Py(),2.0) );
            dpt_ini = fabs( addedParticles[i].MomInit.Pt() - addedParticles[j].MomInit.Pt() );
            dpt_fin = fabs( addedParticles[i].Mom.Pt() - addedParticles[j].Mom.Pt() );

            dp_ini = sqrt( (addedParticles[i].MomInit - addedParticles[j].MomInit).vec2() );
            dp_fin = sqrt( (addedParticles[i].Mom - addedParticles[j].Mom).vec2() );

            rt_init =  addedParticles[i].PosInit.Pt(); // same for both particles

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
//                   partdpt << fabs( sqrt( pow(addedParticles[i].Mom.Px(),2.0) + pow(addedParticles[i].Mom.Py(),2.0) ) - sqrt( pow(addedParticles[j].Mom.Px(),2.0) + pow(addedParticles[j].Mom.Py(),2.0) ) );
//
//                   partdpt.width(25);
//                   partdpt << addedParticles[i].Mom.Px();
//                   partdpt.width(15);
//                   partdpt << addedParticles[i].Mom.Py();
//                   partdpt.width(15);
//                   partdpt << addedParticles[i].Mom.Pz();
//
//                   partdpt.width(25);
//                   partdpt << addedParticles[j].Mom.Px();
//                   partdpt.width(15);
//                   partdpt << addedParticles[j].Mom.Py();
//                   partdpt.width(15);
//                   partdpt << addedParticles[j].Mom.Pz();
//
//                   partdpt.width(25);
//                   partdpt << addedParticles[i].MomInit.Px();
//                   partdpt.width(15);
//                   partdpt << addedParticles[i].MomInit.Py();
//                   partdpt.width(15);
//                   partdpt << addedParticles[i].MomInit.Pz();
//
//                   partdpt.width(25);
//                   partdpt << addedParticles[j].MomInit.Px();
//                   partdpt.width(15);
//                   partdpt << addedParticles[j].MomInit.Py();
//                   partdpt.width(15);
//                   partdpt << addedParticles[j].MomInit.Pz();
//
//                   partdpt.width(25);
//                   partdpt << addedParticles[j].PosInit.X();
//                   partdpt.width(15);
//                   partdpt << addedParticles[j].PosInit.Y();
//                   partdpt.width(15);
//                   partdpt << addedParticles[j].PosInit.Z();
//
//                   partdpt.width(25);
//                   partdpt << rt_init << endl;
//
//                 }


            // azimuthal angle between two charm quarks
            dphi_ini = acos( CosPhi( addedParticles[i].MomInit, addedParticles[j].MomInit ) );
            dphi_fin = acos( CosPhi( addedParticles[i].Mom, addedParticles[j].Mom ) );

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

            if ( addedParticles[i].Mom.Pt2() >= pow( pt_min_assoc, 2.0 ) && addedParticles[j].Mom.Pt2() >= pow( pt_min_assoc, 2.0 ) )
            {
              dphiF22.add( dphi_fin );
              dphiI22.add( dphi_ini );

              if (( addedParticles[i].Mom.Pt2() >= pow( pt_min_trig, 2.0 ) && addedParticles[j].Mom.Pt2() >= pow( pt_min_assoc, 2.0 )) ||
                  ( addedParticles[i].Mom.Pt2() >= pow( pt_min_assoc, 2.0 ) && addedParticles[j].Mom.Pt2() >= pow( pt_min_trig, 2.0 )))
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



            pt_fin_i = addedParticles[i].Mom.Pt();
            pt_fin_j = addedParticles[j].Mom.Pt();

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



  for ( int j = 0;j < addedParticles.size();j++ )
  {
    eta = addedParticles[j].Mom.Pseudorapidity( addedParticles[j].m );

    VectorXYZ V = VectorXYZ( 0.0, 1.0, 0.0 );
    if ( fabs( eta ) <= eta_max )
    {
      dphi_fin = acos( CosPhi(addedParticles[j].Mom, V) );
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
  v2RAA theV2RAA( theConfig, name, filename_prefix, rapidityRanges );
  
  if( theConfig->isStudyNonPromptJpsiInsteadOfElectrons() )
  {
    theV2RAA.setPtBinProperties( 0.0, 30.0, 60 );
    
    theV2RAA.computeFor( bottom, addedParticles, addedParticles.size(), "added", _outputTime, v2jets );
    if( mesonDecay )
      theV2RAA.computeFor( electron_gen, addedPartcl_electron, addedPartcl_electron.size(), "added", _outputTime, v2jets );
  }
  if( studyHQ || studyJpsi )
  {
    theV2RAA.setPtBinProperties( 0.0, 30.0, 60 );
    
    theV2RAA.computeFor( charm, addedParticles, addedParticles.size(), "added", _outputTime, v2jets );
    theV2RAA.computeFor( bottom, addedParticles, addedParticles.size(), "added", _outputTime, v2jets );
    theV2RAA.computeFor( heavy_quark, addedParticles, addedParticles.size(), "added", _outputTime, v2jets );

    if ( name == "initial" || name == "final" )
    {
      if( hadronization_hq )
      {
        theV2RAA.computeFor( dmeson_gen, addedParticlesCopy, addedParticlesCopy.size(), "added", _outputTime, v2jets );
        theV2RAA.computeFor( bmeson_gen, addedParticlesCopy, addedParticlesCopy.size(), "added", _outputTime, v2jets );
      }
      
      if( mesonDecay )
      {
        theV2RAA.computeFor( electron_gen, addedPartcl_electron, addedPartcl_electron.size(), "added", _outputTime, v2jets );
        theV2RAA.computeFor( c_electron, addedPartcl_electron, addedPartcl_electron.size(), "added", _outputTime, v2jets );  // electrons from charm quarks (D mesons actually)
        theV2RAA.computeFor( b_electron, addedPartcl_electron, addedPartcl_electron.size(), "added", _outputTime, v2jets );  // electrons from bottom quarks (B mesons actually)
      }
    }
    
    // also take a look at light parton v2 of background
    theV2RAA.computeFor( gluon, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
    theV2RAA.computeFor( light_quark, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
    
    if ( studyJpsi )
    {
      theV2RAA.setPtBinProperties( 0.0, 15.0, 20 );
      
      theV2RAA.computeFor( jpsi, addedParticles, addedParticles.size(), "added", _outputTime, v2jets );
      theV2RAA.computeFor( jpsi_ini, addedParticles, addedParticles.size(), "added", _outputTime, v2jets );
      theV2RAA.computeFor( jpsi_sec, addedParticles, addedParticles.size(), "added", _outputTime, v2jets );
    }
  }
  else
  {
    theV2RAA.computeFor( gluon, addedParticles, addedParticles.size(), "jets", _outputTime, v2jets );
    theV2RAA.computeFor( light_quark, addedParticles, addedParticles.size(), "jets", _outputTime, v2jets );
    theV2RAA.computeFor( up, addedParticles, addedParticles.size(), "jets", _outputTime, v2jets );
    theV2RAA.computeFor( down, addedParticles, addedParticles.size(), "jets", _outputTime, v2jets );
    theV2RAA.computeFor( strange, addedParticles, addedParticles.size(), "jets", _outputTime, v2jets );
    theV2RAA.computeFor( anti_up, addedParticles, addedParticles.size(), "jets", _outputTime, v2jets );
    theV2RAA.computeFor( anti_down, addedParticles, addedParticles.size(), "jets", _outputTime, v2jets );
    theV2RAA.computeFor( anti_strange, addedParticles, addedParticles.size(), "jets", _outputTime, v2jets );
  }
  
  if( studyBackground )
  {
    theV2RAA.computeFor( gluon, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
    theV2RAA.computeFor( up, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
    theV2RAA.computeFor( down, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
    theV2RAA.computeFor( strange, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
    theV2RAA.computeFor( anti_up, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
    theV2RAA.computeFor( anti_down, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
    theV2RAA.computeFor( anti_strange, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
    theV2RAA.computeFor( light_quark, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
  }
}



v2RAA::v2RAA( config * const c, string name_arg, string filename_prefix_arg, std::vector<analysisRapidityRange> rapidityRanges_arg, const double pt_min_arg, const double pt_max_arg, const int n_g_arg, const double pt_min_background_arg, const double pt_max_background_arg, const int n_g_background_arg ):
    theConfig( c ), name( name_arg ), filename_prefix( filename_prefix_arg ), rapidityRanges( rapidityRanges_arg ), pt_min( pt_min_arg ), pt_max( pt_max_arg ), n_g( n_g_arg ), pt_min_background( pt_min_arg ), pt_max_background( pt_max_arg ), n_g_background( n_g_arg )
{
  eta_bins = rapidityRanges.size();
}




void v2RAA::computeFor( const FLAVOR_TYPE _flavTypeToComputeFor, vector<ParticleOffline>& _particles, const int n_particles, string additionalNameTag, const double _outputTime, const v2Type _v2type )
{
  double eta, pp, pt, v2, xt;
  double sinAlpha, alpha;
  int dummy, flavor, n_bins;
  int alphaIndex;
  
  double _pt_min, _pt_max;

  string filename_v2, filename_v2_summed, filename_v2_tot, filename_yield, filename_pt_angleDependence, type;
  
  type = Particle::getName( _flavTypeToComputeFor );
  
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
    
    // for some scenarios however explicitly the rapidity is measured. So substitute eta by the rapidity:
    if( ( theConfig->isStudyNonPromptJpsiInsteadOfElectrons() &&
          ( _flavTypeToComputeFor == charm ||
            _flavTypeToComputeFor == bottom ||
            _flavTypeToComputeFor == heavy_quark ||
            _flavTypeToComputeFor == c_electron ||
            _flavTypeToComputeFor == b_electron ||
            _flavTypeToComputeFor == electron_gen 
          ) ) ||
        ParticleOffline::mapToGenericFlavorType( _flavTypeToComputeFor ) == dmeson_gen ||
        ParticleOffline::mapToGenericFlavorType( _flavTypeToComputeFor ) == bmeson_gen ||
        ParticleOffline::mapToGenericFlavorType( _flavTypeToComputeFor ) == jpsi
    )
      eta = _particles[i].Mom.Rapidity();

    v2 = ( pow( _particles[i].Mom.Px(), 2.0 ) - pow( _particles[i].Mom.Py(), 2.0 ) ) / pow( pt, 2.0 );

    flavor = _particles[i].FLAVOR;
    
    FLAVOR_TYPE mother_flav;
    if( _flavTypeToComputeFor == c_electron || _flavTypeToComputeFor == b_electron )
    {
      int mother_id = i / theConfig->getNumberElectronStat();
      mother_flav = addedParticlesCopy[mother_id].FLAVOR;
    }

    if( ( pt <= _pt_max && pt > _pt_min ) &&  
        ( ( _flavTypeToComputeFor == flavor ) || 
          ( _flavTypeToComputeFor == ParticleOffline::mapToGenericFlavorType( static_cast<FLAVOR_TYPE>( flavor ) ) ) ||
          ( _flavTypeToComputeFor == charm && ( flavor == anti_charm ) ) ||
          ( _flavTypeToComputeFor == bottom && ( flavor == anti_bottom ) ) ||
          ( _flavTypeToComputeFor == heavy_quark && ( flavor == charm || flavor == bottom || flavor == anti_charm || flavor == anti_bottom ) ) ||
          ( _flavTypeToComputeFor == c_electron && ( ( flavor == electron || flavor == positron ) && ParticleOffline::mapToGenericFlavorType( static_cast<FLAVOR_TYPE>( mother_flav ) ) == dmeson_gen ) ) ||
          ( _flavTypeToComputeFor == b_electron && ( ( flavor == electron || flavor == positron ) && ParticleOffline::mapToGenericFlavorType( static_cast<FLAVOR_TYPE>( mother_flav ) ) == bmeson_gen ) ) ||
          ( _flavTypeToComputeFor == jpsi_ini && ( flavor == jpsi && _particles[i].initially_produced ) ) ||
          ( _flavTypeToComputeFor == jpsi_sec && ( flavor == jpsi && !_particles[i].initially_produced ) )
        ) 
      )
    {
      // individually check for each rapidity range whether this particle needs to be binned
      for ( int yRangeIndex = 0; yRangeIndex < eta_bins; yRangeIndex++ )
      {
        if ( fabs( eta ) >= rapidityRanges[yRangeIndex].yleft && fabs( eta ) <= rapidityRanges[yRangeIndex].yright )
        {
          v2sum[yRangeIndex] += v2;
          NmbInRange[yRangeIndex]++;

          dummy = int(( log( pt ) - log( _pt_min ) ) / d_ln_pt );
          ptBinsV2[yRangeIndex][dummy] += v2;
          ptBinsNmb[yRangeIndex][dummy]++;
          ptBinsAngleDep[yRangeIndex][alphaIndex][dummy]++;
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
  print_v2_tot << "# t = " << _outputTime << endl;
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
  print_v2_summed << "# t = " << _outputTime << endl;
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
  print_yield << "# t = " << _outputTime << endl;
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
      
      double nInBin = double( ptBinsNmb[i][k] ) / theConfig->getTestparticles() / dpt / delta_eta;
      
      if( _v2type == v2jets )
        nInBin = nInBin / theConfig->getNaddedEvents();
      
      if( Particle::mapToGenericFlavorType( _flavTypeToComputeFor ) == jpsi )
        nInBin = nInBin / theConfig->getJpsiTestparticles();
      
      if( Particle::mapToGenericFlavorType( _flavTypeToComputeFor ) == electron_gen )
        nInBin = nInBin / theConfig->getNumberElectronStat();
      
      print_yield.width( 15 );
      print_yield << nInBin;
    }
    print_yield << endl;
  }
  
  
  // print yield for RAA for different angles with respect to the reaction plane
  print_pt_angleDependence << "# " << type << " yield distribution for different angles (alpha) with respect to the reaction plane" << endl;
  print_pt_angleDependence << "# t = " << _outputTime << endl;
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
        print_pt_angleDependence << double( ptBinsAngleDep[i][j][k] ) / theConfig->getTestparticles() / dpt / delta_eta;
      }
      print_pt_angleDependence << endl;
    }    
    print_pt_angleDependence << endl;
    print_pt_angleDependence << endl;
  }
  
    //CMS pt cut for non prompt jpsi: 6.5GeV < pt < 30 GeV, y cut: |y|<2.4
  if( theConfig->isStudyNonPromptJpsiInsteadOfElectrons() && _flavTypeToComputeFor == electron_gen )
  {
    string filename_yield_cmsNonPromptJpsi = filename_prefix + "_nonPromptJpsiCMScuts_yield_" + name;
    fstream print_cmsNonPromptJpsi( filename_yield_cmsNonPromptJpsi.c_str(), ios::out | ios::trunc );
    
    print_cmsNonPromptJpsi << "# yield of non prompt jpsi from B meson decays" << endl;
    print_cmsNonPromptJpsi << "# with same acceptance cuts as for CMS:  6.5GeV < pt < 30 GeV,  |y|<2.4" << endl;
    
    int count_jpsi = 0;
    double delta_eta_cms = 2.0 * 2.4;
    double y, pt;
    for ( int i = 1; i <= n_particles; i++ )
    {
      pt = _particles[i].Mom.Pt();
      y  = _particles[i].Mom.Rapidity();
      
      if( fabs( y ) < 2.4 && pt > 6.5 && pt < 30.0 )
        count_jpsi++;
    }
    
    print_cmsNonPromptJpsi << double( count_jpsi ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getNumberElectronStat() / delta_eta_cms << endl;
  }
  

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
//     y = 0.5*log( (particles_atTimeNow[i].E + particles_atTimeNow[i].PZ) / (particles_atTimeNow[i].E - particles_atTimeNow[i].PZ) );
//     if(fabs(y) <= y_max)
//     {
//       zBins.add( fabs(particles_atTimeNow[i].Z) );
//       xt = sqrt(pow(particles_atTimeNow[i].X,2.0) + pow(particles_atTimeNow[i].Y,2.0));
//       tBins.add(xt);
//     }

    if ( particles_atTimeNow[i].Pos.T() <= tstep[step] )
    {
      y = particles_atTimeNow[i].Mom.Rapidity();
      if ( fabs( y ) <= y_max )
      {
        zBins.add( fabs( particles_atTimeNow[i].Pos.Z() ) );
        tBins.add( particles_atTimeNow[i].Pos.Perp() );
      }



//       zBins.add( fabs(particles_atTimeNow[i].Pos.Z()) );
//       tBins.add( particles_atTimeNow[i].Pos.Perp() );
    }

  }

  zBins.print();
  tBins.print();
}


void analysis::ptDistribution( const FLAVOR_TYPE _flavTypeToComputeFor, vector<ParticleOffline>& _particles, const int n_particles, const int step )
{
  FLAVOR_TYPE genFlavor;
  double pt, y;
  tVecOfArrayOfDoubleVec _ptBins;
  
  // populate temporary _ptBins with temporary "references" to the actual data structures according to the requested flavor
  for ( int yRange = 0; yRange < rapidityRanges.size(); yRange++ )
  {
    switch ( _flavTypeToComputeFor )
    {
      case gluon:
        _ptBins.push_back( ptBins_gluons[yRange] );
        break;
      case light_quark:
        _ptBins.push_back( ptBins_quarks[yRange] );
        break;
      case up:
        _ptBins.push_back( ptBins_ups[yRange] );
        break;
      case down:
        _ptBins.push_back( ptBins_downs[yRange] );
        break;
      case strange:
        _ptBins.push_back( ptBins_stranges[yRange] );
        break;
      case anti_up:
        _ptBins.push_back( ptBins_anti_ups[yRange] );
        break;
      case anti_down:
        _ptBins.push_back( ptBins_anti_downs[yRange] );
        break;  
      case anti_strange:
        _ptBins.push_back( ptBins_anti_stranges[yRange] );
        break;
      case allFlavors:
        _ptBins.push_back(ptBins_all[yRange]);
        break;
      default:
        string errMsg = "error in ptDistribution, flavor not specified";
        throw eAnalysis_error( errMsg );
        break;
    }
  }
  
  // loop over all particles and bin them according to their pt
  for ( int j = 0; j < n_particles; j++ )
  {
    pt = _particles[j].Mom.Pt();
    y = _particles[j].Mom.Rapidity();
    
    // check whether particle has the correct flavor
    genFlavor = ParticleOffline::mapToGenericFlavorType( _particles[j].FLAVOR );
    if ( _flavTypeToComputeFor == allFlavors || _particles[j].FLAVOR == _flavTypeToComputeFor || 
      ( _flavTypeToComputeFor == light_quark && ( genFlavor == light_quark || genFlavor == anti_light_quark ) ) )
    {
      // is it in the possible pt range at all?
      if ( pt < maxPT && pt >= minPT )
      {
        // individually check for each rapidity range whether this particle needs to be binned
        for ( int yRangeIndex = 0; yRangeIndex < rapidityRanges.size(); yRangeIndex++ )
        {
          if ( fabs( y ) >= rapidityRanges[yRangeIndex].yleft && fabs( y ) <= rapidityRanges[yRangeIndex].yright )
          {
            if ( pt == minPT )  // a special case
            {
              ++_ptBins[yRangeIndex][step][0];  // actually bin the pt of the particle
            }
            else
            {
              ++_ptBins[yRangeIndex][step][int(( pt - minPT )/binWidthPT )];  // actually bin the pt of the particle
            }
          }
        }
      }
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
    pt = _particles[j].Mom.Pt();

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
    y = _particles[j].Mom.Rapidity();

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
    Et = _particles[j].Mom.Pt();
    y = _particles[j].Mom.Rapidity();
    
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

  tempEvent.R_proj = addedParticles[jetID].Pos;
  tempEvent.P_proj_out = addedParticles[jetID].Mom;

  vector<jetTrackerSingleEvent> tempVec;
  tempVec.push_back( tempEvent );
  jetTracker.push_back( tempVec );
}


void analysis::addJetEvents_final()
{
  double pt;
  for ( int i = 0; i < addedParticles.size(); i++ )
  {
    pt = addedParticles[i].Mom.Pt();
    if ( pt > jetTracking_PT )
    {
      jetTrackerSingleEvent tempEvent;
      tempEvent.jet_ID_in = -addedParticles[i].unique_id;
      tempEvent.jet_ID_out = -1;
      tempEvent.R_proj = addedParticles[i].Pos;
      tempEvent.P_proj_in = addedParticles[i].Mom;
      
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
  tempEvent.R_proj = addedParticles[jetID].Pos;
  tempEvent.flavor_in = static_cast<int>( addedParticles[jetID].FLAVOR );

  tempEvent.P_proj_in = addedParticles[jetID].Mom;
  tempEvent.P1_in = particles_atTimeNow[partner1].Mom;

  if ( coll_type == c3to2 )
  {
    tempEvent.P2_in = particles_atTimeNow[partner2].Mom;
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
  double pt1 = addedParticles[added_ID].Mom.Pt();
  double pt2 = particles_atTimeNow[ID_2].Mom.Pt();

  double pt3 = -1;
  if ( coll_type == c2to3 )
  {
    pt3 = addedParticles[ID_3].Mom.Pt();
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
    jetTracker[entity_ID].back().P_proj_out = addedParticles[jetID].Mom;
    
    jetTracker[entity_ID].back().flavor_out = static_cast<int>( addedParticles[jetID].FLAVOR );
    jetTracker[entity_ID].back().P1_out = particles_atTimeNow[partner1].Mom;

    if ( coll_type == c2to3 )
    {
      jetTracker[entity_ID].back().P2_out = addedParticles[partner2].Mom;
    }
  }
  else if ( addedParticles[jetID].Mom.Pt2() > pow( jetTracking_PT, 2 ) )
  {
    jetTrackerSingleEvent tempEvent;
    tempEvent.jet_ID_in = -1;
    tempEvent.jet_ID_out = -addedParticles[jetID].unique_id;
    tempEvent.coll_type = production;
    tempEvent.R_proj = addedParticles[jetID].Pos;
    tempEvent.flavor_in = -1;
    tempEvent.flavor_out = static_cast<int>( addedParticles[jetID].FLAVOR );

    tempEvent.P_proj_out = addedParticles[jetID].Mom;
    tempEvent.P1_out = particles_atTimeNow[partner1].Mom;

    if ( coll_type == c2to3 )
    {
      tempEvent.P2_out = addedParticles[partner2].Mom;
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
  string filename = filename_prefix + "_" + name + ".f1";
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
    file << i << sep << addedParticles[i].unique_id << sep << addedParticles[i].cell_id << sep << addedParticles[i].FLAVOR << sep 
         << addedParticles[i].Pos.T() << sep << addedParticles[i].Pos.X() << sep
         << addedParticles[i].Pos.Y() << sep << addedParticles[i].Pos.Z() << sep 
         << addedParticles[i].Mom.E() << sep << addedParticles[i].Mom.Px() << sep 
         << addedParticles[i].Mom.Py() << sep << addedParticles[i].Mom.Pz() << sep 
         << addedParticles[i].md2g << sep << addedParticles[i].md2q << endl;
  }
  file.close();
}



void analysis::writePartclMovie( vector< ParticleOffline >& _particles, const int n_particles, fstream& _oscar, const int step, const int jumpSteps )
{
  const string sep = "  ";
  const int width = 14;
  double cc, dt;
  VectorTXYZ Pos;
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
                 << setw( width ) << _particles[i].Mom.Px() 
                 << setw( width ) << _particles[i].Mom.Py() 
                 << setw( width ) << _particles[i].Mom.Pz() 
                 << setw( width ) << _particles[i].Mom.E() 
                 << setw( width ) << _particles[i].m;

          if ( _particles[i].Pos.T() < time )
          {
            cout << "error in write movie particle data, particles from the past" << endl;
            cout << "timestep=" << time << "   time _particles=" << _particles[i].Pos.T() << endl;
          }
          
          Pos = _particles[i].Pos + _particles[i].Mom * (( time - _particles[i].Pos.T() )/_particles[i].Mom.E() );


          _oscar << setw( width ) << Pos.X() << setw( width ) << Pos.Y() << setw( width ) << Pos.Z() << setw( width ) << Pos.T();

          //       _oscar << "0.0" << sep << "0.0" << sep << "0.0" << sep << "0.0" << sep << "0.0" << sep << "0" << endl;
          //         _oscar << zero << sep << zero << sep << zero << sep << zero << sep << zero << ::std::setw( 10 ) << int(zero) << endl;
          //         _oscar << zero << sep << _particles[i].X_lastInt << sep << _particles[i].Y_lastInt << sep << _particles[i].Z_lastInt << sep << _particles[i].T_lastInt << ::std::setw( 10 ) << int(zero) << endl;
          _oscar << setw( width ) << zero 
                 << setw( width ) << _particles[i].lastInt.X()
                 << setw( width ) << _particles[i].lastInt.Y() 
                 << setw( width ) << _particles[i].lastInt.Z() 
                 << setw( width ) << _particles[i].lastInt.T() 
                 << ::std::setw( 10 ) << -(_particles[i].unique_id) << endl;
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





// heavy quark stuff

void analysis::writeJpsiFugacityOutput( const int step )
{
  const string sep = "  ";

  double fugacity, temp, deltaTemp, enDen;
  int n_jpsi;
  double dr, dz, deta;


  double time;
  if ( step == 0 )
    time = 0.0;
  else if ( step == nTimeSteps )
    return;
  else
    time = tstep[step-1];

  if ( step == 0 )
  {
    // file header
    printJpsiFugacity << "#Jpsi fugacities" << endl;
    printJpsiFugacity << "# time  fugacities, N_ccbar, temp, energy dens. for different boxes" << endl;
    return;
  }

  printJpsiFugacity.width( 10 );
  printJpsiFugacity << time ;
  printJpsiFugacity<< "\t";


  dr = 2.0; //fm
//   dz=1.0; //fm
  deta = 0.5; // spacetime rapidty interval
  dz = time * ( exp( 2.0 * deta ) - 1.0 ) / ( exp( 2.0 * deta ) + 1.0 ); //translated to spatial coordinate z

//   cout << "deta=" << deta << "   dz=" << dz << endl;

  getJpsiFugacity( time, dr, dz, fugacity, n_jpsi, temp, enDen );
  printJpsiFugacity << fugacity;
  printJpsiFugacity << "\t";
  printJpsiFugacity << double(n_jpsi) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles();

  printJpsiFugacity << "\t";
  printJpsiFugacity << temp;
  printJpsiFugacity << "\t";
  printJpsiFugacity << enDen / theConfig->getTestparticles();


  printJpsiFugacity.width( 20 );

  dr = 5.0; //fm
//   dz=1.0; //fm
  deta = 0.5; // spacetime rapidty interval
  dz = time * ( exp( 2.0 * deta ) - 1.0 ) / ( exp( 2.0 * deta ) + 1.0 ); //translated to spatial coordinate z



  getJpsiFugacity( time, dr, dz, fugacity, n_jpsi, temp, enDen );
  printJpsiFugacity << fugacity;
  printJpsiFugacity << "\t";
  printJpsiFugacity << double(n_jpsi) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles();
  printJpsiFugacity << "\t";
  printJpsiFugacity << temp;
  printJpsiFugacity << "\t";
  printJpsiFugacity << enDen / theConfig->getTestparticles();



  printJpsiFugacity.width( 20 );

  dr = 2.0; //fm
//   dz=1.0; //fm
  deta = 1.0; // spacetime rapidty interval
  dz = time * ( exp( 2.0 * deta ) - 1.0 ) / ( exp( 2.0 * deta ) + 1.0 ); //translated to spatial coordinate z

//   cout << "deta=" << deta << "   dz=" << dz << endl;

  getJpsiFugacity( time, dr, dz, fugacity, n_jpsi, temp, enDen );
  printJpsiFugacity << fugacity;
  printJpsiFugacity << "\t";
  printJpsiFugacity << double(n_jpsi) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles();
  printJpsiFugacity << "\t";
  printJpsiFugacity << temp;
  printJpsiFugacity << "\t";
  printJpsiFugacity << enDen / theConfig->getTestparticles();



  printJpsiFugacity.width( 20 );

  dr = 5.0; //fm
//   dz=1.0; //fm
  deta = 1.0; // spacetime rapidty interval
  dz = time * ( exp( 2.0 * deta ) - 1.0 ) / ( exp( 2.0 * deta ) + 1.0 ); //translated to spatial coordinate z

//   cout << "deta=" << deta << "   dz=" << dz << endl;

  getJpsiFugacity( time, dr, dz, fugacity, n_jpsi, temp, enDen );
  printJpsiFugacity << fugacity;
  printJpsiFugacity << "\t";
  printJpsiFugacity << double(n_jpsi) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles();
  printJpsiFugacity << "\t";
  printJpsiFugacity << temp;
  printJpsiFugacity << "\t";
  printJpsiFugacity << enDen / theConfig->getTestparticles();



  printJpsiFugacity << endl;

}

void analysis::getJpsiFugacity( const double time, const double dr, const double dz, double& fugacity, int& n_jpsi, double& temp, double& enDen )
{
  double n_jpsi_equ;

  double V = M_PI * pow( dr, 2.0 ) * 2.0 * dz / pow( 0.197, 3.0 ); // 1/GeV^3

  n_jpsi = 0;

  double e_g = 0.0;
  int n_g = 0;



  for ( int i = 0; i < addedParticles.size(); i++ )
  {
    if ( addedParticles[i].T_creation <= time )
    {
      if (( addedParticles[i].Pos.Perp2() < pow( dr, 2.0 ) )  && ( fabs( addedParticles[i].Pos.Z() ) < dz ) )
      {
        // number of jpsi
        if ( addedParticles[i].FLAVOR == 50 )
        {
          n_jpsi++;
        }
      }
    }
  }
  
  for ( int i = 0; i < particles_atTimeNow.size(); i++ )
  {
    if ( particles_atTimeNow[i].Pos.T() <= time )
    {
      if (( particles_atTimeNow[i].Pos.Perp2() < pow( dr, 2.0 ) )  && ( fabs( particles_atTimeNow[i].Pos.Z() ) < dz ) )
      {
        // compute temp
        if ( particles_atTimeNow[i].FLAVOR == 0 ) // gluon
        {
          e_g += particles_atTimeNow[i].Mom.E();
          n_g++;
        }
      }
    }
  }

  temp = e_g /( 3.0 * n_g ); // GeV
  enDen = e_g /( V * pow( 0.197, 3.0 ) ); // GeV/fm^3


  
  n_jpsi_equ = theInterpolation_nJpsi.getN( temp ) * V;  // GeV^3/GeV^3 = 1


  fugacity = double(n_jpsi)/( n_jpsi_equ * theConfig->getTestparticles() * theConfig->getNaddedEvents() * theConfig->getJpsiTestparticles() );

//   cout << "t=" << time << "  temp=" << temp << "  deltaTemp=" << deltaTemp <<  "  n_jpsi_equ=" << n_jpsi_equ << "  fugacity=" << fugacity << "  n_jpsi=" << n_jpsi << "  V=" << V << endl;

}



void analysis::jpsi_correlations()
{
  double cos_delta_phi, cos_delta_theta, pt, delta_eta, eta, eta_charm;
  string filename;

  filename = filename_prefix + "_jpsi_corr_ptbins_iniJpsi";
  binning ptbins_iniJpsi(filename, 0.0, 6.0, 40);
  filename = filename_prefix + "_jpsi_corr_ptbins_secJpsi";
  binning ptbins_secJpsi(filename, 0.0, 6.0, 40);
  filename = filename_prefix + "_jpsi_corr_phibins_deltaPhiJpsiCharm";
  binning phibins_deltaPhiJpsiCharm(filename, 0.0, M_PI, 30);
  filename = filename_prefix + "_jpsi_corr_phibins_deltaPhiCharmCharm";
  binning phibins_deltaPhiCharmCharm(filename, 0.0, M_PI, 30);
  filename = filename_prefix + "_jpsi_corr_phibins_deltaThetaJpsiCharm";
  binning phibins_deltaThetaJpsiCharm(filename, 0.0, M_PI, 30);
  filename = filename_prefix + "_jpsi_corr_phibins_deltaThetaCharmCharm";
  binning phibins_deltaThetaCharmCharm(filename, 0.0, M_PI, 30);
  filename = filename_prefix + "_jpsi_corr_ptbins_charmFromJpsi";
  binning ptbins_charmFromJpsi(filename, 0.0, 6.0, 40);
  
  filename = filename_prefix + "_jpsi_corr_etabins_detaJpsiCharm";
  binning etabins_detaJpsiCharm(filename, 0.0, 20.0, 70);
  filename = filename_prefix + "_jpsi_corr_etabins_detaCharmCharm";
  binning etabins_detaCharmCharm(filename, 0.0, 20.0, 70);
  filename = filename_prefix + "_jpsi_corr_etabins_detaSecJpsiAllCharm";
  binning etabins_detaSecJpsiAllCharm(filename, 0.0, 20.0, 70);
  filename = filename_prefix + "_jpsi_corr_etabins_detaIniJpsiAllCharm";
  binning etabins_detaIniJpsiAllCharm(filename, 0.0, 20.0, 70);
  filename = filename_prefix + "_jpsi_corr_etabins_detaJpsiCharm_2d";
  binning2d etabins_detaJpsiCharm2d(filename, 0.0, 10.0, 40, 0.0, 10.0, 40);
  
  
  filename = filename_prefix + "_jpsi_corr_eta_test";
  binning eta_test(filename, -7.0, 7.0, 100);
  filename = filename_prefix + "_jpsi_corr_y_test";
  binning y_test(filename, -7.0, 7.0, 100);
  
  
  for ( int i = 0; i < addedParticles.size(); i++ )
  {
    if ( addedParticles[i].FLAVOR == 50 ) // jpsi
    {
      pt = addedParticles[i].Mom.Pt();
      eta = addedParticles[i].Mom.Pseudorapidity( addedParticles[i].m );
      
      eta_test.add(eta);
      y_test.add( addedParticles[i].Mom.Rapidity() );
      
      if( addedParticles[i].initially_produced ) // initial Jpsi
      {
        ptbins_iniJpsi.add(pt);
        
        for ( int j = 0; j < addedParticles.size(); j++ )
        {
          // all charm quarks
          if( addedParticles[j].FLAVOR == 7 || addedParticles[j].FLAVOR == 8 )
          {
            // delta eta
            eta_charm = addedParticles[j].Mom.Pseudorapidity( addedParticles[j].m );
            delta_eta = fabs( eta - eta_charm );
            etabins_detaIniJpsiAllCharm.add( delta_eta );
          }
        }
      }
      else
      {
        ptbins_secJpsi.add(pt);
        
        for ( int j = 0; j < addedParticles.size(); j++ )
        {
          // find partners of charm quarks in Jpsi
          if( ( addedParticles[j].N_EVENT_pp == addedParticles[i].N_EVENT_pp || addedParticles[j].N_EVENT_pp == addedParticles[i].N_EVENT_Cbar ) && i != j )
          {
//             if ( addedParticles[j].FLAVOR != 7 && addedParticles[j].FLAVOR != 8 )
//               cout << "error in analysis::jpsi_correlations() " << addedParticles[j].FLAVOR << "  " << addedParticles[i].FLAVOR << "  " << addedParticles[j].m << "  " << addedParticles[i].m << "  " << addedParticles[j].N_EVENT_pp << "  " << addedParticles[i].N_EVENT_pp << "  " << addedParticles[j].N_EVENT_Cbar << "  " << addedParticles[i].N_EVENT_Cbar << endl;
            
            cos_delta_phi = CosPhi( addedParticles[i].Mom, addedParticles[j].Mom );
            phibins_deltaPhiJpsiCharm.add( acos(cos_delta_phi) );
            
            cos_delta_theta = CosTheta( addedParticles[i].Mom, addedParticles[j].Mom );
            phibins_deltaThetaJpsiCharm.add( acos(cos_delta_theta) );
            
            // delta eta
            eta_charm = addedParticles[j].Mom.Pseudorapidity( addedParticles[j].m );
            delta_eta = fabs( eta - eta_charm );
            etabins_detaJpsiCharm.add( delta_eta );
            etabins_detaJpsiCharm2d.add( eta, eta_charm );
          }
          
          // all charm quarks
          if( addedParticles[j].FLAVOR == 7 || addedParticles[j].FLAVOR == 8 )
          {
            // delta eta
            eta_charm = addedParticles[j].Mom.Pseudorapidity( addedParticles[j].m );
            delta_eta = fabs( eta - eta_charm );
            etabins_detaSecJpsiAllCharm.add( delta_eta );
          }
        }  
      }
    }
    else if ( addedParticles[i].FLAVOR == 7 || addedParticles[i].FLAVOR == 8 )
    {
      if( addedParticles[i].jpsi_dissociation_number != -1 ) // produced in Jpsi dissociation
      {
        eta = addedParticles[i].Mom.Pseudorapidity( addedParticles[i].m );
        
        for ( int j = 0; j < addedParticles.size(); j++ )
        {
          if( addedParticles[j].jpsi_dissociation_number == addedParticles[i].jpsi_dissociation_number && i != j )
          {
            if ( addedParticles[j].FLAVOR != 7 && addedParticles[j].FLAVOR != 8 )
              cout << "error2 in analysis::jpsi_correlations() " << addedParticles[j].FLAVOR << "  " << addedParticles[j].m << "  " << addedParticles[i].FLAVOR << "  " << addedParticles[i].m << "  " << addedParticles[j].jpsi_dissociation_number << endl;
            
            cos_delta_phi = CosPhi( addedParticles[i].Mom, addedParticles[j].Mom );
            phibins_deltaPhiCharmCharm.add( acos(cos_delta_phi) );
      
            cos_delta_theta = CosTheta( addedParticles[i].Mom, addedParticles[j].Mom );
            phibins_deltaThetaCharmCharm.add( acos(cos_delta_theta) );
            
            // delta eta
            eta_charm = addedParticles[j].Mom.Pseudorapidity( addedParticles[j].m );
            delta_eta = fabs( eta - eta_charm );
            etabins_detaCharmCharm.add( delta_eta );
            
            ptbins_charmFromJpsi.add( addedParticles[j].Mom.Pt() );
            ptbins_charmFromJpsi.add( addedParticles[i].Mom.Pt() );
          }
        }
      }
    }
  }
  
  
  ptbins_iniJpsi.print();
  ptbins_secJpsi.print();
  phibins_deltaPhiJpsiCharm.print();
  phibins_deltaPhiCharmCharm.print();
  phibins_deltaThetaJpsiCharm.print();
  phibins_deltaThetaCharmCharm.print();
  ptbins_charmFromJpsi.print();
  etabins_detaJpsiCharm.print();
  etabins_detaCharmCharm.print();
  etabins_detaSecJpsiAllCharm.print();
  etabins_detaIniJpsiAllCharm.print();
  etabins_detaJpsiCharm2d.print();
  eta_test.print();
  y_test.print();
}



void analysis::ini_charm_correlations()
{
  double cos_delta_phi, cos_delta_theta, eta, eta_charm, delta_eta;
  string filename;


  filename = filename_prefix + "_phibins_deltaPhiIniCharm";
  binning phibins_deltaPhiIniCharm(filename, 0.0, M_PI, 30);
  filename = filename_prefix + "_phibins_deltaThetaIniCharm";
  binning phibins_deltaThetaIniCharm(filename, 0.0, M_PI, 30);
  filename = filename_prefix + "_etabins_detaIniCharm";
  binning etabins_detaIniCharm(filename, 0.0, 20.0, 70);

  for ( unsigned int i = 0; i < addedParticles.size(); i++ )
  {
    if ( addedParticles[i].FLAVOR == 7 || addedParticles[i].FLAVOR == 8 )
    {
      eta = addedParticles[i].Mom.Pseudorapidity( addedParticles[i].m );
      
      for ( unsigned int j = i+1; j < addedParticles.size(); j++ )
      {
        if( addedParticles[j].N_EVENT_pp == addedParticles[i].N_EVENT_pp &&  ( addedParticles[j].FLAVOR == 7 || addedParticles[j].FLAVOR == 8 ) )
        {
          cos_delta_phi = CosPhi( addedParticles[i].Mom, addedParticles[j].Mom );
          phibins_deltaPhiIniCharm.add( acos(cos_delta_phi) );
          
          cos_delta_theta = CosTheta( addedParticles[i].Mom, addedParticles[j].Mom );
          phibins_deltaThetaIniCharm.add( acos(cos_delta_theta) );
          
          // delta eta
          eta_charm = addedParticles[j].Mom.Pseudorapidity( addedParticles[j].m );
          delta_eta = fabs( eta - eta_charm );
          etabins_detaIniCharm.add( delta_eta );
        }
      }
    }
  }
  phibins_deltaPhiIniCharm.print();
  phibins_deltaThetaIniCharm.print();
  etabins_detaIniCharm.print();

}


void analysis::writeTempInTube( const int step  )
{
  const string sep = "  ";

  double temp, tempWithQuarks, energyDensity;
  int n_jpsi;
  double dr, dz, deta;


  double time;
  if ( step == 0 )
    time = 0.0;
  else if ( step == nTimeSteps )
    return;
  else
    time = tstep[step-1];

  if ( step == 0 )
  {
    // file header
    printTempInTube << "#temperature" << endl;
    printTempInTube << "# time  temp, tempWithQuarks, energy dens. for different boxes" << endl;
    return;
  }

  printTempInTube.width( 10 );
  printTempInTube << time ;
  printTempInTube<< "\t";


  dr = 2.0; //fm
  deta = 0.5; // spacetime rapidty interval
  dz = time * ( exp( 2.0 * deta ) - 1.0 ) / ( exp( 2.0 * deta ) + 1.0 ); //translated to spatial coordinate z

//   cout << "deta=" << deta << "   dz=" << dz << endl;

  calculateTempInTube( time, dr, dz, temp, tempWithQuarks, energyDensity  );
  printTempInTube << temp;
  printTempInTube << "\t";
  printTempInTube << tempWithQuarks;
  printTempInTube << "\t";
  printTempInTube << energyDensity;
  printTempInTube.width( 20 );
  

  dr = 5.0; //fm
  deta = 0.5; // spacetime rapidty interval
  dz = time * ( exp( 2.0 * deta ) - 1.0 ) / ( exp( 2.0 * deta ) + 1.0 ); //translated to spatial coordinate z

  calculateTempInTube( time, dr, dz, temp, tempWithQuarks, energyDensity  );
  printTempInTube << temp;
  printTempInTube << "\t";
  printTempInTube << tempWithQuarks;
  printTempInTube << "\t";
  printTempInTube << energyDensity;
  printTempInTube.width( 20 );

  dr = 2.0; //fm
  deta = 1.0; // spacetime rapidty interval
  dz = time * ( exp( 2.0 * deta ) - 1.0 ) / ( exp( 2.0 * deta ) + 1.0 ); //translated to spatial coordinate z

  calculateTempInTube( time, dr, dz, temp, tempWithQuarks, energyDensity  );
  printTempInTube << temp;
  printTempInTube << "\t";
  printTempInTube << tempWithQuarks;
  printTempInTube << "\t";
  printTempInTube << energyDensity;
  printTempInTube.width( 20 );
  printTempInTube.width( 20 );


  dr = 5.0; //fm
  deta = 1.0; // spacetime rapidty interval
  dz = time * ( exp( 2.0 * deta ) - 1.0 ) / ( exp( 2.0 * deta ) + 1.0 ); //translated to spatial coordinate z

  calculateTempInTube( time, dr, dz, temp, tempWithQuarks, energyDensity  );
  printTempInTube << temp;
  printTempInTube << "\t";
  printTempInTube << tempWithQuarks;
  printTempInTube << "\t";
  printTempInTube << energyDensity;
  printTempInTube.width( 20 );



  printTempInTube << endl;

}

// writes temperature and velocity of all specified cells in a file, used by Alex Meistrenko as input
void analysis::calculateTempInTube( const double time, const double radius, const double dz, double & temp, double & tempWithQuarks, double & energyDensity  )
{
  int cell_id;
  double pr, XT;
  
  const int minNmbTemp = 30; // minimum number of particles to calculate temperature from


  // total length in z direction
  const double zlength = dz*2.0;

  // volume
  dv = M_PI * pow( radius , 2.0 ) * zlength; // 1/GeV^3

  
  // the following routine is written for several cells -> here we do not need this, just set nCells = 1
  
  nCells = 1;

  numberInCell = new int[nCells]; // number of all particles in cell
  vx_cell = new double[nCells]; // total x-velocity of all particles in cell
  vy_cell = new double[nCells];  
  vz_cell = new double[nCells];
  vr_cell = new double[nCells];
  em_cell = new double[nCells];// total energy of all particles in cell
  prm_cell = new double[nCells];
  pzm_cell = new double[nCells];
  pr2em_cell = new double[nCells];
  pz2em_cell = new double[nCells];
  przem_cell = new double[nCells];
  densn_cell = new double[nCells];
  gama_cell = new double[nCells];
  temp_cell = new double[nCells];
  tempWithQuarks_cell = new double[nCells];


  // set all properties to 0
  for ( int i = 0; i < nCells; i++ )
  {
    numberInCell[i] = 0;
    vx_cell[i] = 0.0; 
    vy_cell[i] = 0.0;  
    vz_cell[i] = 0.0;
    vr_cell[i] = 0.0;
    em_cell[i] = 0.0;
    prm_cell[i] = 0.0;
    pzm_cell[i] = 0.0;
    pr2em_cell[i] = 0.0;
    pz2em_cell[i] = 0.0;
    przem_cell[i] = 0.0;
    densn_cell[i] = 0.0;
    gama_cell[i] = 0.0;
    temp_cell[i] = 0.0;
    tempWithQuarks_cell[i] = 0.0;
  }
  

  // sum over all particles
  for ( int i = 0; i < particles_atTimeNow.size(); i++ )
  {
    if ( FPT_COMP_LE( particles_atTimeNow[i].Pos.T(), time )  && particles_atTimeNow[i].FLAVOR < 7 ) // only gluons and light quarks
    {
      if (( particles_atTimeNow[i].Pos.Perp2() < pow( radius, 2.0 ) )  && ( fabs( particles_atTimeNow[i].Pos.Z() ) < dz ) )
      {
        // set cell id to 0
        cell_id = 0;
        
        ++numberInCell[cell_id];
        
        XT = particles_atTimeNow[i].Pos.Perp();
        if ( XT < 1.0e-5 )
        {
          pr = particles_atTimeNow[i].Mom.Pt();
        }
        else
        {
          pr = ( particles_atTimeNow[i].Mom.Px() * particles_atTimeNow[i].Pos.X()
                 + particles_atTimeNow[i].Mom.Py() * particles_atTimeNow[i].Pos.Y() ) / XT;
        }
        vr_cell[cell_id] += pr / particles_atTimeNow[i].Mom.E();
        vx_cell[cell_id] += particles_atTimeNow[i].Mom.Px() / particles_atTimeNow[i].Mom.E();
        vy_cell[cell_id] += particles_atTimeNow[i].Mom.Py() / particles_atTimeNow[i].Mom.E();
        vz_cell[cell_id] += particles_atTimeNow[i].Mom.Pz() / particles_atTimeNow[i].Mom.E();
        
        em_cell[cell_id] += particles_atTimeNow[i].Mom.E();
        prm_cell[cell_id] += pr;
        pzm_cell[cell_id] += particles_atTimeNow[i].Mom.Pz();
        pr2em_cell[cell_id] += pr * pr / particles_atTimeNow[i].Mom.E();
        pz2em_cell[cell_id] += particles_atTimeNow[i].Mom.Pz() * particles_atTimeNow[i].Mom.Pz() / particles_atTimeNow[i].Mom.E();
        przem_cell[cell_id] += pr * particles_atTimeNow[i].Mom.Pz() / particles_atTimeNow[i].Mom.E();
      }
    }
  }

  // calculate temperature for cell
  for ( int i = 0; i < nCells; i++ )
  {
    // If more than minNmbTemp particles are in the cell the temperature is calculated by them
    if(numberInCell[i] >= minNmbTemp)
    {
      calcTempCell( i );
    }
    else // take also 6 neighbor cells into account
    {
      cout << "error in calculateTempInTube(): to few particles in cell. Number=" << numberInCell[i] << "  time=" << time << endl;
      cout << particles_atTimeNow.size() << "  " << particles_atTimeNow[0].Pos.T() << "  " << particles_atTimeNow[10].Pos.T() << endl;
    }
  }
  
  
  temp = temp_cell[0];
  tempWithQuarks = tempWithQuarks_cell[0];
  energyDensity = em_cell[0]; // GeV/fm^3

  delete[] numberInCell; 
  delete[] vx_cell; 
  delete[] vy_cell;  
  delete[] vz_cell;
  delete[] vr_cell;
  delete[] em_cell;
  delete[] prm_cell;
  delete[] pzm_cell;
  delete[] pr2em_cell;
  delete[] pz2em_cell;
  delete[] przem_cell;
  delete[] densn_cell;
  delete[] gama_cell;
  delete[] temp_cell;
  delete[] tempWithQuarks_cell;
}





// writes temperature and velocity of all specified cells in a file, used by Alex Meistrenko as input
void analysis::writeTempAndVel( const int step  )
{
  int nx, ny, nz, cell_id;
  double time, pr, XT;
 
  const int minNmbTemp = 30; // minimum number of particles to calculate temperature from
//   const int minNmbTemp = 2; // minimum number of particles to calculate temperature from
  const int minNmbTempCell = 2; // minimum number of particles in one cell to calculate a Temperature. If the number is below this value, the cell is taken as empty (no temperature)
  
  binning binNumber("output/number.dat", -0.5, 49.5, 50);
  binning binTempCells("output/tempCells.dat", 0.0, 3., 70);
  binning binTempRings("output/tempRings.dat", 0.0, 2.2, 70);
  
  int count_tempCell = 0;
  int count_tempNeighborCells = 0;
  int count_noTemp = 0;
  
  int IXY = IX * IY;

  if ( step != 0 && step != nTimeSteps )
    time = tstep[step-1];
  else
    return;

  // total length of grid system
  const double xlength = 24.6;
  const double ylength = 24.6;
  const double zlength = 12.3;

  // number of cells in given direction
  // in cascade IX=40   IY=40   IZ=47
  const int nCellsx = 41;
  const int nCellsy = 41;
  const int nCellsz = 41;
  nCells = nCellsx * nCellsy * nCellsz;
  
  dv = xlength * ylength * zlength / nCells; // volume of each cell

//   int numberInCell[nCells]; // number of all particles in cell
//   double vx_cell[nCells]; // total x-velocity of all particles in cell
//   double vy_cell[nCells];  
//   double vz_cell[nCells];
//   double vr_cell[nCells];
//   double em_cell[nCells];// total energy of all particles in cell
//   double prm_cell[nCells];
//   double pzm_cell[nCells];
//   double pr2em_cell[nCells];
//   double pz2em_cell[nCells];
//   double przem_cell[nCells];
//   double densn_cell[nCells];
//   double gama_cell[nCells];
//   double temp_cell[nCells];
//   double tempWithQuarks_cell[nCells];


  numberInCell = new int[nCells]; // number of all particles in cell
  temp_numberInCell = new int[nCells]; // number of all particles in cell
  vx_cell = new double[nCells]; // total x-velocity of all particles in cell
  vy_cell = new double[nCells];  
  vz_cell = new double[nCells];
  vr_cell = new double[nCells];
  em_cell = new double[nCells];// total energy of all particles in cell
  prm_cell = new double[nCells];
  pzm_cell = new double[nCells];
  pr2em_cell = new double[nCells];
  pz2em_cell = new double[nCells];
  przem_cell = new double[nCells];
  densn_cell = new double[nCells];
  gama_cell = new double[nCells];
  temp_cell = new double[nCells];
  tempWithQuarks_cell = new double[nCells];
  
  // for copy. If too few particles in one cell, the particles from the surrounding cells are added from _org. Otherwise one would double or triple add particles from cells
  numberInCell_org = new int[nCells]; // number of all particles in cell
  vx_cell_org = new double[nCells]; // total x-velocity of all particles in cell
  vy_cell_org = new double[nCells];  
  vz_cell_org = new double[nCells];
  vr_cell_org = new double[nCells];
  em_cell_org = new double[nCells];// total energy of all particles in cell
  prm_cell_org = new double[nCells];
  pzm_cell_org = new double[nCells];
  pr2em_cell_org = new double[nCells];
  pz2em_cell_org = new double[nCells];
  przem_cell_org = new double[nCells];
// 

//   for ( int i = 0; i < nCells; i++ )
//   {
//     const int nxny = nCellsx * nCellsy;
//     const int indexZ = i / nxny;
//     const int indexY = (i -  indexZ  * nxny) / nCellsx;
//     const int indexX = i -  indexZ  * nxny - indexY * nCellsx;
//     
//     energy[i] = 0.0;
//     numberInCell[i] = indexX + nCellsx * indexY + nCellsx * nCellsy * indexZ;
//     vx[i] = double(indexX) * xlength / nCellsx;
//     vy[i] = double(indexY) * ylength / nCellsy;
//     vz[i] = double(indexZ) * zlength / nCellsz;
//   }

  // set all properties to 0
  for ( int i = 0; i < nCells; i++ )
  {
    numberInCell[i] = 0;
    temp_numberInCell[i] = 0;
    vx_cell[i] = 0.0; 
    vy_cell[i] = 0.0;  
    vz_cell[i] = 0.0;
    vr_cell[i] = 0.0;
    em_cell[i] = 0.0;
    prm_cell[i] = 0.0;
    pzm_cell[i] = 0.0;
    pr2em_cell[i] = 0.0;
    pz2em_cell[i] = 0.0;
    przem_cell[i] = 0.0;
    densn_cell[i] = 0.0;
    gama_cell[i] = 0.0;
    temp_cell[i] = 0.0;
    tempWithQuarks_cell[i] = 0.0;
    
    numberInCell_org[i] = 0;
    vx_cell_org[i] = 0.0; 
    vy_cell_org[i] = 0.0;  
    vz_cell_org[i] = 0.0;
    vr_cell_org[i] = 0.0;
    em_cell_org[i] = 0.0;
    prm_cell_org[i] = 0.0;
    pzm_cell_org[i] = 0.0;
    pr2em_cell_org[i] = 0.0;
    pz2em_cell_org[i] = 0.0;
    przem_cell_org[i] = 0.0;
  }
  
  // sum over all particles
  for ( int i = 0; i < particles_atTimeNow.size(); i++ )
  {
    if ( FPT_COMP_E( particles_atTimeNow[i].Pos.T(), time ) && particles_atTimeNow[i].FLAVOR < 7 ) // only gluons and light quarks
    {
      // determine cell id
      if ( fabs( particles_atTimeNow[i].Pos.X() - xlength / 2.0 ) < 1.0e-6 )
        nx = nCellsx - 1;
      else
        nx = int(( particles_atTimeNow[i].Pos.X() / xlength + 0.5 ) * nCellsx );

      if ( fabs( particles_atTimeNow[i].Pos.Y() - ylength / 2.0 ) < 1.0e-6 )
        ny = nCellsy - 1;
      else
        ny = int(( particles_atTimeNow[i].Pos.Y() / ylength + 0.5 ) * nCellsy );

      if ( fabs( particles_atTimeNow[i].Pos.Z() - zlength / 2.0 ) < 1.0e-6 )
        nz = nCellsz - 1;
      else
        nz = int(( particles_atTimeNow[i].Pos.Z() / zlength + 0.5 ) * nCellsz );


      if (( nx >= nCellsx ) || ( nx < 0 ) || ( ny >= nCellsy ) || ( ny < 0 ) || ( nz >= nCellsz ) || ( nz < 0 ) )
      {
        cout << "err cell_ID in temp output" << endl;
        cout << particles_atTimeNow[i].Pos.T() << "\t" << particles_atTimeNow[i].Pos.X() << "\t" << particles_atTimeNow[i].Pos.Y();
        cout << "\t" << particles_atTimeNow[i].Pos.Z() << endl;
        cout << nx << "\t" << ny << "\t" << nz << endl;
      }
      else
      {
        cell_id = nx + nCellsx * ny + nCellsx * nCellsy * nz;
        
        ++numberInCell[cell_id];
        
        XT = particles_atTimeNow[i].Pos.Perp();
        if ( XT < 1.0e-5 )
        {
          pr = particles_atTimeNow[i].Mom.Pt();
        }
        else
        {
          pr = ( particles_atTimeNow[i].Mom.Px() * particles_atTimeNow[i].Pos.X()
                 + particles_atTimeNow[i].Mom.Py() * particles_atTimeNow[i].Pos.Y() ) / XT;
        }
        vr_cell[cell_id] += pr / particles_atTimeNow[i].Mom.E();
        vx_cell[cell_id] += particles_atTimeNow[i].Mom.Px() / particles_atTimeNow[i].Mom.E();
        vy_cell[cell_id] += particles_atTimeNow[i].Mom.Py() / particles_atTimeNow[i].Mom.E();
        vz_cell[cell_id] += particles_atTimeNow[i].Mom.Pz() / particles_atTimeNow[i].Mom.E();
        
        em_cell[cell_id] += particles_atTimeNow[i].Mom.E();
        prm_cell[cell_id] += pr;
        pzm_cell[cell_id] += particles_atTimeNow[i].Mom.Pz();
        pr2em_cell[cell_id] += pr * pr / particles_atTimeNow[i].Mom.E();
        pz2em_cell[cell_id] += particles_atTimeNow[i].Mom.Pz() * particles_atTimeNow[i].Mom.Pz() / particles_atTimeNow[i].Mom.E();
        przem_cell[cell_id] += pr * particles_atTimeNow[i].Mom.Pz() / particles_atTimeNow[i].Mom.E();
        
        // temp of particles summed
        tempWithQuarks_cell[cell_id] += particles_atTimeNow[i].temperature;
        if( particles_atTimeNow[i].temperature >= 0.1)
          ++temp_numberInCell[cell_id];
      }
    }
  }
  
  // duplicate properties of cell
  for ( int i = 0; i < nCells; i++ )
  {
    numberInCell_org[i] = numberInCell[i];
    vr_cell_org[i] = vr_cell[i];
    vx_cell_org[i] = vx_cell[i];
    vy_cell_org[i] = vy_cell[i];
    vz_cell_org[i] = vz_cell[i];
    
    em_cell_org[i] = em_cell[i];
    prm_cell_org[i] = prm_cell[i];
    pzm_cell_org[i] = pzm_cell[i];
    pr2em_cell_org[i] = pr2em_cell[i];
    pz2em_cell_org[i] = pz2em_cell[i];
    przem_cell_org[i] = przem_cell[i];
  }

  
  // determine cells which do not have enough particles
  for ( int i = 0; i < nCells; i++ )
  {
    // If less than minNmbTempCell particles are in the cell it is taken as empty
    if(numberInCell[i] >= minNmbTempCell)
    {
      // If less than minNmbTemp particles are in the cell the temperature the 6 neighbor cells are also taken into account
      if(numberInCell[i] < minNmbTemp)
      {
        int cell_id_neighbor;
        
        //neighbors in x direction
        cell_id_neighbor = i-1;
        addNeighborCells( i, cell_id_neighbor );
        cell_id_neighbor = i+1;
        addNeighborCells( i, cell_id_neighbor );
        //neighbors in y direction
        cell_id_neighbor = i+IX;
        addNeighborCells( i, cell_id_neighbor );
        cell_id_neighbor = i-IX;
        addNeighborCells( i, cell_id_neighbor );
        //neighbors in z direction
        cell_id_neighbor = i-IXY;
        addNeighborCells( i, cell_id_neighbor );
        cell_id_neighbor = i+IXY;
        addNeighborCells( i, cell_id_neighbor );
      }
    }
  }
  
  // calculate temperature for each cell
  for ( int i = 0; i < nCells; i++ )
  {
    // If less than minNmbTempCell particles are in the cell it is taken as empty
    if(numberInCell[i] >= minNmbTempCell)
    {
      // If more than minNmbTemp particles are in the cell the temperature is calculated by them
      if(numberInCell[i] >= minNmbTemp)
      {
        calcTempCell( i );
        
        tempWithQuarks_cell[i] = tempWithQuarks_cell[i]/temp_numberInCell[i];

        binTempCells.add(temp_cell[i]);
        count_tempCell++;
      }
      else // take also 6 neighbor cells into account
      {
        count_noTemp++;
        
        temp_cell[i] = 0.0; // not enough particles in surrounding to calculate temperature
        
        vx_cell[i] = 0.0;
        vy_cell[i] = 0.0;
        vz_cell[i] = 0.0;
      }
    }
    else
    {
      if(numberInCell[i] == 0)
      {
        temp_cell[i] = -2.0; // completely empty
        vx_cell[i] = 0.0;
        vy_cell[i] = 0.0;
        vz_cell[i] = 0.0;
      }
      else
      {
        temp_cell[i] = -1.0; // very few particles
        vx_cell[i] = 0.0;
        vy_cell[i] = 0.0;
        vz_cell[i] = 0.0;
      }
    }
  }


  string filename,name;
  stringstream ss;
  ss << time*10;
  filename = filename_prefix + "_tempVel_" + ss.str() + ".dat";

  fstream file( filename.c_str(), ios::out | ios::trunc  );

  if ( step == 1 )
  {
    file << "#temperature and velocity" << endl;
    file << "#simulation parameter:" << endl;
    file << "#time = " << time << endl;
    file << "#testparticles= " << theConfig->getTestparticles() << endl;
    file << "#runtime= " << theConfig->getRuntime() << endl;
    file << "#sqrtS= " << theConfig->getSqrtS() << " GeV" << endl;
    file << "#b= " << theConfig->getImpactParameter() << " fm" << endl;
    file << "#(" << theConfig->getA()  << ") on ("
    << theConfig->getA()  << ")" << endl;
    file << "# T in GeV/fm^3" << endl;
    file << "#" << endl;
    file << "#" << "T" << "\t" << "vx" << "\t" << "vy" << "\t" << "vz" << endl;
  }
  
  for ( int i = 0; i < nCells; i++ )
  {
    file << temp_cell[i] << "\t";
    file << tempWithQuarks_cell[i] << "\t";
    file << vx_cell[i] << "\t";
    file << vy_cell[i] << "\t";
    file << vz_cell[i] << endl;
  }
  
  filename = filename + "_spatial";
  
  fstream file_spatial( filename.c_str(), ios::out | ios::trunc  );
  
  for ( int i = 0; i < nCells; i++ )
  {
    const int nxny = nCellsx * nCellsy;
    const int indexZ = i / nxny;
    const int indexY = (i -  indexZ  * nxny) / nCellsx;
    const int indexX = i -  indexZ  * nxny - indexY * nCellsx;
    
    if(indexY == int(nCellsy/2))
    {
//       file_spatial << i << "\t";
      file_spatial << double(indexZ) * zlength / nCellsz << "\t";
      file_spatial << double(indexX) * xlength / nCellsx << "\t";
      
      file_spatial << temp_cell[i] << "\t";
      file_spatial << tempWithQuarks_cell[i] << "\t";
      file_spatial << vx_cell[i] << "\t";
      file_spatial << vy_cell[i] << "\t";
      file_spatial << vz_cell[i] << endl;
    }
  }

  delete[] numberInCell; 
  delete[] vx_cell; 
  delete[] vy_cell;  
  delete[] vz_cell;
  delete[] vr_cell;
  delete[] em_cell;
  delete[] prm_cell;
  delete[] pzm_cell;
  delete[] pr2em_cell;
  delete[] pz2em_cell;
  delete[] przem_cell;
  delete[] densn_cell;
  delete[] gama_cell;
  delete[] temp_cell;
  delete[] tempWithQuarks_cell;
}

void analysis::calcTempCell( const int cell_id ) 
{
  gama_cell[cell_id] = 1.0;
  vx_cell[cell_id] = vx_cell[cell_id] / numberInCell[cell_id];
  vy_cell[cell_id] = vy_cell[cell_id] / numberInCell[cell_id];
  vz_cell[cell_id] = vz_cell[cell_id] / numberInCell[cell_id];
  vr_cell[cell_id] = vr_cell[cell_id] / numberInCell[cell_id];
  if ( vr_cell[cell_id] > 0.0 ) gama_cell[cell_id] = 1.0 / sqrt( 1.0 - vz_cell[cell_id] * vz_cell[cell_id] - vr_cell[cell_id] * vr_cell[cell_id] );
  else
  {
    vr_cell[cell_id] = 0.0;
    gama_cell[cell_id] = 1.0 / sqrt( 1.0 - vz_cell[cell_id] * vz_cell[cell_id] );
  }
  
  densn_cell[cell_id] = double(numberInCell[cell_id]) / dv / theConfig->getTestparticles() / gama_cell[cell_id];//1/fm^3
  em_cell[cell_id] = ( em_cell[cell_id] - 2.0 * vr_cell[cell_id] * prm_cell[cell_id] - 2.0 * vz_cell[cell_id] * pzm_cell[cell_id] + vr_cell[cell_id] * vr_cell[cell_id] * pr2em_cell[cell_id]
  + vz_cell[cell_id] * vz_cell[cell_id] * pz2em_cell[cell_id] + 2.0 * vr_cell[cell_id] * vz_cell[cell_id] * przem_cell[cell_id] )
  / theConfig->getTestparticles() / dv * gama_cell[cell_id] * gama_cell[cell_id];//GeV/fm^3
  
  temp_cell[cell_id] = em_cell[cell_id] / ( 3.0 * densn_cell[cell_id] );
  // assume thermal equilibrium and additional quark flavor
//   int Nflavor_temp = 3;
//   tempWithQuarks_cell[cell_id] = pow(pi*pi/3.0 / (16.+12.*Nflavor_temp) * em_cell[cell_id]  * pow(0.197,3.0) , 1.0/4.0);
}

void analysis::addNeighborCells( const int cell_id, const int neighborCell_id ) 
{
  if(neighborCell_id >= 0 && neighborCell_id < nCells)
  {
    numberInCell[cell_id] += numberInCell_org[neighborCell_id];
    vr_cell[cell_id] += vr_cell_org[neighborCell_id];
    vx_cell[cell_id] += vx_cell_org[neighborCell_id];
    vy_cell[cell_id] += vy_cell_org[neighborCell_id];
    vz_cell[cell_id] += vz_cell_org[neighborCell_id];
    
    em_cell[cell_id] += em_cell_org[neighborCell_id];
    prm_cell[cell_id] += prm_cell_org[neighborCell_id];
    pzm_cell[cell_id] += pzm_cell_org[neighborCell_id];
    pr2em_cell[cell_id] += pr2em_cell_org[neighborCell_id];
    pz2em_cell[cell_id] += pz2em_cell_org[neighborCell_id];
    przem_cell[cell_id] += przem_cell_org[neighborCell_id];
  }
}






// print dndy, detdy and mean pt of gluons
void analysis::print_dndy(const string subfix )
{
  double y, pt;
  double pt_sum = 0.0;
  
  string filename;

  filename = filename_prefix + "_dndy_" + subfix;
  binning dndy(filename, -10.0, 10.0, 200);
  filename = filename_prefix + "_dedy_" + subfix;
  binningValues dedy(filename, -10.0, 10.0, 200);
  
  filename = filename_prefix + "_dndpt_" + subfix;
  binning dndpt(filename, 0.0, 20.0, 200);
  
  for ( unsigned int i = 0; i < particles_atTimeNow.size(); i++ )
  {
    pt = particles_atTimeNow[i].Mom.Pt();
    y  = particles_atTimeNow[i].Mom.Rapidity();
    
//     if(y < -8.0 || isnan(y))
// //       cout << y << "\t" << particles_atTimeNow[i].Mom.E() << "\t" << particles_atTimeNow[i].Mom.Pz()  << "\t" << pt << "\t" << i << endl;
//       cout << i << "\t";
    
    dndy.add(y);
    dedy.add(y,pt);
    dndpt.add(pt);
    pt_sum += pt;
  }
  
  cout << "mean pt = " << pt_sum/particles_atTimeNow.size() << endl;
  
  dndy.print();
  dedy.print();
  dndpt.print();

}











void analysis::analyseAngleDe()
{
  double costheta; // angle between D meson and e-
  double cosphi; // transverse angle between D meson and e-
  double ept;
  int k_e;
  
  string filename;
  
  
  filename = filename_prefix + "_dmeson_electron_cosTheta";
  binningValues binsCosTheta(filename, 0.0, 10.0, 200);
  filename = filename_prefix + "_dmeson_electron_cosPhiTrans";
  binningValues binsCosPhi(filename, 0.0, 10.0, 200);
  
  for ( unsigned int i = 0; i < addedParticles.size(); i++ )
  {
    if ( addedParticles[i].FLAVOR == 7 || addedParticles[i].FLAVOR == 8 )
    {
      // there are several electrons per dmeson/charm quark
      for(int k = 0; k < theConfig->getNumberElectronStat(); k++)
      {      
        k_e = ( i ) * theConfig->getNumberElectronStat() + k ;
        costheta = CosTheta( addedParticlesCopy[i].Mom, addedPartcl_electron[k_e].Mom );
        
        cosphi = CosPhi( addedParticlesCopy[i].Mom, addedPartcl_electron[k_e].Mom );

        ept = addedPartcl_electron[k_e].Mom.Pt(); // e- pt
        
        binsCosTheta.add(ept, costheta);
        binsCosPhi.add(ept, cosphi);
      }
    }
  }
  
  binsCosTheta.print();
  binsCosPhi.print();
  
}


// void analysis::analyseCharmTestJetEvolution(const int step)
// {
//   double energy_sum = 0.0;
//   int count = 0;
//   for ( int i = 0; i < addedParticles.size(); i++ )
//   {
//     if ( addedParticles[i].FLAVOR == 7 )
//     {
//       count++;
//       energy_sum += addedParticles[i].E;
//     }
//   }
// 
//   charmJetEnergy[step] = energy_sum/count;
//   timestepAnalysed[step] = true;
// }
// 
// void analysis::analyseCharmTestJet()
// {
//   string filename;
//   
//   filename = filename_prefix + "_Ebins";
//   binning Ebins(filename, 0.0, 12.0, 24);
//   filename = filename_prefix + "_xbins";
//   binning xbins(filename, -3.0, 5.0, 70);
//   filename = filename_prefix + "_ybins";
//   binning ybins(filename, -3.0, 3.0, 70);
//   filename = filename_prefix + "_elossbins";
//   binning elossbins(filename, -2.0, 10.0, 70);
//   
//   for ( int i = 0; i < addedParticles.size(); i++ )
//   {
//     if ( addedParticles[i].FLAVOR == 7 )
//     {
//       Ebins.add(addedParticles[i].E);
//       xbins.add(addedParticles[i].X);
//       ybins.add(addedParticles[i].Y);
//       elossbins.add(10.0 - addedParticles[i].E);
//     }
//   }
//   
//   Ebins.print();
//   xbins.print();
//   ybins.print();
//   elossbins.print();
//   
//   filename = filename_prefix + "_jetEvolution";
//   fstream print( filename.c_str(), ios::out | ios::trunc );
// 
//   print << "# time   #mean jet energy" << endl;
//   for ( int i = 0; i < nTimeSteps + 1; i++ )
//   {
//     if ( timestepAnalysed[i] )
//     {
//       print.width( 15 );
//       if ( i == 0 )
//         print << "0";
//       else
//         print << tstep[i-1];
//       print.width( 15 );
//       print << charmJetEnergy[i];
//       print << endl;
//     }
//   }
// 
// }



void analysis::jpsiEvolution( int step )
{
  double pt_min = 0;
  
  if( theConfig->getOutputScheme() == cms_jpsi )
    pt_min = 6.5;
  
  
  int countJpsi_all = 0;
  int countJpsi_ini = 0;
//   int countJpsi_midNormRap = 0;
  int countJpsi_midPseudoRap_all = 0;
  int countJpsi_midPseudoRap_ini = 0;
  int countJpsi_forwardPseudoRap_all = 0;
  int countJpsi_forwardPseudoRap_ini = 0;
  int countJpsi_midNormRap_all = 0;
  int countJpsi_midNormRap_ini = 0;
  int countJpsi_forwardNormRap_all = 0;
  int countJpsi_forwardNormRap_ini = 0;
//   int countJpsi_midSpaceTimeRap = 0;

  for ( unsigned int i = 0; i < addedParticles.size(); i++ )
  {
    if ( addedParticles[i].FLAVOR == 50 )
    {
      double pt = addedParticles[i].Mom.Pt();
      
      countJpsi_all++;
        
      if(addedParticles[i].initially_produced)
        countJpsi_ini++;
      
      if( pt >= pt_min )
      {
        double pseudorap = addedParticles[i].Mom.Pseudorapidity( addedParticles[i].m );
        
        if ( fabs( pseudorap ) >= rapidityRanges[0].yleft && fabs( pseudorap ) <= rapidityRanges[0].yright )
        {
          countJpsi_midPseudoRap_all++;
          if(addedParticles[i].initially_produced)
            countJpsi_midPseudoRap_ini++;
        }
        
        if ( fabs( pseudorap ) >= rapidityRanges[1].yleft && fabs( pseudorap ) <= rapidityRanges[1].yright )
        {
          countJpsi_forwardPseudoRap_all++;
          if(addedParticles[i].initially_produced)
            countJpsi_forwardPseudoRap_ini++;
        }
        
              
        // normal rapidity
        double normrap = addedParticles[i].Mom.Rapidity();
        
        if ( fabs( normrap ) >= rapidityRanges[0].yleft && fabs( normrap ) <= rapidityRanges[0].yright )
        {
          countJpsi_midNormRap_all++;
          if(addedParticles[i].initially_produced)
            countJpsi_midNormRap_ini++;
        }
        
        
        if ( fabs( normrap ) >= rapidityRanges[1].yleft && fabs( normrap ) <= rapidityRanges[1].yright )
        {
          countJpsi_forwardNormRap_all++;
          if(addedParticles[i].initially_produced)
            countJpsi_forwardNormRap_ini++;
        }
          
        
  //       // space time rapidity
  //       double strap = addedParticles[i].Pos.Rapidity();
  //       if( fabs(strap) < midrap )
  //         countJpsi_midSpaceTimeRap++;
      }
    }
  }


  numberJpsi_all_time[step] = countJpsi_all;
  numberJpsi_ini_time[step] = countJpsi_ini;
  numberJpsi_midPseudoRap_all_time[step] = countJpsi_midPseudoRap_all;
  numberJpsi_midPseudoRap_ini_time[step] = countJpsi_midPseudoRap_ini;
  numberJpsi_forwardPseudoRap_all_time[step] = countJpsi_forwardPseudoRap_all;
  numberJpsi_forwardPseudoRap_ini_time[step] = countJpsi_forwardPseudoRap_ini;
  numberJpsi_midNormRap_all_time[step] = countJpsi_midNormRap_all;
  numberJpsi_midNormRap_ini_time[step] = countJpsi_midNormRap_ini;
  numberJpsi_forwardNormRap_all_time[step] = countJpsi_forwardNormRap_all;
  numberJpsi_forwardNormRap_ini_time[step] = countJpsi_forwardNormRap_ini;
//   numberJpsi_midSpaceTimeRap_time[step] = countJpsi_midSpaceTimeRap;
  numberJpsiProd_time[step] = ns_heavy_quarks::jpsicreation;
  numberJpsiDiss_time[step] = ns_heavy_quarks::jpsi_dissociation;
  numberJpsiDissTd_time[step] = ns_heavy_quarks::jpsi_dissociation_from_temperature;
  numberCCbGG_time[step] = ns_heavy_quarks::charmAnnihil;
  timestepAnalysed[step] = true;
}

void analysis::printJpsiEvolution()
{
  string filename = filename_prefix + "_jpsiEvolution";
  fstream print_je( filename.c_str(), ios::out | ios::trunc );
  
  const double delta_y_mid = 2.0 * ( rapidityRanges[0].yright - rapidityRanges[0].yleft );
  const double delta_y_forward = 2.0 * ( rapidityRanges[1].yright - rapidityRanges[1].yleft );

  print_je << "# charm Annihaltion and J/psi processes" << endl;
  print_je << "# time   #J/psis    #J/psi production     # ccb->gg processes" << endl;
  for ( int i = 0; i < nTimeSteps + 1; i++ )
  {
    if ( timestepAnalysed[i] )
    {
      print_je.width( 15 );
      if ( i == 0 )
        print_je << "0";
      else
        print_je << tstep[i-1];
      print_je.width( 15 );
      print_je << double( numberJpsi_all_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles();
      print_je.width( 15 );
      print_je << double( numberJpsi_ini_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles();
      print_je.width( 15 );
      print_je << double( numberJpsiProd_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles();
      print_je.width( 15 );
      print_je << double( numberJpsiDiss_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles();
      print_je.width( 15 );
      print_je << double( numberJpsiDissTd_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles();
      print_je.width( 15 );
      print_je << double( numberCCbGG_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles();

      
      //  / (2.0*1.0) because the rapidity range is 1.0, but for + and . Therefore, times 2.
      print_je.width( 15 );
      print_je << double( numberJpsi_midPseudoRap_all_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles() / delta_y_mid;
      print_je.width( 15 );
      print_je << double( numberJpsi_midPseudoRap_ini_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles() / delta_y_mid;
      print_je.width( 15 );
      print_je << double( numberJpsi_forwardPseudoRap_all_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles() / delta_y_forward;
      print_je.width( 15 );
      print_je << double( numberJpsi_forwardPseudoRap_ini_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles() / delta_y_forward;
      
      print_je.width( 15 );
      print_je << double( numberJpsi_midNormRap_all_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles() / delta_y_mid;
      print_je.width( 15 );
      print_je << double( numberJpsi_midNormRap_ini_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles() / delta_y_mid;
      print_je.width( 15 );
      print_je << double( numberJpsi_forwardNormRap_all_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles() / delta_y_forward;
      print_je.width( 15 );
      print_je << double( numberJpsi_forwardNormRap_ini_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles() / delta_y_forward;
      
//       print_je.width( 15 );
//       print_je << double( numberJpsi_midSpaceTimeRap_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles();
      print_je << endl;
    }
  }
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
  string filename;
  if( studyJets )
    filename = filename_prefix + ".f4";
  else
    filename = "/dev/null";
  fstream file( filename.c_str(), ios::out | ios::trunc );
  
  //---- print header if file is empty ----
  file.seekp( 0, ios::end );
  long size = file.tellp();
  file.seekp( 0, ios::beg );
  if ( size == 0 )
    printHeader( file, jets, end );
  //---------------------------------------
  
  for ( unsigned int jet = 0; jet < jetTracker.size(); jet++ )
  {
    for ( unsigned int event = 0; event < jetTracker[jet].size(); event++ )
    {
      file << jetTracker[jet][event].jet_ID_in << sep << jetTracker[jet][event].jet_ID_out << sep
      << jetTracker[jet][event].coll_type << sep;
      for ( int i = 0; i < 4; i++ )
      {
        file << jetTracker[jet][event].P_proj_in(i) << sep;
      }
      for ( int i = 0; i < 4; i++ )
      {
        file << jetTracker[jet][event].P_proj_out(i) << sep;
      }
      for ( int i = 0; i < 4; i++ )
      {
        file << jetTracker[jet][event].R_proj(i) << sep;
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
  for ( unsigned int i = 1; i < rings.size(); i++ )
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
  for ( unsigned int i = 1; i < rings.size(); i++ )
  {
    rings[i].relocate( rings[i-1].maxRadius, rings[i-1].maxRadius + deltaR );
  }
}



int analysisRingStructure::getIndex( const double _xt ) const
{
  unsigned int index = getIndexPure( _xt );
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
  return getIndex( _particle.Pos.Perp() );
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
