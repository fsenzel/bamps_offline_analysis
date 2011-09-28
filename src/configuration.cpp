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
#include <fstream>
#include <math.h>
#include <vector>

#include "configuration.h"
#include "analysis.h"

using std::cout;
using std::endl;
using std::fstream;
using std::ios;
using std::ifstream;


/** @brief definition of particles vector, defined extern in configuration.h */
std::vector<Particle> ns_casc::particles;
std::vector<Particle> ns_casc::particles_init;
std::vector<Particle> ns_casc::particles_atTimeNow;
std::vector<Particle> ns_casc::particles_atTimeNowCopy;
std::vector<Particle> ns_casc::addedParticles;
std::vector<Particle> ns_casc::addedParticlesCopy;


// global variables which are changed in config.cpp (input parameter)
double Mcharm = 1.5;       // GeV, mass of charm quark
double Mbottom = 4.8;      // GeV, mass of bottom quark
double Kggccb = 1.0;       // K factor for gg -> ccb cross section
double KgQgQ = 1.0;        // K factor for gQ -> gQ cross section
double kappa_gQgQ = 1.0;   // kappa factor for debye mass, usually 0.2 (Peshier,Gossiaux,Aichelin)
bool isCouplingRunningForGQ = false; // true if a running coupling is used for gQ -> gQ 

double dt = 0.1;

bool isIsotropicCrossSec = false;
bool isConstantCrossSec = false;
double constantCrossSec = 10.0; // mb

// interpolation_gQ theIgc;
// interpolation_gQ theIgb;


/**
* The only constructor provided for this class. It needs to be passed the command line arguments,
* from which the name of a possible input file is taken.
*
* Defaults for all parameters are set in case there is no input file or it cannot be read.
*
* @param[in] argc number of command line arguments, passed from the calling process
* @param[in] argv[] command line arguments, passed from the calling process
*/
config::config(const int argc, const char * const argv[])
{
  //--set default values---------------------------------------
  name = "default";
  runtime = 5.0; 
  A = 197.0; 
  Aatomic = 79.0;
  B = 197.0; 
  Batomic = 79.0;
  sqrtS = 200.0;
  impact = 0.0;
  P0 = 1.4;
  testparticles = 35;
  seed = 0;
  jetMfpComputationSwitch = computeMfpDefault;
  interpolationBorder = 50;
  outputForOfflineReconstruction = true;
  detailedParticleOutput = false;
  movieOutput = false;
  freezeOutEnergyDensity = 1.0;
  pathdirOfflineData = "offline_data";
  pythiaParticleFile = "-";
  cgcParticleFile = "-";
  initialStateType = miniJetsInitialState;
  localCluster = true;
  movieOutputBackground = false;
  movieOutputJets = false;
  factor_dt = 0.8;
  //these values are used in case no input file is specified 
  //or in case the input file lacks certain statements
  //-----------------------------------------------------------
  
  switch (argc)
  {
    case 1: 
      cout << "WARNING: Argument specifying input file is missing. Defaults will be used." << endl;
      break;
    case 2:
      inputFile = argv[1];
      readInputfile();
      cout << "reading input file: " << inputFile << endl;
      break;
    default:
      cout << "WARNING: Too many arguments passed. Defaults will be used." << endl;
      break;
  }
  
  cout << "Cascade data folder: " << pathdirOfflineData << endl;
  
  if(dt_specified)
    cout << "dt = " << dt << endl;
  else
    cout << "Get dt for each step from cascade data." << endl;
  
//   if(isConstantCrossSec)
//     cout << "Constant cross section of " << constantCrossSec << " mb for gQ-> gQ." << endl;
//   else
//   {
//     cout << "pQCD cross section for gQ-> gQ with K = " << KgQgQ;
//     if(isCouplingRunningForGQ)
//       cout << " and running coupling";
//     else
//       cout << " and constant coupling";
//     cout << ", kappa_gQ = " << kappa_gQgQ << endl;
//   }
// 
//   if(isIsotropicCrossSec)
//     cout << "Sample angle isotropic." << endl;
// 
//   theIgc.configure(charm, isCouplingRunningForGQ);
//   theIgb.configure(bottom, isCouplingRunningForGQ);
  
  writeInput();
}


/** Returns the length of the name associated with the current job.
* @return name.length()
*/ 
int config::getNameLength() const
{
  return name.length();
}


/**
* Read the input file and obtain the parameters from it. Each parameter in the inputfile (plain text) must be on
* a single line and be specified as follows
*     #<parameter name> value
*
* Where <parameter name> can be: name, time, A(A), Z(A), A(B), Z(B), sqrtS, impact, P0, testparticles, seed, iterMFP, iterPT, outLevel
* The file must be ended by
*     #end
* In case a parameter is not qualified in the input file a default value is specified in the constructor.
* If a parameter is defined multiple times in the input file, the last occurence should be significant. However such experiments
* are not encouraged.
*/
void config::readInputfile()
{  
  extern double Mcharm;       // GeV, mass of charm quark
  extern double Mbottom;       // GeV, mass of charm quark
  extern double Kggccb;       // K factor for gg -> ccb cross section
  extern double KgQgQ;        // K factor for gc -> gc cross section
  extern double kappa_gQgQ;   // kappa factor for debye mass, usually 0.2 (Peshier,Gossiaux,Aichelin)
  extern bool isCouplingRunningForGQ; // true if a running coupling is used for gQ -> gQ 
  extern double dt;
  
  fstream file(inputFile.c_str(), ios::in);
  string content = "";
  
  int tempInitialState = -1;
  int tempJetMfpComputationType = 0;
  
  do
  {
    file >> content;
    
    if (content == "#name")
      file >> name;
    else if (content == "#time")
      file >> runtime;
    else if (content == "#A(A)")
      file >> A;
    else if (content == "#Z(A)")
      file >> Aatomic;
    else if (content == "#A(B)")
      file >> B;
    else if (content == "#Z(B)")
      file >> Batomic;
    else if (content == "#sqrtS")
      file >> sqrtS;
    else if (content == "#impact")
      file >> impact;
    else if (content == "#P0")
      file >> P0;
    else if (content == "#testparticles")
      file >> testparticles;
    else if (content == "#seed")	
      file >> seed;
    else if (content == "#jetMfpComputation")	
      file >> tempJetMfpComputationType;
    else if (content == "#iterpolationBorder")
      file >> interpolationBorder;
    else if (content == "#offlineOutput")
      file >> outputForOfflineReconstruction;
    else if (content == "#pythiaParticleFile")
      file >> pythiaParticleFile;
    else if (content == "#pathdirOfflineData")
      file >> pathdirOfflineData;  
    else if (content == "#cgcParticleFile")
      file >> cgcParticleFile;
    else if (content == "#eCrit")
      file >> freezeOutEnergyDensity;
    else if (content == "#OfflineDataPath")
      file >> pathdirOfflineData;
    else if (content == "#initialState")  // 0 = mini jets;  1 = PYTHIA;  2 = CGC (not implemented yet)
      file >> tempInitialState;
    else if ( content == "#jetMovie")
    {
      int temp;
      file >> temp;
      if ( temp != 0 )
      {
        movieOutputJets = true;
      }
    }
    else if ( content == "#backgroundMovie")
    {
      int temp;
      file >> temp;
      if ( temp != 0 )
      {
        movieOutputBackground = true;
      }
    }
    else if (content == "#localCluster")
    {
      file >> localCluster;
    }
    else if (content == "#dt")
    {
      file >> delta_t;
      dt = delta_t;
      dt_specified = true;
    }
    else if (content == "#isotropicCrossSec")
      file >> isIsotropicCrossSec;
    else if (content == "#constantCrossSec")
    {
      file >> constantCrossSec;
      isConstantCrossSec = true;
    }
    else if ( content == "#nAdded" )
    {
      file >> numberOfParticlesToAdd;
    }
    else if ( content == "#factor_dt" )
    {
      file >> factor_dt;
    }
    else if ( content == "#minPT_added" )
    {
      file >> minimumPT;
    }
    //       else if (content == "#pythiaParticleFileCharm")
      //         file >> pythiaParticleFileCharm;
      //       else if (content == "#pythiaParticleFileCharmTestparticles")
      //         file >> pythiaParticleFileCharmTestparticles;
      
      //       else if (content == "#Mcharm")
      //         file >> Mcharm;
      //       else if (content == "#Mbottom")
      //         file >> Mbottom;
      //       else if (content == "#Kggccb")
      //         file >> Kggccb;
      //       else if (content == "#Kgcgc")
      //         file >> KgQgQ;
      //       else if (content == "#KIniCharm")
      //         file >> KIniCharm;   
      
      //       else if (content == "#kappa_gQgQ")
      //         file >> kappa_gQgQ;  
      //       else if (content == "#Kgcgc")
      //         file >> kappa_gQgQ;  
      //       else if (content == "#asRunGQ")
      //         file >> isCouplingRunningForGQ;   
  }
  while(content != "#end");
  
  if ( tempJetMfpComputationType < 0 || tempJetMfpComputationType > 2 )
  {
    //     string errMsg = "#jetMfpComputation should be specified as 0, 1 or 2";
    //     throw eConfiguration_error( errMsg );
    cout << "  #jetMfpComputation should be specified as 0, 1 or 2" << endl;
    cout << "  using default jetMfpComputation = 0" << endl;
    tempJetMfpComputationType = 0;
  }
  
  if ( tempInitialState >= 0 && tempInitialState < 3 )
  {
    initialStateType = static_cast<INITIAL_STATE_TYPE>( tempInitialState );
  }
  
}

/**
* write out the values of the input parameters to cout
*/
void config::writeInput()
{
  cout << "#name " << name << endl;
  cout << "#time " << runtime << endl;
  cout << "#A(A) " << A << endl;
  cout << "#Z(A) " << Aatomic << endl;
  cout << "#A(B) " << B << endl;
  cout << "#Z(B) " << Batomic << endl;
  cout << "#sqrtS " << sqrtS << endl;
  cout << "#impact " << impact << endl;
  cout << "#P0 " << P0 << endl;
  cout << "#testparticles " << testparticles << endl;
  cout << "#seed " << seed << endl;
  cout << "#iterpolationBorder " << interpolationBorder << endl;
  cout << "#offlineOutput " << outputForOfflineReconstruction << endl;
  cout << "#pythiaParticleFile " << pythiaParticleFile << endl;
  cout << "#cgcParticleFile " << cgcParticleFile << endl;
  cout << "#eCrit " << freezeOutEnergyDensity << endl;
  cout << "#OfflineDataPath " << pathdirOfflineData << endl;
  cout << "#initialState " << initialStateType << endl;
  cout << "#end"<< endl;
  
}
