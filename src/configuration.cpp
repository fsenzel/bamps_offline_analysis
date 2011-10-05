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
#include <sstream>
#include <string>

#include "particle.h"
#include "configuration.h"
#include "analysis.h"
#include "parameter.h"

using std::cout;
using std::endl;
using std::fstream;
using std::ios;
using std::ifstream;

/** @brief definition of particles vector, defined extern in configuration.h */
std::vector<ParticleOffline> ns_casc::particles;
std::vector<ParticleOffline> ns_casc::particles_init;
std::vector<ParticleOffline> ns_casc::particles_atTimeNow;
std::vector<ParticleOffline> ns_casc::particles_atTimeNowCopy;
std::vector<ParticleOffline> ns_casc::addedParticles;
std::vector<ParticleOffline> ns_casc::addedParticlesCopy;

using namespace ns_casc;

double dt = 0.1;

int MaxN, number, numberAdded, IX, IY, IZ;
int cellcut;
double deta;

initParam iparam;


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



void config::readAndPrepareInitialSettings( offlineOutputInterface* const offlineInterface )
{
  boost::shared_ptr< offlineDataSimulationParameters > ptrSimulationData = offlineInterface->readOfflineDataFromArchive< offlineDataSimulationParameters >();
  originalSeed = ptrSimulationData->seed;
  sqrtS = ptrSimulationData->sqrtS;
  impact = ptrSimulationData->impactParameter;
  A = ptrSimulationData->massNumberNucleusA;
  Aatomic = ptrSimulationData->atomicNumberNucleusA;
  B = ptrSimulationData->massNumberNucleusB;
  Batomic = ptrSimulationData->atomicNumberNucleusB;
  testparticles = ptrSimulationData->numberOfTestparticles;
  N_init = ptrSimulationData->initialNumberOfParticles;
  timefirst = ptrSimulationData->firstTimeStep;
  timeshift = ptrSimulationData->timeShift;
  freezeOutEnergyDensity = ptrSimulationData->freezeOutEnergyDensity;
  ringNumber = ptrSimulationData->ringStructureSize;
  centralRingRadius = ptrSimulationData->ringStructureCentralRadius;
  deltaR = ptrSimulationData->ringStructureDeltaR;
  dx = ptrSimulationData->cellSizeDeltaX;
  dy = ptrSimulationData->cellSizeDeltaY;
  transLen = ptrSimulationData->transversalSize;
  IX = ptrSimulationData->gridSizeX;
  IY = ptrSimulationData->gridSizeY;
  IZ = ptrSimulationData->gridSizeZ;
  
  numberOfParticlesToAdd *= testparticles;
  numberAdded = numberOfParticlesToAdd;
  cout << "** add " << numberOfParticlesToAdd << " particles with pt_min = " << minimumPT << endl;
  cout << timefirst << "  " << N_init << endl;

  cellcut = 4;
  //-------------------------//

  number = N_init;
  
  particles.reserve( N_init * 1.8 );
  particles_init.reserve( N_init );
  particles.resize( N_init * 1.8 );
  particles_init.resize( N_init );

  boost::shared_ptr< offlineDataInitialParticles > ptrInitialParticles = offlineInterface->readOfflineDataFromArchive< offlineDataInitialParticles >();
  
  if ( ptrInitialParticles->particleVector->size() != N_init )
  {
    string errMsg = "Wrong particle number";
    throw eConfig_error( errMsg );
  }
  
//   particles_init = *(ptrInitialParticles->particleVector);
  for ( int i = 0; i < ptrInitialParticles->particleVector->size(); i++ )
  {
    particles_init[i].unique_id = (*(ptrInitialParticles->particleVector))[i].unique_id;
    particles_init[i].FLAVOR = (*(ptrInitialParticles->particleVector))[i].FLAVOR;
    particles_init[i].m = (*(ptrInitialParticles->particleVector))[i].m;
    particles_init[i].T = (*(ptrInitialParticles->particleVector))[i].T;
    particles_init[i].X = (*(ptrInitialParticles->particleVector))[i].X;
    particles_init[i].Y = (*(ptrInitialParticles->particleVector))[i].Y;
    particles_init[i].Z = (*(ptrInitialParticles->particleVector))[i].Z;
    particles_init[i].E = (*(ptrInitialParticles->particleVector))[i].E;
    particles_init[i].PX = (*(ptrInitialParticles->particleVector))[i].PX;
    particles_init[i].PY = (*(ptrInitialParticles->particleVector))[i].PY;
    particles_init[i].PZ = (*(ptrInitialParticles->particleVector))[i].PZ;
    particles_init[i].md2g = (*(ptrInitialParticles->particleVector))[i].md2g;
    particles_init[i].md2q = (*(ptrInitialParticles->particleVector))[i].md2q;
  }

  for ( int i = 0; i < particles_init.size(); i++ )
  {   
    particles[i].init = true;
    particles[i].free = false;

    particles[i] = particles_init[i];
    
    if ( particles[i].T < timefirst )
    {
      particles[i].edge = true;
      particles[i].init = false;
    }

    particles[i].PXold = particles_init[i].PX;
    particles[i].PYold = particles_init[i].PY;
    particles[i].PZold = particles_init[i].PZ;
    particles[i].Eold  = particles_init[i].E;
  }

  particles_atTimeNow = particles;
}
