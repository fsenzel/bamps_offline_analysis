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
#include <string>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "config.h"
#include "configuration.h"
#include "particle.h"
#include "parameter.h"


using std::cout;
using std::endl;
using std::fstream;
using std::ios;
using std::string;
using std::vector;
namespace po = boost::program_options;
using namespace ns_casc;

/** @brief definition of particles vector, defined extern in configuration.h */
std::vector<ParticleOffline> ns_casc::particlesEvolving;
std::vector<ParticleOffline> ns_casc::particles_init;
std::vector<ParticleOffline> ns_casc::particles_atTimeNow;
std::vector<ParticleOffline> ns_casc::particles_atTimeNowCopy;
std::vector<ParticleOffline> ns_casc::addedParticles;
std::vector<ParticleOffline> ns_casc::addedParticlesCopy;


double dt = 0.1;

int MaxN, IX, IY, IZ;
int cellcut;
double deta;

initParam iparam;


/**
 * The only constructor provided for this class. Defaults for all parameters are set in case there is
 * no input file or it cannot be read.
 */
config::config() :
 configBase(),
// ---- simulation parameters ---- 
 freezeOutEnergyDensity(0.6),
// ---- initial state options ----
 initialStateType(miniJetsInitialState),
#ifdef LHAPDF_FOUND
 PDFsource( LHAPDF ),
#else
 PDFsource( builtInGRV ),
#endif
 LHAPDFdatasetName("cteq6l"),
 LHAPDFmember(0),
 LHAPDFuseGrid(false),
 nuclearPDFs(false),
 nuclearPDFdatasetName("EPS09"),
 pythiaParticleFile("-"),
 cgcParticleFile("-"),
 P0(1.4),
// ---- output options ----
 outputSwitch_progressLog( true ),
 outputSwitch_detailedParticleOutput(false),
 outputSwitch_movieOutputJets(false),
 outputSwitch_movieOutputBackground(false),
// ---- miscellaneous options ----
 switch_repeatTimesteps(true),
 jetMfpComputationSwitch(computeMfpDefault),
 interpolationBorder(50),
 // ---- offline reconstruction options ----
 pathdirOfflineData("offline_data"),
 originalName("default"),
 originalSeed(0),
 switch_fixed_dt(false),
 fixed_dt(-0.1),
 factor_dt(0.8),
 numberOfParticlesToAdd(0),
 minimumPT(10.0),
//---- program_options groups ----
 initial_state_options("Options and parameters for the initial state used by the BAMPS simulation"),
 offline_options("Offline reconstruction options")
{
  // populate the program_options::options_description groups with names and types of the possible parameters
  // this can add to the options_description groups provided by the base class
  initializeProgramOptions();
}


/**
 * @param[in] argc number of command line arguments, passed from the calling process
 * @param[in] argv[] command line arguments, passed from the calling process
 */ 
void config::readProgramOptions ( const int argc, char* argv[] )
{
  // since parsing the command line and the configuration file only relies on the groups "command_line_options",
  // "pos_options" and "config_file_options" to which the derived class adds but which are members of the base class
  // it is sufficient to call the routine of the base class here
  configBase::readProgramOptions( argc, argv );
}


/**
 * This provides the interface for the actual reading and processing of program options. It needs to be explicitly called
 * after an instance of the configPrototype class has been created.
 * 
 * @param[in] argc number of command line arguments, passed from the calling process
 * @param[in] argv[] command line arguments, passed from the calling process
 */
void config::readAndProcessProgramOptions ( const int argc, char* argv[] )
{
  groupProgramOptions();
  readProgramOptions(argc, argv);
  processProgramOptions();
    
  checkOptionsForSanity();
  
//   if( ParticlePrototype::N_heavy_flavor > 0 )
//   {
//     processHeavyQuarkOptions();
//   }
  
  printUsedConfigurationParameters();
}



/**
 * This is the place where post-processing of values that have been stored in the boost::program_options::variables_map
 * should go. For example conversion from integer to enum types etc.
 */
void config::processProgramOptions()
{
  configBase::processProgramOptions();
   
  if ( switch_fixed_dt )
  {
    if ( !(fixed_dt > 0) )
    {
      string errMsg = "Error: When \"offline.switch_fixed_dt = 1\" is specified, a value (>0) for \"offline.fixed_dt\" needs to be given.";
      throw eConfig_error( errMsg );
    }   
  }
  
  // some special conversions from integer type values to enum values
  if ( vm.count("initial_state.type") )
  {
    if ( vm["initial_state.type"].as<int>() < 3 && vm["initial_state.type"].as<int>() >= 0 )
    {
      initialStateType = static_cast<INITIAL_STATE_TYPE>( vm["initial_state.type"].as<int>() );
    }
    else
    {
      string errMsg = "parameter \"initial_state.type\" out of range";
      throw eConfig_error( errMsg );
    }
  }
  
  if ( vm.count("initial_state.PDFsource") )
  {
    if ( vm["initial_state.PDFsource"].as<unsigned short int>() < 2 && vm["initial_state.PDFsource"].as<unsigned short int>() >= 0 )
    {
      PDFsource = static_cast<PDF_SOURCE_TYPE>( vm["initial_state.PDFsource"].as<unsigned short int>() );
    }
    else
    {
      string errMsg = "parameter \"initial_state.PDFsource\" out of range";
      throw eConfig_error( errMsg );
    }
  }
  
  if ( vm.count("misc.jet_mfp_computation") )
  {
    if ( vm["misc.jet_mfp_computation"].as<int>() < 3 && vm["misc.jet_mfp_computation"].as<int>() >= 0 )
    {
      jetMfpComputationSwitch = static_cast<JET_MFP_COMPUTATION_TYPE>( vm["misc.jet_mfp_computation"].as<int>() );
    }
    else
    {
      string errMsg = "parameter \"misc.jet_mfp_computation\" out of range";
      throw eConfig_error( errMsg );      
    }
  }
}


/**
 * This routine initializes the boost::program_options::options_description objects that give structure to the handling of
 * user provided input.
 * By using a name pattern <group_name>.<option_name> the options can be given in a configuration file in INI format using 
 * group names such as
 * [group_name]
 * option_name1 = value1
 * option_name2 = value2
 */
void config::initializeProgramOptions()
{
  // Add some simulation parameters
  simulation_parameters.add_options()
  ("simulation.e_freeze", po::value<double>( &freezeOutEnergyDensity )->default_value( freezeOutEnergyDensity ), "freeze out energy density (GeV/fm^3)")
  ;
  
  // Group some options related to the initial state
  initial_state_options.add_options()
  ("initial_state.type", po::value<int>()->default_value( static_cast<int>(initialStateType) ), "initial state type (0 = mini-jets, 1 = pythia, 2 = cgc)")
  ("initial_state.PDFsource", po::value<unsigned short int>()->default_value( static_cast<unsigned short int>(PDFsource) ), "which source to use for the PDFs ( 0 = built-in GRV, 1 = PDFs from LHAPDF )")
  ("initial_state.LHAPDFset", po::value<string>( &LHAPDFdatasetName )->default_value( LHAPDFdatasetName ), "name of the LHAPDF data set that should be used")
  ("initial_state.LHAPDFmember", po::value<unsigned short int>( &LHAPDFmember )->default_value( LHAPDFmember ), "which member of the LHAPDF set should be used")
  ("initial_state.LHAPDFgrid", po::value<bool>( &LHAPDFuseGrid )->default_value( LHAPDFuseGrid ), "whether a grid version of the LHAPDF set should be used")
  ("initial_state.nuclearPDF", po::value<bool>( &nuclearPDFs )->default_value( nuclearPDFs ), "whether to use nuclear PDFs (only available together with LHAPDF and mini-jets)")
  ("initial_state.nuclearPDFname", po::value<string>( &nuclearPDFdatasetName )->default_value( nuclearPDFdatasetName ), "name of the nPDF dataset to use (EPS09, EPS09LO, EPS09NLO, EPS08, EKS98)")
  ("initial_state.minijet_P0", po::value<double>( &P0 )->default_value( P0 ), "lower pT cutoff for minijet initial conditions")
  ("initial_state.pythia_file", po::value<string>( &pythiaParticleFile )->default_value( pythiaParticleFile ), "input file providing pythia particle information, needed when initial_state.type = 1")
  ("initial_state.cgc_file", po::value<string>( &cgcParticleFile )->default_value( cgcParticleFile ), "input file providing cgc particle information, needed when initial_state.type = 2")
  ;
  
  // Add some options related to the program output  
  output_options.add_options()
  ("output.progress_log", po::value<bool>( &outputSwitch_progressLog )->default_value( outputSwitch_progressLog ), "write progress information" )
  ("output.particles", po::value<bool>( &outputSwitch_detailedParticleOutput )->default_value( outputSwitch_detailedParticleOutput ), "write detailed particle output")
  ("output.movie_jets", po::value<bool>( &outputSwitch_movieOutputJets )->default_value( outputSwitch_movieOutputJets ), "write movie output for added high-pt particles")
  ("output.movie_medium", po::value<bool>( &outputSwitch_movieOutputBackground )->default_value( outputSwitch_movieOutputBackground ), "write movie output for medium particles")
  ;
  
  // Add some miscellaneous options
  misc_options.add_options()
  ("misc.repeat_timesteps", po::value<bool>( &switch_repeatTimesteps )->default_value( switch_repeatTimesteps ), "repeat timesteps in cases where the probability has been > 1" ) 
  ("misc.interpolation_border", po::value<double>( &interpolationBorder )->default_value( interpolationBorder ), "X where interpolation of MFP is done for E > X*T")
  ("misc.jet_mfp_computation", po::value<int>()->default_value( jetMfpComputationSwitch ), "special treatment for the mean free path of high energy particles")
  ;
  
  // Group offline reconstruction options
  offline_options.add_options()
  ("offline.name", po::value<string>( &originalName )->default_value( originalName ), "name of the original BAMPS run that is to be reconstructed")
  ("offline.offline_data_dir", po::value<string>( &pathdirOfflineData )->default_value( pathdirOfflineData ), "directory from which the \"offline\" data is read")
  ("offline.use_fixed_dt", po::value<bool>( &switch_fixed_dt )->default_value( switch_fixed_dt ), "Indicates whether a fixed dt (provided via fixed_dt) should be used. Use time steps from the original run if not" ) 
  ("offline.fixed_dt", po::value<double>( &fixed_dt )->default_value( fixed_dt ), "fixed dt (time steps) at which the reconstructed medium is \"sampled\" [optional]")
  ("offline.factor_dt", po::value<double>( &factor_dt )->default_value( factor_dt ), "factor with which time steps from the original run should be scaled for use in sampling of the reconstructed medium (should be <1)")
  ("offline.nAdded", po::value<int>( &numberOfParticlesToAdd)->default_value( numberOfParticlesToAdd ), "number of (high-pt) particles that is added on top of the reconstructed medium, using it as a background" )
  ("offline.minPT_added", po::value<double>( &minimumPT )->default_value( minimumPT ), "minimum p_T [GeV] of the added particles")
  ;
}


/**
 * This routine groups the options_description groups into categories that can control which parameters are accessible via
 * the command line or the configuration file, which parameters are visible when requesting detailed help messages etc.
 */
void config::groupProgramOptions()
{
  configBase::groupProgramOptions(); // first add the options already contained in the base class
  
  // Add some groups that are meant to be provided via a configuration file
  config_file_options.add(initial_state_options).add(offline_options);
  
  // Add option groups that are to be printed in a detailed help message
  visible_options.add(initial_state_options).add(offline_options);
}



/**
 * This routine is provided in case future implementations should contain sanity checks for the given values of parameters and
 * options.
 */
void config::checkOptionsForSanity()
{
  configBase::checkOptionsForSanity();  // perform sanity check for base class options
  // sanity checks of parameters and options provided by the derived class can go here  
}



void config::processHeavyQuarkOptions()
{
  configBase::processHeavyQuarkOptions();  
}



/**
 * Print a file in INI format with all current option values. This file could be used a configuration filename
 * in a subsequent run to reproduce the very same option settings.
 * As it is somewhat easier this way, the routine explicitly implements the output of all option groups without calling
 * the corresponding base class routine.
 */
void config::printUsedConfigurationParameters()
{
  string filename = standardOutputDirectoryName + "/" + jobName + "_used_configuration";
  boost::filesystem::path outputFile( filename );
  boost::filesystem::ofstream output( outputFile, ios::trunc );
 
  printOptionsDescriptionToIniFormat( general_options, output );
  printOptionsDescriptionToIniFormat( simulation_parameters, output );
  printOptionsDescriptionToIniFormat( initial_state_options, output );
  printOptionsDescriptionToIniFormat( output_options, output );
  printOptionsDescriptionToIniFormat( misc_options, output );
//   printOptionsDescriptionToIniFormat( heavy_quark_options, output );
  printOptionsDescriptionToIniFormat( offline_options, output );
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
  // This overwrites the value of testparticles set by an input file since it must determined by the offline data. Effectively, specifying Ntest in input file has no consequences, but is present in configbase due to its need to specify in the input file in most BAMPS programs.
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
  N_light_flavors_offline = ptrSimulationData->N_light_flav;
  N_heavy_flavors_offline = ptrSimulationData->N_heavy_flav;
  
  
  Particle::set_N_light_flavor( N_light_flavors_offline );
  Particle::set_N_heavy_flavor( N_heavy_flavors_offline );
  
  cout << "** add " << numberOfParticlesToAdd << " particles with pt_min = " << minimumPT << endl;
  cout << timefirst << "  " << N_init << endl;
  
  cellcut = 4;
  //-------------------------//
  
  particlesEvolving.reserve( N_init * 1.8 );
  particles_init.reserve( N_init );
  particlesEvolving.resize( N_init * 1.8 );
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
    particlesEvolving[i] = particles_init[i];
    
    particlesEvolving[i].init = true;
    particlesEvolving[i].free = false;
    
    if ( particlesEvolving[i].T < timefirst )
    {
      particlesEvolving[i].edge = true;
      particlesEvolving[i].init = false;
    }
    
    particlesEvolving[i].PXold = particles_init[i].PX;
    particlesEvolving[i].PYold = particles_init[i].PY;
    particlesEvolving[i].PZold = particles_init[i].PZ;
    particlesEvolving[i].Eold  = particles_init[i].E;
  }
  
  particles_atTimeNow = particles_init;
}

