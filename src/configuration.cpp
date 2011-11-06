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
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "particle.h"
#include "configuration.h"
#include "analysis.h"
#include "parameter.h"

using std::cout;
using std::endl;
using std::fstream;
using std::ios;
using std::ifstream;
using namespace ns_casc;
namespace po = boost::program_options;

/** @brief definition of particles vector, defined extern in configuration.h */
std::vector<ParticleOffline> ns_casc::particles;
std::vector<ParticleOffline> ns_casc::particles_init;
std::vector<ParticleOffline> ns_casc::particles_atTimeNow;
std::vector<ParticleOffline> ns_casc::particles_atTimeNowCopy;
std::vector<ParticleOffline> ns_casc::addedParticles;
std::vector<ParticleOffline> ns_casc::addedParticlesCopy;


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
config::config( const int argc, char* argv[] )
{
  //--set default values---------------------------------------
  name = "default";
  runtime = 5.0; 
  seed = 0;
  jetMfpComputationSwitch = computeMfpDefault;
  interpolationBorder = 50;
    outputSwitch_detailedParticleOutput = false;
  freezeOutEnergyDensity = 1.0;
  pathdirOfflineData = "offline_data";
  pythiaParticleFile = "-";
  cgcParticleFile = "-";
  initialStateType = miniJetsInitialState;
  localCluster = true;
  movieOutputBackground = false;
  movieOutputJets = false;
  factor_dt = 0.8;
  dt_specified = false;
  //these values are used in case no input file is specified 
  //or in case the input file lacks certain statements
  //-----------------------------------------------------------
  
  // process command line arguments and parameters provided via an input file
  readAndProcessProgramOptions( argc, argv );
  
//   Particle::set_N_light_flavor( N_light_flavors_input );
//   Particle::set_N_heavy_flavor( N_heavy_flavors_input );
//   cout << "N_f = N_f_light_quarks + N_f_heavy_quarks = " << Particle::N_light_flavor << " + " <<   Particle::N_heavy_flavor << endl;
//   
//   checkOptionsForSanity();
//   
//   if( Particle::N_heavy_flavor > 0 )
//   {
//     processHeavyQuarkOptions();
//   }

  printUsedConfigurationParameters();
}


void config::readAndProcessProgramOptions( const int argc, char* argv[] )
{
  // Define some option groups. This gives more structure to the code and also structures the help message that is displayed
  // to the user.
  
  // Usage information. This is meant for the command line.
  po::options_description usage_information("usage information");
  usage_information.add_options()
  ("help,h", "print help message")
  ("detailed-help,d", "print detailed information on all options and parameters");
  
  // Read the name of the configuration file (from the command line)
  po::options_description hidden_options("Hidden options");
  hidden_options.add_options()
  ("config-file,c", po::value< string >(), "configuration file");
  
  // Make "config-file" a positional option, i.e. it can be used without explicitly specifying "--config-file=blabla.txt"
  po::positional_options_description pos_options;
  pos_options.add("config-file", 1);
  
//   // Group some general options  
//   po::options_description general_options("General options");
//   general_options.add_options()
//   ("general.jobname", po::value< string >( &jobName )->default_value( jobName ), "name of the simulation job, assigned to output files")
//   ("general.seed", po::value<long>( &seed )->default_value( seed ), "inital seed for random number generator (0 = pick random seed)")
//   ;
//   
//   // Group some collision parameters
//   po::options_description collision_parameters("Parameters of the heavy ion collision");
//   collision_parameters.add_options()
//   ("collision.A_mass", po::value<double>( &A )->default_value( A ), "mass number of nucleus A" )
//   ("collision.A_atomic", po::value<double>( &Aatomic )->default_value( Aatomic ), "atomic number (number of protons) of nucleus A" )
//   ("collision.B_mass", po::value<double>( &B )->default_value( B ), "mass number of nucleus B" )
//   ("collision.B_atomic", po::value<double>( &Batomic )->default_value( Batomic ), "atomic number (number of protons) of nucleus B" )
//   ("collision.impact", po::value<double>( &impact )->default_value( impact ), "impact parameter (fm)")
//   ("collision.sqrt_s", po::value<double>( &sqrtS )->default_value( sqrtS ), "center of mass energy per nucleon-nucleon pair (GeV)")
//   ;
//   
//   // Group some simulation parameters
//   po::options_description simulation_parameters("Paramters for the BAMPS simulation");
//   simulation_parameters.add_options()
//   ("simulation.time", po::value<double>( &runtime )->default_value( runtime ), "simulated time of the system evolution (fm/c)")
//   ("simulation.N_test", po::value<int>( &testparticles )->default_value( testparticles ), "number of test particles per real parton")
//   ("simulation.e_freeze", po::value<double>( &freezeOutEnergyDensity )->default_value( freezeOutEnergyDensity ), "freeze out energy density (GeV/fm^3)")
//   ("simulation.Nf_light", po::value<int>( &N_light_flavors_input )->default_value( N_light_flavors_input ), "number of light quark flavors")
//   ("simulation.Nf_heavy", po::value<int>( &N_heavy_flavors_input )->default_value( N_heavy_flavors_input ), "number of heavy quark flavors")
//   ;
//   
//   // Group some options related to the initial state
//   po::options_description initial_state_options("Options and parameters for the initial state used by the BAMPS simulation");
//   initial_state_options.add_options()
//   ("initial_state.type", po::value<int>()->default_value( static_cast<int>(initialStateType) ), "initial state type (0 = mini-jets, 1 = pythia, 2 = cgc)")
//   ("initial_state.minijet_P0", po::value<double>( &P0 )->default_value( P0 ), "lower pT cutoff for minijet initial conditions")
//   ("initial_state.pythia_file", po::value<string>( &pythiaParticleFile )->default_value( pythiaParticleFile ), "input file providing pythia particle information, needed when initial_state.type = 1")
//   ("initial_state.cgc_file", po::value<string>( &cgcParticleFile )->default_value( cgcParticleFile ), "input file providing cgc particle information, needed when initial_state.type = 2")
//   ;
//   
//   // Group some options related to the program output  
//   po::options_description output_options("Options for the output generated by the BAMPS simulation");
//   output_options.add_options()
//   ("output.directory", po::value<string>( &standardOutputDirectoryName )->default_value( standardOutputDirectoryName ), "directory to which general output should be written")
//   ("output.offline", po::value<bool>( &outputSwitch_offlineReconstructionData )->default_value( outputSwitch_offlineReconstructionData ), "write data for later reconstruction via \"offline\" routines")
//   ("output.offline_output_dir", po::value<string>( &pathdirOfflineData )->default_value( pathdirOfflineData ), "directory to which the output for \"offline\" data is written (if output.offline true)")
//   ("output.progress_log", po::value<bool>( &outputSwitch_progressLog )->default_value( outputSwitch_progressLog ), "write progress information" )
//   ("output.particles", po::value<bool>( &outputSwitch_detailedParticleOutput )->default_value( outputSwitch_detailedParticleOutput ), "write detailed particle output")
//   ("output.movie", po::value<bool>( &outputSwitch_movieOutput )->default_value( outputSwitch_movieOutput ), "write movie output")
//   ("output.v2", po::value<bool>( &outputSwitch_v2FinalData )->default_value( outputSwitch_v2FinalData ), "v2 output")
//   ("output.v2_intermediate", po::value<bool>( &outputSwitch_v2DataIntermediateSteps )->default_value( outputSwitch_v2DataIntermediateSteps ), "v2 output at intermediate steps")
//   ("output.v2_etabins", po::value<int>( &v2outputEtaBins )->default_value( v2outputEtaBins ), "number of eta bins used for v2 and RAA output")
//   ;
//   
//   // Group miscellaneous options
//   po::options_description misc_options("Miscellaneous options");
//   misc_options.add_options()
//   ("misc.repeat_timesteps", po::value<bool>( &switch_repeatTimesteps )->default_value( switch_repeatTimesteps ), "repeat timesteps in cases where the probability has been > 1" ) 
//   ("misc.interpolation_border", po::value<double>( &interpolationBorder )->default_value( interpolationBorder ), "X where interpolation of MFP is done for E > X*T")
//   ("misc.local_cluster", po::value<bool>( &localCluster )->default_value( localCluster ), "true if job needs to run on local ITP cluster, false otherwise")
//   ("misc.jet_mfp_computation", po::value<int>()->default_value( jetMfpComputationSwitch ), "special treatment for the mean free path of high energy particles")
//   ;
//   
//   // Group heavy quark options
//   po::options_description heavy_quark_options("Heavy quark options");
//   heavy_quark_options.add_options()
//   ("heavy_quark.charm_mass", po::value<double>( &Mcharm_input )->default_value( Mcharm_input ), "charm mass (GeV)")
//   ("heavy_quark.bottom_mass", po::value<double>( &Mbottom_input )->default_value( Mbottom_input ), "bottom mass (GeV)")
//   ("heavy_quark.Kfactor_ggQQbar", po::value<double>( &KggQQbar )->default_value( KggQQbar ), "K factor for process g + g -> Q + Qbar")
//   ("heavy_quark.Kfactor_gQgQ", po::value<double>( &KgQgQ )->default_value( KgQgQ ), "K factor for process g + Q -> g + Q")
//   ("heavy_quark.kappa_gQgQ", po::value<double>( &kappa_gQgQ )->default_value( kappa_gQgQ ), "multiplying factor kappa for the Debye mass in  g + Q -> g + Q processes")
//   ("heavy_quark.running_coupling", po::value<bool>( &couplingRunningForGQ )->default_value( couplingRunningForGQ ), "Running coupling for g+Q -> g+Q processes")
//   ("heavy_quark.iso_xsection", po::value<bool>( &isotropicCrossSecGQ )->default_value( isotropicCrossSecGQ ), "switch for isotropic cross section in g+Q -> g+Q processes")
//   ("heavy_quark.const_xsection", po::value<bool>( &constantCrossSecGQ )->default_value( constantCrossSecGQ ), "switch for constant cross section for g+Q -> g+Q processes")
//   ("heavy_quark.const_xsection_value", po::value<double>( &constantCrossSecValueGQ )->default_value( constantCrossSecValueGQ ), "value for constant g+Q -> g+Q cross section (mb)")
//   ;
  
  // Group options that are meant to be provided via the command line
  po::options_description command_line_options;
  command_line_options.add(usage_information).add(hidden_options);
  
//   // Group options that are meant to be provided via a configuration file
//   po::options_description config_file_options;
//   config_file_options.add(general_options).add(collision_parameters).add(simulation_parameters).add(initial_state_options).add(output_options).add(misc_options).add(heavy_quark_options);
//   
//   // Group options that are to be printed in a detailed help message
//   po::options_description visible_options;
//   visible_options.add(usage_information).add(general_options).add(collision_parameters).add(simulation_parameters).add(initial_state_options).add(output_options).add(misc_options).add(heavy_quark_options);
  
  
  // vm stores all the parameters that are obtained from the command line or the configuration file 
  po::variables_map vm;
  // parse the command line arguments
  po::store( po::command_line_parser(argc, argv).options(command_line_options).positional(pos_options).run(), vm );
  po::notify(vm); // needed to actually populate the map "vm" with all the argument values
  
  // check whether help messages are requested and print them
  if (vm.count("help"))
  {
    cout << "Usage: cascade [options] [configuration file]" << endl;
    cout << usage_information << "\n";
    throw( HELP_MESSAGE_REQUESTED );
  }
  if (vm.count("detailed-help"))
  {
//     cout << "Usage: cascade [options] [configuration file]" << endl;
//     cout << visible_options << "\n";
//     throw( HELP_MESSAGE_REQUESTED );
  }
  
  // check whether a configuration file name has been provided and parse this file if applicable
  if ( vm.count("config-file") )
  {
    boost::filesystem::path configFileName( vm["config-file"].as< string >() );
    if ( boost::filesystem::exists( configFileName ) )
    {
      cout << "Using configuration from: " << vm["config-file"].as< string >() << endl;
      boost::filesystem::ifstream configFile( configFileName );
//       po::store( po::parse_config_file(configFile, config_file_options), vm);
//       notify(vm);
    }
    else
    {
      string errMsg = "Configuration file \"" + vm["config-file"].as< string >() + "\" could not be found.";
      throw eConfig_error( errMsg );
    }
  }
  
//   // some special conversions from integer type values to enum values
//   if ( vm.count("initial_state.type") )
//   {
//     if ( vm["initial_state.type"].as<int>() < 3 && vm["initial_state.type"].as<int>() >= 0 )
//     {
//       initialStateType = static_cast<INITIAL_STATE_TYPE>( vm["initial_state.type"].as<int>() );
//     }
//     else
//     {
//       string errMsg = "parameter \"initial_state.type\" out of range";
//       throw eConfig_error( errMsg );
//     }
//   }
//   
//   if ( vm.count("misc.jet_mfp_computation") )
//   {
//     if ( vm["misc.jet_mfp_computation"].as<int>() < 3 && vm["misc.jet_mfp_computation"].as<int>() >= 0 )
//     {
//       jetMfpComputationSwitch = static_cast<JET_MFP_COMPUTATION_TYPE>( vm["misc.jet_mfp_computation"].as<int>() );
//     }
//     else
//     {
//       string errMsg = "parameter \"misc.jet_mfp_computation\" out of range";
//       throw eConfig_error( errMsg );      
//     }
//   }
  
}


void config::checkOptionsForSanity()
{
  // sanity checks of parameters and options can go here  
}



/**
* write out the values of the input parameters to cout
*/
void config::printUsedConfigurationParameters()
{
  string filename = standardOutputDirectoryName + "/" + jobName + "_used_configuration";
  boost::filesystem::path outputFile( filename );
  boost::filesystem::ofstream output( outputFile, ios::trunc );
  
//   output << "[general]" << endl;
//   output << "jobname = " << jobName << endl;
//   output << "seed = " << seed << endl;
//   output << endl;
//   
//   output << "[collision]" << endl;
//   output << "A_mass = " << A << endl;
//   output << "A_atomic = " << Aatomic << endl;
//   output << "B_mass = " << B << endl;
//   output << "B_atomic = " << Batomic << endl;
//   output << "impact = " << impact << endl;
//   output << "sqrt_s = " << sqrtS << endl;
//   output << endl;
//   
//   output << "[simulation]" << endl;
//   output << "time = " << runtime << endl;
//   output << "N_test = " << testparticles << endl;
//   output << "e_freeze = " << freezeOutEnergyDensity << endl;
//   output << "Nf_light = " << N_light_flavors_input << endl;
//   output << "Nf_heavy = " << N_heavy_flavors_input << endl;
//   output << endl;
//   
//   output << "[initial_state]" << endl;
//   output << "type = " << static_cast<int>(initialStateType) << endl;
//   output << "minijet_P0 = " << P0 << endl;
//   output << "pythia_file = " << pythiaParticleFile << endl;
//   output << "cgc_file = " << cgcParticleFile << endl;
//   output << endl;
//   
//   output << "[output]" << endl;
//   output << "directory = " << standardOutputDirectoryName << endl;
//   output << "offline = " << static_cast<int>(outputSwitch_offlineReconstructionData) << endl;
//   output << "offline_output_dir = " << pathdirOfflineData << endl;
//   output << "progress_log = " << static_cast<int>(outputSwitch_progressLog) << endl;
//   output << "particles = " << static_cast<int>(outputSwitch_detailedParticleOutput) << endl;
//   output << "movie = " << static_cast<int>( outputSwitch_movieOutput ) << endl;
//   output << "v2 = " << static_cast<int>( outputSwitch_v2FinalData ) << endl;
//   output << "v2_intermediate = " << static_cast<int>( outputSwitch_v2DataIntermediateSteps ) << endl;
//   output << "v2_etabins = " << v2outputEtaBins << endl;
//   output << endl;
//   
//   output << "[misc]" << endl;
//   output << "repeat_timesteps = " << static_cast<int>( switch_repeatTimesteps ) << endl;
//   output << "interpolation_border = " << interpolationBorder  << endl;
//   output << "local_cluster = " << static_cast<int>( localCluster ) << endl;
//   output << "jet_mfp_computation = " << static_cast<int>( jetMfpComputationSwitch ) << endl;
//   output << endl;
//  
//   output << "[heavy_quark]" << endl;
//   output << "charm_mass = " << Mcharm_input << endl;
//   output << "bottom_mass = " << Mbottom_input << endl;
//   output << "Kfactor_ggQQbar = " << KggQQbar << endl;
//   output << "Kfactor_gQgQ = " << KgQgQ << endl;
//   output << "kappa_gQgQ = "<< kappa_gQgQ << endl;
//   output << "running_coupling = " << static_cast<int>(couplingRunningForGQ) << endl;
//   output << "iso_xsection = " << static_cast<int>(isotropicCrossSecGQ) << endl;
//   output << "const_xsection = " << static_cast<int>(constantCrossSecGQ) << endl; 
//   output << "const_xsection_value = " << constantCrossSecValueGQ << endl;
//   output << endl;
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
