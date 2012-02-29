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
std::vector<ParticleOffline> ns_casc::particlesEvolving;
std::vector<ParticleOffline> ns_casc::particles_init;
std::vector<ParticleOffline> ns_casc::particles_atTimeNow;
std::vector<ParticleOffline> ns_casc::particles_atTimeNowCopy;
std::vector<ParticleOffline> ns_casc::addedParticles;
std::vector<ParticleOffline> ns_casc::addedParticlesCopy;
// std::vector<ParticleHFelectron> ns_casc::addedPartcl_electron;
std::vector<ParticleOffline> ns_casc::addedPartcl_electron;


double dt = 0.1;

int MaxN, IX, IY, IZ;
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
config::config( const int argc, char* argv[] ) :
// ---- general options ----
 jobName("default"),
 seed(0),
// ---- collision parameters ---- 
 A(197),
 Aatomic(79),
 B(197),
 Batomic(79),
 sqrtS(200),
 impact(0),
// ---- simulation parameters ---- 
 runtime(0.5),
 testparticles(35),
 freezeOutEnergyDensity(0.6),
 N_light_flavors_input(3),
 N_heavy_flavors_input(0), 
 scatt_offlineWithAddedParticles(true),
 scatt_amongOfflineParticles(false),
 scatt_amongAddedParticles(false),
 switchOff_23_32(false),
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
 mcatnloParticleFile("-"),
 P0(1.4),
// ---- output options ----
 outputSwitch_progressLog( true ),
 outputSwitch_detailedParticleOutput(false),
 outputSwitch_movieOutputJets(false),
 outputSwitch_movieOutputBackground(false),
 standardOutputDirectoryName("output"),
 v2RAAoutput(true),
 v2RAAoutputIntermediateSteps(false),
 dndyOutput(false),
 outputScheme(no_output),
// ---- heavy quark options ----
 KggQQbar(1.0),
 KgQgQ(1.0),
 kappa_gQgQ(1.0),
 couplingRunning(false),
 isotropicCrossSecGQ(false),
 constantCrossSecGQ(false),
 constantCrossSecValueGQ(10.0),
 Mcharm_input(1.5),
 Mbottom_input(4.8),
 hadronization_hq(false),
 mesonDecay(false),
 numberElectronStat(20),
 muonsInsteadOfElectrons(false),
 studyNonPromptJpsiInsteadOfElectrons(false),
 N_psi_input(0),
 isotropicCrossSecJpsi(false),
 constantCrossSecJpsi(false),
 constantCrossSecValueJpsi(10.0),
 Mjpsi_input(3.1),
 TdJpsi(0.32),
 jpsi_sigmaAbs(2.8),
 jpsi_agN(0.1),
 shadowing_model(eps08),
 jpsi_formationTime(0.0),
 hqCorrelationsOutput(false),
// ---- miscellaneous options ----
 switch_repeatTimesteps(true),
 jetMfpComputationSwitch(computeMfpDefault),
 interpolationBorder(50),
 localCluster(false),
// ---- offline reconstruction options ----
 pathdirOfflineData("offline_data"),
 originalName("default"),
 originalSeed(0),
 switch_fixed_dt(false),
 fixed_dt(-0.1),
 factor_dt(0.8),
 numberOfParticlesToAdd(0),
 numberOfAddedEvents(1),
 minimumPT(10.0)
{
  // process command line arguments and parameters provided via an input file
  readAndProcessProgramOptions( argc, argv );
  
  // default value of standardOutputDirectoryName is dependent on whether the calculation is performed at the CSC
  // However, it can be changed by hand in the config file
  //create the names of the output files for Fuchs CSC (old cluster) and LOEWE CSC
  char * csc_check_fuchs = getenv("PBS_JOBID");
  char * csc_check_loewe = getenv("SLURM_JOB_ID");
  if( csc_check_fuchs != NULL )
  {
    string jobID( csc_check_fuchs );
    standardOutputDirectoryName = "/local/" + jobID;
  }
  else if( csc_check_loewe != NULL )
  {
    string jobID( csc_check_loewe );
    standardOutputDirectoryName = "/local/" + jobID;
  }
  
  boost::filesystem::path outputDirectory( standardOutputDirectoryName );
  checkAndCreateOutputDirectory( outputDirectory );
  
  Particle::set_N_light_flavor( N_light_flavors_input );
  Particle::set_N_heavy_flavor( N_heavy_flavors_input );
  cout << "N_f of added  = N_f_light_quarks + N_f_heavy_quarks = " << Particle::N_light_flavor << " + " <<   Particle::N_heavy_flavor << endl;
  
  Particle::set_N_psi_states( N_psi_input );
  
  checkOptionsForSanity();
  
  if( Particle::N_heavy_flavor > 0 )
  {
    processHeavyQuarkOptions();
  }

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
  
  // Group some general options  
  po::options_description general_options("General options");
  general_options.add_options()
  ("general.jobname", po::value< string >( &jobName )->default_value( jobName ), "name of the simulation job, assigned to output files")
  ("general.seed", po::value<long>( &seed )->default_value( seed ), "inital seed for random number generator (0 = pick random seed)")
  ;
  
  // Group some simulation parameters
  po::options_description simulation_parameters("Paramters for the BAMPS simulation");
  simulation_parameters.add_options()
  ("simulation.time", po::value<double>( &runtime )->default_value( runtime ), "simulated time of the system evolution (fm/c)")
  ("simulation.Nf_light", po::value<int>( &N_light_flavors_input )->default_value( N_light_flavors_input ), "number of light quark flavors for added particles")
  ("simulation.Nf_heavy", po::value<int>( &N_heavy_flavors_input )->default_value( N_heavy_flavors_input ), "number of heavy quark flavors for added particles")
  ("simulation.scatt_offlineWithAdded", po::value<bool>( &scatt_offlineWithAddedParticles )->default_value( scatt_offlineWithAddedParticles ), "whether offline particles are allowed to scatter with added particles")
  ("simulation.scatt_amongOffline", po::value<bool>( &scatt_amongOfflineParticles )->default_value( scatt_amongOfflineParticles ), "whether offline particles are allowed to scatter with other offline particles")
  ("simulation.scatt_amongAdded", po::value<bool>( &scatt_amongAddedParticles )->default_value( scatt_amongAddedParticles ), "whether added particles are allowed to scatter with other added particles")
  ("simulation.switchOff_23_32", po::value<bool>( &switchOff_23_32 )->default_value( switchOff_23_32 ), "whether 2->3 and 3->2 processed are switched off for added particles")
  
  ;
  
  // Group some options related to the initial state
  po::options_description initial_state_options("Options and parameters for the initial state of added particles used by the BAMPS simulation");
  initial_state_options.add_options()
  ("initial_state.type", po::value<int>()->default_value( static_cast<int>(initialStateType) ), "initial state type (0 = mini-jets, 1 = pythia, 2 = cgc, 3 = mcatnlo)")
  ("initial_state.PDFsource", po::value<unsigned short int>()->default_value( static_cast<unsigned short int>(PDFsource) ), "which source to use for the PDFs ( 0 = built-in GRV, 1 = PDFs from LHAPDF )")
  ("initial_state.LHAPDFset", po::value<string>( &LHAPDFdatasetName )->default_value( LHAPDFdatasetName ), "name of the LHAPDF data set that should be used")
  ("initial_state.LHAPDFmember", po::value<unsigned short int>( &LHAPDFmember )->default_value( LHAPDFmember ), "which member of the LHAPDF set should be used")
  ("initial_state.LHAPDFgrid", po::value<bool>( &LHAPDFuseGrid )->default_value( LHAPDFuseGrid ), "whether a grid version of the LHAPDF set should be used")
  ("initial_state.nuclearPDF", po::value<bool>( &nuclearPDFs )->default_value( nuclearPDFs ), "whether to use nuclear PDFs (only available together with LHAPDF and mini-jets)")
  ("initial_state.nuclearPDFname", po::value<string>( &nuclearPDFdatasetName )->default_value( nuclearPDFdatasetName ), "name of the nPDF dataset to use (EPS09, EPS09LO, EPS09NLO, EPS08, EKS98)")
  ("initial_state.minijet_P0", po::value<double>( &P0 )->default_value( P0 ), "lower pT cutoff for minijet initial conditions")
  ("initial_state.pythia_file", po::value<string>( &pythiaParticleFile )->default_value( pythiaParticleFile ), "input file providing pythia particle information, needed when initial_state.type = 1")
  ("initial_state.cgc_file", po::value<string>( &cgcParticleFile )->default_value( cgcParticleFile ), "input file providing cgc particle information, needed when initial_state.type = 2")
  ("initial_state.mcatnlo_file", po::value<string>( &mcatnloParticleFile )->default_value( mcatnloParticleFile ), "input file providing MC@NLO particle information, needed when initial_state.type = 3")
  ;
  
  // Group some options related to the program output  
  po::options_description output_options("Options for the output generated by the BAMPS simulation");
  output_options.add_options()
  ("output.directory", po::value<string>( &standardOutputDirectoryName )->default_value( standardOutputDirectoryName ), "directory to which general output should be written")
  ("output.progress_log", po::value<bool>( &outputSwitch_progressLog )->default_value( outputSwitch_progressLog ), "write progress information" )
  ("output.particles", po::value<bool>( &outputSwitch_detailedParticleOutput )->default_value( outputSwitch_detailedParticleOutput ), "write detailed particle output")
  ("output.movie_jets", po::value<bool>( &outputSwitch_movieOutputJets )->default_value( outputSwitch_movieOutputJets ), "write movie output for added high-pt particles")
  ("output.movie_medium", po::value<bool>( &outputSwitch_movieOutputBackground )->default_value( outputSwitch_movieOutputBackground ), "write movie output for medium particles")
  ("output.v2RAA", po::value<bool>( &v2RAAoutput )->default_value( v2RAAoutput ), "write v2 and RAA output for added particles")
  ("output.v2RAAoutputIntermediateSteps", po::value<bool>( &v2RAAoutputIntermediateSteps )->default_value( v2RAAoutputIntermediateSteps ), "whether v2 and RAA output are printed at each analyisis time step (otherwise just at beginning and end)")
  ("output.dndyOutput", po::value<bool>( &dndyOutput )->default_value( dndyOutput ), "whether dndy output is written out")
  ("output.outputScheme", po::value<int>()->default_value( static_cast<int>(outputScheme) ), "output scheme id which configures the analysis routines and decides which output is written. The integer for the desired output scheme is given in the OUTPUT_SCHEME enum in configuration.h.")
  ;
  
  // Group heavy quark options
  po::options_description heavy_quark_options("Heavy quark options");
  heavy_quark_options.add_options()
  ("heavy_quark.charm_mass", po::value<double>( &Mcharm_input )->default_value( Mcharm_input ), "charm mass (GeV)")
  ("heavy_quark.bottom_mass", po::value<double>( &Mbottom_input )->default_value( Mbottom_input ), "bottom mass (GeV)")
  ("heavy_quark.Kfactor_ggQQbar", po::value<double>( &KggQQbar )->default_value( KggQQbar ), "K factor for process g + g -> Q + Qbar")
  ("heavy_quark.Kfactor_gQgQ", po::value<double>( &KgQgQ )->default_value( KgQgQ ), "K factor for process g + Q -> g + Q")
  ("heavy_quark.kappa_gQgQ", po::value<double>( &kappa_gQgQ )->default_value( kappa_gQgQ ), "multiplying factor kappa for the Debye mass in  g + Q -> g + Q processes")
  ("heavy_quark.running_coupling", po::value<bool>( &couplingRunning )->default_value( couplingRunning ), "Running coupling for all heavy flavor processes")
  ("heavy_quark.iso_xsection_gQ", po::value<bool>( &isotropicCrossSecGQ )->default_value( isotropicCrossSecGQ ), "switch for isotropic cross section in g+Q -> g+Q processes")
  ("heavy_quark.const_xsection_gQ", po::value<bool>( &constantCrossSecGQ )->default_value( constantCrossSecGQ ), "switch for constant cross section for g+Q -> g+Q processes")
  ("heavy_quark.const_xsection_gQ_value", po::value<double>( &constantCrossSecValueGQ )->default_value( constantCrossSecValueGQ ), "value for constant g+Q -> g+Q cross section (mb)")
  ("heavy_quark.hadronization_hq", po::value<bool>( &hadronization_hq )->default_value( hadronization_hq ), "whether hadronization of heavy quarks to D and B mesons is carried out")
  ("heavy_quark.mesonDecay", po::value<bool>( &mesonDecay )->default_value( mesonDecay ), "whether decay of heavy mesons from charm and bottom to electrons is performed")
  ("heavy_quark.numberElectronStat", po::value<int>( &numberElectronStat )->default_value( numberElectronStat ), "one meson decays to numberElectronStat electrons")
  ("heavy_quark.muonsInsteadOfElectrons", po::value<bool>( &muonsInsteadOfElectrons )->default_value( muonsInsteadOfElectrons ), "decay to muons should be performed instead of electrons")
  ("heavy_quark.studyNonPromptJpsiInsteadOfElectrons", po::value<bool>( &studyNonPromptJpsiInsteadOfElectrons )->default_value( studyNonPromptJpsiInsteadOfElectrons ), "decay of B mesons to non prompt Jpsi should be performed instead to electrons")
  ("heavy_quark.N_psi", po::value<int>( &N_psi_input )->default_value( N_psi_input ), "number of active psi states (1 = J/psi)")
  ("heavy_quark.iso_xsection_jpsi", po::value<bool>( &isotropicCrossSecJpsi )->default_value( isotropicCrossSecJpsi ), "isotropic momentum sampling is employed for process Q + Qbar -> g + J/psi")
  ("heavy_quark.const_xsection_jpsi", po::value<bool>( &constantCrossSecJpsi )->default_value( constantCrossSecJpsi ), "constant cross section is employed for process Q + Qbar -> g + J/psi")
  ("heavy_quark.const_xsection_jpsi_value", po::value<double>( &constantCrossSecValueJpsi )->default_value( constantCrossSecValueJpsi ), "value of constant cross section for process Q + Qbar -> g + J/psi in mb")
  ("heavy_quark.jpsi_mass", po::value<double>( &Mjpsi_input )->default_value( Mjpsi_input ), "J/psi mass (GeV)")
  ("heavy_quark.TdJpsi", po::value<double>( &TdJpsi )->default_value( TdJpsi ), "dissociation temperature of J/psi")
  ("heavy_quark.jpsi_sigmaAbs", po::value<double>( &jpsi_sigmaAbs )->default_value( jpsi_sigmaAbs ), "absorption cross section for initial J/psi in mb")
  ("heavy_quark.jpsi_agN", po::value<double>( &jpsi_agN )->default_value( jpsi_agN ), "parameter for momentum broadening of initial J/psi")
  ("heavy_quark.shadowing_model", po::value<int>( )->default_value( static_cast<int>(shadowing_model) ), "shadowing model used for initial J/psi")
  ("heavy_quark.jpsi_formationTime", po::value<double>( &jpsi_formationTime )->default_value( jpsi_formationTime ), "formation time for initial J/psi in addition to the standard 1/M_T")
  ("heavy_quark.hqCorrelationsOutput", po::value<bool>( &hqCorrelationsOutput )->default_value( hqCorrelationsOutput ), "whether correlation analysis of heavy quark pairs is done")
  ;
  
  // Group miscellaneous options
  po::options_description misc_options("Miscellaneous options");
  misc_options.add_options()
  ("misc.repeat_timesteps", po::value<bool>( &switch_repeatTimesteps )->default_value( switch_repeatTimesteps ), "repeat timesteps in cases where the probability has been > 1" ) 
  ("misc.interpolation_border", po::value<double>( &interpolationBorder )->default_value( interpolationBorder ), "X where interpolation of MFP is done for E > X*T")
  ("misc.local_cluster", po::value<bool>( &localCluster )->default_value( localCluster ), "true if job needs to run on local ITP cluster, false otherwise")
  ("misc.jet_mfp_computation", po::value<int>()->default_value( jetMfpComputationSwitch ), "special treatment for the mean free path of high energy particles")
  ;

  // Group offline reconstruction options
  po::options_description offline_options("Offline reconstruction options");
  offline_options.add_options()
  ("offline.name", po::value<string>( &originalName )->default_value( originalName ), "name of the original BAMPS run that is to be reconstructed")
  ("offline.offline_data_dir", po::value<string>( &pathdirOfflineData )->default_value( pathdirOfflineData ), "directory from which the \"offline\" data is read")
  ("offline.use_fixed_dt", po::value<bool>( &switch_fixed_dt )->default_value( switch_fixed_dt ), "Indicates whether a fixed dt (provided via fixed_dt) should be used. Use time steps from the original run if not" ) 
  ("offline.fixed_dt", po::value<double>( &fixed_dt ), "fixed dt (time steps) at which the reconstructed medium is \"sampled\" [optional]")
  ("offline.factor_dt", po::value<double>( &factor_dt )->default_value( factor_dt ), "factor with which time steps from the original run should be scaled for use in sampling of the reconstructed medium (should be <1)")
  ("offline.nAdded", po::value<int>( &numberOfParticlesToAdd)->default_value( numberOfParticlesToAdd ), "number of (high-pt) particles that is added on top of the reconstructed medium, using it as a background" )
  ("offline.numberOfAddedEvents", po::value<int>( &numberOfAddedEvents)->default_value( numberOfAddedEvents ), "number of heavy ion collision events, set on top of the offline reconstruction. This input is neede if Pythia data files are read in to normalize the total yield of the particles. The number of one such heavy ion collision event includes < number of produced particles in pp > * Ntest * Nbin" )
  ("offline.minPT_added", po::value<double>( &minimumPT )->default_value( minimumPT ), "minimum p_T [GeV] of the added particles. If p_T of added particle falls below this value it does not scatter anymore.")
  ;

  // Group options that are meant to be provided via the command line
  po::options_description command_line_options;
  command_line_options.add(usage_information).add(hidden_options);
  
  // Group options that are meant to be provided via a configuration file
  po::options_description config_file_options;
  config_file_options.add(general_options).add(simulation_parameters).add(initial_state_options).add(output_options).add(heavy_quark_options).add(misc_options).add(offline_options);
  
  // Group options that are to be printed in a detailed help message
  po::options_description visible_options;
  visible_options.add(usage_information).add(general_options).add(simulation_parameters).add(initial_state_options).add(output_options).add(heavy_quark_options).add(misc_options).add(offline_options);
  
  
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
    cout << "Usage: cascade [options] [configuration file]" << endl;
    cout << visible_options << "\n";
    throw( HELP_MESSAGE_REQUESTED );
  }
  
  // check whether a configuration file name has been provided and parse this file if applicable
  if ( vm.count("config-file") )
  {
    boost::filesystem::path configFileName( vm["config-file"].as< string >() );
    if ( boost::filesystem::exists( configFileName ) )
    {
      cout << "Using configuration from: " << vm["config-file"].as< string >() << endl;
      boost::filesystem::ifstream configFile( configFileName );
      po::store( po::parse_config_file(configFile, config_file_options), vm);
      notify(vm);
    }
    else
    {
      string errMsg = "Configuration file \"" + vm["config-file"].as< string >() + "\" could not be found.";
      throw eConfig_error( errMsg );
    }
  }
  
  // automatically set switch_fixed_dt when fixed_dt is set 
  if ( vm.count("offline.fixed_dt") )
  {
    switch_fixed_dt = true;
  }

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
    if ( vm["initial_state.type"].as<int>() < 4 && vm["initial_state.type"].as<int>() >= 0 )
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
  
  if ( vm.count("heavy_quark.shadowing_model") )
  {
    if ( vm["heavy_quark.shadowing_model"].as<int>() < 3 && vm["heavy_quark.shadowing_model"].as<int>() >= 0 )
    {
      shadowing_model = static_cast<shadowModelJpsi>( vm["heavy_quark.shadowing_model"].as<int>() );
    }
    else
    {
      string errMsg = "parameter \"heavy_quark.shadowing_model\" out of range";
      throw eConfig_error( errMsg );
    }
  }
  
  if ( vm.count("output.outputScheme") )
  {
    outputScheme = static_cast<OUTPUT_SCHEME>( vm["output.outputScheme"].as<int>() );
  }
  
}


void config::checkOptionsForSanity()
{
  // sanity checks of parameters and options can go here  
  
  if( ( scatt_amongAddedParticles && Particle::N_psi_states == 0 ) || ( Particle::N_psi_states > 0 && !scatt_amongAddedParticles ) )
  {
//     string errMsg = "Scatterings among added particles and Jpsi not active or vice versa.";
//     throw eConfig_error( errMsg );
  }
  
  if( mesonDecay && !hadronization_hq )
  {
    string errMsg = "Meson decay not possible if heavy quarks cannot hadronize.";
    throw eConfig_error( errMsg );
  }
  
  if( ( muonsInsteadOfElectrons || studyNonPromptJpsiInsteadOfElectrons ) && !mesonDecay )
  {
    string errMsg = "Meson decay to other particles does not make sense if meson decay is switched of.";
    throw eConfig_error( errMsg );
  }
  
  if( ( studyNonPromptJpsiInsteadOfElectrons && outputScheme != cms_hq_nonPromptJpsi ) || ( outputScheme == cms_hq_nonPromptJpsi && !studyNonPromptJpsiInsteadOfElectrons ) )
  {
    string errMsg = "Study Jpsi instead of electrons but no output for this or vice versa.";
    throw eConfig_error( errMsg );
  }
}


void config::processHeavyQuarkOptions()
{
  coupling::set_Nflavor( N_light_flavors_input );
  
  coupling::set_isRunning( couplingRunning );
  Particle::setCharmMass( Mcharm_input );
  cout << "Charm mass: " << Particle::Mcharm;
  if( Particle::N_heavy_flavor > 1 )
  {
    Particle::setBottomMass( Mbottom_input );
    cout << "  bottom mass: " << Particle::Mbottom;
  }
  cout << endl;
  
  if(constantCrossSecGQ)
  {
    cout << "Constant cross section of " << constantCrossSecValueGQ << " mb for gQ-> gQ." << endl;
  }
  else
  {
    cout << "pQCD cross section for gQ-> gQ with K = " << KgQgQ;
    if(couplingRunning)
    {
      cout << " and running coupling";
    }
    else
    {
      cout << " and constant coupling";
    }
    cout << ", kappa_gQ = " << kappa_gQgQ << endl;
  }
  
  if(isotropicCrossSecGQ)
  {
    cout << "Sample angle isotropic for gQ-> gQ." << endl;
  }

  if( mesonDecay )
  {
    if(muonsInsteadOfElectrons)
      cout << "Decay of heavy mesons to muons and not to electrons. " << endl;
    if(studyNonPromptJpsiInsteadOfElectrons)
      cout << "Decay of B mesons to non prompt Jpsi and not to electrons. " << endl;
    if( muonsInsteadOfElectrons && studyNonPromptJpsiInsteadOfElectrons )
    {
      string errMsg = "Decay to muons and non prompt Jpsi not possible.";
      throw eConfig_error( errMsg );
    }
  }


  // study Jpsi
  if( Particle::N_psi_states > 0 )
  {
    Particle::setJpsiMass( Mjpsi_input );
    
    if(constantCrossSecJpsi)
      cout << "Constant cross section of " << constantCrossSecValueJpsi << " mb for ccbar -> Jpsi + g." << endl;
    else
      cout << "Cross section from Peskin for ccbar -> Jpsi + g." << endl;

    if(isotropicCrossSecJpsi)
      cout << "Sample angle isotropic for ccbar -> Jpsi + g." << endl;
    
    cout << "TdJpsi = " << TdJpsi << endl;

    if( jpsi_formationTime != 0.0 )
      cout << "Jpsi formation time: " << jpsi_formationTime << endl;
  }
}


/**
* write out the values of the input parameters to cout
*/
void config::printUsedConfigurationParameters()
{
  string filename = standardOutputDirectoryName + "/" + jobName + "_used_configuration";
  boost::filesystem::path outputFile( filename );
  boost::filesystem::ofstream output( outputFile, ios::trunc );
  
  output << "[general]" << endl;
  output << "jobname = " << jobName << endl;
  output << "seed = " << seed << endl;
  output << endl;
  
  output << "[simulation]" << endl;
  output << "time = " << runtime << endl;
  output << "Nf_light = " << N_light_flavors_input << endl;
  output << "Nf_heavy = " << N_heavy_flavors_input << endl;
  output << "scatt_offlineWithAdded = " << scatt_offlineWithAddedParticles << endl;
  output << "scatt_amongOffline = " << scatt_amongOfflineParticles << endl;
  output << "scatt_amongAdded = " << scatt_amongAddedParticles << endl;
  output << "switchOff_23_32 = " << switchOff_23_32 << endl;
  output << endl;
  
  output << "[initial_state]" << endl;
  output << "type = " << static_cast<int>(initialStateType) << endl;
  output << "PDFsource = " << static_cast<int>(PDFsource) << endl;
  output << "LHAPDFset = " << LHAPDFdatasetName << endl;
  output << "LHAPDFmember = " << LHAPDFmember << endl;
  output << "LHAPDFgrid = " << static_cast<int>(LHAPDFuseGrid) << endl;
  output << "nuclearPDF = " << static_cast<int>(nuclearPDFs) << endl;
  output << "nuclearPDFname = " << nuclearPDFdatasetName << endl;
  output << "minijet_P0 = " << P0 << endl;
  output << "pythia_file = " << pythiaParticleFile << endl;
  output << "cgc_file = " << cgcParticleFile << endl;
  output << "mcatnlo_file = " << mcatnloParticleFile << endl;
  output << endl;
  
  output << "[output]" << endl;
  output << "directory = " << standardOutputDirectoryName << endl;
  output << "progress_log = " << static_cast<int>(outputSwitch_progressLog) << endl;
  output << "particles = " << static_cast<int>(outputSwitch_detailedParticleOutput) << endl;
  output << "movie_jets = " << static_cast<int>( outputSwitch_movieOutputJets ) << endl;
  output << "movie_medium = " << static_cast<int>( outputSwitch_movieOutputBackground ) << endl;
  output << "v2RAA = " << static_cast<int>( v2RAAoutput ) << endl;
  output << "v2RAAoutputIntermediateSteps = " << static_cast<int>( v2RAAoutputIntermediateSteps ) << endl;
  output << "dndyOutput = " << static_cast<int>( dndyOutput ) << endl;
  output << "outputScheme = " << static_cast<int>( outputScheme ) << endl;
  output << endl;
  
  output << "[heavy_quark]" << endl;
  output << "charm_mass = " << Mcharm_input << endl;
  output << "bottom_mass = " << Mbottom_input << endl;
  output << "Kfactor_ggQQbar = " << KggQQbar << endl;
  output << "Kfactor_gQgQ = " << KgQgQ << endl;
  output << "kappa_gQgQ = " << kappa_gQgQ << endl;
  output << "running_coupling = " << couplingRunning << endl;
  output << "iso_xsection_gQ = " << isotropicCrossSecGQ << endl;
  output << "const_xsection_gQ = " << constantCrossSecGQ << endl;
  output << "const_xsection_gQ_value = " << constantCrossSecValueGQ << endl;
  output << "hadronization_hq = " << hadronization_hq << endl;
  output << "mesonDecay = " << mesonDecay << endl;
  output << "numberElectronStat = " << numberElectronStat << endl;
  output << "muonsInsteadOfElectrons = " << muonsInsteadOfElectrons << endl;
  output << "studyNonPromptJpsiInsteadOfElectrons = " << studyNonPromptJpsiInsteadOfElectrons << endl;
  output << "N_psi = " << N_psi_input << endl;
  output << "iso_xsection_jpsi = " << isotropicCrossSecJpsi << endl;
  output << "const_xsection_jpsi = " << constantCrossSecJpsi << endl;
  output << "const_xsection_jpsi_value = " << constantCrossSecValueJpsi << endl;
  output << "jpsi_mass = " << Mjpsi_input << endl;
  output << "TdJpsi = " << TdJpsi << endl;
  output << "jpsi_sigmaAbs = " << jpsi_sigmaAbs << endl;
  output << "jpsi_agN = " << jpsi_agN << endl;
  output << "shadowing_model = " << static_cast<int>( shadowing_model ) << endl;
  output << "jpsi_formationTime = " << jpsi_formationTime << endl;
  output << "hqCorrelationsOutput = " << hqCorrelationsOutput << endl;
  output << endl;
  
  output << "[misc]" << endl;
  output << "repeat_timesteps = " << static_cast<int>( switch_repeatTimesteps ) << endl;
  output << "interpolation_border = " << interpolationBorder  << endl;
  output << "local_cluster = " << static_cast<int>( localCluster ) << endl;
  output << "jet_mfp_computation = " << static_cast<int>( jetMfpComputationSwitch ) << endl;
  output << endl;

  output << "[offline]" << endl;
  output << "name = " << originalName << endl;
  output << "offline_data_dir = " << pathdirOfflineData << endl;
  output << "use_fixed_dt = " << static_cast<int>( switch_fixed_dt ) << endl;
  output << "fixed_dt = " << fixed_dt  << endl;
  output << "factor_dt = " << factor_dt << endl;
  output << "nAdded = " << numberOfParticlesToAdd << endl;
  output << "numberOfAddedEvents = " << numberOfAddedEvents << endl;
  output << "minPT_added = " << minimumPT  << endl;
  output << endl;
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
  N_light_flavors_offline = ptrSimulationData->N_light_flav;
  N_heavy_flavors_offline = ptrSimulationData->N_heavy_flav;
  
  // check if the number of flavors of the offline program is larger than that for added particles. If so set the global particle property accordingly to avoid errors. Particle::N_light_flavor must always be larger/equal the largest actually occuring particle flavor
  if( N_light_flavors_offline > Particle::N_light_flavor )
    Particle::set_N_light_flavor( N_light_flavors_offline );
  if( N_heavy_flavors_offline > Particle::N_heavy_flavor )
    Particle::set_N_heavy_flavor( N_heavy_flavors_offline );
  
  cout << "N_f of medium = N_f_light_quarks + N_f_heavy_quarks = " << N_light_flavors_offline << " + " <<   N_heavy_flavors_offline << endl;
  cout << "N_f overall   = N_f_light_quarks + N_f_heavy_quarks = " << Particle::N_light_flavor << " + " <<   Particle::N_heavy_flavor << endl;
  
  if( N_light_flavors_offline < N_light_flavors_input )
    cout << "There are more light flavors for added particles (" << N_light_flavors_input << ") switched on than for the medium (" << N_light_flavors_offline << "). It is assumed that you know what you do..." << endl;
  
  numberOfParticlesToAdd *= testparticles;
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
