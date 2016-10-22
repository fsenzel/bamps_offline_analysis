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

#include "configBAMPS.h"
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

/** @brief definition of particles vector, declared extern in configuration.h */
std::vector<ParticleOffline> ns_casc::particlesEvolving;
std::vector<ParticleOffline> ns_casc::particles_init;
std::vector<ParticleOffline> ns_casc::particles_atTimeNow;
std::vector<ParticleOffline> ns_casc::particles_atTimeNowCopy;
std::vector<ParticleOffline> ns_casc::addedParticles;
std::vector<ParticleOffline> ns_casc::addedParticlesCopy;
// std::vector<ParticleHFelectron> ns_casc::addedPartcl_electron;
std::vector<ParticleOffline> ns_casc::addedPartcl_electron;
std::vector<ParticleOffline> ns_casc::scatteredMediumParticles;
std::vector<ParticleOffline> ns_casc::noninteractingParticles;
std::vector<ParticleOffline> ns_casc::dileptons;


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
    scatt_offlineWithAddedParticles(true),
    scatt_amongOfflineParticles(false),
    scatt_amongAddedParticles(false),
    scatt_amongBackgroundParticles(false),
    scatt_furtherOfflineParticles( false ),
    jet_tagged( false ),
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
    outputSwitch_scatteredMediumParticlesOutput( false ),
    v2RAAoutput(true),
    v2RAAoutputIntermediateSteps(false),
    dndyOutput(false),
    outputSwitch_photons( false ),
    outputSwitch_QCDparticles( true ),
    outputScheme(no_output),
    analysisTubeRadius( 1.5 ),
    analysisTubedEta( 0.5 ),
// ---- heavy quark options ----
    couplingRunningHeavyQuarksInput(false),
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
    jpsi_testparticles(1),
    hqCorrelationsOutput(false),
// ---- miscellaneous options ----
    switch_repeatTimesteps(true),
    jetMfpComputationSwitch(computeMfpLastTimestep),
    mfpAddedRangeVariation( 100.0 ),
    fixed_mfp_added( 1.0 ),
    onlyMediumEvolution( false ),
    analysisPhotonsTimeCut( 0.0 ),
    analysisPhotonsPTCut( 0.0 ),
    restrictParentPTForPhotonproduction(false),
    minAllowedParentPT( 0.0 ),
    maxAllowedParentPT( 100000.0 ),
    analysisPhotonsNBinsV2(8),
    mfpCellAveraging( false ),
    usedExternalField( 0.0),
//  interpolationBorder(50),
// ---- offline reconstruction options ----
    pathdirOfflineData("offline_data"),
    originalName("default"),
    originalSeed(0),
    switch_fixed_dt(false),
    fixed_dt(-0.1),
    factor_dt(0.8),//WARNING
    numberOfParticlesToAdd(0),
    numberOfAddedEvents(1),
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

    if( Particle::N_heavy_flavor > 0 )
    {
        processHeavyQuarkOptions();
    }

    printUsedConfigurationParameters();
}



/**
 * This is the place where post-processing of values that have been stored in the boost::program_options::variables_map
 * should go. For example conversion from integer to enum types etc.
 */
void config::processProgramOptions()
{
    configBase::processProgramOptions();

    Particle::set_N_psi_states( N_psi_input );

    // check if both initial number of particles and initial number of events are set which does not make sense
    bool ini_number_added_set = false;
    bool ini_number_added_events_set = false;
    if( !( vm["offline.nAdded"].defaulted() || vm["offline.nAdded"].empty() ) )
    {
        ini_number_added_set = true;
    }
    if( !( vm["offline.nAddedEvents"].defaulted() || vm["offline.nAddedEvents"].empty() ) )
    {
        ini_number_added_events_set = true;
    }

    // if both are given and the number of particles is not negative (which means use number of added events), throw error.
    if( ini_number_added_set && ini_number_added_events_set && numberOfParticlesToAdd >= 0 )
    {
        string errMsg = "Both initial number of particles and initial number of events are set. Do not know which one to use.";
        throw eConfig_error( errMsg );
    }
    else if( !ini_number_added_set && ini_number_added_events_set )
    {
        // set initial number of particles to a negative value to use the number of added events instead
        numberOfParticlesToAdd = -1;
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
        if ( vm["initial_state.type"].as<int>() < 8 && vm["initial_state.type"].as<int>() >= 0 )
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
        if ( vm["misc.jet_mfp_computation"].as<int>() < 5 && vm["misc.jet_mfp_computation"].as<int>() >= 0 )
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
        if ( vm["heavy_quark.shadowing_model"].as<int>() < 4 && vm["heavy_quark.shadowing_model"].as<int>() >= 0 )
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

    // check if simulation.jet_tagged is set in input file: if not default value or no value given (latter happens if no input file is given at all)
    // If it is not set and N_heavy_flavors > 0 and no light flavors, set the property jet_tagged = true.
    if ( ( vm["simulation.jet_tagged"].defaulted() || vm["simulation.jet_tagged"].empty() ) &&  N_light_flavors_input < 0 &&
            N_heavy_flavors_input > 0 )
    {
        cout << "Only heavy quarks are switched on. Therefore, set jet_tagged = true." << endl;
        jet_tagged = true;
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
    ("simulation.scatt_offlineWithAdded", po::value<bool>( &scatt_offlineWithAddedParticles )->default_value( scatt_offlineWithAddedParticles ), "whether offline particles are allowed to scatter with added particles")
    ("simulation.scatt_amongOffline", po::value<bool>( &scatt_amongOfflineParticles )->default_value( scatt_amongOfflineParticles ), "whether offline particles are allowed to scatter with other offline particles")
    ("simulation.scatt_amongAdded", po::value<bool>( &scatt_amongAddedParticles )->default_value( scatt_amongAddedParticles ), "whether added particles are allowed to scatter with other added particles")
    ("simulation.scatt_amongBackground", po::value<bool>( &scatt_amongBackgroundParticles )->default_value( scatt_amongBackgroundParticles ), "whether added particles are allowed to scatter with other added particles")
    ("simulation.scatt_furtherOffline", po::value<bool>( &scatt_furtherOfflineParticles )->default_value( scatt_furtherOfflineParticles ), "whether the background can produce photons by MC-sampling only among background particles")
    ("simulation.jet_tagged", po::value<bool>( &jet_tagged )->default_value( jet_tagged ), "whether the added particles are treated as tagges jets")
    ;

    // Group some options related to the initial state
    initial_state_options.add_options()
    ("initial_state.type", po::value<int>()->default_value( static_cast<int>(initialStateType) ), "initial state type (0 = mini-jets, 1 = pythia, 2 = cgc, 3 = mcatnlo, 4 = onlyJpsi)")
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
    ("initial_state.initialPartonPt", po::value<double>( &initialPartonPt )->default_value( initialPartonPt ), "parton pt of fixed initial parton pt")
    ("initial_state.initialPartonFlavor", po::value<int>( &initialPartonFlavor )->default_value( initialPartonFlavor ), "flavor of fixed initial parton pt ( 0 = gg, 1 = u u-bar, 2 = g u )")
    ;

    // Add some options related to the program output
    output_options.add_options()
    ("output.progress_log", po::value<bool>( &outputSwitch_progressLog )->default_value( outputSwitch_progressLog ), "write progress information" )
    ("output.particles", po::value<bool>( &outputSwitch_detailedParticleOutput )->default_value( outputSwitch_detailedParticleOutput ), "write detailed particle output")
    ("output.movie_jets", po::value<bool>( &outputSwitch_movieOutputJets )->default_value( outputSwitch_movieOutputJets ), "write movie output for added high-pt particles")
    ("output.movie_medium", po::value<bool>( &outputSwitch_movieOutputBackground )->default_value( outputSwitch_movieOutputBackground ), "write movie output for medium particles")
    ("output.v2RAA", po::value<bool>( &v2RAAoutput )->default_value( v2RAAoutput ), "write v2 and RAA output for added particles")
    ("output.v2RAAoutputIntermediateSteps", po::value<bool>( &v2RAAoutputIntermediateSteps )->default_value( v2RAAoutputIntermediateSteps ), "whether v2 and RAA output are printed at each analyisis time step (otherwise just at beginning and end)")
    ("output.dndyOutput", po::value<bool>( &dndyOutput )->default_value( dndyOutput ), "whether dndy output is written out")
    ("output.photons", po::value<bool>( &outputSwitch_photons )->default_value( outputSwitch_photons ), "whether output concerning photons is written out")
    ("output.QCDparticles", po::value<bool>( &outputSwitch_QCDparticles )->default_value( outputSwitch_QCDparticles ), "whether output concerning QCD particles (gluons, quarks, etc.) is written out")
    ("output.outputScheme", po::value<int>()->default_value( static_cast<int>(outputScheme) ), "output scheme id which configures the analysis routines and decides which output is written. The integer for the desired output scheme is given in the OUTPUT_SCHEME enum in configuration.h.")
    ("output.scatteredMedium", po::value<bool>( &outputSwitch_scatteredMediumParticlesOutput )->default_value( outputSwitch_scatteredMediumParticlesOutput ), "write scattered medium particles output")
    ("output.radiusAnalysisTube", po::value<double>( &analysisTubeRadius )->default_value( analysisTubeRadius ), "Radius in fm of the central tube for medium analysis.")
    ("output.dEtaAnalysisTube", po::value<double>( &analysisTubedEta )->default_value( analysisTubedEta ), "dEta for the central tube for medium analysis.")
    ;

    // Add heavy quark options
    heavy_quark_options.add_options()
    ("heavy_quark.running_coupling", po::value<bool>( &couplingRunningHeavyQuarksInput )->default_value( couplingRunningHeavyQuarksInput ), "Running coupling for all heavy flavor processes")
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
    ("heavy_quark.shadowing_model", po::value<int>( )->default_value( static_cast<int>(shadowing_model) ), "shadowing model used for initial J/psi: 0 = none, 1 = eks98, 2 = eps08, 3 = eps09")
    ("heavy_quark.jpsi_formationTime", po::value<double>( &jpsi_formationTime )->default_value( jpsi_formationTime ), "formation time for initial J/psi in addition to the standard 1/M_T")
    ("heavy_quark.jpsi_testparticles", po::value<int>( &jpsi_testparticles )->default_value( jpsi_testparticles ), "number of testparticles for J/psi (in addition to the number employed for the rest of added particles)")

    ("heavy_quark.hqCorrelationsOutput", po::value<bool>( &hqCorrelationsOutput )->default_value( hqCorrelationsOutput ), "whether correlation analysis of heavy quark pairs is done")
    ;

    // Add some miscellaneous options
    misc_options.add_options()
    ("misc.repeat_timesteps", po::value<bool>( &switch_repeatTimesteps )->default_value( switch_repeatTimesteps ), "repeat timesteps in cases where the probability has been > 1" )
//   ("misc.interpolation_border", po::value<double>( &interpolationBorder )->default_value( interpolationBorder ), "X where interpolation of MFP is done for E > X*T")
    ("misc.jet_mfp_computation", po::value<int>()->default_value( jetMfpComputationSwitch ), "treatment for the mean free path of added particles ( 0 = computeMfpLastTimestep, 1 = computeMfpIteration, 2 = computeMfpInterpolation, 3 = fixedMfp, 4 = thermalMfpGluon)")
    ("misc.fixed_mfp_added", po::value<double>( &fixed_mfp_added )->default_value( fixed_mfp_added ), "Mean free path of added particles set by hand. Does not depend on energy of particle" )
    ("misc.mfpAddedRangeVariation", po::value<double>( &mfpAddedRangeVariation )->default_value( mfpAddedRangeVariation ), "Range in % in respect to the old mean free path, in which the new value of the mean free path is expected to be" )
    ("misc.analysisPhotonsTimeCut", po::value<double>( &analysisPhotonsTimeCut )->default_value( analysisPhotonsTimeCut ), "in fm/c. At what time the photon production starts. Default 0.0 fm/c." )
    ("misc.analysisPhotonsPTCut", po::value<double>( &analysisPhotonsPTCut )->default_value( analysisPhotonsPTCut ), "in GeV. Above what pT the photon production is counted. Default 0.0 GeV ." )   
    ("misc.restrictParentPTForPhotonproduction", po::value<bool>( &restrictParentPTForPhotonproduction )->default_value( restrictParentPTForPhotonproduction ), "if true, restrict the highest parent PT for Photon-Production to lie between minAllowedParentPT and maxAllowedParentPT." )
    ("misc.minAllowedParentPT", po::value<double>( &minAllowedParentPT )->default_value( minAllowedParentPT ), "in GeV. Highest PT of the two parents must be greater than this parameter." )     
    ("misc.maxAllowedParentPT", po::value<double>( &maxAllowedParentPT )->default_value( maxAllowedParentPT ), "in GeV. Highest PT of the two parents must be smaller than this parameter." )
    ("misc.mfpCellAveraging", po::value<bool>( &mfpCellAveraging )->default_value( mfpCellAveraging ), "If the specific mean free path-calculation for 23-photonproduction should be averaged over 5 cells." )
    ("misc.analysisPhotonsNBinsV2", po::value<int>( &analysisPhotonsNBinsV2 )->default_value( analysisPhotonsNBinsV2 ), "Number of Bins for the V2-Pt Analysis.")
    ("misc.usedExternalField", po::value<double>( &usedExternalField )->default_value( usedExternalField ), "Externel Bfield used for online runs. Given as eB/mpi^2." )
;


    // Group offline reconstruction options
    offline_options.add_options()
    ("offline.name", po::value<string>( &originalName )->default_value( originalName ), "name of the original BAMPS run that is to be reconstructed")
    ("offline.offline_data_dir", po::value<string>( &pathdirOfflineData )->default_value( pathdirOfflineData ), "directory from which the \"offline\" data is read")
    ("offline.use_fixed_dt", po::value<bool>( &switch_fixed_dt )->default_value( switch_fixed_dt ), "Indicates whether a fixed dt (provided via fixed_dt) should be used. Use time steps from the original run if not" )
    ("offline.fixed_dt", po::value<double>( &fixed_dt )->default_value( fixed_dt ), "fixed dt (time steps) at which the reconstructed medium is \"sampled\" [optional]")
    ("offline.factor_dt", po::value<double>( &factor_dt )->default_value( factor_dt ), "factor with which time steps from the original run should be scaled for use in sampling of the reconstructed medium (should be <1)")
    ("offline.nAdded", po::value<int>( &numberOfParticlesToAdd)->default_value( numberOfParticlesToAdd ), "number of (high-pt) particles that is added on top of the reconstructed medium, using it as a background. Do not use together with nAddedEvents." )
    ("offline.nAddedEvents", po::value<int>( &numberOfAddedEvents)->default_value( numberOfAddedEvents ), "number of heavy ion collision events, set on top of the offline reconstruction. Do not use together with nAdded. This input is needed if Pythia data files are read in to normalize the total yield of the particles. The number of one such heavy ion collision event includes < number of produced particles in pp > * Ntest * Nbin" )
    ("offline.minPT_added", po::value<double>( &minimumPT )->default_value( minimumPT ), "minimum p_T [GeV] of the added particles. If p_T of added particle falls below this value it does not scatter anymore.")
    ("offline.onlyMediumEvolution", po::value<bool>( &onlyMediumEvolution )->default_value( onlyMediumEvolution ), "whether to do only medium evolution without scatterings or rate and debye mass calculations.")

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
    if( ( scatt_amongAddedParticles && Particle::N_psi_states == 0 ) || ( Particle::N_psi_states > 0 && !scatt_amongAddedParticles ) )
    {
        string errMsg = "Scatterings among added particles and Jpsi not active or vice versa.";
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

    if( ( studyNonPromptJpsiInsteadOfElectrons && ( outputScheme != cms_hq_nonPromptJpsi && outputScheme != alice_hq_nonPromptJpsi ) ) || ( ( outputScheme == cms_hq_nonPromptJpsi || outputScheme == alice_hq_nonPromptJpsi ) && !studyNonPromptJpsiInsteadOfElectrons ) )
    {
        string errMsg = "Study Jpsi instead of electrons but no output for this or vice versa.";
        throw eConfig_error( errMsg );
    }

    if( initialStateType == onlyJpsiInitialState && Particle::N_psi_states == 0 )
    {
        string errMsg = "Only Jpsi in initial state, but N_psi_states = 0.";
        throw eConfig_error( errMsg );
    }

    // check if misc.fixed_mfp_added is set in input file: if not default value or no value given (latter happens if no input file is given at all)
    if ( !( vm["misc.fixed_mfp_added"].defaulted() || vm["misc.fixed_mfp_added"].empty() ) && jetMfpComputationSwitch != fixedMfp )
    {
        string errMsg = "Option fixed_mfp_added can only be set if jetMfpComputationSwitch is set to fixedMfp.";
        throw eConfig_error( errMsg );
    }

    if ( N_heavy_flavors_input > 0 && !jet_tagged )
    {
        string errMsg = "If heavy quarks are involved, jets must be tagged.";
        throw eConfig_error( errMsg );
    }

    if ( scatt_furtherOfflineParticles && !jet_tagged )
    {
        string errMsg = "If recoiled particles are further evolved, jets must be tagged.";
        throw eConfig_error( errMsg );
    }
}

void config::processHeavyQuarkOptions()
{
    configBase::processHeavyQuarkOptions();

    // see if running coupling for all particles and running coupling for heavy quarks are both set simultaneously
    if( ( !( vm["simulation.running_coupling"].defaulted() || vm["simulation.running_coupling"].empty() ) ) &&
            ( !( vm["heavy_quark.running_coupling"].defaulted() || vm["heavy_quark.running_coupling"].empty() ) ) )
    {
        if( couplingRunning != couplingRunningHeavyQuarksInput )
        {
            string errMsg = "Different coupling scheme for heavy and light particles is not implemented.";
            throw eConfig_error( errMsg );
        }
    }
    // if it is only set for heavy quarks and no light quarks are present, set it globally.
    else if ( !( vm["heavy_quark.running_coupling"].defaulted() || vm["heavy_quark.running_coupling"].empty() ) )
    {
        if( N_light_flavors_input < 0 )
        {
            cout << "Although the coupling scheme has just been set for heavy quarks, it is adapted for all particles since there are no other particles, N_f_light = " << N_light_flavors_input << "." << endl;
            couplingRunning = couplingRunningHeavyQuarksInput;
            coupling::set_isRunning( couplingRunning );
        }
        else
        {
            string errMsg = "Heavy flavor coupling is set, but also light partons are present.";
            throw eConfig_error( errMsg );
        }
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
    printOptionsDescriptionToIniFormat( parameters23, output );
    printOptionsDescriptionToIniFormat( misc_options, output );
    if( Particle::N_heavy_flavor > 0 )
        printOptionsDescriptionToIniFormat( heavy_quark_options, output );
    printOptionsDescriptionToIniFormat( offline_options, output );
    printOptionsDescriptionToIniFormat( crossSection_options, output);
}


void config::readAndPrepareInitialSettings( offlineOutputInterface* const offlineInterface )
{
    double freezeOutEnergyDensity_offline;

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
    freezeOutEnergyDensity_offline = ptrSimulationData->freezeOutEnergyDensity;
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

    if( freezeOutEnergyDensity_offline > freezeOutEnergyDensity )
    {
        string errMsg = "Freeze-out energy density of offline data is larger than the input value of the simulation.";
        throw eConfig_error( errMsg );
    }

    if( numberOfParticlesToAdd >= 0 )
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
    int maxID = 0;
    for ( int i = 0; i < ptrInitialParticles->particleVector->size(); i++ )
    {
        particles_init[i].unique_id = (*(ptrInitialParticles->particleVector))[i].unique_id;
        if( particles_init[i].unique_id > maxID )
            maxID = particles_init[i].unique_id;
        particles_init[i].FLAVOR = (*(ptrInitialParticles->particleVector))[i].FLAVOR;
        particles_init[i].m = (*(ptrInitialParticles->particleVector))[i].m;
        particles_init[i].Pos = (*(ptrInitialParticles->particleVector))[i].Pos;
        particles_init[i].Mom = (*(ptrInitialParticles->particleVector))[i].Mom;
        particles_init[i].md2g = (*(ptrInitialParticles->particleVector))[i].md2g;
        particles_init[i].md2q = (*(ptrInitialParticles->particleVector))[i].md2q;
        particles_init[i].isAlreadyInAddedParticles.resize( static_cast< int >( numberOfParticlesToAdd / 2 ), false );
    }
    Particle::unique_id_counter = maxID + 1;

    for ( int i = 0; i < particles_init.size(); i++ )
    {
        particlesEvolving[i] = particles_init[i];

        particlesEvolving[i].init = true;
        particlesEvolving[i].free = false;

        if ( particlesEvolving[i].Pos.T() < timefirst )
        {
            particlesEvolving[i].edge = true;
            particlesEvolving[i].init = false;
        }

        particlesEvolving[i].Old  = particles_init[i].Mom;
    }

    particles_atTimeNow = particles_init;
}

