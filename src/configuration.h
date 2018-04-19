//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


/** @file
 * @brief This file provides global objects and a configuration interface
 */


#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <string>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <stdint.h>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "configurationbase.h"
#include "particleOffline.h"
#include "offlineoutput.h"
#include "interpolation_iniJpsi.h"
#include "initialmodel.h"


/** @brief Enumeration type for PDF sources */
enum PDF_SOURCE_TYPE { builtInGRV, LHAPDF };

/** @brief Enumeration type for different variants of computing the mean free path of added particles */
enum JET_MFP_COMPUTATION_TYPE { computeMfpLastTimestep, computeMfpIteration, computeMfpInterpolation, fixedMfp, thermalMfpGluon };

/** @brief Enumeration type for different output schemes to decide which kind of output is printed 
 * set output studies according to the output scheme here in analysis::handle_output_studies( OUTPUT_SCHEME _outputScheme )
*/
enum OUTPUT_SCHEME { no_output = 0,  
light_parton_phenix = 10,
light_parton_lhc = 20,
// heavy quarks:
phenix_hq_electrons = 101,
alice_hq_electrons = 111,
alice_hq_muons = 112,
alice_hq_dmesons = 113,
cms_hq_nonPromptJpsi = 114,
alice_hq_nonPromptJpsi = 115,
phenix_jpsi = 131,
alice_jpsi = 141,
cms_jpsi = 151,
hq_corr = 181,
// jets:
background_jets = 161,
central_densities = 171
};


/** @brief This extends the namespace ns_casc (see globalsettings.h) to hold a global particle vector */
namespace ns_casc
{  
  /** @brief Vector to hold all particles in the simulation, GLOBAL 
   * 
   * Declared as extern here. There MUST be ONE definition in a cpp-file that includes this header in the global 
   * scope (it's in configuration.cpp).
   * All files that include this header can then access the vector.
   */
  extern std::vector<ParticleOffline> particlesEvolving;
  extern std::vector<ParticleOffline> particles_init;
  extern std::vector<ParticleOffline> particles_atTimeNow;
  extern std::vector<ParticleOffline> particles_atTimeNowCopy;
  extern std::vector<ParticleOffline> addedParticles;
  extern std::vector<ParticleOffline> addedParticlesCopy;
//   extern std::vector<ParticleHFelectron> addedPartcl_electron;
  extern std::vector<ParticleOffline> addedPartcl_electron;
  extern std::vector<ParticleOffline> scatteredMediumParticles;

}
//--------------------------------------------------------//

using std::string;


extern double dt;


/**
 * @brief interface for runtime configuration
 *
 * This class provides a basic interface for settings made by the user at runtime. It is based on boost::program_options 
 * and can handle command line arguments as well as INI-type configuration files.
 * 
 * It derives from and extends the base class configBase that is provided in the BAMPS lib. configBase provides the basic
 * functionality, some helper routines (e.g. printing of used parameters) and also some basic boost::program_options::options_description groups
 * and options. The structure of configBase allows this derived class to extend the functionality by 
 *  a) adding to existing boost::program_options::options_description objects
 *  b) defining new boost::program_options::options_description objects and adding those to either the command line or the 
 *    configuration file option groups
 * 
 * This is done in config::initializeProgramOptions() and in config::groupProgramOptions(). In principle these two steps 
 * would be sufficient as they eventually only add to the groups configBase::command_line_options and configBase::config_file_options
 * which are members of the base class. Therefore reading and processing - including the extended options - could then be
 * handled by the base class.
 * However, in order to provide the flexibility to add more sophisticated checking, processing, etc. to the extended options
 * this class explicitly re-implements the processing routines of the base class configBase. Mostly by simply calling 
 * the appropriate member routines of the base class. 
 * 
 * The program flow is therefore identical to the one descriped in the documentation of configBase
 */
class config : public configBase
{
 public:
  /** @brief Constructor that internally reads the provided input file */
  config();
  /** @brief Standard constructor */
  ~config() {};
  
  /** @brief processes command line arguments and settings provided via a configuration file */
  void readAndProcessProgramOptions(const int argc, char* argv[]);
  
  /** ---- simulation parameters ---- */ 
  /** @brief Interface for config::freezeOutEnergyDensity */
  double getFreezeOutEnergyDensity() const { return freezeOutEnergyDensity; }
  
  /** @brief Interface for config::scatt_offlineWithAddedParticles */
  bool isScatt_offlineWithAddedParticles() const {return scatt_offlineWithAddedParticles;}
  
  /** @brief Interface for config::scatt_amongOfflineParticles */
  bool isScatt_amongOfflineParticles() const {return scatt_amongOfflineParticles;}
  
  /** @brief Interface for config::scatt_amongAddedParticles */
  bool isScatt_amongAddedParticles() const {return scatt_amongAddedParticles;}

  /** @brief Interface for config::scatt_furtherOfflineParticles */
  bool isScatt_furtherOfflineParticles() const {return scatt_furtherOfflineParticles;}

  /** @brief Interface for config::N_light_flavors_input */
  int getNlightFlavorsAdded() const {return N_light_flavors_input;}
  
  /** @brief Interface for config::N_heavy_flavors_input */
  int getNheavyFlavorsAdded() const {return N_heavy_flavors_input;}
  
  /** @brief Interface for config::jet_tagged */
  bool isJetTagged() const {return jet_tagged;}
  /** ------------------------------- */

  /** ---- initial state options ---- */ 
  /** @brief Interface for config::initialStateType */
  INITIAL_STATE_TYPE getInitialStateType() const { return initialStateType; }
  
  /** @brief Interface for config::PDFsource */
  PDF_SOURCE_TYPE getPDFsource() const { return PDFsource; }
  
  /** @brief Interface for config::LHAPDFdatasetName */
  string getLHAPDFdatasetName() const { return LHAPDFdatasetName; }
  
  /** @brief Interface for config::LHAPDFmember */
  unsigned short int getLHAPDFmember() const { return LHAPDFmember; }
  
  /** @brief Interface for config::LHAPDFuseGrid */
  bool getLHAPDFuseGrid() const { return LHAPDFuseGrid; }
  
  /** @brief Interface for config::nuclearPDFs */
  bool useNuclearPDFs() const { return nuclearPDFs; }
  
  /** @brief Interface for config::nuclearPDFdatasetName */
  string getNuclearPDFdatasetName() const { return nuclearPDFdatasetName; }
  
  /** @brief Interface for config::pythiaParticleFile */
  string getPythiaParticleFile() const { return pythiaParticleFile; }

  /** @brief Interface for config::cgcParticleFile */
  string getCgcParticleFile() const { return cgcParticleFile; }
  
  /** @brief Interface for config::mcatnloParticleFileCharm */
  string getMcatnloParticleFile() const {return mcatnloParticleFile;}
  
  /** @brief Interface for config::P0 */
  double getPtCutoff() const { return P0; }
  
  /** @brief Interface for config::initialPartonFlavor */
  FLAVOR_TYPE getInitialPartonFlavor() const { return initialPartonFlavor; };

  /** ------------------------------- */

  /** -------- output options ------- */   
  /** @brief Interface for config::outputSwitch_progressLog */
  bool doOutput_progressLog() const { return outputSwitch_progressLog; }  
  
  /** @brief Interface for config::outputSwitch_detailedParticleOutput */
  bool doOutput_detailedParticleOutput() const { return outputSwitch_detailedParticleOutput; }
  
  /** @brief Interface for config::outputSwitch_movieOutput */
  bool doOutput_movieOutputJets() const {return outputSwitch_movieOutputJets;}
  
  /** @brief Interface for config::outputSwitch_movieOutput */
  bool doOutput_movieOutputBackground() const {return outputSwitch_movieOutputBackground;}

  /** @brief Interface for config::outputSwitch_scatteredMediumParticlesOutput */
  bool doOutput_scatteredMediumParticles() const {return outputSwitch_scatteredMediumParticlesOutput;}

  /** @brief Interface for config::v2RAAoutput */
  bool isV2RAAoutput() const {return  v2RAAoutput;}
  
  /** @brief Interface for config::v2RAAoutputIntermediateSteps */
  bool isV2RAAoutputIntermediateSteps() const {return v2RAAoutputIntermediateSteps;}
  
  /** @brief Interface for config::dndyOutput */
  bool isDndyOutput() const {return  dndyOutput;}
  
  /** @brief Interface for config::outputScheme */
  OUTPUT_SCHEME getOutputScheme() const {return  outputScheme;}
  /** ------------------------------- */
  
  /** -------- heavy quark options ------- */ 
  // hadronization of open heavy flavor
  /** @brief Interface for config::hadronization_hq */
  bool isHadronizationHQ() const {return hadronization_hq;}
  
  /** @brief Interface for config::mesonDecay */
  bool isMesonDecay() const {return mesonDecay;}
  
  /** @brief Interface for config::numberElectronStat */
  int getNumberElectronStat() const {return  numberElectronStat;}
  
  /** @brief Interface for config::muonsInsteadOfElectrons */
  bool isMuonsInsteadOfElectrons() const {return muonsInsteadOfElectrons;}
  
  /** @brief Interface for config::studyNonPromptJpsiInsteadOfElectrons */
  bool isStudyNonPromptJpsiInsteadOfElectrons() const {return studyNonPromptJpsiInsteadOfElectrons;}
  
  
  // Jpsi
  /** @brief Interface for config::isotropicCrossSecJpsi */
  bool isIsotropicCrossSecJpsi() const {return  isotropicCrossSecJpsi;}
  
  /** @brief Interface for config::constantCrossSecJpsi */
  bool isConstantCrossSecJpsi() const {return  constantCrossSecJpsi;}
  
  /** @brief Interface for config::constantCrossSecValueJpsi */
  double getConstantCrossSecValueJpsi() const {return  constantCrossSecValueJpsi; }
  
  /** @brief Interface for config::TdJpsi */
  double getTdJpsi() const {return  TdJpsi;}
  
  /** @brief Interface for config::jpsi_sigmaAbs */
  double getSigmaAbs() const {return jpsi_sigmaAbs;} 
  
  /** @brief Interface for config::jpsi_agN */
  double getJpsiagN() const {return jpsi_agN;} 
  
  /** @brief Interface for config::shadowing_model */
  shadowModelJpsi getShadowingModel() const  {return shadowing_model;} 

  /** @brief Interface for config::jpsi_formationTime */
  double getJpsiFormationTime() const  {return jpsi_formationTime;} 
  
  /** @brief Interface for config::jpsi_testparticles */
  int getJpsiTestparticles() const  {return jpsi_testparticles;} 
  
  
  // heavy quark output
  /** @brief Interface for config::hqCorrelationsOutput */
  bool isHqCorrelationsOutput() const {return hqCorrelationsOutput;}
  /** ------------------------------------ */

  /** -------- miscellaneous options ------- */ 
  bool repeatTimesteps() const { return switch_repeatTimesteps; }

  /** @brief Interface for config::jetMfpComputationSwitch */
  JET_MFP_COMPUTATION_TYPE getJetMfpComputationType() const {return jetMfpComputationSwitch;}
  
  /** @brief Interface for config::mfpAddedRangeVariation */
  double getMfpAddedRangeVariation() const {return mfpAddedRangeVariation;}
  
  /** @brief Interface for config::mfp_added */
  double getFixedMfpAdded() const {return fixed_mfp_added;}
  
//   /** @brief Interface for config::interpolationBorder */
//   double getMFPInterpolationBorder() const {return interpolationBorder;}
//   
  /** @brief Interface for config::tau_min */
  double getTauMin() const { return tau_min; }

  /** @brief Interface for config::minNumberForTemperature */
  double getMinNumberForTemperature() const { return minNumberForTemperature; }

  /** ------------------------------------ */
  
  /** -------- offline reconstruction options ------- */ 
  /** @brief Interface for config::originalName */
  string getOriginalName() const {return originalName;} 
  
  /** @brief Interface for config::pathdirOfflineData */
  string getPathdirOfflineData() const { return pathdirOfflineData; }
  
  /** @brief Interface for config::pathdirOfflineData, returns pointer to char */
  const char* getPathdirOfflineDataChar() const { return pathdirOfflineData.c_str(); }
  
  /** @brief Interface for config::fixed_dt */
  double getFixed_dt() const {return fixed_dt;}
  /** @brief Interface for config::switch_fixed_dt */
  bool useFixed_dt() const {return switch_fixed_dt;}
  /** @brief Interface for config::factor_dt */
  double getFactor_dt() const { return factor_dt; }
  
  /** @brief Interface for config::numberOfParticlesToAdd */
  int getNumberOfParticlesToAdd() const { return numberOfParticlesToAdd; }
  /** @brief Interface for config::numberOfAddedEvents */
  int getNaddedEvents() const { return numberOfAddedEvents; }
  /** @brief Interface for config::minimumPT */
  double getMinimumPT() const { return minimumPT; }
  
  /** @brief Interface for config::ringNumber */
  int getRingNumber() const { return ringNumber; }
  /** @brief Interface for config::centralRingRadius */
  double getCentralRingRadius() const { return centralRingRadius; }
  /** @brief Interface for config::deltaR */
  double getDeltaR() const { return deltaR; }
  
  /** @brief Interface for config::N_init */
  int getN_init() const { return N_init; }
  /** @brief Interface for config::dx */
  double get_dx() const { return dx; }
  /** @brief Interface for config::dy */
  double get_dy() const { return dy; }
  /** @brief Interface for config::transLen */
  double getTransLen() const { return transLen; }
  
  /** @brief Interface for config::timeshift */
  double getTimeshift() const {return timeshift;}
  /** @brief Interface for config::timefirst */
  double getTimefirst() const {return timefirst;}
  
  // collision parameters
  /** @brief Interface for config::A */
  double getA() const { return A; }
  /** @brief Interface for config::Aatomic */
  double getAatomic() const { return Aatomic; }
  /** @brief Interface for config::B */
  double getB() const { return B; }
  /** @brief Interface for config::Batomic */
  double getBatomic() const { return Batomic; }
  /** @brief Interface for config::sqrtS */
  double getSqrtS() const { return sqrtS; }
  /** @brief Interface for config::impact */
  double getImpactParameter() const { return impact; }  
  /** ----------------------------------------------- */ 
  
  /** ----- auxiliary routines ----- */
  void readAndPrepareInitialSettings( offlineOutputInterface*const _offlineInterface );
  /** ------------------------------ */

  /** Set seed chosen during runtime */
  void setSeed ( uint32_t _seed ){ seed = _seed; };
  
 protected:
   /** ----- auxiliary routines ----- */
   /** @brief Sort the options into groups */   
   void groupProgramOptions();
   
   /** @brief Actually define all the options */
   void initializeProgramOptions();
   
   /** @brief Read command line arguments and settings provided via a configuration file and via the commmand line */
   void readProgramOptions(const int argc, char* argv[]);
   
   /** @brief Do some processing on program options that have previously been read via configPrototype::readProgramOptions */
   void processProgramOptions();
   
   /** @brief Print a complete configuration file using all current parameter values */
   void printUsedConfigurationParameters();   
   
   /** @brief Do some checks on user-provided options and parameters */
   void checkOptionsForSanity();
   
   /** @brief Some processing of heavy quark options */
   void processHeavyQuarkOptions();
   /** ------------------------------ */
   
   /** ------ boost::program_options objects ------- */ 
   // base class provides:
   //    po::options_description command_line_options;
   //    po::options_description config_file_options;
   //    po::options_description visible_options;
   //    po::positional_options_description pos_options;
   //    po::variables_map vm;
   //    po::options_description usage_information;
   //    po::options_description hidden_options;
   // 
   //    po::options_description general_options;
   //    po::options_description simulation_parameters;
   //    po::options_description output_options;
   //    po::options_description misc_options;
   //    po::options_description heavy_quark_options;
   
   po::options_description initial_state_options;   
   po::options_description offline_options;
   /** ------ boost::program_options objects ------- */ 

   
  /** ------ general options ------- */  
  // base class provides:
  // string jobname;
  // long seed;
  /** ------------------------------ */
  
  /** ---- simulation parameters ---- */ 
  // base class provides:
  // double runtime  
  // int testparticles
  // int N_light_flavors_input
  // int N_heavy_flavors_input
  // bool couplingRunning
  
  /** @brief Energy density for freeze out (in GeV/fm^3)
   * Particles in regions with energy densities below this threshold will be freely streaming
   */
  double freezeOutEnergyDensity;
  
  /** @brief Whether offline particles are allowed to scatter with added particles */
  bool scatt_offlineWithAddedParticles;
  
  /** @brief Whether offline particles are allowed to scatter with other offline particles */
  bool scatt_amongOfflineParticles;
  
  /** @brief Whether added particles are allowed to scatter with other added particles */
  bool scatt_amongAddedParticles;
  
  /** @brief Whether the added particles are treated as tagges jets */
  bool jet_tagged;
  /** ------------------------------- */

  /** @brief Whether scattered offline particles are allowed to scatter again with other added particles */
  bool scatt_furtherOfflineParticles;
  /** ------------------------------- */

  /** ---- initial state options ---- */ 
  /** @brief Which type of initial state to use */
  INITIAL_STATE_TYPE initialStateType;
  
  /** @brief Which source to use for the parton distribution functions
   * 0 = built-in GRV parametrization
   * 1 = PDF parametrizations provided by the LHAPDF library
   */
  PDF_SOURCE_TYPE PDFsource;
  
  /** @brief Name of the LHAPDF data set that should be used */
  string LHAPDFdatasetName;
  
  /** @brief Which member of the given LHAPDF data set should be used */
  unsigned short int LHAPDFmember;
  
  /** @brief Whether to use the grid version of the given LHAPDF set */
  bool LHAPDFuseGrid;
  
  /** @brief Whether to use nPDFs (only available in combination with LHAPDF and minijets) */
  bool nuclearPDFs;
  
  /** @brief Name of the nPDF set that is to be used */
  string nuclearPDFdatasetName;
                     
  /** @brief Relative or full path (including filename) of PYTHIA output file with particle data
   * Declare as "-" if particle momenta should be sampled via glauber method
   */
  string pythiaParticleFile;
  
  /** @brief Relative or full path (including filename) of color glass condensate output file with particle data
   */
  string cgcParticleFile; 
  
  /** @brief Relative or full path (including filename) of MC@NLO output file with particle data
   */
  string mcatnloParticleFile;
  
  /** @brief Lower PT-cutoff [GeV] used for spectrum of initial conditions */
  double P0;  
  
  /** @brief Flavor of initial parton pair */
  FLAVOR_TYPE initialPartonFlavor;
  
  /** ------------------------------- */
  
  /** -------- output options ------- */ 
  // provided by base class:
  // string standardOutputDirectoryName
  
  /** @brief Specify whether progress output should be written to a file */
  bool outputSwitch_progressLog;
  
  /** @brief Specify whether detailed particle output should be written */
  bool outputSwitch_detailedParticleOutput;
  
  /** @brief Specify whether movie output should be written (for added jet particles) */
  bool outputSwitch_movieOutputJets;
  
  /** @brief Specify whether movie output should be written (for the reconstructed background) */
  bool outputSwitch_movieOutputBackground;

  /** @brief Specify whether scattered medium particles output should be written */
  bool outputSwitch_scatteredMediumParticlesOutput;
  
  /** @brief Whether v2 and RAA output are printed */
  bool v2RAAoutput;
  
  /** @brief Whether v2 and RAA output are printed at each analyisis time step (otherwise just at beginning and end) */
  bool v2RAAoutputIntermediateSteps;
  
  /** @brief Whether dndy output is written out */
  bool dndyOutput;
  
  /** @brief Output schemes to decide which kind of output is printed */
  OUTPUT_SCHEME outputScheme;
  /** ------------------------------- */ 

  /** -------- heavy quark options ------- */ 
  // provided by base class:
  // open heavy flavor
  //   double KggQQbar
  //   double KgQgQ
  //   double kappa_gQgQ
  //   bool isotropicCrossSecGQ
  //   bool constantCrossSecGQ
  //   double constantCrossSecValueGQ
  //   double Mcharm_input
  //   double Mbottom_input
  
  /** @brief Whether a running coupling is employed for heavy quark processes. Do not use here. Instead set running coupling also for light particles under misc. This option here is only to be compatible with older inputfiles and does nothing else to set the coupling globally. */
  bool couplingRunningHeavyQuarksInput;
  
  // hadronization of open heavy flavor
  /** @brief Whether hadronization of heavy quarks to D and B mesons is carried out */
  bool hadronization_hq;
  
  /** @brief Whether decay of heavy mesons from charm and bottom to electrons is performed */
  bool mesonDecay;
  
  /** @brief To improve statistics one meson decays to numberElectronStat electrons */
  int numberElectronStat;
  
  /** @brief Whether decay to muons should be performed instead of electrons */
  bool muonsInsteadOfElectrons;
  
  /** @brief Whether decay of B mesons to non prompt Jpsi should be performed instead to electrons */
  bool studyNonPromptJpsiInsteadOfElectrons;
  
  
  // Jpsi
  /** @brief How many psi states are allowed */
  int N_psi_input;
  
  /** @brief Whether an isotropic momentum sampling is employed for process Q + Qbar -> g + J/psi */
  bool isotropicCrossSecJpsi;
  
  /** @brief Whether a constant cross section is employed for process Q + Qbar -> g + J/psi */
  bool constantCrossSecJpsi;
  
  /** @brief Value of constant cross section for process Q + Qbar -> g + J/psi in mb */
  double constantCrossSecValueJpsi;
  
  /** @brief Mass of J/psi */
  double Mjpsi_input;

  /** @brief Dissociation temperature of J/psi */
  double TdJpsi;
  
  /** @brief Absorption cross section for initial J/psi in mb */
  double jpsi_sigmaAbs;
  
  /** @brief Parameter for momentum broadening of initial J/psi */
  double jpsi_agN;
  
  /** @brief Shadowing model used for initial J/psi */
  shadowModelJpsi shadowing_model;

  /** @brief Formation time for initial J/psi in addition to the standard 1/M_T */
  double jpsi_formationTime;
  
  /** @brief Number of testparticles for J/psi (in addition to the number employed for the rest of added particles) */
  int jpsi_testparticles;
  
  
  
  // heavy quark output
  /** @brief Whether correlation analysis of heavy quark pairs is done */
  bool hqCorrelationsOutput;
  /** ------------------------------------ */
  
  /** -------- miscellaneous options ------- */ 
  // provided by base class:

  /** @brief How to compute the mean free path high energy particles?
   *
   * 0 = computeMfpLastTimestep = use mean free path from last time step if available
   * 1 = computeMfpIteration = iterative computation every time step for every jet
   * 2 = computeMfpInterpolation = interpolate mean free path from tables
   * 3 = fixedMfp = set fixed mean free path by hand
   * 4 = thermalMfpGluon = set mean free path of emitted gluon which is assumed to be thermal
   */
  JET_MFP_COMPUTATION_TYPE jetMfpComputationSwitch;
  
  /** @brief Range in % in respect to the old mean free path, in which the new value of the mean free path is expected to be. 
   * For new value the region [oldvalue * (1 - mfpAddedRangeVariation/100) , oldvalue * (1 + mfpAddedRangeVariation/100) ] is tested. Used in iterate_mfp_bisection(),  
   */
  double mfpAddedRangeVariation;

  /** @brief Fixed mean free path which is used if jetMfpComputationSwitch == fixedMfp.  */
  double fixed_mfp_added; // fm
  
//   /** @brief X where interpolation of MFP is done for E > X*T */
//   double interpolationBorder; 
  
  /** @brief Whether timesteps are repeated in cases where the probability has been > 1 */
  bool switch_repeatTimesteps;

  /** @brief Proper time after which partons are allowed to scatter / radiate */
  double tau_min;

  /** @brief Minimum number necessary for calculating temperature */
  double minNumberForTemperature;
  
  /** ------------------------------- */

  /** -------- offline reconstruction options ------- */ 
  /** @brief Path to which the output needed for offline analysis is written */
  string pathdirOfflineData;
  
  /** @brief The name of the BAMPS job that is processed (input file names) */
  string originalName;
  
  /** @brief Initial seed used in the full BAMPS run */
  uint32_t originalSeed;
  
  /** @brief A fixed dt (time steps) at which the reconstructed medium is "sampled" [optional] */
  double fixed_dt;
  /** @brief Indicates whether a fixed dt (provided via fixed_dt) should be used. Use time steps from the original run if not. */
  bool switch_fixed_dt;
  /** @brief Factor with which time steps from the original run should be scaled for use in sampling of the reconstructed medium (should be <1) */
  double factor_dt;
  
  /** @brief Number of (high-pt) particles that is added on top of the reconstructed medium, using it as a background. This is needed for mini-jet initial conditions. However, to do a correct normalisation for particles yields also set numberOfAddedEvents */
  int numberOfParticlesToAdd;
  /** @brief Number of heavy ion collision events, set on top of the offline reconstruction. 
   *  This input is needed if Pythia data files are read in to normalize the total yield of the particles. 
   * The number of particles in one such heavy ion collision event is equal to < number of produced particles 
   * in pp > * Ntest * Nbin. Consequently, 1 would mean that one adds as many particles as there are in a heavy 
   * ion collision times Ntest, making it analogously to a standard BAMPS simulation.
   * 
   * To use the number of events to sample the initial number of particles instead of the given number of particles via numberOfParticlesToAdd set numberOfParticlesToAdd to a negative value.
  */
  int numberOfAddedEvents;
  /** @brief Minimum p_T [GeV] of the added particles */
  double minimumPT;
  
  
  // the following parameters are read at runtime from the offline data recorded by the original run, 
  // see config::readAndPrepareInitialSettings
  // testparticles is defined in configBase but overwritten in readAndPrepareInitialSettings.
  /** @brief centralRingRadius as used in the original run */
  double centralRingRadius;
  /** @brief deltaR as used in the original run */
  double deltaR;
  /** @brief ringnumber as used in the original run */
  int ringNumber;
  /** @brief dx as used in the original run */
  double dx;
  /** @brief dy as used in the original run */
  double dy;
  /** @brief transversal length as used in the original run */
  double transLen;
  /** @brief initial number of particles as used in the original run */
  int N_init;
  /** @brief first time step as used in the original run */
  double timefirst;
  /** @brief "time shift" as used in the original run */
  double timeshift;
  /** @brief number of active light flavors in the offline run */
  int N_light_flavors_offline;
  /** @brief number of active heavy flavors in the offline run */
  int N_heavy_flavors_offline;
  /** @brief Mass number of nucleus A  */
  double A;  
  /** @brief Atomic number, i.e. number of protons, of nucleus A */
  double Aatomic;
  /** @brief Mass number of nucleus B  */
  double B;
  /** @brief Atomic number, i.e. number of protons, of nucleus B */
  double Batomic;
  /** @brief Center of mass energy per NN pair in GeV */
  double sqrtS;   
  /** @brief Impact parameter in fm */
  double impact;  
  /** ----------------------------------------------- */
};

#endif
