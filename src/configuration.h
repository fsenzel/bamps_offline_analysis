//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------
//---------------------------------------------


/** @file
 * @brief This file provides global objects and a configuration interface
 */


#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <string>
#include <vector>
#include <stdexcept>
#include <stdint.h>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "globalsettings.h"
#include "particle.h"
#include "offlineoutput.h"
#include "interpolation_iniJpsi.h"

#define HELP_MESSAGE_REQUESTED 1


/** @brief Enumeration type for possible initial state models */
enum INITIAL_STATE_TYPE { miniJetsInitialState, pythiaInitialState, cgcInitialState, mcatnloInitialState };

/** @brief Enumeration type for different variants of computing the mean free path of high-pt particles */
enum JET_MFP_COMPUTATION_TYPE { computeMfpDefault, computeMfpIteration, computeMfpInterpolation };

/** @brief Enumeration type for different output schemes to decide which kind of output is printed 
 * set output studies according to the output scheme here in analysis::handle_output_studies( OUTPUT_SCHEME _outputScheme )
*/
enum OUTPUT_SCHEME { no_output = 0, 
// heavy quarks:
phenix_hq_electrons = 101,
alice_hq_electrons = 111,
alice_hq_muons = 112,
alice_hq_dmesons = 113,
cms_hq_nonPromptJpsi = 114,
phenix_jpsi = 131
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
}
//--------------------------------------------------------//

using std::string;


extern double dt;


/**
 * @brief Interface for settings made at runtime
 *
 * This class provides an interface for settings made by the user at runtime. These need to be specified in a configuration
 * file whose name is passed to the program as the first (and only) command line parameter. This file is read when the constructor
 * of the class config is called.
 */
class config
{
 public:
  /** @brief Constructor that internally reads the provided input file */
  config(const int argc, char* argv[]);
  /** @brief Standard constructor */
  ~config() {};

  /** @brief processes command line arguments and settings provided via a configuration file */
  void readAndProcessProgramOptions(const int argc, char* argv[]);
  
  /** ------ general options ------- */  
  /** @brief Interface for config::jobname */
  string getJobName() const { return jobName; }
  
  /** @brief Interface for config::seed */
  long getSeed() const { return seed; }
  /** ------------------------------- */
  
  /** ---- collision parameters ---- */ 
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
  /** ------------------------------- */
  
  /** ---- simulation parameters ---- */ 
  /** @brief Interface for config::runtime */
  double getRuntime() const { return runtime; }
  
  /** @brief Interface for config::testparticles */
  int getTestparticles() const { return testparticles; }
  
  /** @brief Interface for config::freezeOutEnergyDensity */
  double getFreezeOutEnergyDensity() const { return freezeOutEnergyDensity; }
  
  /** @brief Interface for config::scatt_offlineWithAddedParticles */
  bool isScatt_offlineWithAddedParticles() const {return scatt_offlineWithAddedParticles;}
  
  /** @brief Interface for config::scatt_amongOfflineParticles */
  bool isScatt_amongOfflineParticles() const {return scatt_amongOfflineParticles;}
  
  /** @brief Interface for config::scatt_amongAddedParticles */
  bool isScatt_amongAddedParticles() const {return scatt_amongAddedParticles;}
  
  /** @brief Interface for config::N_light_flavors_input */
  int getNlightFlavorsAdded() const {return N_light_flavors_input;}
  
  /** @brief Interface for config::N_heavy_flavors_input */
  int getNheavyFlavorsAdded() const {return N_heavy_flavors_input;}
  
  /** @brief Interface for config::switchOff_23_32 */
  bool isSwitchOff_23_32() const {return switchOff_23_32;}
  /** ------------------------------- */
  
  /** ---- initial state options ---- */ 
  /** @brief Interface for config::initialStateType */
  INITIAL_STATE_TYPE getInitialStateType() const { return initialStateType; }
  
  /** @brief Interface for config::pythiaParticleFile */
  string getPythiaParticleFile() const { return pythiaParticleFile; }

  /** @brief Interface for config::cgcParticleFile */
  string getCgcParticleFile() const { return cgcParticleFile; }
  
  /** @brief Interface for config::mcatnloParticleFileCharm */
  string getMcatnloParticleFile() const {return mcatnloParticleFile;}
  
  /** @brief Interface for config::P0 */
  double getPtCutoff() const { return P0; }
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
 
  /** @brief Interface for config::standardOutputDirectoryName */
  string getStandardOutputDirectoryName() const { return standardOutputDirectoryName; }
  
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
  // open heavy flavor
  /** @brief Interface for config::KggQQbar */
  double getKggQQb() const {return KggQQbar;}
  
  /** @brief Interface for config::KgQgQ */
  double getKgQgQ() const {return KgQgQ;}
  
  /** @brief Interface for config::kappa_gQgQ */
  double getKappa_gQgQ() const {return kappa_gQgQ;}
  
  /** @brief Interface for config::couplingRunning */
  bool isCouplingRunning() const {return couplingRunning;}
  
  /** @brief Interface for config::isotropicCrossSecGQ */
  bool isIsotropicCrossSecGQ() const {return isotropicCrossSecGQ;}
  
  /** @brief Interface for config::constantCrossSecGQ */
  bool isConstantCrossSecGQ() const {return constantCrossSecGQ;}
  
  /** @brief Interface for config::constantCrossSecValueGQ */
  double getConstantCrossSecValueGQ() const {return constantCrossSecValueGQ; }


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
  
  
  // heavy quark output
  /** @brief Interface for config::hqCorrelationsOutput */
  bool isHqCorrelationsOutput() const {return hqCorrelationsOutput;}
  /** ------------------------------------ */

  /** -------- miscellaneous options ------- */ 
  bool repeatTimesteps() const { return switch_repeatTimesteps; }
  
  /** @brief Interface for config::interpolationBorder */
  double getMFPInterpolationBorder() const {return interpolationBorder;}

  /** @brief Interface for config::localCluster */
  bool isLocalCluster() const {return localCluster;}
  
  /** @brief Interface for config::jetMfpComputationSwitch */
  JET_MFP_COMPUTATION_TYPE getJetMfpComputationType() const {return jetMfpComputationSwitch;}
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
  /** ----------------------------------------------- */ 
  
  /** ----- auxiliary routines ----- */
  void readAndPrepareInitialSettings( offlineOutputInterface*const _offlineInterface );
  /** ------------------------------ */

 private:
  /** ----- auxiliary routines ----- */  
  /** @brief Write out all input parameters */
  void printUsedConfigurationParameters();
  
  /** @brief Process some settings needed for heavy quark runs */
  void processHeavyQuarkOptions();
  
  /** @brief Do some checks on user-provided options and parameters */
  void checkOptionsForSanity();
  
  /** @brief Create output directory if necessary */
  void checkAndCreateOutputDirectory( boost::filesystem::path& _dir );
  /** ------------------------------ */
  
  /** ------ general options ------- */  
  /** @brief The name of the current job - assigned to output files */
  string jobName;
  
  /** @brief Initial seed for the random number generator. seed = 0 forces the random generation of a seed */
  long seed;
  /** ------------------------------ */
  
  /** ---- collision parameters ---- */  
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
  /** ------------------------------ */
  
  /** ---- simulation parameters ---- */ 
  /** @brief Total simulated runtime in fm/c */
  double runtime;  

  /** @brief Number of testparticles per real particle */
  int testparticles;
  
  /** @brief Energy density for freeze out (in GeV/fm^3)
   * Particles in regions with energy densities below this threshold will be freely streaming
   */
  double freezeOutEnergyDensity;
  
  /** @brief number of active light flavors for added particles 
   * ( -1: no gluons nor light quarks, 0: only gluons, 1: including up, 2: including down, 3: including strange)
   */
  int N_light_flavors_input;
  /** @brief number of active heavy flavors for added particles 
   * ( 0: no charm and bottom, 1: only charm, 2: charm and bottom), see particleprototype.h, global parameter 
   */
  int N_heavy_flavors_input;
  
  /** @brief Whether offline particles are allowed to scatter with added particles */
  bool scatt_offlineWithAddedParticles;
  
  /** @brief Whether offline particles are allowed to scatter with other offline particles */
  bool scatt_amongOfflineParticles;
  
  /** @brief Whether added particles are allowed to scatter with other added particles */
  bool scatt_amongAddedParticles;
  
  /** @brief Whether 2->3 and 3->2 processed are switched off for added particles */
  bool switchOff_23_32;
  /** ------------------------------- */
  
  /** ---- initial state options ---- */ 
  /** @brief Which type of initial state to use */
  INITIAL_STATE_TYPE initialStateType;
                     
  /** @brief Relative or full path (including filename) of PYTHIA output file with particle data
   * Declare as "-" if particle momenta should be sampled via glauber method
   */
  string pythiaParticleFile;
  
  /** @brief Relative or full path (including filename) of color glass condensate output file with particle data
   * Declare as "-" if particle momenta should be sampled via glauber method
   */
  string cgcParticleFile; 
  
  /** @brief Relative or full path (including filename) of MC@NLO output file with particle data
   * Declare as "-" if particle momenta should be sampled via glauber method
   */
  string mcatnloParticleFile;
  
  /** @brief Lower PT-cutoff [GeV] used for minijet initial conditions */
  double P0;  
  /** ------------------------------- */
  
  /** -------- output options ------- */ 
  /** @brief Specify whether progress output should be written to a file */
  bool outputSwitch_progressLog;
  
  /** @brief Specify whether detailed particle output should be written */
  bool outputSwitch_detailedParticleOutput;
  
  /** @brief Specify whether movie output should be written (for added jet particles) */
  bool outputSwitch_movieOutputJets;

  /** @brief Specify whether movie output should be written (for the reconstructed background) */
  bool outputSwitch_movieOutputBackground;
  
  /** @brief Directory to which general output should be written */
  string standardOutputDirectoryName;
  
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
  // open heavy flavor
  /** @brief K factor for process g + g -> Q + Qbar */
  double KggQQbar;
  
  /** @brief K factor for process g + Q -> g + Q */
  double KgQgQ;
  
  /** @brief Kappa for Debye screening for process g + Q -> g + Q, usually 0.2 (Peshier,Gossiaux) */
  double kappa_gQgQ;
  
  /** @brief Whether a running coupling is employed for all process, for which running coupling is implemented (currently all processes involving heavy quarks) */
  bool couplingRunning;
  
  /** @brief Whether an isotropic momentum sampling is employed for process g + Q -> g + Q */
  bool isotropicCrossSecGQ;
  
  /** @brief Whether a constant cross section is employed for process g + Q -> g + Q */
  bool constantCrossSecGQ;
  
  /** @brief Value of constant cross section for process g + Q -> g + Q in mb */
  double constantCrossSecValueGQ;
  
  /** @brief Mass of charm quarks */
  double Mcharm_input;
  
  /** @brief Mass of bottom quarks */
  double Mbottom_input;
  
  
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
  
  
  // heavy quark output
  /** @brief Whether correlation analysis of heavy quark pairs is done */
  bool hqCorrelationsOutput;
  /** ------------------------------------ */
  
  /** -------- miscellaneous options ------- */ 
  /** @brief X where interpolation of MFP is done for E > X*T */
  double interpolationBorder; 
  
  /** @brief whether jobs should run on local itp cluster or CSC or the like
   * chose false for CSC and true for local queuing system, like housewifes, dwarfs etc. 
   * local system accesses personal filesystem, no need for transfering data from /local/ to working environment 
   */
  bool localCluster;
  
   /** @brief How to compute the mean free path high energy particles?
   *
   * 0 = computeMfpDefault = default, i.e. no special treatment
   * 1 = computeMfpIteration = iterative computation
   * 2 = computeMfpInterpolation = use tabulated mfp data and interpolation functions
   */
  JET_MFP_COMPUTATION_TYPE jetMfpComputationSwitch;
  
  /** @brief Whether timesteps are repeated in cases where the probability has been > 1 */
  bool switch_repeatTimesteps;
  /** ------------------------------------ */
  
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
  /** @brief Number of heavy ion collision events, set on top of the offline reconstruction. This input is neede if Pythia data files are read in to normalize the total yield of the particles. The number of one such heavy ion collision event includes < number of produced particles in pp > * Ntest * Nbin */
  int numberOfAddedEvents;
  /** @brief Minimum p_T [GeV] of the added particles */
  double minimumPT; 
  

  // the following parameters are read at runtime from the offline data recorded by the original run, 
  // see config::readAndPrepareInitialSettings
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
  /** ----------------------------------------------- */
};

inline void config::checkAndCreateOutputDirectory(boost::filesystem::path& _dir)
{
  if ( boost::filesystem::exists( _dir ) )
  {
    if ( boost::filesystem::is_directory( _dir ) )
    {
      return;
    }
    else
    {
      boost::filesystem::path renamePath( _dir.string() + ".backup" );
      std::cout << "File with name " << _dir.string() << " blocks the creation of an output folder for offline reconstruction." << std::endl;
      std::cout << "It is renamed to " << renamePath.string() << std::endl;
      boost::filesystem::rename( _dir, renamePath );
      boost::filesystem::create_directory( _dir );       
    }
  }
  else
  {
    std::cout << "Creating output folder: " << _dir.string() << std::endl;
    boost::filesystem::create_directory( _dir );    
  }
}

/** @brief exception class for handling unexpected critical behaviour within the configuration of the BAMPS run  */
class eConfig_error : public std::runtime_error
{
  public:
    explicit eConfig_error(const std::string& what) : std::runtime_error(what) {};
    
    virtual ~eConfig_error() throw() {};
};


#endif
