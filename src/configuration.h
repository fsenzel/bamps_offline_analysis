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

#include "globalsettings.h"
#include "particle.h"
#include "offlineoutput.h"

#define HELP_MESSAGE_REQUESTED 1


/** @brief Enumeration type for possible initial state models */
enum INITIAL_STATE_TYPE { miniJetsInitialState, pythiaInitialState, cgcInitialState };

/** @brief Enumeration type for different variants of computing the mean free path of high-pt particles */
enum JET_MFP_COMPUTATION_TYPE { computeMfpDefault, computeMfpIteration, computeMfpInterpolation };


/** @brief This extends the namespace ns_casc (see globalsettings.h) to hold a global particle vector */
namespace ns_casc
{  
  /** @brief Vector to hold all particles in the simulation, GLOBAL 
  *
  * Declared as extern here. There MUST be ONE definition in a cpp-file that includes this header in the global 
  * scope (it's in configuration.cpp).
  * All files that include this header can then access the vector.
  */
  extern std::vector<ParticleOffline> particles;
  extern std::vector<ParticleOffline> particles_init;
  extern std::vector<ParticleOffline> particles_atTimeNow;
  extern std::vector<ParticleOffline> particles_atTimeNowCopy;
  extern std::vector<ParticleOffline> addedParticles;
  extern std::vector<ParticleOffline> addedParticlesCopy;
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
  /** ------------------------------- */
  
  /** ---- initial state options ---- */ 
  /** @brief Interface for config::initialStateType */
  INITIAL_STATE_TYPE getInitialStateType() const { return initialStateType; }
  
  /** @brief Interface for config::pythiaParticleFile */
  string getPythiaParticleFile() const { return pythiaParticleFile; }

  /** @brief Interface for config::cgcParticleFile */
  string getCgcParticleFile() const { return cgcParticleFile; }
  
  /** @brief Interface for config::P0 */
  double getPtCutoff() const { return P0; }
  /** ------------------------------- */
  
  /** -------- output options ------- */   
  /** @brief Interface for config::outputSwitch_progressLog */
  bool doOutput_progressLog() const { return outputSwitch_progressLog; }  

  /** @brief Interface for config::outputSwitch_detailedParticleOutput */
  bool doOutput_detailedParticleOutput() const { return outputSwitch_detailedParticleOutput; }
  
  /** @brief Interface for config::outputSwitch_movieOutput */
  bool doOutput_movieOutput() const {return outputSwitch_movieOutput;}
 
  /** @brief Interface for config::standardOutputDirectoryName */
  string getStandardOutputDirectoryName() const { return standardOutputDirectoryName; }
  /** ------------------------------- */
  
 /** @brief Interface for config::name */
  string getName() const {return name;}
  /** @brief Interface for config::interpolationSwitch */
  JET_MFP_COMPUTATION_TYPE getJetMfpComputationType() const {return jetMfpComputationSwitch;}
  /** @brief Interface for config::interpolationSwitch */
  double getMFPInterpolationBorder() const {return interpolationBorder;}


  /** @brief Interface for config::pathdirOfflineData */
  string getPathdirOfflineData() const { return pathdirOfflineData; }
  /** @brief Interface for config::pathdirOfflineData, returns pointer to char */
  const char* getPathdirOfflineDataChar() const { return pathdirOfflineData.c_str(); }
  
  /** Stuff special to offline analysis */
  void readAndPrepareInitialSettings(offlineOutputInterface*const _offlineInterface);
  
  double getMinimumPT() const { return minimumPT; }
  int getNumberOfParticlesToAdd() const { return numberOfParticlesToAdd; }
  bool movieOutputJets;
  bool movieOutputBackground;
  bool isLocalCluster() const {return localCluster;}
  double getTimeshift() const {return timeshift;}
  int getRingNumber() const { return ringNumber; }
  double getCentralRingRadius() const { return centralRingRadius; }
  double getDeltaR() const { return deltaR; }
  double getDt() const {return delta_t;}
  bool DtSpecified() const {return dt_specified;}
  double getFactor_dt() const { return factor_dt; }
  double getTimefirst() const {return timefirst;}
  
  int getN_init() const { return N_init; }
  double get_dx() const { return dx; }
  double get_dy() const { return dy; }
  double getTransLen() const { return transLen; }
 

 private:
  /** ----- auxiliary routines ----- */  
  /** @brief Write out all input parameters */
  void printUsedConfigurationParameters();
  
  /** @brief Do some checks on user-provided options and parameters */
  void checkOptionsForSanity();
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
  
  /** @brief Lower PT-cutoff [GeV] used for minijet initial conditions */
  double P0;  
  /** ------------------------------- */
  
  /** -------- output options ------- */ 
  /** @brief Specify whether progress output should be written to a file */
  bool outputSwitch_progressLog;
  
  /** @brief Specify whether detailed particle output should be written */
  bool outputSwitch_detailedParticleOutput;
  
  /** @brief Specify whether movie output should be written */
  bool outputSwitch_movieOutput;
  
  /** @brief Directory to which general output should be written */
  string standardOutputDirectoryName;
  /** ------------------------------- */
  
  
  
  
  /** @brief Path to which the output needed for offline analysis is written */
  string pathdirOfflineData;


  /** @brief How to compute the mean free path high energy particles?
   *
   * 0 = computeMfpDefault = default, i.e. no special treatment
   * 1 = computeMfpIteration = iterative computation
   * 2 = computeMfpInterpolation = use tabulated mfp data and interpolation functions
   */
  JET_MFP_COMPUTATION_TYPE jetMfpComputationSwitch;
  /** @brief X where interpolation of MFP is done for E > X*T */
  double interpolationBorder; 
  
  
  /** @brief The name of the BAMPS job that is processed (input file names) */
  string name;

  /** @brief Directory to which general output should be written */
  string standardOutputDirectoryName;
  
  
  /** Stuff special to the offline analysis */
  /** @brief Initial seed used in the full BAMPS run */
  uint32_t originalSeed;
  
  double delta_t;   // if dt is specified in input file, use this delta t. If not get it from cascade data.
  bool dt_specified;
  
  double centralRingRadius, deltaR;
  int ringNumber;
  double timeshift;
  int numberOfParticlesToAdd;
  double minimumPT; 
  double factor_dt;
  bool localCluster; // chose false for CSC and true for local queuing system, like housewifes, dwarfs etc., local system accesses personal filesystem, no need for transfering data from /local/ to working environment 
  double dx;
  double dy;
  double transLen;
  int N_init;
  double timefirst;
};



/** @brief exception class for handling unexpected critical behaviour within the configuration of the BAMPS run  */
class eConfig_error : public std::runtime_error
{
  public:
    explicit eConfig_error(const std::string& what) : std::runtime_error(what) {};
    
    virtual ~eConfig_error() throw() {};
};


#endif
