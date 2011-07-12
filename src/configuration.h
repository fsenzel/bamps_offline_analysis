//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL: file:///home/fochler/svnrepos/cascade_full/main/trunk/src/configuration.h $
//$LastChangedDate: 2008-02-29 15:50:28 +0100 (Fri, 29 Feb 2008) $
//$LastChangedRevision: 34 $
//$LastChangedBy: fochler $
//---------------------------------------------
//---------------------------------------------


#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <string>
#include <stdexcept>
#include <vector>
#include <stdint.h>

#include "particle.h"


/** @brief namespace for global cascade paramters, objects etc. */
namespace ns_casc
{
  //-- global parameters -----------------------------//
  /** @brief (lambda_QCD)^2 in GeV^2, global parameter*/
  const double lambda2 = 0.04;
  
  /** @brief number of colors, global parameter */
  const int Ncolor = 3;
  
  /** @brief number of active flavors, see particle.h, global parameter */
  const int Nflavor = 3;
  
  /** @brief K-factor */
  const double K = 2.0;                

  /** @brief lower cutoff for some t-chanel processes */
  const double tcut = -0.1;

  
  const int first = 2; 
  
  /**
  * @brief auxiliary variable, denotes "a large number"
  *
  * yeah, I know that is's not really infinity, but what the heck
  */
  const double infinity = 300000;
  //--------------------------------------------------//
  
 
  /** @brief vector to hold all particles in the simulation, GLOBAL 
  *
  * Declared as extern here. There MUST be ONE definition in a cpp-file that includes this header in the global 
  * scope (it's in configuration.cpp).
  * All files that include this header can then access the vector.
  */
  extern std::vector<Particle> particles;
  extern std::vector<Particle> particles_init;
  extern std::vector<Particle> particles_atTimeNow;
  extern std::vector<Particle> particles_atTimeNowCopy;
  extern std::vector<Particle> addedParticles;
  extern std::vector<Particle> addedParticlesCopy;
}


extern double dt;




using std::string;

class config
{
 public:
  config(const int argc, const char * const argv[]);
  ~config() {};
  
  void setting();

  string getName() const {return name;}
  string getOutputName() const {return outputName;}
//   string getPythiaParticleFileCharm() const {return pythiaParticleFileCharm;}
//   int getPythiaParticleFileCharmTestparticles() const {return pythiaParticleFileCharmTestparticles;}
  int getNameLength() const;
  double getRuntime() const {return runtime;}
  double getDt() const {return delta_t;}
  bool DtSpecified() const {return dt_specified;}
  double getSqrtS() const {return sqrtS;}
  int getTestparticles() const {return testparticles;}
  double getImpactParameter() const {return impact;}
  int getOutputLevel() const {return outputLevel;}
  bool isLocalCluster() const {return localCluster;}
//   double getKIniCharm() const {return KIniCharm;}
  double getTimefirst() const {return timefirst;}
  double getA() const {return A;}
  double getAatomic() const {return Aatomic;}
  double getB() const {return B;}
  double getBatomic() const {return Batomic;}
  double getPtCutoff() const {return P0;}
  double getTimeshift() const {return timeshift;}
  string getPathdirCascadeData() const {return pathdirCascadeData;}
  uint32_t getSeed() const {return seed;}
  double getMinimumPT() const { return minimumPT; }
  int getNumberOfParticlesToAdd() const { return numberOfParticlesToAdd; }
  double getFreezeOutEnergyDensity() const { return freezeOutEnergyDensity; }
  double getFactor_dt() const { return factor_dt; }
  int getN_init() const { return N_init; }
  int getRingNumber() const { return ringNumber; }
  double getCentralRingRadius() const { return centralRingRadius; }
  double getDeltaR() const { return deltaR; }
  double get_dx() const { return dx; }
  double get_dy() const { return dy; }
  double getTransLen() const { return transLen; }

  bool movieOutputJets;
  bool movieOutputBackground;
  
  string getPythiaParticleFile() const { return pythiaParticleFile; }

 private:
  void readInputfile();

  string inputFile, pathdirCascadeData;

  double centralRingRadius, deltaR;
  int ringNumber;
  double timefirst,timeshift;
  double A;          //mass number of nucleus A  
  double Aatomic;    //atomic number, i.e. number of protons, of nucleus A   
  double B;          //mass number of nucleus B
  double Batomic;    //atomic number of nucleus B 
  
  double sqrtS;      //c.m. energy per NN pair, GeV
  double impact;     //impact parameter in fm
  double P0;         //lower PT-cutoff, GeV
  int testparticles; //number of testparticles per real particle
  
  double dx;
  double dy;
  double transLen;
  
  int getJobID();
  
  bool dt_specified;
  
  int numberOfParticlesToAdd;
  double minimumPT; 
  double freezeOutEnergyDensity;
  double factor_dt;

  double runtime;    //total simulated runtime in fm/c
  
  double delta_t;   // if dt is specified in input file, use this delta t. If not get it from cascade data.

  
  uint32_t originalSeed;    //initial seed used in the full BAMPS run
  uint32_t seed;            //seed used in the offline run
  
  int N_init;

  int outputLevel;   //1 = minimal (30M), 2 = with <name>_step?.f1 (110M), 3 = with reconstruction files (71M)
                     //4 = all (151M), all sizes are approx. for a run to 5.2 fm/c, uncompressed
                     //5 = only charm quark data
                     //6 = movie
                     //7 = v2 and R_AA
                     //8 = hydro: T_mu,nu
  
  
  string name;        //name of the BAMPS job that is processed (input file names)
  string outputName;  //name of this job - assigned to output files
  
  string pythiaParticleFile; //filepath to PYTHIA output file with particle data,  declare as "-" if particle momenta should be sampled via glauber method
//   string pythiaParticleFileCharm; // filepath to PYTHIA output file with particle data (can contain both charm quarks and gluons, the latter one being rejected) if using minijet or CGC for initial conditions, since both cannot sample charm quarks.  Declare as "-" if PYTHIA is used as initial model (contains already charm quarks). Usually one can use same PYTHIA particle files as for the case where PYTHIA is model for initial conditions
//   
//   int pythiaParticleFileCharmTestparticles; // number of testparticles used in particle data file, if it differs from global testparticle number take that into account
//   
//   string cgcParticleFile; //filepath to color glass condensate output file with particle data,  declare as "-" if particle momenta should be sampled via glauber method
//   
//   double KIniCharm; // K factor for initial charm for pythia charm file, K=2 means twice as much initial charm
//                     // This is different from the definition in the cascade. In cascade the initial charm from pythia is multiplied by this factor. Here we already increased the initial charm in Pythia with this factor and use it here only for the right scaling of the output. 
  
  bool localCluster; // chose false for CSC and true for local queuing system, like housewifes, dwarfs etc.
                                // local system accesses personal filesystem, no need for transfering data from /local/ to working environment
                                
};


/** @brief exception class for handling unexpected critical behaviour within the configuration  */
class eConfig_error : public std::runtime_error
{
public:
  explicit eConfig_error( const std::string& what ) : std::runtime_error( what ) {};

  virtual ~eConfig_error() throw() {};
};


#endif
