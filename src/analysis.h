#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <stdint.h>
#include <time.h>
#include <fstream>
#include <vector>
#include <math.h>

#include "configuration.h"
#include "particle.h"
#include "ringstructure.h"

using std::vector;
using std::fstream;


class config;

enum jetTrackerCollType {initial_jet, final_jet, c2to2, c2to3, c3to2, production};
enum anaType {ptSpectrum, ptSpectrumSoft, rapidityDistribution, quarkNumbers, all, initial, final, jets};
enum v2Type {v2jets, v2background};

const string sep = "\t";


class analysisRingContainer
{
  public:
    analysisRingContainer() : minRadius( 0 ), maxRadius( 0 ), deltaR( 0 ), lambdaGluon(0), lambdaQuark(0), collectedGluon(0), collectedQuark(0) {};
    analysisRingContainer( const double _minR, const double _maxR ) : minRadius( _minR ), maxRadius( _maxR ), deltaR( _maxR - _minR ), lambdaGluon(0), lambdaQuark(0), collectedGluon(0), collectedQuark(0) {};
    ~analysisRingContainer() {};
    
    analysisRingContainer& operator+=( const analysisRingContainer& rhs );   
    
    void clear() { lambdaGluon = 0; lambdaQuark = 0; collectedGluon = 0; collectedQuark = 0; }
    void relocate( const double _minR, const double _maxR )
    {
      minRadius = _minR;
      maxRadius = _maxR;
      deltaR = _maxR - _minR;
      clear();
    }
    
    double getVolume( const double _dz ) const
    {
      return ( M_PI * ( pow( maxRadius, 2 ) - pow( minRadius, 2 ) ) * _dz );   //fm^3
    }
    
    double minRadius;
    double maxRadius;
    double deltaR;
    
    double lambdaGluon;
    double lambdaQuark;
    int collectedGluon;
    int collectedQuark;
};



class analysisRingStructure
{
  public:
    analysisRingStructure() : numberOfRings( 0 ), centralRingRadius( 0 ), totalRadius( 0 ), deltaR( 0 ) { rings.resize(0); }
    analysisRingStructure( const int _nRings, const double _centralRadius, const double _deltaR );
    ~analysisRingStructure() {};
    
    void resize( const int _nRings, const double _centralRadius, const double _deltaR );
    
    analysisRingContainer& getRing( const double _xt ) { return rings[ getIndex( _xt ) ]; }
    analysisRingContainer& operator[]( const int _index );
    
    analysisRingStructure& operator+=( analysisRingStructure& rhs );    
    
    int getIndex( const double _xt ) const;
    int getIndexPure( const double _xt ) const;
    int getIndex( const Particle& _particle ) const;
    int size() const { return numberOfRings; }
    double getCentralRadius() const { return centralRingRadius; }
    double getDeltaR() const { return deltaR; }
    
    void clear() 
    { 
      for( int i = 0; i < numberOfRings; i++ ) 
      { 
        rings[i].clear(); 
      }
    }  
    
  private:
    std::vector<analysisRingContainer> rings;
    
    int numberOfRings;
    double centralRingRadius;
    double totalRadius;
    double deltaR;
};




class jetTrackerSingleEvent
{
  public:
    jetTrackerSingleEvent();
    ~jetTrackerSingleEvent();
    
    double R_proj[4];
    
    double P_proj_in[4];
    double P_proj_out[4];
    
    double P1_in[4], P2_in[4];
    double P1_out[4], P2_out[4];
    
    int flavor_in;
    int flavor_out;
    
    jetTrackerCollType coll_type;
    double lambda;
    double xSection;
    int cell_ID;
    int jet_ID_in, jet_ID_out;
};


class analysis
{
public:
  analysis(config * const c);  
  ~analysis();
  
  void initialOutput();
  void intermediateOutput(const int nn);
  void finalOutput( const double _stoptime );
  void movieOutput( const int step, const int jumpSteps );
  void movieOutputMedium( const int step, const int jumpSteps );
  
  void mfpJetsOutput( const int step, const int jumpSteps );  
  
  void printCentralDensities( const double _time );
  
  void writePartclMovie( vector< Particle >& _particles, const int n_particles, fstream& _oscar, const int step, const int jumpSteps );
  
  void collectPtData( const int step );
  void collectPtDataInitial();
  void collectYData( const int step );
  void collectYDataInitial();
  void collectEtData( const int step );
  void collectEtDataInitial();
  
  void setSeed( uint32_t _s ) { seed = _s; }


  int addJetEvent_in( const int ID_1, const int ID_2, const int ID_3, const jetTrackerCollType coll_type,
  const double cross_section, const int cell_ID, const double lambda );
  void addJetEvent_initial( const int jetID );
  void addJetEvents_final();
  void addJetEvent_out( const int entity_ID, const int added_ID, const int ID_2, const int ID_3, const jetTrackerCollType coll_type );
  void removeJetEvent_in( const int entity_ID );
  void exchangeJetID( const int oldID, const int newID );
  void makeJetTrackerCopy();
  void restoreJetTracker();
  
  double getJetTracking_PT() const {return jetTracking_PT;}
  
  void volumeMidrap(const int ) const;
  
  double tstep[65];
  double tstep_movie[500];
  
  analysisRingStructure rings;
  ringStructure centralRingsCopyFromCascade;

private:

  config * const theConfig ;
  
  uint32_t seed;
  time_t start;
  
  string filename_prefix;
  
  void writePartclMovie( const int step ) const;
  void yDistribution( const FLAVOR_TYPE _flavTypeToComputeFor, vector< Particle >& _particles, const int n_particles, const int step );
  void transverseEnergyDistribution( const FLAVOR_TYPE _flavTypeToComputeFor, vector< Particle >& _particles, const int n_particles, const int step );
  void ptDistribution( const FLAVOR_TYPE _flavTypeToComputeFor, vector< Particle >& _particles, const int n_particles, const int step );
  void ptSoftDistribution( const FLAVOR_TYPE _flavTypeToComputeFor, vector< Particle >& _particles, const int n_particles, const int step );
  void particleOutput( const int step );
  void jetTrackerOutput();
  
  
  
  void printHeader( fstream & f, const anaType mode, const time_t end );
  
  void computeV2RAA( string, const double _outputTime );
  void onePartclCorrelations();
  void twoPartclCorrelations();
  void printPtSpectra( const FLAVOR_TYPE _flavTypeToComputeFor);
  void printSoftPtSpectra( const FLAVOR_TYPE _flavTypeToComputeFor );
  void printYDistribution();
  
  
  int nTimeSteps;
  int nTimeSteps_movie;
  
  std::vector< std::vector<jetTrackerSingleEvent> > jetTracker;
  std::vector< std::vector<jetTrackerSingleEvent> > jetTracker_copy;
  
  
  std::vector<double> *yBins_gluon, *yBins_up, *yBins_down, *yBins_strange, *yBins_anti_up, *yBins_anti_down, *yBins_anti_strange;
  std::vector<double> *transverseEnergyGluons, *transverseEnergyQuarks, *transverseEnergyAntiQuarks;
  std::vector<double> yBinLabels;
  double minY, maxY, binWidthY;
  int numberBinsY;  
  
  std::vector<double> *ptBinsDY1_gluons, *ptBinsDY2_gluons, *ptBinsDY3_gluons, *ptBinsDY4_gluons, *ptBinsDY5_gluons, *ptBinsAll_gluons;
  std::vector<double> *ptBinsDY1_quarks, *ptBinsDY2_quarks, *ptBinsDY3_quarks, *ptBinsDY4_quarks, *ptBinsDY5_quarks, *ptBinsAll_quarks;
  std::vector<double> *ptBinsDY1_ups, *ptBinsDY2_ups, *ptBinsDY3_ups, *ptBinsDY4_ups, *ptBinsDY5_ups, *ptBinsAll_ups;
  std::vector<double> *ptBinsDY1_downs, *ptBinsDY2_downs, *ptBinsDY3_downs, *ptBinsDY4_downs, *ptBinsDY5_downs, *ptBinsAll_downs;
  std::vector<double> *ptBinsDY1_stranges, *ptBinsDY2_stranges, *ptBinsDY3_stranges, *ptBinsDY4_stranges, *ptBinsDY5_stranges, *ptBinsAll_stranges;
  std::vector<double> *ptBinsDY1_anti_ups, *ptBinsDY2_anti_ups, *ptBinsDY3_anti_ups, *ptBinsDY4_anti_ups, *ptBinsDY5_anti_ups, *ptBinsAll_anti_ups;
  std::vector<double> *ptBinsDY1_anti_downs, *ptBinsDY2_anti_downs, *ptBinsDY3_anti_downs, *ptBinsDY4_anti_downs, *ptBinsDY5_anti_downs, *ptBinsAll_anti_downs;
  std::vector<double> *ptBinsDY1_anti_stranges, *ptBinsDY2_anti_stranges, *ptBinsDY3_anti_stranges, *ptBinsDY4_anti_stranges, *ptBinsDY5_anti_stranges, *ptBinsAll_anti_stranges;
  std::vector<double> *ptBinsDY1_all, *ptBinsDY2_all, *ptBinsDY3_all, *ptBinsDY4_all, *ptBinsDY5_all, *ptBinsAll_all;
  std::vector<double> ptBinLabels;
  double minPT, maxPT, binWidthPT;
  double maxPTSoft;
  int numberBinsPT;
  std::vector<double> *ptBinsSoftAll_gluons, *ptBinsSoftAll_quarks, *ptBinsSoftAll_all;
  std::vector<double> *ptBinsSoftAll_ups, *ptBinsSoftAll_downs, *ptBinsSoftAll_stranges, *ptBinsSoftAll_anti_ups, *ptBinsSoftAll_anti_downs, *ptBinsSoftAll_anti_stranges;
  std::vector<double> ptSoftBinLabels;
  double binWidthSoftPT;
  int numberBinsSoftPT;
  
  double jetTracking_PT;  
  
  fstream oscarBackground;
  fstream oscarJets;
  fstream centralDensitiesOutputFile;
  fstream mfpJetsOutputFile;
};




class v2RAA
{
public:
  v2RAA( config * const c, string name_arg, string filename_prefix_arg );
  
  void computeFor( const FLAVOR_TYPE _flavTypeToComputeFor, std::vector< Particle >& _particles, const int n_particles, string additionalNameTag, const double _outputTime, const v2Type _v2type );

private:
  
  config * const theConfig ;
  
  double pt_min,pt_max;  
  int n_c,n_b,n_g,eta_bins; 
  
  double pt_min_background, pt_max_background;
  int n_g_background;
  
  string name,filename_prefix;
  

};



/** @brief exception class for handling unexpected critical behaviour within simulations of heavy ion collisions  */
class eAnalysis_error : public std::runtime_error
{
public:
  explicit eAnalysis_error( const std::string& what ) : std::runtime_error( what ) {};

  virtual ~eAnalysis_error() throw() {};
};




#endif

