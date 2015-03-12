//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <fstream>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <boost/concept_check.hpp>

#include "configuration.h"
#include "particle.h"
#include "ringstructure.h"
#include "interpolation_nJpsi.h"

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
    int getIndex( const ParticleOffline& _particle ) const;
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


class analysisRapidityRange
{
  public:
    analysisRapidityRange() : yleft(0), yright(0) {};
    analysisRapidityRange( const double _left, const double _right ) : yleft(_left), yright(_right) {};
  
    void reset( const double _left, const double _right ) { yleft = _left; yright = _right; }
    double ydelta() const { return (yright - yleft); }
    
    double yleft;
    double yright;
};


class jetTrackerSingleEvent
{
public:
  jetTrackerSingleEvent() : 
    R_proj(), 
    P_proj_in(),
    P_proj_out(), 
    P1_in(), P2_in(), P1_out(), P2_out(),
    lambda(-1), 
    xSection(-1), 
    cell_ID(-1)
  { };
  ~jetTrackerSingleEvent() {};
  
  VectorTXYZ R_proj;
  
  VectorEPxPyPz P_proj_in;
  VectorEPxPyPz P_proj_out;
  
  VectorEPxPyPz P1_in;
  VectorEPxPyPz P2_in;
  VectorEPxPyPz P1_out;
  VectorEPxPyPz P2_out;
  
  int flavor_in;
  int flavor_out;
  
  jetTrackerCollType coll_type;
  double lambda;
  double xSection;
  int cell_ID;
  int jet_ID_in, jet_ID_out;
};



typedef boost::shared_array< std::vector<double> > tArrayOfDoubleVec;
typedef std::vector< tArrayOfDoubleVec > tVecOfArrayOfDoubleVec;

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
  
  void writePartclMovie( vector< ParticleOffline >& _particles, const int n_particles, fstream& _oscar, const int step, const int jumpSteps );
  
  void collectPtData( const int step );
  void collectPtDataInitial();
  void collectYData( const int step );
  void collectYDataInitial();
  void collectEtData( const int step );
  void collectEtDataInitial();
  
  void setSeed( uint32_t _s ) { seed = _s; }
  uint32_t getSeed( ) const { return seed; }


  int addJetEvent_in( const int ID_1, const int ID_2, const int ID_3, 
                      const jetTrackerCollType coll_type,
                      const double cross_section, const int cell_ID, 
                      const double lambda );
  void addJetEvent_initial( const int jetID );
  void addJetEvents_final();
  void addJetEvent_out( const int entity_ID, 
                        const int added_ID, const int ID_2, const int ID_3, 
                        const jetTrackerCollType coll_type );
  void removeJetEvent_in( const int entity_ID );
  void exchangeJetID( const int oldID, const int newID );
  void makeJetTrackerCopy();
  void restoreJetTracker();
  
  double getJetTracking_PT() const {return jetTracking_PT;}
  
  void registerProgressInformationForOutput( 
           const double _time, 
           const double _dt, 
           const int _nAddedParticles,
           const int _nMediumParticles, 
           const int _nColl, 
           const int _nColl22,
           const int _nColl23, 
           const int _nColl32 );
  void volumeMidrap(const int ) const;
  
private:
  config * const theConfig ;

public:

  double tstep[120];
  double tstep_movie[500];
  
  analysisRingStructure rings;
  ringStructure centralRingsCopyFromCascade;

  //--------------------//
  //for jet background analysis
  void scatteredMediumParticlesOutput( const int step );
  void mediumParticlesOutput( const int step );
  //--------------------//
private:

  uint32_t seed;
  time_t start;
  
  string filename_prefix;
  
  
  OUTPUT_SCHEME outputScheme;
  
  void handle_output_studies( OUTPUT_SCHEME _outputScheme );
  
  // switches for studies:
  bool studyHQ;
  bool studyJpsi;
  bool studyTempInTube;
  bool studyParticleOutput;
  bool studyTempAndVelocity;
  bool studyPtSpectra;
  bool studyEtSpectra;
  bool studyYDistribution;
  bool studyJets;
  bool studyCentralDensity;
  bool studyBackground;
  bool studyScatteredMediumParticles;
  
  
  
  void writePartclMovie( const int step ) const;
  void yDistribution( const FLAVOR_TYPE _flavTypeToComputeFor, vector< ParticleOffline >& _particles, const int n_particles, const int step );
  void transverseEnergyDistribution( const FLAVOR_TYPE _flavTypeToComputeFor, vector< ParticleOffline >& _particles, const int n_particles, const int step );
  void ptDistribution( const FLAVOR_TYPE _flavTypeToComputeFor, vector< ParticleOffline >& _particles, const int n_particles, const int step );
  void ptSoftDistribution( const FLAVOR_TYPE _flavTypeToComputeFor, vector< ParticleOffline >& _particles, const int n_particles, const int step );
  void particleOutput( const int step );
  void jetTrackerOutput();
  
  
  
  void printHeader( fstream & f, const anaType mode, const time_t end );
  
  void computeV2RAA( string, const double _outputTime );

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

  std::vector<analysisRapidityRange> rapidityRanges;
  
  double minPT, maxPT, binWidthPT;
  int numberBinsPT;
  std::vector<double> ptBinLabels;
  tVecOfArrayOfDoubleVec ptBins_gluons;
  tVecOfArrayOfDoubleVec ptBins_quarks;
  tVecOfArrayOfDoubleVec ptBins_ups;
  tVecOfArrayOfDoubleVec ptBins_downs;
  tVecOfArrayOfDoubleVec ptBins_stranges;
  tVecOfArrayOfDoubleVec ptBins_anti_ups;
  tVecOfArrayOfDoubleVec ptBins_anti_downs;
  tVecOfArrayOfDoubleVec ptBins_anti_stranges;
  tVecOfArrayOfDoubleVec ptBins_all;

  double maxPTSoft;
  double binWidthSoftPT;
  int numberBinsSoftPT;
  std::vector<double> *ptBinsSoftAll_gluons, *ptBinsSoftAll_quarks, *ptBinsSoftAll_all;
  std::vector<double> *ptBinsSoftAll_ups, *ptBinsSoftAll_downs, *ptBinsSoftAll_stranges, *ptBinsSoftAll_anti_ups, *ptBinsSoftAll_anti_downs, *ptBinsSoftAll_anti_stranges;
  std::vector<double> ptSoftBinLabels;
  
  double jetTracking_PT;  
  
  fstream oscarBackground;
  fstream oscarJets;
  fstream centralDensitiesOutputFile;
  fstream mfpJetsOutputFile;
  fstream progressLogFile;
  
  
  
  void calcTempCell( const int cell_id );
  void writeTempAndVel( const int step  );
  void addNeighborCells( const int cell_id, const int neighborCell_id );
  void writeTempInTube( const int step  );
  void calculateTempInTube( const double time, const double radius, const double dz, double & temp, double & tempWithQuarks, double & energyDensity  );
  void print_dndy(const string subfix );
  
  bool v2output;
  bool v2outputIntermediateSteps;
  bool dndyOutput;
  
  
  int *numberInCell; 
  int *temp_numberInCell; 
  double *vx_cell; 
  double *vy_cell;  
  double *vz_cell;
  double *vr_cell;
  double *em_cell;
  double *prm_cell;
  double *pzm_cell;
  double *pr2em_cell;
  double *pz2em_cell;
  double *przem_cell;
  double *densn_cell;
  double *gama_cell;
  double *temp_cell;
  double *tempWithQuarks_cell;
  double dv;
  int nCells;
  
  int *numberInCell_org;
  double *vx_cell_org; 
  double *vy_cell_org;  
  double *vz_cell_org;
  double *vr_cell_org;
  double *em_cell_org;
  double *prm_cell_org;
  double *pzm_cell_org;
  double *pr2em_cell_org;
  double *pz2em_cell_org;
  double *przem_cell_org;
  
  
  
  // heavy flavor stuff
  
  fstream printJpsiFugacity;
  fstream printTempInTube;
  
  void onePartclCorrelations();
  void twoPartclCorrelations();
  void jpsiEvolution(int step);
  void printJpsiEvolution();
  void analyseAngleDe();
  void writeJpsiFugacityOutput(const int step);
  void getJpsiFugacity(const double time, const double dr, const double dz, double& fugacity, int& n_charm, double& temp, double& deltaTemp);
  void jpsi_correlations();
  void ini_charm_correlations();
  
  interpolation_nJpsi theInterpolation_nJpsi;

  bool particleCorrelationsOutput;
  bool hadronization_hq;
  bool mesonDecay;

  int *numberJpsi_all_time;
  int *numberJpsi_ini_time;  

  int *numberJpsi_midPseudoRap_all_time;
  int *numberJpsi_midPseudoRap_ini_time;
  int *numberJpsi_forwardPseudoRap_all_time;
  int *numberJpsi_forwardPseudoRap_ini_time;
  int *numberJpsi_midNormRap_all_time;
  int *numberJpsi_midNormRap_ini_time;
  int *numberJpsi_forwardNormRap_all_time;
  int *numberJpsi_forwardNormRap_ini_time;
//   int *numberJpsi_midSpaceTimeRap_time;
  int *numberJpsiProd_time;
  int *numberCCbGG_time;
  bool *timestepAnalysed;

  int* numberJpsiDiss_time;
  int* numberJpsiDissTd_time;

  vector< vector < int > > showerParticlesInEvent;
  
};




class v2RAA
{
public:
  v2RAA( config * const c, string name_arg, string filename_prefix_arg, std::vector<analysisRapidityRange> rapidityRanges_arg, const double pt_min_arg = 5.0, const double pt_max_arg = 40.0, const int n_g_arg = 35, const double pt_min_background_arg = 0.0, const double pt_max_background_arg = 5.0, const int n_g_background_arg = 25 );
  
  void setPtBinProperties( const double pt_min_arg, const double pt_max_arg, const int n_g_arg, const double pt_min_background_arg = 0.0, const double pt_max_background_arg = 5.0, const int n_g_background_arg = 25 )
  {
    pt_min = pt_min_arg;
    pt_max = pt_max_arg;
    n_g = n_g_arg;
    pt_min_background = pt_min_background_arg;
    pt_max_background = pt_max_background_arg;
    n_g_background = n_g_background_arg;
  }
  
  void computeFor( const FLAVOR_TYPE _flavTypeToComputeFor, std::vector< ParticleOffline >& _particles, const int n_particles, string additionalNameTag, const double _outputTime, const v2Type _v2type );

private:
  
  config * const theConfig ;
  
  double pt_min,pt_max;  
  int n_c,n_b,n_g,eta_bins; 
  
  double pt_min_background, pt_max_background;
  int n_g_background;
  
  string name,filename_prefix;
  
  std::vector<analysisRapidityRange> rapidityRanges;
  

};



/** @brief exception class for handling unexpected critical behaviour within simulations of heavy ion collisions  */
class eAnalysis_error : public std::runtime_error
{
public:
  explicit eAnalysis_error( const std::string& what ) : std::runtime_error( what ) {};

  virtual ~eAnalysis_error() throw() {};
};




#endif

