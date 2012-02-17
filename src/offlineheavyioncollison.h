//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------
//---------------------------------------------


#ifndef OFFLINEHEAVYIONCOLLISON_H
#define OFFLINEHEAVYIONCOLLISON_H


#include <fstream>
#include <stdexcept>
#include <vector>
#include <algorithm>

#include "configuration.h"
#include "analysis.h"
#include "interpolation23.h"
#include "interpolation22.h" 
#include "woodsaxon.h"
#include "ringstructure.h"
#include "cellcontainer.h"
#include "scattering23.h"
#include "scattering22.h"
#include "scattering32.h"
#include "offlineoutput.h"


class offlineHeavyIonCollision
{
public:
  offlineHeavyIonCollision( config* const _config, offlineOutputInterface* const _offlineInterface );
  ~offlineHeavyIonCollision();

  void init();
  void mainFramework( analysis& aa );

  double evolveMedium( const double evolveToTime, bool& _endOfDataFiles );
  void scattering( const double nexttime, bool& again, analysis& aa );
  void cell_ID( double time );

  void scatterEdgeParticles( std::list< int >& _offlineParticleList, std::list< int >& _addedParticleList, const double nexttime );

private:
  /** @brief Pointer to the interface for the output of data needed for later offline reconstruction */
  offlineOutputInterface* offlineInterface;
  
  config * const theConfig;
  interpolation23 theI23;
  
  /** @brief  interpolation22 object that provides access to tabulated values for the cross section of all 2->2 processes with running coupling */
  interpolation22 theI22;

  //--parameters taken from config--------------------------------------------
  double stop;       //total simulated runtime in fm/c
  double A;          //mass number of nucleus A
  double Aatomic;    //atomic number, i.e. number of protons, of nucleus A
  double B;          //mass number of nucleus B
  double Batomic;    //atomic number of nucleus B
  double sqrtS;      //c.m. energy per NN pair, GeV
  double P0;         //lower PT-cutoff, GeV
  double Bimp;       //impact parameter in fm
  int testpartcl;    //number of testparticles per real particle
  //--------------------------------------------------------------------------

  ringStructure rings;
  
  std::vector< std::vector<double> > rateGluons;
  std::vector< std::vector<double> > rateQuarks;
  std::vector< std::vector<double> > rateAntiQuarks;
  std::vector< std::vector<double> > rateGluons_prev;
  std::vector< std::vector<double> > rateQuarks_prev;
  std::vector< std::vector<double> > rateAntiQuarks_prev;

  int currentNumber;

  double stoptime;
  double stoptime_last;
  
  int numberEvolvingParticles;

  WoodSaxon WoodSaxonParameter;


  void scatt2223_offlineWithAddedParticles( cellContainer& _cell, std::vector< int >& _allParticlesList, std::vector< int >& _gluonList,
                                     cellContainer& _cellAdded, std::vector< int >& _allParticlesListAdded, std::vector< int >& _gluonListAdded,
                                     const double scaleFactor, bool& again, analysis& aa, const double nexttime, analysisRingStructure& _analysisRings );
  
  void scatt22_amongAddedParticles( cellContainer& _cellAdded, std::vector< int >& _allParticlesListAdded, const double scaleFactor, bool& again, const double nexttime );

  void scatt32_offlineWithAddedParticles( cellContainer& _cell, std::vector< int >& _allParticlesList, std::vector< int >& _gluonList,
                                   cellContainer& _cellAdded, std::vector< int >& _allParticlesListAdded, std::vector< int >& _gluonListAdded,
                                   int& n32, bool& again, analysis& aa, const double nexttime );

  int scatt32_offlineWithAddedParticles_utility( scattering32& scatt32_obj, std::list< int >& _cellMembersAdded, std::vector< int >& _allParticlesListAdded, std::vector< int >& _gluonListAdded, const int iscat, const int jscat, const int kscat, int& n32, const double picked_ratio, const double nexttime  );
  int scatt23_offlineWithAddedParticles_utility( scattering23& scatt23_obj, cellContainer& _cell, int iscat, const int jscat, bool& again, const double nexttime );
  void scatt22_offlineWithAddedParticles_utility( scattering22& scatt22_obj, const int iscat, const int jscat, int& typ, const double nexttime );
  void scatt22_amongAddedParticles_utility( scattering22& scatt22_obj, std::list< int >& _cellMembersAdded, std::vector< int >& _allParticlesListAdded, const int iscat, const int jscat, int& typ, const double nexttime );
  
  /** @brief Goes through all added particles and decays J/psi if the temperature of the surrounding medium is larger than the dissocation temperature */
  void jpsi_dissociation_td( const double time );

  double iterateMFP( std::vector< int >& _allParticlesList, std::vector< int >& _gluonList, const int jetID, const double dt, const double dv );
  
  void removeDeadParticles( analysis& _aa );
  
  int binomial( const int N, const int k ) const;
};


/** @brief exception class for handling unexpected critical behaviour within simulations of heavy ion collisions  */
class eHIC_error : public std::runtime_error
{
public:
  explicit eHIC_error( const std::string& what ) : std::runtime_error( what ) {};

  virtual ~eHIC_error() throw() {};
};



/** @brief Utility function that finds an element in a vector and subsequently removes it
*
* Usage: removeElementFromVector<int> ( vector, element)
*
* @param[in,out] _vec Vector from which _elementToRemove needs to be erased
* @param[in] _elementToRemove The element that should be removed
*/
template <class T>
void removeElementFromVector( std::vector<T>& _vec, const T _elementToRemove )
{
  typename std::vector<T>::iterator findIter;
  findIter = std::find( _vec.begin(), _vec.end(), _elementToRemove );
  
  if ( findIter != _vec.end() )
  {
    _vec.erase( findIter );
  }
  else
  {
    std::string errMsg = "Removal of element from vector failed. Unrecoverable error.";
    throw eHIC_error( errMsg );
  }
}
#endif // OFFLINEHEAVYIONCOLLISON_H
