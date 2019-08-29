//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


#ifndef CELLCONTAINER_H
#define CELLCONTAINER_H

#include <list>
#include <stdexcept>

#include "ratesmanager.h"
#include "coordinateBins.h"
#include "particle.h"
#include "globalsettings.h"


class cornerCoordinates
{
public:
  cornerCoordinates();
  cornerCoordinates( const double _x_min, const double _x_max, const double _y_min, const double _y_max, const int _etaIndex );
  ~cornerCoordinates() {};

  void setCorners( const double _x_min, const double _x_max, const double _y_min, const double _y_max, const int _etaIndex );
  void clear() { setCorners( 0, 0, 0, 0, 0 ); }
  double getVolume( const coordinateEtaBins& _etaBins, const double _time ) const;


  double x_min;
  double x_max;
  double y_min;
  double y_max;
//   double z_min;
//   double z_max;
//   double eta_min;
//   double eta_max;
  int etaIndex;
};


class cellContainer
{
public:
  cellContainer();
  ~cellContainer() {};

  std::list<int> particleList;

  ratesManager rates;

  void clear();
  void resetStoredValues();
  bool empty() const { return particleList.empty(); }
  int size() const { return particleList.size(); }
  void setCoordinates( const int _index, const double _dx, const int _nx, const double _sizeX, const double _dy, const int _ny, const double _sizeY );
  void prepareAverages();
  void writeAveragesToParticle( Particle& _particle ) const;
  
  void getClusterInformation( double &T, VectorTXYZ &beta, const double minNumber, const double Ntest, const bool calcTviaEdens );
  void addParticleToCluster( Particle particleToAdd );
  void addVolumeToCluster( const double dV ) { volume_cluster += dV; };

  int index;
  cornerCoordinates corner;

  bool averagesPrepared;

  int nCollectedAll2223;
  int nCollected22;
  int nCollected23;

  double alpha_s_22;
  double alpha_s_23;

  double md2g_scaled_22;
  double md2q_scaled_22;
  double md2g_scaled_23;
  double md2q_scaled_23;

  double sigma_22;
  double sigma_23;

  double lambdaScaled;

//  AMY
//  Attention: These quantities are not for single cell but cell + neighbour cells = cluster.
  double number_gluon_cluster;
  double number_quark_cluster;
  double inverseE_gluons_cluster;
  double inverseE_quarks_cluster;
  VectorEPxPyPz p_cluster;
  double E_cluster;
  double px_cluster;
  double py_cluster;
  double pz_cluster;
  double pxx_cluster;
  double pyy_cluster;
  double pzz_cluster;
  double pxy_cluster;
  double pxz_cluster;
  double pyz_cluster;
  double pz2overE2_cluster;
  double pt2overpz2_cluster;
  double volume_cluster;

private:

};


/** @brief exception class for handling unexpected critical behaviour within cell objects  */
class eCell_error : public std::runtime_error
{
public:
  explicit eCell_error( const std::string& what ) : std::runtime_error( what ) {};

  virtual ~eCell_error() throw() {};
};

#endif // CELLCONTAINER_H
// kate: indent-mode cstyle; space-indent on; indent-width 0;
