//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------
//---------------------------------------------


#include <math.h>
#include <string>


#include "cellcontainer.h"



cellContainer::cellContainer() :
  averagesPrepared( 0 ),
  nCollectedAll2223( 0 ),
  nCollected22( 0 ),
  nCollected23( 0 ),
  alpha_s_22( 0 ),
  alpha_s_23( 0 ),
  md2g_scaled_22( 0 ),
  md2q_scaled_22( 0 ),
  md2g_scaled_23( 0 ),
  md2q_scaled_23( 0 ),
  sigma_22( 0 ),
  sigma_23( 0 ),
  lambdaScaled( 0 ),
  number_gluon_cluster( 0.0 ),
  number_quark_cluster( 0.0 ),
  inverseE_gluons_cluster( 0.0 ),
  inverseE_quarks_cluster( 0.0 ),
  p_cluster( VectorEPxPyPz( 0.0, 0.0, 0.0, 0.0 ) )
{
  particleList.clear();
}



void cellContainer::clear()
{
  particleList.clear();
  rates.clear();
  corner.clear();

  nCollectedAll2223 = 0;
  nCollected22 = 0;
  nCollected23 = 0;
  alpha_s_22 = 0;
  alpha_s_23 = 0;
  sigma_22 = 0;
  sigma_23 = 0;
  md2g_scaled_22 = 0;
  md2q_scaled_22 = 0;
  md2g_scaled_23 = 0;
  md2q_scaled_23 = 0;
  lambdaScaled = 0;
  number_gluon_cluster =  0.0;
  number_quark_cluster =  0.0;
  inverseE_gluons_cluster =  0.0;
  inverseE_quarks_cluster =  0.0;
  p_cluster = VectorEPxPyPz( 0.0, 0.0, 0.0, 0.0 );
  averagesPrepared = false;
}


void cellContainer::resetStoredValues()
{
  nCollectedAll2223 = 0;
  nCollected22 = 0;
  nCollected23 = 0;
  alpha_s_22 = 0;
  alpha_s_23 = 0;
  sigma_22 = 0;
  sigma_23 = 0;
  md2g_scaled_22 = 0;
  md2q_scaled_22 = 0;
  md2g_scaled_23 = 0;
  md2q_scaled_23 = 0;
  lambdaScaled = 0;
  number_gluon_cluster =  0.0;
  number_quark_cluster =  0.0;
  inverseE_gluons_cluster =  0.0;
  inverseE_quarks_cluster =  0.0;
  p_cluster = VectorEPxPyPz( 0.0, 0.0, 0.0, 0.0 );

  rates.clear();
  averagesPrepared = false;
}



void cellContainer::setCoordinates( const int _index, const double _dx, const int _nx, const double _sizeX, const double _dy, const int _ny, const double _sizeY )
{
  const int nxny = _nx * _ny;
  const int indexEta = _index / nxny;
  const int indexY = _index - ( _index / nxny ) * nxny;
  const int indexX = indexY - ( indexY / _nx ) * _nx;

  double leftY = -( _sizeY / 2.0 ) + _dy * ( indexY / _nx );
  double leftX = -( _sizeX / 2.0 ) + _dx * indexX;

  corner.setCorners( leftX, leftX + _dx, leftY, leftY + _dy, indexEta );
  index = _index;
}



void cellContainer::prepareAverages()
{
  if ( !averagesPrepared )
  {
    averagesPrepared = true;
    if ( nCollectedAll2223 > 0 )
    {
      sigma_22 /= static_cast<double>( nCollectedAll2223 );
      sigma_23 /= static_cast<double>( nCollectedAll2223 );
    }
    else
    {
      sigma_22 = 0;
      sigma_23 = 0;
    }

    if ( nCollected22 > 0 )
    {
      md2g_scaled_22 /= static_cast<double>( nCollected22 );
      md2q_scaled_22 /= static_cast<double>( nCollected22 );
      alpha_s_22 /= static_cast<double>( nCollected22 );
    }
    else
    {
      md2g_scaled_22 = 0;
      md2q_scaled_22 = 0;
      alpha_s_22 = 0;
    }

    if ( nCollected23 > 0 )
    {
      md2g_scaled_23 /= static_cast<double>( nCollected23 );
      md2q_scaled_23 /= static_cast<double>( nCollected23 );
      alpha_s_23 /= static_cast<double>( nCollected23 );
      lambdaScaled /= static_cast<double>( nCollected23 );
    }
    else
    {
      md2g_scaled_23 = 0;
      md2q_scaled_23 = 0;
      alpha_s_23 = 0;
      lambdaScaled = 0;
    }
  }
  else
  {
    std::string errMsg = "prepareAverages called for cell that has already been averaged";
    throw eCell_error( errMsg );
  }
}

// The following routine should not be used in full/offlineAnalysis !!!

void cellContainer::writeAveragesToParticle( Particle& _particle ) const
{
  if ( !averagesPrepared )
  {
    std::string errMsg = "writeAveragesToParticle(..) called without prior call to prepareAverages()";
    throw eCell_error( errMsg );
  }
  
  _particle.cs22 = sigma_22;               //1/GeV^2
  _particle.cs23 = sigma_23;               //1/GeV^2
  _particle.md2g_scaled_22 = md2g_scaled_22;
  _particle.md2q_scaled_22 = md2q_scaled_22;
  _particle.md2g_scaled_23 = md2g_scaled_23;
  _particle.md2q_scaled_23 = md2q_scaled_23;
  _particle.as22 = alpha_s_22;
  _particle.as23 = alpha_s_23;
  _particle.lambda_scaled = lambdaScaled;
  
}

void cellContainer::addParticleToCluster( Particle particleToAdd )
{
  p_cluster += particleToAdd.Mom;

  if ( particleToAdd.FLAVOR == gluon )
  {
    inverseE_gluons_cluster += 1.0 / particleToAdd.Mom.E();
    ++number_gluon_cluster;
  }
  else
  {
    inverseE_quarks_cluster += 1.0 / particleToAdd.Mom.E();
    ++number_quark_cluster;
  }
}

void cellContainer::getClusterInformation( double &T, VectorTXYZ &beta, const double minNumber )
{
  const double number_cluster = number_gluon_cluster + number_quark_cluster;

  if( number_cluster >= minNumber )
  {
    beta = p_cluster.NormalizeToE();

    const double gamma = 1.0 / sqrt( 1 - beta.vec2() );

    // T* as calculated in JHEP 0301 (2003) 030, eqs. (1.6), (1.7)
    const double Cf = ( 4.0 / 3.0 );
    const double Ca = 3.0;
    const double I = 0.5 * ( Cf * number_quark_cluster + Ca * number_gluon_cluster ) / gamma; // unitless, gamma comes from boost of volume in partice density
    const double J = Cf * inverseE_quarks_cluster + Ca * inverseE_gluons_cluster; // 1 / GeV
   
    T = I / J;
  }
  else
  {
    beta = VectorTXYZ( 1.0, 0.0, 0.0, 0.0 );
    T = 0.0;
  }
}


cornerCoordinates::cornerCoordinates() :
    x_min( 0 ),
    x_max( 0 ),
    y_min( 0 ),
    y_max( 0 ),
    etaIndex( -1 )
{

}


cornerCoordinates::cornerCoordinates( const double _x_min, const double _x_max, const double _y_min, const double _y_max, const int _etaIndex ) :
    x_min( _x_min ),
    x_max( _x_max ),
    y_min( _y_min ),
    y_max( _y_max ),
    etaIndex( _etaIndex )
{
}



void cornerCoordinates::setCorners( const double _x_min, const double _x_max, const double _y_min, const double _y_max,  const int _etaIndex )
{
  x_min = _x_min;
  x_max = _x_max;
  y_min = _y_min;
  y_max = _y_max;
  etaIndex = _etaIndex;
}


double cornerCoordinates::getVolume( const coordinateEtaBins& _etaBins, const double _time ) const
{
  double deltaZ = _time * ( tanh( _etaBins[etaIndex].right ) - tanh( _etaBins[etaIndex].left ) );
  return (( x_max - x_min ) * ( y_max - y_min ) * deltaZ );
}

// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;
