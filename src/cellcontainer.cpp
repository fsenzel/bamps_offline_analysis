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
  p_cluster( VectorEPxPyPz( 0.0, 0.0, 0.0, 0.0 ) ),
  E_cluster( 0.0 ),
  px_cluster( 0.0 ),
  py_cluster( 0.0 ),
  pz_cluster( 0.0 ),
  pxx_cluster( 0.0 ),
  pyy_cluster( 0.0 ),
  pzz_cluster( 0.0 ),
  pxy_cluster( 0.0 ),
  pxz_cluster( 0.0 ),
  pyz_cluster( 0.0 ),
  volume_cluster( 0.0 )
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
  E_cluster = 0.0;
  px_cluster = 0.0;
  py_cluster = 0.0;
  pz_cluster = 0.0;
  pxx_cluster =  0.0;
  pyy_cluster =  0.0;
  pzz_cluster =  0.0;
  pxy_cluster =  0.0;
  pxz_cluster =  0.0;
  pyz_cluster =  0.0;
  volume_cluster = 0.0;
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
  E_cluster = 0.0;
  px_cluster = 0.0;
  py_cluster = 0.0;
  pz_cluster = 0.0;
  pxx_cluster =  0.0;
  pyy_cluster =  0.0;
  pzz_cluster =  0.0;
  pxy_cluster =  0.0;
  pxz_cluster =  0.0;
  pyz_cluster =  0.0;
  volume_cluster = 0.0;

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
  E_cluster += particleToAdd.Mom.E();
  px_cluster += particleToAdd.Mom.Px();
  py_cluster += particleToAdd.Mom.Py();
  pz_cluster += particleToAdd.Mom.Pz();
  pxx_cluster += particleToAdd.Mom.Px() * particleToAdd.Mom.Px() /  particleToAdd.Mom.E();
  pyy_cluster += particleToAdd.Mom.Py() * particleToAdd.Mom.Py() /  particleToAdd.Mom.E();
  pzz_cluster += particleToAdd.Mom.Pz() * particleToAdd.Mom.Pz() /  particleToAdd.Mom.E();
  pxy_cluster += particleToAdd.Mom.Px() * particleToAdd.Mom.Py() /  particleToAdd.Mom.E();
  pxz_cluster += particleToAdd.Mom.Px() * particleToAdd.Mom.Pz() /  particleToAdd.Mom.E();
  pyz_cluster += particleToAdd.Mom.Py() * particleToAdd.Mom.Pz() /  particleToAdd.Mom.E();

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

void cellContainer::getClusterInformation( double &T, VectorTXYZ &beta, const double minNumber, const double Ntest, const bool calcTviaEdens )
{
  const double number_cluster = number_gluon_cluster + number_quark_cluster;

  if( number_cluster >= minNumber )
  {
    // Attention: p_cluster is now normalized to E and becomes beta.
    beta = p_cluster.NormalizeToE();

    const double gamma = 1.0 / sqrt( 1 - beta.vec2() );

    if( !calcTviaEdens )
    {
      // T* as calculated in JHEP 0301 (2003) 030, eqs. (1.6), (1.7)
      const double Cf = ( 4.0 / 3.0 );
      const double Ca = 3.0;
      const double I = 0.5 * ( Cf * number_quark_cluster + Ca * number_gluon_cluster ) / gamma; // unitless, gamma comes from boost of volume in partice density
      const double J = Cf * inverseE_quarks_cluster + Ca * inverseE_gluons_cluster; // 1 / GeV
     
      T = I / J;
    }
    else
    {
      // Temperature calculated via the quartic root of the energy density
      
      const double gG = 2 * ( pow( ns_casc::Ncolor, 2 ) - 1 ); //16
      const double gQ = 2.0 * ns_casc::Ncolor * Particle::N_light_flavor; //18

      const double T00 = E_cluster / ( Ntest * volume_cluster );
      const double T01 = px_cluster / ( Ntest * volume_cluster );
      const double T02 = py_cluster / ( Ntest * volume_cluster );
      const double T03 = pz_cluster / ( Ntest * volume_cluster );
      const double T11 = pxx_cluster / ( Ntest * volume_cluster );
      const double T22 = pyy_cluster / ( Ntest * volume_cluster );
      const double T33 = pzz_cluster / ( Ntest * volume_cluster );
      const double T12 = pxy_cluster / ( Ntest * volume_cluster );
      const double T13 = pxz_cluster / ( Ntest * volume_cluster );
      const double T23 = pyz_cluster / ( Ntest * volume_cluster );    
    
      const double L00 = gamma;
      const double L01 = -gamma * beta.X();
      const double L02 = -gamma * beta.Y();
      const double L03 = -gamma * beta.Z();
      
      //LRF energy density
      const double energyDensity = L00*L00*T00       +   L01*L01*T11     +   L02*L02*T22   +   L03*L03*T33
                                  + 2*L00*L01*T01   +   2*L00*L02*T02   +   2*L00*L03*T03
                                  + 2*L01*L02*T12   +   2*L01*L03*T13   +   2*L02*L03*T23;//GeV/fm^3
      
      T = sqrt( sqrt( energyDensity * pow( 0.197, 3.0 ) * M_PI * M_PI / ( 3.0 * ( gG + 2.0 * gQ ) ) ) );
    }
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
