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

#include "ringcontainer.h"
#include "particle.h"
#include "configuration.h"
#include "FPT_compare.h"


ringContainer::ringContainer() : minRadius ( 0 ), maxRadius ( 0 ), deltaR ( 0 ), numberOfParticles ( 0 ), numberOfGluons ( 0 ), numberOfQuarks ( 0 ),
  numberOfActiveParticles ( 0 ), md2g ( 0 ), md2q ( 0 ),
  v_x ( 0 ), v_y ( 0 ), v_z ( 0 ), v_r ( 0 ), inverseE_gluons ( 0 ), inverseE_quarks ( 0 ), E ( 0 ), p_z ( 0 ), p_t ( 0 ),
  gamma ( 0 ), particleDensity ( 0 ), gluonDensity ( 0 ), quarkDensity ( 0 ), energyDensity ( 0 ), averagesPrepared ( false ), volume ( 0 ),
  pr2_over_E ( 0 ), pz2_over_E ( 0 ), pr_pz_over_E ( 0 ), rates(), numberOfCollectedRateObjects ( 0 )
{
  rates.normalizeRates();
}



ringContainer::ringContainer ( const double _minR, const double _maxR ) : minRadius ( _minR ), maxRadius ( _maxR ), deltaR ( _maxR - _minR ),
  numberOfParticles ( 0 ), numberOfActiveParticles ( 0 ), numberOfGluons ( 0 ), numberOfQuarks ( 0 ), md2g ( 0 ), md2q ( 0 ),
  v_x ( 0 ), v_y ( 0 ), v_z ( 0 ), v_r ( 0 ), inverseE_gluons ( 0 ), inverseE_quarks ( 0 ), E ( 0 ), p_z ( 0 ), p_t ( 0 ),
  gamma ( 0 ), particleDensity ( 0 ), gluonDensity ( 0 ), quarkDensity ( 0 ), energyDensity ( 0 ), averagesPrepared ( false ), volume ( 0 ),
  pr2_over_E ( 0 ), pz2_over_E ( 0 ), pr_pz_over_E ( 0 ), rates(), numberOfCollectedRateObjects ( 0 )
{
  rates.normalizeRates();
}


void ringContainer::clear()
{
  averagesPrepared = false;
  numberOfParticles = 0;
  numberOfActiveParticles = 0;
  numberOfGluons = 0;
  numberOfQuarks = 0;
  particleDensity = 0;
  gluonDensity = 0;
  quarkDensity = 0;
  md2g = 0;
  md2q = 0;
  v_x = 0;
  v_y = 0;
  v_z = 0;
  v_r = 0;
  inverseE_gluons = 0;
  inverseE_quarks = 0;
  E = 0;
  p_r = 0;
  p_z = 0;
  p_t = 0;
  pr2_over_E = 0;
  pz2_over_E = 0;
  pr_pz_over_E = 0;

  rates.clear();
  rates.normalizeRates();
  numberOfCollectedRateObjects = 0;
}



void ringContainer::addParticle ( const ParticleOffline& _particle )
{
  ++numberOfParticles;
  if ( _particle.FLAVOR == gluon )
  {
    ++numberOfGluons;
  }
  else
  {
    ++numberOfQuarks;
  }

  if ( !_particle.free || FPT_COMP_G ( _particle.rate, 0 ) )
  {
    ++numberOfActiveParticles;
  }


  double xt = sqrt ( pow ( _particle.X, 2 ) + pow ( _particle.Y, 2 ) );

  E += _particle.E;
  if ( _particle.FLAVOR == gluon )
  {
    inverseE_gluons += 1 / _particle.E;
  }
  else
  {
    inverseE_quarks += 1 / _particle.E;
  }
  v_z += _particle.PZ / _particle.E;

  double pr;
  if ( xt < 1.0e-5 )
  {
    pr = sqrt ( pow ( _particle.PX, 2 ) + pow ( _particle.PY, 2 ) );
    v_x += _particle.PX / _particle.E;
    v_y += _particle.PY / _particle.E;
  }
  else
  {
    pr = ( _particle.PX * _particle.X + _particle.PY * _particle.Y ) / xt;
    v_x += _particle.PX * _particle.X / ( _particle.E * xt );
    v_y += _particle.PY * _particle.Y / ( _particle.E * xt );
  }
  v_r += pr / _particle.E;

  p_r += pr;
  p_z += _particle.PZ;
  p_t += sqrt ( _particle.PX * _particle.PX + _particle.PY * _particle.PY );
  pr2_over_E += pow ( pr, 2 ) / _particle.E;
  pz2_over_E += pow ( _particle.PZ, 2 ) / _particle.E;
  pr_pz_over_E += pr * _particle.PZ / _particle.E;
}



void ringContainer::addParticleInFormGeom ( const ParticleOffline& _particle, const double _time )
{
  ++numberOfParticles;
  if ( _particle.FLAVOR == gluon )
  {
    ++numberOfGluons;
  }
  else
  {
    ++numberOfQuarks;
  }

  double Eold = sqrt ( pow ( _particle.PXold, 2 ) + pow ( _particle.PYold, 2 ) + pow ( _particle.PZold, 2 ) + pow ( _particle.m, 2 ) );
  double cc = ( _time - _particle.T ) / Eold;
  double zz = _particle.Z + _particle.PZold * cc;

  double xx = _particle.X + _particle.PXold * cc;
  double yy = _particle.Y + _particle.PYold * cc;
  double xt = sqrt ( pow ( xx, 2 ) + pow ( yy, 2 ) );

  E += Eold;
  if ( _particle.FLAVOR == gluon )
  {
    inverseE_gluons += 1 / Eold;
  }
  else
  {
    inverseE_quarks += 1 / Eold;
  }
  v_z += _particle.PZold / Eold;

  double pr;
  if ( xt < 1.0e-5 )
  {
    pr = sqrt ( pow ( _particle.PXold, 2 ) + pow ( _particle.PYold, 2 ) );
    v_x += _particle.PXold / Eold;
    v_y += _particle.PYold / Eold;
  }
  else
  {
    pr = ( _particle.PXold * xx + _particle.PYold * yy ) / xt;
    v_x += _particle.PXold * xx / ( Eold * xt );
    v_y += _particle.PYold * yy / ( Eold * xt );
  }
  v_r += pr / Eold;

  p_r += pr;
  p_z += _particle.PZold;
  p_t += sqrt ( _particle.PXold * _particle.PXold + _particle.PYold * _particle.PYold );
  pr2_over_E += pow ( pr, 2 ) / Eold;
  pz2_over_E += pow ( _particle.PZold, 2 ) / Eold;
  pr_pz_over_E += pr * _particle.PZold / Eold;
}







void ringContainer::addRates ( const ParticleOffline& _particle )
{
  rates.addParticleBasedRates ( _particle, GeV );;
  ++numberOfCollectedRateObjects;
}


void ringContainer::relocate ( const double _minR, const double _maxR )
{
  minRadius = _minR;
  maxRadius = _maxR;
  deltaR = _maxR - _minR;
  clear();
}


double ringContainer::getVolume ( const double _dz ) const
{
  return ( M_PI * ( pow ( maxRadius, 2 ) - pow ( minRadius, 2 ) ) * _dz ); //fm^3
}



double ringContainer::getEnergyDensity() const
{
  if ( !averagesPrepared )
  {
    std::string errMsg = "ringContainer: averaged quantity requested without prior call to prepareAverages()";
    throw eRingContainer_error ( errMsg );
  }

  return energyDensity;
}




double ringContainer::getParticleDensity() const
{
  if ( !averagesPrepared )
  {
    std::string errMsg = "ringContainer: averaged quantity requested without prior call to prepareAverages()";
    throw eRingContainer_error ( errMsg );
  }

  return particleDensity;
}



double ringContainer::getGluonDensity() const
{
  if ( !averagesPrepared )
  {
    std::string errMsg = "ringContainer: averaged quantity requested without prior call to prepareAverages()";
    throw eRingContainer_error ( errMsg );
  }

  return gluonDensity;
}



double ringContainer::getQuarkDensity() const
{
  if ( !averagesPrepared )
  {
    std::string errMsg = "ringContainer: averaged quantity requested without prior call to prepareAverages()";
    throw eRingContainer_error ( errMsg );
  }

  return quarkDensity;
}



double ringContainer::getGamma() const
{
  if ( !averagesPrepared )
  {
    std::string errMsg = "ringContainer: averaged quantity requested without prior call to prepareAverages()";
    throw eRingContainer_error ( errMsg );
  }

  return gamma;
}



double ringContainer::getAveraged_md2g() const
{
  if ( !averagesPrepared )
  {
    std::string errMsg = "ringContainer: averaged quantity requested without prior call to prepareAverages()";
    throw eRingContainer_error ( errMsg );
  }

  return md2g;
}



double ringContainer::getAveraged_md2q() const
{
  if ( !averagesPrepared )
  {
    std::string errMsg = "ringContainer: averaged quantity requested without prior call to prepareAverages()";
    throw eRingContainer_error ( errMsg );
  }

  return md2q;
}



void ringContainer::prepareAverages ( const double _dz, const int _Ntest )
{
  volume = getVolume ( _dz );
  averagesPrepared = true;

  if ( numberOfCollectedRateObjects > 0 )
  {
    rates.prepareParticleBasedAverages();
  }

  if ( numberOfParticles > 0 )
  {
    double gG = 2 * ( pow ( ns_casc::Ncolor, 2 ) - 1 );
    double gQ = 2.0 * ns_casc::Ncolor * Particle::N_light_flavor;

    double invEg = inverseE_gluons / ( gG * _Ntest );
    double invEq;
    if ( Particle::N_light_flavor == 0 )
    {
      invEq = 0;
    }
    else
    {
      invEq = inverseE_quarks / ( 2.0 * gQ * _Ntest );
    }
    md2g = pow ( 0.197, 3 ) * 16 * M_PI / volume * ( ns_casc::Ncolor * invEg + Particle::N_light_flavor * invEq );
    md2q = pow ( 0.197, 3 ) * 2 * M_PI / volume * 8.0 / 3.0 * ( invEg + invEq );

    gamma = 1 / sqrt ( 1.0 - pow ( getAveraged_v_z(), 2 ) - pow ( getAveraged_v_r(), 2 ) );
    particleDensity = numberOfParticles / ( _Ntest * volume * gamma );
    gluonDensity = numberOfGluons / ( _Ntest * volume * gamma );
    quarkDensity = numberOfQuarks / ( _Ntest * volume * gamma );

    energyDensity = ( E - ( 2 * getAveraged_v_r() * p_r ) - ( 2 * getAveraged_v_z() * p_z )
                      + ( pow ( getAveraged_v_r(), 2 ) * pr2_over_E ) + ( pow ( getAveraged_v_z(), 2 ) * pz2_over_E )
                      + ( 2 * getAveraged_v_r() * getAveraged_v_z() * pr_pz_over_E ) ) * pow ( gamma, 2 ) / ( _Ntest * volume );               //GeV/fm^3
  }
}


double ringContainer::getEffectiveTemperature() const
{
  if ( !averagesPrepared )
  {
    std::string errMsg = "ringContainer: averaged quantity requested without prior call to prepareAverages()";
    throw eRingContainer_error ( errMsg );
  }

  return energyDensity / ( 3 * particleDensity );
}




double ringContainer::transformEnergyToComovingFrame ( double _P[4] ) const
{
  if ( !averagesPrepared )
  {
    std::string errMsg = "ringContainer: averaged quantity requested without prior call to prepareAverages()";
    throw eRingContainer_error ( errMsg );
  }

  double Edash = 0;
  Edash = gamma * ( _P[0] - ( getAveraged_v_x() * _P[1] + getAveraged_v_y() * _P[2] + getAveraged_v_z() * _P[3] ) );

  return Edash;
}


// kate: indent-mode cstyle; indent-width 2; replace-tabs on; 
