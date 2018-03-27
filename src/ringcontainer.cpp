//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------
//---------------------------------------------

// at revision 912, this is identical to full/branches/vector4D/src/ringcontainer.cpp

#include <math.h>

#include "ringcontainer.h"
#include "particle.h"
#include "configuration.h"
#include "FPT_compare.h"


ringContainer::ringContainer() :
    minRadius( 0 ), maxRadius( 0 ), deltaR( 0 ),
    rates(),
    averagesPrepared( false ),
    numberOfParticles( 0 ),
    numberOfGluons( 0 ), numberOfQuarks( 0 ),
    numberOfActiveParticles( 0 ),
    md2g( 0 ), md2q( 0 ),
    v_x( 0 ), v_y( 0 ), v_z( 0 ), v_r( 0 ),
    E( 0 ), inverseE_gluons( 0 ), inverseE_quarks( 0 ),
    p_z( 0 ), p_r( 0 ), p_t( 0 ),
    pr2_over_E( 0 ), pz2_over_E( 0 ), pr_pz_over_E( 0 ),
    gamma( 0 ), energyDensity( 0 ), particleDensity( 0 ),
    gluonDensity( 0 ), quarkDensity( 0 ), volume( 0 ),
    numberOfCollectedRateObjects( 0 )
{
    rates.normalizeRates();
}



ringContainer::ringContainer( const double _minR, const double _maxR ) :
    minRadius( _minR ), maxRadius( _maxR ), deltaR( _maxR - _minR ),
    rates(),
    averagesPrepared( false ),
    numberOfParticles( 0 ),
    numberOfGluons( 0 ), numberOfQuarks( 0 ),
    numberOfActiveParticles( 0 ),
    md2g( 0 ), md2q( 0 ),
    v_x( 0 ), v_y( 0 ), v_z( 0 ), v_r( 0 ),
    E( 0 ), inverseE_gluons( 0 ), inverseE_quarks( 0 ),
    p_z( 0 ), p_r( 0 ), p_t( 0 ),
    pr2_over_E( 0 ), pz2_over_E( 0 ), pr_pz_over_E( 0 ),
    gamma( 0 ), energyDensity( 0 ), particleDensity( 0 ),
    gluonDensity( 0 ), quarkDensity( 0 ), volume( 0 ),
    numberOfCollectedRateObjects( 0 )
{
    rates.normalizeRates();
}


void ringContainer::clear()
{
    averagesPrepared = false;
    numberOfParticles = 0;
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



void ringContainer::addParticle( const Particle& _particle )
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

    if ( !_particle.free || FPT_COMP_G( _particle.rate22, 0 ) )
    {
        ++numberOfActiveParticles;
    }


    double xt = _particle.Pos.Perp();
    double oneE = 1 / _particle.Mom.E();


    E += _particle.Mom.E();
    if ( _particle.FLAVOR == gluon )
    {
        inverseE_gluons += oneE;
    }
    else
    {
        inverseE_quarks += oneE;
    }
    v_z += _particle.Mom.Pz() * oneE;

    double pr;
    if ( xt < 1.0e-5 )
    {
        pr = _particle.Mom.Pt();
        v_x += _particle.Mom.Px() * oneE;
        v_y += _particle.Mom.Py() * oneE;
    }
    else
    {
        //Multiply Unit vector in the direction pointing from 0 to the position of the particle scalar with the momentum to project the transverse component out
        double h_x =  _particle.Mom.Px() * _particle.Pos.X() / xt;
        double h_y =  _particle.Mom.Py() * _particle.Pos.Y() / xt;
        pr = h_x + h_y;
        v_x += h_x * oneE;
        v_y += h_y * oneE;
    }
    v_r += pr * oneE;

    p_r += pr;
    p_z += _particle.Mom.Pz();
    p_t += _particle.Mom.Perp();
    pr2_over_E += pow( pr, 2 ) * oneE;
    pz2_over_E += pow( _particle.Mom.Pz(), 2 ) * oneE;
    pr_pz_over_E += pr * _particle.Mom.Pz() * oneE;
}



void ringContainer::addParticleInFormGeom( const Particle& _particle, const double _time )
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

    double Eold = _particle.Old.E();
    double oneE = 1 / Eold;

    double cc = ( _time - _particle.Pos.T() ) * oneE;

    double xx = _particle.Pos.X() + _particle.Old.Px() * cc;
    double yy = _particle.Pos.Y() + _particle.Old.Py() * cc;
    double xt = sqrt( pow( xx, 2 ) + pow( yy, 2 ) );

    E += Eold;
    if ( _particle.FLAVOR == gluon )
    {
        inverseE_gluons += oneE;
    }
    else
    {
        inverseE_quarks += oneE;
    }
    v_z += _particle.Old.Pz() * oneE;

    double pr;
    if ( xt < 1.0e-5 )
    {
        pr = _particle.Old.Pt();
        v_x += _particle.Old.Px() * oneE;
        v_y += _particle.Old.Py() * oneE;
    }
    else
    {
        double h_x = _particle.Old.Px() * xx / xt;
        double h_y = _particle.Old.Py() * yy / xt;
        pr = h_x + h_y;
        v_x += h_x * oneE;
        v_y += h_y * oneE;
    }
    v_r += pr * oneE;

    p_r += pr;
    p_z += _particle.Mom.Pz();
    p_t += _particle.Mom.Perp();
    pr2_over_E += pow( pr, 2 ) * oneE;
    pz2_over_E += pow( _particle.Old.Pz(), 2 ) * oneE;
    pr_pz_over_E += pr * _particle.Old.Pz() * oneE;
}







void ringContainer::addRates( const Particle& _particle )
{
    rates.addParticleBasedRates( _particle, GeV );;
    ++numberOfCollectedRateObjects;
}


void ringContainer::relocate( const double _minR, const double _maxR )
{
    minRadius = _minR;
    maxRadius = _maxR;
    deltaR = _maxR - _minR;
    clear();
}

//WARNING
//This is the ring volume in LAB frame. When this is used, be aware that the computed quantity is a LAB-frame quantity (if no other boost-parameter is added.)
double ringContainer::getVolume( const double _dz ) const
{
    return ( M_PI * ( pow( maxRadius, 2 ) - pow( minRadius, 2 ) ) * _dz );   //fm^3
}



double ringContainer::getEnergyDensity() const
{
    if ( !averagesPrepared )
    {
        std::string errMsg = "ringContainer: averaged quantity requested without prior call to prepareAverages()";
        throw eRingContainer_error( errMsg );
    }

    return energyDensity;
}




double ringContainer::getParticleDensity() const
{
    if ( !averagesPrepared )
    {
        std::string errMsg = "ringContainer: averaged quantity requested without prior call to prepareAverages()";
        throw eRingContainer_error( errMsg );
    }

    return particleDensity;
}



double ringContainer::getGluonDensity() const
{
    if ( !averagesPrepared )
    {
        std::string errMsg = "ringContainer: averaged quantity requested without prior call to prepareAverages()";
        throw eRingContainer_error( errMsg );
    }

    return gluonDensity;
}



double ringContainer::getQuarkDensity() const
{
    if ( !averagesPrepared )
    {
        std::string errMsg = "ringContainer: averaged quantity requested without prior call to prepareAverages()";
        throw eRingContainer_error( errMsg );
    }

    return quarkDensity;
}



double ringContainer::getGamma() const
{
    if ( !averagesPrepared )
    {
        std::string errMsg = "ringContainer: averaged quantity requested without prior call to prepareAverages()";
        throw eRingContainer_error( errMsg );
    }

    return gamma;
}



double ringContainer::getAveraged_md2g() const
{
    if ( !averagesPrepared )
    {
        std::string errMsg = "ringContainer: averaged quantity requested without prior call to prepareAverages()";
        throw eRingContainer_error( errMsg );
    }

    return md2g;
}



double ringContainer::getAveraged_md2q() const
{
    if ( !averagesPrepared )
    {
        std::string errMsg = "ringContainer: averaged quantity requested without prior call to prepareAverages()";
        throw eRingContainer_error( errMsg );
    }

    return md2q;
}



void ringContainer::prepareAverages( const double _dz, const int _Ntest )
{
    volume = getVolume( _dz );
    averagesPrepared = true;

    if ( numberOfCollectedRateObjects > 0 )
    {
        rates.prepareParticleBasedAverages();
    }

    if ( numberOfParticles > 0 )
    {
        double gG = 2 * ( pow( ns_casc::Ncolor, 2 ) - 1 );
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
        gamma = 1 / sqrt( 1.0 - pow( getAveraged_v_z(), 2 ) - pow( getAveraged_v_r(), 2 ) );
        
        md2g = pow( 0.197, 3 ) * 16 * M_PI / (volume/gamma) * ( ns_casc::Ncolor * invEg + Particle::N_light_flavor * invEq );
        md2q = pow( 0.197, 3 ) * 2 * M_PI / (volume/gamma) * 8.0 / 3.0 * ( invEg + invEq );

        //all w/o alpha_s
        //md2g = 16 * pi * int (d^3p/(2pi)^3/E) (n_c f_g + n_f f_q)
        //     = 16 * M_PI * 1/(V*gamma) ( ns_casc::Ncolor * invEg + Particle::N_light_flavor * invEq )
        
       
        // Boost-Gamma factor correct?
        particleDensity = numberOfParticles / ( _Ntest * volume * gamma );
        gluonDensity = numberOfGluons / ( _Ntest * volume * gamma );
        quarkDensity = numberOfQuarks / ( _Ntest * volume * gamma );

        energyDensity = ( E - ( 2 * getAveraged_v_r() * p_r ) - ( 2 * getAveraged_v_z() * p_z )
                          + ( pow( getAveraged_v_r(), 2 ) * pr2_over_E ) + ( pow( getAveraged_v_z(), 2 ) * pz2_over_E )
                          + ( 2 * getAveraged_v_r() * getAveraged_v_z() * pr_pz_over_E ) ) * pow( gamma, 2 ) / ( _Ntest * volume );                //GeV/fm^3
                          
        //cout << "Quark fugacity " <<  quarkDensity/(12/pow(M_PI,2.0)*pow((energyDensity / (3 * particleDensity)),3.0))  << endl;                                        
        //std::cout << "sqrt(Debye-Mass^2*pi/8/(Nc+Nf)): "<< sqrt(md2g*M_PI/8.0/(Particle::N_light_flavor + ns_casc::Ncolor)) << "\t"<< (energyDensity / (3 * particleDensity)) << "\t" << sqrt(md2g*M_PI/8.0/(Particle::N_light_flavor + ns_casc::Ncolor)) / ( (energyDensity / (3 * particleDensity)) ) << "\t" << particleDensity/md2g << std::endl;                  
 
        // Obtain the effective temperature like in AMY, https://arxiv.org/abs/hep-ph/0209353v3, Eq. (1.7) and (1.6) and (A9)
        // Boost-Gamma factor canceled out
        double IGluon =     0.5 *  ( (numberOfGluons / ( _Ntest * volume) )  /  gG ) * pow( 0.197, 3 );   //GeV^3
        double IQuark =     0.5 *  ( (numberOfQuarks / ( _Ntest * volume) )  / ( ns_casc::Ncolor * Particle::N_light_flavor * 2 * 2 ) ) * pow( 0.197, 3 ); //GeV^3
            
        double JGluon =     pow( 0.197, 3 )/ (volume) * invEg; //GeV^2
        double JQuark =     pow( 0.197, 3 )/ (volume) * invEq; //GeV^2
        
        TStarGEV =   (IGluon + IQuark)/(JGluon + JQuark); 

        //std::cout << "Temp/GeV = " << TStarGEV <<  std::endl;
        
    }
}

/**
 * @brief Calculate effective temperature by e/3n = (second over first moment of f)
 * 
 * 
 */
double ringContainer::getEffectiveTemperature() const
{
    if ( !averagesPrepared )
    {
        std::string errMsg = "ringContainer: averaged quantity requested without prior call to prepareAverages()";
        throw eRingContainer_error( errMsg );
    }

    return energyDensity / (3 * particleDensity);
}


/**
 * @brief Calculate effective temperature by n/mD^2 = (first over zeroth moment of f)
 * 
 * 
 */
double ringContainer::getEffectiveTemperatureFirstOverZerothMoment()
{
  //This corresponds to https://arxiv.org/abs/hep-ph/0209353v3, Eq. (1.7) and (1.6) and (A9)
  // AMY_mathcal_I = in equil. similar to n     ~ T^3
  // AMY_mathcal_J = in equil. similar to md^2  ~ T^2
  
  return TStarGEV;
}



double ringContainer::transformEnergyToComovingFrame(VectorEPxPyPz & P) const
{
    if ( !averagesPrepared )
    {
        std::string errMsg = "ringContainer: averaged quantity requested without prior call to prepareAverages()";
        throw eRingContainer_error( errMsg );
    }

    double Edash = 0;
    if ( numberOfParticles > 0 )
    {
        Edash = gamma * ( P.E() - ( v_x * P.Px() + v_y * P.Py() + v_z * P.Pz() ) / numberOfParticles );
    }

    return Edash;
}


// kate: indent-mode cstyle; replace-tabs on; 
