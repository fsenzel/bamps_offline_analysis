//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#ifndef MFPFORHEAVYIONCOLLISION_H
#define MFPFORHEAVYIONCOLLISION_H

#include <vector>

#include "mfp_data.h"
#include "particle.h"
#include "configuration.h"
#include "unit_enum.h"


class mfpForHeavyIonCollision
{
public:
    mfpForHeavyIonCollision( config* const _c );
    ~mfpForHeavyIonCollision();

    double getMeanFreePath( const double _E, const FLAVOR_TYPE _F, const double _T, const double _ngTest, const double _nqTest, const UNIT_TYPE _unit ) const;
    double getMeanFreePathForPhotonBremsstrahlung(const double _E, const FLAVOR_TYPE _F, const double _T, const double _ngTest, const double _nqTest, const UNIT_TYPE _unit) const;
    void loadData();


private:
    double fugacityDependence( const double _gluonFugacity, const double _quarkFugacity ) const;
    double getGluonDensity( const double T, const double fugacity = 1 ) const ;
    double getQuarkDensity( const double T, const double fugacity = 1 ) const;

    int getStartIndexForTemperatureValues( const double _T, const unsigned int _nValues ) const;

    /** @brief interpolation routine */
    void polint( const double xa[], const double ya[], const int n, const double x, double *y, double *dy ) const;

    std::vector<mfp_data> mfpData;
    std::vector<double> temperaturesForMfpData;

    double fitParameterGluon;
    double fitParameterQuark;
    double fugacityDependenceEvaluatedAtFugacityOne;

    bool data_loaded;

    config * const theConfig;
};


/** @brief exception class for handling unexpected behaviour when reading the files containing the MPF data */
class eMFP_heavy_ion_error : public std::runtime_error
{
public:
    explicit eMFP_heavy_ion_error( const std::string& what ) : std::runtime_error( what ) {};

    virtual ~eMFP_heavy_ion_error() throw() {};
};

#endif // MFPFORHEAVYIONCOLLISION_H
