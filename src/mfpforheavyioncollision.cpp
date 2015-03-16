#include "mfpforheavyioncollision.h"

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <math.h>

#include "mfp_data.h"
#include "unit_enum.h"

using namespace std;

mfpForHeavyIonCollision::mfpForHeavyIonCollision( config*const _c ) :
    fitParameterGluon(-0.569485),
    fitParameterQuark(-1.07081),
    data_loaded( false ),
    theConfig( _c )
{
    //just to be sure
    temperaturesForMfpData.clear();
    mfpData.clear();
}




mfpForHeavyIonCollision::~mfpForHeavyIonCollision()
{
    temperaturesForMfpData.clear();
    mfpData.clear();
}


void mfpForHeavyIonCollision::loadData()
{
    //write temperatures for which the mfp data should be read from tables into a vector
    // T in GeV
    temperaturesForMfpData.push_back( 0.1 );
    temperaturesForMfpData.push_back( 0.2 );
    temperaturesForMfpData.push_back( 0.3 );
    temperaturesForMfpData.push_back( 0.4 );
    temperaturesForMfpData.push_back( 0.6 );
    temperaturesForMfpData.push_back( 0.8 );
    temperaturesForMfpData.push_back( 1.0 );
    temperaturesForMfpData.push_back( 1.5 );

    //read temperature from the above created vector, initialize corresponding mfp_data object and push it into vector
    for ( unsigned int i = 0; i < temperaturesForMfpData.size(); i++ )
    {
        stringstream ssT;
        ssT << int( temperaturesForMfpData[i] * 1000 + 0.001 ); // convert from GeV to MeV and write to stringstream
        mfp_data tempMFP( ssT.str() );

        mfpData.push_back( tempMFP );
    }

    fugacityDependenceEvaluatedAtFugacityOne = fugacityDependence( 1, 1 );

    data_loaded = true;
}


double mfpForHeavyIonCollision::getMeanFreePath(const double _E, const FLAVOR_TYPE _F, const double _T, const double _ngTest, const double _nqTest, const UNIT_TYPE _unit ) const
{
    if( !data_loaded )
    {
        string errMsg = "MFP data not loaded yet.";
        throw eMFP_heavy_ion_error( errMsg );
    }

    const double gluonFugacity = _ngTest / ( getGluonDensity( _T ) );
    double quarkFugacity;
    if( Particle::N_light_flavor != 0 )
        quarkFugacity = _nqTest / ( getQuarkDensity( _T ) );
    else
        quarkFugacity = 0.0;

    const double scaleForFugacity = fugacityDependence( gluonFugacity, quarkFugacity ) / fugacityDependenceEvaluatedAtFugacityOne;

    const int nInterpolValues = 4;
    int startIndex = getStartIndexForTemperatureValues( _T, nInterpolValues );

    double xa[nInterpolValues+1], ya[nInterpolValues+1];

    for ( int i = 0; i < nInterpolValues; i++ )
    {
        xa[i+1] = temperaturesForMfpData[ startIndex + i ];
        ya[i+1] = mfpData[ startIndex + i ].getLambda( _E, _F );  //1/GeV
    }

    double interpolResult = 0, interpolError = 0;
    polint( xa, ya, nInterpolValues, _T, &interpolResult, &interpolError );

    double conversionFactor = 1;
    switch( _unit )
    {
    case GeV:
        conversionFactor = 1;
        break;
    case fm:
        conversionFactor = 0.197;
        break;
    }
    interpolResult *= conversionFactor;

    interpolResult *= scaleForFugacity;

    return interpolResult;
}



double mfpForHeavyIonCollision::fugacityDependence( const double _gluonFugacity, const double _quarkFugacity ) const
{
    double termGluon = sqrt( 2 * _gluonFugacity - fitParameterGluon * pow( _gluonFugacity, 2 ) );
    double termQuark = sqrt( 2 * _quarkFugacity - fitParameterQuark * pow( _quarkFugacity, 2 ) );
    return ( 1 / ( termGluon + termQuark ) );
}




double mfpForHeavyIonCollision::getGluonDensity( const double T, const double fugacity ) const
{
    return fugacity * 16 * pow( T, 3.0 ) / pow( M_PI, 2.0 ) / pow( 0.197, 3.0 );
}



double mfpForHeavyIonCollision::getQuarkDensity( const double T, const double fugacity ) const
{
    return fugacity * 12 * Particle::N_light_flavor * pow( T, 3.0 ) / pow( M_PI, 2.0 ) / pow( 0.197, 3.0 );
}



int mfpForHeavyIonCollision::getStartIndexForTemperatureValues(const double _T, const unsigned int _nValues) const
{
    if ( _nValues > temperaturesForMfpData.size() )
    {
        string errMsg = "requested number of interpolation points larger than temperature array size";
        throw eMFP_heavy_ion_error( errMsg );
    }

    unsigned int nn = 0;
    while ( temperaturesForMfpData[nn] < _T && nn < temperaturesForMfpData.size() )
    {
        ++nn;
    }

    nn = nn - _nValues / 2;
    if ( nn < 0 )
    {
        nn = 0;
    }
    else if ( nn + _nValues > temperaturesForMfpData.size() )
    {
        while ( nn + _nValues > temperaturesForMfpData.size() )
        {
            --nn;
        }
    }

    return nn;
}




/**
* Routine that linearly interpolates. Attention: Offset = 1 !
*
* @param[in] xa[] array of x-values
* @param[in] ya[] array of y-values
* @param[in] n number of x- and y-values used for the interpolation
* @param[in] x value at which the interpolated y should be evaluated
* @param[out] y result
* @param[out] dy some sort of error
*/
void mfpForHeavyIonCollision::polint(const double xa[], const double ya[], const int n, const double x, double* y, double* dy) const
{

    int ns = 1;
    double den, dif, dift, ho, hp, w;
    double *c, *d;

    dif = fabs( x - xa[1] );

    c = new double[n+1];
    d = new double[n+1];

    for ( int i = 1; i <= n; i++ )
    {
        if (( dift = fabs( x - xa[i] ) ) < dif )
        {
            ns = i;
            dif = dift;
        }
        c[i] = ya[i];
        d[i] = ya[i];
    }

    *y = ya[ns--];

    for ( int m = 1; m < n; m++ )
    {
        for ( int i = 1; i <= n - m; i++ )
        {
            ho = xa[i] - x;
            hp = xa[i+m] - x;
            w = c[i+1] - d[i];
            if (( den = ho - hp ) == 0.0 )
            {
                cout << "Error in routine polint" << endl;
                cout << x << endl;
                cout << xa[1] << "\t" << ya[1] << endl;
                cout << xa[2] << "\t" << ya[2] << endl;
                cout << xa[3] << "\t" << ya[3] << endl;
                cout << xa[4] << "\t" << ya[4] << endl;
            }
            den = w / den;
            d[i] = hp * den;
            c[i] = ho * den;
        }

        *y += ( *dy = ( 2 * ns < ( n - m ) ? c[ns+1] : d[ns--] ) );
    }

    delete[] c;
    delete[] d;
}
