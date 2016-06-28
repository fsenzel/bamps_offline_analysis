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
    lambdaScaled( 0 )
{
    SpecificMFP.assign(3,0);
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
    SpecificMFP.assign(3,0);

    averagesPrepared = false;
}


void cellContainer::resetStoredValues()
{
    SpecificMFP.assign(3,0);
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
    rates.clear();
    averagesPrepared = false;
}

/**
 * Routine to store (collect) a probability for specified 2->2 interaction process
 *
 * @param[in] _specificInteractionType Type of the rate to be stored (0: qq -> qq, 1: qqbar->qqbar + qqbar->qqbarDash, 2: qqbarDash->qqbarDash, qq'->qq')
 * @param[in] _P Probability to be stored
 */
void cellContainer::addSpecificRate(int _specificInteractionType, const double _P )
{
  rates.addSpecific(_specificInteractionType, _P);
}

void cellContainer::compute22MeanFreePath()
{
  SpecificMFP[0] = rates.getLambdaSpecific( 0, fm);
  SpecificMFP[1] = rates.getLambdaSpecific( 1, fm);
  SpecificMFP[2] = rates.getLambdaSpecific( 2, fm);
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

// kate: indent-mode cstyle; replace-tabs on; 
