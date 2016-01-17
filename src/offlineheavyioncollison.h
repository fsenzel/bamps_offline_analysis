//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


#ifndef OFFLINEHEAVYIONCOLLISON_H
#define OFFLINEHEAVYIONCOLLISON_H

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <vector>


#include "analysis.h"
#include "cellcontainer.h"
#include "configuration.h"
#include "interpolation22.h"
#include "interpolation23.h"
#include "mfpforheavyioncollision.h"
#include "offlineoutput.h"
#include "ringstructure.h"
#include "scattering22.h"
#include "scattering23.h"
#include "scattering32.h"
#include "woodsaxon.h"


class offlineHeavyIonCollision
{
public:
    /**
     * @brief Constructor
     */
    offlineHeavyIonCollision( config* const _config, offlineOutputInterface* const _offlineInterface, analysis* const _analysis );

    /**
     * @brief Destructor
     */
    ~offlineHeavyIonCollision();

    /**
     * @brief Initialize particle momenta, positions, etc.
     */
    void initialize();

    /**
     * @brief The main framework of the simulation
     */
    void mainFramework();

    double evolveMedium( const double evolveToTime, bool& _endOfDataFiles );


private:
    /**
     * @brief A pointer to the global config object that encapsulates
     * user defined settings etc.
     */
    config* const theConfig;


    analysis* const theAnalysis;

    /** @brief  interpolation23 object that provides access to tabulated values for the cross section of all 2->3 processes with running coupling */
    interpolation23 theI23_massless;
    interpolation23 theI23_charm_m1;
    interpolation23 theI23_charm_m2;
    interpolation23 theI23_bottom_m1;
    interpolation23 theI23_bottom_m2;
    interpolation23 theI23_photons;

    /** @brief  interpolation22 object that provides access to tabulated values for the cross section of all 2->2 processes with running coupling */
    interpolation22 theI22;

    //---- parameters that are copied from the config object for convenience ----
    /** @brief Runtime to simulate [in fm/c] */
    double stop;
    /** @brief Mass number of nucleus A */
    double A;
    /** @brief Atomic number, i.e. number of protons, of nucleus A */
    double Aatomic;
    /** @brief Mass number of nucleus A */
    double B;
    /** @brief Atomic number, i.e. number of protons, of nucleus B */
    double Batomic;
    /** @brief Center of momentum energy per NN pair [GeV] */
    double sqrtS;
    /** @brief Lower pT cutoff [GeV] for mini-jet initial conditions, see config::initialStateType */
    double P0;
    /** @brief Impact parameter [in fm] */
    double Bimp;
    /** @brief Number of test particles per real particle */
    int testpartcl;
    //--------------------------------------------------------------------------

    /**
     * @brief Total Number of produced photons
     *
     */
    int totalPhotonNumber;
    int numberB, numberA;
    double averageS_smallv2;
    double averageS_Bigv2;

    
    
    /**
     * @brief Accumulated number of sampling errors in
     * scattering23::getMomenta23
     *
     * Can give hints on the validity of the envelope function. Should
     * be "small".
     */
    int nGet23Errors;

    /**
     * @brief Accumulated number of sampling errors in
     * scattering32::getMomenta32
     *
     * Only needed for rejection sampling, 0 for Metropolis sampling per
     * construction. For rejection sampling: can give hints on the
     * validity of the envelope function. Should be "small".
     */
    int nGet32Errors;

    /**
     * @brief The rings that are used for averaging several variables
     * (Debye masses, rates, etc.)
     */
    ringStructure rings;

    /** @brief Pointer to the interface for the output of data needed for later offline reconstruction */
    offlineOutputInterface* offlineInterface;

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


    void cell_ID( double time );

    void scattering( const double nexttime, bool& again );

    void scatterEdgeParticles( std::list< int >& _offlineParticleList, std::list< int >& _addedParticleList, const double nexttime );


    void scatt2223_offlineWithAddedParticles( cellContainer& _cell, std::vector< int >& _allParticlesList, std::vector< int >& _gluonList,
            cellContainer& _cellAdded, std::vector< int >& _allParticlesListAdded, std::vector< int >& _gluonListAdded,
            const double scaleFactor, bool& again, const double nexttime, analysisRingStructure& _analysisRings );

    void scatt22_amongAddedParticles( cellContainer& _cellAdded, std::vector< int >& _allParticlesListAdded, const double scaleFactor, bool& again, const double nexttime );

    void scatt22_amongBackgroundParticles_photons( cellContainer& _cellAdded, std::vector< int >& _allParticlesListAdded, const double scaleFactor, bool& again, const double nexttime );
    void scatt23_amongBackgroundParticles_photons( cellContainer& _cellAdded, std::vector< int >& _allParticlesListAdded, const double scaleFactor, bool& again, const double nexttime );

    void scatt32_offlineWithAddedParticles( cellContainer& _cell, std::vector< int >& _allParticlesList, std::vector< int >& _gluonList,
                                            cellContainer& _cellAdded, std::vector< int >& _allParticlesListAdded, std::vector< int >& _gluonListAdded,
                                            int& n32, bool& again, const double nexttime );

    int scatt32_offlineWithAddedParticles_utility( scattering32& scatt32_obj, std::list< int >& _cellMembersAdded, std::vector< int >& _allParticlesListAdded, std::vector< int >& _gluonListAdded, const int iscat, const int jscat, const int kscat, int& n32, const double nexttime  );
    int scatt23_offlineWithAddedParticles_utility( scattering23& scatt23_obj, cellContainer& _cell, int iscat, const int jscat, bool& again, const double nexttime );
    void scatt22_offlineWithAddedParticles_utility( scattering22& scatt22_obj, std::list< int >& _cellMembersAdded, std::vector< int >& _allParticlesListAdded, const int iscat, const int jscat, int& typ, const double nexttime );
    void scatt22_amongAddedParticles_utility( scattering22& scatt22_obj, std::list< int >& _cellMembersAdded, std::vector< int >& _allParticlesListAdded, const int iscat, const int jscat, int& typ, const double nexttime );
    void scatt22_amongBackgroundParticles_photons_utility( scattering22& scatt22_obj, std::list< int >& _cellMembersAdded, std::vector< int >& _allParticlesListAdded, const int iscat, const int jscat, int& typ, const double nexttime );
    void scatt23_amongBackgroundParticles_photons_utility( scattering23& scatt23_obj, std::list< int >& _cellMembersAdded, std::vector< int >& _allParticlesListAdded, const int iscat, const int jscat, int& typ, const double nexttime );

    /** @brief Goes through all added particles and decays J/psi if the temperature of the surrounding medium is larger than the dissocation temperature */
    void jpsi_dissociation_td( const double time );

    double iterateMFP( std::vector< int >& _allParticlesList, std::vector< int >& _gluonList, const int jetID, const double dt, const double dv );
    double iterate_mfp_bisection( std::vector< int >& _allParticlesList, std::vector< int >& _gluonList, const int jetID, const double dt, const double dv, const double lambda_old );

    /**
     * @brief A pointer to a \p mfpForHeavyIonCollision object that
     * provides interpolated mean free path data for added particles
     *
     * This is only needed when \sa config::jetMfpComputationSwitch is set
     * to \sa JET_MFP_COMPUTATION_TYPE::computeMfpInterpolation or JET_MFP_COMPUTATION_TYPE::thermalMfpGluon.
     */
    mfpForHeavyIonCollision theMFP;

    void removeDeadParticles();

    int binomial( const int N, const int k ) const;

    double addVelocities( const double vx_1, const double vy_1, const double vz_1, const double vx_2, const double vy_2, const double vz_2 );
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
