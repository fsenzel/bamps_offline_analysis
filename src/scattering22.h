//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------
#include "particle.h"

/** @file 
 * @brief definitions for the scattering22 class
 */

/** @author Oliver Fochler */

#ifndef SCATTERING22_H
#define SCATTERING22_H


//#define alpha_s(s) 12.0*pi/((33.0-2.0*Nflavor)*log(s/lambda2))
#define alpha_s(s) 0.3


/**
 * @brief Provides routines and data structures for 2->2 scatterings.
 *
 * The class scattering22 encapsulates all routines and data structures needed for 2->2 scattering processes. It provides methods
 * to compute the cross section for a given particle pair (#getXSection22), sample the outgoing momenta (#getMomenta22) and to
 * asigne these momenta to the outgoing particles (#setNewMomenta22).
 *
 * A scattering22 object can be used for more than one particle doublet, just call #setParameter with the appropriate parameters
 * before the evaluation of each particle pair.
 */
class scattering22
{
  public:
    /** @brief Standard constructor. When used, setParameter needs to be called prior to other methods! */
    scattering22();
    /** @brief Constructor taking the same arguments as setParamter */
    scattering22(const double P1arg[], const double P2arg[], const FLAVOR_TYPE F1_arg, const FLAVOR_TYPE F2_arg, 
                  const double s_arg, const double md2g_arg, const double md2q_arg);
    /** @brief standard destructor */
    ~scattering22();
    
    /** @brief Sets the internal parameter needed for a specific particle pair. */    
    void setParameter(const double P1arg[], const double P2arg[], const FLAVOR_TYPE F1_arg, const FLAVOR_TYPE F2_arg,
                      const double s_arg, const double md2g_arg, const double md2q_arg);
    
    /** @brief Returns the cross section for the given particle pair (set with #setParameter or in the constructor). Wrapper */
    double getXSection22() const { int temp; return getXSection22( temp ); }
    /** @brief Returns the cross section for the given particle pair (set with #setParameter or in the constructor). */
    double getXSection22( int& initialStateIndex ) const;  
    
    /** @brief Returns the ELASTIC cross section for the given particle pair (set with #setParameter or in the constructor). Wrapper */
    double getXSectionElastic() const { int temp; return getXSectionElastic( temp ); }
    /** @brief Returns the ELASTIC cross section for the given particle pair (set with #setParameter or in the constructor). */
    double getXSectionElastic( int& initialStateIndex ) const;  
    
    
    /** @brief computes the transport cross section */
    double getTransportXSection22() const;
          
    /** @brief samples new momenta */
    void getMomenta22( double& PT2, int& typ, FLAVOR_TYPE & F1arg, FLAVOR_TYPE & F2arg ) const;
    /** @brief samples new momenta for ELASTIC processes*/
    void getMomentaElastic( double& PT2, int& typ ) const;    
    /** @brief sample new momenta isotropically */
    void getMomenta22_isotropic(double& PT2, int& typ) const;
    /** @brief sets new momenta for outgoing particles */
    void setNewMomenta22(double P1[4], double P2[4], const double R1[4], const double R2[4], const double PT2);
    
  
  private:
    void rotation(const double[], const double[], const double[], double[], double[]) const;
    /** @brief static function for direct access to the g+g cross section without need for a scattering22 object */
    static double xSectionGluons(const FLAVOR_TYPE F1_arg, const FLAVOR_TYPE F2_arg, const double s_arg, const double md2g_arg, const double md2q_arg);
    
    /** @brief Boost velocity from lab framet to CMS of particle pair.*/
    double beta_vec[4];
    
    /** @brief Original momentum vector P1 in the lab frame. */
    double * P1;
    /** @brief Original momentum vector P2 in the lab frame. */  
    double * P2;
    
    /** @brief Momentum vector P1 boosted to the CMS of the colliding particle triplet.*/
    double * P1cm;
    /** @brief Momentum vector P2 boosted to the CMS of the colliding particle triplet.*/
    double * P2cm;
    
    /** @brief flavor of particle 1 */
    FLAVOR_TYPE F1;
    /** @brief flavor of particle 2 */
    FLAVOR_TYPE F2;
    
    /** @brief gluonic debye mass (squared) divided by alpha_s */
    double md2_gluon;
    /** @brief debye mass (squared) for quarks divided by alpha_s */
    double md2_quark;
    /** @brief mandelstam variable s for the given particle pair */
    double s;
};


#endif
