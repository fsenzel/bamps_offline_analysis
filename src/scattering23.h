//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL: file:///home/fochler/svnrepos/cascade_full/main/branches/restructure_branch/src/scattering23.h $
//$LastChangedDate: 2010-07-06 16:24:24 +0200 (Tue, 06 Jul 2010) $
//$LastChangedRevision: 115 $
//$LastChangedBy: fochler $
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------

/** @file 
* @brief Declarations for class integrand23.
*/

#ifndef SCATTERING23_H
#define SCATTERING23_H

#include "interpolation23.h"
#include "particle.h"

/** @brief Class encapsulating the calculation of 2->3 scatterings */
class scattering23
{
  public:
    /** @brief Basic constructor. scattering23::setParameter MUST be used when using this constructor. */
    scattering23(const interpolation23 * const theI23_arg);
    /** @brief Constructor setting all necessary parameters. */
    scattering23(const interpolation23 * const theI23_arg, const double vx, const double vy, const double vz, const double P1arg[], const double P2arg[], 
                 const FLAVOR_TYPE F1_arg, const FLAVOR_TYPE F2_arg,
                 const double sqrtS_arg, const double md2g_scaled_arg, const double lambda_scaled_arg);
    /** @brief Destructor */
    ~scattering23();
    
    /** @brief Sets the internal parameter needed for a specific particle doublet. */
    double setParameter(const double vx, const double vy, const double vz, const double P1arg[], const double P2arg[],
                        const FLAVOR_TYPE F1_arg, const FLAVOR_TYPE F2_arg,
                        const double sqrtS_arg, const double md2g_scaled_arg, const double lambda_scaled_arg);
    
                        
    /** @brief Returns the integral over the 2->3 matrix element. Wrapper */
    double getIntegral23() const { int temp; return getIntegral23( temp ); }
    /** @brief Returns the integral over the 2->3 matrix element. */
    double getIntegral23( int& initialStateIndex ) const;  
    
    /** @brief Samples new momenta.*/ 
    int getMomenta23( double& pt1, double& pt3, double& y, double& phi, double& pz1, int& typ, FLAVOR_TYPE& F1arg, FLAVOR_TYPE& F2arg) const;
    
    /** @brief Sets new momenta.*/
    void setNewMomenta23(double P1[4], double P2[4], double P3[4], const double R1[4], const double R2[4],
                     const double PT1, const double PT3,const double y3, const double phi, const double PZ1);
    /** @brief Utility routine giving access to cos_theta. */
    double getCosTheta() const {return cos_theta;}
    
    
  private:
    /** @brief Computes the maximal range in rapidity y */
    bool maximal_Y_range_23(double& ymin, double& ymax) const;
    void rotation(const double[], const double[], const double[], double[], double[]) const;
    
    /** @brief Velocity with which the computational cell moves (if any). */
    double cell_beta_vec[4];
    /** @brief Boost velocity from lab frame to CMS. */
    double beta_vec[4];
    /** @brief absolute value of #beta_vec */
    double beta; //absolute value of beta_vec
    
    /** @brief cos(theta) where theta is the angle between the boost velocity vector and the axis of the CMS system */
    double cos_theta;
    /** @brief angle between the boost velocity vector and the axis of the CMS system */
    double theta;
    
    /** @brief Original momentum of particle 1 in lab frame. */
    double * P1;  
    /** @brief Original momentum of particle 2 in lab frame. */
    double * P2;
    
    /** @brief Momentum of particle 1 in rest frame of the computational cell. */
    double * P1cell;  
    /** @brief Momentum of particle 2 in rest frame of the computational cell. */
    double * P2cell;
    
    /** @brief Momentum of particle 1 in CM frame of the scattering. */
    double * P1cm;
    /** @brief Momentum of particle 2 in CM frame of the scattering. */
    double * P2cm;
    
    /** @brief flavor of particle 1 */
    FLAVOR_TYPE F1;
    /** @brief flavor of particle 2 */
    FLAVOR_TYPE F2;
    
    /** @brief Debye mass squared, scaled by 1/s */
    double md2g_scaled;
    /** @brief Center of mass energy, sqrt(s) */
    double sqrtS;
    /** @brief Mean free path, scaled by sqrt(s) */
    double lambda_scaled;
    
    /** @brief Pointer to an #interpolation23 object that is used for interpolating the integrated 2->3 matrix element */
    const interpolation23 * const theI23;
};


/** @brief exception class for handling unexpected critical behaviour within 2->3 routines  */
class eScatt23_error : public std::runtime_error
{
  public:
    explicit eScatt23_error(const std::string& what) : std::runtime_error(what) {};

    virtual ~eScatt23_error() throw() {};
};


#endif
