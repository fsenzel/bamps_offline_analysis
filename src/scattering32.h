//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL: file:///home/fochler/svnrepos/cascade_full/main/branches/restructure_branch/src/scattering32.h $
//$LastChangedDate: 2010-07-06 16:24:24 +0200 (Tue, 06 Jul 2010) $
//$LastChangedRevision: 115 $
//$LastChangedBy: fochler $
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------

/** @file 
 * @brief Declarations for the class scattering32.
 */


#ifndef SCATTERING32_H
#define SCATTERING32_H

#include "integrand32.h"
#include "particle.h"


using std::fstream;


/**
 * @brief Provides routines and data structures for 3->2 scatterings.
 * @author Oliver Fochler
 *
 * The class scattering32 encapsulates all routines and data structures needed for 3->2 scattering processes. It provides methods to
 * calculate the integral over the matrix element I32, #getIntegral32_vegas (using numerical integration) or #getIntegral32_fast 
 * (using an estimation procedure). Sampling of new momenta after a 3->2 collision is provided by #getMomenta32 and #setNewMomenta32
 * then sets the particle momenta to the new values.
 *
 * A scattering32 object can be used for more than one particle triplet, just call #setParameter with the appropriate parameters
 * before the evaluation of each triplet.
 */
class scattering32
{
  public:
    /** @brief Standard constructor. When used, setParameter needs to be called prior to other methods! */
    scattering32();
    /** @brief Constructor taking the same arguments as setParamter */
    scattering32(const double vx, const double vy, const double vz,
                 const double P1arg[], const double P2arg[], const double P3arg[],
                 const FLAVOR_TYPE F1_arg, const FLAVOR_TYPE F2_arg, const FLAVOR_TYPE F3_arg,
                 const double sqrtS_arg, const double md2g_scaled_arg, const double lambda_scaled_arg);
    /** @brief standard destructor */
    ~scattering32();
    
    /** @brief Sets the internal parameter needed for a specific particle triplet. */
    double setParameter(const double vx, const double vy, const double vz, 
                        const double P1arg[], const double P2arg[], const double P3arg[],
                        const FLAVOR_TYPE F1_arg, const FLAVOR_TYPE F2_arg, const FLAVOR_TYPE F3_arg,
                        const double sqrtS_arg, const double md2g_scaled_arg, const double lambda_scaled_arg);
    
    
    /** @brief Returns I32 using numerical integration. */
    double getIntegral32_vegas( int& initialStateIndex );
    /** @brief Returns I32 using numerical integration. Wrapper.*/
    double getIntegral32_vegas() { int temp; return getIntegral32_vegas( temp ); }
    /** @brief Returns I32 using an estimation procedure. */
    double getIntegral32_fast( int& initialStateIndex );
    /** @brief Returns I32 using an estimation procedure. Wrapper. */    
    double getIntegral32_fast() { int temp; return getIntegral32_fast( temp ); }
    
    /** @brief Samples new momenta. Wrapper method.*/ 
    int getMomenta32( double& u, double& phi, const double picked_ratio, int& typ, FLAVOR_TYPE& F1arg, FLAVOR_TYPE& F2arg);   
    /** @brief Sets new momenta.*/
    void setNewMomenta32(double P1_arg[], double P2_arg[], const double u, const double phi) const;
    
    
    /** @brief Returns the number indicating the actual order of P1, P2, P3. */
    int getOrder() const;    
    /** @brief Returns #total_ratio*/
    double getRatio() const;    
    /** @brief Returns #collision*/
    bool getCollisionStatus() const;
        
    
  private:
    /** @brief Transforms the boost vector #beta_vec to the reference frame given by P1_arg and P3_arg*/
    void rotate_beta(const double P1_arg[], const double P3_arg[], double beta_new[]);  
    
    /** @brief Get new momenta via the rejection method. Called by wrapper getMomenta32(double& u, double& phi)*/
    int getMomenta32_rejection(double& u, double& phi) const;
    /** @brief Get new momenta via the metropolis algorithm. Called by wrapper getMomenta32(double& u, double& phi)*/
    int getMomenta32_metropolis(double& u_arg, double& phi_arg) const;
    /** @brief Get new momenta when using estimation method and #direct_estimate = false. Called by wrapper getMomenta32(double& u, double& phi)*/
    int getMomenta32_estimate_rejection(double& u, double& phi, const double picked_ratio);
    
    /** @brief Returns the matrix element at given u and phi*/
    double getMatrixElement(const double u, const double phi) const;
    
    /** @brief Returns the estimated I32.*/
    double get_I32estimate(const double E1, const double E3, const double cosgamma_arg) const;
    /** @brief Returns the ratio of the "real" I32 to the estimated I32.*/
    double get_I32estimate_ratio(const double E1, const double E3, const double cosgamma_arg, const int N_arg);
    /** @brief Utility method needed by #get_I32estimate.*/
    double H(const double a, const double u) const;
    
    /** @brief Function object representing the integrand.
     *
     * integrand32 overloads the () operator, can be called like a normal function.
     */
    integrand32 theIntegrand;
    
    /** @brief Collective velocity of the cell. */
    double cell_beta_vec[4];
    /** @brief Boost velocity from cell rest frame to CMS of particle triplet.*/
    double beta_vec[4];
    /** @brief Absolute value of #beta_vec */
    double beta_abs;
    /**
     * @brief Boost vector #beta_vec rotated into possible reference frames given by p1 and p3.
     *         
     * Array (7 entries but index 0 is not in use) of vectors for holding the vector beta (#beta_vec) rotated 
     * into the reference frame given by p1 and p3 (or rather: the 6 possible choices for p1 and p3) - see notes.
     */
    double rotated_beta[7][4];
    
    /** @brief Original momentum vector P1 in the lab frame. */
    double * P1;
    /** @brief Original momentum vector P2 in the lab frame. */  
    double * P2;
    /** @brief Original momentum vector P3 in the lab frame. */
    double * P3;
    
    /** @brief Momentum vector P1 boosted to the cell frame.*/
    double * P1cell;
    /** @brief Momentum vector P2 boosted to the cell frame.*/  
    double * P2cell;
    /** @brief Momentum vector P3 boosted to the cell frame.*/
    double * P3cell;
    
    /** @brief Momentum vector P1 boosted to the CMS of the colliding particle triplet.*/
    double * P1cm;
    /** @brief Momentum vector P2 boosted to the CMS of the colliding particle triplet.*/
    double * P2cm;
    /** @brief Momentum vector P3 boosted to the CMS of the colliding particle triplet.*/
    double * P3cm;
    
    /** @brief flavor of particle 1 */
    FLAVOR_TYPE F1;
    /** @brief flavor of particle 2 */
    FLAVOR_TYPE F2;
    /** @brief flavor of particle 3 */
    FLAVOR_TYPE F3;
    
    
    /** @brief Debye mass squared, scaled by 1/s */
    double md2g_scaled;
    /** @brief Center of mass energy, sqrt(s) */
    double sqrtS;
    /** @brief Mean free path, scaled by sqrt(s) */
    double lambda_scaled;
    
    /**
     * @brief Determines whether the estimation procedure affects the sampling of new momenta or not.
     *
     * Default: true
     *
     * When using the estimation procedure #getIntegral32_fast for calculating I32, the ratio of the estimated result to the
     * real result is estimated by means of #get_I32estimate_ratio. New momenta are sampled according to the full matrix element
     * as usual.
     *
     *
     * Optional: false
     *
     * When set to false, the return value of #getIntegral32_fast is the estimated value for I32 without any correction by the
     * ratio to the expected true result. Momenta are then sampled according to the estimated matrix element and an additional
     * decision whether the collision takes place or not is made in #getMomenta32_estimate_rejection. The rate calculated 
     * directly from I32 in the main program, however, has to be corrected by the estimated ratio. This is done by
     * calling #getRatio.
     */
    bool direct_estimate;
    /** 
     * @brief Whether a collision really takes place. Use when #direct_estimate = false. 
     *
     * Whether the collision takes place or not, set in #getMomenta32_estimate_rejection. 
     * Only in use when #direct_estimate is set to false.
     */
    bool collision;
     /** 
     * @brief Ratio of real to estimated result, needed to correct the collision rate in the calling procedure. 
     *
     * As default this is set to 1, i.e. the rate calculated in the calling program directly from I32 is not changed.
     * When #direct_estimate is set to false, this rate needs to be corrected. Then #total_ratio returns (via #getRatio)
     * the ratio with which the rate needs to be corrected.
      */
    double total_ratio;
    
    /** @brief Energy of the particle selected by #getIntegral32_vegas or #getIntegral32_fast to be particle 1. */
    double E1_selected;
    /** @brief Energy of the particle selected by #getIntegral32_vegas or #getIntegral32_fast to be particle 3. */
    double E3_selected;
    /** @brief cos(gamma) for the selected pair (p1,p3). */
    double cos_gamma;
    /** @brief Indicates the selected pair (p1,p3) */      
    int N;  
};

#endif
