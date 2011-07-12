//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL: file:///home/fochler/svnrepos/cascade_full/main/branches/restructure_branch/src/interpolation23.h $
//$LastChangedDate: 2010-07-06 16:24:24 +0200 (Tue, 06 Jul 2010) $
//$LastChangedRevision: 115 $
//$LastChangedBy: fochler $
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------

/** @file 
 * @brief Declarations for the class interpolation23.
 */


#ifndef INTERPOLATION23_H
#define INTERPOLATION23_H

/**
@author Oliver Fochler
 */

#include <vector>
#include <stdexcept>


enum I23_ERROR_CODE {I23_SUCCESS, I23_READOUT_ERROR, I23_FILE_ERROR, I23_SIZE_ERROR}; 

    
class interpolation23
{
  public:
    interpolation23();
    ~interpolation23();
    
    /** @brief find I23 via interpolation from tabulated values */
    double getI23(const double, const double, const double, const double) const;
    /** @brief returns true when interpolation can be performed with the given parameters, false otherwise */
    bool canInterpolate(const double, const double, const double, const double) const;
    
    
  private:
    /** @brief This routine reads the I23 tables from a file */
    I23_ERROR_CODE readTables();
    int get_index(const int, const int, const int, const int) const;
    bool inRange(const double, const double, const double, const double) const;
    int getShift(const std::vector<double>&) const;
    
    /** @brief get I23 using linear interpolation */
    double getI23_linearInterpolation(const double, const double, const double, const double) const;
    
    /** @brief get I23 using polynomial interpolation */
    double getI23_polyInterpolation(const double, const double, const double, const double) const;
    
    /** @brief get I23 using spline interpolation */
    double getI23_splineInterpolation(const double, const double, const double, const double) const;
    
    /** @brief polynomial interpolation and extrapolation as found in numerical recipes */
    void polynomialInterpolation(const std::vector<double>& xa, const std::vector<double>& ya, const int n, const double x, 
                                 double& y, double& dy) const;
    
    /** @brief the data structure holding the tabulated values */
    std::vector<double> I23data;        
    
    
    /** @brief number of entries in a = ln(md2) direction */
    int n_i;
    /** @brief number of entries in b = ln(lambda) direction */
    int n_j;
    /** @brief number of entries in c = ln(gamma) direction for lower values of ln(gamma) = 0 .. 3.4 */
    int n_k_low;
    /** @brief number of entries in c = ln(gamma) direction for higher values of ln(gamma) = 4.4 .. 7.4 */
    int n_k_high;
    /** @brief total number of entries in c = ln(gamma) direction */
    int n_k;   
    /** @brief number of entries in d = cos(theta) direction */
    int n_l;
    
    /** @brief total number of entries to be expected in the input file */
    int expected_entries;
    
    /** @brief spacing of tabulated a values */
    double delta_a;
    /** @brief spacing of tabulated b values */
    double delta_b;
    /** @brief spacing of tabulated c values for small c */
    double delta_c_low;
    /** @brief spacing of tabulated c values for large c*/
    double delta_c_high;
    /** @brief spacing of tabulated d values */
    double delta_d;
    
    /** @brief lowest tabulated value in a-direction */
    double a_start;
    /** @brief lowest tabulated value in b-direction */
    double b_start;
    /** @brief lowest tabulated value in c-direction */
    double c_start;
    /** @brief lowest tabulated value in d-direction */
    double d_start;
    
    /** @brief boundary between small and high c-values with a different spacing */
    double c_separator;
    
    /** @brief a 'very low' negative number */
    int minus_infinity;
};



/** @brief exception class for handling unexpected behaviour when reading the file containing the I23 tables */
class eI23_read_error : public std::runtime_error
{
  public:
    explicit eI23_read_error(const std::string& what) : std::runtime_error(what) {};

    virtual ~eI23_read_error() throw() {};
};

#endif
