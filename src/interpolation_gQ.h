//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL: file:///home/fochler/svnrepos/cascade_box/main/branches/quark_branch/src/interpolation_gc.h $
//$LastChangedDate: 2009-06-23 12:06:07 +0200 (Tue, 23 Jun 2009) $
//$LastChangedRevision: 90 $
//$LastChangedBy: fochler $
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------

/** @file 
 * @brief Declarations for the class interpolation_gQ.
 */


#ifndef INTERPOLATIONGC_H
#define INTERPOLATIONGC_H

/**
@author Oliver Fochler
 */

#include <vector>
#include <stdexcept>
#include "configuration.h"


enum IgQ_ERROR_CODE {IgQ_SUCCESS, IgQ_READOUT_ERROR, IgQ_FILE_ERROR, IgQ_SIZE_ERROR}; 

    
class interpolation_gQ
{
  public:
    interpolation_gQ();
    ~interpolation_gQ();
    
    void configure(const partclType, const bool);
    
    /** @brief find IgQ via interpolation from tabulated values */
    double getIgQ(const double, const double) const;
    
    
  private:
    /** @brief This routine reads the IgQ tables from a file */
    IgQ_ERROR_CODE readTables();
    int get_index(const int, const int) const;
//     bool inRange(const double, const double, const double, const double) const;
//     int getShift(const std::vector<double>&) const;
    
    /** @brief get IgQ using linear interpolation */
    double getIgQ_linearInterpolation(const double, const double) const;
    

    
    /** @brief the data structure holding the tabulated values */
    std::vector<double> IgQdata;        
    
    
    /** @brief number of entries in a = ln(md2) direction */
    int n_i;
    /** @brief number of entries in b = ln(s) direction */
    int n_j;
    
    /** @brief total number of entries to be expected in the input file */
    int expected_entries;
    
    /** @brief spacing of tabulated a values */
    double delta_a;
    /** @brief spacing of tabulated b values */
    double delta_b;

    
    /** @brief lowest tabulated value in a-direction */
    double a_start;
    /** @brief lowest tabulated value in b-direction */
    double b_start;
    
    double Mass;
    bool asRun;
    partclType species;
};



/** @brief exception class for handling unexpected behaviour when reading the file containing the IgQ tables */
class eIgQ_read_error : public std::runtime_error
{
  public:
    explicit eIgQ_read_error(const std::string& what) : std::runtime_error(what) {};

    virtual ~eIgQ_read_error() throw() {};
};

#endif
