//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL: file:///home/fochler/svnrepos/cascade_box/main/branches/quark_branch/src/binning.h $
//$LastChangedDate: 2008-02-22 15:28:55 +0100 (Fri, 22 Feb 2008) $
//$LastChangedRevision: 21 $
//$LastChangedBy: fochler $
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------


#ifndef BINNING_H
#define BINNING_H

#include <vector>
#include <string>
using std::vector;
using std::string;

/**
@author Oliver Fochler
*/

enum binning_type {default_binning, dynamic, fixed_number, fixed_min_max_number, fixed_width};


class binning
{
  public:
    binning();
    binning(const string);
    
    binning(const double, const double, const int);
    binning(const string, const double, const double, const int);
    
    binning(const int);
    binning(const string, const int);
    
    binning(const double);
    binning(const string, const double);
    
    ~binning();
    
    void setWidth(const double);
    void setNumber(const int);
    void setFilename(const string name) {filename = name;}
    void setMinMaxN(const double x_min, const double x_max, const int n);
    void setMinMaxWidth(const double x_min, const double x_max, const double width_arg);
    
    double getWidth() const {return binWidth;}
    int getNBins() const {return nBins;}
    
    double getBinLabel(const int i);
    double getBin(const int i);
    int getBinRaw(const int i);
    
    void add(const double);
    void print();
    
  private:
    void setFilename() {setFilename("defaultBinning.dat");}
    void finish();
    
    vector<double> theList;
    vector<int> theBins;
    
    bool finished;
    
    double x_min, x_max, binWidth;
    int nBins;
    int binnedEvents;
    
    binning_type btype;
    
    string filename;    
};

#endif
