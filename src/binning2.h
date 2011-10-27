//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------
//---------------------------------------------


#ifndef BINNING2_H
#define BINNING2_H

#include "binning.h"

#include <vector>
#include <string>
using std::vector;
using std::string;



// enum binning_type {default_binning, dynamic, fixed_number, fixed_min_max_number, fixed_width};


class binning2d
{
  public:
//     binning();
//     binning(const string);
//     
//     binning(const double, const double, const int);
    binning2d(const string name, const double min_x_arg, const double max_x_arg, const int n_x, const double min_y_arg, const double max_y_arg, const int n_y);
    
//     binning(const int);
//     binning(const string, const int);
//     
//     binning(const double);
//     binning(const string, const double);
//     
//     ~binning();
    
//     void setWidth(const double);
//     void setNumber(const int);
    void setFilename(const string name) {filename = name;}
//     void setMinMaxN(const double x_min, const double x_max, const int n);
//     void setMinMaxWidth(const double x_min, const double x_max, const double width_arg);
    
//     double getWidth() const {return binWidth;}
//     int getNBins() const {return nBins;}
    
//     double getBinLabel(const int i);
//     double getBin(const int i);
//     int getBinRaw(const int i);
    
    void add(const double x, const double y);
    void print();
    
  private:
    void setFilename() {setFilename("defaultBinning.dat");}
//     void finish();
    
//     vector<double> theList;
//     vector<int> theBins;
    vector<vector<int> > theBins;
    
    bool finished;
    
    double x_min, x_max, y_min, y_max, binWidth_x, binWidth_y;
    int nBins_x, nBins_y;
    int binnedEvents;
    
    binning_type btype;
    
    string filename;    
};


class binningValues
{
  public:
//     binning();
//     binning(const string);
//     
//     binning(const double, const double, const int);
    binningValues(const string name, const double min_arg, const double max_arg, const int n);
    
//     binning(const int);
//     binning(const string, const int);
//     
//     binning(const double);
//     binning(const string, const double);
//     
//     ~binning();
    
//     void setWidth(const double);
//     void setNumber(const int);
    void setFilename(const string name) {filename = name;}
//     void setMinMaxN(const double x_min, const double x_max, const int n);
//     void setMinMaxWidth(const double x_min, const double x_max, const double width_arg);
    
//     double getWidth() const {return binWidth;}
//     int getNBins() const {return nBins;}
    
//     double getBinLabel(const int i);
//     double getBin(const int i);
//     int getBinRaw(const int i);
    
    void add(const double x, const double v);
    void print();
    
  private:
    void setFilename() {setFilename("defaultBinning.dat");}
//     void finish();
    
//     vector<double> theList;
//     vector<int> theBins;
    vector<double > theBins;
    vector<int > NumberInBins;
    
    bool finished;
    
    double x_min, x_max, binWidth;
    int nBins;
    int binnedEvents;
    
    binning_type btype;
    
    string filename;    
};


class binningValues2d
{
  public:
//     binning();
//     binning(const string);
//     
//     binning(const double, const double, const int);
    binningValues2d(const string name, const double min_x_arg, const double max_x_arg, const int n_x, const double min_y_arg, const double max_y_arg, const int n_y);
    
//     binning(const int);
//     binning(const string, const int);
//     
//     binning(const double);
//     binning(const string, const double);
//     
//     ~binning();
    
//     void setWidth(const double);
//     void setNumber(const int);
    void setFilename(const string name) {filename = name;}
//     void setMinMaxN(const double x_min, const double x_max, const int n);
//     void setMinMaxWidth(const double x_min, const double x_max, const double width_arg);
    
//     double getWidth() const {return binWidth;}
//     int getNBins() const {return nBins;}
    
//     double getBinLabel(const int i);
//     double getBin(const int i);
//     int getBinRaw(const int i);
    
    void add(const double x, const double y, const double v);
    void print();
    
  private:
    void setFilename() {setFilename("defaultBinning.dat");}
    void finish();
    
//     vector<double> theList;
//     vector<int> theBins;
    vector<vector<double> > theBins;
    vector<vector<int> > NumberInBins;
    
    bool finished;
    
    double x_min, x_max, y_min, y_max, binWidth_x, binWidth_y;
    int nBins_x, nBins_y;
    int binnedEvents;
    
    binning_type btype;
    
    string filename;    
};

#endif
