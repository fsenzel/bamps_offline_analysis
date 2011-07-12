#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <vector>
using std::vector;


/**
	@author Oliver Fochler <fochler@th.physik.uni-frankfurt.de>
 */
class interpolation
{
  public:
    interpolation();
    interpolation(const vector<double>&, const vector<double>&);
    interpolation(const vector<double>&, const vector<double>&, const double, const double);
    ~interpolation();
    
    void spline_construction( const vector<double>&, const vector<double>&, const double S_initial = 1e20, const double S_final = 1e20);
    double spline_interpolation(const vector<double>&, const vector<double>&, const double, const int index = -1) const;
    
    
  private:
    vector<double> S;  // the second derivatives

    int bisection(const vector<double>&, const double) const;
    
};

#endif
