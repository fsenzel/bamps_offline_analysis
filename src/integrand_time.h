//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL: file:///home/fochler/svnrepos/cascade_full/main/branches/LHC_initial_conditions/src/integrand_time.h $
//$LastChangedDate: 2010-12-14 11:09:56 +0100 (Tue, 14 Dec 2010) $
//$LastChangedRevision: 241 $
//$LastChangedBy: fochler $
//---------------------------------------------
//---------------------------------------------


#ifndef INTEGRAND_TIME_H
#define INTEGRAND_TIME_H

#include "configuration.h"
#include "vegas.h"
#include "woodsaxon.h"

class integrand_time : public integrand
{
  public:
    integrand_time() {};
    integrand_time(const double b, const double t) : bImp(b),time(t) {};
    ~integrand_time() {};
  
    /** @brief Overloaded operator() that makes integrand23 a function object - for use with CUBA-vegas */
    void operator()(const int*, const double [], const int*, double []) const;
    /** @brief Overloaded operator() that makes integrand23 a function object - for use with NR-vegas */
    double operator()(const double [], double) const;
  
    void setB(const double b) { bImp = b; }
    void setTime(const double t) { time = t; }
    void setWoodSaxonParameter( const WoodSaxon& _w ) { woodSaxonParameter = _w; }
  
  private:
    double bImp;
    double time;
    WoodSaxon woodSaxonParameter;
};


#endif