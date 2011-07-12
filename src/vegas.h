//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL: file:///home/fochler/svnrepos/cascade_full/main/branches/restructure_branch/src/vegas.h $
//$LastChangedDate: 2009-07-20 11:50:59 +0200 (Mon, 20 Jul 2009) $
//$LastChangedRevision: 77 $
//$LastChangedBy: fochler $
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------

#ifndef VEGAS_H
#define VEGAS_H


// integrand is a prototype for a function object that can be passed to the vegas routine.
// The operator () is overloaded as a pure virtual function, i.e. an object of integrand cannot be created directly.
// Rather one must derive a class from integrand that implements operator()(...), create an object of this class and
// pass it to vegas.
// For example
//                class integrand23 : public integrand
//                {  
//                 public:
//                   double operator()(double[], double) const;
//     
//                  <...more stuff...>
//                };
// The overloaded operator() allows for calls such as theIntegrand(a,b) where theIntegrand is of type integrand23. Just as if theIntegrand
// was a normal function.

class integrand
{
  public:
    virtual void operator()(const int *ndim, const double xx[], const int *ncomp, double ff[]) const = 0;
    virtual double operator()(const double [], double) const = 0;   // pure virtual (denoted by =0), must be overridden in derived class
    virtual ~integrand() {};
};


void vegas(int ndim, integrand& fxn, double *tgral, double *sd, double *chi2a);  // takes a function object derived from integrand as argument
void rebin(double rc, int nd, double r[], double xin[], double xi[]);

#endif
