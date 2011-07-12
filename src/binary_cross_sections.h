//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------


#ifndef BINARY_XSECTIONS
#define BINARY_XSECTIONS

//#define alpha_s(s) 12.0*M_PI/((33.0-2.0*Nflavor)*log(s/lambda2))
#define alpha_s(s) 0.3

#include <math.h>
#include <iostream>
#include "particle.h"
#include "configuration.h"
#include "random.h"


/** @brief Prototype for the computation of binary cross sections and the sampling of momentum transfers */
class xsection_generic
{
  public:
    /**
    * @brief Constructor
    * @param[in] _s Mandelstam s
    * @param[in] _md2g Gluon debye mass squared
    * @param[in] _md2q Quark debye mass squared
    */
    xsection_generic(const double _s, const double _md2g, const double _md2q)
      : s(_s), md2g(_md2g), md2q(_md2q) {};
    
    /**
    * @brief Compute the total cross section
    * @return Total cross section in 1/GeV^2
    */
    double totalCrossSection() const {};

    /**
    * @brief Compute the total cross section without prefactors - needed for sampling
    * @return Total cross section without prefactors
    */
    double totalCrossSection_noPreFactors() const {};
    
    /**
    * @brief Inverse of the integral over the differential cross section without prefactors - needed for sampling
    *
    * Let the integral over the differential cross section yield Integrate[f(x'),{x',0,x}] = A(x) then this function
    * returns x(A) (without any prefactors as in totalCrossSection_noPreFactors() ). A uniformly sampled A (in [0, A_max])
    * then immediately yields x sampled according to f(x). Here A_max corresponds to the total cross section (without prefactors).
    *
    * @return Inverse integral over the differential cross section
    */
    double inverseIntegral_noPreFactors( const double _A ) const {};
    
    /** @brief Mandelstam s */
    double s;
    
  protected:
    /** @brief Gluon debye mass squared */
    double md2g;
    /** @brief Quark debye mass squared */
    double md2q;
};



/** @brief Binary cross section and sampling of momentum transfers for gg -> gg processes */
class xsection_gg_gg : public xsection_generic
{
  public:
    xsection_gg_gg(const double _s, const double _md2g, const double _md2q)
    : xsection_generic( _s, _md2g, _md2q ) {};
    
    double totalCrossSection() const
    {
      return ( 9 * M_PI * pow( alpha_s(s), 2 ) * s / ( 2 * md2g * ( 4*md2g + s) ) );
//       return 0;
    }
    
    double totalCrossSection_noPreFactors() const  { return ( s / ( 4*md2g + s) ); }
    double inverseIntegral_noPreFactors( const double _A ) const  { return ( _A * md2g / (1 - _A) ); }
    
  private:
};



/** @brief Binary cross section and sampling of momentum transfers for qqbar -> gg processes */
class xsection_qqbar_gg : public xsection_generic
{
  public:
    xsection_qqbar_gg(const double _s, const double _md2g, const double _md2q)
    : xsection_generic( _s, _md2g, _md2q ) {};
    
    double totalCrossSection() const
    {
      return ( 32 * M_PI * pow( alpha_s(s), 2 ) / ( 27 * s ) * log( 1 + s / (4*md2q) ) ) ;
//       return 0;
    }
    
    double totalCrossSection_noPreFactors() const  { return log( 1 + s / (4*md2q) ) ; }
    double inverseIntegral_noPreFactors( const double _A ) const  { return ( md2q * (exp(_A) - 1) ); }
    
  private:
};



/** @brief Binary cross section and sampling of momentum transfers for qqbar -> qqbar processes */
class xsection_qqbar_qqbar : public xsection_generic
{
  public:
    xsection_qqbar_qqbar(const double _s, const double _md2g, const double _md2q)
    : xsection_generic( _s, _md2g, _md2q ) {};
    
    double totalCrossSection() const
    {
      return ( 8 * M_PI * pow( alpha_s(s), 2 ) * s / ( 9 * md2g * (4 * md2g + s) ) ) ;
//       return 0;
    }
    
    double totalCrossSection_noPreFactors() const  { return  ( s / (4 * md2g + s) ) ; }
    double inverseIntegral_noPreFactors( const double _A ) const  { return ( _A * md2g / (1 - _A) ); }
    
  private:
};



/** @brief Binary cross section and sampling of momentum transfers for qqbar -> q'qbar' processes */
class xsection_qqbar_qqbarDash : public xsection_generic
{
  public:
    xsection_qqbar_qqbarDash(const double _s, const double _md2g, const double _md2q)
    : xsection_generic( _s, _md2g, _md2q ) {};
    
    double totalCrossSection() const
    {
      return ( (ns_casc::Nflavor - 1) * 8 * M_PI * pow( alpha_s(s), 2 ) / ( 27 * s ) ) ;
//       return 0;
    }
    
    double totalCrossSection_noPreFactors() const  { return (2 * s / 3); }
    double inverseIntegral_noPreFactors( const double _A ) const
    {
      return (-( s/2 + pow(s,2)/(2.*pow(-6*_A*pow(s,2) + 2*pow(s,3) + sqrt(36*pow(_A,2)*pow(s,4) - 24*_A*pow(s,5) + 5*pow(s,6)), 0.3333333))
      - pow(-6*_A*pow(s,2) + 2*pow(s,3) + sqrt(36*pow(_A,2)*pow(s,4) - 24*_A*pow(s,5) + 5*pow(s,6)), 0.33333333)/2. ));
    }
    
  private:
};



/** @brief Binary cross section and sampling of momentum transfers for qg -> qg processes */
class xsection_qg_qg : public xsection_generic
{
  public:
    xsection_qg_qg(const double _s, const double _md2g, const double _md2q)
    : xsection_generic( _s, _md2g, _md2q ) {};
    
    double totalCrossSection() const
    {
      return ( 2 * M_PI * pow( alpha_s(s), 2 ) * s / ( md2g * ( 4*md2g + s) )  ) ;
//       return 0;
    }
    
    double totalCrossSection_noPreFactors() const  { return ( s /  (4*md2g + s) ); }
    double inverseIntegral_noPreFactors( const double _A ) const  { return ( _A * md2g / (1 - _A)  ); }
    
  private:
};



/** @brief Binary cross section and sampling of momentum transfers for gg -> qqbar processes */
class xsection_gg_qqbar : public xsection_generic
{
  public:
    xsection_gg_qqbar(const double _s, const double _md2g, const double _md2q)
    : xsection_generic( _s, _md2g, _md2q ) {};
    
    double totalCrossSection() const
    {
      return ( ns_casc::Nflavor * M_PI * pow( alpha_s(s), 2 ) / ( 3 * s ) * log( 1 + s / (4*md2q) ) ) ;
//       return 0;
    }
    
    double totalCrossSection_noPreFactors() const  { return ( log( 1 + s / (4*md2q) ) ); }
    double inverseIntegral_noPreFactors( const double _A ) const  { return ( md2q * (exp(_A) - 1) ); }
    
  private:
};



/** @brief Binary cross section and sampling of momentum transfers for qq -> qq (qbarqbar -> qbarqbar) processes */
class xsection_qq_qq : public xsection_generic
{
  public:
    xsection_qq_qq(const double _s, const double _md2g, const double _md2q)
    : xsection_generic( _s, _md2g, _md2q ) {};
    
    double totalCrossSection() const
    {
      return ( 8 * M_PI * pow( alpha_s(s), 2 ) * s / ( 9 * md2g * ( 4*md2g + s) ) ) ;
//       return 0;
    }
    
    double totalCrossSection_noPreFactors() const  { return ( s /  (4*md2g + s) ); }
    double inverseIntegral_noPreFactors( const double _A ) const  { return ( _A * md2g / (1 - _A) ); }
    
  private:
};



/** @brief Binary cross section and sampling of momentum transfers for qq' -> qq' processes */
class xsection_qqdash_qqdash : public xsection_generic
{
  public:
    xsection_qqdash_qqdash(const double _s, const double _md2g, const double _md2q)
    : xsection_generic( _s, _md2g, _md2q ) {};
    
    double totalCrossSection() const
    {
      return ( 8 * M_PI * pow( alpha_s(s), 2 ) * s / ( 9 * md2g * (4 * md2g + s) ) ) ;
//       return 0;
    }
    
    double totalCrossSection_noPreFactors() const  { return  ( s / (4 * md2g + s) ) ; }
    double inverseIntegral_noPreFactors( const double _A ) const  { return ( _A * md2g / (1 - _A) ); }
    
  private:
};




/** @brief Default type trait for sampling in mandelstam t (instead of qt^2) 
*
* This template serves the purpose to determine whether a class derived from xsection_generic
* actually samples in mandelstam t or directly in the transverse momentum transfer qt^2 (the default).
* In order to change the default value for a certain process type, a specialized version of this template
* needs to be created.
*/
template<typename T>
struct returns_mandelstam_t { static const bool value = false; };

/** @brief Specialized type trait for sampling in mandelstam t for xsection_qqbar_qqbarDash */ 
template<>
struct returns_mandelstam_t<xsection_qqbar_qqbarDash> { static const bool value = true; };


/** @brief Sample the transverse momentum transfer qt^2*/
template<class T>
inline double sampleBinaryPT2( T& _xsectionObj )
{
  // if the cross section object for this process returns mandelstam t instead of qt^2, we need to convert it
  if ( returns_mandelstam_t<T>::value )
  {
    double t = _xsectionObj.inverseIntegral_noPreFactors( ran2() * _xsectionObj.totalCrossSection_noPreFactors() );
    return ( -( pow(t,2)/ _xsectionObj.s + t ) );
  }
  else
  {
    return _xsectionObj.inverseIntegral_noPreFactors( ran2() * _xsectionObj.totalCrossSection_noPreFactors() );
  }
}


inline void sampleFlavor(FLAVOR_TYPE & F1, FLAVOR_TYPE & F2, FLAVOR_TYPE exclude = gluon )
{
  int flav;
  do {
    flav = static_cast<int>( ns_casc::Nflavor * ran2() ) + 1;
    if ( flav > ns_casc::Nflavor ) 
      flav = ns_casc::Nflavor;
    F2 = static_cast<FLAVOR_TYPE>( 2 * flav );
    F1 = static_cast<FLAVOR_TYPE>( 2 * flav - 1 );
  } while ( F1 == exclude || F2 == exclude );
}


#endif
