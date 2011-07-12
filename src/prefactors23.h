//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------


#ifndef PREFACTORS23
#define PREFACTORS23

//#define alpha_s(s) 12.0*M_PI/((33.0-2.0*Nflavor)*log(s/lambda2))
#define alpha_s(s) 0.3

#include <math.h>
#include <iostream>
#include "particle.h"
#include "configuration.h"
#include "random.h"


/** @brief Prototype for the management of prefactors (relative to gg -> ggg) associated with 2 -> 3 interactions */
class prefactor23_generic
{
  public:
    /**
    * @brief Constructor
    */
    prefactor23_generic( const double _sF23, const double _sF32 ) : symmetryFactor_23( _sF23 ), symmetryFactor_32( _sF32 ) {};
    
    /**
    * @brief Return the associated prefactor
    * @return The prefactor associated with the process in question, relative to the gg -> ggg case
    */
    double prefactor() const {};

    
    /**
    * @brief Return the associated symmetry factor for the 2 -> 3 process
    * @return The symmetry factor associated with the process in question, relativ to the gg -> ggg case
    */
    double symFactor23() const { return ( symmetryFactor_23 ); }
    
    
    /**
    * @brief Return the associated symmetry factor for the 3 -> 2 process
    * @return The symmetry factor associated with the process in question, relativ to the ggg -> gg case
    */
    double symFactor32() const { return ( symmetryFactor_32 ); }
    
    
  protected:
    double symmetryFactor_23;
    double symmetryFactor_32;
};



/** @brief Prefactor for g + g -> g+ g + g processes, 1 per definition */
class prefactor23_gg_ggg : public prefactor23_generic
{
  public:
    prefactor23_gg_ggg() : prefactor23_generic( 1, 1 ) {};
    
    double prefactor() const
    {
//       return 0;
      return 1;
    }
    
  private:
};



/** @brief Prefactor for g + g -> q + qbar + g processes */
class prefactor23_gg_qqbarg : public prefactor23_generic
{
  public:
    prefactor23_gg_qqbarg() : prefactor23_generic( ((9 * ns_casc::Nflavor) / 8), 1 ) {};
    
    double prefactor() const
    {
      return 0 ;
    }
    
  private:
};



/** @brief Prefactor for q + qbar -> g + g + g processes */
class prefactor23_qqbar_ggg : public prefactor23_generic
{
  public:
    prefactor23_qqbar_ggg() : prefactor23_generic( (8 / (9 * ns_casc::Nflavor)), 2 ) {};
    
    double prefactor() const
    {
      return 0 ;
    }
        
  private:
};



/** @brief Prefactor for q + qbar -> q + qbar + g processes */
class prefactor23_qqbar_qqbarg : public prefactor23_generic
{
  public:
    prefactor23_qqbar_qqbarg() : prefactor23_generic( 2, 2 ) {};
    
    double prefactor() const
    {
      if ( ns_casc::Nflavor > 0 )
      {
//         return 0 ;
//         return 1 ;
        return ( 16.0 / 81.0 );
      }
      else
      {
        return 0;
      }
    }
    
  private:
};



/** @brief Prefactor for q + qbar -> q' + qbar' + g processes */
class prefactor23_qqbar_qqbarDashg : public prefactor23_generic
{
  public:
    prefactor23_qqbar_qqbarDashg() : prefactor23_generic( 2, 2 ) {};
    
    double prefactor() const
    {
      return 0 ;
    }
    
  private:
};



/** @brief Prefactor for q + g -> q + g + g (and qbar + g -> qbar + g + g) processes */
class prefactor23_qg_qgg : public prefactor23_generic
{
  public:
    prefactor23_qg_qgg() : prefactor23_generic( 2, 2 ) {};
    
    double prefactor() const
    {
      if ( ns_casc::Nflavor > 0 )
      {
//         return 0 ;
//         return 1 ;
        return ( 4.0 / 9.0 );
      }
      else
      {
        return 0;
      }
    }

  private:
};



/** @brief Prefactor for q + q -> q + q + g (qbar + qbar -> qbar + qbar + g) processes */
class prefactor23_qq_qqg : public prefactor23_generic
{
  public:
    prefactor23_qq_qqg() : prefactor23_generic( 1, 1 ) {};
    
    double prefactor() const
    {
      if ( ns_casc::Nflavor > 0 )
      {
//         return 0 ;
//         return 1 ;
        return ( 32.0 / 81.0 );
      }
      else
      {
        return 0;
      }
    }
     
  private:
};



/** @brief Prefactor for q + q' -> q + q' + g processes */
class prefactor23_qqdash_qqdashg : public prefactor23_generic
{
  public:
    prefactor23_qqdash_qqdashg() : prefactor23_generic( 2, 2 ) {};
    
    double prefactor() const
    {
      if ( ns_casc::Nflavor > 0 )
      {
//         return 0 ;
//         return 1 ;
        return ( 16.0 / 81.0 );
      }
      else
      {
        return 0;
      }
    }
    
    
  private:
};



inline void sampleFlavor23(FLAVOR_TYPE & F1, FLAVOR_TYPE & F2, FLAVOR_TYPE exclude = gluon )
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
