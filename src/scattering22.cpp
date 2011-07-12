//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL: file:///home/fochler/svnrepos/cascade_full/main/branches/restructure_branch/src/scattering22.cpp $
//$LastChangedDate: 2010-07-23 17:05:12 +0200 (Fri, 23 Jul 2010) $
//$LastChangedRevision: 153 $
//$LastChangedBy: fochler $
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------

#include <iostream>
#include <math.h>
#include <algorithm>

#include "scattering22.h"
#include "lorentz.h"
#include "configuration.h"
#include "random.h"
#include "binary_cross_sections.h"

#define alpha_s(s) 0.3

using std::cout;
using std::endl;


scattering22::scattering22()
    : P1( NULL ), P2( NULL ), P1cm( NULL ), P2cm( NULL ), F1( gluon ), F2( gluon ), md2_gluon( 0 ), md2_quark( 0 ), s( 0 )
{
}


/**
 * Call this constructor when creating a scattering22 object with already know properties of a particle pair.
 * Alternatively the parameters can be set by calling #scattering22::setParameter.
 *
 * @param[in] P1arg[] 4-momentum vector of ingoing particle 1
 * @param[in] P2arg[] 4-momentum vector of ingoing particle 2
 * @param[in] F1_arg flavor of ingoing particle 1
 * @param[in] F2_arg flavor of ingoing particle 1
 * @param[in] s_arg mandelstam s
 * @param[in] md2g_arg gluon debye mass squared divided by alpha_s
 * @param[in] md2q_arg quark debye mass squared divided by alpha_s
 */
scattering22::scattering22( const double P1arg[], const double P2arg[], const FLAVOR_TYPE F1_arg, const FLAVOR_TYPE F2_arg,
                            const double s_arg, const double md2g_arg, const double md2q_arg )
    : F1( F1_arg ), F2( F2_arg ), s( s_arg ), md2_gluon( md2g_arg ), md2_quark( md2q_arg )
{
  P1 = new double[4];
  P2 = new double[4];
  P1cm = new double[4];
  P2cm = new double[4];

  for ( int i = 0; i <= 3; i++ )
  {
    P1[i] = P1arg[i];
    P2[i] = P2arg[i];
  }

  double totE = P1[0] + P2[0];
  for ( int i = 1; i <= 3; i++ )
    beta_vec[i] = ( P1[i] + P2[i] ) / totE;
  //beta = sqrt( pow(beta_vec[1],2.0) + pow(beta_vec[2],2.0) + pow(beta_vec[3],2.0) );

  lorentz( beta_vec, P1, P1cm );
  lorentz( beta_vec, P2, P2cm );
}



scattering22::~scattering22()
{
  delete[] P1;
  P1 = NULL;
  delete[] P2;
  P2 = NULL;
  delete[] P1cm;
  P1cm = NULL;
  delete[] P2cm;
  P2cm = NULL;
}


/**
 * Used to set the parameters. With this an existing scattering22 object can be re-used for a new particle pair.
 *
 * @param[in] P1arg[] 4-momentum vector of ingoing particle 1
 * @param[in] P2arg[] 4-momentum vector of ingoing particle 2
 * @param[in] F1_arg flavor of ingoing particle 1
 * @param[in] F2_arg flavor of ingoing particle 1
 * @param[in] s_arg mandelstam s
 * @param[in] md2g_arg gluon debye mass squared divided by alpha_s
 * @param[in] md2q_arg quark debye mass squared divided by alpha_s
 */
void scattering22::setParameter( const double P1arg[], const double P2arg[], const FLAVOR_TYPE F1_arg, const FLAVOR_TYPE F2_arg,
                                 const double s_arg, const double md2g_arg, const double md2q_arg )
{
  F1 = F1_arg;
  F2 = F2_arg;
  s = s_arg;
  md2_gluon = md2g_arg;
  md2_quark = md2q_arg;

  delete[] P1;
  P1 = NULL;
  delete[] P2;
  P2 = NULL;
  delete[] P1cm;
  P1cm = NULL;
  delete[] P2cm;
  P2cm = NULL;

  P1 = new double[4];
  P2 = new double[4];
  P1cm = new double[4];
  P2cm = new double[4];

  for ( int i = 0; i <= 3; i++ )
  {
    P1[i] = P1arg[i];
    P2[i] = P2arg[i];
  }

  double totE = P1[0] + P2[0];
  for ( int i = 1; i <= 3; i++ )
    beta_vec[i] = ( P1[i] + P2[i] ) / totE;
  //beta = sqrt( pow(beta_vec[1],2.0) + pow(beta_vec[2],2.0) + pow(beta_vec[3],2.0) );

  lorentz( beta_vec, P1, P1cm );
  lorentz( beta_vec, P2, P2cm );
}


/**
 * Get the cross section for 2->2 processes
 *
 * @param[out] initialStateIndex Integer flag set according to the initial state that has been processed.
 * @return cross section in units 1/GeV^2
 */
double scattering22::getXSection22( int& initialStateIndex ) const
{
  if ( s < ( 1.1 * ns_casc::lambda2 ) )
  {
    return 0;
  }

  // sort F1 and F2 such that comparisons below are easier
  unsigned int _F1 = std::min( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );
  unsigned int _F2 = std::max( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );
  double cs = 0;

  if (( _F1 + _F2 ) == 0 ) // gg -> gg, gg -> qqbar
  {
    xsection_gg_gg csObj1( s, md2_gluon, md2_quark );
    xsection_gg_qqbar csObj2( s, md2_gluon, md2_quark );

    initialStateIndex = 0;

    cs = csObj1.totalCrossSection() + csObj2.totalCrossSection();
  }
  else if ( _F1 == _F2 )  // qq -> qq, qbarqbar -> qbarqbar
  {
    xsection_qq_qq csObj1( s, md2_gluon, md2_quark );

    initialStateIndex = 1;

    cs = csObj1.totalCrossSection();
  }
  else if (( _F1 * _F2 ) == 0 ) // gq -> gq, gqbar -> gqbar
  {
    xsection_qg_qg csObj1( s, md2_gluon, md2_quark );

    initialStateIndex = 2;

    cs = csObj1.totalCrossSection();
  }
  else if (( _F2 - _F1 ) == 1 && ( _F2 % 2 ) == 0 )  // qqbar -> qqbar, qqbar -> q'qbar', qqbar -> gg
  {
    xsection_qqbar_qqbar csObj1( s, md2_gluon, md2_quark );
    xsection_qqbar_qqbarDash csObj2( s, md2_gluon, md2_quark );
    xsection_qqbar_gg csObj3( s, md2_gluon, md2_quark );

    initialStateIndex = 3;

    cs = csObj1.totalCrossSection() + csObj2.totalCrossSection() + csObj3.totalCrossSection();
  }
  else // qq' -> qq', qqbar' -> qqbar'
  {
    xsection_qqdash_qqdash csObj1( s, md2_gluon, md2_quark );

    initialStateIndex = 4;

    cs = csObj1.totalCrossSection();
  }

  return cs;
}



/**
 * Get the cross section for ELASTIC 2->2 processes
 *
 * @param[out] initialStateIndex Integer flag set according to the initial state that has been processed.
 * @return cross section in units 1/GeV^2
 */
double scattering22::getXSectionElastic( int& initialStateIndex ) const
{
  if ( s < ( 1.1 * ns_casc::lambda2 ) )
  {
    return 0;
  }

  // sort F1 and F2 such that comparisons below are easier
  unsigned int _F1 = std::min( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );
  unsigned int _F2 = std::max( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );
  double cs = 0;

  if (( _F1 + _F2 ) == 0 ) // gg -> gg
  {
    xsection_gg_gg csObj1( s, md2_gluon, md2_quark );
    initialStateIndex = 0;
    cs = csObj1.totalCrossSection();
  }
  else if ( _F1 == _F2 )  // qq -> qq, qbarqbar -> qbarqbar
  {
    xsection_qq_qq csObj1( s, md2_gluon, md2_quark );
    initialStateIndex = 1;
    cs = csObj1.totalCrossSection();
  }
  else if (( _F1 * _F2 ) == 0 ) // gq -> gq, gqbar -> gqbar
  {
    xsection_qg_qg csObj1( s, md2_gluon, md2_quark );
    initialStateIndex = 2;
    cs = csObj1.totalCrossSection();
  }
  else if (( _F2 - _F1 ) == 1 && ( _F2 % 2 ) == 0 )  // qqbar -> qqbar, qqbar -> q'qbar', qqbar -> gg
  {
    xsection_qqbar_qqbar csObj1( s, md2_gluon, md2_quark );
    initialStateIndex = 3;
    cs = csObj1.totalCrossSection();
  }
  else // qq' -> qq', qqbar' -> qqbar'
  {
    xsection_qqdash_qqdash csObj1( s, md2_gluon, md2_quark );
    initialStateIndex = 4;
    cs = csObj1.totalCrossSection();
  }

  return cs;
}












/**
 * This routine samples the transverse momentum transfer of a given collision according to the differential
 * cross section.
 *
 * @param[out] PT2 transverse momentum transfer squared
 * @param[out] typ the type of collision, used for analysis purposes (e.g. 221 for gg->gg)
 * @param[out] F1arg Sampled flavor of outgoing particle 1
 * @param[out] F2arg Sampled flavor of outgoing particle 2
 */
void scattering22::getMomenta22( double& PT2, int& typ, FLAVOR_TYPE& F1arg, FLAVOR_TYPE& F2arg ) const
{
  F1arg = F1;
  F2arg = F2;
  // sort F1 and F2 such that comparisons below are easier
  // convert to unsigned int to allow for some math that speeds up the decisions
  unsigned int _F1 = std::min( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );
  unsigned int _F2 = std::max( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );

  if (( _F1 + _F2 ) == 0 ) // gg -> gg, gg -> qqbar
  {
    xsection_gg_gg csObj1( s, md2_gluon, md2_quark );
    xsection_gg_qqbar csObj2( s, md2_gluon, md2_quark );
    double cs1 = csObj1.totalCrossSection();
    double cs2 = csObj2.totalCrossSection();

    if ( ran2() < cs1 / ( cs1 + cs2 ) )
    {
      typ = 221; // gg -> gg
      PT2 = sampleBinaryPT2( csObj1 );
    }
    else
    {
      typ = 222; // gg -> qqbar
      PT2 = sampleBinaryPT2( csObj2 );
      sampleFlavor( F1arg, F2arg );
    }
  }
  else if ( _F1 == _F2 )  // qq -> qq, qbarqbar -> qbarqbar
  {
    xsection_qq_qq csObj1( s, md2_gluon, md2_quark );
    typ = 227;
    PT2 = sampleBinaryPT2( csObj1 );
  }
  else if (( _F1 * _F2 ) == 0 ) // gq -> gq, gqbar -> gqbar
  {
    xsection_qg_qg csObj1( s, md2_gluon, md2_quark );
    typ = 223;
    PT2 = sampleBinaryPT2( csObj1 );
  }
  else if (( _F2 - _F1 ) == 1 && ( _F2 % 2 ) == 0 )  // qqbar -> qqbar, qqbar -> q'qbar', qqbar -> gg
  {
    xsection_qqbar_qqbar csObj1( s, md2_gluon, md2_quark );
    xsection_qqbar_qqbarDash csObj2( s, md2_gluon, md2_quark );
    xsection_qqbar_gg csObj3( s, md2_gluon, md2_quark );
    double cs1 = csObj1.totalCrossSection();
    double cs2 = csObj2.totalCrossSection();
    double cs3 = csObj3.totalCrossSection();

    double select = ran2() * ( cs1 + cs2 + cs3 );

    if ( select < cs1 )  // qqbar -> qqbar
    {
      typ = 224;
      PT2 = sampleBinaryPT2( csObj1 );
    }
    else if ( select < ( cs1 + cs2 ) )  // qqbar -> q'qbar'
    {
      typ = 225;
      PT2 = sampleBinaryPT2( csObj2 );
      sampleFlavor( F1arg, F2arg, F1 );  // sample flavor excluding F1 (or anti-F1 = F2)
    }
    else  // qqbar -> gg
    {
      typ = 226;
      PT2 = sampleBinaryPT2( csObj3 );
      F1arg = gluon;
      F2arg = gluon;
    }
  }
  else // qq' -> qq', qqbar' -> qqbar'
  {
    typ = 228; // qq' -> qq', qqbar' -> qqbar'

    xsection_qqdash_qqdash csObj1( s, md2_gluon, md2_quark );
    PT2 = sampleBinaryPT2( csObj1 );
  }
}




/**
 * This routine samples the transverse momentum transfer of a given ELASTIC collision according to the differential
 * cross section.
 *
 * @param[out] PT2 transverse momentum transfer squared
 * @param[out] typ the type of collision, used for analysis purposes (e.g. 221 for gg->gg)
 */
void scattering22::getMomentaElastic( double& PT2, int& typ ) const
{
  // sort F1 and F2 such that comparisons below are easier
  // convert to unsigned int to allow for some math that speeds up the decisions
  unsigned int _F1 = std::min( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );
  unsigned int _F2 = std::max( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );

  if (( _F1 + _F2 ) == 0 ) // gg -> gg
  {
    xsection_gg_gg csObj1( s, md2_gluon, md2_quark );
    typ = 221; // gg -> gg
    PT2 = sampleBinaryPT2( csObj1 );
  }
  else if ( _F1 == _F2 )  // qq -> qq, qbarqbar -> qbarqbar
  {
    xsection_qq_qq csObj1( s, md2_gluon, md2_quark );
    typ = 227;
    PT2 = sampleBinaryPT2( csObj1 );
  }
  else if (( _F1 * _F2 ) == 0 ) // gq -> gq, gqbar -> gqbar
  {
    xsection_qg_qg csObj1( s, md2_gluon, md2_quark );
    typ = 223;
    PT2 = sampleBinaryPT2( csObj1 );
  }
  else if (( _F2 - _F1 ) == 1 && ( _F2 % 2 ) == 0 )  // qqbar -> qqbar
  {
    xsection_qqbar_qqbar csObj1( s, md2_gluon, md2_quark );
    typ = 224;
    PT2 = sampleBinaryPT2( csObj1 );
  }
  else // qq' -> qq', qqbar' -> qqbar'
  {
    typ = 228; // qq' -> qq', qqbar' -> qqbar'
    xsection_qqdash_qqdash csObj1( s, md2_gluon, md2_quark );
    PT2 = sampleBinaryPT2( csObj1 );
  }
}









void scattering22::getMomenta22_isotropic( double& PT2, int& typ ) const
{
  double theta_min = 0;
  double theta_max = M_PI / 2;
  double theta;

  if (( F1 + F2 ) == 0 ) // gluon,gluon -> gluon,gluon
  {
    typ = 221;// gg->gg
    theta = theta_min + ran2() * ( theta_max - theta_min );
    PT2 = s / 4 * pow( sin( theta ), 2.0 );
  }
  else
    cout << "Purley gluonic plasma expected when using the isotropic 2->2 routines." << endl;
}



/**
 * Samples new momentum vectors for the outgoing particles, using a given PT2.
 *
 * @param[out] P1[] Momentum vector of outgoing particle 1
 * @param[out] P2[] Momentum vector of outgoing particle 2
 * @param[in] R1[] Space-time vector of ingoing particle 1
 * @param[in] R2[] Space-time vector of ingoing particle 2
 * @param[in] PT2 Transverse momentum transfer squared, obtained via #scattering22::getMomenta22
 */
void scattering22::setNewMomenta22( double P1[4], double P2[4], const double R1[4], const double R2[4], const double PT2 )
{
  double PT, PZ, c;
  double R1cm[4], R2cm[4];
  double PP[4], TT[4];

  lorentz( beta_vec, R1, R1cm );
  lorentz( beta_vec, R2, R2cm );

  rotation( P1cm, R1cm, R2cm, PP, TT );

  PT = sqrt( PT2 );
  PZ = sqrt( s / 4 - PT2 );

  c = 0.0;
  for ( int i = 1; i <= 3; i++ )
  {
    P1cm[i] = PT * TT[i] + PZ * PP[i];
    P2cm[i] = -P1cm[i];
    c += P1cm[i] * P1cm[i];
  }
  P1cm[0] = P2cm[0] = sqrt( c );

  double inv_beta[4];
  for ( int j = 1;j <= 3;j++ )
  {
    inv_beta[j] = -beta_vec[j];
  }
    
  lorentz( inv_beta, P1cm, P1 );
  lorentz( inv_beta, P2cm, P2 );
}


/**
 * Auxiliary routine. Determines a direction for the transverse momtentum vector in the azimuthal plane.
 * Used by #scattering22::setNewMomenta22
 *
 * @param[in] P[] Momentum vector of ingoing particle 1
 * @param[in] R1[] Space-time vector of ingoing particle 1
 * @param[in] R2[] Space-time vector of ingoing particle 2
 * @param[out] TT[] Direction in the azimuthal plane between the two particles.
 */
void scattering22::rotation( const double P[4], const double R1[4], const double R2[4], double PP[4], double TT[4] ) const
{

  double min = 1.0, phi, sinus, cosinus, c1, c2, c3;
  double transv[4];
  int mini = -1;


  c1 = c2 = 0.0;
  for ( int i = 1;i <= 3;i++ )
  {
    c1 += P[i] * P[i];
    c2 += P[i] * ( R2[i] - R1[i] );
  }

  c3 = 0.0;
  for ( int j = 1;j <= 3;j++ )
  {
    PP[j] = P[j] / sqrt( c1 );
    TT[j] = c1 * ( R1[j] - R2[j] ) + c2 * P[j];
    c3 += TT[j] * TT[j];
  }
  c3 = sqrt( c3 );

  // r2-r1 // p1 => TT[] = 0
  if ( c3 < 1.0e-8 )
  {
    for ( int i = 1;i <= 3;i++ )
    {
      if ( fabs( PP[i] ) < min )
      {
        min = fabs( PP[i] );
        mini = i;
      }
    }
    switch ( mini )
    {
    case 1:
      TT[1] = 0.0;
      TT[2] = PP[3];
      TT[3] = -PP[2];
      break;
    case 2:
      TT[1] = -PP[3];
      TT[2] = 0.0;
      TT[3] = PP[1];
      break;
    case 3:
      TT[1] = PP[2];
      TT[2] = -PP[1];
      TT[3] = 0.0;
      break;
    default:
      cout << "Error in rotation()" << endl;
    }
  }
  //----------------------

  phi = 2.0 * M_PI * ran2();
  sinus = sin( phi );
  cosinus = cos( phi );
  transv[1] = TT[1] * cosinus + ( PP[2] * TT[3] - PP[3] * TT[2] ) * sinus;
  transv[2] = TT[2] * cosinus + ( PP[3] * TT[1] - PP[1] * TT[3] ) * sinus;
  transv[3] = TT[3] * cosinus + ( PP[1] * TT[2] - PP[2] * TT[1] ) * sinus;

  c3 = 0.0;
  for ( int j = 1;j <= 3;j++ )
    c3 += transv[j] * transv[j];
  c3 = sqrt( c3 );

  for ( int i = 1;i <= 3;i++ )
    TT[i] = transv[i] / c3;
}


/**
 * Computes the transport cross section in the case of gluon gluon scattering.
 *
 * @return transport cross section in units 1/GeV^2
 */
double scattering22::getTransportXSection22() const
{
  double as = 0.3;
  double prefactor = 36 * M_PI * pow( as, 2 ) / 2;
  double A = log(( s + 4 * md2_gluon ) / ( 4 * md2_gluon ) );
  double B = s / ( s + 4 * md2_gluon );

  return prefactor / s * ( A - B );
}


/*
// older/alternative implementations for the 2->2 cross section
// not adapted to the scattering22 structure yet
//
double crosection(int F1, int F2, const double M1, const double M2, const  double ss) //mb
{
  double factor,as,M,cs;
  int tmp;

  factor = 1/2.5682;// dimension transform from 1/GeV^2 to mb
  as = alpha_s(ss);
  M=2.0*(M1*M1+M2*M2);
  if(F1 < F2){
    tmp=F1;
    F1=F2;
    F2=tmp;
  }

  if((F1+F2) == 0){ // gluon-gluon-scattering
    cs=9.0*M_PI*as*as/2.0/ss*(17.0/6.0+6.0*tcut/ss+pow(tcut/ss,2.0)
        +2.0/3.0*pow(tcut/ss,3.0)-2.0*(ss/tcut+ss/(ss+tcut))
        -2.0*log(-ss/tcut-1.0));
  }
  else if((F1*F2) == 0){ // gluon-quark-scattering
    cs=M_PI*as*as/ss*((11.0*ss-2.0*M)*(ss+2.0*tcut)/9.0/ss/ss
        +2.0/9.0*(ss+M)/ss*log((ss+tcut-M/2.0)/(-tcut-M/2.0))
        -(1.0+pow((ss-M)/ss,2.0))*(ss/(tcut+M/2.0)+ss/(ss+tcut-M/2.0)));
  }
  else if(F1 == F2){ // identical quark-quark-scattering
    cs=M_PI*as*as/ss*(8.0/9.0*(ss+2.0*tcut)/ss
        -16.0*((ss-M)/9.0/ss+ss/27.0/(ss-M))*log((ss+tcut-M/2.0)/(-tcut-M/2.0))
        -8.0/9.0*(1.0+pow((ss-M)/ss,2.0))*(ss/(tcut+M/2.0)+ss/(ss+tcut-M/2.0)));
  }
  else if(((F1-F2) == 1) && (((F1/2)*2) == F1)){ // quark-antiquark-scattering
    cs=2.0/27.0*M_PI*as*as/ss*(4.0+3.0*pow(M/ss,2.0)
        +(4.0+6.0*pow(M/ss,2.0))*tcut/ss
        -4.0*(ss-M)*(2.0*ss+M)/ss/ss*log((ss+tcut-M/2.0)/(-tcut-M/2.0))
        -6.0*(1.0+pow((ss-M)/ss,2.0))*(ss/(tcut+M/2.0)+ss/(ss+tcut-M/2.0)));
  }
  else { // quark-quark-scattering
    cs=4.0/9.0*M_PI*as*as/ss*(1.0+2.0*tcut/ss
        -2.0*(ss-M)/ss*log((ss+tcut-M/2.0)/(-tcut-M/2.0))
        -(1.0+pow((ss-M)/ss,2.0))*(ss/(tcut+M/2.0)+ss/(ss+tcut-M/2.0)));
  }

  return cs*factor;
}



double crosection1(int F1, int F2, const double s, const double md2) // dimension 1/GeV�
{

  double as,cs,c;
  int tmp;

  as=alpha_s(s);
  //as=0.3;

  if(F2 < F1)
  {
    tmp=F1;
    F1=F2;
    F2=tmp;
  }

  c=s/4.0+md2;

  if((F1+F2) == 0){ // gluon,gluon -> gluon,gluon; quark,antiquark
    cs=M_PI*as*as*(9.0/8.0*s/(md2*c)+Nflavor/(3.0*s)*log(c/md2));
  }
  else if((F1*F2) == 0){ // gluon,quark -> gluon,quark
    cs=M_PI/2.0*as*as*s/(md2*c);
  }
  else if(((F2-F1) == 1) && (((F2/2)*2) == F2)){
  // quark,antiquark -> quark,antiquark; quark',antiquark'; gluon,gluon
    cs=M_PI*as*as/27.0*(6.0*s/(md2*c)+(Nflavor-1.0)*8.0/s+32.0/s*log(c/md2));
  }
  else { // quark,quark -> quark,quark; quark,quark' -> quark,quark'
    cs=2.0/9.0*M_PI*as*as*s/(md2*c);
  }

  return cs;
}

*/
// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on; 
