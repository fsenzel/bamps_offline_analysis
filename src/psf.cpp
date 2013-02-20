//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------
//---------------------------------------------


#include <math.h>
#include "psf.h"

//
// parameters of the Glück-Reya-Vogt structure function
// for a proton     Z.Phys C 67, 433-447(1995)
//



#define param(a,s) (a[0]+a[1]*s+a[2]*s*s)
#define param1(a,s) (a[0]+a[1]*s+a[2]*s*s+a[3]*s*s*s)
#define param2(a,s) (a[0]+a[1]*sqrt(s)+a[2]*s)
#define B(a,x) a[0]*pow(double(x),double(a[1]))*(1.0+a[2]*pow(double(x),double(a[3]))+a[4]*x+a[5]*pow(double(x),1.5))*pow(double(1.0-x),double(a[6]))
#define C(a,s,x) (pow(double(x),double(a[0]))*(a[1]+a[2]*x+a[3]*x*x)*pow(double(log(1.0/x)),double(a[4]))+pow(double(s),double(a[5]))*exp(-a[6]+sqrt(fabs(a[7]*pow(double(s),double(a[8]))*log(1.0/x)))))*pow(double(1.0-x),double(a[9]))
#define D(a,s,x) pow(double(s),double(a[0]))/pow(double(log(1.0/x)),double(a[1]))*(1.0+a[6]*sqrt(x)+a[7]*x)*pow(double(1.0-x),double(a[2]))*exp(-a[3]+sqrt(fabs(a[4]*pow(double(s),double(a[5]))*log(1.0/x))))




void xpGRV ( double Q2, double x1, double x2, int nuc1, int nuc2, double F1[9], double F2[9] )
{
// nuc1 and nuc2 indicate the nucleons are either porton(1) or neutron(-1).

  const double mu2 = 0.23;
  const double lambda2 = 0.232 * 0.232; //GeVÂ²
  double s, tmp;
  double xpp[6][10];
  static double set_u[7][3] = {{ 2.284, 0.802, 0.055},
    { 0.590, -0.024, 0.000},
    { -0.499, -0.138, -0.076},
    { 0.131, 0.063, 0.000},
    { 0.213, 2.669, -0.728},
    { 8.854, -9.135, 1.979},
    { 2.997, 0.753, -0.076}
  },
  set_d[7][3] = {{ 0.371, 0.083, 0.039},
    { 0.376, 0.000, 0.000},
    { -0.509, 3.310, -1.248},
    { 0.486, 0.062, 0.000},
    { 12.41, -10.52, 2.267},
    { 6.373, -6.208, 1.418},
    { 3.691, 0.799, -0.071}
  },
  set_du[7][3] = {{ 0.082, 0.014, 0.008},
    { 0.409, -0.005, 0.000},
    { -38.07, 36.13, -0.656},
    { 0.799, 0.071, 0.000},
    { 90.31, -74.15, 7.645},
    { 0.000, 0.000, 0.000},
    { 7.486, 1.217, -0.159}
  },
  set_G[9][3] = {{ 1.742, -0.930, 0.000},
    { 7.486, -2.185, 0.000},
    { 16.69, -22.74, 5.779},
    { -25.59, 29.71, -7.296},
    { 0.000, 0.000, -0.399},
    { 0.524, 0.000, 0.000},
    { 0.807, 2.005, 0.000},
    { 3.841, 0.316, 0.000},
    { 1.088, 0.000, 0.000}
  },
  set_GD[4] = { 2.792, 2.215, 0.422, -0.104},
  set_ud[10][3] = {{ 0.410, -0.232, 0.000},
    { 0.890, -0.140, 0.000},
    { -0.981, 0.000, 0.000},
    { 0.320, 0.683, 0.000},
    { 0.534, -0.457, 0.000},
    { 1.451, 0.000, 0.000},
    { 4.119, 1.713, 0.000},
    { 0.682, 2.978, 0.000},
    { 0.271, 0.000, 0.000},
    { 4.752, 1.164, 0.286}
  },
  set_s[6][3] = {{ 0.914, 0.000, 0.000},
    { 1.798, -0.596, 0.000},
    { 6.379, -0.350, 0.142},
    { 3.981, 1.638, 0.000},
    { 6.402, 0.000, 0.000},
    { 0.577, 0.000, 0.000}
  },
  set_sAB[2][3] = {{ -5.548, 3.669, -0.616},
    { 18.92, -16.73, 5.168}
  };

  s = log ( log ( Q2 / lambda2 ) / log ( mu2 / lambda2 ) );

  //parameters for xu_v, xd_v, x(dbar-ubar)
  for ( int i = 0; i < 7; i++ )
  {
    xpp[0][i] = param ( set_u[i], s );
    xpp[1][i] = param ( set_d[i], s );
    xpp[2][i] = param ( set_du[i], s );
  }

  //parameters for xG
  for ( int j = 0; j < 9; j++ )
  {
    xpp[3][j] = param ( set_G[j], s );
  }
  xpp[3][9] = param1 ( set_GD, s );

  //parameters for x(ubar+dbar)
  for ( int i = 0; i < 10; i++ )
  {
    xpp[4][i] = param ( set_ud[i], s );
  }

  //parameters for xs
  for ( int j = 0; j < 6; j++ )
  {
    xpp[5][j] = param ( set_s[j], s );
  }
  for ( int i = 0; i < 2; i++ )
  {
    xpp[5][i + 6] = param2 ( set_sAB[i], s );
  }

  // [0]g,[1]u,[2]ubar,[3]d,[4]dbar,[5]s,[6]sbar,[7]c,[8]cbar

  F1[2] = ( -B ( xpp[2], x1 ) + C ( xpp[4], s, x1 ) ) / 2.0;
  F1[4] = ( B ( xpp[2], x1 ) + C ( xpp[4], s, x1 ) ) / 2.0;
  F1[1] = B ( xpp[0], x1 );
  F1[3] = B ( xpp[1], x1 );
  if ( nuc1 == -1 )
  {
    tmp = F1[1];
    F1[1] = F1[3];
    F1[3] = tmp;
  }
  F1[1] = F1[1] + F1[2];
  F1[3] = F1[3] + F1[4];
  F1[5] = F1[6] = D ( xpp[5], s, x1 );
  F1[7] = F1[8] = 0.0;
  F1[0] = C ( xpp[3], s, x1 );

  F2[2] = ( -B ( xpp[2], x2 ) + C ( xpp[4], s, x2 ) ) / 2.0;
  F2[4] = ( B ( xpp[2], x2 ) + C ( xpp[4], s, x2 ) ) / 2.0;
  F2[1] = B ( xpp[0], x2 );
  F2[3] = B ( xpp[1], x2 );
  if ( nuc2 == -1 )
  {
    tmp = F2[1];
    F2[1] = F2[3];
    F2[3] = tmp;
  }
  F2[1] = F2[1] + F2[2];
  F2[3] = F2[3] + F2[4];
  F2[5] = F2[6] = D ( xpp[5], s, x2 );
  F2[7] = F2[8] = 0.0;
  F2[0] = C ( xpp[3], s, x2 );
}





//
// parameters of the Duke-Owens structure function set 1.
// Phys.Rev. D 30 (1984) 49.
//
#define param(a,s) (a[0]+a[1]*s+a[2]*s*s)
#define A(a,x) a[0]*pow(double(x),double(a[1]))*pow(double(1.-x),double(a[2]))*(1.+a[3]*x+a[4]*x*x+a[5]*pow(double(x),3))

void xpDO ( double Q2, double x1, double x2, int nuc1, int nuc2, double F1[9], double F2[9] )
{
// nuc1 and nuc2 indicate the nucleons are either porton(1) or neutron(-1).

  const double Q02 = 4.0; //GeVÂ²
  const double lambda2 = 0.04; //GeVÂ²
  double s, tmp;
  double xpp[5][6];
  static double set1_ud[3][3] = {{ 0.419, 0.004, -0.007},
    { 3.460, 0.724, -0.066},
    { 4.400, -4.860, 1.330}
  },
  set1_d[3][3] = {{ 0.763, -0.237, 0.026},
    { 4.000, 0.627, -0.019},
    { 0.000, -0.421, 0.033}
  },
  set1_S[6][3] = {{ 1.265, -1.132, 0.293},
    { 0.000, -0.372, -0.029},
    { 8.050, 1.590, -0.153},
    { 0.000, 6.310, -0.273},
    { 0.000, -10.50, -3.170},
    { 0.000, 14.70, 9.800}
  },
  set1_c[6][3] = {{ 0.000, 0.135, -0.075},
    { -0.036, -0.222, -0.058},
    { 6.350, 3.260, -0.909},
    { 0.000, -3.030, 1.500},
    { 0.000, 17.40, -11.30},
    { 0.000, -17.90, 15.60}
  },
  set1_G[6][3] = {{ 1.560, -1.710, 0.638},
    { 0.000, -0.949, 0.325},
    { 6.000, 1.440, -1.050},
    { 9.000, -7.190, 0.255},
    { 0.000, -16.50, 10.90},
    { 0.000, 15.30, -10.10}
  };

  s = log ( log ( Q2 / lambda2 ) / log ( Q02 / lambda2 ) );

  //parameters for x(u_v+d_v)
  for ( int i = 1; i < 4; i++ )
  {
    xpp[0][i] = param ( set1_ud[i - 1], s );
  }
  xpp[0][0] = 3. / ( beta ( xpp[0][1], xpp[0][2] + 1. ) * ( 1. + xpp[0][3] * xpp[0][1] /
                     ( xpp[0][1] + xpp[0][2] + 1. ) ) );
  xpp[0][4] = 0.;
  xpp[0][5] = 0.;

  //parameters for xd_v
  for ( int j = 1; j < 4; j++ )
  {
    xpp[1][j] = param ( set1_d[j - 1], s );
  }
  xpp[1][0] = 1. / ( beta ( xpp[1][1], xpp[1][2] + 1. ) * ( 1. + xpp[1][3] * xpp[1][1] /
                     ( xpp[1][1] + xpp[1][2] + 1. ) ) );
  xpp[1][4] = 0.;
  xpp[1][5] = 0.;

  //parameters for xS, S=2(u_s+d_s+s_s)
  for ( int i = 0; i < 6; i++ )
  {
    xpp[2][i] = param ( set1_S[i], s );
  }

  //parameters for xc
  for ( int j = 0; j < 6; j++ )
  {
    xpp[3][j] = param ( set1_c[j], s );
  }

  //parameters for xG
  for ( int i = 0; i < 6; i++ )
  {
    xpp[4][i] = param ( set1_G[i], s );
  }

  // [0]g,[1]u,[2]ubar,[3]d,[4]dbar,[5]s,[6]sbar,[7]c,[8]cbar
  F1[2] = F1[4] = F1[5] = F1[6] = A ( xpp[2], x1 ) / 6.;
  F1[3] = A ( xpp[1], x1 ) + F1[4];
  F1[1] = A ( xpp[0], x1 ) - F1[3] + 2.*F1[2];
  if ( nuc1 == -1 )
  {
    tmp = F1[1];
    F1[1] = F1[3];
    F1[3] = tmp;
  }
  F1[7] = F1[8] = A ( xpp[3], x1 );
  F1[0] = A ( xpp[4], x1 );

  F2[2] = F2[4] = F2[5] = F2[6] = A ( xpp[2], x2 ) / 6.;
  F2[3] = A ( xpp[1], x2 ) + F2[4];
  F2[1] = A ( xpp[0], x2 ) - F2[3] + 2.*F2[2];
  if ( nuc2 == -1 )
  {
    tmp = F2[1];
    F2[1] = F2[3];
    F2[3] = tmp;
  }
  F2[7] = F2[8] = A ( xpp[3], x2 );
  F2[0] = A ( xpp[4], x2 );
}


// Euler-Beta Function, for x>0 and y>0.
double beta ( double x, double y )
{
  return double ( exp ( gammln ( x ) + gammln ( y ) - gammln ( x + y ) ) );
}

//Gamma Function, returns the value ln(gamma(xx)) for xx>0.
float gammln ( float xx )
{

  double x, y, tmp, ser;
  static double cof[6] = {76.18009172947146, -86.50532032941677,
                          24.01409824083091, -1.231739572450155,
                          0.1208650973866179e-2, -0.5395239384953e-5
                         };

  y = x = xx;
  tmp = x + 5.5;
  tmp -= ( x + 0.5 ) * log ( tmp );
  ser = 1.000000000190015;
  for ( int j = 0; j <= 5; j++ )
    ser += cof[j] / ++y;
  return -tmp + log ( 2.5066282746310005 * ser / x );

}
