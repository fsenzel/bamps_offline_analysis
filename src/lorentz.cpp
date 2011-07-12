//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL: file:///home/fochler/svnrepos/cascade_full/main/trunk/src/lorentz.cpp $
//$LastChangedDate: 2008-02-22 14:53:42 +0100 (Fri, 22 Feb 2008) $
//$LastChangedRevision: 33 $
//$LastChangedBy: fochler $
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------


#include <iostream>
#include <math.h>

#include "lorentz.h"
#include "FPT_compare.h"


// Lorentz transformation
void lorentz(const double beta[4], const double A[4], double LA[4])
{
  double beta2,cc,c0,c1,c2,c3,gama;

  beta2=0.0;
  for(int i=1;i<=3;i++)
    beta2 += beta[i]*beta[i];

  if(beta2 < 1.0e-10)
  {
    for(int k=0;k<=3;k++) 
      LA[k]=A[k];
    return;
  }

  if(beta2 > (1.0-1.0e-10)) 
  {
    std::cout << beta2  <<"  error in lorentz()" << std::endl;
    std::cout << A[0] << " " << A[1] << " " << A[2] << " " << A[3] << std::endl;
    std::cout << beta[1] << " " << beta[2] << " "  << beta[3] << "  " << std::endl;
  }
  else
  {
    gama=1.0/sqrt(1.0-beta2);
      //cout << beta2 << endl;
      //cout << A[0] << " " << A[1] << " " << A[2] << " " << A[3] << endl;
  }

  LA[0]=gama*(A[0]-beta[1]*A[1]-beta[2]*A[2]-beta[3]*A[3]);
  for(int j=1;j<=3;j++){
    cc=(gama-1.0)*beta[j]/beta2;
    c0=-gama*beta[j];
    c1=cc*beta[1];
    c2=cc*beta[2];
    c3=cc*beta[3];
    LA[j]=c0*A[0]+c1*A[1]+c2*A[2]+c3*A[3]+A[j];
  }
}
