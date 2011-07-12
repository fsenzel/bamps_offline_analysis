//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL: file:///home/fochler/svnrepos/cascade_full/main/branches/restructure_branch/src/psf.h $
//$LastChangedDate: 2008-02-05 15:11:34 +0100 (Tue, 05 Feb 2008) $
//$LastChangedRevision: 18 $
//$LastChangedBy: fochler $
//---------------------------------------------
//---------------------------------------------


#ifndef PSF_H
#define PSF_H


void xpGRV(double Q2, double x1, double x2, int nuc1, int nuc2, double F1[9], double F2[9]);

// parameters of the Duke-Owens structure function set 1.  
void xpDO(double Q2, double x1, double x2, int nuc1, int nuc2, double F1[9], double F2[9]);

// Euler-Beta Function, for x>0 and y>0.  
double beta(double x, double y);

//Gamma Function, returns the value ln(gamma(xx)) for xx>0.
float gammln(float xx);


#endif
