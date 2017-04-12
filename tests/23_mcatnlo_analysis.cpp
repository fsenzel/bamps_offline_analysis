//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/trunk/tests/21_checkinit.cpp $
//$LastChangedDate: 2014-03-21 13:29:33 +0100 (Fri, 21 Mar 2014) $
//$LastChangedRevision: 1641 $
//$LastChangedBy: gallmei $
//---------------------------------------------
//---------------------------------------------

#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
// #include <stdlib.h>
// #include <exception>
// #include <stdint.h>
// #include <execinfo.h>
// #include <signal.h>
#include <vector>

#include "hadronization_hq.h"
#include "mesonDecay.h"
#include "particleOffline.h"
#include "analysis.h"

using namespace std;

vector<ParticleOffline> partcl;
vector<ParticleOffline> partcl_dmeson;
vector<ParticleOffline> partcl_electron;

int nmb_electron_per_hq;
const int nmb_electron_per_hq_charm = 1;
const int nmb_electron_per_hq_bottom = 1;
bool muonsInsteadOfElectrons_global = false;
bool nonPromptJpsiInsteadOfElectrons_global = false;

// const double pp_inelast_cross_section = 42.0; //mb RHIC
const double pp_inelast_cross_section = 64.0; //mb // LHC 2.76 TeV
// const double pp_inelast_cross_section = 69.4; //mb // LHC 7 TeV, http://arxiv.org/abs/1104.0326

// MC@NLO (values stem from MC_NLO.log
// charm

// RHIC @ 200 GeV
// const double ccb_cross_section = 0.2318E+09; //pb for pairs
// const double ccb_cross_section = 5.34E+008; //pb for pairs, M=1.3  renormalisation=factorization scale = 0.48
// const double ccb_cross_section = 5.25E+008; //pb for pairs, M=1.5  renormalisation=factorization scale = 0.35
// const double ccb_cross_section = 0.4936E+09; //pb  2.000D+02     1.300D+00 0.50 0.50 0.50 0.4936D+09 2.D+06
// const double ccb_cross_section = 1.15E+009; //pb  2.000D+02     1.300D+00 0.50 0.50 0.50 0.1148D+10 3.D+06 EPS08
// const double ccb_cross_section = 0.5399E+09; //pb   2.000D+02     1.300D+00 0.67 0.67 0.67 0.5399D+09 1.D+06
// const double ccb_cross_section = 0.5183E+09; //pb   2.000D+02     1.300D+00 0.67 0.67 0.67 0.5183D+09 1.D+06 EPS08
// const double ccb_cross_section = 0.5747E+09; //pb    2.000D+02     1.300D+00 0.65 0.65 0.65 0.5747D+09 2.D+06
// const double ccb_cross_section = 0.5546E+09; //pb   2.000D+02     1.300D+00 0.65 0.65 0.65 0.5546D+09 2.D+06 EPS08
// const double ccb_cross_section = 0.5556E+09; //pb    2.000D+02     1.300D+00 0.65 0.65 0.65 0.5556D+09 2.D+06 EKS98

// LHC @ 2.76 TeV
// const double ccb_cross_section = 0.2574E+10; //pb   2.760D+03     1.300D+00 0.65 0.65 0.65 0.2574D+10 1.D+07
// const double ccb_cross_section = 0.1549E+10; //pb   2.760D+03     1.300D+00 0.65 0.65 0.65 0.1549D+10 1.D+07 EPS08
// const double ccb_cross_section = 0.1552E+10; //pb   2.760D+03     1.300D+00 0.65 0.65 0.65 0.1552D+10 1.D+07 EKS98
// const double ccb_cross_section = 0.2708E+10; //pb   2.760D+03     1.300D+00 1.00 1.00 1.00 0.2708D+10 1.D+07
// const double ccb_cross_section = 0.4103E+09; //pb   2.760D+03     1.300D+00 1.00 1.00 1.00 0.4103D+09 5.D+06 shadowing??
// const double ccb_cross_section = 0.4363E+10; //pb   2.760D+03     1.300D+00 0.50 0.50 0.50 0.4363D+10 3.D+07
const double ccb_cross_section = 0.1049E+11; //pb   2.760D+03     1.300D+00 0.40 0.40 0.40 0.1049D+11 7.D+07

// LHC @ 7 TeV
// const double ccb_cross_section = 0.4224E+10; //pb    7.000D+03     1.300D+00 0.65 0.65 0.65 0.4224D+10 2.D+07
// const double ccb_cross_section = 0.5192E+10; //pb     7.000D+03     1.300D+00 1.00 1.00 1.00 0.5192D+10 2.D+07
// const double ccb_cross_section = 0.5624E+10; //pb      7.000D+03     1.300D+00 1.20 1.20 1.20 0.5624D+10 2.D+07
// const double ccb_cross_section = 0.6026E+10; //pb    7.000D+03     1.300D+00 1.50 1.50 1.50 0.6026D+10 2.D+07
// const double ccb_cross_section = 0.4400E+10; //pb    7.000D+03     1.300D+00 0.80 0.80 0.80 0.4400D+10 2.D+07

// PHENIX
// const double ccb_cross_section = 0.544E+09; //pb for pairs

// MC@NLO
// bottom

// RHIC @ 200 GeV
// const double bbb_cross_section = 0.1731E+07; //pb for pairs
// const double bbb_cross_section = 3.27E+006; //pb for pairs, M=4.7  renormalisation=factorization scale = 0.35
// const double bbb_cross_section = 3.65E+006; //pb for pairs, M=4.6   renormalisation=factorization scale = 0.35   
// const double bbb_cross_section = 3.06E+006; //pb for pairs, M=4.6   renormalisation=factorization scale = 0.48  
// const double bbb_cross_section = 0.2992E+07; //pb 2.000D+02     4.600D+00 0.50 0.50 0.50 0.2992D+07 4.D+03
// const double bbb_cross_section = 0.3367E+07; //pb  2.000D+02     4.600D+00 0.40 0.40 0.40 0.3367D+07 5.D+03
// const double bbb_cross_section = 0.3358E+07; //pb  2.000D+02        4.600D+00 0.40 0.40 0.40 0.3358D+07 5.D+03 EPS08
// const double bbb_cross_section = 0.3308E+07; //pb   2.000D+02     4.600D+00 0.40 0.40 0.40 0.3308D+07 5.D+03 EKS98


// LHC @ 2.76 TeV
// const double bbb_cross_section = 0.1196E+09; //pb   2.760D+03     4.600D+00 0.40 0.40 0.40 0.1196D+09 2.D+05
// const double bbb_cross_section = 0.3822E+08; //pb   2.760D+03     4.600D+00 0.40 0.40 0.40 0.3822D+08 1.D+05 EPS08
// const double bbb_cross_section = 0.3808E+08; //pb   2.760D+03     4.600D+00 0.40 0.40 0.40 0.3808D+08 1.D+05 EKS98
// const double bbb_cross_section = 0.1005E+09; //pb   2.760D+03     4.600D+00 1.00 1.00 1.00 0.1005D+09 2.D+05
// const double bbb_cross_section = 0.7563E+07; //pb   2.760D+03     4.600D+00 1.00 1.00 1.00 0.7563D+07 4.D+04 EPS08


// LHC @ 7 TeV
const double bbb_cross_section = 0.2783E+09; //pb    7.000D+03     4.600D+00 0.40 0.40 0.40 0.2783D+09 6.D+05
// const double bbb_cross_section = 0.2615E+09; //pb    7.000D+03     4.600D+00 1.00 1.00 1.00 0.2615D+09 5.D+05
// const double bbb_cross_section = 0.2518E+09; //pb     7.000D+03     4.600D+00 1.50 1.50 1.50 0.2518D+09 5.D+05
// const double bbb_cross_section = 0.2469E+09; //pb   7.000D+03     4.600D+00 2.00 2.00 2.00 0.2469D+09 5.D+05

//scaled like ccb of phenix
// const double bbb_cross_section = 0.35E+07; //pb for pairs

unsigned int CountLines(std::istream& in)
{
    return std::count(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>(), '\n');
}

void binning_for_plots( string filename_prefix, const vector<ParticleOffline> partcl, const FLAVOR_TYPE ana_type, const FLAVOR_TYPE hq_flavor, const double nmb_QQb);

/**
 * @brief The main routine of BAMPS
 */
int main( int argc, char *argv[] )
{
  if ( !( argc == 3 ) )
  {
    cout << "Usage: ./22_mcatnlo_analysis <flavor> <path> <output_filename_prefix>" << endl;
  }
  {    
    FLAVOR_TYPE flavor = static_cast<FLAVOR_TYPE>( atoi( argv[1] ) );
    string path = argv[2];
    string filename_prefix = argv[3];
    
    stringstream filename;
    
    double nmb_QQb;    
    // number of heavy quark pairs in one pp collision
    if( flavor == charm )
    {
        nmb_QQb = ( ccb_cross_section*pow(10.0,-12) )/ ( pp_inelast_cross_section*pow(10.0,-3) ) ;
        nmb_electron_per_hq = nmb_electron_per_hq_charm;
        filename << path << "/charm_2.dat";
    }
    else if( flavor == bottom )
    {
        nmb_QQb = ( bbb_cross_section*pow(10.0,-12) )/ ( pp_inelast_cross_section*pow(10.0,-3) ) ;
        nmb_electron_per_hq = nmb_electron_per_hq_bottom;
        filename << path << "/bottom_2.dat";
    }
    
    cout << "file: " << filename.str() << endl;
    
    cout << "QQbar pairs heavy ion = " << nmb_QQb * 1200 << endl;
    
    // number of single hq + antihq in one pp collision
    nmb_QQb = nmb_QQb * 2;

    // read file
    ifstream read(filename.str());
    
    // count number of particles in file
    int numberInInputFile = CountLines(read);
    int number_electrons = numberInInputFile * nmb_electron_per_hq;
    
    cout << "number of particles in inputfile: " << numberInInputFile << endl;
    
    // move cursor to beginning of file
    read.seekg (0, ios::beg);

    partcl.resize( numberInInputFile );
    partcl_dmeson.resize( numberInInputFile );
    partcl_electron.resize( numberInInputFile );
    
    double read_dummy;
    for (int i = 1; i <= numberInInputFile; i++)
    {
      // structure of file
      // ID, hard event, flavour  Energy  Px  Py  Pz  M, flavor pdg,, flavor herwig, status
      int flav;
      read >> read_dummy >> read_dummy >> flav >> partcl[i].Mom.E() >> partcl[i].Mom.Px() >> partcl[i].Mom.Py() >> partcl[i].Mom.Pz() >> partcl[i].m >> read_dummy >> read_dummy >> read_dummy;
        
      if(flav == 4)
        partcl[i].FLAVOR = charm; //charm
      else if(flav == -4)
        partcl[i].FLAVOR = anti_charm; // anti charm
      else if(flav == 5)
        partcl[i].FLAVOR = bottom; // bottom
      else if(flav == -5)
        partcl[i].FLAVOR = anti_bottom; // anti bottom
      else
        cout << "wrong flavor: " << flav << endl;
    }
    read.close();

    hadronization_hq ppHadronization;
    partcl_dmeson = ppHadronization.heavyQuarkFragmentation( partcl );

//     mesonDecay ppMesonDecay( nmb_electron_per_hq, muonsInsteadOfElectrons_global, nonPromptJpsiInsteadOfElectrons_global );
//     partcl_electron = ppMesonDecay.decayToElectronsPythia( partcl_dmeson );
  
    binning_for_plots( filename_prefix, partcl, flavor, flavor, nmb_QQb );
    binning_for_plots( filename_prefix, partcl_dmeson, dmeson_gen, flavor, nmb_QQb );
    binning_for_plots( filename_prefix, partcl_electron, electron, flavor, nmb_QQb );
  
    return EXIT_SUCCESS;
  }
}


void binning_for_plots( string filename_prefix, const vector<ParticleOffline> partcl, const FLAVOR_TYPE ana_type, const FLAVOR_TYPE hq_flavor, const double nmb_QQb)
{
  string name,filename;
  
  if( hq_flavor == charm)
  {
    if( ana_type == charm)
      name = "charm";
    else if( ana_type == dmeson_gen)
      name = "dmeson";
    else if( ana_type == electron)
      name = "c_electron";
    else
      cout << "error" << endl;
  }
  else if( hq_flavor == bottom)
  {
    if( ana_type == bottom)
      name = "bottom";
    else if( ana_type == dmeson_gen)
      name = "bmeson";
    else if( ana_type == electron)
      name = "b_electron";
    else
      cout << "error" << endl;
  }
  else
    cout << "error" << endl;
  

  const int n = 150; // number of Bins
  int dummy;
  double dy,dpt,y_out,pt_out, pt;
  double y_min, y_max,y;
  
  // compute dN/dY
  y_min = -10.0;
  y_max = 10.0;
  
  dy = (y_max-y_min)/n;
  int yBins[n+1]; // bins of t
  
  for(int j=0;j<n+1;j++)
    yBins[j]=0;
  
  for(int i = 0; i < partcl.size(); i++)
  {
    y = 0.5*log( ( partcl[i].Mom.E()+partcl[i].Mom.Pz() )/( partcl[i].Mom.E()-partcl[i].Mom.Pz() ));
    if (y <= y_max && y >= y_min)
    {
      dummy = int((y-y_min)/dy);
      yBins[dummy] += 1;
    }
  } 
  
  // dndy normalized to one if integrated, multiplied with number of charm pairs
  filename = "output/"+filename_prefix+"_"+name+"_dndy.dat";
  ofstream print1(filename.c_str());

  for(int k=0;k<n+1;k++)
  {
    print1.width(15);
    y_out = double(k)*dy+y_min;
    print1 << y_out;
    print1.width(15);
    print1 << double(yBins[k])/partcl.size()/dy*nmb_QQb << endl;
  }

  // compute dE_t/dY
  y_min = -10.0;
  y_max = 10.0;
  
  dy = (y_max-y_min)/n;
  double eBins[n+1]; // bins of t
  
  for(int j=0;j<n+1;j++)
    eBins[j]=0.0;
  
  for(int i = 0; i < partcl.size(); i++)
  {
    y = 0.5*log( ( partcl[i].Mom.E()+partcl[i].Mom.Pz() )/( partcl[i].Mom.E()-partcl[i].Mom.Pz() ));
    if (y <= y_max && y >= y_min)
    {
      dummy = int((y-y_min)/dy);
      eBins[dummy] += sqrt(partcl[i].Mom.Px2() + partcl[i].Mom.Py2() + partcl[i].m * partcl[i].m);
    }
  } 
  
  filename = "output/"+filename_prefix+"_"+name+"_dedy.dat";
  ofstream print1e(filename.c_str());

  for(int k=0;k<n+1;k++)
  {
    print1e.width(15);
    y_out = double(k)*dy+y_min;
    print1e << y_out;
    print1e.width(15);
    print1e << double(eBins[k])/partcl.size()/dy*nmb_QQb << endl;
  }
  
  // compute dN/dpt/dY
  y_min = -0.5;
  y_max = 0.5;
  double pt_min = 0.0;
  double pt_max = 30.0;

  dpt = (pt_max-pt_min)/n;
  dy = y_max-y_min;
  
  //hard
  int ptBins[n+1]; // Bins of t
  
  for(int j=0;j<n+1;j++)
    ptBins[j]=0;
  
  for(int i = 0; i < partcl.size(); i++)
  {
    y = 0.5*log( ( partcl[i].Mom.E()+partcl[i].Mom.Pz() )/( partcl[i].Mom.E()-partcl[i].Mom.Pz() ));
    pt = partcl[i].Mom.Perp();
    
    if (y<=y_max && y>=y_min && pt<=pt_max && pt>=pt_min)
    {
      dummy = int((pt-pt_min)/dpt);
      ptBins[dummy] += 1;
    }
  } 
  
  filename = "output/"+filename_prefix+"_"+name+"_dndptdy.dat";
  ofstream print2(filename.c_str());

  for(int k=0;k<n+1;k++)
  {
    print2.width(15);
    pt_out = double(k)*dpt+pt_min+dpt/2.0;
    print2 << pt_out;
    print2.width(15);
    print2 << 1.0/pt_out*double(ptBins[k])/partcl.size()/dpt/dy*nmb_QQb << endl;
  }
  
  // compute dN/dpt/dY in forward rapidity
  y_min = 2.5;
  y_max = 4.0;
  pt_min = 0.0;
  pt_max = 30.0;
  
  dpt = (pt_max-pt_min)/n;
  dy = 2.0 * (y_max-y_min);
  
  for(int j=0;j<n+1;j++)
    ptBins[j]=0;
  
  for(int i = 0; i < partcl.size(); i++)
  {
    y = 0.5*log( ( partcl[i].Mom.E()+partcl[i].Mom.Pz() )/( partcl[i].Mom.E()-partcl[i].Mom.Pz() ));
    pt = partcl[i].Mom.Perp();
    
    if (fabs(y)<=y_max && fabs(y)>=y_min && pt<=pt_max && pt>=pt_min)
    {
      dummy = int((pt-pt_min)/dpt);
      ptBins[dummy] += 1;
    }
  } 
  
  filename = "output/"+filename_prefix+"_"+name+"_dndptdy_forwardRap.dat";
  ofstream printForwardRap(filename.c_str());

  for(int k=0;k<n+1;k++)
  {
    printForwardRap.width(15);
    pt_out = double(k)*dpt+pt_min+dpt/2.0;
    printForwardRap << pt_out;
    printForwardRap.width(15);
    printForwardRap << 1.0/pt_out*double(ptBins[k])/partcl.size()/dpt/dy*nmb_QQb << endl;
  }

  
  // compute dN/dpt/d eta
  double eta_min = -0.5;
  double eta_max = 0.5;
  pt_min = 0.0;
  pt_max = 30.0;
  
  dpt = (pt_max-pt_min)/n;
  double deta = eta_max-eta_min;
  double eta,pp;
  
  for(int j=0;j<n+1;j++)
    ptBins[j]=0;
  
  for(int i = 0; i < partcl.size(); i++)
  {
    pp = sqrt( pow( partcl[i].Mom.E(), 2.0 ) - pow( partcl[i].m, 2.0 ) );
    // pseudorapidity
    eta = 0.5 * log(( pp + partcl[i].Mom.Pz() ) / ( pp - partcl[i].Mom.Pz() ) );
    pt = partcl[i].Mom.Perp();
    
    if (eta<=eta_max && eta>=eta_min && pt<=pt_max && pt>=pt_min)
    {
      dummy = int((pt-pt_min)/dpt);
      ptBins[dummy] += 1;
    }
  } 
  
  filename = "output/"+filename_prefix+"_"+name+"_dndptdeta.dat";
  ofstream print3(filename.c_str());

  for(int k=0;k<n+1;k++)
  {
    print3.width(15);
    pt_out = double(k)*dpt+pt_min+dpt/2.0;
    print3 << pt_out;
    print3.width(15);
    print3 << 1.0/pt_out*double(ptBins[k])/partcl.size()/dpt/deta*nmb_QQb << endl;
  }
  
  // do CMS output for non-prompt Jpsi
  if( nonPromptJpsiInsteadOfElectrons_global )
  {
    // compute dN/dpt/dY in forward rapidity
    y_min = 0.0;
    y_max = 1.2;
    pt_min = 0.0;
    pt_max = 30.0;
    
    dpt = (pt_max-pt_min)/n;
    dy = 2.0 * (y_max-y_min);
    
    for(int j=0;j<n+1;j++)
      ptBins[j]=0;
    
    for(int i = 0; i < partcl.size(); i++)
    {
      y = 0.5*log( ( partcl[i].Mom.E()+partcl[i].Mom.Pz() )/( partcl[i].Mom.E()-partcl[i].Mom.Pz() ));
      pt = partcl[i].Mom.Perp();
      
      if (fabs(y)<=y_max && fabs(y)>=y_min && pt<=pt_max && pt>=pt_min)
      {
        dummy = int((pt-pt_min)/dpt);
        ptBins[dummy] += 1;
      }
    } 
    
    filename = "output/"+filename_prefix+"_"+name+"_dndptdy_nonPromptJpsi_0_12.dat";
    ofstream printNPJ_012(filename.c_str());

    for(int k=0;k<n+1;k++)
    {
      printNPJ_012.width(15);
      pt_out = double(k)*dpt+pt_min+dpt/2.0;
      printNPJ_012 << pt_out;
      printNPJ_012.width(15);
      printNPJ_012 << 1.0/pt_out*double(ptBins[k])/partcl.size()/dpt/dy*nmb_QQb << endl;
    }
    
    // compute dN/dpt/dY in forward rapidity
    y_min = 1.2;
    y_max = 1.6;
    pt_min = 0.0;
    pt_max = 30.0;
    
    dpt = (pt_max-pt_min)/n;
    dy = 2.0 * (y_max-y_min);
    
    for(int j=0;j<n+1;j++)
      ptBins[j]=0;
    
    for(int i = 0; i < partcl.size(); i++)
    {
      y = 0.5*log( ( partcl[i].Mom.E()+partcl[i].Mom.Pz() )/( partcl[i].Mom.E()-partcl[i].Mom.Pz() ));
      pt = partcl[i].Mom.Perp();
     
      if (fabs(y)<=y_max && fabs(y)>=y_min && pt<=pt_max && pt>=pt_min)
      {
        dummy = int((pt-pt_min)/dpt);
        ptBins[dummy] += 1;
      }
    } 
    
    filename = "output/"+filename_prefix+"_"+name+"_dndptdy_nonPromptJpsi_12_16.dat";
    ofstream printNPJ_1216(filename.c_str());

    for(int k=0;k<n+1;k++)
    {
      printNPJ_1216.width(15);
      pt_out = double(k)*dpt+pt_min+dpt/2.0;
      printNPJ_1216 << pt_out;
      printNPJ_1216.width(15);
      printNPJ_1216 << 1.0/pt_out*double(ptBins[k])/partcl.size()/dpt/dy*nmb_QQb << endl;
    }
    
    
    // compute dN/dpt/dY in forward rapidity
    y_min = 1.6;
    y_max = 2.4;
    pt_min = 0.0;
    pt_max = 30.0;
    
    dpt = (pt_max-pt_min)/n;
    dy = 2.0 * (y_max-y_min);
    
    for(int j=0;j<n+1;j++)
      ptBins[j]=0;
    
    for(int i = 0; i < partcl.size(); i++)
    {
      y = 0.5*log( ( partcl[i].Mom.E()+partcl[i].Mom.Pz() )/( partcl[i].Mom.E()-partcl[i].Mom.Pz() ));
      pt = partcl[i].Mom.Perp();
     
      if (fabs(y)<=y_max && fabs(y)>=y_min && pt<=pt_max && pt>=pt_min)
      {
        dummy = int((pt-pt_min)/dpt);
        ptBins[dummy] += 1;
      }
    } 
    
    filename = "output/"+filename_prefix+"_"+name+"_dndptdy_nonPromptJpsi_16_24.dat";
    ofstream printNPJ_1624(filename.c_str());

    for(int k=0;k<n+1;k++)
    {
      printNPJ_1624.width(15);
      pt_out = double(k)*dpt+pt_min+dpt/2.0;
      printNPJ_1624 << pt_out;
      printNPJ_1624.width(15);
      printNPJ_1624 << 1.0/pt_out*double(ptBins[k])/partcl.size()/dpt/dy*nmb_QQb << endl;
    }
  }
}
