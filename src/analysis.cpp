//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------
//---------------------------------------------


#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <sstream>
#include <string>
#include <stdio.h> // for getenv()
#include <stdlib.h> // for getenv()

#include "analysis.h"
#include "binning.h"
#include "binning2.h"
#include "FPT_compare.h"
#include "random.h"

using namespace std;
using namespace ns_casc;

extern int IX, IY;

namespace ns_heavy_quarks
{
  extern int jpsi_dissociation;
  extern int jpsi_dissociation_from_temperature;
  extern int jpsicreation;
  extern int charmAnnihil;
}


analysis::analysis( config* const c ):
    theConfig( c ),
    rings( c->getRingNumber(), c->getCentralRingRadius(), c->getDeltaR() ),
    centralRingsCopyFromCascade( c->getRingNumber(), c->getCentralRingRadius(), c->getDeltaR() ),
    filename_prefix( c->getStandardOutputDirectoryName() + "/" + c->getJobName() ),
    v2output( c->isV2RAAoutput() ),
    dndyOutput( c->isDndyOutput() ),
    particleCorrelationsOutput( false ),
    hadronization_hq( c->isHadronizationHQ() ),
    mesonDecay( c->isMesonDecay() ),
//     charmTestJet( theConfig->isCharmTestJet() ),
    v2outputIntermediateSteps( c->isV2RAAoutputIntermediateSteps() ),
    studyParticleOutput( c->doOutput_detailedParticleOutput() ),
    outputScheme( c->getOutputScheme() )
{
  //---get time and date---
  time( &start );
  //-----------------------
  
  cout << "Output Scheme: " << outputScheme << endl;
  handle_output_studies( outputScheme );
   
  // to create the tsteps use such a bash script:
  // for (( I=50; $I <= 80; I++ )); do   J=$(echo "scale=1; ($I+1)/10" | bc); echo "tstep[$I]=$J;" ; done
  //for (( I=0; $I <= 800; I++ )); do   J=$(echo "scale=2; ($I+15)/100" | bc); echo "tstep[$I]=$J;" ; done
  
    
  if( studyJpsi ) // jpsi evolution: more timesteps
  {
    tstep[0]=.005;
    tstep[1]=.006; // choose to be small that no analysis is performed, just start for all files at 0.2 in next step
    tstep[2]=.20;
    tstep[3]=.25;
    tstep[4]=.30;
    tstep[5]=.35;
    tstep[6]=.40;
    tstep[7]=.45;
    tstep[8]=.50;
    tstep[9]=.55;
    tstep[10]=.60;
    tstep[11]=.65;
    tstep[12]=.70;
    tstep[13]=.75;
    tstep[14]=.80;
    tstep[15]=.85;
    tstep[16]=.90;
    tstep[17]=.95;
    tstep[18]=1.00;
    tstep[19]=1.1;
    tstep[20]=1.2;
    tstep[21]=1.3;
    tstep[22]=1.4;
    tstep[23]=1.5;
    tstep[24]=1.6;
    tstep[25]=1.7;
    tstep[26]=1.8;
    tstep[27]=1.9;
    tstep[28]=2.0;
    tstep[29]=2.1;
    tstep[30]=2.2;
    tstep[31]=2.3;
    tstep[32]=2.4;
    tstep[33]=2.5;
    tstep[34]=2.6;
    tstep[35]=2.7;
    tstep[36]=2.8;
    tstep[37]=2.9;
    tstep[38]=3.0;
    tstep[39]=3.3;
    tstep[40]=3.7;
    tstep[41]=4.0;
    tstep[42]=4.3;
    tstep[43]=4.7;
    tstep[44]=5.0;
    tstep[45]=5.3;
    tstep[46]=5.7;
    tstep[47]=6.0;
    tstep[48]=6.3;
    tstep[49]=6.7;
    tstep[50]=7.0;
    tstep[51]=7.3;
    tstep[52]=7.7;
    tstep[53]=8.0;
    tstep[54]=8.5;
    tstep[55]=9.0;
    tstep[56]=9.5;
    tstep[57]=10.0;
    tstep[58]=10.5;                                                                                                                                                                     
    tstep[59]=11.0;                                                                                                                                                                     
    tstep[60]=11.5;                                                                                                                                                                     
    tstep[61]=12.0;                                                                                                                                                                     
    tstep[62]=12.5;                                                                                                                                                                     
    tstep[63]=13.0;                                                                                                                                                                     
    tstep[64]=13.5;                                                                                                                                                                     
    tstep[65]=14.0;                                                                                                                                                                     
    tstep[66]=14.5;                                                                                                                                                                     
    tstep[67]=15.0;                                                                                                                                                                     
    tstep[68]=15.5;                                                                                                                                                                     
    tstep[69]=16.0;                                                                                                                                                                     
    tstep[70]=16.5;                                                                                                                                                                     
    tstep[71]=17.0;                                                                                                                                                                     
    tstep[72]=17.5;                                                                                                                                                                     
    tstep[73]=18.0;                                                                                                                                                                     
    tstep[74]=18.5;                                                                                                                                                                     
    tstep[75]=19.0;                                                                                                                                                                     
    tstep[76]=19.5;                                                                                                                                                                     
    tstep[77]=20.0;
    tstep[78] = infinity; //fm/c
    nTimeSteps = 79;
    
    if( theConfig->getRuntime() > 20.0 )
    {
      tstep[78]=21;
      tstep[79]=22;
      tstep[80]=23;
      tstep[81]=24;
      tstep[82]=25;
      tstep[83]=26;
      tstep[84]=27;
      tstep[85]=28;
      tstep[86]=29;
      tstep[87]=30;
      tstep[88]=31;
      tstep[89]=32;
      tstep[90]=33;
      tstep[91]=34;
      tstep[92]=35;
      tstep[93]=36;
      tstep[94]=37;
      tstep[95]=38;
      tstep[96]=39;
      tstep[97]=40;
      tstep[98]=42;
      tstep[99]=44;
      tstep[100]=46;
      tstep[101]=48;
      tstep[102]=50;
      tstep[103]=52;
      tstep[104]=54;
      tstep[105]=56;
      tstep[106]=58;
      tstep[107]=60;
      tstep[108]=62;
      tstep[109]=64;
      tstep[110]=66;
      tstep[111]=68;
      tstep[112]=70;
      tstep[113]=72;
      tstep[114]=74;
      tstep[115]=76;
      tstep[116]=78;
      tstep[117]=80;
      tstep[118] = infinity; //fm/c
      nTimeSteps = 119;
    }
  }
  else if (studyPhotons)
  {
    tstep[0]=.15;
    tstep[1]=.2;
    tstep[2]=.3;
    tstep[3]=.4;
    tstep[4]=.5;
    tstep[5]=.6;
    tstep[6]=.7;
    tstep[7]=.8;
    tstep[8]=.9;
    tstep[9]=1.0;
    tstep[10] = 2.0;
    tstep[11] = 3.0;
    tstep[12] = 4.0;
    tstep[13] = 5.0;
    tstep[14] = 6.0;
    tstep[15] = 7.0;
    tstep[16] = 8.0;
    tstep[17] = 9.0;
    tstep[18] = 10.0;
    tstep[19] = infinity;
    nTimeSteps = 20;

    /* Paper 2016
    tstep[0] = 0.2;
    tstep[1] = 0.4;
    tstep[2] = 0.6;
    tstep[3] = 0.8;
    tstep[4] = 1.0;
    tstep[5] = 2.0;
    tstep[6] = 3.0;
    tstep[7] = 4.0;
    tstep[8] = 5.0;
    tstep[9] = 6.0;
    tstep[10] = 7.0;
    tstep[11] = 8.0;
    tstep[12] = 9.0;
    tstep[13] = 10.0;
    tstep[14] = infinity;
    nTimeSteps = 15;*/
  } 
  else if (studyDileptons)
  {
    tstep[0] = 0.2;
    tstep[1] = 0.4;
    tstep[2] = 0.6;
    tstep[3] = 0.8;
    tstep[4] = 1.0;
    tstep[5] = 2.0;
    tstep[6] = 3.0;
    tstep[7] = 4.0;
    tstep[8] = 5.0;
    tstep[9] = 6.0;
    tstep[10] = 7.0;
    tstep[11] = 8.0;
    tstep[12] = 9.0;
    tstep[13] = 10.0;
    tstep[14] = infinity;
    nTimeSteps = 15;
  } 
  else if (studyEtSpectra)
  {
      tstep[0] =1.;      //fm/c
      tstep[1] = 2.;      //fm/c
      tstep[2] = 3.;      //fm/c
      tstep[3] = 4.;      //fm/c
      tstep[4] = 5.;      //fm/c
      tstep[5] = 6.;      //fm/c
      tstep[6] = 7.;      //fm/c
      tstep[7] = infinity; //fm/cn
      nTimeSteps = 8;
  }
  else if (studyNumberOfBackgroundQuarks)
  {
    tstep[0]=.15;
    tstep[1]=.2;
    tstep[2]=.3;
    tstep[3]=.4;
    tstep[4]=.5;
    tstep[5]=.6;
    tstep[6]=.7;
    tstep[7]=.8;
    tstep[8]=.9;
    tstep[9]=1.0;
    tstep[10]=1.1;
    tstep[11]=1.2;
    tstep[12]=1.3;
    tstep[13]=1.4;
    tstep[14]=1.5;
    tstep[15]=1.6;
    tstep[16]=1.7;
    tstep[17]=1.8;
    tstep[18]=1.9;
    tstep[19]=2.0;
    tstep[20]=2.1;
    tstep[21]=2.2;
    tstep[22]=2.3;
    tstep[23]=2.4;
    tstep[24]=2.5;
    tstep[25]=2.6;
    tstep[26]=2.7;
    tstep[27]=2.8;
    tstep[28]=2.9;
    tstep[29]=3.0;
    tstep[30]=3.1;
    tstep[31]=3.2;
    tstep[32]=3.3;
    tstep[33]=3.4;
    tstep[34]=3.5;
    tstep[35]=3.6;
    tstep[36]=3.7;
    tstep[37]=3.8;
    tstep[38]=3.9;
    tstep[39]=4.0;
    tstep[40]=4.1;
    tstep[41]=4.2;
    tstep[42]=4.3;
    tstep[43]=4.4;
    tstep[44]=4.5;
    tstep[45]=4.6;
    tstep[46]=4.7;
    tstep[47]=4.8;
    tstep[48]=4.9;
    tstep[49]=5.0;
    tstep[50]=5.1;
    tstep[51]=5.2;
    tstep[52]=5.3;
    tstep[53]=5.4;
    tstep[54]=5.5;
    tstep[55]=5.6;
    tstep[56]=5.7;
    tstep[57]=5.8;
    tstep[58]=5.9;
    tstep[59]=6.0;
    tstep[60]=6.1;
    tstep[61]=6.2;
    tstep[62]=6.3;
    tstep[63]=6.4;
    tstep[64]=6.5;
    tstep[65]=6.6;
    tstep[66]=6.7;
    tstep[67]=6.8;
    tstep[68]=6.9;
    tstep[69]=7.0;
    tstep[70]=7.1;
    tstep[71]=7.2;
    tstep[72]=7.3;
    tstep[73]=7.4;
    tstep[74]=7.5;
    tstep[75]=7.6;
    tstep[76]=7.7;
    tstep[77]=7.8;
    tstep[78]=7.9;
    tstep[79]=8.0;
    tstep[80]=8.1;
    tstep[81]=infinity;
    nTimeSteps = 82;
  }
  else if(studySpatialPhotons) 
  {
    tstep[0]=.15;
    tstep[1]=.2;
    tstep[2]=.3;
    tstep[3]=.4;
    tstep[4]=.5;
    tstep[5]=.6;
    tstep[6]=.7;
    tstep[7]=.8;
    tstep[8]=.9;
    tstep[9]=1.0;
    tstep[10]=1.1;
    tstep[11]=1.2;
    tstep[12]=1.3;
    tstep[13]=1.4;
    tstep[14]=1.5;
    tstep[15]=1.6;
    tstep[16]=1.7;
    tstep[17]=1.8;
    tstep[18]=1.9;
    tstep[19]=2.0;
    tstep[20]=2.1;
    tstep[21]=2.2;
    tstep[22]=2.3;
    tstep[23]=2.4;
    tstep[24]=2.5;
    tstep[25]=2.6;
    tstep[26]=2.7;
    tstep[27]=2.8;
    tstep[28]=2.9;
    tstep[29]=3.0;
    tstep[30]=3.1;
    tstep[31]=3.2;
    tstep[32]=3.3;
    tstep[33]=3.4;
    tstep[34]=3.5;
    tstep[35]=3.6;
    tstep[36]=3.7;
    tstep[37]=3.8;
    tstep[38]=3.9;
    tstep[39]=4.0;
    tstep[40]=4.1;
    tstep[41]=4.2;
    tstep[42]=4.3;
    tstep[43]=4.4;
    tstep[44]=4.5;
    tstep[45]=4.6;
    tstep[46]=4.7;
    tstep[47]=4.8;
    tstep[48]=4.9;
    tstep[49]=5.0;
    tstep[50]=5.1;
    tstep[51]=5.2;
    tstep[52]=5.3;
    tstep[53]=5.4;
    tstep[54]=5.5;
    tstep[55]=5.6;
    tstep[56]=5.7;
    tstep[57]=5.8;
    tstep[58]=5.9;
    tstep[59]=6.0;
    tstep[60]=6.1;
    tstep[61]=6.2;
    tstep[62]=6.3;
    tstep[63]=6.4;
    tstep[64]=6.5;
    tstep[65]=6.6;
    tstep[66]=6.7;
    tstep[67]=6.8;
    tstep[68]=6.9;
    tstep[69]=7.0;
    tstep[70]=7.1;
    tstep[71]=7.2;
    tstep[72]=7.3;
    tstep[73]=7.4;
    tstep[74]=7.5;
    tstep[75]=7.6;
    tstep[76]=7.7;
    tstep[77]=7.8;
    tstep[78]=7.9;
    tstep[79]=8.0;
    tstep[80]=8.1;
    tstep[81]=infinity;
    nTimeSteps = 82;
  }
  else if(studyTempCustom) 
  {
    tstep[0]=.15;
    tstep[1]=.2;
    tstep[2]=.3;
    tstep[3]=.4;
    tstep[4]=.5;
    tstep[5]=.6;
    tstep[6]=.7;
    tstep[7]=.8;
    tstep[8]=.9;
    tstep[9]=1.0;
    tstep[10]=1.1;
    tstep[11]=1.2;
    tstep[12]=1.3;
    tstep[13]=1.4;
    tstep[14]=1.5;
    tstep[15]=1.6;
    tstep[16]=1.7;
    tstep[17]=1.8;
    tstep[18]=1.9;
    tstep[19]=2.0;
    tstep[20]=2.1;
    tstep[21]=2.2;
    tstep[22]=2.3;
    tstep[23]=2.4;
    tstep[24]=2.5;
    tstep[25]=2.6;
    tstep[26]=2.7;
    tstep[27]=2.8;
    tstep[28]=2.9;
    tstep[29]=3.0;
    tstep[30]=3.1;
    tstep[31]=3.2;
    tstep[32]=3.3;
    tstep[33]=3.4;
    tstep[34]=3.5;
    tstep[35]=3.6;
    tstep[36]=3.7;
    tstep[37]=3.8;
    tstep[38]=3.9;
    tstep[39]=4.0;
    tstep[40]=4.1;
    tstep[41]=4.2;
    tstep[42]=4.3;
    tstep[43]=4.4;
    tstep[44]=4.5;
    tstep[45]=4.6;
    tstep[46]=4.7;
    tstep[47]=4.8;
    tstep[48]=4.9;
    tstep[49]=5.0;
    tstep[50]=5.1;
    tstep[51]=5.2;
    tstep[52]=5.3;
    tstep[53]=5.4;
    tstep[54]=5.5;
    tstep[55]=5.6;
    tstep[56]=5.7;
    tstep[57]=5.8;
    tstep[58]=5.9;
    tstep[59]=6.0;
    tstep[60]=6.1;
    tstep[61]=6.2;
    tstep[62]=6.3;
    tstep[63]=6.4;
    tstep[64]=6.5;
    tstep[65]=6.6;
    tstep[66]=6.7;
    tstep[67]=6.8;
    tstep[68]=6.9;
    tstep[69]=7.0;
    tstep[70]=7.1;
    tstep[71]=7.2;
    tstep[72]=7.3;
    tstep[73]=7.4;
    tstep[74]=7.5;
    tstep[75]=7.6;
    tstep[76]=7.7;
    tstep[77]=7.8;
    tstep[78]=7.9;
    tstep[79]=8.0;
    tstep[80]=8.1;
    tstep[81]=infinity;
    nTimeSteps = 82;
  }
  else if(studyThermalisation||studyPartons) 
  { 
    tstep[0]=.15;
    tstep[1]=.2;
    tstep[2]=.3;
    tstep[3]=.4;
    tstep[4]=.5;
    tstep[5]=.6;
    tstep[6]=.7;
    tstep[7]=.8;
    tstep[8]=.9;
    tstep[9]=1.0;
    tstep[10]=1.1;
    tstep[11]=1.2;
    tstep[12]=1.3;
    tstep[13]=1.4;
    tstep[14]=1.5;
    tstep[15]=1.6;
    tstep[16]=1.7;
    tstep[17]=1.8;
    tstep[18]=1.9;
    tstep[19]=2.0;
    tstep[20]=2.1;
    tstep[21]=2.2;
    tstep[22]=2.3;
    tstep[23]=2.4;
    tstep[24]=2.5;
    tstep[25]=2.6;
    tstep[26]=2.7;
    tstep[27]=2.8;
    tstep[28]=2.9;
    tstep[29]=3.0;
    tstep[30]=3.1;
    tstep[31]=3.2;
    tstep[32]=3.3;
    tstep[33]=3.4;
    tstep[34]=3.5;
    tstep[35]=3.6;
    tstep[36]=3.7;
    tstep[37]=3.8;
    tstep[38]=3.9;
    tstep[39]=4.0;
    tstep[40]=4.1;
    tstep[41]=4.2;
    tstep[42]=4.3;
    tstep[43]=4.4;
    tstep[44]=4.5;
    tstep[45]=4.6;
    tstep[46]=4.7;
    tstep[47]=4.8;
    tstep[48]=4.9;
    tstep[49]=5.0;
    tstep[50]=5.1;
    tstep[51]=5.2;
    tstep[52]=5.3;
    tstep[53]=5.4;
    tstep[54]=5.5;
    tstep[55]=5.6;
    tstep[56]=5.7;
    tstep[57]=5.8;
    tstep[58]=5.9;
    tstep[59]=6.0;
    tstep[60]=6.1;
    tstep[61]=6.2;
    tstep[62]=6.3;
    tstep[63]=6.4;
    tstep[64]=6.5;
    tstep[65]=6.6;
    tstep[66]=6.7;
    tstep[67]=6.8;
    tstep[68]=6.9;
    tstep[69]=7.0;
    tstep[70]=7.1;
    tstep[71]=7.2;
    tstep[72]=7.3;
    tstep[73]=7.4;
    tstep[74]=7.5;
    tstep[75]=7.6;
    tstep[76]=7.7;
    tstep[77]=7.8;
    tstep[78]=7.9;
    tstep[79]=8.0;
    tstep[80]=8.1;
    tstep[81]=infinity;
    nTimeSteps = 82;    
    /*tstep[0]=.15;
    tstep[1]=.16;
    tstep[2]=.17;
    tstep[3]=.18;
    tstep[4]=.19;
    tstep[5]=.20;
    tstep[6]=.21;
    tstep[7]=.22;
    tstep[8]=.23;
    tstep[9]=.24;
    tstep[10]=.25;
    tstep[11]=.26;
    tstep[12]=.27;
    tstep[13]=.28;
    tstep[14]=.29;
    tstep[15]=.30;
    tstep[16]=.31;
    tstep[17]=.32;
    tstep[18]=.33;
    tstep[19]=.34;
    tstep[20]=.35;
    tstep[21]=.36;
    tstep[22]=.37;
    tstep[23]=.38;
    tstep[24]=.39;
    tstep[25]=.40;
    tstep[26]=.41;
    tstep[27]=.42;
    tstep[28]=.43;
    tstep[29]=.44;
    tstep[30]=.45;
    tstep[31]=.46;
    tstep[32]=.47;
    tstep[33]=.48;
    tstep[34]=.49;
    tstep[35]=.50;
    tstep[36]=.51;
    tstep[37]=.52;
    tstep[38]=.53;
    tstep[39]=.54;
    tstep[40]=.55;
    tstep[41]=.56;
    tstep[42]=.57;
    tstep[43]=.58;
    tstep[44]=.59;
    tstep[45]=.60;
    tstep[46]=.61;
    tstep[47]=.62;
    tstep[48]=.63;
    tstep[49]=.64;
    tstep[50]=.65;
    tstep[51]=.66;
    tstep[52]=.67;
    tstep[53]=.68;
    tstep[54]=.69;
    tstep[55]=.70;
    tstep[56]=.71;
    tstep[57]=.72;
    tstep[58]=.73;
    tstep[59]=.74;
    tstep[60]=.75;
    tstep[61]=.76;
    tstep[62]=.77;
    tstep[63]=.78;
    tstep[64]=.79;
    tstep[65]=.80;
    tstep[66]=.81;
    tstep[67]=.82;
    tstep[68]=.83;
    tstep[69]=.84;
    tstep[70]=.85;
    tstep[71]=.86;
    tstep[72]=.87;
    tstep[73]=.88;
    tstep[74]=.89;
    tstep[75]=.90;
    tstep[76]=.91;
    tstep[77]=.92;
    tstep[78]=.93;
    tstep[79]=.94;
    tstep[80]=.95;
    tstep[81]=.96;
    tstep[82]=.97;
    tstep[83]=.98;
    tstep[84]=.99;
    tstep[85]=1.00;
    tstep[86]=1.01;
    tstep[87]=1.02;
    tstep[88]=1.03;
    tstep[89]=1.04;
    tstep[90]=1.05;
    tstep[91]=1.06;
    tstep[92]=1.07;
    tstep[93]=1.08;
    tstep[94]=1.09;
    tstep[95]=1.10;
    tstep[96]=1.11;
    tstep[97]=1.12;
    tstep[98]=1.13;
    tstep[99]=1.14;
    tstep[100]=1.15;
    tstep[101]=1.16;
    tstep[102]=1.17;
    tstep[103]=1.18;
    tstep[104]=1.19;
    tstep[105]=1.20;
    tstep[106]=1.21;
    tstep[107]=1.22;
    tstep[108]=1.23;
    tstep[109]=1.24;
    tstep[110]=1.25;
    tstep[111]=1.26;
    tstep[112]=1.27;
    tstep[113]=1.28;
    tstep[114]=1.29;
    tstep[115]=1.30;
    tstep[116]=1.31;
    tstep[117]=1.32;
    tstep[118]=1.33;
    tstep[119]=1.34;
    tstep[120]=1.35;
    tstep[121]=1.36;
    tstep[122]=1.37;
    tstep[123]=1.38;
    tstep[124]=1.39;
    tstep[125]=1.40;
    tstep[126]=1.41;
    tstep[127]=1.42;
    tstep[128]=1.43;
    tstep[129]=1.44;
    tstep[130]=1.45;
    tstep[131]=1.46;
    tstep[132]=1.47;
    tstep[133]=1.48;
    tstep[134]=1.49;
    tstep[135]=1.50;
    tstep[136]=1.51;
    tstep[137]=1.52;
    tstep[138]=1.53;
    tstep[139]=1.54;
    tstep[140]=1.55;
    tstep[141]=1.56;
    tstep[142]=1.57;
    tstep[143]=1.58;
    tstep[144]=1.59;
    tstep[145]=1.60;
    tstep[146]=1.61;
    tstep[147]=1.62;
    tstep[148]=1.63;
    tstep[149]=1.64;
    tstep[150]=1.65;
    tstep[151]=1.66;
    tstep[152]=1.67;
    tstep[153]=1.68;
    tstep[154]=1.69;
    tstep[155]=1.70;
    tstep[156]=1.71;
    tstep[157]=1.72;
    tstep[158]=1.73;
    tstep[159]=1.74;
    tstep[160]=1.75;
    tstep[161]=1.76;
    tstep[162]=1.77;
    tstep[163]=1.78;
    tstep[164]=1.79;
    tstep[165]=1.80;
    tstep[166]=1.81;
    tstep[167]=1.82;
    tstep[168]=1.83;
    tstep[169]=1.84;
    tstep[170]=1.85;
    tstep[171]=1.86;
    tstep[172]=1.87;
    tstep[173]=1.88;
    tstep[174]=1.89;
    tstep[175]=1.90;
    tstep[176]=1.91;
    tstep[177]=1.92;
    tstep[178]=1.93;
    tstep[179]=1.94;
    tstep[180]=1.95;
    tstep[181]=1.96;
    tstep[182]=1.97;
    tstep[183]=1.98;
    tstep[184]=1.99;
    tstep[185]=2.00;
    tstep[186]=2.01;
    tstep[187]=2.02;
    tstep[188]=2.03;
    tstep[189]=2.04;
    tstep[190]=2.05;
    tstep[191]=2.06;
    tstep[192]=2.07;
    tstep[193]=2.08;
    tstep[194]=2.09;
    tstep[195]=2.10;
    tstep[196]=2.11;
    tstep[197]=2.12;
    tstep[198]=2.13;
    tstep[199]=2.14;
    tstep[200]=2.15;
    tstep[201]=2.16;
    tstep[202]=2.17;
    tstep[203]=2.18;
    tstep[204]=2.19;
    tstep[205]=2.20;
    tstep[206]=2.21;
    tstep[207]=2.22;
    tstep[208]=2.23;
    tstep[209]=2.24;
    tstep[210]=2.25;
    tstep[211]=2.26;
    tstep[212]=2.27;
    tstep[213]=2.28;
    tstep[214]=2.29;
    tstep[215]=2.30;
    tstep[216]=2.31;
    tstep[217]=2.32;
    tstep[218]=2.33;
    tstep[219]=2.34;
    tstep[220]=2.35;
    tstep[221]=2.36;
    tstep[222]=2.37;
    tstep[223]=2.38;
    tstep[224]=2.39;
    tstep[225]=2.40;
    tstep[226]=2.41;
    tstep[227]=2.42;
    tstep[228]=2.43;
    tstep[229]=2.44;
    tstep[230]=2.45;
    tstep[231]=2.46;
    tstep[232]=2.47;
    tstep[233]=2.48;
    tstep[234]=2.49;
    tstep[235]=2.50;
    tstep[236]=2.51;
    tstep[237]=2.52;
    tstep[238]=2.53;
    tstep[239]=2.54;
    tstep[240]=2.55;
    tstep[241]=2.56;
    tstep[242]=2.57;
    tstep[243]=2.58;
    tstep[244]=2.59;
    tstep[245]=2.60;
    tstep[246]=2.61;
    tstep[247]=2.62;
    tstep[248]=2.63;
    tstep[249]=2.64;
    tstep[250]=2.65;
    tstep[251]=2.66;
    tstep[252]=2.67;
    tstep[253]=2.68;
    tstep[254]=2.69;
    tstep[255]=2.70;
    tstep[256]=2.71;
    tstep[257]=2.72;
    tstep[258]=2.73;
    tstep[259]=2.74;
    tstep[260]=2.75;
    tstep[261]=2.76;
    tstep[262]=2.77;
    tstep[263]=2.78;
    tstep[264]=2.79;
    tstep[265]=2.80;
    tstep[266]=2.81;
    tstep[267]=2.82;
    tstep[268]=2.83;
    tstep[269]=2.84;
    tstep[270]=2.85;
    tstep[271]=2.86;
    tstep[272]=2.87;
    tstep[273]=2.88;
    tstep[274]=2.89;
    tstep[275]=2.90;
    tstep[276]=2.91;
    tstep[277]=2.92;
    tstep[278]=2.93;
    tstep[279]=2.94;
    tstep[280]=2.95;
    tstep[281]=2.96;
    tstep[282]=2.97;
    tstep[283]=2.98;
    tstep[284]=2.99;
    tstep[285]=3.00;
    tstep[286]=3.01;
    tstep[287]=3.02;
    tstep[288]=3.03;
    tstep[289]=3.04;
    tstep[290]=3.05;
    tstep[291]=3.06;
    tstep[292]=3.07;
    tstep[293]=3.08;
    tstep[294]=3.09;
    tstep[295]=3.10;
    tstep[296]=3.11;
    tstep[297]=3.12;
    tstep[298]=3.13;
    tstep[299]=3.14;
    tstep[300]=3.15;
    tstep[301]=3.16;
    tstep[302]=3.17;
    tstep[303]=3.18;
    tstep[304]=3.19;
    tstep[305]=3.20;
    tstep[306]=3.21;
    tstep[307]=3.22;
    tstep[308]=3.23;
    tstep[309]=3.24;
    tstep[310]=3.25;
    tstep[311]=3.26;
    tstep[312]=3.27;
    tstep[313]=3.28;
    tstep[314]=3.29;
    tstep[315]=3.30;
    tstep[316]=3.31;
    tstep[317]=3.32;
    tstep[318]=3.33;
    tstep[319]=3.34;
    tstep[320]=3.35;
    tstep[321]=3.36;
    tstep[322]=3.37;
    tstep[323]=3.38;
    tstep[324]=3.39;
    tstep[325]=3.40;
    tstep[326]=3.41;
    tstep[327]=3.42;
    tstep[328]=3.43;
    tstep[329]=3.44;
    tstep[330]=3.45;
    tstep[331]=3.46;
    tstep[332]=3.47;
    tstep[333]=3.48;
    tstep[334]=3.49;
    tstep[335]=3.50;
    tstep[336]=3.51;
    tstep[337]=3.52;
    tstep[338]=3.53;
    tstep[339]=3.54;
    tstep[340]=3.55;
    tstep[341]=3.56;
    tstep[342]=3.57;
    tstep[343]=3.58;
    tstep[344]=3.59;
    tstep[345]=3.60;
    tstep[346]=3.61;
    tstep[347]=3.62;
    tstep[348]=3.63;
    tstep[349]=3.64;
    tstep[350]=3.65;
    tstep[351]=3.66;
    tstep[352]=3.67;
    tstep[353]=3.68;
    tstep[354]=3.69;
    tstep[355]=3.70;
    tstep[356]=3.71;
    tstep[357]=3.72;
    tstep[358]=3.73;
    tstep[359]=3.74;
    tstep[360]=3.75;
    tstep[361]=3.76;
    tstep[362]=3.77;
    tstep[363]=3.78;
    tstep[364]=3.79;
    tstep[365]=3.80;
    tstep[366]=3.81;
    tstep[367]=3.82;
    tstep[368]=3.83;
    tstep[369]=3.84;
    tstep[370]=3.85;
    tstep[371]=3.86;
    tstep[372]=3.87;
    tstep[373]=3.88;
    tstep[374]=3.89;
    tstep[375]=3.90;
    tstep[376]=3.91;
    tstep[377]=3.92;
    tstep[378]=3.93;
    tstep[379]=3.94;
    tstep[380]=3.95;
    tstep[381]=3.96;
    tstep[382]=3.97;
    tstep[383]=3.98;
    tstep[384]=3.99;
    tstep[385]=4.00;
    tstep[386]=4.01;
    tstep[387]=4.02;
    tstep[388]=4.03;
    tstep[389]=4.04;
    tstep[390]=4.05;
    tstep[391]=4.06;
    tstep[392]=4.07;
    tstep[393]=4.08;
    tstep[394]=4.09;
    tstep[395]=4.10;
    tstep[396]=4.11;
    tstep[397]=4.12;
    tstep[398]=4.13;
    tstep[399]=4.14;
    tstep[400]=4.15;
    tstep[401]=4.16;
    tstep[402]=4.17;
    tstep[403]=4.18;
    tstep[404]=4.19;
    tstep[405]=4.20;
    tstep[406]=4.21;
    tstep[407]=4.22;
    tstep[408]=4.23;
    tstep[409]=4.24;
    tstep[410]=4.25;
    tstep[411]=4.26;
    tstep[412]=4.27;
    tstep[413]=4.28;
    tstep[414]=4.29;
    tstep[415]=4.30;
    tstep[416]=4.31;
    tstep[417]=4.32;
    tstep[418]=4.33;
    tstep[419]=4.34;
    tstep[420]=4.35;
    tstep[421]=4.36;
    tstep[422]=4.37;
    tstep[423]=4.38;
    tstep[424]=4.39;
    tstep[425]=4.40;
    tstep[426]=4.41;
    tstep[427]=4.42;
    tstep[428]=4.43;
    tstep[429]=4.44;
    tstep[430]=4.45;
    tstep[431]=4.46;
    tstep[432]=4.47;
    tstep[433]=4.48;
    tstep[434]=4.49;
    tstep[435]=4.50;
    tstep[436]=4.51;
    tstep[437]=4.52;
    tstep[438]=4.53;
    tstep[439]=4.54;
    tstep[440]=4.55;
    tstep[441]=4.56;
    tstep[442]=4.57;
    tstep[443]=4.58;
    tstep[444]=4.59;
    tstep[445]=4.60;
    tstep[446]=4.61;
    tstep[447]=4.62;
    tstep[448]=4.63;
    tstep[449]=4.64;
    tstep[450]=4.65;
    tstep[451]=4.66;
    tstep[452]=4.67;
    tstep[453]=4.68;
    tstep[454]=4.69;
    tstep[455]=4.70;
    tstep[456]=4.71;
    tstep[457]=4.72;
    tstep[458]=4.73;
    tstep[459]=4.74;
    tstep[460]=4.75;
    tstep[461]=4.76;
    tstep[462]=4.77;
    tstep[463]=4.78;
    tstep[464]=4.79;
    tstep[465]=4.80;
    tstep[466]=4.81;
    tstep[467]=4.82;
    tstep[468]=4.83;
    tstep[469]=4.84;
    tstep[470]=4.85;
    tstep[471]=4.86;
    tstep[472]=4.87;
    tstep[473]=4.88;
    tstep[474]=4.89;
    tstep[475]=4.90;
    tstep[476]=4.91;
    tstep[477]=4.92;
    tstep[478]=4.93;
    tstep[479]=4.94;
    tstep[480]=4.95;
    tstep[481]=4.96;
    tstep[482]=4.97;
    tstep[483]=4.98;
    tstep[484]=4.99;
    tstep[485]=5.00;
    tstep[486]=5.01;
    tstep[487]=5.02;
    tstep[488]=5.03;
    tstep[489]=5.04;
    tstep[490]=5.05;
    tstep[491]=5.06;
    tstep[492]=5.07;
    tstep[493]=5.08;
    tstep[494]=5.09;
    tstep[495]=5.10;
    tstep[496]=5.11;
    tstep[497]=5.12;
    tstep[498]=5.13;
    tstep[499]=5.14;
    tstep[500]=5.15;
    tstep[501]=5.16;
    tstep[502]=5.17;
    tstep[503]=5.18;
    tstep[504]=5.19;
    tstep[505]=5.20;
    tstep[506]=5.21;
    tstep[507]=5.22;
    tstep[508]=5.23;
    tstep[509]=5.24;
    tstep[510]=5.25;
    tstep[511]=5.26;
    tstep[512]=5.27;
    tstep[513]=5.28;
    tstep[514]=5.29;
    tstep[515]=5.30;
    tstep[516]=5.31;
    tstep[517]=5.32;
    tstep[518]=5.33;
    tstep[519]=5.34;
    tstep[520]=5.35;
    tstep[521]=5.36;
    tstep[522]=5.37;
    tstep[523]=5.38;
    tstep[524]=5.39;
    tstep[525]=5.40;
    tstep[526]=5.41;
    tstep[527]=5.42;
    tstep[528]=5.43;
    tstep[529]=5.44;
    tstep[530]=5.45;
    tstep[531]=5.46;
    tstep[532]=5.47;
    tstep[533]=5.48;
    tstep[534]=5.49;
    tstep[535]=5.50;
    tstep[536]=5.51;
    tstep[537]=5.52;
    tstep[538]=5.53;
    tstep[539]=5.54;
    tstep[540]=5.55;
    tstep[541]=5.56;
    tstep[542]=5.57;
    tstep[543]=5.58;
    tstep[544]=5.59;
    tstep[545]=5.60;
    tstep[546]=5.61;
    tstep[547]=5.62;
    tstep[548]=5.63;
    tstep[549]=5.64;
    tstep[550]=5.65;
    tstep[551]=5.66;
    tstep[552]=5.67;
    tstep[553]=5.68;
    tstep[554]=5.69;
    tstep[555]=5.70;
    tstep[556]=5.71;
    tstep[557]=5.72;
    tstep[558]=5.73;
    tstep[559]=5.74;
    tstep[560]=5.75;
    tstep[561]=5.76;
    tstep[562]=5.77;
    tstep[563]=5.78;
    tstep[564]=5.79;
    tstep[565]=5.80;
    tstep[566]=5.81;
    tstep[567]=5.82;
    tstep[568]=5.83;
    tstep[569]=5.84;
    tstep[570]=5.85;
    tstep[571]=5.86;
    tstep[572]=5.87;
    tstep[573]=5.88;
    tstep[574]=5.89;
    tstep[575]=5.90;
    tstep[576]=5.91;
    tstep[577]=5.92;
    tstep[578]=5.93;
    tstep[579]=5.94;
    tstep[580]=5.95;
    tstep[581]=5.96;
    tstep[582]=5.97;
    tstep[583]=5.98;
    tstep[584]=5.99;
    tstep[585]=6.00;
    tstep[586]=6.01;
    tstep[587]=6.02;
    tstep[588]=6.03;
    tstep[589]=6.04;
    tstep[590]=6.05;
    tstep[591]=6.06;
    tstep[592]=6.07;
    tstep[593]=6.08;
    tstep[594]=6.09;
    tstep[595]=6.10;
    tstep[596]=6.11;
    tstep[597]=6.12;
    tstep[598]=6.13;
    tstep[599]=6.14;
    tstep[600]=6.15;
    tstep[601]=6.16;
    tstep[602]=6.17;
    tstep[603]=6.18;
    tstep[604]=6.19;
    tstep[605]=6.20;
    tstep[606]=6.21;
    tstep[607]=6.22;
    tstep[608]=6.23;
    tstep[609]=6.24;
    tstep[610]=6.25;
    tstep[611]=6.26;
    tstep[612]=6.27;
    tstep[613]=6.28;
    tstep[614]=6.29;
    tstep[615]=6.30;
    tstep[616]=6.31;
    tstep[617]=6.32;
    tstep[618]=6.33;
    tstep[619]=6.34;
    tstep[620]=6.35;
    tstep[621]=6.36;
    tstep[622]=6.37;
    tstep[623]=6.38;
    tstep[624]=6.39;
    tstep[625]=6.40;
    tstep[626]=6.41;
    tstep[627]=6.42;
    tstep[628]=6.43;
    tstep[629]=6.44;
    tstep[630]=6.45;
    tstep[631]=6.46;
    tstep[632]=6.47;
    tstep[633]=6.48;
    tstep[634]=6.49;
    tstep[635]=6.50;
    tstep[636]=6.51;
    tstep[637]=6.52;
    tstep[638]=6.53;
    tstep[639]=6.54;
    tstep[640]=6.55;
    tstep[641]=6.56;
    tstep[642]=6.57;
    tstep[643]=6.58;
    tstep[644]=6.59;
    tstep[645]=6.60;
    tstep[646]=6.61;
    tstep[647]=6.62;
    tstep[648]=6.63;
    tstep[649]=6.64;
    tstep[650]=6.65;
    tstep[651]=6.66;
    tstep[652]=6.67;
    tstep[653]=6.68;
    tstep[654]=6.69;
    tstep[655]=6.70;
    tstep[656]=6.71;
    tstep[657]=6.72;
    tstep[658]=6.73;
    tstep[659]=6.74;
    tstep[660]=6.75;
    tstep[661]=6.76;
    tstep[662]=6.77;
    tstep[663]=6.78;
    tstep[664]=6.79;
    tstep[665]=6.80;
    tstep[666]=6.81;
    tstep[667]=6.82;
    tstep[668]=6.83;
    tstep[669]=6.84;
    tstep[670]=6.85;
    tstep[671]=6.86;
    tstep[672]=6.87;
    tstep[673]=6.88;
    tstep[674]=6.89;
    tstep[675]=6.90;
    tstep[676]=6.91;
    tstep[677]=6.92;
    tstep[678]=6.93;
    tstep[679]=6.94;
    tstep[680]=6.95;
    tstep[681]=6.96;
    tstep[682]=6.97;
    tstep[683]=6.98;
    tstep[684]=6.99;
    tstep[685]=7.00;
    tstep[686]=7.01;
    tstep[687]=7.02;
    tstep[688]=7.03;
    tstep[689]=7.04;
    tstep[690]=7.05;
    tstep[691]=7.06;
    tstep[692]=7.07;
    tstep[693]=7.08;
    tstep[694]=7.09;
    tstep[695]=7.10;
    tstep[696]=7.11;
    tstep[697]=7.12;
    tstep[698]=7.13;
    tstep[699]=7.14;
    tstep[700]=7.15;
    tstep[701]=7.16;
    tstep[702]=7.17;
    tstep[703]=7.18;
    tstep[704]=7.19;
    tstep[705]=7.20;
    tstep[706]=7.21;
    tstep[707]=7.22;
    tstep[708]=7.23;
    tstep[709]=7.24;
    tstep[710]=7.25;
    tstep[711]=7.26;
    tstep[712]=7.27;
    tstep[713]=7.28;
    tstep[714]=7.29;
    tstep[715]=7.30;
    tstep[716]=7.31;
    tstep[717]=7.32;
    tstep[718]=7.33;
    tstep[719]=7.34;
    tstep[720]=7.35;
    tstep[721]=7.36;
    tstep[722]=7.37;
    tstep[723]=7.38;
    tstep[724]=7.39;
    tstep[725]=7.40;
    tstep[726]=7.41;
    tstep[727]=7.42;
    tstep[728]=7.43;
    tstep[729]=7.44;
    tstep[730]=7.45;
    tstep[731]=7.46;
    tstep[732]=7.47;
    tstep[733]=7.48;
    tstep[734]=7.49;
    tstep[735]=7.50;
    tstep[736]=7.51;
    tstep[737]=7.52;
    tstep[738]=7.53;
    tstep[739]=7.54;
    tstep[740]=7.55;
    tstep[741]=7.56;
    tstep[742]=7.57;
    tstep[743]=7.58;
    tstep[744]=7.59;
    tstep[745]=7.60;
    tstep[746]=7.61;
    tstep[747]=7.62;
    tstep[748]=7.63;
    tstep[749]=7.64;
    tstep[750]=7.65;
    tstep[751]=7.66;
    tstep[752]=7.67;
    tstep[753]=7.68;
    tstep[754]=7.69;
    tstep[755]=7.70;
    tstep[756]=7.71;
    tstep[757]=7.72;
    tstep[758]=7.73;
    tstep[759]=7.74;
    tstep[760]=7.75;
    tstep[761]=7.76;
    tstep[762]=7.77;
    tstep[763]=7.78;
    tstep[764]=7.79;
    tstep[765]=7.80;
    tstep[766]=7.81;
    tstep[767]=7.82;
    tstep[768]=7.83;
    tstep[769]=7.84;
    tstep[770]=7.85;
    tstep[771]=7.86;
    tstep[772]=7.87;
    tstep[773]=7.88;
    tstep[774]=7.89;
    tstep[775]=7.90;
    tstep[776]=7.91;
    tstep[777]=7.92;
    tstep[778]=7.93;
    tstep[779]=7.94;
    tstep[780]=7.95;
    tstep[781]=7.96;
    tstep[782]=7.97;
    tstep[783]=7.98;
    tstep[784]=7.99;
    tstep[785]=8.00;
    tstep[786]=8.01;
    tstep[787]=8.02;
    tstep[788]=8.03;
    tstep[789]=8.04;
    tstep[790]=8.05;
    tstep[791]=8.06;
    tstep[792]=8.07;
    tstep[793]=8.08;
    tstep[794]=8.09;
    tstep[795]=8.10;
    tstep[796]=8.11;
    tstep[797]=8.12;
    tstep[798]=8.13;
    tstep[799]=8.14;
    tstep[800]=8.15;
    tstep[801]=infinity;
    nTimeSteps = 802;*/
  }
  else
  {
    tstep[0] = 0.1;      //fm/c
    tstep[1] = 0.5;      //fm/c
    tstep[2] = 1.0;      //fm/c
    tstep[3] = 1.5;      //fm/c
    tstep[4] = 2.0;      //fm/c
    tstep[5] = 2.5;      //fm/c
    tstep[6] = 3.0;      //fm/c
    tstep[7] = 3.5;      //fm/c
    tstep[8] = 4.0;      //fm/c
    tstep[9] = 4.5;      //fm/c
    tstep[10] = 5.0;      //fm/c
    tstep[11] = 5.5;      //fm/c
    tstep[12] = 6.0;      //fm/c
    tstep[13] = 6.5;      //fm/c
    tstep[14] = 7.0;      //fm/c
    tstep[15] = 7.5;      //fm/c
    tstep[16] = 8.0;      //fm/c
    tstep[17] = 9.0;      //fm/c
    tstep[18] = 10.0;
    tstep[19] = infinity; //fm/c
    nTimeSteps = 20;
  }
  //--------------------------------------
  
    
  if (studyNumberOfBackgroundQuarks)
  {
    NumberOfQuarks.resize(nTimeSteps+1);
    NumberOfAntiquarks.resize(nTimeSteps+1);
    NumberOfGluons.resize(nTimeSteps+1);
  
    for (int i=0;i<nTimeSteps-1;i++)
    {
      NumberOfQuarks[i]=0;
      NumberOfAntiquarks[i]=0;
      NumberOfGluons[i]=0;
    }
  }
  
  //---- times for output of data --------
  int tempCount_tstep = 0;
  tstep_movie[tempCount_tstep] = 0.1;
  do
  {
    ++tempCount_tstep;
    tstep_movie[tempCount_tstep] = tstep_movie[tempCount_tstep - 1] + 0.1;
  }
  while( tstep_movie[tempCount_tstep] <= 10.0 );
  
  ++tempCount_tstep;
  tstep_movie[tempCount_tstep] = infinity;
  
  nTimeSteps_movie = tempCount_tstep + 1;
  //--------------------------------------
 
  
  //--------------------------------------
  string name_fug, name_temp, mfpName, centralDensitiesName, oscarName_jet, oscarName_background;
  if ( studyJpsi )
    name_fug = filename_prefix + "_jpsi_fugacity";
  else
    name_fug = "/dev/null";
  
  
  
  if( studyTempInTube )
    name_temp = filename_prefix + "_temperature";
  else
    name_temp = "/dev/null";
  
  if( studyJets )
    mfpName = filename_prefix + "_mfp_jets";
  else
    mfpName = "/dev/null";
  
  if( studyCentralDensity )
    centralDensitiesName = filename_prefix + "_central_density";
  else
    centralDensitiesName = "/dev/null";
  
  // movie output
  if( theConfig->doOutput_movieOutputBackground() )
    oscarName_background = filename_prefix + "_background.oscar";
  else
    oscarName_background = "/dev/null";
  
  if( theConfig->doOutput_movieOutputJets() )
    oscarName_jet = filename_prefix + "_jets.oscar";
  else
    oscarName_jet = "/dev/null";

  printJpsiFugacity.open( name_fug.c_str(), ios::out );
  printTempInTube.open( name_temp.c_str(), ios::out );
  mfpJetsOutputFile.open( mfpName.c_str(), ios::out );
  centralDensitiesOutputFile.open( centralDensitiesName.c_str(), ios::out );
  oscarBackground.open( oscarName_background.c_str(), ios::out );
  oscarJets.open( oscarName_jet.c_str(), ios::out );
  //--------------------------------------

  jetTracking_PT = 10.0;
  
  

  if( studyPtSpectra )
  {
    //---- initialisation of PT-binning ----
    minPT = 1.4;
    maxPTSoft = 3.0;
    maxPT = 34.4;
    binWidthPT = 1.0;
    numberBinsPT = int(( maxPT - minPT + 0.001 ) / binWidthPT );
    
    tArrayOfDoubleVec tmpArray;
    for ( unsigned int i = 0; i < rapidityRanges.size(); i++ )
    {
      tmpArray.reset( new vector<double>[nTimeSteps+2] );     //+2 because of initial and final timesteps
      ptBins_gluons.push_back( tmpArray );
      tmpArray.reset( new vector<double>[nTimeSteps+2] );     //+2 because of initial and final timesteps
      ptBins_quarks.push_back( tmpArray );
      tmpArray.reset( new vector<double>[nTimeSteps+2] );     //+2 because of initial and final timesteps
      ptBins_ups.push_back( tmpArray );
      tmpArray.reset( new vector<double>[nTimeSteps+2] );     //+2 because of initial and final timesteps
      ptBins_downs.push_back( tmpArray );
      tmpArray.reset( new vector<double>[nTimeSteps+2] );     //+2 because of initial and final timesteps
      ptBins_stranges.push_back( tmpArray );
      tmpArray.reset( new vector<double>[nTimeSteps+2] );     //+2 because of initial and final timesteps
      ptBins_anti_ups.push_back( tmpArray );
      tmpArray.reset( new vector<double>[nTimeSteps+2] );     //+2 because of initial and final timesteps
      ptBins_anti_downs.push_back( tmpArray );
      tmpArray.reset( new vector<double>[nTimeSteps+2] );     //+2 because of initial and final timesteps
      ptBins_anti_stranges.push_back( tmpArray );
      tmpArray.reset( new vector<double>[nTimeSteps+2] );     //+2 because of initial and final timesteps
      ptBins_all.push_back( tmpArray );
    }

    for ( unsigned int ny = 0; ny < rapidityRanges.size(); ny++ )
    {
      for ( int nt = 0; nt < nTimeSteps + 2; nt++ )
      {
        ptBins_gluons[ny][nt].resize( numberBinsPT + 1, 0 );
        ptBins_quarks[ny][nt].resize( numberBinsPT + 1, 0 );
        ptBins_ups[ny][nt].resize( numberBinsPT + 1, 0 );
        ptBins_downs[ny][nt].resize( numberBinsPT + 1, 0 );
        ptBins_stranges[ny][nt].resize( numberBinsPT + 1, 0 );
        ptBins_anti_ups[ny][nt].resize( numberBinsPT + 1, 0 );
        ptBins_anti_downs[ny][nt].resize( numberBinsPT + 1, 0 );
        ptBins_anti_stranges[ny][nt].resize( numberBinsPT + 1, 0 );
        ptBins_all[ny][nt].resize( numberBinsPT + 1, 0 );
      }
    }
    
    //write the bin labels
    for ( int i = 0; i < numberBinsPT; i++ )
    {
      ptBinLabels.push_back( minPT + ( i * binWidthPT ) + ( binWidthPT / 2 ) );
    }
    cout << "number of bins: " << numberBinsPT << "  binWidth: " << binWidthPT << endl;
    //---- initialisation of PT-binning ----
    
    
    //---- initialisation of softPT-binning ----
    binWidthSoftPT = 0.1;
    numberBinsSoftPT = int(( maxPTSoft + 0.001 ) / binWidthSoftPT );
    ptBinsSoftAll_gluons = new vector<double>[nTimeSteps+1];
    ptBinsSoftAll_quarks = new vector<double>[nTimeSteps+1];
    ptBinsSoftAll_ups = new vector<double>[nTimeSteps+1];
    ptBinsSoftAll_downs = new vector<double>[nTimeSteps+1];
    ptBinsSoftAll_stranges = new vector<double>[nTimeSteps+1];
    ptBinsSoftAll_anti_ups = new vector<double>[nTimeSteps+1];
    ptBinsSoftAll_anti_downs = new vector<double>[nTimeSteps+1];
    ptBinsSoftAll_anti_stranges = new vector<double>[nTimeSteps+1];
    ptBinsSoftAll_all = new vector<double>[nTimeSteps+1];
    for ( int j = 0; j < nTimeSteps + 1; j++ )          //+1 because of initial timestep
      for ( int i = 0; i < numberBinsSoftPT; i++ )
      {
        ptBinsSoftAll_gluons[j].push_back( 0 );
        ptBinsSoftAll_quarks[j].push_back( 0 );
        ptBinsSoftAll_ups[j].push_back( 0 );
        ptBinsSoftAll_downs[j].push_back( 0 );
        ptBinsSoftAll_stranges[j].push_back( 0 );
        ptBinsSoftAll_anti_ups[j].push_back( 0 );
        ptBinsSoftAll_anti_downs[j].push_back( 0 );
        ptBinsSoftAll_anti_stranges[j].push_back( 0 );
        ptBinsSoftAll_all[j].push_back( 0 );
      }
    //write the bin labels
    for ( int i = 0; i < numberBinsSoftPT; i++ )
    {
      ptSoftBinLabels.push_back(( i * binWidthSoftPT ) + ( binWidthSoftPT / 2 ) );
    }
    //--------------------------------------
  }
  
  if( studyYDistribution || studyEtSpectra )
  {
    //------ initialisation of rapidity binning ------
    minY = -6.0;
    maxY = 6.0;
    binWidthY = 0.1;
    numberBinsY = int(( maxY - minY + 0.00001 ) / binWidthY );
      
    yBins_gluon = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
    yBins_up = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
    yBins_down = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
    yBins_strange = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
    yBins_anti_up = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
    yBins_anti_down = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
    yBins_anti_strange = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
    
    for ( int j = 0; j < nTimeSteps + 2; j++ )          //+2 because of initial and final timesteps
    {
      for ( int i = 0; i <= numberBinsY; i++ )
      {
        yBins_gluon[j].push_back( 0 );
        yBins_up[j].push_back( 0 );
        yBins_down[j].push_back( 0 );
        yBins_strange[j].push_back( 0 );
        yBins_anti_up[j].push_back( 0 );
        yBins_anti_down[j].push_back( 0 );
        yBins_anti_strange[j].push_back( 0 );
      }
    }
    
    //write the bin labels
    for ( int i = 0; i < numberBinsY; i++ )
    {
      yBinLabels.push_back( minY + ( i * binWidthY ) + ( binWidthY / 2 ) );
    }
    cout << "number of y-bins: " << numberBinsY << "  binWidth: " << binWidthY << endl;
    //--------------------------------------------------
  }

  if( studyEtSpectra )
  {
    //------ initialisation of transverse energy binning ------
    //------ use the same bins as for rapidity binning ------
    transverseEnergyGluons = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
    transverseEnergyQuarks = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
    transverseEnergyAntiQuarks = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
    
    for ( int j = 0; j < nTimeSteps + 2; j++ )          //+2 because of initial and final timesteps
    {
      for ( int i = 0; i <= numberBinsY; i++ )
      {
        transverseEnergyGluons[j].push_back( 0 );
        transverseEnergyQuarks[j].push_back( 0 );
        transverseEnergyAntiQuarks[j].push_back( 0 );
      }
    }
    //--------------------------------------------------
  }

  if(studyMFP)
  {
    cout << " !!! Only computes inverse rates (specific mean free path) and prints them out!!!\n 23-Photonproduction switched off, never mind the input-switch!!!" << endl;
    
    string full_filename = theConfig->getStandardOutputDirectoryName() + "/" +theConfig->getJobName() + "_MFPaveraged" +  ".dat";
    cout << "filename: " <<  full_filename << endl;
    cout << endl;
    fstream outfile( full_filename.c_str(), ios::out | ios::trunc);
    outfile << "#AverageInverseRates" << endl;
    outfile << "#t" << "\t" << "#Lambda1"<< "\t" << "#Lambda2"<< "\t" << "#Lambda3"<<endl;
    outfile.close();     
  }
  
  if( studyJpsi )
  {
    numberJpsi_all_time = new int[nTimeSteps+1];
    numberJpsi_ini_time = new int[nTimeSteps+1];
    numberJpsi_midPseudoRap_all_time = new int[nTimeSteps+1];
    numberJpsi_midPseudoRap_ini_time = new int[nTimeSteps+1];
    numberJpsi_forwardPseudoRap_all_time = new int[nTimeSteps+1];
    numberJpsi_forwardPseudoRap_ini_time = new int[nTimeSteps+1];
    numberJpsi_midNormRap_all_time = new int[nTimeSteps+1];
    numberJpsi_midNormRap_ini_time = new int[nTimeSteps+1];
    numberJpsi_forwardNormRap_all_time = new int[nTimeSteps+1];
    numberJpsi_forwardNormRap_ini_time = new int[nTimeSteps+1];
  //   numberJpsi_midSpaceTimeRap_time = new int[nTimeSteps+1];
    numberJpsiProd_time = new int[nTimeSteps+1];
    numberJpsiDiss_time = new int[nTimeSteps+1];
    numberJpsiDissTd_time = new int[nTimeSteps+1];
    numberCCbGG_time = new int[nTimeSteps+1];
    timestepAnalysed = new bool[nTimeSteps+1];
  //   charmJetEnergy = new double[nTimeSteps+1];

    for ( int i = 0; i < nTimeSteps + 1; i++ )
    {
      numberJpsi_all_time[i] = numberJpsi_ini_time[i] = numberJpsiDiss_time[i] = numberJpsiDissTd_time[i] = 
      numberJpsi_midPseudoRap_all_time[i] = numberJpsi_midPseudoRap_ini_time[i] = numberJpsi_forwardPseudoRap_all_time[i] = numberJpsi_forwardPseudoRap_ini_time[i] = 
      numberJpsi_midNormRap_all_time[i] = numberJpsi_midNormRap_ini_time[i] = numberJpsi_forwardNormRap_all_time[i] = numberJpsi_forwardNormRap_ini_time[i] = 
      numberJpsiProd_time[i] = numberCCbGG_time[i] = 0;
      timestepAnalysed[i] = false;
  //     charmJetEnergy[i] = 0.0;
    }
    
    theInterpolation_nJpsi.configure();
  }
  
  if( theConfig->doOutput_progressLog() )
  {
    string filename = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + "_progressLog";
    progressLogFile.open( filename.c_str(), ios::out | ios::trunc );
  }
  showerParticlesInEvent.resize( static_cast<int>( theConfig->getNumberOfParticlesToAdd() / 2 ) );
  
  //PhotonEnergyBins configuration, to avoid unneccessary errors.   
  PhotondNOverTwoPiptdydptBin.setMinMaxN( 0.0, 3.0, 100 );
  cellV2.setMinMaxN(-1.0,1.0,20);
  cellV2Weighted.setMinMaxN(-1.0,1.0,20);
  
  if( studyPhotons )
  { 
    //PhotonEnergyBins configuration   
    PhotondNOverTwoPiptdydptBin.setMinMaxN( 0.0, 3.0, 100 );
    cellV2.setMinMaxN(-1.0,1.0,20);
    cellV2Weighted.setMinMaxN(-1.0,1.0,20);
    
    //---- initialisation of PT-binning ----
    minPTPhotons = 0.0;
    maxPTPhotons = 5.0;
    binWidthPTPhotons = 0.03;
    numberBinsPTPhotons = int(( maxPTPhotons - minPTPhotons + 0.001 ) / binWidthPTPhotons );
    
    tArrayOfDoubleVec tmpArray2;
    for ( unsigned int i = 0; i < rapidityRanges.size(); i++ )
    { 
      tmpArray2.reset( new vector<double>[nTimeSteps+2] );     //+2 because of initial and final timesteps
      ptBins_photons.push_back( tmpArray2 ); 
    } 
    for ( unsigned int ny = 0; ny < rapidityRanges.size(); ny++ )
    {
      ptBins_photons[ny][0].resize( numberBinsPTPhotons + 1, 0 );
    }    
    //write the bin labels
    for ( int i = 0; i < numberBinsPTPhotons; i++ )
    {
      ptBinLabelsPhotons.push_back( minPT + ( i * binWidthPTPhotons ) + ( binWidthPTPhotons / 2 ) );
    }
    cout << "Photon-binning: Number of bins: " << numberBinsPTPhotons << "  binWidth: " << binWidthPTPhotons << endl;
    //---- initialisation of PT-binning ----
    
  }
  if(studyThermalisation)
  {
    cout << "Thermalization analysis" << endl;
    
    radiusAnalysisTube = theConfig->getAnalysisTubeRadius();
    dEtaAnalysisTube = theConfig->getAnalysisTubedEta();
    
    cout << "Analysis in Tube: r= " << radiusAnalysisTube << "   dEta = +- " << dEtaAnalysisTube << endl;
    
        //PhotonEnergyBins configuration   
    partonEnergies.setMinMaxN( 0.0, 3.0, 20 );
    quarkEnergies.setMinMaxN( 0.0, 3.0, 20 );
    gluonEnergies.setMinMaxN( 0.0, 3.0, 20 );
  }
}


analysis::~analysis()
{
  //---- get time and calculate real time needed by simulation ----
  time_t end;
  time( &end );
  int secs = ( int )difftime( end, start );
  int hours = secs / 3600, minutes = ( secs % 3600 ) / 60, seconds = secs - 3600 * hours - 60 * minutes;
  cout << hours << "h" << minutes << "m" << seconds << "s" <<  endl;
  //---------------------------------------------------------------

  oscarBackground.close();
  oscarJets.close();
  jetTrackerOutput();
}


void analysis::handle_output_studies( OUTPUT_SCHEME _outputScheme )
{
  // By default switches are set off (please do not edit here, but define output scheme and change below, for more information see https://th.physik.uni-frankfurt.de/~bamps/cgi-bin/trac/wiki/RestrictedWiki/AnalysisOutputScheme):
  studyHQ = false; // heavy quarks
  studyJpsi = false; // jpsi
  studyTempInTube = false; // temperature in a centred tube in heavy ion collision
  studyTempAndVelocity = false; // write out temperature and velocity for each cell
  studyPtSpectra = false; // transverse momentum pt spectra (partly the same as in v2RAA class)
  studyEtSpectra = false; // transverse energy spectra
  studyYDistribution = false; // rapidity distribution
  studyJets = false; // study jets
  studyCentralDensity = false; // density in central part of collision
  studyBackground = false; // print also properties like v2, RAA of background
  studyScatteredMediumParticles = false; // print scattered medium particles
  studyPhotons = false; // transverse momentum spectra of photons
  studyNumberOfBackgroundQuarks = false; // study Number Of Background quarks and antiquarks and gluons
  studySpatialPhotons = false;
  studyTempAndVelocity = false;
  studyTempCustom = false;
  studyMFP = false; //study the specific mean free path
  studyDileptons = false;
  studyThermalisation = false;
  studyPartons = false;
  
  //---- defining standard rapidity ranges ----
  // only use positiv ranges since the investigated collision systems usually are symmetric in +-y and we therefore only compare the absolute value of y
  analysisRapidityRange yRange;
  yRange.reset( 0, 0.5 );
  rapidityRanges.push_back(yRange);
  yRange.reset( 0, 0.8 );
  rapidityRanges.push_back(yRange);
  yRange.reset( 0, 1.0 );
  rapidityRanges.push_back(yRange);
  yRange.reset( 0, 1.5 );
  rapidityRanges.push_back(yRange);
  yRange.reset( 0, 2.0 );
  rapidityRanges.push_back(yRange);
  yRange.reset( 0, 2.5 );
  rapidityRanges.push_back(yRange);  
  yRange.reset( 0, infinity );
  rapidityRanges.push_back(yRange);  
  //---- defining rapidity ranges ----

  
  // add a new case for your outpute scheme which you can create in configuration.h
  switch ( _outputScheme )
  {
    case studySpecificMFP:
      studyMFP = true;
      break;
    case studyOnlyEtSpectra:  
      studyEtSpectra = true;
      break;
    case studyNumberOfParticles:
      studyNumberOfBackgroundQuarks = true;
      break;
    case SpatialPhotons: 
      studySpatialPhotons = true;
      break; 
    case temperatureCustom:
      studyTempCustom = true;
      break;
    case thermalization:
      studyThermalisation = true;
      break;
    case temperature_velocity_cells:
      studyTempAndVelocity = true;
      break;
    case phenix_hq_electrons:
      studyHQ = true;
      
      rapidityRanges.clear();
      yRange.reset( 0, 0.35 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.5 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 1.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 2.0 );
      rapidityRanges.push_back(yRange);
      break;
    case alice_hq_electrons:
      studyHQ = true;
      
      rapidityRanges.clear();
      yRange.reset( 0, 0.8 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 1.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 2.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.5 );
      rapidityRanges.push_back(yRange);
      break;
    case alice_hq_muons:
      studyHQ = true;

      rapidityRanges.clear();
      yRange.reset( 2.5, 4.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.35 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.5 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.8 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 1.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 2.0 );
      rapidityRanges.push_back(yRange);
      break;
    case alice_hq_dmesons:
      studyHQ = true;
      
      rapidityRanges.clear();
      yRange.reset( 0, 0.5 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.8 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 1.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 2.0 );
      rapidityRanges.push_back(yRange);
      break;
    case cms_hq_nonPromptJpsi:
      studyHQ = true;
      
      rapidityRanges.clear();
      yRange.reset( 0, 2.4 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.35 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.5 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 1.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 2.0 );
      rapidityRanges.push_back(yRange);
      break;
    case alice_hq_nonPromptJpsi:
      studyHQ = true;
      
      rapidityRanges.clear();
      yRange.reset( 0, 0.9 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.35 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.5 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 1.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 2.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 2.4 );
      rapidityRanges.push_back(yRange);
      break;
    case phenix_jpsi:
      studyJpsi = true;
      studyHQ = true;
      studyTempInTube = true;

      rapidityRanges.clear();
      yRange.reset( 0, 0.35 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 1.2, 2.2 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.5 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 1.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 2.0 );
      rapidityRanges.push_back(yRange);
      break;
    case alice_jpsi:
      studyJpsi = true;
      studyHQ = true;
      studyTempInTube = true;

      rapidityRanges.clear();
      yRange.reset( 0, 0.9 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 2.5, 4.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.35 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.5 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 1.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 2.0 );
      rapidityRanges.push_back(yRange);
      break;
    case cms_jpsi:
      studyJpsi = true;
      studyHQ = true;
      studyTempInTube = true;

      rapidityRanges.clear();
      yRange.reset( 0, 2.4 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.9 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 2.5, 4.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.35 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.5 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 1.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 2.0 );
      rapidityRanges.push_back(yRange);
      break;
    case background_jets:
      studyScatteredMediumParticles = true;
      break;
      
    case light_parton_lhc:      
      rapidityRanges.clear();
      yRange.reset( 0, 0.8 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.5 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 1.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 2.0 );
      rapidityRanges.push_back(yRange);
      break;
      
    case light_parton_phenix:      
      rapidityRanges.clear();
      yRange.reset( 0, 0.35 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.5 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.8 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 1.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 2.0 );
      rapidityRanges.push_back(yRange);
      break;
    
    case central_densities:
      studyCentralDensity = true;
      break;
      
    case photons:     
      rapidityRanges.clear();
      yRange.reset( 0, 0.35 );//rapidityRange[0] // RHIC
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.7 );//1                 // LHC
      rapidityRanges.push_back(yRange);
      studyPhotons=true;
      break;
    
    case photons_plus_background:  
      rapidityRanges.clear();
      yRange.reset( 0, 0.35 );//rapidityRange[0] // RHIC
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.7 );//1                 // LHC
      rapidityRanges.push_back(yRange);
      studyBackground = true;
      studyPhotons=true;
      break;    
      
    case dileptonStudies:     
      rapidityRanges.clear();
      yRange.reset( 0, 0.5 );//rapidityRange[0] // RHIC
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.7 );//1                 // LHC
      rapidityRanges.push_back(yRange);
      studyDileptons=true;
      break;    
    
    case only_background:
      rapidityRanges.clear();
      yRange.reset( 0, 0.35 );//rapidityRange[0] // RHIC
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.5 );//1                 // LHC
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 1.0 );//1                 // LHC
      rapidityRanges.push_back(yRange);
      studyPartons=true;
      break;      
      
    default:
      break;
  } 
}






void analysis::collectPtDataInitial()
{
  if( studyPtSpectra )
  {
    ptDistribution( gluon, addedParticles, addedParticles.size(), 0 );
    ptDistribution( light_quark, addedParticles, addedParticles.size(), 0 );
    ptDistribution( allFlavors, addedParticles, addedParticles.size(), 0 );
    ptSoftDistribution( gluon, particles_atTimeNow, particles_atTimeNow.size(), 0 );
    ptSoftDistribution( light_quark, particles_atTimeNow, particles_atTimeNow.size(), 0 );
    ptSoftDistribution( allFlavors, particles_atTimeNow, particles_atTimeNow.size(), 0 );
  }
  
}



void analysis::collectPtData( const int step )
{
  //are issued starting with 0
  //0 is reserved for initial output, thus step is incremented
  if( studyPtSpectra )
  {
    ptDistribution( gluon, addedParticles, addedParticles.size(), step + 1 );
    ptDistribution( light_quark, addedParticles, addedParticles.size(), step + 1 );
    ptDistribution( allFlavors, addedParticles, addedParticles.size(), step + 1 );
    ptSoftDistribution( gluon, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
    ptSoftDistribution( light_quark, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
    ptSoftDistribution( allFlavors, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
  }
  
  if( studyPhotons)
  {  
    
  }
}


void analysis::collectYDataInitial()
{
  if( studyYDistribution )
    collectYData( -1 );
}


void analysis::collectYData( const int step )
{
  //are issued starting with 0
  //0 is reserved for initial output, thus step is incremented
  if( studyYDistribution )
  {
    yDistribution( gluon, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
    yDistribution( up, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
    yDistribution( down, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
    yDistribution( strange, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
    yDistribution( anti_up, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
    yDistribution( anti_down, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
    yDistribution( anti_strange, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
  }
}


void analysis::collectEtDataInitial()
{
  if( studyEtSpectra )
  {
    transverseEnergyDistribution( gluon, particles_atTimeNow, particles_atTimeNow.size(),0 );
    transverseEnergyDistribution( light_quark, particles_atTimeNow, particles_atTimeNow.size(),0 );
    transverseEnergyDistribution( anti_light_quark, particles_atTimeNow, particles_atTimeNow.size(), 0 );
    printEtDistribution(0);
  }
}


void analysis::collectEtData( const int step )
{
  //are issued starting with 0
  //0 is reserved for initial output, thus step is incremented
  //if( studyEtSpectra )
  //{
  //  transverseEnergyDistribution( gluon, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
  //  transverseEnergyDistribution( light_quark, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
  //  transverseEnergyDistribution( anti_light_quark, particles_atTimeNow, particles_atTimeNow.size(), step + 1 );
  //  printEtDistribution(step+1);
 // }
}

void analysis::printEtDistribution(const int step )
{
  string filename,timename;
  stringstream ss,sst;

  if ( step == 0 )
  {
    timename = "initial";
  }
  else if ( step == nTimeSteps + 1 )
  { 
    timename = "final";
  }
  else
  {
    ss << step;
    timename = "timestep_" + ss.str()+"_time_";
    sst << tstep[step];
    timename += sst.str() + "fm";
  }

  //creates filename, for example: "./output/run32_step1.f1", "./output/run32_initial.f1" or the likes
  filename = filename_prefix + "_dEtdy" + "_" + timename + ".dat";
  fstream file( filename.c_str(), ios::out | ios::trunc );

  //---- print header if file is empty ----
  time_t end;
  time( &end );

  file.seekp( 0, ios::end );
  long size = file.tellp();
  file.seekp( 0, ios::beg );
  if ( size == 0 )
    printHeader( file, dEtdy, end );
  
  file << "#Et_Gluons\t#EtQuarks\t#EtQuarks/dy\t#EtGluons/dy" << endl;
 for ( int i = 0; i <= numberBinsY; i++ )
  {
    file <<  minY + (  (i*binWidthY )  +  (  (i+1)*binWidthY ) )/2.;
    file << "\t";
    file <<   transverseEnergyGluons[step][i];
    file << "\t";
    file <<   transverseEnergyQuarks[step][i];
    file << "\t";
    file <<   transverseEnergyGluons[step][i]/binWidthY;
    file << "\t";
    file <<   transverseEnergyQuarks[step][i]/binWidthY;
    file << endl;
  }
  

 
 file.close();
 
 
 
  
  
}





/** Bin the Pt distribution for photons with two methods. For double-checking and some tests.
 * 
 * @param[in] PTofThisSinglePhoton PT value [GeV] of single photon
 * @param[in] etaOfThisSinglePhoton eta value of single photon
 * @param[in] EofSinglePhoton E value [GeV] of single photon
 */
void analysis::PtDistributionPhotons( const double PTofThisSinglePhoton, const double etaOfThisSinglePhoton, const double EofSinglePhoton )
{
  //Binning-Class:
  if ( fabs( etaOfThisSinglePhoton ) >= rapidityRanges[2].yleft && fabs( etaOfThisSinglePhoton ) <= rapidityRanges[2].yright )
  {  
    PhotondNOverTwoPiptdydptBin.add(PTofThisSinglePhoton); //PTofThisSinglePhoton/(2*M_PI*DeltaRapidity*PTofThisSinglePhoton ) 
  }

  //New Method:
  double pt, y;
 
  pt = PTofThisSinglePhoton;
  y = etaOfThisSinglePhoton;
    
  // is it in the possible pt range at all?
  if ( pt < maxPTPhotons && pt >= minPTPhotons )
  {
    // individually check for each rapidity range whether this particle needs to be binned
    for ( int yRangeIndex = 0; yRangeIndex < rapidityRanges.size(); yRangeIndex++ )
    {
      if ( fabs( y ) >= rapidityRanges[yRangeIndex].yleft && fabs( y ) <= rapidityRanges[yRangeIndex].yright )
      {
        if ( pt == minPTPhotons )  // a special case
        {
          ++ptBins_photons[yRangeIndex][0][0];  // actually bin the pt of the particle
        }
        else
        {
          ++ptBins_photons[yRangeIndex][0][int(( pt - minPTPhotons )/binWidthPTPhotons )];  // actually bin the pt of the particle
        }
      }
    }
  }  
}

/** Bin the cellV2 distribution
 * 
 * @param[in] cellV2ofThisCell
 */
void analysis::cellV2Distribution( const double cellV2ofThisCell, const unsigned int weight)
{ 
  cellV2.add(cellV2ofThisCell); 
  if(weight>0)
  {
    for(unsigned int i = 1;i<weight;i++)
    {
      cellV2Weighted.add(cellV2ofThisCell);
    }
  }  
}


void analysis::printCellV2Distribution( const double nexttime, unsigned int timestepCount)
{
  time_t end;
  time( &end );
  
  string timename,timeInFm;
  stringstream ss,ss2;

  ss << timestepCount;
  ss2 << nexttime;
  timename = ss.str();
  timeInFm = ss2.str();
  
  //NORMAL, not weighted:  
  string filename1 = filename_prefix + "_CellV2Distribution_time_" + timename ;
  fstream file1( filename1.c_str(), ios::out | ios::trunc );

  //print header if file is empty
  file1.seekp( 0, ios::end );
  long size = file1.tellp();
  file1.seekp( 0, ios::beg );
  if ( size == 0 )
    printHeader( file1, cellV2DistributionHeader,end );

  file1 << "#time in fm:" << endl;
  file1 << "#4321\t" << timeInFm << endl;
  file1 << endl;
  file1 << "#Mean cell v2" << endl;
  file1 << cellV2.getMean()<<"\t0.\t1234\t1234" << endl;
  file1 << cellV2.getMean()<<"\t10000\t1234\t1234" << endl;
  file1 << endl;
  cellV2.print(file1);
  file1.close();
  cellV2.deleteBins();

  //Weighhted with Photon number produced in cell:   
  string filename2 = filename_prefix + "_CellV2DistributionWeightedWithPhotons_time_" + timename ;
  fstream file2( filename2.c_str(), ios::out | ios::trunc );

  file2 << "#Weighted wit number of photons produced in cell" << endl;
  
  //print header if file is empty
  file2.seekp( 0, ios::end );
  long size2 = file2.tellp();
  file2.seekp( 0, ios::beg );
  if ( size2 == 0 )
    printHeader( file2, cellV2DistributionHeader,end );

  file2 << "#time in fm:" << endl;
  file2 << "#4321\t" << timeInFm << endl;
  file2 << endl;
  file2 << "#Mean cell v2" << endl;
  file2 << cellV2Weighted.getMean()<<"\t0.\t1234\t1234" << endl;
  file2 << cellV2Weighted.getMean()<<"\t10000\t1234\t1234" << endl;
  file2 << endl;
  cellV2Weighted.print(file2);
  file2.close();
  cellV2Weighted.deleteBins();
  
}


/** Photonspectra output
 *  Prints the Spectra from the function PtDistributionPhotons out. For Double-Checking and Tests.
 *  Generates two files: 
 *    ...filename_prefix + "_PhotondNOverTwoPiptdydpt_Old
 *    ...filename_prefix + "_PhotondNOverTwoPiptdydpt_New
 */
void analysis::photonSpectrumOutput()
{
  time_t end;
  time( &end );
  
  cout << "testp: " << theConfig->getTestparticles() << endl;
  
  //****** Old Method:
  string filename1 = filename_prefix + "_PhotondNOverTwoPiptdydpt_Old";
  fstream file1( filename1.c_str(), ios::out | ios::trunc );

  //print header if file is empty
  file1.seekp( 0, ios::end );
  long size = file1.tellp();
  file1.seekp( 0, ios::beg );
  if ( size == 0 )
    printHeader( file1, photonPtDist,end );

  
  double DeltaRapidity = fabs(rapidityRanges[2].yright - rapidityRanges[2].yleft)*2.0;
  
  //WARNING: Divide manually by K-Factor later, if not 1!!!!!!!!!!!!!!!!!!!!!
  PhotondNOverTwoPiptdydptBin.print(file1,1.0/(2.0*M_PI*DeltaRapidity*theConfig->getTestparticles()));
  file1.close();

  //***** New Method:
  
  
  string filename2 = filename_prefix + "_PhotondNOverTwoPiptdydpt_New";
  fstream file2( filename2.c_str(), ios::out | ios::trunc );

  //---- print header if file is empty ----
  //file2.seekp( 0, ios::end );
  //long size2 = file2.tellp();
  //file2.seekp( 0, ios::beg );
  //if ( size2 == 0 )
  //  printHeader( file2, ptSpectrum, end );
  //---------------------------------------

  /*cout << "Rapidity ranges" << endl;
  cout << rapidityRanges[0].yright << " " << rapidityRanges[0].yleft << endl;
  cout << rapidityRanges[1].yright << " " << rapidityRanges[1].yleft << endl;  
  cout << rapidityRanges[2].yright << " " << rapidityRanges[2].yleft << endl;
  cout << rapidityRanges[3].yright << " " << rapidityRanges[3].yleft << endl;*/
  
  //------------- the actual output ---------------
  for ( unsigned int yRangeIndex = 0; yRangeIndex < rapidityRanges.size(); yRangeIndex++ )
  {
    file2 << "#y in +- [" << rapidityRanges[yRangeIndex].yleft << ", " << rapidityRanges[yRangeIndex].yright << "]" << sep; 
  } 
  file2 << endl;
 
  for ( int i = 0; i < numberBinsPTPhotons; i++ )
  {
    file2 << ptBinLabelsPhotons[i] << sep;
    for ( unsigned int yRangeIndex = 0; yRangeIndex < rapidityRanges.size(); yRangeIndex++ )
    {
      DeltaRapidity = fabs(rapidityRanges[yRangeIndex].yright - rapidityRanges[yRangeIndex].yleft)*2.0;
      //WARNING: Divide manually by K-Factor later, if not 1!!!!!!!!!!!!!!!!!!!!!
      file2 <<  ptBins_photons[yRangeIndex][0][i]/(2.0*M_PI*DeltaRapidity*binWidthPTPhotons*ptBinLabelsPhotons[i]*theConfig->getTestparticles())  << sep;
    } 
    file2 << endl;
  }
  file2 << endl << endl;
  //-------------------------------------------------

  file2.close();
   
}


void analysis::initialOutput()
{
  if ( v2output )
    computeV2RAA( "initial", 0  );

  if ( studyJpsi ) // to consider charm annihaltion is just useful if added particles can scatter
  { 
    jpsiEvolution( 0 );
    ini_charm_correlations();
    writeJpsiFugacityOutput( 0 );
  }
  
  if( studyTempInTube)
     writeTempInTube( 0 );
  
  if ( dndyOutput )
    print_dndy( "initial" );
  
  if ( studyParticleOutput )
    particleOutput( 0 );
  
  if (studyTempCustom) //Moritz
  {
    writeTempCustom(0);   
    writeTempInTube(0);
    writeCustomTube(0);
  }    
  
  if(studyThermalisation)
  {
    cout << "Thermalisation analysis initial timestep " << endl;
    writeCustomTube(0);
  }
  
  if(studyPartons)
  {
    computeV2RAA( "initial", 0  );    
  }
  
  //if ( studyPhotons ){ do nothing because no initial photons there :-)}
  
  
}



void analysis::intermediateOutput( const int nn )
{
  string name;
  stringstream ss;

  ss << tstep[nn];
  name = ss.str() + "fm";

 if( studyEtSpectra)
  { 
    cout << "Analyse Et Spectra at analysis timestep " << nn << " at time " <<  tstep[nn] << endl;
    transverseEnergyDistribution( gluon, particles_atTimeNow, particles_atTimeNow.size(), nn+1 );
    transverseEnergyDistribution( light_quark, particles_atTimeNow, particles_atTimeNow.size(), nn+1 );
    printEtDistribution(nn+1);
  }
  
  if( studyNumberOfBackgroundQuarks )
  {
    saveNumberOfMediumParticles(nn);
  }
  if(studyPartons)
  {
     computeV2RAA( name, tstep[nn] );    
  }
  
  if ( v2output && v2outputIntermediateSteps )
    computeV2RAA( name, tstep[nn] );

  if ( studyJpsi ) // to consider charm annihaltion is just useful if added particles can scatter
  {
    jpsiEvolution( nn + 1 );
    writeJpsiFugacityOutput( nn + 1 );
  }
  
  if( studyTempInTube)
     writeTempInTube( nn + 1 ); 
  
  if ( studyTempAndVelocity ) // hydro
    writeTempAndVel( nn + 1  );
  
  if (studyTempCustom) //Moritz
  {
    writeTempCustom(nn+1);   
    writeTempInTube(nn+1);
    writeCustomTube(nn+1);
  }  
  
  if(studyThermalisation)
  {
    //cout << "Thermalisation analysis at analysis timestep " << nn << " at time " <<  tstep[nn] << endl;
    writeCustomTube(nn);
  }
  
//   if ( charmTestJet )
//    analyseCharmTestJetEvolution( nn + 1 );
  
  if ( dndyOutput )
    print_dndy( name );

  if (studySpatialPhotons)
  {
    writePhotonSpaceProfile(nn+1);   
  }
  
}



void analysis::finalOutput( const double _stoptime )
{
 int step = _stoptime;
 
  if( studyNumberOfBackgroundQuarks )
  {
    printNumberOfMediumParticles();
  }
 
 if( studyEtSpectra)
  { 
    cout << "Analyse Et Spectra at final timestep (maybe the same output as the last analysis timestep.) " << nTimeSteps +1 << endl;
    transverseEnergyDistribution( gluon, particles_atTimeNow, particles_atTimeNow.size(), nTimeSteps +1 );
    transverseEnergyDistribution( light_quark, particles_atTimeNow, particles_atTimeNow.size(), nTimeSteps +1);
    printEtDistribution(nTimeSteps +1);
  }
  
  if( studyPtSpectra )
  {  
    printPtSpectra( gluon );
    printPtSpectra( light_quark );
    printPtSpectra( allFlavors );
    printSoftPtSpectra( gluon );
    printSoftPtSpectra( light_quark );
    printSoftPtSpectra( allFlavors );
  }
  
  if ( studyParticleOutput )
    particleOutput( nTimeSteps );
  
  if( studyYDistribution )
    printYDistribution();
  
  if ( v2output )
    computeV2RAA( "final", _stoptime );

  if ( studyPartons )
    computeV2RAA( "final", _stoptime );  
  
  if ( particleCorrelationsOutput )
  {
    onePartclCorrelations();
    twoPartclCorrelations();
  }

  if ( studyJpsi ) // to consider charm annihaltion is just useful if added particles can scatter
  {
    printJpsiEvolution();
//     jpsi_correlations();
  }
  
//   if ( charmTestJet )
//    analyseCharmTestJet();
  
  if( hadronization_hq && mesonDecay )
    analyseAngleDe();
  
  if ( dndyOutput )
    print_dndy( "final" );
  if( studyScatteredMediumParticles )
  {
    scatteredMediumParticlesOutput( nTimeSteps );
    mediumParticlesOutput( nTimeSteps );
  }
}



void analysis::movieOutput( const int step, const int jumpSteps )
{
  writePartclMovie( addedParticles, addedParticles.size(), oscarJets, step, jumpSteps );
//   writePartclMovie( particles_atTimeNow, particles_atTimeNow.size(), oscarBackground, step );
}


void analysis::movieOutputMedium( const int step, const int jumpSteps )
{
  writePartclMovie( particles_atTimeNow, particles_atTimeNow.size(), oscarBackground, step, jumpSteps );
}



void analysis::printYDistribution()
{
  time_t end;
  time( &end );

  string filename = filename_prefix + "_" + "_yDistribution";
  fstream file( filename.c_str(), ios::out | ios::trunc );

  //---- print header if file is empty ----
  file.seekp( 0, ios::end );
  long size = file.tellp();
  file.seekp( 0, ios::beg );
  if ( size == 0 )
    printHeader( file, rapidityDistribution, end );
  //---------------------------------------
  
  
  for ( int j = 0; j <= nTimeSteps; j++ )
  {
    if ( tstep[j-1] <= theConfig->getRuntime() || j == 0 || j == nTimeSteps )
    {
      if ( j == 0 )
      {
        file << "#t = 0" << endl;
      }
      else
      {
        file << "#t = " << tstep[j-1] << endl;
      }
      for ( int i = 0; i < numberBinsY; i++ )
      {
        file << yBinLabels[i] << sep << yBins_gluon[j][i] << sep << yBins_up[j][i] << sep << yBins_down[j][i] 
        << sep << yBins_strange[j][i] << sep << yBins_anti_up[j][i] << sep << yBins_anti_down[j][i] << sep 
        << yBins_anti_strange[j][i] << endl; 
      }
      file << endl << endl;
    }
  }
  
  file.close();
  
  
  //-------- quark numbers summed over all rapidities
  
  filename =  filename_prefix + "_" + "_quarkNumbers";
  file.open( filename.c_str(), ios::out | ios::trunc );
  
  //---- print header if file is empty ----
  file.seekp( 0, ios::end );
  size = file.tellp();
  file.seekp( 0, ios::beg );
  if ( size == 0 )
    printHeader( file, quarkNumbers, end );
  //---------------------------------------
  
  double nGluons, nUps, nDowns, nStranges, nAntiUps, nAntiDowns, nAntiStranges;
  
  for ( int j = 0; j <= nTimeSteps; j++ )
  {
    nGluons = nUps = nDowns = nStranges = nAntiUps = nAntiDowns = nAntiStranges = 0;
    if ( tstep[j-1] <= theConfig->getRuntime() || j == 0 || j == nTimeSteps )
    {
      if ( j == 0 )
      {
        file << 0 << sep;
      }
      else
      {
        file << tstep[j-1] << sep;
      }
      
      for ( int i = 0; i < numberBinsY; i++ )
      {
        nGluons += yBins_gluon[j][i];
        nUps += yBins_up[j][i];
        nDowns += yBins_down[j][i];
        nStranges += yBins_strange[j][i];
        nAntiUps += yBins_anti_up[j][i];
        nAntiDowns += yBins_anti_down[j][i];
        nAntiStranges += yBins_anti_strange[j][i];
      }
      file << nGluons << sep << nUps << sep << nDowns << sep << nStranges << sep << nAntiUps << sep << nAntiDowns << sep << nAntiStranges << endl;    
    }
  }
  
  
  file.close();
}




void analysis::printPtSpectra( const FLAVOR_TYPE _flavTypeToComputeFor )
{
  time_t end;
  time( &end );

  // populate temporary _ptBins with temporary "references" to the actual data structures according to the requested flavor
  tVecOfArrayOfDoubleVec _ptBins;
  for ( unsigned int yRange = 0; yRange < rapidityRanges.size(); yRange++ )
  {
    switch ( _flavTypeToComputeFor )
    {
      case gluon:
        _ptBins.push_back( ptBins_gluons[yRange] );
        break;
      case light_quark:
        _ptBins.push_back( ptBins_quarks[yRange] );
        break;
      case up:
        _ptBins.push_back( ptBins_ups[yRange] );
        break;
      case down:
        _ptBins.push_back( ptBins_downs[yRange] );
        break;
      case strange:
        _ptBins.push_back( ptBins_stranges[yRange] );
        break;
      case anti_up:
        _ptBins.push_back( ptBins_anti_ups[yRange] );
        break;
      case anti_down:
        _ptBins.push_back( ptBins_anti_downs[yRange] );
        break;  
      case anti_strange:
        _ptBins.push_back( ptBins_anti_stranges[yRange] );
        break;
      case allFlavors:
        _ptBins.push_back(ptBins_all[yRange]);
        break;
      default:
        string errMsg = "error in ptDistribution, flavor not specified";
        throw eAnalysis_error( errMsg );
        break;
    }
  }

  string type;
  if ( _flavTypeToComputeFor == gluon )
  {
    type = "gluon";
  }
  else if ( _flavTypeToComputeFor == light_quark || _flavTypeToComputeFor == anti_light_quark )
  {
    type = "quark";
  }
  if ( _flavTypeToComputeFor == up )
  {
    type = "up";
  }
  if ( _flavTypeToComputeFor == down )
  {
    type = "down";
  }
  if ( _flavTypeToComputeFor == strange )
  {
    type = "strange";
  }
  if ( _flavTypeToComputeFor == anti_up )
  {
    type = "anti_up";
  }
  if ( _flavTypeToComputeFor == anti_down )
  {
    type = "anti_down";
  }
  if ( _flavTypeToComputeFor == anti_strange )
  {
    type = "anti_strange";
  }
  else if ( _flavTypeToComputeFor == allFlavors )
  {
    type = "allFlavors";
  }

  string filename = filename_prefix + "_" + type + "_spectra.f2";
  fstream file( filename.c_str(), ios::out | ios::trunc );

  //---- print header if file is empty ----
  file.seekp( 0, ios::end );
  long size = file.tellp();
  file.seekp( 0, ios::beg );
  if ( size == 0 )
    printHeader( file, ptSpectrum, end );
  //---------------------------------------

  //------------- the actual output ---------------
  for ( unsigned int yRangeIndex = 0; yRangeIndex < rapidityRanges.size(); yRangeIndex++ )
  {
    file << "#y in +- [" << rapidityRanges[yRangeIndex].yleft << ", " << rapidityRanges[yRangeIndex].yright << "]" << endl;
    for ( int i = 0; i < numberBinsPT; i++ )
    {
      file << ptBinLabels[i] << sep;
      for ( int j = 0; j <= nTimeSteps; j++ )
      {
        if ( j == 0 || j == nTimeSteps )
          file << _ptBins[yRangeIndex][j][i] << sep;
        else if ( tstep[j-1] <= theConfig->getRuntime() )
          file << _ptBins[yRangeIndex][j][i] << sep;
      }
      file << endl;
    }
    file << endl << endl;
  }
  //-------------------------------------------------

  file.close();
}



void analysis::printSoftPtSpectra( const FLAVOR_TYPE _flavTypeToComputeFor )
{
  time_t end;
  time( &end );

  vector<double> * _ptBinsSoftAll;

  if ( _flavTypeToComputeFor == gluon )
  {
    _ptBinsSoftAll = ptBinsSoftAll_gluons;
  }
  else if ( _flavTypeToComputeFor == light_quark )
  {
    _ptBinsSoftAll = ptBinsSoftAll_quarks;
  }
  else if ( _flavTypeToComputeFor == up )
  {
    _ptBinsSoftAll = ptBinsSoftAll_ups;
  }
  else if ( _flavTypeToComputeFor == down )
  {
    _ptBinsSoftAll = ptBinsSoftAll_downs;
  }
  else if ( _flavTypeToComputeFor == strange )
  {
    _ptBinsSoftAll = ptBinsSoftAll_stranges;
  }
  else if ( _flavTypeToComputeFor == anti_up )
  {
    _ptBinsSoftAll = ptBinsSoftAll_anti_ups;
  }
  else if ( _flavTypeToComputeFor == anti_down )
  {
    _ptBinsSoftAll = ptBinsSoftAll_anti_downs;
  }
  else if ( _flavTypeToComputeFor == anti_strange )
  {
    _ptBinsSoftAll = ptBinsSoftAll_anti_stranges;
  }
  else if ( _flavTypeToComputeFor == allFlavors )
  {
    _ptBinsSoftAll = ptBinsSoftAll_all;
  }
  else
  {
    string errMsg = "error in printSoftPtSpectra, flavor not specified";
    throw eAnalysis_error( errMsg );
  }

  string type;
  if ( _flavTypeToComputeFor == gluon )
  {
    type = "gluon";
  }
  else if ( _flavTypeToComputeFor == light_quark || _flavTypeToComputeFor == anti_light_quark )
  {
    type = "quark";
  }
  if ( _flavTypeToComputeFor == up )
  {
    type = "up";
  }
  if ( _flavTypeToComputeFor == down )
  {
    type = "down";
  }
  if ( _flavTypeToComputeFor == strange )
  {
    type = "strange";
  }
  if ( _flavTypeToComputeFor == anti_up )
  {
    type = "anti_up";
  }
  if ( _flavTypeToComputeFor == anti_down )
  {
    type = "anti_down";
  }
  if ( _flavTypeToComputeFor == anti_strange )
  {
    type = "anti_strange";
  }
  else if ( _flavTypeToComputeFor == allFlavors )
  {
    type = "allFlavors";
  }

  string filename = filename_prefix + "_" + type + "_soft_spectra.f2";
  fstream file( filename.c_str(), ios::out | ios::trunc );

  //---- print header if file is empty ----
  file.seekp( 0, ios::end );
  long size = file.tellp();
  file.seekp( 0, ios::beg );
  if ( size == 0 )
    printHeader( file, ptSpectrum, end );
  //---------------------------------------


  //---------------------- y in [-inf,inf] ---------------------
  file << "#y in [-inf,inf]" << endl;
  for ( int i = 0;i < numberBinsSoftPT;i++ )
  {
    file << ptSoftBinLabels[i] << sep;
    for ( int j = 0; j <= nTimeSteps; j++ )
    {
      if ( j == 0 || j == nTimeSteps )
        file << _ptBinsSoftAll[j][i] << sep;
      else if ( tstep[j-1] <= theConfig->getRuntime() )
        file << _ptBinsSoftAll[j][i] << sep;
    }
    file << endl;
  }
  //------------------------------------------------------------

  file.close();
  //---------------------------------------------------------------

}




void analysis::onePartclCorrelations()
{
  const double eta_max = 1.0;
//  const double eta_max = 10000.0;//

  string filename;
  double pt_fin, pt_init, eta, dpt_vec, dpt_scal, rt_init, dE;


  filename = filename_prefix + "_ptIni_ptFin_1correl";
  binning2d ptIniptFin( filename, 0.0, 11.1, 50, 0.0, 11.1, 50 );

  filename = filename_prefix + "_ptIni_dpt_1correl";
  binning2d ptInidpt( filename, 0.0, 11.1, 50, 0.0, 11.1, 50 );

  filename = filename_prefix + "_posIni_dpt_1correl";
  binning2d posInidpt( filename, 0.0, 8.6, 50, 0.0, 10.1, 70 );

  filename = filename_prefix + "_dx_dpt_1correl";
  binning2d dxdpt( filename, 0.0, 6.1, 50, 0.0, 8.1, 70 );

  filename = filename_prefix + "_dx_dE_1correl";
  binning2d dxdE( filename, 0.0, 6.1, 50, -2.1, 8.1, 70 );

  filename = filename_prefix + "_dx_dE_precise_1correl";
  binning2d dxdE_precise( filename, 0.0, 6.1, 50, 0.0, 0.21, 50 );

  filename = filename_prefix + "_posIni_ptIni_ptFinAv_1correl";
  binningValues2d posIniptIniptFin( filename, 0.0, 8.1, 35, 0.0, 15.0, 50 );

  filename = filename_prefix + "_posIni_ptIni_dptAv_1correl";
  binningValues2d posIniptInidpt( filename, 0.0, 8.1, 35, 0.0, 15.0, 50 );

//   filename = filename_prefix + "_posIni_ptIni_dptptIAv_1correl";
  filename = "/dev/null";
  binningValues2d posIniptInidptptI( filename, 0.0, 8.1, 35, 0.0, 15.0, 50 );

  filename = filename_prefix + "_ptFin_posIni";
  binningValues ptFinposIni( filename, 0.0, 15.0, 30 );
  filename = filename_prefix + "_ptIni_posIni";
  binningValues ptIniposIni( filename, 0.0, 15.0, 30 );


//   filename = filename_prefix + "_ptFin_ptFin_test";
  filename = "/dev/null";
  binningValues2d ptFinptFin( filename, 0.0, 11.1, 50, 0.0, 11.1, 50 );



  double q_hat_scal_sum = 0.0, q_hat_vec_sum = 0.0;
  int count_q_hat_sum = 0;

  for ( unsigned int i = 0;i < addedParticles.size();i++ )
  {
//     if(addedParticles[i].T <= time)
//     {
    eta = addedParticles[i].Mom.Pseudorapidity( addedParticles[i].m );

    if ( fabs( eta ) <= eta_max )
    {
      pt_init = addedParticles[i].MomInit.Pt();
      pt_fin  = addedParticles[i].Mom.Pt();

      dpt_vec = (addedParticles[i].MomInit - addedParticles[i].Mom).Pt();

      dpt_scal = pt_init - pt_fin;
      rt_init  = addedParticles[i].PosInit.Pt();

      dE = addedParticles[i].MomInit.E() - addedParticles[i].Mom.E();

      ptIniptFin.add( pt_init, pt_fin );
      ptInidpt.add( pt_init, dpt_vec );
      posInidpt.add( rt_init, dpt_vec );
      dxdpt.add( addedParticles[i].X_traveled, dpt_vec );
      dxdE.add( addedParticles[i].X_traveled, dE );
      dxdE_precise.add( addedParticles[i].X_traveled, dE );

      posIniptIniptFin.add( rt_init, pt_init, pt_fin );
      posIniptInidpt.add( rt_init, pt_init, dpt_scal );
      posIniptInidptptI.add( rt_init, pt_init, dpt_scal / pt_init );
      ptFinposIni.add( pt_fin, rt_init );
      ptIniposIni.add( pt_init, rt_init );

      ptFinptFin.add( pt_init, pt_fin, dpt_scal );

      if ( addedParticles[i].X_traveled > 0.0 )
      {
        q_hat_scal_sum += pow( dpt_scal, 2.0 ) / addedParticles[i].X_traveled;
        q_hat_vec_sum += pow( dpt_vec, 2.0 ) / addedParticles[i].X_traveled;
        count_q_hat_sum++;
      }
    }
//     }
  }

  ptIniptFin.print();
  ptInidpt.print();
  posInidpt.print();
  dxdpt.print();
  dxdE.print();
  dxdE_precise.print();
  posIniptIniptFin.print();
  posIniptInidpt.print();
  posIniptInidptptI.print();
  ptFinposIni.print();
  ptIniposIni.print();

  ptFinptFin.print();

  filename = filename_prefix + "_q_hat";
  fstream print_q_hat( filename.c_str(), ios::out | ios::trunc );

//   cout << "Final output: q_hat = " << q_hat_scal_sum/count_q_hat_sum << endl;
  print_q_hat << q_hat_scal_sum / count_q_hat_sum << endl << q_hat_vec_sum / count_q_hat_sum << endl;

}

void analysis::twoPartclCorrelations()
{
  const double eta_max = 1.0;
//  const double eta_max = 10000.0; //!

  const double pt_min_trig = 4.0; //GeV, cut on trigger particle, which is considered for angle correlations dphi
  const double pt_min_assoc = 2.0; //GeV, cut on associated particle, which is considered for angle correlations dphi
  const double pt_min_dNdr = 5.0;//GeV, cut on particle for dN/dr_t

  string filename;
  double eta, rt_init;
  double dpt_ini, dpt_fin, dp_ini, dp_fin, dphi_fin, dphi_ini, scal_prod, length1, length2, pt_fin_i, pt_fin_j;

  filename = filename_prefix + "_dpt_partons";
  fstream partdpt( filename.c_str(), ios::out | ios::app );


  filename = filename_prefix + "_dptIni_dptFin_2correl";
  binning2d dptIdptF2( filename, 0.0, 8.1, 60, 0.0, 10.5, 70 );

  filename = filename_prefix + "_dptIni_dptFin_posIniAv_2correl";
  binningValues2d dptIdptFposI2( filename, 0.0, 8.1, 60, 0.0, 10.5, 70 );

//   filename = filename_prefix + "_dpIni_dpFin_2correl";
//   binning2d dpIdpF2(filename, 0.0, 10.5, 40, 0.0, 10.5, 40);

  filename = filename_prefix + "_posIni_dptFin_2correl";
  binning2d posIdptF2( filename, 0.0, 8.1, 35, 0.0, 10.5, 80 );

  filename = filename_prefix + "_posIni_dptIni_2correl";
  binning2d posIdptI2( filename, 0.0, 8.1, 35, 0.0, 10.5, 80 );



  filename = filename_prefix + "_dphiFin_4_2_2correl";
  binning dphiF42( filename, 0.0, M_PI, 80 );
  filename = filename_prefix + "_dphiIni_4_2_2correl";
  binning dphiI42( filename, 0.0, M_PI, 80 );
  filename = filename_prefix + "_dphiFin_2_2_2correl";
  binning dphiF22( filename, 0.0, M_PI, 80 );
  filename = filename_prefix + "_dphiIni_2_2_2correl";
  binning dphiI22( filename, 0.0, M_PI, 80 );
  filename = filename_prefix + "_dphiFin_all_2correl";
  binning dphiFall( filename, 0.0, M_PI, 80 );
  filename = filename_prefix + "_dphiIni_all_2correl";
  binning dphiIall( filename, 0.0, M_PI, 80 );
  filename = filename_prefix + "_dphiFin_rt3_2correl";
  binning dphiFrt3( filename, 0.0, M_PI, 80 );
  filename = filename_prefix + "_dphiIni_rt3_2correl";
  binning dphiIrt3( filename, 0.0, M_PI, 80 );

  filename = filename_prefix + "_dphiFin_posIni_all_2correl";
  binning2d dphiFposIniall( filename, 0.0, M_PI, 50, 0.0, 8.1, 35 );
  filename = filename_prefix + "_dphiIni_posIni_all_2correl";
  binning2d dphiIposIniall( filename, 0.0, M_PI, 50, 0.0, 8.1, 35 );

  filename = filename_prefix + "_dphi_iso";
  binning dphiiso( filename, 0.0, M_PI, 80 );

  filename = filename_prefix + "_dN_drInit_all";
  binning dNdrI_all( filename, 0.0, 15.0, 60 );
  filename = filename_prefix + "_dN_drInit_pt5";
  binning dNdrI_pt5( filename, 0.0, 15.0, 60 );
  filename = filename_prefix + "_dN_drInit_ptbar5";
  binning dNdrI_ptbar5( filename, 0.0, 15.0, 60 );
  filename = filename_prefix + "_dN_drInit_ptbar5_dpt02";
  binning dNdrI_ptbar5_dpt02( filename, 0.0, 15.0, 60 );

  // 2 particle correlations, charm and anti-charm which are produced in same reaction
  for ( unsigned int i = 0;i < addedParticles.size();i++ )
  {
//     if(addedParticles[i].T <= time)
//     {
    eta = addedParticles[i].Mom.Pseudorapidity( addedParticles[i].m );

    if ( fabs( eta ) <= eta_max )
    {
      for ( unsigned int j = i + 1;j < addedParticles.size();j++ )
      {
//         if ( addedParticles[i].N_EVENT_pp == addedParticles[j].N_EVENT_pp )
        if ( true )
        {

          //           if(addedParticles[i].T <= time)
          //           {
          eta = addedParticles[j].Mom.Pseudorapidity( addedParticles[j].m );

          if ( fabs( eta ) <= eta_max )
          {
//                 dpt_ini = sqrt( pow(addedParticles[i].MomInit.Px()-addedParticles[j].MomInit.Px(),2.0) + pow(addedParticles[i].MomInit.Py()-addedParticles[j].MomInit.Py(),2.0) );
//                 dpt_fin = sqrt( pow(addedParticles[i].Mom.Px()-addedParticles[j].Mom.Px(),2.0) + pow(addedParticles[i].Mom.Py()-addedParticles[j].Mom.Py(),2.0) );
            dpt_ini = fabs( addedParticles[i].MomInit.Pt() - addedParticles[j].MomInit.Pt() );
            dpt_fin = fabs( addedParticles[i].Mom.Pt() - addedParticles[j].Mom.Pt() );

            dp_ini = sqrt( (addedParticles[i].MomInit - addedParticles[j].MomInit).vec2() );
            dp_fin = sqrt( (addedParticles[i].Mom - addedParticles[j].Mom).vec2() );

            rt_init =  addedParticles[i].PosInit.Pt(); // same for both particles

            dptIdptF2.add( dpt_ini, dpt_fin );
            dptIdptFposI2.add( dpt_ini, dpt_fin, rt_init );
//                 dpIdpF2.add(dp_ini,dp_fin);
            posIdptF2.add( rt_init, dpt_fin );
            posIdptI2.add( rt_init, dpt_ini );



//                 if(dpt_fin > 7.0)
//                 {
//                   partdpt.width(15);
//                   partdpt << dpt_fin;
//                   partdpt.width(15);
//                   partdpt << fabs( sqrt( pow(addedParticles[i].Mom.Px(),2.0) + pow(addedParticles[i].Mom.Py(),2.0) ) - sqrt( pow(addedParticles[j].Mom.Px(),2.0) + pow(addedParticles[j].Mom.Py(),2.0) ) );
//
//                   partdpt.width(25);
//                   partdpt << addedParticles[i].Mom.Px();
//                   partdpt.width(15);
//                   partdpt << addedParticles[i].Mom.Py();
//                   partdpt.width(15);
//                   partdpt << addedParticles[i].Mom.Pz();
//
//                   partdpt.width(25);
//                   partdpt << addedParticles[j].Mom.Px();
//                   partdpt.width(15);
//                   partdpt << addedParticles[j].Mom.Py();
//                   partdpt.width(15);
//                   partdpt << addedParticles[j].Mom.Pz();
//
//                   partdpt.width(25);
//                   partdpt << addedParticles[i].MomInit.Px();
//                   partdpt.width(15);
//                   partdpt << addedParticles[i].MomInit.Py();
//                   partdpt.width(15);
//                   partdpt << addedParticles[i].MomInit.Pz();
//
//                   partdpt.width(25);
//                   partdpt << addedParticles[j].MomInit.Px();
//                   partdpt.width(15);
//                   partdpt << addedParticles[j].MomInit.Py();
//                   partdpt.width(15);
//                   partdpt << addedParticles[j].MomInit.Pz();
//
//                   partdpt.width(25);
//                   partdpt << addedParticles[j].PosInit.X();
//                   partdpt.width(15);
//                   partdpt << addedParticles[j].PosInit.Y();
//                   partdpt.width(15);
//                   partdpt << addedParticles[j].PosInit.Z();
//
//                   partdpt.width(25);
//                   partdpt << rt_init << endl;
//
//                 }


            // azimuthal angle between two charm quarks
            dphi_ini = acos( CosPhi( addedParticles[i].MomInit, addedParticles[j].MomInit ) );
            dphi_fin = acos( CosPhi( addedParticles[i].Mom, addedParticles[j].Mom ) );

            // total angle between two charm quarks
//                 scal_prod = addedParticles[i].PX_init*addedParticles[j].PX_init + addedParticles[i].PY_init*addedParticles[j].PY_init + addedParticles[i].PZ_init*addedParticles[j].PZ_init;
//                 length1 = sqrt( pow(addedParticles[i].PX_init,2.0) + pow(addedParticles[i].PY_init,2.0) + pow(addedParticles[i].PZ_init,2.0)  );
//                 length2 = sqrt( pow(addedParticles[j].PX_init,2.0) + pow(addedParticles[j].PY_init,2.0) + pow(addedParticles[j].PZ_init,2.0)  );
//                 dphi_ini = acos( scal_prod / length1 / length2 );
//
//                 scal_prod = addedParticles[i].PX*addedParticles[j].PX + addedParticles[i].PY*addedParticles[j].PY + addedParticles[i].PZ*addedParticles[j].PZ;
//                 length1 = sqrt( pow(addedParticles[i].PX,2.0) + pow(addedParticles[i].PY,2.0) + pow(addedParticles[i].PZ,2.0)  );
//                 length2 = sqrt( pow(addedParticles[j].PX,2.0) + pow(addedParticles[j].PY,2.0) + pow(addedParticles[j].PZ,2.0)  );
//                 dphi_fin = acos( scal_prod / length1 / length2 );

//                 dphi_fin = cos(dphi_fin);
//                 dphi_ini = cos(dphi_ini);

            dphiFall.add( dphi_fin );
            dphiIall.add( dphi_ini );

            if ( addedParticles[i].Mom.Pt2() >= pow( pt_min_assoc, 2.0 ) && addedParticles[j].Mom.Pt2() >= pow( pt_min_assoc, 2.0 ) )
            {
              dphiF22.add( dphi_fin );
              dphiI22.add( dphi_ini );

              if (( addedParticles[i].Mom.Pt2() >= pow( pt_min_trig, 2.0 ) && addedParticles[j].Mom.Pt2() >= pow( pt_min_assoc, 2.0 )) ||
                  ( addedParticles[i].Mom.Pt2() >= pow( pt_min_assoc, 2.0 ) && addedParticles[j].Mom.Pt2() >= pow( pt_min_trig, 2.0 )))
              {
                dphiF42.add( dphi_fin );
                dphiI42.add( dphi_ini );
              }
            }

            if ( rt_init )
            {
              dphiFrt3.add( dphi_fin );
              dphiIrt3.add( dphi_ini );
            }

            dphiFposIniall.add( dphi_fin, rt_init );
            dphiIposIniall.add( dphi_ini, rt_init );



            pt_fin_i = addedParticles[i].Mom.Pt();
            pt_fin_j = addedParticles[j].Mom.Pt();

            dNdrI_all.add( rt_init );
            if ( pt_fin_i > pt_min_dNdr )
              dNdrI_pt5.add( rt_init );
            if (( pt_fin_i + pt_fin_j ) / 2.0 > pt_min_dNdr )
              dNdrI_ptbar5.add( rt_init );
            if ((( pt_fin_i + pt_fin_j ) / 2.0 > pt_min_dNdr ) && ( dpt_fin / (( pt_fin_i + pt_fin_j ) / 2.0 ) < 0.2 ) )
              dNdrI_ptbar5_dpt02.add( rt_init );

          }
          //           }
          break;
        }
      }

    }
//     }
  }



  for ( int j = 0;j < addedParticles.size();j++ )
  {
    eta = addedParticles[j].Mom.Pseudorapidity( addedParticles[j].m );

    VectorXYZ V = VectorXYZ( 0.0, 1.0, 0.0 );
    if ( fabs( eta ) <= eta_max )
    {
      dphi_fin = acos( CosPhi(addedParticles[j].Mom, V) );
      dphiiso.add( dphi_fin );
    }
  }

  dptIdptF2.print();
  dptIdptFposI2.print();
//   dpIdpF2.print();
  posIdptF2.print();
  posIdptI2.print();
  dphiF42.print();
  dphiI42.print();
  dphiF22.print();
  dphiI22.print();
  dphiFall.print();
  dphiIall.print();
  dphiiso.print();
  dphiFposIniall.print();
  dphiIposIniall.print();
  dphiFrt3.print();
  dphiIrt3.print();

  dNdrI_all.print();
  dNdrI_pt5.print();
  dNdrI_ptbar5.print();
  dNdrI_ptbar5_dpt02.print();

}




// compute v2 of gluons,quarks,photons and charm quarks
void analysis::computeV2RAA( string name, const double _outputTime )
{
  v2RAA theV2RAA( theConfig, name, filename_prefix, rapidityRanges );
  
  const double pt_min_v2RAA = theConfig->getMinimumPT();
  double pt_max_v2RAA, nbins_v2RAA;
  if( pt_min_v2RAA < 35 )
  {
    pt_max_v2RAA = 55.0;
    nbins_v2RAA = 67; // good for pt_min = 0
    if( FPT_COMP_GE( pt_min_v2RAA, 6.0 ) )
      nbins_v2RAA = 25;
    else if( FPT_COMP_GE( pt_min_v2RAA, 3.0 ) )
      nbins_v2RAA = 34;
  }
  else
  {
    pt_max_v2RAA = 150.0;
    nbins_v2RAA = 60;
    
    if( FPT_COMP_GE( pt_min_v2RAA, 70.0 ) )
      nbins_v2RAA = 52;
    else if( FPT_COMP_GE( pt_min_v2RAA, 90.0 ) )
      nbins_v2RAA = 45;
  }
  if (theConfig->doOutput_QCDparticles())
  {
    if( theConfig->isStudyNonPromptJpsiInsteadOfElectrons() )
    {
      theV2RAA.setPtBinProperties( pt_min_v2RAA, pt_max_v2RAA, nbins_v2RAA );
      
      theV2RAA.computeFor( bottom, addedParticles, addedParticles.size(), "added", _outputTime, v2jets );
      if( mesonDecay )
        theV2RAA.computeFor( electron_gen, addedPartcl_electron, addedPartcl_electron.size(), "added", _outputTime, v2jets );
    }
    if( studyHQ || studyJpsi )
    {
      theV2RAA.setPtBinProperties( pt_min_v2RAA, pt_max_v2RAA, nbins_v2RAA );
      
      theV2RAA.computeFor( charm, addedParticles, addedParticles.size(), "added", _outputTime, v2jets );
      theV2RAA.computeFor( bottom, addedParticles, addedParticles.size(), "added", _outputTime, v2jets );
      theV2RAA.computeFor( heavy_quark, addedParticles, addedParticles.size(), "added", _outputTime, v2jets );

      if ( name == "initial" || name == "final" )
      {
        if( hadronization_hq )
        {
          theV2RAA.computeFor( dmeson_gen, addedParticlesCopy, addedParticlesCopy.size(), "added", _outputTime, v2jets );
          theV2RAA.computeFor( bmeson_gen, addedParticlesCopy, addedParticlesCopy.size(), "added", _outputTime, v2jets );
        }
        
        if( mesonDecay )
        {
          theV2RAA.computeFor( electron_gen, addedPartcl_electron, addedPartcl_electron.size(), "added", _outputTime, v2jets );
          theV2RAA.computeFor( c_electron, addedPartcl_electron, addedPartcl_electron.size(), "added", _outputTime, v2jets );  // electrons from charm quarks (D mesons actually)
          theV2RAA.computeFor( b_electron, addedPartcl_electron, addedPartcl_electron.size(), "added", _outputTime, v2jets );  // electrons from bottom quarks (B mesons actually)
        }
      }
      
      // also take a look at light parton v2 of background
      theV2RAA.computeFor( gluon, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
      theV2RAA.computeFor( light_quark, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
      
      if ( studyJpsi )
      {
        theV2RAA.setPtBinProperties( 0.0, 15.0, 20 );
        
        theV2RAA.computeFor( jpsi, addedParticles, addedParticles.size(), "added", _outputTime, v2jets );
        theV2RAA.computeFor( jpsi_ini, addedParticles, addedParticles.size(), "added", _outputTime, v2jets );
        theV2RAA.computeFor( jpsi_sec, addedParticles, addedParticles.size(), "added", _outputTime, v2jets );
      }
    }
    else if ( theConfig->getPtCutoff() > 80.0 )
    {
      theV2RAA.setPtBinProperties( 0.8*theConfig->getPtCutoff(), 3.0*theConfig->getPtCutoff(), static_cast< int >( 2.2*theConfig->getPtCutoff() ) );
      
      theV2RAA.computeFor( gluon, addedParticles, addedParticles.size(), "jets", _outputTime, v2jets );
      theV2RAA.computeFor( light_quark, addedParticles, addedParticles.size(), "jets", _outputTime, v2jets );
    }
    else 
    {
      theV2RAA.setPtBinProperties( pt_min_v2RAA, pt_max_v2RAA, nbins_v2RAA );
      theV2RAA.computeFor( gluon, addedParticles, addedParticles.size(), "jets", _outputTime, v2jets );
      theV2RAA.computeFor( light_quark, addedParticles, addedParticles.size(), "jets", _outputTime, v2jets );
      theV2RAA.computeFor( up, addedParticles, addedParticles.size(), "jets", _outputTime, v2jets );
      theV2RAA.computeFor( down, addedParticles, addedParticles.size(), "jets", _outputTime, v2jets );
      theV2RAA.computeFor( strange, addedParticles, addedParticles.size(), "jets", _outputTime, v2jets );
      theV2RAA.computeFor( anti_up, addedParticles, addedParticles.size(), "jets", _outputTime, v2jets );
      theV2RAA.computeFor( anti_down, addedParticles, addedParticles.size(), "jets", _outputTime, v2jets );
      theV2RAA.computeFor( anti_strange, addedParticles, addedParticles.size(), "jets", _outputTime, v2jets );
    }
    
    if( studyBackground )
    {
      theV2RAA.computeFor( gluon, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
      theV2RAA.computeFor( up, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
      theV2RAA.computeFor( down, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
      theV2RAA.computeFor( strange, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
      theV2RAA.computeFor( anti_up, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
      theV2RAA.computeFor( anti_down, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
      theV2RAA.computeFor( anti_strange, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
      theV2RAA.computeFor( light_quark, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
    }
  }
  if ( studyPartons )
  {
    cout << "Analyse partons." << endl;
    theV2RAA.setPtBinProperties( 0.0, 6.0, 8, 0.0, 6.0, 8 );
    theV2RAA.computeFor( light_quark, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
    theV2RAA.computeFor( gluon, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
    theV2RAA.computeFor( allFlavors, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
  }
  if( studyPhotons )
  {
    cout << "analyse Photons. Size of Photon vector: " << noninteractingParticles.size() <<  endl;
    //if(theConfig->v2_bigger > theConfig->v2_smaller) cout << "bigger! " << theConfig->v2_bigger - theConfig->v2_smaller << endl;
    //if(theConfig->v2_bigger < theConfig->v2_smaller) cout << "smaller! " << -theConfig->v2_bigger + theConfig->v2_smaller << endl; 
    //cout << "# Number of colliding pairs with average positive v2: " << theConfig->countPositiveV2 << endl;
    //cout << "# Number of colliding pairs with average negative v2: " << theConfig->countNegativeV2 << endl;     
    theV2RAA.setPtBinProperties( 0.0, 6.0, 70, 0.0, 6.0, 70 );
    theV2RAA.computeFor( photon, noninteractingParticles, noninteractingParticles.size(), "LO", _outputTime, v2background );
    //theV2RAA.computeFor( light_quark, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
    //theV2RAA.computeFor( gluon, particles_atTimeNow, particles_atTimeNow.size(), "background", _outputTime, v2background );
  }
  
  if( studyDileptons )
  {
    cout << "analyse Dileptons. Size of Dilepton pair vector: " << dileptons.size() <<  endl;
    computeAndPrintDileptonSpectra(name, _outputTime );
  }

}

void analysis::computeAndPrintDileptonSpectra(string name, const double _outputTime )
{
  binningInRanges DileptonInvMassSeveralPT;
  binningInRanges DileptonPTSeveralInvMass;  
  binningInRanges DileptonInvMassSeveralTimes;
  binningInRanges DileptonInvMassCuts;
  
  double pGEV, invMassGeV, productionTimeFm;
  
  DileptonInvMassSeveralPT.setNumberRanges(5);
  DileptonInvMassSeveralPT.setMinMaxN(0.0, 4.0, 50);
  DileptonInvMassSeveralPT.setRangevalues(0,0.0,10000.0);//in pt full inclusive spectrum
  DileptonInvMassSeveralPT.setRangevalues(1,0.0,1.0);
  DileptonInvMassSeveralPT.setRangevalues(2,1.0,2.0);
  DileptonInvMassSeveralPT.setRangevalues(3,2.0,10.0);
  DileptonInvMassSeveralPT.setRangevalues(4,0.0,3.0);
  
  DileptonInvMassSeveralTimes.setNumberRanges(5);//in pt full inclusive spectrum
  DileptonInvMassSeveralTimes.setMinMaxN(0.0, 4.0, 50);
  DileptonInvMassSeveralTimes.setRangevalues(0,0.5,10.0); // Timeslots 0.5-10fm
  DileptonInvMassSeveralTimes.setRangevalues(1,1.0,10.0); // Timeslots 1-10fm
  DileptonInvMassSeveralTimes.setRangevalues(2,2.0,10.0); // Timeslots 2-10fm
  DileptonInvMassSeveralTimes.setRangevalues(3,3.0,10.0); // Timeslots 3-10fm
  DileptonInvMassSeveralTimes.setRangevalues(4,0.0,10.0); // Full inclusive

  DileptonPTSeveralInvMass.setNumberRanges(5);
  DileptonPTSeveralInvMass.setMinMaxN(0.0, 4.0, 50);
  DileptonPTSeveralInvMass.setRangevalues(0,0.0,10000.0);//full inclusive spectrum
  DileptonPTSeveralInvMass.setRangevalues(1,0.0,1.0);
  DileptonPTSeveralInvMass.setRangevalues(2,1.0,2.0);
  DileptonPTSeveralInvMass.setRangevalues(3,2.0,10.0);
  DileptonPTSeveralInvMass.setRangevalues(4,0.0,3.0);
  
  DileptonInvMassCuts.setNumberRanges(2);
  DileptonInvMassCuts.setMinMaxN(0.0, 4.0, 50);
  DileptonInvMassCuts.setRangevalues(0,0.2,1000.0); // PHENIX Run 10 acceptance
  DileptonInvMassCuts.setRangevalues(1,0.0,1000.0); // PHENIX Run 10 acceptance
  
  for ( unsigned int j = 0; j < dileptons.size(); j++ )
  {
    pGEV = sqrt(dileptons[j].Mom.vec2());
    invMassGeV = dileptons[j].m;
    productionTimeFm = dileptons[j].production_time;
    
    DileptonInvMassSeveralPT.add(pGEV,invMassGeV);
    DileptonInvMassSeveralTimes.add(productionTimeFm,invMassGeV);
    DileptonPTSeveralInvMass.add(invMassGeV,pGEV);
    
    //Experimental Cuts
    // y-cut
    double SingleElectronYCut = 0.35; 
    if (abs(dileptons[j].dilepton_y_min)<SingleElectronYCut && abs(dileptons[j].dilepton_y_max)<SingleElectronYCut )
    {
      DileptonInvMassCuts.add(dileptons[j].dilepton_pt_min,invMassGeV);
    }
 
  }
  
  string filename_dRdMTimeslots = filename_prefix + "_" + "dilepton" + "_dR_dM_Time" + name;  
  fstream file3( filename_dRdMTimeslots, ios::app | ios::out );
  //print header if file is empty
  file3.seekp( 0, ios::end );
  long size3 = file3.tellp();
  file3.seekp( 0, ios::beg );
  if ( size3 == 0 )
    printHeader( file3, dileptondRdMTimeDist, _outputTime );
  DileptonInvMassSeveralTimes.print(file3);
  
  string filename_dRdM = filename_prefix + "_" + "dilepton" + "_dR_dM_" + name;  
  fstream file( filename_dRdM, ios::app | ios::out );
  //print header if file is empty
  file.seekp( 0, ios::end );
  long size = file.tellp();
  file.seekp( 0, ios::beg );
  if ( size == 0 )
    printHeader( file, dileptondRdMDist, _outputTime );
  DileptonInvMassSeveralPT.print(file);

  string filename_dRdP = filename_prefix + "_" + "dilepton" + "_dR_dP_" + name;    
  fstream file2( filename_dRdP, ios::app | ios::out );
  //print header if file is empty
  file2.seekp( 0, ios::end );
  long size2 = file2.tellp();
  file2.seekp( 0, ios::beg );
  if ( size2 == 0 )
    printHeader( file2, dileptondRdPDist, _outputTime );
  DileptonPTSeveralInvMass.print(file2);
  
  string filename_dNdMCuts = filename_prefix + "_" + "dileptonCuts" + "_dN_dM_" + name;    
  fstream file4( filename_dNdMCuts, ios::app | ios::out );
  //print header if file is empty
  file4.seekp( 0, ios::end );
  long size4 = file4.tellp();
  file4.seekp( 0, ios::beg );
  if ( size4 == 0 )
    printHeader( file4, dileptondNdMCuts, _outputTime );
  DileptonInvMassCuts.print(file4);

}         
          
          



v2RAA::v2RAA( config * const c, string name_arg, string filename_prefix_arg, std::vector<analysisRapidityRange> rapidityRanges_arg, const double pt_min_arg, const double pt_max_arg, const int n_g_arg, const double pt_min_background_arg, const double pt_max_background_arg, const int n_g_background_arg ):
    theConfig( c ), name( name_arg ), filename_prefix( filename_prefix_arg ), rapidityRanges( rapidityRanges_arg ), pt_min( pt_min_arg ), pt_max( pt_max_arg ), n_g( n_g_arg ), pt_min_background( pt_min_arg ), pt_max_background( pt_max_arg ), n_g_background( n_g_arg )
{
  eta_bins = rapidityRanges.size();
  
  //Photon Angle Bins config
  PhotonNumberVsAngleBin.setMinMaxN( 0.0 , 90, 100 );
  
  studyInitialCutOffEffect = true;
  lower_time_cutoff_for_v2 = theConfig->getAnalysisPhotonsTimeCut();
  lower_pt_cutoff_for_v2 = theConfig->getAnalysisPhotonsPTCut();
}




void v2RAA::computeFor( const FLAVOR_TYPE _flavTypeToComputeFor, vector<ParticleOffline>& _particles, const int n_particles, string additionalNameTag, const double _outputTime, const v2Type _v2type )
{
  double eta, pp, pt, v2, xt;
  double sinAlpha, alpha;
  int dummy, flavor, n_bins;
  int alphaIndex;
  int numberOfPhotonV2PtBins = theConfig->getAnalysisPhotonsNBinsV2();
  double _pt_min, _pt_max;

  string filename_v2, filename_v2_summed, filename_v2_tot, filename_yield, filename_pt_angleDependence, type, filename_dRdM, filename_dRdP;
  
  type = Particle::getName( _flavTypeToComputeFor );
  
  if ( _v2type == v2background )
  {
    n_bins = n_g_background;
    _pt_max = pt_max_background;
    _pt_min = pt_min_background;
  }
  else
  {
    n_bins = n_g;
    _pt_max = pt_max;
    _pt_min = pt_min;
  }
  
  // avoid problem with binning pt logaritmitically: cannot deal with pt = 0
  if( _pt_min < 0.1 )
    _pt_min = 0.1;
  
  const double d_ln_pt = ( log( _pt_max ) - log( _pt_min ) ) / n_bins;
  double d_ln_pt_v2=0.;
  double _pt_min_v2=0.;
  double _pt_max_v2=0.;
  
  double v2sum[eta_bins];
  int NmbInRange[eta_bins];
  double v2sumInitialCutOff[eta_bins];
  int NmbInRangeInitialCutOff[eta_bins];
  double v2sumPtCutOff[eta_bins];
  int NmbInRangePtCutOff[eta_bins];   
  int NmbInnerRegion = 0;
  for ( int j = 0;j < eta_bins;j++ )
  {
    v2sum[j] = 0.0;
    NmbInRange[j] = 0;
    v2sumInitialCutOff[j] = 0.0;
    NmbInRangeInitialCutOff[j] = 0;
    v2sumPtCutOff[j] = 0.0;
    NmbInRangePtCutOff[j] = 0;       
  }

  int v2pt_binnumber;
  v2pt_binnumber = 10;   
  _pt_min_v2 = 0.1;
  _pt_max_v2 = 3.;
  d_ln_pt_v2 = ( log( _pt_max_v2 ) - log( _pt_min_v2 ) ) / v2pt_binnumber;
  
  if(_flavTypeToComputeFor==photon)
  {
    v2pt_binnumber = numberOfPhotonV2PtBins; //Paper 2016: 8  
    _pt_min_v2 = 0.1;
    _pt_max_v2 = 6.0;
    d_ln_pt_v2 = ( log( _pt_max_v2 ) - log( _pt_min_v2 ) ) / v2pt_binnumber;
  }
  else
  {
    v2pt_binnumber =n_bins;
  }
  double ptBinsV2[eta_bins][v2pt_binnumber+1];
  double ptBinsV2InitialCutOff[eta_bins][v2pt_binnumber+1];  
  int ptBinsNmbV2[eta_bins][v2pt_binnumber+1];
  int ptBinsNmbInitialCutOffV2[eta_bins][v2pt_binnumber+1];
  
  int ptBinsNmb[eta_bins][n_bins+1];
  int ptBinsNmbCompton[eta_bins][n_bins+1];
  int ptBinsNmbAnnihilation[eta_bins][n_bins+1];
  
  int ptBinsNmbDileptonsdRdM[eta_bins][n_bins+1];
  int ptBinsNmbDileptonsdRdP[eta_bins][n_bins+1];
  
  int ptBinsNmbInitialCutOff[eta_bins][n_bins+1];  
  int ptBinsInnerRegion[n_bins+1];
  
  for ( int j = 0;j < n_bins + 1;j++ )
  {
    ptBinsInnerRegion[j] = 0;
    for ( int i = 0;i < eta_bins;i++ )
    {
      ptBinsNmb[i][j] = 0.0;
      ptBinsNmbCompton[i][j] = 0.0;
      ptBinsNmbAnnihilation[i][j] = 0.0;
      ptBinsNmbInitialCutOff[i][j] = 0.0;
      ptBinsNmbDileptonsdRdM[i][j] = 0.0;
      ptBinsNmbDileptonsdRdP[i][j] = 0.0;
    }
  }
  for ( int j = 0;j < v2pt_binnumber + 1;j++ )
  {
    for ( int i = 0;i < eta_bins;i++ )
    {
      ptBinsV2[i][j] = 0.0;
      ptBinsNmbV2[i][j] = 0.0;

      ptBinsV2InitialCutOff[i][j] = 0.0;
      ptBinsNmbInitialCutOffV2[i][j] =0.0;
    }
  }
  
  
  
  const double deltaAlpha = 30; // degrees
  const int nAlphaBins = 3;  // 90° / 15°
  double ptBinsAngleDep[eta_bins][nAlphaBins][n_bins+1];
  for ( int i = 0; i < eta_bins; i++ )
  {
    for ( int j = 0; j < nAlphaBins; j++ )
    {
      for ( int k = 0; k < n_bins + 1; k++ )
      {
        ptBinsAngleDep[i][j][k] = 0;
      }
    }
  }
  
  

  
  
  // compute angle and bin it
  double angle;
  for ( int i = 0; i < n_particles; i++ )
  {
   //angle = acos(_particles[i].Mom.Px()/_particles[i].Mom.E())
   //_particles[i].Mom.Phi_x()*180.0/M_PI;
   angle = acos(fabs(_particles[i].Mom.Px())/_particles[i].Mom.Pt())*180.0/M_PI;
   PhotonNumberVsAngleBin.add(angle);
  }  
  //print angle distribution
  string filename_angle = filename_prefix + "_" + type + "_" + additionalNameTag + "_angle_distribution_" + name;
  fstream print_angle_distribution( filename_angle.c_str(), ios::out | ios::trunc );

  print_angle_distribution << "# Number of particles emitted vs angle to x-Axis " << type << endl;
  print_angle_distribution << "# t = " << _outputTime << endl;
  print_angle_distribution << "# angle| \t number of photons| \t number of photons/binWidth| ..."  << endl;
 
  PhotonNumberVsAngleBin.print(print_angle_distribution);
  
    
  // compute v2 and bin it into pt bins
  for ( int i = 0; i < n_particles; i++ )
  {
    pt = _particles[i].Mom.Pt();
    xt = _particles[i].Pos.Pt();

    sinAlpha = _particles[i].Mom.Py() / pt;
    alpha = asin( fabs( sinAlpha ) );
    alpha = alpha * 180 / M_PI;
    
    alphaIndex = static_cast<int>( alpha / deltaAlpha );
    if ( alphaIndex >= nAlphaBins )
    {
      alphaIndex = nAlphaBins - 1;
    }
        
    // for most particle species we are interested in the pseudorapidity (for massless particle it does not matter anyhow)
    eta = _particles[i].Mom.Pseudorapidity(_particles[i].m);
    
    // for some scenarios however explicitly the rapidity is measured. So substitute eta by the rapidity:
    if( ( theConfig->isStudyNonPromptJpsiInsteadOfElectrons() &&
          ( _flavTypeToComputeFor == charm ||
            _flavTypeToComputeFor == bottom ||
            _flavTypeToComputeFor == heavy_quark ||
            _flavTypeToComputeFor == c_electron ||
            _flavTypeToComputeFor == b_electron ||
            _flavTypeToComputeFor == electron_gen 
          ) ) ||
        ParticleOffline::mapToGenericFlavorType( _flavTypeToComputeFor ) == dmeson_gen ||
        ParticleOffline::mapToGenericFlavorType( _flavTypeToComputeFor ) == bmeson_gen ||
        ParticleOffline::mapToGenericFlavorType( _flavTypeToComputeFor ) == jpsi
    )
      eta = _particles[i].Mom.Rapidity();

    v2 = ( pow( _particles[i].Mom.Px(), 2.0 ) - pow( _particles[i].Mom.Py(), 2.0 ) ) / pow( pt, 2.0 );

    flavor = _particles[i].FLAVOR;
    
    FLAVOR_TYPE mother_flav;
    if( _flavTypeToComputeFor == c_electron || _flavTypeToComputeFor == b_electron )
    {
      int mother_id = i / theConfig->getNumberElectronStat();
      mother_flav = addedParticlesCopy[mother_id].FLAVOR;
    }

    FLAVOR_TYPE genFlavor = ParticleOffline::mapToGenericFlavorType( static_cast<FLAVOR_TYPE>( flavor ) );
    
    if( ( _flavTypeToComputeFor == flavor ) || 
          ( _flavTypeToComputeFor == allFlavors ) ||
          ( _flavTypeToComputeFor == light_quark && ( genFlavor == light_quark || genFlavor == anti_light_quark ) ) ||
          ( _flavTypeToComputeFor == genFlavor ) ||
          ( _flavTypeToComputeFor == charm && ( flavor == anti_charm ) ) ||
          ( _flavTypeToComputeFor == bottom && ( flavor == anti_bottom ) ) ||
          ( _flavTypeToComputeFor == heavy_quark && ( flavor == charm || flavor == bottom || flavor == anti_charm || flavor == anti_bottom ) ) ||
          ( _flavTypeToComputeFor == c_electron && ( ( flavor == electron || flavor == positron ) && ParticleOffline::mapToGenericFlavorType( static_cast<FLAVOR_TYPE>( mother_flav ) ) == dmeson_gen ) ) ||
          ( _flavTypeToComputeFor == b_electron && ( ( flavor == electron || flavor == positron ) && ParticleOffline::mapToGenericFlavorType( static_cast<FLAVOR_TYPE>( mother_flav ) ) == bmeson_gen ) ) ||
          ( _flavTypeToComputeFor == jpsi_ini && ( flavor == jpsi && _particles[i].initially_produced ) ) ||
          ( _flavTypeToComputeFor == jpsi_sec && ( flavor == jpsi && !_particles[i].initially_produced ) )
      ) 
    {
      // individually check for each rapidity range whether this particle needs to be binned
      for ( int yRangeIndex = 0; yRangeIndex < eta_bins; yRangeIndex++ )
      {
        if ( fabs( eta ) >= rapidityRanges[yRangeIndex].yleft && fabs( eta ) <= rapidityRanges[yRangeIndex].yright )
        {
          if (_particles[i].production_time > lower_time_cutoff_for_v2 )//TODO Include check if photon or not
          {
             v2sumInitialCutOff[yRangeIndex] += v2;
             NmbInRangeInitialCutOff[yRangeIndex]++;
          }
          if (pt > lower_pt_cutoff_for_v2 )//TODO Include check if photon or not
          {
             v2sumPtCutOff[yRangeIndex] += v2;
             NmbInRangePtCutOff[yRangeIndex]++;
          }
          v2sum[yRangeIndex] += v2;
          NmbInRange[yRangeIndex]++;
          
          // Differential V2, has different number of pt bins
          if ( pt <= _pt_max_v2 && pt > _pt_min_v2 )
          {
            dummy = int(( log( pt ) - log( _pt_min_v2 ) ) / d_ln_pt_v2 );
            if (_particles[i].production_time > lower_time_cutoff_for_v2 )
            {            
              ptBinsV2InitialCutOff[yRangeIndex][dummy] += v2;
              ptBinsNmbInitialCutOffV2[yRangeIndex][dummy]++;
            }           
            ptBinsV2[yRangeIndex][dummy] += v2;     
            ptBinsNmbV2[yRangeIndex][dummy]++;
          }
          
          
          // Yield
          if ( pt <= _pt_max && pt > _pt_min )
          {
            dummy = int(( log( pt ) - log( _pt_min ) ) / d_ln_pt );
            if (_particles[i].production_time > lower_time_cutoff_for_v2 )//TODO Include check if photon or not
            {            
              ptBinsNmbInitialCutOff[yRangeIndex][dummy]++;
            }

            ptBinsNmb[yRangeIndex][dummy]++;
            ptBinsAngleDep[yRangeIndex][alphaIndex][dummy]++;
            
            if(_particles[i].production_mechanism == OnlyCompton)
            {
              ptBinsNmbCompton[yRangeIndex][dummy]++;
            }
            if(_particles[i].production_mechanism == OnlyAnnihilation)
            {
              ptBinsNmbAnnihilation[yRangeIndex][dummy]++;
            }                       
          }          
        }
      }
    }
  }

  cout << "binning done" << endl;
  
  int binMax = 0;
  int binMin = n_particles;
  for ( int k = 0; k < n_bins + 1; k++ )
  {
    if ( ptBinsNmb[0][k] > binMax )
    {
      binMax = ptBinsNmb[0][k];
    }
    if ( ptBinsNmb[0][k] < binMin && ptBinsNmb[0][k] != 0 )
    {
      binMin = ptBinsNmb[0][k];
    }
  }

  // TEST:
  // Number of average v2 of colliding particle pairs 
  /*if(_flavTypeToComputeFor == photon )
  {
    string filename_average_initial_v2 = filename_prefix + "_" + type + "_" + additionalNameTag + "_collidingPairsInitialV2_" + name;
  
    //cout << "# Number of colliding pairs with average positive v2: " << theConfig->countPositiveV2 << endl;
    //cout << "# Number of colliding pairs with average negative v2: " << theConfig->countNegativeV2 << endl;  
    //print angle distribution
    
    fstream print_average_initial_v2( filename_average_initial_v2.c_str(), ios::out | ios::trunc );
    print_average_initial_v2 << "#Number of colliding pairs with positive <v2> | Number of colliding pairs with negative <v2>"  << endl;
    print_average_initial_v2 << theConfig->countPositiveV2 << '\t' << theConfig->countNegativeV2 << endl;   
  }*/
  
  
  
  
  double pt_out;

  filename_v2_summed = filename_prefix + "_" + type + "_" + additionalNameTag + "_v2_pt_summed_" + name;
  filename_v2_tot = filename_prefix + "_" + type + "_" + additionalNameTag + "_v2_tot_" + name;
  filename_yield = filename_prefix + "_" + type + "_" + additionalNameTag + "_yield_pt_" + name;

  filename_pt_angleDependence = filename_prefix + "_" + type + "_" + additionalNameTag + "_pt_angular_dependence_" + name;

  fstream print_v2_summed( filename_v2_summed.c_str(), ios::out | ios::trunc );
  fstream print_v2_tot( filename_v2_tot.c_str(), ios::out | ios::trunc );
  fstream print_yield( filename_yield.c_str(), ios::out | ios::trunc );
  fstream print_pt_angleDependence( filename_pt_angleDependence.c_str(), ios::out | ios::trunc );

//   cout << "total v2 of " << type << " in Y=+-0.35 = " << v2sum[0]/NmbInRange[0] << endl;

  // print total v2
  print_v2_tot << "# total v2 of " << type << endl;
  print_v2_tot << "# t = " << _outputTime << endl;
  print_v2_tot << "# bin statistics for 0.35 mid-rapidity:  Avg per bin=" << double( NmbInRange[0] ) / n_c << "   Min=" << binMin << "   Max=" << binMax << endl;
  //print_v2_tot << "# Number of colliding pairs with average positive v2: " << theConfig->countPositiveV2 << endl;
  //print_v2_tot << "# Number of colliding pairs with average negative v2: " << theConfig->countNegativeV2 << endl;  
  print_v2_tot << "#pt_min\t#total v2 for different rapidity bins" << endl;
  //This is the trigger info for Moritz Analysis scripts
  print_v2_tot << _pt_min;
  
  print_v2_tot.width( 15 );
  for ( int i = 0;i < eta_bins;i++ )
  {
    if ( NmbInRange[i] > 0 )
    {
      print_v2_tot <<  v2sum[i] / NmbInRange[i];
    }
    else
    {
      print_v2_tot <<  0;
    }
    print_v2_tot.width( 15 );
  }
  print_v2_tot << endl;
  print_v2_tot << "[#v2_sum\t#number in range], for different rapidity bins" << endl;
  
  for ( int i = 0;i < eta_bins;i++ )
  {
    print_v2_tot <<  v2sum[i];
    print_v2_tot.width( 15 );
    print_v2_tot <<  NmbInRange[i];
    print_v2_tot.width( 15 );
  }
  print_v2_tot << endl;

  /** WARNING This output is not necesary
  print_v2_tot << _pt_max;
  print_v2_tot.width( 15 );
  for ( int i = 0;i < eta_bins;i++ )
  {
    if ( NmbInRange[i] > 0 )
    {
      print_v2_tot <<  v2sum[i] / NmbInRange[i];
    }
    else
    {
      print_v2_tot <<  0;
    }
    print_v2_tot.width( 15 );
  }
  for ( int i = 0;i < eta_bins;i++ )
  {
    print_v2_tot <<  v2sum[i];
    print_v2_tot.width( 15 );
    print_v2_tot <<  NmbInRange[i];
    print_v2_tot.width( 15 );
  }
  print_v2_tot << endl;
  */
  
  print_v2_tot << "# Cut off initial timesteps! Only count from time " << lower_time_cutoff_for_v2 << " fm" << endl;
  print_v2_tot << "#1234\t#average integrated v2, for all rapidity bins" << endl; 
  print_v2_tot << 1234;
  print_v2_tot.width( 15 );
  
  for ( int i = 0;i < eta_bins;i++ )
  {
    if ( NmbInRange[i] > 0 )
    {
      print_v2_tot <<  v2sumInitialCutOff[i] / NmbInRangeInitialCutOff[i];
    }
    else
    {
      print_v2_tot <<  0;
    }
    print_v2_tot << '\t' ;
  } 
  print_v2_tot << endl;
  
  print_v2_tot << "[#sum\t#number], for all rapidity bins, for lower time cutoff" << endl; 
  for ( int i = 0;i < eta_bins;i++ )
  {
    print_v2_tot <<  v2sumInitialCutOff[i] << '\t';
    print_v2_tot <<  NmbInRangeInitialCutOff[i] << '\t';
  }
  print_v2_tot << endl;
  print_v2_tot << "# Cut off pt: " << lower_pt_cutoff_for_v2 << " GeV" << endl;
  print_v2_tot << "#4321\t#average integrated v2, for all rapidity bins" << endl; 
  print_v2_tot << 4321;
  print_v2_tot.width( 15 );
  for ( int i = 0;i < eta_bins;i++ )
  {
    if ( NmbInRange[i] > 0 )
    {
      print_v2_tot <<  v2sumPtCutOff[i] / NmbInRangePtCutOff[i];
    }
    else
    {
      print_v2_tot <<  0;
    }
    print_v2_tot << '\t' ;
  } 
  print_v2_tot << endl;
  print_v2_tot << "[#sum\t#number], for all rapidity bins, for lower pt cutoff" << endl; 
  for ( int i = 0;i < eta_bins;i++ )
  {
    print_v2_tot <<  v2sumPtCutOff[i] << '\t';
    print_v2_tot <<  NmbInRangePtCutOff[i] << '\t';
  }
  print_v2_tot << endl;
  print_v2_tot << "#5678\t#Total number of photons:" << endl;
  print_v2_tot << "5678" << '\t' <<  noninteractingParticles.size() << endl;

  

  //V2 over pt
  // print summed output, v2 is not computed, but summed v2 and the number in one bin
  print_v2_summed << "# summed v2 of " << type << endl;
  print_v2_summed << "# t = " << _outputTime << endl;
  print_v2_summed << "# Total summed photon v2" << "\t" << theConfig->v2average_debug << endl;
  print_v2_summed << "#";
  print_v2_summed.width( 14 );
  print_v2_summed << "#pt\t[#summed v_2\t#number in bin],for different rapidity bins" << endl;
  

  for ( int k = 0;k < v2pt_binnumber + 1;k++ )
  {
    print_v2_summed.width( 15 );
    pt_out = exp( double( k ) * d_ln_pt_v2 + log( _pt_min_v2 ) + d_ln_pt_v2 / 2.0 ); // +dpt/2.0 to shift value in the middle of the intervall, important for small number of bins
    print_v2_summed << pt_out;
    for ( int i = 0;i < eta_bins;i++ )
    {
      print_v2_summed.width( 15 );
      print_v2_summed << ptBinsV2[i][k];
      print_v2_summed.width( 10 );
      print_v2_summed << ptBinsNmbV2[i][k];
      if(studyInitialCutOffEffect)
      {
        print_v2_summed.width( 20 );
        print_v2_summed << ptBinsV2InitialCutOff[i][k];      
        print_v2_summed.width( 10 );
        print_v2_summed << ptBinsNmbInitialCutOffV2[i][k];       
      }
    }
    print_v2_summed << endl;
  }

  cout << "differential v2 printout done" << endl;
  
  //*****
  //***** Print the PT-Spectrum
  //*****
  // print yield for RAA
  print_yield << "# " << type << " yield distribution" << endl;
  print_yield << "# t = " << _outputTime << " fm" << endl;
  
  print_yield << "#pt     \t";
  for ( int i = 0;i < eta_bins;i++ )
  {
    const double delta_eta =  ( rapidityRanges[i].yright - rapidityRanges[i].yleft );
    print_yield << "#Tot_Rap_" << delta_eta << "\t#Compt_Rap_"<<delta_eta<<"\t#Annih_Rap_"<<delta_eta;   
  }
  print_yield << endl;
  for ( int k = 0;k < n_bins + 1;k++ )
  {
    //print_yield.width( 15 );
    pt_out = exp( double( k ) * d_ln_pt + log( _pt_min ) + d_ln_pt / 2.0 ); // +dpt/2.0 to shift value in the middle of the intervall, important for small number of bins
    
    print_yield << pt_out;
    
    const double dpt = exp( double( k ) * d_ln_pt + log( _pt_min ) + d_ln_pt ) - exp( double( k ) * d_ln_pt + log( _pt_min ) );
    
    for ( int i = 0;i < eta_bins;i++ )
    {
      const double delta_eta = 2.0 * ( rapidityRanges[i].yright - rapidityRanges[i].yleft );
      
      //WARNING: Moritz likes it to divide also by pt_out and 2pi
      //WARNING: DIVIDE MANUALLY BY K_FACTOR!!!!
      //double nInBin = double( ptBinsNmb[i][k] ) / (double(theConfig->getkFactorEMProcesses())*double(theConfig->getTestparticles())* double(dpt) * double(delta_eta) * double(pt_out) *(2.0*M_PI));
      double nInBin = double( ptBinsNmb[i][k] ) / (double(theConfig->getTestparticles())* double(dpt) * double(delta_eta) * double(pt_out) *(2.0*M_PI));
      
      if( _v2type == v2jets )
        nInBin = nInBin / theConfig->getNaddedEvents();
      
      if( Particle::mapToGenericFlavorType( _flavTypeToComputeFor ) == jpsi )
        nInBin = nInBin / theConfig->getJpsiTestparticles();
      
      if( Particle::mapToGenericFlavorType( _flavTypeToComputeFor ) == electron_gen )
        nInBin = nInBin / theConfig->getNumberElectronStat();
      
      print_yield.width( 15 );
      print_yield << nInBin;
      
      if( theConfig->doScattering_22_photons())
      {
        double nInBinCompton = double( ptBinsNmbCompton[i][k] ) / (double(theConfig->getTestparticles())* double(dpt) * double(delta_eta) * double(pt_out) *(2.0*M_PI));
        double nInBinAnnihilation = double( ptBinsNmbAnnihilation[i][k] ) / (double(theConfig->getTestparticles())* double(dpt) * double(delta_eta) * double(pt_out) *(2.0*M_PI));
        print_yield.width( 15 );
        print_yield << nInBinCompton;  
        print_yield.width( 15 );
        print_yield << nInBinAnnihilation;  
      }

      //DEBUG:
      /*if (pt_out >0.5 && pt_out < 0.7 && delta_eta==0.7)
      {
        //cout << "k-Faktor: " << theConfig->getkFactorEMProcesses()<<endl;
        cout << "pt=" << pt_out << "\t dN/2piptdptdy=" << nInBin << endl;
      }*/
      
      
    }
    print_yield << endl;
  }
  
  
  // print yield for RAA for different angles with respect to the reaction plane
  print_pt_angleDependence << "# " << type << " yield distribution for different angles (alpha) with respect to the reaction plane" << endl;
  print_pt_angleDependence << "# t = " << _outputTime << endl;
  print_pt_angleDependence << "#";
  print_pt_angleDependence.width( 14 );
  print_pt_angleDependence << "pt       yield for different rapidity bins" << endl;
  for ( int j = 0; j < nAlphaBins; j++ )
  {
    print_pt_angleDependence << "#alpha in [ " << j * deltaAlpha << "°, " << (j+1)*deltaAlpha << "° ] "<< endl;
    for ( int k = 0;k < n_bins + 1;k++ )
    {
      print_pt_angleDependence.width( 15 );
      pt_out = exp( double( k ) * d_ln_pt + log( _pt_min ) + d_ln_pt / 2.0 ); // +dpt/2.0 to shift value in the middle of the intervall, important for small number of bins
      print_pt_angleDependence << pt_out;
      const double dpt = exp( double( k ) * d_ln_pt + log( _pt_min ) + d_ln_pt ) - exp( double( k ) * d_ln_pt + log( _pt_min ) );
      for ( int i = 0;i < eta_bins;i++ )
      {
        const double delta_eta = 2.0 * ( rapidityRanges[i].yright - rapidityRanges[i].yleft );
        
        print_pt_angleDependence.width( 15 );
        print_pt_angleDependence << double( ptBinsAngleDep[i][j][k] ) / theConfig->getTestparticles() / dpt / delta_eta;
      }
      print_pt_angleDependence << endl;
    }    
    print_pt_angleDependence << endl;
    print_pt_angleDependence << endl;
  }
  
    //CMS pt cut for non prompt jpsi: 6.5GeV < pt < 30 GeV, y cut: |y|<2.4
  if( theConfig->isStudyNonPromptJpsiInsteadOfElectrons() && _flavTypeToComputeFor == electron_gen && ( theConfig->getOutputScheme() == cms_hq_nonPromptJpsi || theConfig->getOutputScheme() == alice_hq_nonPromptJpsi ) )
  {
    string filename_yield_nonPromptJpsi;
    string text;
    
    double eta_jpsi, pt_min_jpsi, pt_max_jpsi;
    
    if( theConfig->getOutputScheme() == cms_hq_nonPromptJpsi )
    {
      filename_yield_nonPromptJpsi = filename_prefix + "_nonPromptJpsiCMScuts_yield_" + name;
      text = "# with same acceptance cuts as for CMS:  6.5GeV < pt < 30 GeV,  |y|<2.4";
      eta_jpsi = 2.4;
      pt_min_jpsi = 6.5;
      pt_max_jpsi = 30.0;
    }
    else if( theConfig->getOutputScheme() == alice_hq_nonPromptJpsi )
    {
      filename_yield_nonPromptJpsi = filename_prefix + "_nonPromptJpsiALICEcuts_yield_" + name;
      text = "# with same acceptance cuts as for ALICE:  2GeV < pt < 30 GeV,  |y|<0.9";
      eta_jpsi = 0.9;
      pt_min_jpsi = 2.0;
      pt_max_jpsi = 30.0;
    }

    fstream print_nonPromptJpsi( filename_yield_nonPromptJpsi.c_str(), ios::out | ios::trunc );
    
    print_nonPromptJpsi << "# yield of non prompt jpsi from B meson decays" << endl;
    print_nonPromptJpsi << text << endl;
    print_nonPromptJpsi << "# <pt>     dN/dy" << endl;
    
    int count_jpsi = 0;
    double pt_sum = 0;
    double delta_eta_jpsi = 2.0 * eta_jpsi;
    double y, pt;
    for ( int i = 1; i <= n_particles; i++ )
    {
      pt = _particles[i].Mom.Pt();
      y  = _particles[i].Mom.Rapidity();
      
      if( fabs( y ) < eta_jpsi && pt > pt_min_jpsi && pt < pt_max_jpsi )
      {
        count_jpsi++;
        pt_sum += pt;
      }
    }
    
    print_nonPromptJpsi << pt_sum / double( count_jpsi ) << "\t";
    
    print_nonPromptJpsi << double( count_jpsi ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getNumberElectronStat() / delta_eta_jpsi << endl;
  }
  
}




void analysis::volumeMidrap( const int step ) const
{
  string name;
  stringstream ss;

  ss << step;
  name = "step" + ss.str();


  string filename_z = filename_prefix + "_zMidrap_" + name;
  string filename_t = filename_prefix + "_tMidrap_" + name;

  double y, xt;
  const double y_max = 0.5;

  int numberOfBins = 80;

  binning zBins( filename_z, numberOfBins );
  binning tBins( filename_t, numberOfBins );

  for ( int i = 0;i < particles_atTimeNow.size();i++ )
  {
//     y = 0.5*log( (particles_atTimeNow[i].E + particles_atTimeNow[i].PZ) / (particles_atTimeNow[i].E - particles_atTimeNow[i].PZ) );
//     if(fabs(y) <= y_max)
//     {
//       zBins.add( fabs(particles_atTimeNow[i].Z) );
//       xt = sqrt(pow(particles_atTimeNow[i].X,2.0) + pow(particles_atTimeNow[i].Y,2.0));
//       tBins.add(xt);
//     }

    if ( particles_atTimeNow[i].Pos.T() <= tstep[step] )
    {
      y = particles_atTimeNow[i].Mom.Rapidity();
      if ( fabs( y ) <= y_max )
      {
        zBins.add( fabs( particles_atTimeNow[i].Pos.Z() ) );
        tBins.add( particles_atTimeNow[i].Pos.Perp() );
      }



//       zBins.add( fabs(particles_atTimeNow[i].Pos.Z()) );
//       tBins.add( particles_atTimeNow[i].Pos.Perp() );
    }

  }

  zBins.print();
  tBins.print();
}


void analysis::ptDistribution( const FLAVOR_TYPE _flavTypeToComputeFor, vector<ParticleOffline>& _particles, const int n_particles, const int step )
{
  FLAVOR_TYPE genFlavor;
  double pt, y;
  tVecOfArrayOfDoubleVec _ptBins;
  
  // populate temporary _ptBins with temporary "references" to the actual data structures according to the requested flavor
  for ( int yRange = 0; yRange < rapidityRanges.size(); yRange++ )
  {
    switch ( _flavTypeToComputeFor )
    {
      case gluon:
        _ptBins.push_back( ptBins_gluons[yRange] );
        break;
      case light_quark:
        _ptBins.push_back( ptBins_quarks[yRange] );
        break;
      case up:
        _ptBins.push_back( ptBins_ups[yRange] );
        break;
      case down:
        _ptBins.push_back( ptBins_downs[yRange] );
        break;
      case strange:
        _ptBins.push_back( ptBins_stranges[yRange] );
        break;
      case anti_up:
        _ptBins.push_back( ptBins_anti_ups[yRange] );
        break;
      case anti_down:
        _ptBins.push_back( ptBins_anti_downs[yRange] );
        break;  
      case anti_strange:
        _ptBins.push_back( ptBins_anti_stranges[yRange] );
        break;
      case allFlavors:
        _ptBins.push_back(ptBins_all[yRange]);
        break;
      default:
        string errMsg = "error in ptDistribution, flavor not specified";
        throw eAnalysis_error( errMsg );
        break;
    }
  }
  
  // loop over all particles and bin them according to their pt
  for ( int j = 0; j < n_particles; j++ )
  {
    pt = _particles[j].Mom.Pt();
    y = _particles[j].Mom.Rapidity();
    
    // check whether particle has the correct flavor
    genFlavor = ParticleOffline::mapToGenericFlavorType( _particles[j].FLAVOR );
    if ( _flavTypeToComputeFor == allFlavors || _particles[j].FLAVOR == _flavTypeToComputeFor || 
      ( _flavTypeToComputeFor == light_quark && ( genFlavor == light_quark || genFlavor == anti_light_quark ) ) )
    {
      // is it in the possible pt range at all?
      if ( pt < maxPT && pt >= minPT )
      {
        // individually check for each rapidity range whether this particle needs to be binned
        for ( int yRangeIndex = 0; yRangeIndex < rapidityRanges.size(); yRangeIndex++ )
        {
          if ( fabs( y ) >= rapidityRanges[yRangeIndex].yleft && fabs( y ) <= rapidityRanges[yRangeIndex].yright )
          {
            if ( pt == minPT )  // a special case
            {
              ++_ptBins[yRangeIndex][step][0];  // actually bin the pt of the particle
            }
            else
            {
              ++_ptBins[yRangeIndex][step][int(( pt - minPT )/binWidthPT )];  // actually bin the pt of the particle
            }
          }
        }
      }
    }
  }
  
}



void analysis::ptSoftDistribution( const FLAVOR_TYPE _flavTypeToComputeFor, vector<ParticleOffline>& _particles, const int n_particles, const int step )
{
  double pt;

  vector<double>* _ptBinsSoftAll;

  if ( _flavTypeToComputeFor == gluon )
  {
    _ptBinsSoftAll = ptBinsSoftAll_gluons;
  }
  else if ( _flavTypeToComputeFor == light_quark )
  {
    _ptBinsSoftAll = ptBinsSoftAll_quarks;
  }
  else if ( _flavTypeToComputeFor == up )
  {
    _ptBinsSoftAll = ptBinsSoftAll_ups;
  }
  else if ( _flavTypeToComputeFor == down )
  {
    _ptBinsSoftAll = ptBinsSoftAll_downs;
  }
  else if ( _flavTypeToComputeFor == strange )
  {
    _ptBinsSoftAll = ptBinsSoftAll_stranges;
  }
  else if ( _flavTypeToComputeFor == anti_up )
  {
    _ptBinsSoftAll = ptBinsSoftAll_anti_ups;
  }
  else if ( _flavTypeToComputeFor == anti_down )
  {
    _ptBinsSoftAll = ptBinsSoftAll_anti_downs;
  }
  else if ( _flavTypeToComputeFor == anti_strange )
  {
    _ptBinsSoftAll = ptBinsSoftAll_anti_stranges;
  }
  else if ( _flavTypeToComputeFor == allFlavors )
  {
    _ptBinsSoftAll = ptBinsSoftAll_all;
  }
  else
  {
    string errMsg = "error in ptSoftDistribution, flavor not specified";
    throw eAnalysis_error( errMsg );
  }

  FLAVOR_TYPE genFlavor;
  for ( int j = 0; j < n_particles; j++ )
  {
    pt = _particles[j].Mom.Pt();

    genFlavor = ParticleOffline::mapToGenericFlavorType( _particles[j].FLAVOR );
    if ( _flavTypeToComputeFor == allFlavors || _particles[j].FLAVOR == _flavTypeToComputeFor || 
      ( _flavTypeToComputeFor == light_quark && ( genFlavor == light_quark || genFlavor == anti_light_quark ) ) )
    {
      //------------------------ y in [-inf,inf] -----------------------
      if ( pt <= maxPTSoft && pt >= 0 )
      {
        if ( pt == maxPTSoft )
        {
          ++_ptBinsSoftAll[step][numberBinsSoftPT - 1];
        }
        else
        {
          if ( pt == 0 )
            ++_ptBinsSoftAll[step][0];
          else
            ++_ptBinsSoftAll[step][int( pt/binWidthSoftPT )];
        }
      }
      //----------------------------------------------------------------
    }
  }
}



void analysis::yDistribution( const FLAVOR_TYPE _flavTypeToComputeFor, vector<ParticleOffline>& _particles, const int n_particles, const int step )
{
  double y;
  vector<double> * _yBins;


  if ( _flavTypeToComputeFor == gluon )
  {
    _yBins = yBins_gluon;
  }
  else if ( _flavTypeToComputeFor == up )
  {
    _yBins = yBins_up;
  }
  else if ( _flavTypeToComputeFor == down )
  {
    _yBins = yBins_down;
  }
  else if ( _flavTypeToComputeFor == strange )
  {
    _yBins = yBins_strange;
  }
  else if ( _flavTypeToComputeFor == anti_up )
  {
    _yBins = yBins_anti_up;
  }
  else if ( _flavTypeToComputeFor == anti_down )
  {
    _yBins = yBins_anti_down;
  }
  else if ( _flavTypeToComputeFor == anti_strange )
  {
    _yBins = yBins_anti_strange;
  }
  else
  {
    string errMsg = "error in yDistribution, flavor not specified";
    throw eAnalysis_error( errMsg );
  }

  for ( int j = 0; j < n_particles; j++ )
  {
    y = _particles[j].Mom.Rapidity();

    if ( _particles[j].FLAVOR == _flavTypeToComputeFor )
    {
      if ( y <= maxY && y >= minY)
      {
        if ( y == maxY )
        {
          ++_yBins[step][numberBinsPT - 1];
        }
        else
        {
          if ( y == minY )
            ++_yBins[step][0];
          else
            ++_yBins[step][int(( y - minY )/binWidthY )];
        }
      }
      //----------------------------------------------------------------
    }
  }
}




void analysis::transverseEnergyDistribution( const FLAVOR_TYPE _flavTypeToComputeFor, vector<ParticleOffline>& _particles, const int n_particles, const int step )
{
  double Et, y;
  vector<double> * _EtBins;
  
  if ( _flavTypeToComputeFor == gluon )
  {
    _EtBins = transverseEnergyGluons;
  }
  else if ( _flavTypeToComputeFor == light_quark )
  {
    _EtBins = transverseEnergyQuarks;
  }
  else if ( _flavTypeToComputeFor == anti_light_quark )
  {
    _EtBins = transverseEnergyAntiQuarks;
  }
  else
  {
    string errMsg = "error in yDistribution, flavor not specified";
    throw eAnalysis_error( errMsg );
  }
  
  for ( int j = 0; j < n_particles; j++ )
  {
    Et = _particles[j].Mom.Pt();
    y = _particles[j].Mom.Rapidity();
    
    FLAVOR_TYPE genFlavor = ParticleOffline::mapToGenericFlavorType( _particles[j].FLAVOR );
    if ( (  genFlavor == _flavTypeToComputeFor ) ||
      ( _flavTypeToComputeFor == light_quark && ( genFlavor == light_quark || genFlavor == anti_light_quark ) ) )
    {
      if ( y <= maxY && y >= minY)
      {
        if ( y == maxY )
        {
          _EtBins[step][numberBinsPT - 1] += Et;
        }
        else
        {
          if ( y == minY )
            _EtBins[step][0] += Et;
          else
            _EtBins[step][int(( y - minY )/binWidthY )] += Et;
        }
      }
      //----------------------------------------------------------------
    }
  }
}




void analysis::removeJetEvent_in( const int entity_ID )
{
  jetTracker[entity_ID].pop_back();
}


void analysis::makeJetTrackerCopy()
{
  jetTracker_copy = jetTracker;
}


void analysis::restoreJetTracker()
{
  jetTracker = jetTracker_copy;
}


void analysis::exchangeJetID( const int oldID, const int newID )
{
  int entity_index = 0;
  while ( entity_index < jetTracker.size() && jetTracker[entity_index].back().jet_ID_out != -addedParticles[oldID].unique_id )
  {
    entity_index++;
  }

  if ( entity_index < jetTracker.size() )
  {
    jetTracker[entity_index].back().jet_ID_out = -addedParticles[newID].unique_id;
  }
  else
  {
    cout << "error in exchangeJetID()" << endl;
  }

}



void analysis::addJetEvent_initial( const int jetID )
{
  jetTrackerSingleEvent tempEvent;
  tempEvent.jet_ID_in = -1;
  tempEvent.jet_ID_out = -addedParticles[jetID].unique_id;
  tempEvent.coll_type = initial_jet;
  
  tempEvent.flavor_in = -1;
  tempEvent.flavor_out = static_cast<int>( addedParticles[jetID].FLAVOR );

  tempEvent.R_proj = addedParticles[jetID].Pos;
  tempEvent.P_proj_out = addedParticles[jetID].Mom;

  vector<jetTrackerSingleEvent> tempVec;
  tempVec.push_back( tempEvent );
  jetTracker.push_back( tempVec );
}


void analysis::addJetEvents_final()
{
  double pt;
  for ( int i = 0; i < addedParticles.size(); i++ )
  {
    pt = addedParticles[i].Mom.Pt();
    if ( pt > jetTracking_PT )
    {
      jetTrackerSingleEvent tempEvent;
      tempEvent.jet_ID_in = -addedParticles[i].unique_id;
      tempEvent.jet_ID_out = -1;
      tempEvent.R_proj = addedParticles[i].Pos;
      tempEvent.P_proj_in = addedParticles[i].Mom;
      
      tempEvent.flavor_out = -1;
      tempEvent.flavor_in = static_cast<int>( addedParticles[i].FLAVOR );
      
      tempEvent.coll_type = final_jet;

      int entity_index = 0;
      while ( entity_index < jetTracker.size() && jetTracker[entity_index].back().jet_ID_out != -addedParticles[i].unique_id )
      {
        entity_index++;
      }

      if ( entity_index < jetTracker.size() )
      {
        jetTracker[entity_index].push_back( tempEvent );
      }
    }
  }

}



int analysis::addJetEvent_in( const int ID_1, const int ID_2, const int added_ID, const jetTrackerCollType coll_type,
                              const double cross_section, const int cell_ID, const double lambda )
{
//   double jet_pt = sqrt( pow( addedParticles[added_ID].PX, 2.0 ) + pow( addedParticles[added_ID].PY, 2.0 ) );
//   double pt1 = sqrt( pow( particles_atTimeNow[ID_1].PX, 2.0 ) + pow( particles_atTimeNow[ID_1].PY, 2.0 ) );
//   double pt2 = -1;
//   if ( ID_2 > 0 )
//   {
//     pt2 = sqrt( pow( particles_atTimeNow[ID_2].PX, 2.0 ) + pow( particles_atTimeNow[ID_2].PY, 2.0 ) );
//   }

  int jetID = added_ID;
  int partner1 = ID_1;
  int partner2 = ID_2;


  jetTrackerSingleEvent tempEvent;
  tempEvent.jet_ID_in = -addedParticles[jetID].unique_id;
  tempEvent.R_proj = addedParticles[jetID].Pos;
  tempEvent.flavor_in = static_cast<int>( addedParticles[jetID].FLAVOR );

  tempEvent.P_proj_in = addedParticles[jetID].Mom;
  tempEvent.P1_in = particles_atTimeNow[partner1].Mom;

  if ( coll_type == c3to2 )
  {
    tempEvent.P2_in = particles_atTimeNow[partner2].Mom;
  }

  tempEvent.lambda = lambda;
  tempEvent.xSection = cross_section;
  tempEvent.cell_ID = cell_ID;
  tempEvent.coll_type = coll_type;

  int entity_index = 0;
  while ( entity_index < jetTracker.size() && jetTracker[entity_index].back().jet_ID_out != -addedParticles[jetID].unique_id )
  {
    entity_index++;
  }

  if ( entity_index < jetTracker.size() )
  {
    jetTracker[entity_index].push_back( tempEvent );
  }
  else
  {
    vector< jetTrackerSingleEvent > tempEntity;
    tempEntity.push_back( tempEvent );
    jetTracker.push_back( tempEntity );
  }

  return entity_index;
}


void analysis::addJetEvent_out( const int entity_ID, const int added_ID, const int ID_2, const int ID_3, const jetTrackerCollType coll_type )
{
  double pt1 = addedParticles[added_ID].Mom.Pt();
  double pt2 = particles_atTimeNow[ID_2].Mom.Pt();

  double pt3 = -1;
  if ( coll_type == c2to3 )
  {
    pt3 = addedParticles[ID_3].Mom.Pt();
  }

  int jetID;
  int partner1 = -1, partner2 = -1;

  jetID = added_ID;
  partner1 = ID_2;
  partner2 = ID_3;
  
  if ( pt3 > pt1 )
  {
    jetID = ID_3;
    partner2 = added_ID;
  }

  if ( entity_ID != -1 )
  {
    jetTracker[entity_ID].back().jet_ID_out = -addedParticles[jetID].unique_id;
    jetTracker[entity_ID].back().P_proj_out = addedParticles[jetID].Mom;
    
    jetTracker[entity_ID].back().flavor_out = static_cast<int>( addedParticles[jetID].FLAVOR );
    jetTracker[entity_ID].back().P1_out = particles_atTimeNow[partner1].Mom;

    if ( coll_type == c2to3 )
    {
      jetTracker[entity_ID].back().P2_out = addedParticles[partner2].Mom;
    }
  }
  else if ( addedParticles[jetID].Mom.Pt2() > pow( jetTracking_PT, 2 ) )
  {
    jetTrackerSingleEvent tempEvent;
    tempEvent.jet_ID_in = -1;
    tempEvent.jet_ID_out = -addedParticles[jetID].unique_id;
    tempEvent.coll_type = production;
    tempEvent.R_proj = addedParticles[jetID].Pos;
    tempEvent.flavor_in = -1;
    tempEvent.flavor_out = static_cast<int>( addedParticles[jetID].FLAVOR );

    tempEvent.P_proj_out = addedParticles[jetID].Mom;
    tempEvent.P1_out = particles_atTimeNow[partner1].Mom;

    if ( coll_type == c2to3 )
    {
      tempEvent.P2_out = addedParticles[partner2].Mom;
    }

    vector< jetTrackerSingleEvent > tempEntity;
    tempEntity.push_back( tempEvent );
    jetTracker.push_back( tempEntity );
  }
}



void analysis::particleOutput( const int step )
{
  string name;
  stringstream ss;

  if ( step == 0 )
    name = "initial";
  else if ( step == nTimeSteps )
    name = "final";
  else
  {
    ss << step;
    name = "step" + ss.str();
  }

  //creates filename, for example: "./output/run32_step1.f1", "./output/run32_initial.f1" or the likes
  string filename = filename_prefix + "_" + name + ".f1";
  fstream file( filename.c_str(), ios::out | ios::trunc );

  //---- print header if file is empty ----
  time_t end;
  time( &end );

  file.seekp( 0, ios::end );
  long size = file.tellp();
  file.seekp( 0, ios::beg );
  if ( size == 0 )
    printHeader( file, all, end );
  //---------------------------------------

  for ( int i = 0; i < addedParticles.size(); i++ )
  {
    file << i << sep << addedParticles[i].unique_id << sep << addedParticles[i].cell_id << sep << addedParticles[i].FLAVOR << sep 
         << addedParticles[i].Pos.T() << sep << addedParticles[i].Pos.X() << sep
         << addedParticles[i].Pos.Y() << sep << addedParticles[i].Pos.Z() << sep 
         << addedParticles[i].Mom.E() << sep << addedParticles[i].Mom.Px() << sep 
         << addedParticles[i].Mom.Py() << sep << addedParticles[i].Mom.Pz() << sep 
      //         << addedParticles[i].md2g << sep << addedParticles[i].md2q << sep 
      //         << addedParticles[i].N_EVENT_pp << endl;
         << addedParticles[i].md2g << sep 
         << addedParticles[i].X_traveled << sep 
         << addedParticles[i].N_EVENT_pp << endl;
  }
  file.close();
  
  // if heavy quarks fragment to D mesons write also the D mesons to file
  if( theConfig->isHadronizationHQ() )
  {
    //creates filename, for example: "./output/run32_step1.f1", "./output/run32_initial.f1" or the likes
    string filename_meson = filename_prefix + "_" + name + "_meson.f1";
    fstream file_meson( filename_meson.c_str(), ios::out | ios::trunc );

    //---- print header if file is empty ----
    time_t end;
    time( &end );

    file_meson.seekp( 0, ios::end );
    long size = file_meson.tellp();
    file_meson.seekp( 0, ios::beg );
    if ( size == 0 )
      printHeader( file_meson, all, end );
    //---------------------------------------

    for ( int i = 0; i < addedParticlesCopy.size(); i++ )
    {
      file_meson << i << sep << addedParticlesCopy[i].unique_id << sep << addedParticlesCopy[i].cell_id << sep << addedParticlesCopy[i].FLAVOR << sep 
                 << addedParticlesCopy[i].Pos.T() << sep << addedParticlesCopy[i].Pos.X() << sep
                 << addedParticlesCopy[i].Pos.Y() << sep << addedParticlesCopy[i].Pos.Z() << sep 
                 << addedParticlesCopy[i].Mom.E()  << sep << addedParticlesCopy[i].Mom.Px() << sep 
                 << addedParticlesCopy[i].Mom.Py() << sep << addedParticlesCopy[i].Mom.Pz() << sep 
        //                 << addedParticlesCopy[i].md2g << sep << addedParticlesCopy[i].md2q << sep << addedParticlesCopy[i].N_EVENT_pp << endl;
                 << addedParticlesCopy[i].md2g << sep << addedParticlesCopy[i].X_traveled << sep << addedParticlesCopy[i].N_EVENT_pp << endl;
    }
    file_meson.close();
  }
}



void analysis::writePartclMovie( vector< ParticleOffline >& _particles, const int n_particles, fstream& _oscar, const int step, const int jumpSteps )
{
  const string sep = "  ";
  const int width = 14;
  double cc, dt;
  VectorTXYZ Pos;
  const double zero = 0.0;

  double time;

  const int nOutput = 1; // every nOutput-th particle is writen to particle file, eg nOutput=2 every second particle, nOutput=1 for every particle
  int nCount_selected = nOutput; // dummy variable for counting

  if ( step == 0 )
  {
    time = 0.0;
  }
  else if ( step == nTimeSteps_movie - 1 )
  {
    return;
  }
  else
  {
    time = tstep_movie[step];
  }

  // determine number of timesteps, nTimeSteps is too large if runtime is shorter than last timestep
  int numberOfTimeSteps = 0;
  while ( theConfig->getRuntime() >= tstep_movie[numberOfTimeSteps] && numberOfTimeSteps < nTimeSteps_movie - 1 )
  {
    numberOfTimeSteps++;
  }
  numberOfTimeSteps -= jumpSteps;

  int numberActiveParticles = 0;
  for ( int i = 0; i < n_particles; i++ )
  {
    if ( _particles[i].T_creation <= time )
    {
      if ( nCount_selected != nOutput ) // do not write particle to file, increase count by one
      {
        nCount_selected++;
      }
      else // write particle to file
      {
        numberActiveParticles++;
        nCount_selected = 1;
      }
    }
  }
  nCount_selected = nOutput; // reset dummy variable for counting



  // <<---------------------------------------------------
  // oscar output

  if ( step == 0 ) // header, no particles yet
  {
    // file header
    _oscar << "OSC1997A    " << endl;
    _oscar << "final_id_p_x" << endl;

    // code name etc.
    // code_name version (aproj, zproj)+(atarg, ztarg) refframe, ebeam, ntestpart
    _oscar << ::std::scientific << ::std::uppercase
    << "BAMPS        0.2    (" << ::std::setw( 3 ) << int( theConfig->getA() ) << "," << ::std::setw( 6 )
    <<  int( theConfig->getAatomic() ) << ")+(" << ::std::setw( 3 ) << int( theConfig->getB() ) << ","
    << ::std::setw( 6 ) << int( theConfig->getBatomic() ) << ")  CMS " << ::std::setprecision( 4 )
    << theConfig->getSqrtS() << "  " << ::std::setw( 8 ) << theConfig->getTestparticles() << endl;
  }
  else // all other timesteps
  {
    // event header for this timestep
    // event npart bimp phi #timesteps timestep
    _oscar << "          1  " << ::std::setw( 10 ) << numberActiveParticles << "  " << ::std::resetiosflags( ::std::ios::scientific )
    << ::std::setw( 8 ) << ::std::fixed << ::std::setprecision( 3 ) << theConfig->getImpactParameter() << "  "
    << ::std::setw( 8 ) << ::std::fixed << ::std::setprecision( 3 ) << 0.000 << "  "  << ::std::setw( 4 )
    << numberOfTimeSteps << "  " << ::std::setw( 4 ) << step << endl;


    // write particle data
    // ipart, id, px, py, pz, p0, mass, x, y, z, t, flag#1, x_f, y_f, z_f, t_f, iflag#1
    for ( int i = 0; i < n_particles; i++ )
    {
      if ( _particles[i].T_creation <= time )
      {
        // find out if particle should be written to file (in case nOutput is not 1)
        if ( nCount_selected != nOutput ) // do not write particle to file, increase count by one
        {
          nCount_selected++;
        }
        else // write particle to file
        {
          _oscar << ::std::setw( 11 ) << i+1;
          _oscar << ::std::setw( 12 ) << ParticleOffline::mapToPDGCodes( _particles[i].FLAVOR ); //PDG group codes
          nCount_selected = 1;

          _oscar << ::std::scientific << ::std::uppercase << ::std::setprecision( 6 )
                 << setw( width ) << _particles[i].Mom.Px() 
                 << setw( width ) << _particles[i].Mom.Py() 
                 << setw( width ) << _particles[i].Mom.Pz() 
                 << setw( width ) << _particles[i].Mom.E() 
                 << setw( width ) << _particles[i].m;

          if ( _particles[i].Pos.T() < time )
          {
            cout << "error in write movie particle data, particles from the past" << endl;
            cout << "timestep=" << time << "   time _particles=" << _particles[i].Pos.T() << endl;
          }
          
          Pos = _particles[i].Pos + _particles[i].Mom * (( time - _particles[i].Pos.T() )/_particles[i].Mom.E() );


          _oscar << setw( width ) << Pos.X() << setw( width ) << Pos.Y() << setw( width ) << Pos.Z() << setw( width ) << Pos.T();

          //       _oscar << "0.0" << sep << "0.0" << sep << "0.0" << sep << "0.0" << sep << "0.0" << sep << "0" << endl;
          //         _oscar << zero << sep << zero << sep << zero << sep << zero << sep << zero << ::std::setw( 10 ) << int(zero) << endl;
          //         _oscar << zero << sep << _particles[i].X_lastInt << sep << _particles[i].Y_lastInt << sep << _particles[i].Z_lastInt << sep << _particles[i].T_lastInt << ::std::setw( 10 ) << int(zero) << endl;
          _oscar << setw( width ) << zero 
                 << setw( width ) << _particles[i].lastInt.X()
                 << setw( width ) << _particles[i].lastInt.Y() 
                 << setw( width ) << _particles[i].lastInt.Z() 
                 << setw( width ) << _particles[i].lastInt.T() 
                 << ::std::setw( 10 ) << -(_particles[i].unique_id) << endl;
          // --------------------------------------------------->>
        }
      }
    }
  }
}






void analysis::printHeader( fstream & f, const anaType mode, const time_t end )
{
  switch ( mode )
  {
  case cellV2DistributionHeader:
    f << "#V2-distribution of cells of background particles "<<endl;
    break;
  case numbOfPartcles:
    f << "#Number of Quarks and Gluons" << endl;
    break;
  case dEtdy:
  f << "#Et-spectrum" << endl;
  break;
  case ptSpectrum:
    f << "#pt-spectra " << endl;
    break;
  case ptSpectrumSoft:
    f << "#soft pt-spectra " << endl;
    break;
  case jets:
    f << "#jet tracking informatino " << endl;
    break;
  case all:
    f << "#information on particle locations, momenta etc. at step indicated by filename" << endl;
    f << "#initial = 0 fm/c" << endl;
    for ( int i = 0;i < nTimeSteps - 1;i++ )
      f << "#step " << ( i + 1 ) << " = " << tstep[i] << " fm/c" << endl;
    f << "#final = " << theConfig->getRuntime() << " fm/c" << endl;
    break;
  case rapidityDistribution:
    f << "#rapidity distribution " << endl;
    f << "#one block per timestep " << endl;
    break;
  case quarkNumbers:  
    f << "#quark numbers summed over all rapidities" << endl;
    break;
  case photonPtDist:
    f << "#spectrum of produced photons - ";
    break;
  case dileptondRdMDist:
    f << "#spectrum of produced dileptons dR/dM - ";
    break;    
  case dileptondRdPDist:
    f << "#spectrum of produced dileptons dR/dP - ";
    break;    
  case dileptondRdMTimeDist:
    f << "#spectrum of produced dileptons dR/dM produced in different time slots ";
    break;  
  case dileptondNdMCuts:
    f << "#Invariant Mass Yield with experimental cuts applied.";
  default:
    f << "#undefined ";
  }

  f << "#start: " << ctime( &start );
  f << "#end: " << ctime( &end );
  f << "#" << endl;
  f << "#simulation parameter:" << endl;
  f << "#testparticles= " << theConfig->getTestparticles() << endl;
  f << "#runtime= " << theConfig->getRuntime() << endl;
  f << "#sqrtS= " << theConfig->getSqrtS() << " GeV" << endl;
  f << "#P0= " << theConfig->getPtCutoff() << " GeV" << endl;
  f << "#b= " << theConfig->getImpactParameter() << " fm" << endl;
  f << "#(" << theConfig->getA() << "," << theConfig->getAatomic() << ") on ("
  << theConfig->getB() << "," << theConfig->getBatomic() << ")" << endl;
  f << "#seed for random generator ran2(): " << seed << endl;
  f << "#" << endl;

  stringstream ss;

  switch ( mode )
  {
  case ptSpectrum:
    f << "#numbers NOT yet corrected for width of bins and number of testparticles!" << endl;
    f << "#binWidth= " << binWidthPT << endl;
    f << "#testparticles= " << theConfig->getTestparticles() << endl;
    f << "#" << endl;
    f << "#pt [GeV]" << sep;
    for ( int j = 0; j <= nTimeSteps; j++ )
    {
      if ( j == 0 )
        f << "initial" << sep;
      else if ( j == nTimeSteps )
        f << "final = " << theConfig->getRuntime() << " fm/c";
      else if ( tstep[j-1] <= theConfig->getRuntime() )
      {
        ss << tstep[j-1];
        f << ss.str() << " fm/c" << sep;
        ss.str( "" );
      }
    }
    f << endl;
    break;
  case ptSpectrumSoft:
    f << "#numbers NOT yet corrected for width of bins and number of testparticles!" << endl;
    f << "#binWidth= " << binWidthSoftPT << endl;
    f << "#testparticles= " << theConfig->getTestparticles() << endl;
    f << "#" << endl;
    f << "#pt [GeV]" << sep;
    for ( int j = 0; j <= nTimeSteps; j++ )
    {
      if ( j == 0 )
        f << "initial" << sep;
      else if ( j == nTimeSteps )
        f << "final = " << theConfig->getRuntime() << " fm/c";
      else if ( tstep[j-1] <= theConfig->getRuntime() )
      {
        ss << tstep[j-1];
        f << ss.str() << " fm/c" << sep;
        ss.str( "" );
      }
    }
    f << endl;
    break;
  case jets:
    f << "#jetID_in" << sep << "jetID_out" << sep << "coll_type" << sep
    << "P_in[0]" << sep << "P_in[1]" << sep << "P_in[2]" << sep  << "P_in[3]" << sep
    << "P_out[0]" << sep << "P_out[1]" << sep << "P_out[2]" << sep  << "P_out[3]" << sep
    << "R[0]" << sep << "R[1]" << sep << "R[2]" << sep << "R[3]" << sep
    << "xSection" << sep << "lambda" << endl;
    break;
  case rapidityDistribution:
    f << "# y   gluon   up   down   strange   anti-up   anti-down   anti-strange" << endl;
    break;
  case quarkNumbers:
    f << "# t   gluon   up   down   strange   anti-up   anti-down   anti-strange" << endl;
    break;
  case all:
    f << "#ID" << sep << "uniqueID" << sep << "Cell ID" << sep << "Flavor" << sep << "t" << sep << "x" << sep << "y" << sep << "z" << sep
    << "E"  << sep << "Px"  << sep << "Py"  << sep << "Pz" << sep << "md2g/alpha_s [GeV^2]" << sep
    << "md2q/alpha_s [GeV^2]" << sep << "Init. eventID" << endl;
    break;
  case photonPtDist:
    f << "#photons per energy bin (divided by bin width)" << endl;
    f << "# bin width = " << PhotondNOverTwoPiptdydptBin.getWidth() << endl;    
    break;
  default:
    f << endl;
  }
  f << "#" << endl;
}



void analysis::mfpJetsOutput( const int step, const int jumpSteps )
{
  const string sep = "\t ";
  double time;

  if ( step == 0 )
  {
    time = 0.0;
  }
  else if ( step == nTimeSteps_movie - 1 )
  {
    return;
  }
  else
  {
    time = tstep_movie[step];
  }
  
  mfpJetsOutputFile << time;
  for ( int i = 0; i < rings.size(); i++ )
  {
    if ( rings[i].collectedGluon != 0 )
    {
      mfpJetsOutputFile << sep << rings[i].lambdaGluon / rings[i].collectedGluon;
    }
    else
    {
      mfpJetsOutputFile << sep << 0; 
    }
    if ( rings[i].collectedQuark != 0 )
    {
      mfpJetsOutputFile << sep << rings[i].lambdaQuark / rings[i].collectedQuark;
    }
    else
    {
      mfpJetsOutputFile << sep << 0; 
    }
    mfpJetsOutputFile << sep << rings[i].collectedGluon << sep << rings[i].collectedQuark;  
  }
  mfpJetsOutputFile << endl;
  
  rings.clear();

}



void analysis::printCentralDensities(const double _time)
{
  const string sep = "\t ";
  
  centralDensitiesOutputFile << _time << sep << centralRingsCopyFromCascade.y_left << sep << centralRingsCopyFromCascade.y_right; 
  for ( int i = 0; i < centralRingsCopyFromCascade.size(); i++ )
  {
    centralDensitiesOutputFile << sep << centralRingsCopyFromCascade[i].getEnergyDensity() << sep 
    << centralRingsCopyFromCascade[i].getGluonDensity() << sep << centralRingsCopyFromCascade[i].getQuarkDensity() << sep
    << centralRingsCopyFromCascade[i].getEffectiveTemperature() << sep << centralRingsCopyFromCascade[i].getAveraged_md2g() << sep 
    << centralRingsCopyFromCascade[i].getAveraged_md2q(); 
  }
  centralDensitiesOutputFile << endl;

  centralRingsCopyFromCascade.clear();
}





// heavy quark stuff

void analysis::writeJpsiFugacityOutput( const int step )
{
  const string sep = "  ";

  double fugacity, temp, deltaTemp, enDen;
  int n_jpsi;
  double dr, dz, deta;


  double time;
  if ( step == 0 )
    time = 0.0;
  else if ( step == nTimeSteps )
    return;
  else
    time = tstep[step-1];

  if ( step == 0 )
  {
    // file header
    printJpsiFugacity << "#Jpsi fugacities" << endl;
    printJpsiFugacity << "# time  fugacities, N_ccbar, temp, energy dens. for different boxes" << endl;
    return;
  }

  printJpsiFugacity.width( 10 );
  printJpsiFugacity << time ;
  printJpsiFugacity<< "\t";


  dr = 2.0; //fm
//   dz=1.0; //fm
  deta = 0.5; // spacetime rapidty interval
  dz = time * ( exp( 2.0 * deta ) - 1.0 ) / ( exp( 2.0 * deta ) + 1.0 ); //translated to spatial coordinate z

//   cout << "deta=" << deta << "   dz=" << dz << endl;

  getJpsiFugacity( time, dr, dz, fugacity, n_jpsi, temp, enDen );
  printJpsiFugacity << fugacity;
  printJpsiFugacity << "\t";
  printJpsiFugacity << double(n_jpsi) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles();

  printJpsiFugacity << "\t";
  printJpsiFugacity << temp;
  printJpsiFugacity << "\t";
  printJpsiFugacity << enDen / theConfig->getTestparticles();


  printJpsiFugacity.width( 20 );

  dr = 5.0; //fm
//   dz=1.0; //fm
  deta = 0.5; // spacetime rapidty interval
  dz = time * ( exp( 2.0 * deta ) - 1.0 ) / ( exp( 2.0 * deta ) + 1.0 ); //translated to spatial coordinate z



  getJpsiFugacity( time, dr, dz, fugacity, n_jpsi, temp, enDen );
  printJpsiFugacity << fugacity;
  printJpsiFugacity << "\t";
  printJpsiFugacity << double(n_jpsi) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles();
  printJpsiFugacity << "\t";
  printJpsiFugacity << temp;
  printJpsiFugacity << "\t";
  printJpsiFugacity << enDen / theConfig->getTestparticles();



  printJpsiFugacity.width( 20 );

  dr = 2.0; //fm
//   dz=1.0; //fm
  deta = 1.0; // spacetime rapidty interval
  dz = time * ( exp( 2.0 * deta ) - 1.0 ) / ( exp( 2.0 * deta ) + 1.0 ); //translated to spatial coordinate z

//   cout << "deta=" << deta << "   dz=" << dz << endl;

  getJpsiFugacity( time, dr, dz, fugacity, n_jpsi, temp, enDen );
  printJpsiFugacity << fugacity;
  printJpsiFugacity << "\t";
  printJpsiFugacity << double(n_jpsi) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles();
  printJpsiFugacity << "\t";
  printJpsiFugacity << temp;
  printJpsiFugacity << "\t";
  printJpsiFugacity << enDen / theConfig->getTestparticles();



  printJpsiFugacity.width( 20 );

  dr = 5.0; //fm
//   dz=1.0; //fm
  deta = 1.0; // spacetime rapidty interval
  dz = time * ( exp( 2.0 * deta ) - 1.0 ) / ( exp( 2.0 * deta ) + 1.0 ); //translated to spatial coordinate z

//   cout << "deta=" << deta << "   dz=" << dz << endl;

  getJpsiFugacity( time, dr, dz, fugacity, n_jpsi, temp, enDen );
  printJpsiFugacity << fugacity;
  printJpsiFugacity << "\t";
  printJpsiFugacity << double(n_jpsi) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles();
  printJpsiFugacity << "\t";
  printJpsiFugacity << temp;
  printJpsiFugacity << "\t";
  printJpsiFugacity << enDen / theConfig->getTestparticles();



  printJpsiFugacity << endl;

}

void analysis::getJpsiFugacity( const double time, const double dr, const double dz, double& fugacity, int& n_jpsi, double& temp, double& enDen )
{
  double n_jpsi_equ;

  double V = M_PI * pow( dr, 2.0 ) * 2.0 * dz / pow( 0.197, 3.0 ); // 1/GeV^3

  n_jpsi = 0;

  double e_g = 0.0;
  int n_g = 0;



  for ( int i = 0; i < addedParticles.size(); i++ )
  {
    if ( addedParticles[i].T_creation <= time )
    {
      if (( addedParticles[i].Pos.Perp2() < pow( dr, 2.0 ) )  && ( fabs( addedParticles[i].Pos.Z() ) < dz ) )
      {
        // number of jpsi
        if ( addedParticles[i].FLAVOR == 50 )
        {
          n_jpsi++;
        }
      }
    }
  }
  
  for ( int i = 0; i < particles_atTimeNow.size(); i++ )
  {
    if ( particles_atTimeNow[i].Pos.T() <= time )
    {
      if (( particles_atTimeNow[i].Pos.Perp2() < pow( dr, 2.0 ) )  && ( fabs( particles_atTimeNow[i].Pos.Z() ) < dz ) )
      {
        // compute temp
        if ( particles_atTimeNow[i].FLAVOR == 0 ) // gluon
        {
          e_g += particles_atTimeNow[i].Mom.E();
          n_g++;
        }
      }
    }
  }

  temp = e_g /( 3.0 * n_g ); // GeV
  enDen = e_g /( V * pow( 0.197, 3.0 ) ); // GeV/fm^3


  
  n_jpsi_equ = theInterpolation_nJpsi.getN( temp ) * V;  // GeV^3/GeV^3 = 1


  fugacity = double(n_jpsi)/( n_jpsi_equ * theConfig->getTestparticles() * theConfig->getNaddedEvents() * theConfig->getJpsiTestparticles() );

//   cout << "t=" << time << "  temp=" << temp << "  deltaTemp=" << deltaTemp <<  "  n_jpsi_equ=" << n_jpsi_equ << "  fugacity=" << fugacity << "  n_jpsi=" << n_jpsi << "  V=" << V << endl;

}



void analysis::jpsi_correlations()
{
  double cos_delta_phi, cos_delta_theta, pt, delta_eta, eta, eta_charm;
  string filename;

  filename = filename_prefix + "_jpsi_corr_ptbins_iniJpsi";
  binning ptbins_iniJpsi(filename, 0.0, 6.0, 40);
  filename = filename_prefix + "_jpsi_corr_ptbins_secJpsi";
  binning ptbins_secJpsi(filename, 0.0, 6.0, 40);
  filename = filename_prefix + "_jpsi_corr_phibins_deltaPhiJpsiCharm";
  binning phibins_deltaPhiJpsiCharm(filename, 0.0, M_PI, 30);
  filename = filename_prefix + "_jpsi_corr_phibins_deltaPhiCharmCharm";
  binning phibins_deltaPhiCharmCharm(filename, 0.0, M_PI, 30);
  filename = filename_prefix + "_jpsi_corr_phibins_deltaThetaJpsiCharm";
  binning phibins_deltaThetaJpsiCharm(filename, 0.0, M_PI, 30);
  filename = filename_prefix + "_jpsi_corr_phibins_deltaThetaCharmCharm";
  binning phibins_deltaThetaCharmCharm(filename, 0.0, M_PI, 30);
  filename = filename_prefix + "_jpsi_corr_ptbins_charmFromJpsi";
  binning ptbins_charmFromJpsi(filename, 0.0, 6.0, 40);
  
  filename = filename_prefix + "_jpsi_corr_etabins_detaJpsiCharm";
  binning etabins_detaJpsiCharm(filename, 0.0, 20.0, 70);
  filename = filename_prefix + "_jpsi_corr_etabins_detaCharmCharm";
  binning etabins_detaCharmCharm(filename, 0.0, 20.0, 70);
  filename = filename_prefix + "_jpsi_corr_etabins_detaSecJpsiAllCharm";
  binning etabins_detaSecJpsiAllCharm(filename, 0.0, 20.0, 70);
  filename = filename_prefix + "_jpsi_corr_etabins_detaIniJpsiAllCharm";
  binning etabins_detaIniJpsiAllCharm(filename, 0.0, 20.0, 70);
  filename = filename_prefix + "_jpsi_corr_etabins_detaJpsiCharm_2d";
  binning2d etabins_detaJpsiCharm2d(filename, 0.0, 10.0, 40, 0.0, 10.0, 40);
  
  
  filename = filename_prefix + "_jpsi_corr_eta_test";
  binning eta_test(filename, -7.0, 7.0, 100);
  filename = filename_prefix + "_jpsi_corr_y_test";
  binning y_test(filename, -7.0, 7.0, 100);
  
  
  for ( int i = 0; i < addedParticles.size(); i++ )
  {
    if ( addedParticles[i].FLAVOR == 50 ) // jpsi
    {
      pt = addedParticles[i].Mom.Pt();
      eta = addedParticles[i].Mom.Pseudorapidity( addedParticles[i].m );
      
      eta_test.add(eta);
      y_test.add( addedParticles[i].Mom.Rapidity() );
      
      if( addedParticles[i].initially_produced ) // initial Jpsi
      {
        ptbins_iniJpsi.add(pt);
        
        for ( int j = 0; j < addedParticles.size(); j++ )
        {
          // all charm quarks
          if( addedParticles[j].FLAVOR == 7 || addedParticles[j].FLAVOR == 8 )
          {
            // delta eta
            eta_charm = addedParticles[j].Mom.Pseudorapidity( addedParticles[j].m );
            delta_eta = fabs( eta - eta_charm );
            etabins_detaIniJpsiAllCharm.add( delta_eta );
          }
        }
      }
      else
      {
        ptbins_secJpsi.add(pt);
        
        for ( int j = 0; j < addedParticles.size(); j++ )
        {
          // find partners of charm quarks in Jpsi
          if( ( addedParticles[j].N_EVENT_pp == addedParticles[i].N_EVENT_pp || addedParticles[j].N_EVENT_pp == addedParticles[i].N_EVENT_Cbar ) && i != j )
          {
//             if ( addedParticles[j].FLAVOR != 7 && addedParticles[j].FLAVOR != 8 )
//               cout << "error in analysis::jpsi_correlations() " << addedParticles[j].FLAVOR << "  " << addedParticles[i].FLAVOR << "  " << addedParticles[j].m << "  " << addedParticles[i].m << "  " << addedParticles[j].N_EVENT_pp << "  " << addedParticles[i].N_EVENT_pp << "  " << addedParticles[j].N_EVENT_Cbar << "  " << addedParticles[i].N_EVENT_Cbar << endl;
            
            cos_delta_phi = CosPhi( addedParticles[i].Mom, addedParticles[j].Mom );
            phibins_deltaPhiJpsiCharm.add( acos(cos_delta_phi) );
            
            cos_delta_theta = CosTheta( addedParticles[i].Mom, addedParticles[j].Mom );
            phibins_deltaThetaJpsiCharm.add( acos(cos_delta_theta) );
            
            // delta eta
            eta_charm = addedParticles[j].Mom.Pseudorapidity( addedParticles[j].m );
            delta_eta = fabs( eta - eta_charm );
            etabins_detaJpsiCharm.add( delta_eta );
            etabins_detaJpsiCharm2d.add( eta, eta_charm );
          }
          
          // all charm quarks
          if( addedParticles[j].FLAVOR == 7 || addedParticles[j].FLAVOR == 8 )
          {
            // delta eta
            eta_charm = addedParticles[j].Mom.Pseudorapidity( addedParticles[j].m );
            delta_eta = fabs( eta - eta_charm );
            etabins_detaSecJpsiAllCharm.add( delta_eta );
          }
        }  
      }
    }
    else if ( addedParticles[i].FLAVOR == 7 || addedParticles[i].FLAVOR == 8 )
    {
      if( addedParticles[i].jpsi_dissociation_number != -1 ) // produced in Jpsi dissociation
      {
        eta = addedParticles[i].Mom.Pseudorapidity( addedParticles[i].m );
        
        for ( int j = 0; j < addedParticles.size(); j++ )
        {
          if( addedParticles[j].jpsi_dissociation_number == addedParticles[i].jpsi_dissociation_number && i != j )
          {
            if ( addedParticles[j].FLAVOR != 7 && addedParticles[j].FLAVOR != 8 )
              cout << "error2 in analysis::jpsi_correlations() " << addedParticles[j].FLAVOR << "  " << addedParticles[j].m << "  " << addedParticles[i].FLAVOR << "  " << addedParticles[i].m << "  " << addedParticles[j].jpsi_dissociation_number << endl;
            
            cos_delta_phi = CosPhi( addedParticles[i].Mom, addedParticles[j].Mom );
            phibins_deltaPhiCharmCharm.add( acos(cos_delta_phi) );
      
            cos_delta_theta = CosTheta( addedParticles[i].Mom, addedParticles[j].Mom );
            phibins_deltaThetaCharmCharm.add( acos(cos_delta_theta) );
            
            // delta eta
            eta_charm = addedParticles[j].Mom.Pseudorapidity( addedParticles[j].m );
            delta_eta = fabs( eta - eta_charm );
            etabins_detaCharmCharm.add( delta_eta );
            
            ptbins_charmFromJpsi.add( addedParticles[j].Mom.Pt() );
            ptbins_charmFromJpsi.add( addedParticles[i].Mom.Pt() );
          }
        }
      }
    }
  }
  
  
  ptbins_iniJpsi.print();
  ptbins_secJpsi.print();
  phibins_deltaPhiJpsiCharm.print();
  phibins_deltaPhiCharmCharm.print();
  phibins_deltaThetaJpsiCharm.print();
  phibins_deltaThetaCharmCharm.print();
  ptbins_charmFromJpsi.print();
  etabins_detaJpsiCharm.print();
  etabins_detaCharmCharm.print();
  etabins_detaSecJpsiAllCharm.print();
  etabins_detaIniJpsiAllCharm.print();
  etabins_detaJpsiCharm2d.print();
  eta_test.print();
  y_test.print();
}



void analysis::ini_charm_correlations()
{
  double cos_delta_phi, cos_delta_theta, eta, eta_charm, delta_eta;
  string filename;


  filename = filename_prefix + "_phibins_deltaPhiIniCharm";
  binning phibins_deltaPhiIniCharm(filename, 0.0, M_PI, 30);
  filename = filename_prefix + "_phibins_deltaThetaIniCharm";
  binning phibins_deltaThetaIniCharm(filename, 0.0, M_PI, 30);
  filename = filename_prefix + "_etabins_detaIniCharm";
  binning etabins_detaIniCharm(filename, 0.0, 20.0, 70);

  for ( unsigned int i = 0; i < addedParticles.size(); i++ )
  {
    if ( addedParticles[i].FLAVOR == 7 || addedParticles[i].FLAVOR == 8 )
    {
      eta = addedParticles[i].Mom.Pseudorapidity( addedParticles[i].m );
      
      for ( unsigned int j = i+1; j < addedParticles.size(); j++ )
      {
        if( addedParticles[j].N_EVENT_pp == addedParticles[i].N_EVENT_pp &&  ( addedParticles[j].FLAVOR == 7 || addedParticles[j].FLAVOR == 8 ) )
        {
          cos_delta_phi = CosPhi( addedParticles[i].Mom, addedParticles[j].Mom );
          phibins_deltaPhiIniCharm.add( acos(cos_delta_phi) );
          
          cos_delta_theta = CosTheta( addedParticles[i].Mom, addedParticles[j].Mom );
          phibins_deltaThetaIniCharm.add( acos(cos_delta_theta) );
          
          // delta eta
          eta_charm = addedParticles[j].Mom.Pseudorapidity( addedParticles[j].m );
          delta_eta = fabs( eta - eta_charm );
          etabins_detaIniCharm.add( delta_eta );
        }
      }
    }
  }
  phibins_deltaPhiIniCharm.print();
  phibins_deltaThetaIniCharm.print();
  etabins_detaIniCharm.print();

}


void analysis::writeTempInTube( const int step  )
{
  const string sep = "  ";

  double temp, tempWithQuarks, energyDensity, fugacityGluons, fugacityQuarks, densityGluons, densityQuarks;
  int n_jpsi;
  double dr, dz, deta;

  string filename = filename_prefix + "_TempFugInTube" + ".dat";
  //filename = filename + "_spatial";
  
  fstream printTempInTube( filename.c_str(), ios::out | ios::app  );
  

  double time;
  if ( step == 0 )
    time = 0.0;
  else if ( step == nTimeSteps )
    return;
  else
    time = tstep[step-1];

  if ( step == 0 )
  {
    // file header
    //printTempInTube << "#temperature" << endl;
    printTempInTube << "# time\ttemp,\ttempWithQuarks,\tenergy dens.,\tfugGluons,\tfugacityQuarks,\tdensGluons,\tdensQuarks." << endl;
    return;
  }

  //printTempInTube.width( 10 );
  printTempInTube << time ;
  printTempInTube<< "\t";


  dr = 1.5; //fm
  deta = 0.5; // spacetime rapidty interval
  dz = time * ( exp( 2.0 * deta ) - 1.0 ) / ( exp( 2.0 * deta ) + 1.0 ); //translated to spatial coordinate z

//   cout << "deta=" << deta << "   dz=" << dz << endl;

  calculateTempInTube( time, dr, dz, temp, tempWithQuarks, energyDensity, fugacityGluons, fugacityQuarks, densityGluons, densityQuarks );
  printTempInTube << temp;
  printTempInTube << "\t";
  printTempInTube << tempWithQuarks;
  printTempInTube << "\t";
  printTempInTube << energyDensity;
  printTempInTube << "\t";
  printTempInTube << fugacityGluons;
  printTempInTube << "\t";
  printTempInTube << fugacityQuarks;
  printTempInTube << "\t";
  printTempInTube << densityGluons;
  printTempInTube << "\t";
  printTempInTube << densityQuarks;
  //printTempInTube.width( 20 );
  
/*
  dr = 5.0; //fm
  deta = 0.5; // spacetime rapidty interval
  dz = time * ( exp( 2.0 * deta ) - 1.0 ) / ( exp( 2.0 * deta ) + 1.0 ); //translated to spatial coordinate z

  calculateTempInTube( time, dr, dz, temp, tempWithQuarks, energyDensity, fugacityGluons, fugacityQuarks, densityGluons, densityQuarks );
  printTempInTube << temp;
  printTempInTube << "\t";
  printTempInTube << tempWithQuarks;
  printTempInTube << "\t";
  printTempInTube << energyDensity;
  printTempInTube << "\t";
  printTempInTube << fugacityGluons;
  printTempInTube << "\t";
  printTempInTube << fugacityQuarks;
  printTempInTube << "\t";
  printTempInTube << densityGluons;
  printTempInTube << "\t";
  printTempInTube << densityQuarks;
  printTempInTube.width( 20 );

  dr = 2.0; //fm
  deta = 1.0; // spacetime rapidty interval
  dz = time * ( exp( 2.0 * deta ) - 1.0 ) / ( exp( 2.0 * deta ) + 1.0 ); //translated to spatial coordinate z

  calculateTempInTube( time, dr, dz, temp, tempWithQuarks, energyDensity, fugacityGluons, fugacityQuarks, densityGluons, densityQuarks );
  printTempInTube << temp;
  printTempInTube << "\t";
  printTempInTube << tempWithQuarks;
  printTempInTube << "\t";
  printTempInTube << energyDensity;
  printTempInTube << "\t";
  printTempInTube << fugacityGluons;
  printTempInTube << "\t";
  printTempInTube << fugacityQuarks;
  printTempInTube << "\t";
  printTempInTube << densityGluons;
  printTempInTube << "\t";
  printTempInTube << densityQuarks;
  printTempInTube.width( 20 );


  dr = 5.0; //fm
  deta = 1.0; // spacetime rapidty interval
  dz = time * ( exp( 2.0 * deta ) - 1.0 ) / ( exp( 2.0 * deta ) + 1.0 ); //translated to spatial coordinate z

  calculateTempInTube( time, dr, dz, temp, tempWithQuarks, energyDensity, fugacityGluons, fugacityQuarks, densityGluons, densityQuarks );
  printTempInTube << temp;
  printTempInTube << "\t";
  printTempInTube << tempWithQuarks;
  printTempInTube << "\t";
  printTempInTube << energyDensity;
  printTempInTube << "\t";
  printTempInTube << fugacityGluons;
  printTempInTube << "\t";
  printTempInTube << fugacityQuarks;
  printTempInTube << "\t";
  printTempInTube << densityGluons;
  printTempInTube << "\t";
  printTempInTube << densityQuarks;
  printTempInTube.width( 20 );
*/


  printTempInTube << endl;

}

void analysis::writeCustomTube(const int step)
{
  const string sep = "  ";

  double totEnergy, EnergyG, EnergyQ, totalPT;
  int NumberG, NumberQ, totNumber ;
  int n_jpsi;
  double dr, dz, deta;
  double v2 = 0.0;
  double v2g,v2q;
  double v2_a,v2_b,v2_c;
  double v2g_a,v2g_b,v2g_c;
  double v2q_a,v2q_b,v2q_c;

  double rLarmor,rLarmor_a,rLarmor_b,rLarmor_c;
  double rLarmorg,rLarmorg_a,rLarmorg_b,rLarmorg_c;
  double rLarmorq,rLarmorq_a,rLarmorq_b,rLarmorq_c;

  double IsoX,IsoY,IsoZ;
  double IsoXg,IsoYg,IsoZg;
  double IsoXq,IsoYq,IsoZq;

  
  string filename = filename_prefix + "_EnergiesInTube" + ".dat";
  string filename2 = filename_prefix + "_thermalization" + ".dat";
  //filename = filename + "_spatial";
  
  fstream printTempInTube( filename.c_str(), ios::out | ios::app  );
  fstream printThermalization( filename2.c_str(), ios::out | ios::app  );

  double time;
  if ( step == 0 )
    time = 0.0;
  else if ( step == nTimeSteps )
    return;
  else
    time = tstep[step];

  if ( step == 0 )
  {
    // file header
    //printTempInTube << "#temperature" << endl;
    printThermalization << "#[TIME-1\t#IsoX-2\t#IsoY-3\t#IsoZ-4\t#IsoXq-5\t#IsoYq-6\t#IsoZq-7\t#IsoXg-8\t#IsoYg-9\t#IsoZg-10\t#v2Parton-11\t#v2Quark-12\t#v2Glue-13\t#rLarmor-14\t#rLarmorq-15\t#rLarmorg-16\t]#x#[r,eta:Input,#r10+eta1,#r100+eta1,#r100+eta20]" << endl;

    printTempInTube << "#time\tNumberQuarks\tNumberGluons\tNumberTotal\tEnergyQuarks\tEnergyGluons\tEnergyTotal\tmeanEnergy\tmeanEnergyQuarks\tmeanEnergyGluons" << endl;
    return;
  }

  //printTempInTube.width( 10 );
  printTempInTube << time ;
  printTempInTube<< "\t";
  printThermalization << time ;
  printThermalization << "\t";

  dr = radiusAnalysisTube; //fm
  deta = dEtaAnalysisTube; // spacetime rapidty interval
  dz = time * ( exp( 2.0 * deta ) - 1.0 ) / ( exp( 2.0 * deta ) + 1.0 ); //translated to spatial coordinate z
  double dzFull = time * ( exp( 2.0 * 20 ) - 1.0 ) / ( exp( 2.0 * 20 ) + 1.0 ); //translated to spatial coordinate z 
  double dzMid = time * ( exp( 2.0 * 1.0 ) - 1.0 ) / ( exp( 2.0 * 1.0 ) + 1.0 ); //translated to spatial coordinate z 
  
  
  
//   cout << "deta=" << deta << "   dz=" << dz << endl;
//   cout << "time in alaysis=" << time << endl;
//   calculateTubeCustom( time, 10, dzMid, totEnergy, EnergyG, EnergyQ, NumberG, NumberQ, totNumber, totalPT, IsoX_a, IsoY_a, IsoZ_a, v2_a, v2g_a, v2q_a, rLarmor_a, rLarmorg_a, rLarmorq_a );
//   calculateTubeCustom( time, 100, dzMid, totEnergy, EnergyG, EnergyQ, NumberG, NumberQ, totNumber, totalPT, IsoX_b, IsoY_b, IsoZ_b, v2_b, v2g_b, v2q_b, rLarmor_b, rLarmorg_b, rLarmorq_b );
//   calculateTubeCustom( time, 100, dzFull, totEnergy, EnergyG, EnergyQ, NumberG, NumberQ, totNumber, totalPT, IsoX_c, IsoY_c, IsoZ_c, v2_c, v2g_c, v2q_c, rLarmor_c, rLarmorg_c, rLarmorq_c );
  
  
  printTempInTube << NumberQ;
  printTempInTube << "\t";
  printTempInTube << NumberG;
  printTempInTube << "\t";
  printTempInTube << totNumber;
  printTempInTube << "\t";
  printTempInTube << EnergyQ;
  printTempInTube << "\t";
  printTempInTube << EnergyG;
  printTempInTube << "\t";
  printTempInTube << totEnergy;
  printTempInTube << "\t";
  printTempInTube << totEnergy/totNumber;  
  printTempInTube << "\t";
  printTempInTube << EnergyQ/NumberQ;  
  printTempInTube << "\t";
  printTempInTube << EnergyG/NumberG;  
  printTempInTube << "\t";
  printTempInTube << totalPT/totNumber;   
  printTempInTube << endl;

  calculateTubeCustom( time, dr, dz, totEnergy, EnergyG, EnergyQ, NumberG, NumberQ, totNumber, totalPT, IsoX, IsoY, IsoZ, IsoXq, IsoYq, IsoZq, IsoXg, IsoYg, IsoZg, v2, v2g, v2q, rLarmor, rLarmorg, rLarmorq );
  printThermalization << IsoX;
  printThermalization << "\t";
  printThermalization << IsoY;
  printThermalization << "\t";  
  printThermalization << IsoZ;
  printThermalization << "\t";  
  printThermalization << IsoXq;
  printThermalization << "\t";
  printThermalization << IsoYq;
  printThermalization << "\t";  
  printThermalization << IsoZq;
  printThermalization << "\t";  
  printThermalization << IsoXg;
  printThermalization << "\t";
  printThermalization << IsoYg;
  printThermalization << "\t";  
  printThermalization << IsoZg;
  printThermalization << "\t";  
  printThermalization << v2;
  printThermalization << "\t"; 
  printThermalization << v2q;
  printThermalization << "\t";  
  printThermalization << v2g;
  printThermalization << "\t";  
  printThermalization << rLarmor;
  printThermalization << "\t";
  printThermalization << rLarmorq;
  printThermalization << "\t";
  printThermalization << rLarmorg;
  printThermalization << "\t";
  
  calculateTubeCustom( time, 10, dzMid, totEnergy, EnergyG, EnergyQ, NumberG, NumberQ, totNumber, totalPT, IsoX, IsoY, IsoZ, IsoXq, IsoYq, IsoZq, IsoXg, IsoYg, IsoZg, v2, v2g, v2q, rLarmor, rLarmorg, rLarmorq );
  printThermalization << IsoX;
  printThermalization << "\t";
  printThermalization << IsoY;
  printThermalization << "\t";  
  printThermalization << IsoZ;
  printThermalization << "\t";  
  printThermalization << IsoXq;
  printThermalization << "\t";
  printThermalization << IsoYq;
  printThermalization << "\t";  
  printThermalization << IsoZq;
  printThermalization << "\t";  
  printThermalization << IsoXg;
  printThermalization << "\t";
  printThermalization << IsoYg;
  printThermalization << "\t";  
  printThermalization << IsoZg;
  printThermalization << "\t";  
  printThermalization << v2;
  printThermalization << "\t"; 
  printThermalization << v2q;
  printThermalization << "\t";  
  printThermalization << v2g;
  printThermalization << "\t"; 
  printThermalization << rLarmor;
  printThermalization << "\t";
  printThermalization << rLarmorq;
  printThermalization << "\t";
  printThermalization << rLarmorg;
  printThermalization << "\t";

  calculateTubeCustom( time, 100, dzMid, totEnergy, EnergyG, EnergyQ, NumberG, NumberQ, totNumber, totalPT, IsoX, IsoY, IsoZ, IsoXq, IsoYq, IsoZq, IsoXg, IsoYg, IsoZg, v2, v2g, v2q, rLarmor, rLarmorg, rLarmorq );
  printThermalization << IsoX;
  printThermalization << "\t";
  printThermalization << IsoY;
  printThermalization << "\t";  
  printThermalization << IsoZ;
  printThermalization << "\t";  
  printThermalization << IsoXq;
  printThermalization << "\t";
  printThermalization << IsoYq;
  printThermalization << "\t";  
  printThermalization << IsoZq;
  printThermalization << "\t";  
  printThermalization << IsoXg;
  printThermalization << "\t";
  printThermalization << IsoYg;
  printThermalization << "\t";  
  printThermalization << IsoZg;
  printThermalization << "\t";  
  printThermalization << v2;
  printThermalization << "\t"; 
  printThermalization << v2q;
  printThermalization << "\t";  
  printThermalization << v2g;
  printThermalization << "\t"; 
  printThermalization << rLarmor;
  printThermalization << "\t";
  printThermalization << rLarmorq;
  printThermalization << "\t";
  printThermalization << rLarmorg;
  printThermalization << "\t";
 
  calculateTubeCustom( time, 100, dzFull, totEnergy, EnergyG, EnergyQ, NumberG, NumberQ, totNumber, totalPT, IsoX, IsoY, IsoZ, IsoXq, IsoYq, IsoZq, IsoXg, IsoYg, IsoZg, v2, v2g, v2q, rLarmor, rLarmorg, rLarmorq );
  printThermalization << IsoX;
  printThermalization << "\t";
  printThermalization << IsoY;
  printThermalization << "\t";  
  printThermalization << IsoZ;
  printThermalization << "\t";  
  printThermalization << IsoXq;
  printThermalization << "\t";
  printThermalization << IsoYq;
  printThermalization << "\t";  
  printThermalization << IsoZq;
  printThermalization << "\t";  
  printThermalization << IsoXg;
  printThermalization << "\t";
  printThermalization << IsoYg;
  printThermalization << "\t";  
  printThermalization << IsoZg;
  printThermalization << "\t";  
  printThermalization << v2;
  printThermalization << "\t"; 
  printThermalization << v2q;
  printThermalization << "\t";  
  printThermalization << v2g;
  printThermalization << "\t";
  printThermalization << rLarmor;
  printThermalization << "\t";
  printThermalization << rLarmorq;
  printThermalization << "\t";
  printThermalization << rLarmorg;
  printThermalization << "\t";
  printThermalization << endl; 
  
/*  double precision=0.1;
  if( (((1.-precision)<IsoX*3.)&&(IsoX*3.<(1.+precision)))&&(((1.-precision)<IsoY*3.)&&(IsoY*3.<(1.+precision)))&&(((1.-precision)<IsoZ*3.)&&(IsoZ*3.<(1.+precision))) )
  {
    printThermalization << "1";
  }
  else
  {
    printThermalization << "0";
  }
  printThermalization << endl;
*/  

  stringstream ss;
  ss << time;
  
  string filename3 = filename_prefix + "_partonSpectrum_" + ss.str();
  fstream printPartonSpectrum( filename3.c_str(), ios::out | ios::trunc  );
  partonEnergies.print(printPartonSpectrum);
  
  string filename4 = filename_prefix + "_gluonSpectrum_" + ss.str();
  fstream printGluonSpectrum( filename4.c_str(), ios::out | ios::trunc  );
  gluonEnergies.print(printGluonSpectrum);  
  
  string filename5 = filename_prefix + "_quarkSpectrum_" + ss.str();
  fstream printQuarkSpectrum( filename5.c_str(), ios::out | ios::trunc );
  quarkEnergies.print(printQuarkSpectrum);  
  
  cout << "Medium Analyse in Custom Tube done for time " << time << endl;
  
}

void analysis::calculateTubeCustom(const double time, const double radius, const double dz, double& totalEnergy,double& totalEnergyGluons,double& totalEnergyQuarks, int & totalNumberGluons, int & totalNumberQuarks, int & totalNumber, double & totalPT, double& IsoX, double& IsoY, double& IsoZ,double& IsoXq, double& IsoYq, double& IsoZq,double& IsoXg, double& IsoYg, double& IsoZg, double& v2,double& v2g,double& v2q, double& rLarmor, double& rLarmorg , double& rLarmorq  )
{
  int cell_id;
  double pr, XT;
  
  const int minNmbTemp = 30; // minimum number of particles to calculate temperature from

  // total length in z direction
  const double zlength = dz*2.0;

  // volume
  dv = M_PI * pow( radius , 2.0 ) * zlength; // fm^3

  totalNumber=0; // number of all particles in cell
  totalNumberGluons = 0;
  totalNumberQuarks = 0;
  totalEnergy =0.;
  totalEnergyQuarks =0.;
  totalEnergyGluons =0.;
  totalPT=0.;
  v2=0.;
  v2g=0.;
  v2q=0.;
  double totIsoX=0.;
  double totIsoY=0.;
  double totIsoZ=0.;
  rLarmorg = 0.0;
  rLarmorq = 0.0;  
  rLarmor = 0.0;  
  IsoX = 0.;
  IsoY = 0.;
  IsoZ = 0.;
  IsoXg = 0.;
  IsoYg = 0.;
  IsoZg = 0.;
  IsoXq = 0.;
  IsoYq = 0.;
  IsoZq = 0.;  
  //
  //&& ( fabs( particles_atTimeNow[i].Mom.Rapidity() ) < 0.5 )
  // sum over all particles
  // && ( fabs( particles_atTimeNow[i].Mom.Rapidity() ) < 0.1 ) 
  for ( int i = 0; i < particles_atTimeNow.size(); i++ )
  {
    if( (  particles_atTimeNow[i].FLAVOR < 7 ) && FPT_COMP_LE( particles_atTimeNow[i].Pos.T(), time ) ) // only gluons and light quarks
    {
      if (( particles_atTimeNow[i].Pos.Perp2() < pow( radius, 2.0 ) )  && ( fabs( particles_atTimeNow[i].Pos.Z() ) < dz ))
      {
        totIsoX += pow(particles_atTimeNow[i].Mom.Px()/particles_atTimeNow[i].Mom.E(),2.0);
        totIsoY += pow(particles_atTimeNow[i].Mom.Py()/particles_atTimeNow[i].Mom.E(),2.0);
        totIsoZ += pow(particles_atTimeNow[i].Mom.Pz()/particles_atTimeNow[i].Mom.E(),2.0);
        ++totalNumber;
        totalEnergy += particles_atTimeNow[i].Mom.E();
        totalPT += particles_atTimeNow[i].Mom.Pt();
        partonEnergies.add(particles_atTimeNow[i].Mom.E());
        v2 += particles_atTimeNow[i].Mom.FlowV2();
        rLarmor +=  0.197*sqrt(pow(particles_atTimeNow[i].Mom.Pz(),2.0)+pow(particles_atTimeNow[i].Mom.Px(),2.0))/(theConfig->getUsedExternalField()*pow(0.14,2.0));
        if( particles_atTimeNow[i].FLAVOR == gluon )
        {
          IsoXg += pow(particles_atTimeNow[i].Mom.Px()/particles_atTimeNow[i].Mom.E(),2.0);
          IsoYg += pow(particles_atTimeNow[i].Mom.Py()/particles_atTimeNow[i].Mom.E(),2.0);
          IsoZg += pow(particles_atTimeNow[i].Mom.Pz()/particles_atTimeNow[i].Mom.E(),2.0);
          rLarmorg +=  0.197*sqrt(pow(particles_atTimeNow[i].Mom.Pz(),2.0)+pow(particles_atTimeNow[i].Mom.Px(),2.0))/(theConfig->getUsedExternalField()*pow(0.14,2.0));
          totalEnergyGluons += particles_atTimeNow[i].Mom.E();
          totalNumberGluons++;
          gluonEnergies.add(particles_atTimeNow[i].Mom.E());
          v2g+= particles_atTimeNow[i].Mom.FlowV2();
        }
        else
        {
          IsoXq += pow(particles_atTimeNow[i].Mom.Px()/particles_atTimeNow[i].Mom.E(),2.0);
          IsoYq += pow(particles_atTimeNow[i].Mom.Py()/particles_atTimeNow[i].Mom.E(),2.0);
          IsoZq += pow(particles_atTimeNow[i].Mom.Pz()/particles_atTimeNow[i].Mom.E(),2.0);
          rLarmorq +=  0.197*sqrt(pow(particles_atTimeNow[i].Mom.Pz(),2.0)+pow(particles_atTimeNow[i].Mom.Px(),2.0))/(theConfig->getUsedExternalField()*pow(0.14,2.0));
          totalEnergyQuarks += particles_atTimeNow[i].Mom.E();
          totalNumberQuarks++;
          quarkEnergies.add(particles_atTimeNow[i].Mom.E());
          v2q +=  particles_atTimeNow[i].Mom.FlowV2();
        }
      }
    }
  }
  IsoX = totIsoX/totalNumber;
  IsoY = totIsoY/totalNumber;
  IsoZ = totIsoZ/totalNumber;
  IsoXg /= totalNumberGluons;
  IsoYg /= totalNumberGluons;
  IsoZg /= totalNumberGluons;
  IsoXq /= totalNumberQuarks;
  IsoYq /= totalNumberQuarks;
  IsoZq /= totalNumberQuarks;
  v2  /= totalNumber;
  v2g /= totalNumberGluons;
  v2q /= totalNumberQuarks;
  rLarmor  /= totalNumber;
  rLarmorg /= totalNumberGluons;
  rLarmorq /= totalNumberQuarks;
}


// writes temperature and velocity of all specified cells in a file, used by Alex Meistrenko as input
void analysis::calculateTempInTube( const double time, const double radius, const double dz, double & temp, double & tempWithQuarks, double & energyDensity, double & fugacityGluons, double & fugacityQuarks, double & densityGluons, double & densityQuarks  )
{  
  int cell_id;
  double pr, XT;
  
  const int minNmbTemp = 30; // minimum number of particles to calculate temperature from


  // total length in z direction
  const double zlength = dz*2.0;

  // volume
  dv = M_PI * pow( radius , 2.0 ) * zlength; // fm^3

  
  // the following routine is written for several cells -> here we do not need this, just set nCells = 1
  
  nCells = 1;

  numberInCell = new int[nCells]; // number of all particles in cell
  vx_cell = new double[nCells]; // total x-velocity of all particles in cell
  vy_cell = new double[nCells];  
  vz_cell = new double[nCells];
  vr_cell = new double[nCells];
  em_cell = new double[nCells];// total energy of all particles in cell
  prm_cell = new double[nCells];
  pzm_cell = new double[nCells];
  pr2em_cell = new double[nCells];
  pz2em_cell = new double[nCells];
  przem_cell = new double[nCells];
  densn_cell = new double[nCells];
  gama_cell = new double[nCells];
  temp_cell = new double[nCells];
  tempWithQuarks_cell = new double[nCells];
  
  int numberGluons = 0;
  int numberQuarks = 0;
  
  double numberGluonsEquil, numberQuarksEquil;
  


  // set all properties to 0
  for ( int i = 0; i < nCells; i++ )
  {
    numberInCell[i] = 0;
    vx_cell[i] = 0.0; 
    vy_cell[i] = 0.0;  
    vz_cell[i] = 0.0;
    vr_cell[i] = 0.0;
    em_cell[i] = 0.0;
    prm_cell[i] = 0.0;
    pzm_cell[i] = 0.0;
    pr2em_cell[i] = 0.0;
    pz2em_cell[i] = 0.0;
    przem_cell[i] = 0.0;
    densn_cell[i] = 0.0;
    gama_cell[i] = 0.0;
    temp_cell[i] = 0.0;
    tempWithQuarks_cell[i] = 0.0;
  }
  

  // sum over all particles
  for ( int i = 0; i < particles_atTimeNow.size(); i++ )
  {
    if ( FPT_COMP_LE( particles_atTimeNow[i].Pos.T(), time )  && particles_atTimeNow[i].FLAVOR < 7 ) // only gluons and light quarks
    {
      if (( particles_atTimeNow[i].Pos.Perp2() < pow( radius, 2.0 ) )  && ( fabs( particles_atTimeNow[i].Pos.Z() ) < dz ) )
      {
        // set cell id to 0
        cell_id = 0;
        
        ++numberInCell[cell_id];
        if( particles_atTimeNow[i].FLAVOR == gluon )
          numberGluons++;
        else
          numberQuarks++;
        
        XT = particles_atTimeNow[i].Pos.Perp();
        if ( XT < 1.0e-5 )
        {
          pr = particles_atTimeNow[i].Mom.Pt();
        }
        else
        {
          pr = ( particles_atTimeNow[i].Mom.Px() * particles_atTimeNow[i].Pos.X()
                 + particles_atTimeNow[i].Mom.Py() * particles_atTimeNow[i].Pos.Y() ) / XT;
        }
        vr_cell[cell_id] += pr / particles_atTimeNow[i].Mom.E();
        vx_cell[cell_id] += particles_atTimeNow[i].Mom.Px() / particles_atTimeNow[i].Mom.E();
        vy_cell[cell_id] += particles_atTimeNow[i].Mom.Py() / particles_atTimeNow[i].Mom.E();
        vz_cell[cell_id] += particles_atTimeNow[i].Mom.Pz() / particles_atTimeNow[i].Mom.E();
        
        em_cell[cell_id] += particles_atTimeNow[i].Mom.E();
        prm_cell[cell_id] += pr;
        pzm_cell[cell_id] += particles_atTimeNow[i].Mom.Pz();
        pr2em_cell[cell_id] += pr * pr / particles_atTimeNow[i].Mom.E();
        pz2em_cell[cell_id] += particles_atTimeNow[i].Mom.Pz() * particles_atTimeNow[i].Mom.Pz() / particles_atTimeNow[i].Mom.E();
        przem_cell[cell_id] += pr * particles_atTimeNow[i].Mom.Pz() / particles_atTimeNow[i].Mom.E();
      }
    }
  }

  // calculate temperature for cell
  for ( int i = 0; i < nCells; i++ )
  {
    // If more than minNmbTemp particles are in the cell the temperature is calculated by them
    if(numberInCell[i] >= minNmbTemp)
    {
      calcTempCell( i );
    }
    else // take also 6 neighbor cells into account
    {
      cout << "error in calculateTempInTube(): to few particles in cell. Number=" << numberInCell[i] << "  time=" << time << endl;
      cout << particles_atTimeNow.size() << "  " << particles_atTimeNow[0].Pos.T() << "  " << particles_atTimeNow[10].Pos.T() << endl;
    }
  }
  
  
  temp = temp_cell[0];
  tempWithQuarks = tempWithQuarks_cell[0];
  energyDensity = em_cell[0]; // GeV/fm^3
  
  numberGluonsEquil = 16.0 * pow( temp , 3.0 ) / pow( M_PI , 2.0) * dv / pow( 0.197 , 3.0 ) * theConfig->getTestparticles();
  numberQuarksEquil = 12.0 * Particle::N_light_flavor * pow( temp , 3.0 ) / pow( M_PI , 2.0) * dv / pow( 0.197 , 3.0 ) * theConfig->getTestparticles();
  
  fugacityGluons = double(numberGluons) / numberGluonsEquil;
  fugacityQuarks = double(numberQuarks) / numberQuarksEquil;
  
  densityGluons = double(numberGluons) / dv / theConfig->getTestparticles(); // 1/fm^3
  densityQuarks = double(numberQuarks) / dv / theConfig->getTestparticles(); // 1/fm^3

  delete[] numberInCell; 
  delete[] vx_cell; 
  delete[] vy_cell;  
  delete[] vz_cell;
  delete[] vr_cell;
  delete[] em_cell;
  delete[] prm_cell;
  delete[] pzm_cell;
  delete[] pr2em_cell;
  delete[] pz2em_cell;
  delete[] przem_cell;
  delete[] densn_cell;
  delete[] gama_cell;
  delete[] temp_cell;
  delete[] tempWithQuarks_cell;
}


// writes spatial profile of photons
void analysis::writePhotonSpaceProfile( const int step  )
{
  string filename,name;
  stringstream ss;
  double time,time_from_step_before;
  
  
  if ( step == 1 )
  {  
    time = tstep[step-1];
    time_from_step_before = 0.0;
  }else if ( step > 1 && step != nTimeSteps )
  {
    time = tstep[step-1];
    time_from_step_before = tstep[step-2];
  }
  else
    return;
  
  ss << step;
  filename = filename_prefix + "_phSpatial_" + ss.str() + ".dat";
  //filename = filename + "_spatial";
  
  fstream file_spatial( filename.c_str(), ios::out | ios::trunc  );
  
  int IXY = IX * IY;
  // total length of grid system
  const double xlength = 11;
  const double ylength = 11;
  // number of cells in given direction
  // in cascade IX=40   IY=40   IZ=47
  const int nCellsx = 200;
  const int nCellsy = 200;  
  const int zlength = 47;
  const int nCellsz = 41;  
  int cell_id,cell_id2D,nx,ny,nz;
  const int nCells = nCellsx * nCellsy * nCellsz;
  const int nCells2D = nCellsx * nCellsy;
  double dv = xlength * ylength * zlength / nCells; // volume of each cell
  double dx = xlength / nCellsx;
  double dy = ylength / nCellsy;
  double v2 = 0.;
  
  numberInCell = new int[nCells]; // number of all particles in cell
  numberInCell2D = new int[nCells2D]; // number of all particles in cell
  photonNumberInCell2D = new int[nCells2D]; // number of all photons in cell
  photonV2SumInCell2D = new double[nCells2D]; // Sum of v2 in Cell
  
    // set all properties to 0
  for ( int i = 0; i < nCells; i++ )
  {
    numberInCell[i] = 0;
  }
  for ( int i = 0; i < nCells2D; i++ )
  {
    numberInCell2D[i] = 0;
    photonNumberInCell2D[i]=0;
    photonV2SumInCell2D[i]=0.0;
  }
  
  //******************************
  //1.) Background Particle Number
  //******************************
  /*for ( int i = 0; i < particles_atTimeNow.size(); i++ )
  {
    // determine cell id
    if ( fabs( particles_atTimeNow[i].Pos.X() - xlength / 2.0 ) < 1.0e-6 )
      nx = nCellsx - 1;
    else
      nx = int(( particles_atTimeNow[i].Pos.X() / xlength + 0.5 ) * nCellsx );

    if ( fabs( particles_atTimeNow[i].Pos.Y() - ylength / 2.0 ) < 1.0e-6 )
      ny = nCellsy - 1;
    else
      ny = int(( particles_atTimeNow[i].Pos.Y() / ylength + 0.5 ) * nCellsy );

    if ( fabs( particles_atTimeNow[i].Pos.Z() - zlength / 2.0 ) < 1.0e-6 )
      nz = nCellsz - 1;
    else
      nz = int(( particles_atTimeNow[i].Pos.Z() / zlength + 0.5 ) * nCellsz );


    if (( nx >= nCellsx ) || ( nx < 0 ) || ( ny >= nCellsy ) || ( ny < 0 ) || ( nz >= nCellsz ) || ( nz < 0 ) )
    {
      cout << "err cell_ID in temp output" << endl;
      cout << particles_atTimeNow[i].Pos.T() << "\t" << particles_atTimeNow[i].Pos.X() << "\t" << particles_atTimeNow[i].Pos.Y();
      cout << "\t" << particles_atTimeNow[i].Pos.Z() << endl;
      cout << nx << "\t" << ny << "\t" << nz << endl;
    }
    else
    {
      cell_id = nx + nCellsx * ny + nCellsx * nCellsy * nz;
      cell_id2D = nx + nCellsx * ny;
      ++numberInCell[cell_id];
      ++numberInCell2D[cell_id2D];
    }
  }
  */
  
  //******************************
  //2.) Photon Number
  //******************************
  for ( int i = 0; i < noninteractingParticles.size(); i++ )
  {
    if (noninteractingParticles[i].Mom.Rapidity() > 0.35 )
    {
      continue;
    }
    if (noninteractingParticles[i].production_time < time_from_step_before || noninteractingParticles[i].production_time > time )
    {
      continue;
    }    
    // determine cell id
    if ( fabs( noninteractingParticles[i].Pos.X() - xlength / 2.0 ) < 1.0e-6 )
      nx = nCellsx - 1;
    else
      nx = int(( noninteractingParticles[i].Pos.X() / xlength + 0.5 ) * nCellsx );

    if ( fabs( noninteractingParticles[i].Pos.Y() - ylength / 2.0 ) < 1.0e-6 )
      ny = nCellsy - 1;
    else
      ny = int(( noninteractingParticles[i].Pos.Y() / ylength + 0.5 ) * nCellsy );

    //if ( fabs( noninteractingParticles[i].Pos.Z() - zlength / 2.0 ) < 1.0e-6 )
    //  nz = nCellsz - 1;
    //else
    //  nz = int(( noninteractingParticles[i].Pos.Z() / zlength + 0.5 ) * nCellsz );


    if (( nx >= nCellsx ) || ( nx < 0 ) || ( ny >= nCellsy ) || ( ny < 0 ) )
    {
      cout << "err cell_ID in temp output" << endl;
      cout << particles_atTimeNow[i].Pos.T() << "\t" << particles_atTimeNow[i].Pos.X() << "\t" << particles_atTimeNow[i].Pos.Y();
      cout << "\t" << particles_atTimeNow[i].Pos.Z() << endl;
      cout << nx << "\t" << ny << "\t" << nz << endl;
    }
    else
    {
      cell_id2D = nx + nCellsx * ny;
      v2 = (pow(noninteractingParticles[i].Mom.Px(),2.0)-pow(noninteractingParticles[i].Mom.Py(),2.0))/pow(noninteractingParticles[i].Mom.Pt(),2.0);
      photonV2SumInCell2D[cell_id2D] += v2;
      ++photonNumberInCell2D[cell_id2D];
    }
  }
  
  file_spatial << "#x\t#y\t#PhNum\t#PhV2"<< endl;
  for(int i = 0; i < nCellsx;i++)
  {
    for(int j = 0; j < nCellsy; j++ )
    {
      int cellID2Dhere = i + nCellsx * j ;
      if(photonNumberInCell2D[cellID2Dhere]>0 && photonV2SumInCell2D[cellID2Dhere]>0)
      {
        file_spatial << dx*i << "\t" << dy*j << "\t" << photonNumberInCell2D[cellID2Dhere] << "\t" << 1.0 << endl;
      }else if (photonNumberInCell2D[cellID2Dhere]>0 && photonV2SumInCell2D[cellID2Dhere]<0)
      {
        file_spatial << dx*i << "\t" << dy*j << "\t" << photonNumberInCell2D[cellID2Dhere] << "\t" << -1.0 << endl;        
      }else
      {
        file_spatial << dx*i << "\t" << dy*j << "\t" << photonNumberInCell2D[cellID2Dhere] << "\t" << "0" << endl;          
      }
    }
  }
  file_spatial.close();
}

// writes temperature and velocity of all specified cells in a file, used by Moritz
void analysis::writeTempCustom( const int step  )
{
  
  int nx, ny, nz, cell_id;
  double time, pr, XT;
 
  const int minNmbTemp = 50; // minimum number of particles to calculate temperature from
//   const int minNmbTemp = 2; // minimum number of particles to calculate temperature from
  const int minNmbTempCell = 10; // minimum number of particles in one cell to calculate a Temperature. If the number is below this value, the cell is taken as empty (no temperature)
  
  binning binNumber("output/number.dat", -0.5, 49.5, 50);
  binning binTempCells("output/tempCells.dat", 0.0, 3., 70);
  binning binTempRings("output/tempRings.dat", 0.0, 2.2, 70);
  
  int count_tempCell = 0;
  int count_tempNeighborCells = 0;
  int count_noTemp = 0;
  
  int IXY = IX * IY;

  if ( step != 0 && step != nTimeSteps )
    time = tstep[step-1];
  else
    return;
  
 
  // total length of grid system
  const double xlength = 25;//24.6;
  const double ylength = 25;//24.6;
  const double zlength = 13;//12.3;

  // number of cells in given direction
  // in cascade IX=40   IY=40   IZ=47
  const int nCellsx = 10;
  const int nCellsy = 10;
  const int nCellsz = 50;
  nCells = nCellsx * nCellsy * nCellsz;
  
  dv = xlength * ylength * zlength / nCells; // volume of each cell

//   int numberInCell[nCells]; // number of all particles in cell
//   double vx_cell[nCells]; // total x-velocity of all particles in cell
//   double vy_cell[nCells];  
//   double vz_cell[nCells];
//   double vr_cell[nCells];
//   double em_cell[nCells];// total energy of all particles in cell
//   double prm_cell[nCells];
//   double pzm_cell[nCells];
//   double pr2em_cell[nCells];
//   double pz2em_cell[nCells];
//   double przem_cell[nCells];
//   double densn_cell[nCells];
//   double gama_cell[nCells];
//   double temp_cell[nCells];
//   double tempWithQuarks_cell[nCells];


  numberInCell = new int[nCells]; // number of all particles in cell
  numberInCellQuarks = new int[nCells]; // number of all particles in cell
  numberInCellGluons = new int[nCells]; // number of all particles in cell
  fugacityGluons = new double[nCells]; 
  fugacityQuarks = new double[nCells]; 
  
  
  temp_numberInCell = new int[nCells]; // number of all particles in cell
  vx_cell = new double[nCells]; // total x-velocity of all particles in cell
  vy_cell = new double[nCells];  
  vz_cell = new double[nCells];
  vr_cell = new double[nCells];
  em_cell = new double[nCells];// total energy of all particles in cell
  prm_cell = new double[nCells];
  pzm_cell = new double[nCells];
  pr2em_cell = new double[nCells];
  pz2em_cell = new double[nCells];
  przem_cell = new double[nCells];
  densn_cell = new double[nCells];
  gama_cell = new double[nCells];
  temp_cell = new double[nCells];
  tempWithQuarks_cell = new double[nCells];
  
  // for copy. If too few particles in one cell, the particles from the surrounding cells are added from _org. Otherwise one would double or triple add particles from cells
  numberInCell_org = new int[nCells]; // number of all particles in cell
  vx_cell_org = new double[nCells]; // total x-velocity of all particles in cell
  vy_cell_org = new double[nCells];  
  vz_cell_org = new double[nCells];
  vr_cell_org = new double[nCells];
  em_cell_org = new double[nCells];// total energy of all particles in cell
  prm_cell_org = new double[nCells];
  pzm_cell_org = new double[nCells];
  pr2em_cell_org = new double[nCells];
  pz2em_cell_org = new double[nCells];
  przem_cell_org = new double[nCells];
// 

//   for ( int i = 0; i < nCells; i++ )
//   {
//     const int nxny = nCellsx * nCellsy;
//     const int indexZ = i / nxny;
//     const int indexY = (i -  indexZ  * nxny) / nCellsx;
//     const int indexX = i -  indexZ  * nxny - indexY * nCellsx;
//     
//     energy[i] = 0.0;
//     numberInCell[i] = indexX + nCellsx * indexY + nCellsx * nCellsy * indexZ;
//     vx[i] = double(indexX) * xlength / nCellsx;
//     vy[i] = double(indexY) * ylength / nCellsy;
//     vz[i] = double(indexZ) * zlength / nCellsz;
//   }

  // set all properties to 0
  for ( int i = 0; i < nCells; i++ )
  {
    numberInCell[i] = 0;
    numberInCellQuarks[i] = 0;
    numberInCellGluons[i] = 0;
    fugacityQuarks[i]=0;
    fugacityGluons[i]=0;
    temp_numberInCell[i] = 0;
    vx_cell[i] = 0.0; 
    vy_cell[i] = 0.0;  
    vz_cell[i] = 0.0;
    vr_cell[i] = 0.0;
    em_cell[i] = 0.0;
    prm_cell[i] = 0.0;
    pzm_cell[i] = 0.0;
    pr2em_cell[i] = 0.0;
    pz2em_cell[i] = 0.0;
    przem_cell[i] = 0.0;
    densn_cell[i] = 0.0;
    gama_cell[i] = 0.0;
    temp_cell[i] = 0.0;
    tempWithQuarks_cell[i] = 0.0;
    
    numberInCell_org[i] = 0;
    vx_cell_org[i] = 0.0; 
    vy_cell_org[i] = 0.0;  
    vz_cell_org[i] = 0.0;
    vr_cell_org[i] = 0.0;
    em_cell_org[i] = 0.0;
    prm_cell_org[i] = 0.0;
    pzm_cell_org[i] = 0.0;
    pr2em_cell_org[i] = 0.0;
    pz2em_cell_org[i] = 0.0;
    przem_cell_org[i] = 0.0;
  }
  
  // sum over all particles
  for ( int i = 0; i < particles_atTimeNow.size(); i++ )
  {
    if ( FPT_COMP_E( particles_atTimeNow[i].Pos.T(), time ) && particles_atTimeNow[i].FLAVOR < 7 ) // only gluons and light quarks
    {
      // determine cell id
      if ( fabs( particles_atTimeNow[i].Pos.X() - xlength / 2.0 ) < 1.0e-6 )
        nx = nCellsx - 1;
      else
        nx = int(( particles_atTimeNow[i].Pos.X() / xlength + 0.5 ) * nCellsx );

      if ( fabs( particles_atTimeNow[i].Pos.Y() - ylength / 2.0 ) < 1.0e-6 )
        ny = nCellsy - 1;
      else
        ny = int(( particles_atTimeNow[i].Pos.Y() / ylength + 0.5 ) * nCellsy );

      if ( fabs( particles_atTimeNow[i].Pos.Z() - zlength / 2.0 ) < 1.0e-6 )
        nz = nCellsz - 1;
      else
        nz = int(( particles_atTimeNow[i].Pos.Z() / zlength + 0.5 ) * nCellsz );


      if (( nx >= nCellsx ) || ( nx < 0 ) || ( ny >= nCellsy ) || ( ny < 0 ) || ( nz >= nCellsz ) || ( nz < 0 ) )
      {
        /*cout << "err cell_ID in temp output" << endl;
        cout << particles_atTimeNow[i].Pos.T() << "\t" << particles_atTimeNow[i].Pos.X() << "\t" << particles_atTimeNow[i].Pos.Y();
        cout << "\t" << particles_atTimeNow[i].Pos.Z() << endl;
        cout << nx << "\t" << ny << "\t" << nz << endl;*/
      }
      else
      {
        cell_id = nx + nCellsx * ny + nCellsx * nCellsy * nz;
        
        if(particles_atTimeNow[i].FLAVOR != gluon )
        {
          ++numberInCellQuarks[cell_id];
        }
        if(particles_atTimeNow[i].FLAVOR == gluon)
        {
          ++numberInCellGluons[cell_id];
        }        
        ++numberInCell[cell_id];
        
        XT = particles_atTimeNow[i].Pos.Perp();
        if ( XT < 1.0e-5 )
        {
          pr = particles_atTimeNow[i].Mom.Pt();
        }
        else
        {
          pr = ( particles_atTimeNow[i].Mom.Px() * particles_atTimeNow[i].Pos.X()
                 + particles_atTimeNow[i].Mom.Py() * particles_atTimeNow[i].Pos.Y() ) / XT;
        }
        vr_cell[cell_id] += pr / particles_atTimeNow[i].Mom.E();
        vx_cell[cell_id] += particles_atTimeNow[i].Mom.Px() / particles_atTimeNow[i].Mom.E();
        vy_cell[cell_id] += particles_atTimeNow[i].Mom.Py() / particles_atTimeNow[i].Mom.E();
        vz_cell[cell_id] += particles_atTimeNow[i].Mom.Pz() / particles_atTimeNow[i].Mom.E();
        
        em_cell[cell_id] += particles_atTimeNow[i].Mom.E();
        prm_cell[cell_id] += pr;
        pzm_cell[cell_id] += particles_atTimeNow[i].Mom.Pz();
        pr2em_cell[cell_id] += pr * pr / particles_atTimeNow[i].Mom.E();
        pz2em_cell[cell_id] += particles_atTimeNow[i].Mom.Pz() * particles_atTimeNow[i].Mom.Pz() / particles_atTimeNow[i].Mom.E();
        przem_cell[cell_id] += pr * particles_atTimeNow[i].Mom.Pz() / particles_atTimeNow[i].Mom.E();
        
        // temp of particles summed
        tempWithQuarks_cell[cell_id] += particles_atTimeNow[i].temperature;
        if( particles_atTimeNow[i].temperature >= 0.1)
          ++temp_numberInCell[cell_id];
      }
    }
  }
  
  // duplicate properties of cell
  for ( int i = 0; i < nCells; i++ )
  {
    numberInCell_org[i] = numberInCell[i];
    vr_cell_org[i] = vr_cell[i];
    vx_cell_org[i] = vx_cell[i];
    vy_cell_org[i] = vy_cell[i];
    vz_cell_org[i] = vz_cell[i];
    
    em_cell_org[i] = em_cell[i];
    prm_cell_org[i] = prm_cell[i];
    pzm_cell_org[i] = pzm_cell[i];
    pr2em_cell_org[i] = pr2em_cell[i];
    pz2em_cell_org[i] = pz2em_cell[i];
    przem_cell_org[i] = przem_cell[i];
  }

  
  // determine cells which do not have enough particles
  for ( int i = 0; i < nCells; i++ )
  {
    // If less than minNmbTempCell particles are in the cell it is taken as empty
    if(numberInCell[i] >= minNmbTempCell)
    {
      // If less than minNmbTemp particles are in the cell the temperature the 6 neighbor cells are also taken into account
      if(numberInCell[i] < minNmbTemp)
      {
        int cell_id_neighbor;
        
        //neighbors in x direction
        cell_id_neighbor = i-1;
        addNeighborCells( i, cell_id_neighbor );
        cell_id_neighbor = i+1;
        addNeighborCells( i, cell_id_neighbor );
        //neighbors in y direction
        cell_id_neighbor = i+IX;
        addNeighborCells( i, cell_id_neighbor );
        cell_id_neighbor = i-IX;
        addNeighborCells( i, cell_id_neighbor );
        //neighbors in z direction
        cell_id_neighbor = i-IXY;
        addNeighborCells( i, cell_id_neighbor );
        cell_id_neighbor = i+IXY;
        addNeighborCells( i, cell_id_neighbor );
      }
    }
  }
  
  // calculate temperature for each cell
  for ( int i = 0; i < nCells; i++ )
  {
    // If less than minNmbTempCell particles are in the cell it is taken as empty
    if(numberInCell[i] >= minNmbTempCell)
    {
      // If more than minNmbTemp particles are in the cell the temperature is calculated by them
      if(numberInCell[i] >= minNmbTemp)
      {
        calcTempCell( i );
        
        tempWithQuarks_cell[i] = tempWithQuarks_cell[i]/temp_numberInCell[i];

        binTempCells.add(temp_cell[i]);
        count_tempCell++;
      }
      else // take also 6 neighbor cells into account
      {
        count_noTemp++;
        
        temp_cell[i] = 0.0; // not enough particles in surrounding to calculate temperature
        
        vx_cell[i] = 0.0;
        vy_cell[i] = 0.0;
        vz_cell[i] = 0.0;
      }
    }
    else
    {
      if(numberInCell[i] == 0)
      {
        temp_cell[i] = -2.0; // completely empty
        vx_cell[i] = 0.0;
        vy_cell[i] = 0.0;
        vz_cell[i] = 0.0;
      }
      else
      {
        temp_cell[i] = -1.0; // very few particles
        vx_cell[i] = 0.0;
        vy_cell[i] = 0.0;
        vz_cell[i] = 0.0;
      }
    }
  }

  double avFQ = 0.;
  double avFG = 0.;
  double avT = 0.;
  for(int i=1;i<=nCells;i++)
  {
    avFQ += fugacityQuarks[i];  
    avFG += fugacityGluons[i];
    avT += temp_cell[i];
  }
  avFQ /=count_tempCell; 
  avFG /=count_tempCell;
  avT /=count_tempCell;
  
  
  string filename,name;
  stringstream ss;
  ss << time*10;
  //+ ss.str() 
  filename = filename_prefix + "_central" + ".dat";
  fstream file_central( filename.c_str(), ios::out | ios::app  );
  
  const int middle_cell_id = nCellsx/2 + nCellsx * nCellsy/2 +nCellsx*nCellsy*nCellsz/2;
  cout  << time <<'\t'<< "T   " << '\t' << "FUgQ   " << '\t'<< "FugG  " << '\t' << "nQUARKS  " << '\t'<< "nGluons  "  << '\t' << "nTot   " << endl;
  cout  << time <<'\t'<< temp_cell[middle_cell_id] << '\t' << fugacityQuarks[middle_cell_id] << '\t'<<fugacityGluons[middle_cell_id] << '\t' << numberInCellQuarks[middle_cell_id] << '\t'<< numberInCellGluons[middle_cell_id]<< '\t' << numberInCell[middle_cell_id] << '\t'<< numberInCellQuarks[middle_cell_id]+numberInCellGluons[middle_cell_id]  << endl;

  file_central << time <<'\t'<<avT<<'\t'<<avFQ<<'\t'<<avFG<<'\t'<< temp_cell[middle_cell_id] << '\t' << fugacityQuarks[middle_cell_id] << '\t'<<fugacityGluons[middle_cell_id] << '\t' << numberInCellQuarks[middle_cell_id] << '\t'<< numberInCellGluons[middle_cell_id]<< '\t' << numberInCell[middle_cell_id] << endl;

  
  delete[] numberInCell; 
  delete[] vx_cell; 
  delete[] vy_cell;  
  delete[] vz_cell;
  delete[] vr_cell;
  delete[] em_cell;
  delete[] prm_cell;
  delete[] pzm_cell;
  delete[] pr2em_cell;
  delete[] pz2em_cell;
  delete[] przem_cell;
  delete[] densn_cell;
  delete[] gama_cell;
  delete[] temp_cell;
  delete[] tempWithQuarks_cell;
}



// writes temperature and velocity of all specified cells in a file, used by Alex Meistrenko as input
void analysis::writeTempAndVel( const int step  )
{
  int nx, ny, nz, cell_id;
  double time, pr, XT;
 
  const int minNmbTemp = 30; // minimum number of particles to calculate temperature from
//   const int minNmbTemp = 2; // minimum number of particles to calculate temperature from
  const int minNmbTempCell = 2; // minimum number of particles in one cell to calculate a Temperature. If the number is below this value, the cell is taken as empty (no temperature)
  
  binning binNumber("output/number.dat", -0.5, 49.5, 50);
  binning binTempCells("output/tempCells.dat", 0.0, 3., 70);
  binning binTempRings("output/tempRings.dat", 0.0, 2.2, 70);
  
  int count_tempCell = 0;
  int count_tempNeighborCells = 0;
  int count_noTemp = 0;
  
  int IXY = IX * IY;

  if ( step != 0 && step != nTimeSteps )
    time = tstep[step-1];
  else
    return;

  // total length of grid system
  const double xlength = 25;//24.6;
  const double ylength = 25;//24.6;
  const double zlength = 25;//12.3;

  // number of cells in given direction
  // in cascade IX=40   IY=40   IZ=47
  const int nCellsx = 10;
  const int nCellsy = 10;
  const int nCellsz = 10;
  nCells = nCellsx * nCellsy * nCellsz;
  
  dv = xlength * ylength * zlength / nCells; // volume of each cell

//   int numberInCell[nCells]; // number of all particles in cell
//   double vx_cell[nCells]; // total x-velocity of all particles in cell
//   double vy_cell[nCells];  
//   double vz_cell[nCells];
//   double vr_cell[nCells];
//   double em_cell[nCells];// total energy of all particles in cell
//   double prm_cell[nCells];
//   double pzm_cell[nCells];
//   double pr2em_cell[nCells];
//   double pz2em_cell[nCells];
//   double przem_cell[nCells];
//   double densn_cell[nCells];
//   double gama_cell[nCells];
//   double temp_cell[nCells];
//   double tempWithQuarks_cell[nCells];


  numberInCell = new int[nCells]; // number of all particles in cell
  temp_numberInCell = new int[nCells]; // number of all particles in cell
  vx_cell = new double[nCells]; // total x-velocity of all particles in cell
  vy_cell = new double[nCells];  
  vz_cell = new double[nCells];
  vr_cell = new double[nCells];
  em_cell = new double[nCells];// total energy of all particles in cell
  prm_cell = new double[nCells];
  pzm_cell = new double[nCells];
  pr2em_cell = new double[nCells];
  pz2em_cell = new double[nCells];
  przem_cell = new double[nCells];
  densn_cell = new double[nCells];
  gama_cell = new double[nCells];
  temp_cell = new double[nCells];
  tempWithQuarks_cell = new double[nCells];
  
  // for copy. If too few particles in one cell, the particles from the surrounding cells are added from _org. Otherwise one would double or triple add particles from cells
  numberInCell_org = new int[nCells]; // number of all particles in cell
  vx_cell_org = new double[nCells]; // total x-velocity of all particles in cell
  vy_cell_org = new double[nCells];  
  vz_cell_org = new double[nCells];
  vr_cell_org = new double[nCells];
  em_cell_org = new double[nCells];// total energy of all particles in cell
  prm_cell_org = new double[nCells];
  pzm_cell_org = new double[nCells];
  pr2em_cell_org = new double[nCells];
  pz2em_cell_org = new double[nCells];
  przem_cell_org = new double[nCells];
// 

//   for ( int i = 0; i < nCells; i++ )
//   {
//     const int nxny = nCellsx * nCellsy;
//     const int indexZ = i / nxny;
//     const int indexY = (i -  indexZ  * nxny) / nCellsx;
//     const int indexX = i -  indexZ  * nxny - indexY * nCellsx;
//     
//     energy[i] = 0.0;
//     numberInCell[i] = indexX + nCellsx * indexY + nCellsx * nCellsy * indexZ;
//     vx[i] = double(indexX) * xlength / nCellsx;
//     vy[i] = double(indexY) * ylength / nCellsy;
//     vz[i] = double(indexZ) * zlength / nCellsz;
//   }

  // set all properties to 0
  for ( int i = 0; i < nCells; i++ )
  {
    numberInCell[i] = 0;
    temp_numberInCell[i] = 0;
    vx_cell[i] = 0.0; 
    vy_cell[i] = 0.0;  
    vz_cell[i] = 0.0;
    vr_cell[i] = 0.0;
    em_cell[i] = 0.0;
    prm_cell[i] = 0.0;
    pzm_cell[i] = 0.0;
    pr2em_cell[i] = 0.0;
    pz2em_cell[i] = 0.0;
    przem_cell[i] = 0.0;
    densn_cell[i] = 0.0;
    gama_cell[i] = 0.0;
    temp_cell[i] = 0.0;
    tempWithQuarks_cell[i] = 0.0;
    
    numberInCell_org[i] = 0;
    vx_cell_org[i] = 0.0; 
    vy_cell_org[i] = 0.0;  
    vz_cell_org[i] = 0.0;
    vr_cell_org[i] = 0.0;
    em_cell_org[i] = 0.0;
    prm_cell_org[i] = 0.0;
    pzm_cell_org[i] = 0.0;
    pr2em_cell_org[i] = 0.0;
    pz2em_cell_org[i] = 0.0;
    przem_cell_org[i] = 0.0;
  }
  
  // sum over all particles
  for ( int i = 0; i < particles_atTimeNow.size(); i++ )
  {
    if ( FPT_COMP_E( particles_atTimeNow[i].Pos.T(), time ) && particles_atTimeNow[i].FLAVOR < 7 ) // only gluons and light quarks
    {
      // determine cell id
      if ( fabs( particles_atTimeNow[i].Pos.X() - xlength / 2.0 ) < 1.0e-6 )
        nx = nCellsx - 1;
      else
        nx = int(( particles_atTimeNow[i].Pos.X() / xlength + 0.5 ) * nCellsx );

      if ( fabs( particles_atTimeNow[i].Pos.Y() - ylength / 2.0 ) < 1.0e-6 )
        ny = nCellsy - 1;
      else
        ny = int(( particles_atTimeNow[i].Pos.Y() / ylength + 0.5 ) * nCellsy );

      if ( fabs( particles_atTimeNow[i].Pos.Z() - zlength / 2.0 ) < 1.0e-6 )
        nz = nCellsz - 1;
      else
        nz = int(( particles_atTimeNow[i].Pos.Z() / zlength + 0.5 ) * nCellsz );


      if (( nx >= nCellsx ) || ( nx < 0 ) || ( ny >= nCellsy ) || ( ny < 0 ) || ( nz >= nCellsz ) || ( nz < 0 ) )
      {
        cout << "err cell_ID in temp output" << endl;
        cout << particles_atTimeNow[i].Pos.T() << "\t" << particles_atTimeNow[i].Pos.X() << "\t" << particles_atTimeNow[i].Pos.Y();
        cout << "\t" << particles_atTimeNow[i].Pos.Z() << endl;
        cout << nx << "\t" << ny << "\t" << nz << endl;
      }
      else
      {
        cell_id = nx + nCellsx * ny + nCellsx * nCellsy * nz;
        
        ++numberInCell[cell_id];
        
        XT = particles_atTimeNow[i].Pos.Perp();
        if ( XT < 1.0e-5 )
        {
          pr = particles_atTimeNow[i].Mom.Pt();
        }
        else
        {
          pr = ( particles_atTimeNow[i].Mom.Px() * particles_atTimeNow[i].Pos.X()
                 + particles_atTimeNow[i].Mom.Py() * particles_atTimeNow[i].Pos.Y() ) / XT;
        }
        vr_cell[cell_id] += pr / particles_atTimeNow[i].Mom.E();
        vx_cell[cell_id] += particles_atTimeNow[i].Mom.Px() / particles_atTimeNow[i].Mom.E();
        vy_cell[cell_id] += particles_atTimeNow[i].Mom.Py() / particles_atTimeNow[i].Mom.E();
        vz_cell[cell_id] += particles_atTimeNow[i].Mom.Pz() / particles_atTimeNow[i].Mom.E();
        
        em_cell[cell_id] += particles_atTimeNow[i].Mom.E();
        prm_cell[cell_id] += pr;
        pzm_cell[cell_id] += particles_atTimeNow[i].Mom.Pz();
        pr2em_cell[cell_id] += pr * pr / particles_atTimeNow[i].Mom.E();
        pz2em_cell[cell_id] += particles_atTimeNow[i].Mom.Pz() * particles_atTimeNow[i].Mom.Pz() / particles_atTimeNow[i].Mom.E();
        przem_cell[cell_id] += pr * particles_atTimeNow[i].Mom.Pz() / particles_atTimeNow[i].Mom.E();
        
        // temp of particles summed
        tempWithQuarks_cell[cell_id] += particles_atTimeNow[i].temperature;
        if( particles_atTimeNow[i].temperature >= 0.1)
          ++temp_numberInCell[cell_id];
      }
    }
  }
  
  // duplicate properties of cell
  for ( int i = 0; i < nCells; i++ )
  {
    numberInCell_org[i] = numberInCell[i];
    vr_cell_org[i] = vr_cell[i];
    vx_cell_org[i] = vx_cell[i];
    vy_cell_org[i] = vy_cell[i];
    vz_cell_org[i] = vz_cell[i];
    
    em_cell_org[i] = em_cell[i];
    prm_cell_org[i] = prm_cell[i];
    pzm_cell_org[i] = pzm_cell[i];
    pr2em_cell_org[i] = pr2em_cell[i];
    pz2em_cell_org[i] = pz2em_cell[i];
    przem_cell_org[i] = przem_cell[i];
  }

  
  // determine cells which do not have enough particles
  for ( int i = 0; i < nCells; i++ )
  {
    // If less than minNmbTempCell particles are in the cell it is taken as empty
    if(numberInCell[i] >= minNmbTempCell)
    {
      // If less than minNmbTemp particles are in the cell the temperature the 6 neighbor cells are also taken into account
      if(numberInCell[i] < minNmbTemp)
      {
        int cell_id_neighbor;
        
        //neighbors in x direction
        cell_id_neighbor = i-1;
        addNeighborCells( i, cell_id_neighbor );
        cell_id_neighbor = i+1;
        addNeighborCells( i, cell_id_neighbor );
        //neighbors in y direction
        cell_id_neighbor = i+IX;
        addNeighborCells( i, cell_id_neighbor );
        cell_id_neighbor = i-IX;
        addNeighborCells( i, cell_id_neighbor );
        //neighbors in z direction
        cell_id_neighbor = i-IXY;
        addNeighborCells( i, cell_id_neighbor );
        cell_id_neighbor = i+IXY;
        addNeighborCells( i, cell_id_neighbor );
      }
    }
  }
  
  // calculate temperature for each cell
  for ( int i = 0; i < nCells; i++ )
  {
    // If less than minNmbTempCell particles are in the cell it is taken as empty
    if(numberInCell[i] >= minNmbTempCell)
    {
      // If more than minNmbTemp particles are in the cell the temperature is calculated by them
      if(numberInCell[i] >= minNmbTemp)
      {
        calcTempCell( i );
        
        tempWithQuarks_cell[i] = tempWithQuarks_cell[i]/temp_numberInCell[i];

        binTempCells.add(temp_cell[i]);
        count_tempCell++;
      }
      else // take also 6 neighbor cells into account
      {
        count_noTemp++;
        
        temp_cell[i] = 0.0; // not enough particles in surrounding to calculate temperature
        
        vx_cell[i] = 0.0;
        vy_cell[i] = 0.0;
        vz_cell[i] = 0.0;
      }
    }
    else
    {
      if(numberInCell[i] == 0)
      {
        temp_cell[i] = -2.0; // completely empty
        vx_cell[i] = 0.0;
        vy_cell[i] = 0.0;
        vz_cell[i] = 0.0;
      }
      else
      {
        temp_cell[i] = -1.0; // very few particles
        vx_cell[i] = 0.0;
        vy_cell[i] = 0.0;
        vz_cell[i] = 0.0;
      }
    }
  }


  string filename,name;
  stringstream ss;
  ss << time*10;
  // 1) Cell Info
  filename = filename_prefix + "_tempVel_" + ss.str() + ".dat";

  fstream file( filename.c_str(), ios::out | ios::trunc  );

  if ( step == 1 )
  {
    file << "#temperature and velocity" << endl;
    file << "#simulation parameter:" << endl;
    file << "#time = " << time << endl;
    file << "#testparticles= " << theConfig->getTestparticles() << endl;
    file << "#runtime= " << theConfig->getRuntime() << endl;
    file << "#sqrtS= " << theConfig->getSqrtS() << " GeV" << endl;
    file << "#b= " << theConfig->getImpactParameter() << " fm" << endl;
    file << "#(" << theConfig->getA()  << ") on ("
    << theConfig->getA()  << ")" << endl;
    file << "# T in GeV/fm^3" << endl;
    file << "#" << endl;
    file << "#" << "T" << "\t" << "vx" << "\t" << "vy" << "\t" << "vz" << endl;
  }
  
  for ( int i = 0; i < nCells; i++ )
  {
    file << temp_cell[i] << "\t";
    file << tempWithQuarks_cell[i] << "\t";
    file << vx_cell[i] << "\t";
    file << vy_cell[i] << "\t";
    file << vz_cell[i] << endl;
  }
  
  // 2) Spatial info about cells
  filename = filename + "_spatial";
  
  fstream file_spatial( filename.c_str(), ios::out | ios::trunc  );
  
  for ( int i = 0; i < nCells; i++ )
  {
    const int nxny = nCellsx * nCellsy;
    const int indexZ = i / nxny; //Integer
    const int indexY = (i -  indexZ  * nxny) / nCellsx;
    const int indexX = i -  indexZ  * nxny - indexY * nCellsx;
    
    if(indexY == int(nCellsy/2))
    {
//       file_spatial << i << "\t";
      file_spatial << double(indexZ) * zlength / nCellsz << "\t";
      file_spatial << double(indexX) * xlength / nCellsx << "\t";

      file_spatial << temp_cell[i] << "\t";
      file_spatial << tempWithQuarks_cell[i] << "\t";
      file_spatial << vx_cell[i] << "\t";
      file_spatial << vy_cell[i] << "\t";
      file_spatial << vz_cell[i] << endl;
    }
  }

  // 2) Central info 
  filename = filename + "_central";
  fstream file_central( filename.c_str(), ios::out | ios::trunc  );
  
  const int middle_cell_id = nCellsx/2 + nCellsx * nCellsy/2 +nCellsx*nCellsy*nCellsz/2;
  file_central << temp_cell[middle_cell_id] << endl;

  
  delete[] numberInCell; 
  delete[] vx_cell; 
  delete[] vy_cell;  
  delete[] vz_cell;
  delete[] vr_cell;
  delete[] em_cell;
  delete[] prm_cell;
  delete[] pzm_cell;
  delete[] pr2em_cell;
  delete[] pz2em_cell;
  delete[] przem_cell;
  delete[] densn_cell;
  delete[] gama_cell;
  delete[] temp_cell;
  delete[] tempWithQuarks_cell;
}

void analysis::calcTempCell( const int cell_id ) 
{
  gama_cell[cell_id] = 1.0;
  vx_cell[cell_id] = vx_cell[cell_id] / numberInCell[cell_id];
  vy_cell[cell_id] = vy_cell[cell_id] / numberInCell[cell_id];
  vz_cell[cell_id] = vz_cell[cell_id] / numberInCell[cell_id];
  vr_cell[cell_id] = vr_cell[cell_id] / numberInCell[cell_id];
  if ( vr_cell[cell_id] > 0.0 ) gama_cell[cell_id] = 1.0 / sqrt( 1.0 - vz_cell[cell_id] * vz_cell[cell_id] - vr_cell[cell_id] * vr_cell[cell_id] );
  else
  {
    vr_cell[cell_id] = 0.0;
    gama_cell[cell_id] = 1.0 / sqrt( 1.0 - vz_cell[cell_id] * vz_cell[cell_id] );
  }
  
  densn_cell[cell_id] = double(numberInCell[cell_id]) / dv / theConfig->getTestparticles() / gama_cell[cell_id];//1/fm^3
  em_cell[cell_id] = ( em_cell[cell_id] - 2.0 * vr_cell[cell_id] * prm_cell[cell_id] - 2.0 * vz_cell[cell_id] * pzm_cell[cell_id] + vr_cell[cell_id] * vr_cell[cell_id] * pr2em_cell[cell_id]
  + vz_cell[cell_id] * vz_cell[cell_id] * pz2em_cell[cell_id] + 2.0 * vr_cell[cell_id] * vz_cell[cell_id] * przem_cell[cell_id] )
  / theConfig->getTestparticles() / dv * gama_cell[cell_id] * gama_cell[cell_id];//GeV/fm^3
  
  temp_cell[cell_id] = em_cell[cell_id] / ( 3.0 * densn_cell[cell_id] );
  
  double neqQ = 3.0*2.0*3.0*2.0/pow(M_PI,2.0)*pow(temp_cell[cell_id],3.0)/(pow(0.197,3.0));
  double neqG = 16.0/pow(M_PI,2.0)*pow(temp_cell[cell_id],3.0)/(pow(0.197,3.0));
  fugacityQuarks[cell_id] = (double(numberInCellQuarks[cell_id])/ dv/ theConfig->getTestparticles()* gama_cell[cell_id] )/neqQ;
  fugacityGluons[cell_id] = (double(numberInCellGluons[cell_id])/ dv/ theConfig->getTestparticles()* gama_cell[cell_id] )/neqG;
  // assume thermal equilibrium and additional quark flavor* gama_cell[cell_id]
//   int Nflavor_temp = 3;
//   tempWithQuarks_cell[cell_id] = pow(pi*pi/3.0 / (16.+12.*Nflavor_temp) * em_cell[cell_id]  * pow(0.197,3.0) , 1.0/4.0);
}

void analysis::addNeighborCells( const int cell_id, const int neighborCell_id ) 
{
  if(neighborCell_id >= 0 && neighborCell_id < nCells)
  {
    numberInCell[cell_id] += numberInCell_org[neighborCell_id];
    vr_cell[cell_id] += vr_cell_org[neighborCell_id];
    vx_cell[cell_id] += vx_cell_org[neighborCell_id];
    vy_cell[cell_id] += vy_cell_org[neighborCell_id];
    vz_cell[cell_id] += vz_cell_org[neighborCell_id];
    
    em_cell[cell_id] += em_cell_org[neighborCell_id];
    prm_cell[cell_id] += prm_cell_org[neighborCell_id];
    pzm_cell[cell_id] += pzm_cell_org[neighborCell_id];
    pr2em_cell[cell_id] += pr2em_cell_org[neighborCell_id];
    pz2em_cell[cell_id] += pz2em_cell_org[neighborCell_id];
    przem_cell[cell_id] += przem_cell_org[neighborCell_id];
  }
}






// print dndy, detdy and mean pt of gluons
void analysis::print_dndy(const string subfix )
{
  double y, pt;
  double pt_sum = 0.0;
  
  string filename;

  filename = filename_prefix + "_dndy_" + subfix;
  binning dndy(filename, -10.0, 10.0, 200);
  filename = filename_prefix + "_dedy_" + subfix;
  binningValues dedy(filename, -10.0, 10.0, 200);
  
  filename = filename_prefix + "_dndpt_" + subfix;
  binning dndpt(filename, 0.0, 20.0, 200);
  
  for ( unsigned int i = 0; i < particles_atTimeNow.size(); i++ )
  {
    pt = particles_atTimeNow[i].Mom.Pt();
    y  = particles_atTimeNow[i].Mom.Rapidity();
    
//     if(y < -8.0 || isnan(y))
// //       cout << y << "\t" << particles_atTimeNow[i].Mom.E() << "\t" << particles_atTimeNow[i].Mom.Pz()  << "\t" << pt << "\t" << i << endl;
//       cout << i << "\t";
    
    dndy.add(y);
    dedy.add(y,pt);
    dndpt.add(pt);
    pt_sum += pt;
  }
  
  cout << "mean pt = " << pt_sum/particles_atTimeNow.size() << endl;
  
  dndy.print();
  dedy.print();
  dndpt.print();

}











void analysis::analyseAngleDe()
{
  double costheta; // angle between D meson and e-
  double cosphi; // transverse angle between D meson and e-
  double ept;
  int k_e;
  
  string filename;
  
  
  filename = filename_prefix + "_dmeson_electron_cosTheta";
  binningValues binsCosTheta(filename, 0.0, 10.0, 200);
  filename = filename_prefix + "_dmeson_electron_cosPhiTrans";
  binningValues binsCosPhi(filename, 0.0, 10.0, 200);
  
  for ( unsigned int i = 0; i < addedParticles.size(); i++ )
  {
    if ( addedParticles[i].FLAVOR == 7 || addedParticles[i].FLAVOR == 8 )
    {
      // there are several electrons per dmeson/charm quark
      for(int k = 0; k < theConfig->getNumberElectronStat(); k++)
      {      
        k_e = ( i ) * theConfig->getNumberElectronStat() + k ;
        costheta = CosTheta( addedParticlesCopy[i].Mom, addedPartcl_electron[k_e].Mom );
        
        cosphi = CosPhi( addedParticlesCopy[i].Mom, addedPartcl_electron[k_e].Mom );

        ept = addedPartcl_electron[k_e].Mom.Pt(); // e- pt
        
        binsCosTheta.add(ept, costheta);
        binsCosPhi.add(ept, cosphi);
      }
    }
  }
  
  binsCosTheta.print();
  binsCosPhi.print();
  
}


// void analysis::analyseCharmTestJetEvolution(const int step)
// {
//   double energy_sum = 0.0;
//   int count = 0;
//   for ( int i = 0; i < addedParticles.size(); i++ )
//   {
//     if ( addedParticles[i].FLAVOR == 7 )
//     {
//       count++;
//       energy_sum += addedParticles[i].E;
//     }
//   }
// 
//   charmJetEnergy[step] = energy_sum/count;
//   timestepAnalysed[step] = true;
// }
// 
// void analysis::analyseCharmTestJet()
// {
//   string filename;
//   
//   filename = filename_prefix + "_Ebins";
//   binning Ebins(filename, 0.0, 12.0, 24);
//   filename = filename_prefix + "_xbins";
//   binning xbins(filename, -3.0, 5.0, 70);
//   filename = filename_prefix + "_ybins";
//   binning ybins(filename, -3.0, 3.0, 70);
//   filename = filename_prefix + "_elossbins";
//   binning elossbins(filename, -2.0, 10.0, 70);
//   
//   for ( int i = 0; i < addedParticles.size(); i++ )
//   {
//     if ( addedParticles[i].FLAVOR == 7 )
//     {
//       Ebins.add(addedParticles[i].E);
//       xbins.add(addedParticles[i].X);
//       ybins.add(addedParticles[i].Y);
//       elossbins.add(10.0 - addedParticles[i].E);
//     }
//   }
//   
//   Ebins.print();
//   xbins.print();
//   ybins.print();
//   elossbins.print();
//   
//   filename = filename_prefix + "_jetEvolution";
//   fstream print( filename.c_str(), ios::out | ios::trunc );
// 
//   print << "# time   #mean jet energy" << endl;
//   for ( int i = 0; i < nTimeSteps + 1; i++ )
//   {
//     if ( timestepAnalysed[i] )
//     {
//       print.width( 15 );
//       if ( i == 0 )
//         print << "0";
//       else
//         print << tstep[i-1];
//       print.width( 15 );
//       print << charmJetEnergy[i];
//       print << endl;
//     }
//   }
// 
// }



void analysis::jpsiEvolution( int step )
{
  double pt_min = 0;
  
  if( theConfig->getOutputScheme() == cms_jpsi )
    pt_min = 6.5;
  
  
  int countJpsi_all = 0;
  int countJpsi_ini = 0;
//   int countJpsi_midNormRap = 0;
  int countJpsi_midPseudoRap_all = 0;
  int countJpsi_midPseudoRap_ini = 0;
  int countJpsi_forwardPseudoRap_all = 0;
  int countJpsi_forwardPseudoRap_ini = 0;
  int countJpsi_midNormRap_all = 0;
  int countJpsi_midNormRap_ini = 0;
  int countJpsi_forwardNormRap_all = 0;
  int countJpsi_forwardNormRap_ini = 0;
//   int countJpsi_midSpaceTimeRap = 0;

  for ( unsigned int i = 0; i < addedParticles.size(); i++ )
  {
    if ( addedParticles[i].FLAVOR == 50 )
    {
      double pt = addedParticles[i].Mom.Pt();
      
      countJpsi_all++;
        
      if(addedParticles[i].initially_produced)
        countJpsi_ini++;
      
      if( pt >= pt_min )
      {
        double pseudorap = addedParticles[i].Mom.Pseudorapidity( addedParticles[i].m );
        
        if ( fabs( pseudorap ) >= rapidityRanges[0].yleft && fabs( pseudorap ) <= rapidityRanges[0].yright )
        {
          countJpsi_midPseudoRap_all++;
          if(addedParticles[i].initially_produced)
            countJpsi_midPseudoRap_ini++;
        }
        
        if ( fabs( pseudorap ) >= rapidityRanges[1].yleft && fabs( pseudorap ) <= rapidityRanges[1].yright )
        {
          countJpsi_forwardPseudoRap_all++;
          if(addedParticles[i].initially_produced)
            countJpsi_forwardPseudoRap_ini++;
        }
        
              
        // normal rapidity
        double normrap = addedParticles[i].Mom.Rapidity();
        
        if ( fabs( normrap ) >= rapidityRanges[0].yleft && fabs( normrap ) <= rapidityRanges[0].yright )
        {
          countJpsi_midNormRap_all++;
          if(addedParticles[i].initially_produced)
            countJpsi_midNormRap_ini++;
        }
        
        
        if ( fabs( normrap ) >= rapidityRanges[1].yleft && fabs( normrap ) <= rapidityRanges[1].yright )
        {
          countJpsi_forwardNormRap_all++;
          if(addedParticles[i].initially_produced)
            countJpsi_forwardNormRap_ini++;
        }
          
        
  //       // space time rapidity
  //       double strap = addedParticles[i].Pos.Rapidity();
  //       if( fabs(strap) < midrap )
  //         countJpsi_midSpaceTimeRap++;
      }
    }
  }


  numberJpsi_all_time[step] = countJpsi_all;
  numberJpsi_ini_time[step] = countJpsi_ini;
  numberJpsi_midPseudoRap_all_time[step] = countJpsi_midPseudoRap_all;
  numberJpsi_midPseudoRap_ini_time[step] = countJpsi_midPseudoRap_ini;
  numberJpsi_forwardPseudoRap_all_time[step] = countJpsi_forwardPseudoRap_all;
  numberJpsi_forwardPseudoRap_ini_time[step] = countJpsi_forwardPseudoRap_ini;
  numberJpsi_midNormRap_all_time[step] = countJpsi_midNormRap_all;
  numberJpsi_midNormRap_ini_time[step] = countJpsi_midNormRap_ini;
  numberJpsi_forwardNormRap_all_time[step] = countJpsi_forwardNormRap_all;
  numberJpsi_forwardNormRap_ini_time[step] = countJpsi_forwardNormRap_ini;
//   numberJpsi_midSpaceTimeRap_time[step] = countJpsi_midSpaceTimeRap;
  numberJpsiProd_time[step] = ns_heavy_quarks::jpsicreation;
  numberJpsiDiss_time[step] = ns_heavy_quarks::jpsi_dissociation;
  numberJpsiDissTd_time[step] = ns_heavy_quarks::jpsi_dissociation_from_temperature;
  numberCCbGG_time[step] = ns_heavy_quarks::charmAnnihil;
  timestepAnalysed[step] = true;
}

void analysis::printJpsiEvolution()
{
  string filename = filename_prefix + "_jpsiEvolution";
  fstream print_je( filename.c_str(), ios::out | ios::trunc );
  
  const double delta_y_mid = 2.0 * ( rapidityRanges[0].yright - rapidityRanges[0].yleft );
  const double delta_y_forward = 2.0 * ( rapidityRanges[1].yright - rapidityRanges[1].yleft );

  print_je << "# charm Annihaltion and J/psi processes" << endl;
  print_je << "# time   #J/psis    #J/psi production     # ccb->gg processes" << endl;
  for ( int i = 0; i < nTimeSteps + 1; i++ )
  {
    if ( timestepAnalysed[i] )
    {
      print_je.width( 15 );
      if ( i == 0 )
        print_je << "0";
      else
        print_je << tstep[i-1];
      print_je.width( 15 );
      print_je << double( numberJpsi_all_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles();
      print_je.width( 15 );
      print_je << double( numberJpsi_ini_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles();
      print_je.width( 15 );
      print_je << double( numberJpsiProd_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles();
      print_je.width( 15 );
      print_je << double( numberJpsiDiss_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles();
      print_je.width( 15 );
      print_je << double( numberJpsiDissTd_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles();
      print_je.width( 15 );
      print_je << double( numberCCbGG_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles();

      
      //  / (2.0*1.0) because the rapidity range is 1.0, but for + and . Therefore, times 2.
      print_je.width( 15 );
      print_je << double( numberJpsi_midPseudoRap_all_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles() / delta_y_mid;
      print_je.width( 15 );
      print_je << double( numberJpsi_midPseudoRap_ini_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles() / delta_y_mid;
      print_je.width( 15 );
      print_je << double( numberJpsi_forwardPseudoRap_all_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles() / delta_y_forward;
      print_je.width( 15 );
      print_je << double( numberJpsi_forwardPseudoRap_ini_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles() / delta_y_forward;
      
      print_je.width( 15 );
      print_je << double( numberJpsi_midNormRap_all_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles() / delta_y_mid;
      print_je.width( 15 );
      print_je << double( numberJpsi_midNormRap_ini_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles() / delta_y_mid;
      print_je.width( 15 );
      print_je << double( numberJpsi_forwardNormRap_all_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles() / delta_y_forward;
      print_je.width( 15 );
      print_je << double( numberJpsi_forwardNormRap_ini_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles() / delta_y_forward;
      
//       print_je.width( 15 );
//       print_je << double( numberJpsi_midSpaceTimeRap_time[i] ) / theConfig->getTestparticles() / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles();
      print_je << endl;
    }
  }
}


void analysis::registerProgressInformationForOutput( const double _time, const double _dt, const int _nAddedParticles, const int _nMediumParticles, const int _nColl, const int _nColl22, const int _nColl23, const int _nColl32 )
{
  if( progressLogFile.good() )
  {
    string sep = "\t";
    progressLogFile << _time << sep;
    progressLogFile << _dt << sep;
    progressLogFile << _nAddedParticles << sep;
    progressLogFile << _nMediumParticles << sep;
    progressLogFile << _nColl << sep;
    progressLogFile << _nColl22 << sep;
    progressLogFile << _nColl23 << sep;
    progressLogFile << _nColl32 << endl;
  }
  else
  {
    string errMsg = "error when attempting to write progress information log";
    throw eAnalysis_error( errMsg );
  }
}


void analysis::scatteredMediumParticlesOutput( const int step )
{
  string name;
  stringstream ss;

  if( step == 0 )
    name = "initial";
  else
    if( step == nTimeSteps )
      name = "final";
    else
    {
      ss << step;
      name = "step" + ss.str();
    }

  //creates filename, for example: "./output/run32_step1.f1", "./output/run32_initial.f1" or the likes
  string filename = filename_prefix + "_" + name + "_scatteredMediumParticles.f1";
  fstream file( filename.c_str(), ios::out | ios::trunc );

  //---- print header if file is empty ----
  time_t end;
  time( &end );

  file.seekp( 0, ios::end );
  long size = file.tellp();
  file.seekp( 0, ios::beg );
  if( size == 0 )
    printHeader( file, all, end );
  //---------------------------------------

  for( int i = 0; i < scatteredMediumParticles.size(); i++ )
  {
    file << i << sep << scatteredMediumParticles[i].unique_id << sep << scatteredMediumParticles[i].cell_id 
         << sep << scatteredMediumParticles[i].FLAVOR << sep 
         << scatteredMediumParticles[i].Pos.T() << sep << scatteredMediumParticles[i].Pos.X() << sep
         << scatteredMediumParticles[i].Pos.Y() << sep << scatteredMediumParticles[i].Pos.Z() << sep 
         << scatteredMediumParticles[i].Mom.E() << sep << scatteredMediumParticles[i].Mom.Px() << sep 
         << scatteredMediumParticles[i].Mom.Py() << sep << scatteredMediumParticles[i].Mom.Pz() << sep 
         << scatteredMediumParticles[i].md2g << sep << scatteredMediumParticles[i].md2q << sep 
         << scatteredMediumParticles[i].N_EVENT_pp << endl;
  }
  file.close();
}

void analysis::saveNumberOfMediumParticles( const int step)
{
  for( int i = 0; i < particles_atTimeNow.size(); i++ )
  {
    if(fabs(particles_atTimeNow[i].Pos.Rapidity() )<0.1 )
    {
      FLAVOR_TYPE genFlavor = ParticleOffline::mapToGenericFlavorType( static_cast<FLAVOR_TYPE>( particles_atTimeNow[i].FLAVOR ) );

      switch(genFlavor)
      {
        case light_quark:
          NumberOfQuarks[step]++;
          break;
        case anti_light_quark:
          NumberOfAntiquarks[step]++;
          break;
        case gluon:
          NumberOfGluons[step]++;
          break;
      }  
    }
  }
}

void analysis::printNumberOfMediumParticles()
{  
   //creates filename
  string filename = filename_prefix + "_" + "_NumberOfMediumParticles.dat";
  fstream file( filename.c_str(), ios::out | ios::trunc ); 
  //---- print header if file is empty ----
  time_t end;
  time( &end );

  file.seekp( 0, ios::end );
  long size = file.tellp();
  file.seekp( 0, ios::beg );
  if( size == 0 )
    printHeader( file, numbOfPartcles, end );
  //---------------------------------------  
  file << "#All number divided by Ntest = " << theConfig->getTestparticles() << endl;
  file << "#Quarks\t#nAntiquarks\t#nGluons" << endl;
  for (int i=0;i<nTimeSteps-1;i++)
  {
    if(NumberOfQuarks[i] >0 || NumberOfAntiquarks[i]>0 || NumberOfGluons[i]>0)
    {
      file << tstep[i] << "\t" << NumberOfQuarks[i]/theConfig->getTestparticles();
      file << "\t" << NumberOfAntiquarks[i]/theConfig->getTestparticles();
      file << "\t" << NumberOfGluons[i]/theConfig->getTestparticles();
      file << endl;
    }    
  }
  
  file.close();  
}


void analysis::mediumParticlesOutput( const int step )
{
  string name;
  stringstream ss;

  if( step == 0 )
    name = "initial";
  else
    if( step == nTimeSteps )
      name = "final";
    else
    {
      ss << step;
      name = "step" + ss.str();
    }

  //creates filename, for example: "./output/run32_step1.f1", "./output/run32_initial.f1" or the likes
  string filename = filename_prefix + "_" + name + "_allMediumParticles.f1";
  fstream file( filename.c_str(), ios::out | ios::trunc );

  //---- print header if file is empty ----
  time_t end;
  time( &end );

  file.seekp( 0, ios::end );
  long size = file.tellp();
  file.seekp( 0, ios::beg );
  if( size == 0 )
    printHeader( file, all, end );
  //---------------------------------------

  for( int i = 0; i < particles_atTimeNow.size(); i++ )
  {
    file << i << sep << particles_atTimeNow[i].unique_id << sep << particles_atTimeNow[i].cell_id << sep 
         << particles_atTimeNow[i].FLAVOR << sep 
         << particles_atTimeNow[i].Pos.T() << sep << particles_atTimeNow[i].Pos.X() << sep
         << particles_atTimeNow[i].Pos.Y() << sep << particles_atTimeNow[i].Pos.Z() << sep 
         << particles_atTimeNow[i].Mom.E() << sep << particles_atTimeNow[i].Mom.Px() << sep 
         << particles_atTimeNow[i].Mom.Py() << sep << particles_atTimeNow[i].Mom.Pz() << sep 
         << particles_atTimeNow[i].md2g << sep << particles_atTimeNow[i].md2q << sep 
         << particles_atTimeNow[i].N_EVENT_pp << endl;
  }
  file.close();
}


void analysis::jetTrackerOutput()
{
  //---- get time and calculate real time needed by simulation ----
  time_t end;
  time( &end );
  int secs = ( int )difftime( end, start );
  int hours = secs / 3600, minutes = ( secs % 3600 ) / 60, seconds = secs - 3600 * hours - 60 * minutes;
  cout << hours << "h" << minutes << "m" << seconds << "s" <<  endl;
  //---------------------------------------------------------------
  
  
  //---- write output file for pt-spectra -------------------------
  string filename;
  if( studyJets )
    filename = filename_prefix + ".f4";
  else
    filename = "/dev/null";
  fstream file( filename.c_str(), ios::out | ios::trunc );
  
  //---- print header if file is empty ----
  file.seekp( 0, ios::end );
  long size = file.tellp();
  file.seekp( 0, ios::beg );
  if ( size == 0 )
    printHeader( file, jets, end );
  //---------------------------------------
  
  for ( unsigned int jet = 0; jet < jetTracker.size(); jet++ )
  {
    for ( unsigned int event = 0; event < jetTracker[jet].size(); event++ )
    {
      file << jetTracker[jet][event].jet_ID_in << sep << jetTracker[jet][event].jet_ID_out << sep
      << jetTracker[jet][event].coll_type << sep;
      for ( int i = 0; i < 4; i++ )
      {
        file << jetTracker[jet][event].P_proj_in(i) << sep;
      }
      for ( int i = 0; i < 4; i++ )
      {
        file << jetTracker[jet][event].P_proj_out(i) << sep;
      }
      for ( int i = 0; i < 4; i++ )
      {
        file << jetTracker[jet][event].R_proj(i) << sep;
      }
      file << jetTracker[jet][event].xSection << sep << jetTracker[jet][event].lambda 
      << sep << jetTracker[jet][event].flavor_in << sep << jetTracker[jet][event].flavor_out << endl;
    }
    file << endl << endl;
  }
  
}


analysisRingStructure::analysisRingStructure( const int _nRings, const double _centralRadius, const double _deltaR ) : numberOfRings( _nRings ),
    centralRingRadius( _centralRadius ), deltaR( _deltaR )
{
  rings.resize( _nRings );

  rings[0].relocate( 0, _centralRadius );

  totalRadius = _centralRadius + ( _nRings - 1 ) * deltaR;
  for ( unsigned int i = 1; i < rings.size(); i++ )
  {
    rings[i].relocate( rings[i-1].maxRadius, rings[i-1].maxRadius + deltaR );
  }
}



void analysisRingStructure::resize( const int _nRings, const double _centralRadius, const double _deltaR )
{
  numberOfRings = _nRings;
  centralRingRadius = _centralRadius;
  deltaR = _deltaR;
  
  rings.clear();
  rings.resize( _nRings );

  rings[0].relocate( 0, _centralRadius );

  totalRadius = _centralRadius + ( _nRings - 1 ) * deltaR;
  for ( unsigned int i = 1; i < rings.size(); i++ )
  {
    rings[i].relocate( rings[i-1].maxRadius, rings[i-1].maxRadius + deltaR );
  }
}



int analysisRingStructure::getIndex( const double _xt ) const
{
  unsigned int index = getIndexPure( _xt );
  if ( index >= rings.size() )
  {
    return ( static_cast<int>( rings.size() ) - 1 );
  }
  else
  {
    return index;
  }
}



int analysisRingStructure::getIndexPure( const double _xt ) const
{
  if ( _xt < 0 )
  {
    std::string errMsg = "transverse position xt < 0";
    throw eAnalysis_error( errMsg );
  }
  
  if ( _xt < centralRingRadius )
  {
    return 0;
  }
  else
  {
//     if ( _xt > totalRadius )
//     {
//       std::string errMsg = "transverse position xt > R (R = total radius of ring structure)";
//       throw eAnalysis_error( errMsg );
//     }
    
    int index = static_cast<int>( ( _xt - centralRingRadius ) / deltaR ) + 1;  // +1 since index 0 is for the central ring
    return index;
  }
}





int analysisRingStructure::getIndex( const ParticleOffline& _particle ) const
{
  return getIndex( _particle.Pos.Perp() );
}



analysisRingContainer& analysisRingStructure::operator[]( const int index )
{
  if ( index < 0 || index >= numberOfRings )
  {
    std::string errMsg = "index out of range in analysisRingStructure";
    throw eAnalysis_error( errMsg );
  }
  
  return rings[ index ];
}



analysisRingStructure& analysisRingStructure::operator+=( analysisRingStructure& rhs )
{
  if ( !( rhs.size() == ( *this ).size() ) )
  {
    std::string errMsg = "analysisRingStructure::operator+= only works for equally sized structures";
    throw eAnalysis_error( errMsg );
  }

  for ( int i = 0; i < numberOfRings; i++ )
  {
    rings[i] += rhs[i];
  }

  return ( *this );
}



analysisRingContainer& analysisRingContainer::operator+=( const analysisRingContainer & rhs )
{
  lambdaGluon += rhs.lambdaGluon;
  lambdaQuark += rhs.lambdaQuark;
  collectedGluon += rhs.collectedGluon;
  collectedQuark += rhs.collectedQuark;
  
  return ( *this );
}



// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;
