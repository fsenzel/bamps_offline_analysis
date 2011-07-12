#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "particle.h"
#include "configuration.h"
#include "parameter.h"

using namespace std;
using namespace ns_casc;

int MaxN, number, numberAdded, IX, IY, IZ;
// double as;

// int nmd2,nlam,nI23;
// double *md2a,*lam,*I23,*I23t;
int cellcut;
double deta;

initParam iparam;

void config::setting()
{
  ifstream readinpara;

  string filename;

  //-------------------------//
  // variables from cascade

  double dummy;

  filename = pathdirCascadeData + "/" + getName() + "_parameters.dat";
  readinpara.open( filename.c_str() );
  if ( !readinpara.good() )
  {
    string errMsg = "Error in reading " + filename;
    throw eConfig_error( errMsg );
  }
  else
  {
    readinpara >> originalSeed >> sqrtS >> impact >> A >> Aatomic >> B >> Batomic >> dummy 
      >> testparticles >> N_init >> timefirst >> timeshift >> freezeOutEnergyDensity
      >> ringNumber  >> centralRingRadius >> deltaR >> dx >> dy >> transLen >> IX >> IY >> IZ >> deta;
//     readinpara >> originalSeed >> sqrtS >> impact >> A >> Aatomic >> B >> Batomic >> P0 
//       >> testparticles >> N_init >> timefirst >> timeshift >> freezeOutEnergyDensity
//       >> ringNumber  >> centralRingRadius >> deltaR >> dx >> dy >> transLen >> IX >> IY >> IZ >> deta;
  }

  numberOfParticlesToAdd *= testparticles;
  numberAdded = numberOfParticlesToAdd;
  cout << "** add " << numberOfParticlesToAdd << " particles with pt_min = " << minimumPT << endl;
  cout << timefirst << "  " << N_init << endl;

  
  cellcut = 4;
  //-------------------------//

  number = N_init;
  
  particles.reserve( N_init * 1.8 );
  particles_init.reserve( N_init );
  
  particles.resize( N_init * 1.8 );
  particles_init.resize( N_init );

  filename = pathdirCascadeData + "/" + getName() + "_initial.f1";
  fstream readinp( filename.c_str(), ios::in );
  if ( !readinp.good() )
  {
    string errMsg = "Error in reading " + filename;
    throw eConfig_error( errMsg );
  }

//   double dummy;
  int FLAVOR;
  string line;
  stringstream inStream;
  int i = 0;
  while( !readinp.eof() )
  {
    getline( readinp, line );
    //only lines which are non empty and do not contain a # are taken into account
    if (line.find("#",0) == string::npos && line.length() != 0)    
    {
      inStream.clear();
      inStream.str(line);
      inStream >> dummy >> particles_init[i].unique_id >> dummy >> FLAVOR 
      >> particles_init[i].T >> particles_init[i].X >> particles_init[i].Y >> particles_init[i].Z
      >> particles_init[i].E >> particles_init[i].PX >> particles_init[i].PY >> particles_init[i].PZ
      >> particles_init[i].md2g >> particles_init[i].md2q;
      
      inStream.str("");
      particles_init[i].FLAVOR = static_cast<FLAVOR_TYPE>( FLAVOR );
      ++i;
    }
  }
  readinp.close();

  if ( particles_init.size() != N_init )
  {
    string errMsg = "Wrong particle number";
    throw eConfig_error( errMsg );
  }

  for ( int i = 0; i < particles_init.size(); i++ )
  {
    // compute energy from momenta again, to get higher precision than the numbers from file (otherwise energy conservation is not fullfilled to the last decimal place
    particles_init[i].E = sqrt( pow( particles_init[i].PX, 2.0 ) + pow( particles_init[i].PY, 2.0 ) + pow( particles_init[i].PZ, 2.0 ) + pow( particles_init[i].m, 2.0 ) );
    
    particles[i].init = true;
    particles[i].free = false;

    particles[i] = particles_init[i];
    
    if ( particles[i].T < timefirst )
    {
      particles[i].edge = true;
      particles[i].init = false;
    }

    particles[i].PXold = particles_init[i].PX;
    particles[i].PYold = particles_init[i].PY;
    particles[i].PZold = particles_init[i].PZ;
    particles[i].Eold  = particles_init[i].E;
  }

  particles_atTimeNow = particles;
}


// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on; 
