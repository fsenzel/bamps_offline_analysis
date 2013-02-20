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
#include <string>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <list>
#include <algorithm> // for std::count

#include "configuration.h"
#include "initialmodel_cgc.h"
#include "random.h"
#include "particle.h"


using std::cout;
using std::endl;


initialModel_CGC::initialModel_CGC ( const config& _config ) : initialModel()
{
  filename_cgcParticleFile = _config.getCgcParticleFile();
  cout << "CGC particle data file: " << filename_cgcParticleFile << endl;
  std::ifstream countCGCParticles ( filename_cgcParticleFile.c_str() );
  if ( countCGCParticles.good() )
  {
    numberOfParticlesToGenerate = std::count ( std::istreambuf_iterator<char> ( countCGCParticles ), std::istreambuf_iterator<char>(), '\n' ); // counts number of lines in PYTHIA particle data file
  }
  else
    cout << "Error at opening CGC particle data file." << endl;

  countCGCParticles.close();
}



void initialModel_CGC::populateParticleVector ( std::vector< Particle >& _particles )
{
  /**
  * Reserve memory for the Particle vector. Due to possible particle creation this size should be chosen somewhat larger
  * than the initial number of particles. Doing this, push_back operations to add new particles won't lead to internal
  * memory re-allocation in the vector, which could possibly be time consuming.
  * However push_back operations are always possible, even when the reserved size is used up (as long as there is physical
  * memory available of course). And modern day compilers probably optimize to the extent where it doesn't really matter -
  * but still, it's good practice.
  */
  _particles.reserve ( static_cast<int> ( numberOfParticlesToGenerate * 1.2 ) );
  /**
  * Now the particle vector is re-sized to hold the initial number of particles (this is NOT done when reserving memory!).
  * The particles are initialised with the standard constructor, actual attributes such as momentum etc. MUST be set later!
  */
  _particles.resize ( numberOfParticlesToGenerate );

  //----
  const double tau0 = 0.1; // intitial time 1/1.4 GeV^-1 = 0.14 fm (~1.4 GeV ist saturation momentum)
  //----

  std::ifstream cgcParticles ( filename_cgcParticleFile.c_str() );
  char buf[250];
  double rx, ry, y, pt, phi, eta;
  int n = 0;

  do
  {
    cgcParticles.getline ( buf, 200 );
    std::sscanf ( buf, "%lf %lf %lf %lf ", &rx, &ry, &y, &pt );
    if ( sqrt ( rx * rx + ry * ry ) <= 100.0 ) // just to check if rx, ry are from actual data
    {
      n++;

      _particles[n].T = tau0 * cosh ( y );
      _particles[n].X = rx;
      _particles[n].Y = ry;
      _particles[n].Z = tau0 * sinh ( y );

      _particles[n].PZ = pt * sinh ( y );
      phi = 2.0 * M_PI * ran2();
      _particles[n].PX = pt * cos ( phi );
      _particles[n].PY = pt * sin ( phi );
      _particles[n].E = sqrt ( _particles[n].PX * _particles[n].PX + _particles[n].PY * _particles[n].PY + _particles[n].PZ * _particles[n].PZ );

      _particles[n].m = 0.0;
      _particles[n].FLAVOR = gluon;

    }
  }
  while ( ! ( cgcParticles.eof() ) );

  if ( n != numberOfParticlesToGenerate )
    cout << "Error in initialModel_CGC::populateParticleVector, wrong number of particles" << endl;
}


