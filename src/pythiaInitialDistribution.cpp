//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL: svn+ssh://fochler@th.physik.uni-frankfurt.de/home/fochler/svnrepos/cascade_full/main/trunk/src/minijets.cpp $
//$LastChangedDate: 2010-10-22 14:29:14 +0200 (Fri, 22 Oct 2010) $
//$LastChangedRevision: 217 $
//$LastChangedBy: fochler $
//---------------------------------------------
//---------------------------------------------


#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <list>
#include <algorithm> // for std::count

#include <sstream>
#include <list>

#include "configuration.h"
#include "pythiaInitialDistribution.h"
#include "psf.h"
#include "random.h"
#include "particle.h"
#include "FPT_compare.h"
#include "integrand_time.h"

using std::cout;
using std::endl;
using namespace ns_casc;


pythiaInitialDistribution::pythiaInitialDistribution( const config& _config, WoodSaxon& _WoodSaxonParameter, STORED_TABLE_USAGE _storedTableUsage ) :
    filename_samplingData_collisionTimes( "data/samplT.dat" ),
    gamma( 1 ),     // set in computeWoodSaxonParameter()
    velocity( 0 ),  // set in computeWoodSaxonParameter()
    RA0( 0 ),       // set in computeWoodSaxonParameter()
    RA( 0 ),        // set in computeWoodSaxonParameter()
    dA( 0 ),        // set in computeWoodSaxonParameter()
    n0A( 0 ),       // set in computeWoodSaxonParameter()
    nEntries_collisionTimes( 0 ),
    numberOfParticlesToGenerate( 0 ),
    numberOfTestparticles( _config.getTestparticles() ),
    A( _config.getA() ),
    Aatomic( _config.getAatomic() ),
    B( _config.getB() ),
    Batomic( _config.getBatomic() )
{
  impactParameter = _config.getImpactParameter();
  sqrtS_perNN = _config.getSqrtS();
  filename_pythiaParticleFile = _config.getPythiaParticleFile();
  computeWoodSaxonParameters( _config, _WoodSaxonParameter );
  WoodSaxonParameter = _WoodSaxonParameter;


  setDataFilesProperties( _storedTableUsage );


  cout << "PYTHIA particle data file: " << filename_pythiaParticleFile << endl;

  std::ifstream countPythiaParticles( filename_pythiaParticleFile.c_str() );
  if ( countPythiaParticles.good() )
  {
    numberOfParticlesToGenerate = std::count( std::istreambuf_iterator<char>( countPythiaParticles ), std::istreambuf_iterator<char>(), '\n' ); // counts number of lines in PYTHIA particle data file
  }
  else
  {
    cout << "Error at opening PYTHIA particle data file. Glauber model will be used to compute particle number." << endl;
  }

  countPythiaParticles.close();
}



void pythiaInitialDistribution::setDataFilesProperties( STORED_TABLE_USAGE _storedTableUsage )
{
  if ( _storedTableUsage == useStoredTables )
  {
    cout << "======= Using stored data sets for sampling of initial state =======" << endl;


    if ( sqrtS_perNN == 5500.0 ) // LHC, 5.5 TeV
    {
      // for p0=3.3 GeV
      nEntries_collisionTimes = 127;

      if ( impactParameter == 0.0 )
      {
        filename_samplingData_collisionTimes = "data/samplT_b0.dat";//b=0 fm
      }
      else
      {
        cout << "======= Generating data sets for sampling of initial state =======" << endl;
        generateSamplingDataSets();
        cout << "======= data sets generated =======" << endl;

        filename_samplingData_collisionTimes = "data/samplT.dat";
      }

    }
    else if ( sqrtS_perNN == 2760.0 ) // LHC, 2.76 TeV
    {
      // for p0=2.9 GeV
      nEntries_collisionTimes = 129;

      if ( impactParameter == 0.0 )
      {
        filename_samplingData_collisionTimes = "data/samplT_b0.dat";//b=0 fm
      }
      else
      {
        cout << "======= Generating data sets for sampling of initial state =======" << endl;
        generateSamplingDataSets();
        cout << "======= data sets generated =======" << endl;

        filename_samplingData_collisionTimes = "data/samplT.dat";
      }
    }
    else if ( sqrtS_perNN == 200.0 ) // RHIC 200 MeV
    {
      nEntries_collisionTimes = 143;

      if ( impactParameter == 0.0 )
        filename_samplingData_collisionTimes = "data/samplT_b0.dat";//b=0 fm
      else if ( impactParameter == 3.3 )
        filename_samplingData_collisionTimes = "data/samplT_b3_3.dat";//b=3.3 fm
      else if ( impactParameter == 4.6 )
        filename_samplingData_collisionTimes = "data/samplT_b4_6.dat";//b=4.6 fm
      else if ( impactParameter == 5.8 )
        filename_samplingData_collisionTimes = "data/samplT_b5_8.dat";//b=5.8 fm
      else if ( impactParameter == 8.6 )
        filename_samplingData_collisionTimes = "data/samplT_b8_6.dat";//b=8.6 fm
      else if ( impactParameter == 8.2 )
        filename_samplingData_collisionTimes = "data/samplT_b8_2.dat";//b=8.2 fm
      else if ( impactParameter == 10.3 )
        filename_samplingData_collisionTimes = "data/samplT_b10_3.dat";//b=10.3 fm
      else
      {
        cout << "======= Generating data sets for sampling of initial state =======" << endl;
        generateSamplingDataSets();
        cout << "======= data sets generated =======" << endl;

        filename_samplingData_collisionTimes = "data/samplT.dat";
      }
    }
    else
    {
      cout << "======= Generating data sets for sampling of initial state =======" << endl;
      generateSamplingDataSets();
      cout << "======= data sets generated =======" << endl;

      filename_samplingData_collisionTimes = "data/samplT.dat";
    }
  }
  else
  {
    cout << "======= Generating data sets for sampling of initial state (mini jets) =======" << endl;
    generateSamplingDataSets();
    cout << "======= data sets generated =======" << endl;

    filename_samplingData_collisionTimes = "data/samplT.dat";
  }

  if ( !samplingDataSetsExist() )
  {
    cout << "======= Generating data sets for sampling of initial state =======" << endl;
    generateSamplingDataSets();
    cout << "======= data sets generated =======" << endl;

    filename_samplingData_collisionTimes = "data/samplT.dat";
  }

}



void pythiaInitialDistribution::generateSamplingDataSets()
{
  double Tab = generateTimeDistribution( nEntries_collisionTimes );
  cout << "++++  Tab = " << Tab << "1/mb" << endl;
}



void pythiaInitialDistribution::populateParticleVector( std::vector< Particle >& _particles )
{
  /**
  * Reserve memory for the Particle vector. Due to possible particle creation this size should be chosen somewhat larger
  * than the initial number of particles. Doing this, push_back operations to add new particles won't lead to internal
  * memory re-allocation in the vector, which could possibly be time consuming.
  * However push_back operations are always possible, even when the reserved size is used up (as long as there is physical
  * memory available of course). And modern day compilers probably optimize to the extent where it doesn't really matter -
  * but still, it's good practice.
  */
  _particles.reserve( static_cast<int>( numberOfParticlesToGenerate * 1.2 ) );
  /**
  * Now the particle vector is re-sized to hold the initial number of particles (this is NOT done when reserving memory!).
  * The particles are initialised with the standard constructor, actual attributes such as momentum etc. MUST be set later!
  */
  _particles.resize( numberOfParticlesToGenerate );

  sampleMomenta( _particles );
  samplePositions( _particles );
}



void pythiaInitialDistribution::populateParticleVector(std::vector< Particle >& _particles, const int _numberOfParticlesToGenerate, const double _minimumPT)
{
  minimumPT = _minimumPT;
  numberOfParticlesToGenerate = _numberOfParticlesToGenerate;
  _particles.reserve( numberOfParticlesToGenerate * 1.5 );
  _particles.resize( numberOfParticlesToGenerate );
  sampleMomenta( _particles );
  samplePositions( _particles );
}





// same as for miniJets
void pythiaInitialDistribution::computeWoodSaxonParameters( const config& _config, WoodSaxon& _WoodSaxonParameter )
{
  const double mproton = 0.938272;
  const double neps = 1.0e-3;
  const double A = _config.getA();
  const double b = _config.getImpactParameter();

  gamma = _config.getSqrtS() / ( 2 * mproton );
  velocity = sqrt( 1.0 - 1 / pow( gamma, 2 ) );

  RA = 1.12 * pow( A, 1.0 / 3.0 ) - 0.86 * pow( A, -1.0 / 3.0 );
  dA = 0.54;
  n0A = 3 * A / ( 4 * M_PI * pow( RA, 3.0 ) * ( 1.0 + pow(( M_PI * dA / RA ), 2 ) ) );

  RA0 = dA * log( gamma * n0A / neps - 1.0 ) + RA;  //radius where density has dropped to neps = 10^-3

  if ( b > ( 2 * RA0 ) )
  {
    std::string errMsg = "Impact parameter b too large. b > 2 R_A0";
    throw ePythia_error( errMsg );
  }

  _WoodSaxonParameter.dA = dA;
  _WoodSaxonParameter.RA = RA;
  _WoodSaxonParameter.RA0 = RA0;
  _WoodSaxonParameter.gamma = gamma;
  _WoodSaxonParameter.n0A = n0A;
  _WoodSaxonParameter.velocity = velocity;


  cout << "long. extension = " << 2 * RA0 / gamma << " fm" << endl;
  cout << "trans. radius = " << RA << " fm" << endl;
  cout << "impact parameter = " << b << " fm" << endl;
  cout << "overlap time = " << sqrt( pow( 2 * RA0, 2 ) - pow( b, 2 ) ) / gamma / velocity << " fm/c" << endl;
}



void pythiaInitialDistribution::sampleMomenta( std::vector< Particle >& _particles )
{

  std::ifstream readPythiaParticles( filename_pythiaParticleFile.c_str() );
  cout << "Read particle momentum from PYTHIA data file." << endl;
  for ( int i = 0; i < numberOfParticlesToGenerate; i++ )
  {
    int flavTemp;
    // structure of file
    // number of pythia event  is hard?   flavour  energy  px  py  pz  m
    readPythiaParticles >> _particles[i].N_EVENT >> _particles[i].HARD >> flavTemp >> _particles[i].E >> _particles[i].PX >> _particles[i].PY >> _particles[i].PZ >> _particles[i].m;
    _particles[i].FLAVOR = static_cast<FLAVOR_TYPE>( flavTemp );
    
    if ( flavTemp <= 6 )
    {
      _particles[i].m = 0;
    }
    
    // to avoid rounding errors compute energy from momenta and mass E^2=p^2+m^2
    _particles[i].E = sqrt( pow( _particles[i].PX, 2.0 ) + pow( _particles[i].PY, 2.0 ) + pow( _particles[i].PZ, 2.0 ) + pow( _particles[i].m, 2.0 ) );
  }
}




void pythiaInitialDistribution::samplePositions( std::vector< Particle >& _particles )
{
  cout << "Sample particle positions for particles from PYTHIA." << endl;
  // particles are already read from file in momentum()
  int nmb_of_events = _particles.back().N_EVENT;
  int event_tmp = 0;
  double x_tmp, y_tmp, z_tmp, t_tmp;

  bool soft_event[nmb_of_events + 1]; // array to store which events are soft and which not (true if event is soft)
  for ( int i = 1; i <= nmb_of_events; i++ )
  {
    soft_event[i] = false;
  }
  int partcl_soft_event[nmb_of_events + 1]; // array to store one particle number of each event to obtain the position of the event afterwards for soft particles

  // sample positions for binary collisions/events where hard partons are produced. Sample also if there are also soft particles at this event (reason: scaling behavior of soft and hard partons differ)
  for ( int j = 0; j < numberOfParticlesToGenerate; j++ )
  {
    if ( _particles[j].HARD )
    {
      if ( _particles[j].N_EVENT != event_tmp ) // first hard particle of an event
      {
        sample_T_one_partcl( _particles, j );
        sample_XYZ_one_partcl( _particles, j, soft_event[_particles[j].N_EVENT] );

        x_tmp = _particles[j].X;
        y_tmp = _particles[j].Y;
        z_tmp = _particles[j].Z;
        t_tmp = _particles[j].T;

        partcl_soft_event[_particles[j].N_EVENT] = j;

        event_tmp = _particles[j].N_EVENT;
      }
      else // all other hard particles of same event get same positions
      {
        _particles[j].X = x_tmp;
        _particles[j].Y = y_tmp;
        _particles[j].Z = z_tmp;
        _particles[j].T = t_tmp;
      }
    }
  }

  // start insert for non vanishing b=Bimp
  if ( impactParameter != 0.0 ) // soft particles are deleted for b!=0, because treatment of these is written for b=0 and would cause an error
  {
    std::vector<Particle>::iterator iIt;
    //!! delete particles
    for ( iIt = _particles.begin(); iIt != _particles.end(); )
    {
      if ( ( *iIt ).HARD == false )
      {
        iIt = _particles.erase( iIt );
      }
      else
      {
        ++iIt;
      }
    }
    cout << "Soft particles have been deleted since b!=0 and this would cause an error!!!" << endl;
  }
  // end of insert
  else
  {
    int sum = 0;
    int event = 1;
    // get positions for soft partons from events which have besides hard also soft partons
    int soft_break = 0;
    for ( int j = 0; j < numberOfParticlesToGenerate; j++ )
    {
      if ( !_particles[j].HARD )
      {
        sum++;
      }
    }
    cout <<  "Soft particles:" << sum << endl;

    sum = 0;
    for ( int i = 1; i <= nmb_of_events; i++ )
    {
      if ( soft_event[i] )
      {
        sum++;
      }
    }
    if ( sum == 0 )
    {
      string errMsg = "ERROR: No soft events! This will result in an infinite loop!";
      throw ePythia_error( errMsg );
    }
    cout << "Soft events: " << sum << endl;

    for ( int s = 0; s <  numberOfParticlesToGenerate; s++ )
    {
      if ( !_particles[s].HARD )
      {
        while ( !soft_event[event] )
        {
          if ( event < nmb_of_events )
          {
            event++;
          }
          else
          {
            soft_break = s;   // there are more soft particles than events for it
            event = 1;
          }
        }

        // event is soft
        _particles[s].X = _particles[partcl_soft_event[event]].X;
        _particles[s].Y = _particles[partcl_soft_event[event]].Y;
        _particles[s].Z = _particles[partcl_soft_event[event]].Z;
        _particles[s].T = _particles[partcl_soft_event[event]].T;
        event++;
      }
    }
  }
}



//---------------------------------//
// sampling of the collision times //
//---------------------------------//
// sample collision times for only one particle with id=number
void pythiaInitialDistribution::sample_T_one_partcl( std::vector< Particle >& _particles, const int number ) const
{
  const int ord = 4;
  int nn;
  double *xx, *yy;
  double sample_y, t, dt;
  double xa[ord+1], ya[ord+1];

  std::fstream read( filename_samplingData_collisionTimes.c_str(), std::ios::in );

  xx = new double[nEntries_collisionTimes];
  yy = new double[nEntries_collisionTimes];
  for ( int i = 0;i < nEntries_collisionTimes;i++ )
  {
    read >> xx[i] >> yy[i];
  }

  sample_y = double( ran2() );

  nn = 0;
  while ( yy[++nn] < sample_y );
  nn = nn - ord / 2;
  if ( nn < 0 )
  {
    nn = 0;
  }
  else if (( nn + ord - 1 > nEntries_collisionTimes - 1 ) || ( yy[nn+ord-1] == 1.0 ) )
  {
    while ( yy[++nn] < 1.0 );
    nn = nn - ord + 1;
  }

  for ( int i = 1;i <= ord;i++ )
  {
    ya[i] = yy[nn+i-1];
    xa[i] = xx[nn+i-1];
  }

  polint( ya, xa, ord, sample_y, &t, &dt );
  _particles[number].T = t;

  read.close();

  delete[] xx;
  delete[] yy;
}



//---------------------------------------------------------------//
// sampling of the position of one parton at the given time//
//---------------------------------------------------------------//
// sample position for only one particle with id=number
void pythiaInitialDistribution::sample_XYZ_one_partcl( std::vector< Particle >& _particles, const int number, bool& soft ) const
{
  double X, Y, Z, T;
  double vt, gv, gvt, tc1, tc2, tc3, xx, bA, bB, zA, zB, fds, fd;
  double zmax, zmin, xmax, xmin, ymax, ymin, max;
  double c1, c2;

  double p_soft, densityA_max, L_z;
  double sigma = 40.0; //  p+p cross section in mb
  sigma = sigma * 0.1; //in fm^2
  double m = 2.0; // L_z is z component of radius where densityA is only 1/m-th of maximum value
//   double m = 2.718; // L_z is z component of radius where densityA is only 1/m-th of maximum value

  T = _particles[number].T;
  vt = velocity * fabs( T );
  gv = gamma * velocity;
  gvt = gv * fabs( T );

  c1 = sqrt(( RA0 + impactParameter ) * fabs( RA0 - impactParameter ) );
  tc1 = ( RA0 - c1 ) / 2.0 / gv;
  tc2 = sqrt( impactParameter * ( 2.0 * RA0 - impactParameter ) ) / 2.0 / gv;
  tc3 = ( RA0 + c1 ) / 2.0 / gv;

  xx = impactParameter / 2.0 + gvt * sqrt( RA0 * RA0 / ( impactParameter * impactParameter / 4.0 + gvt * gvt ) - 1.0 );

  // z-axis in overlap region
  if (( RA0 > impactParameter ) && ( fabs( T ) >= tc1 ) && ( fabs( T ) <= tc3 ) ) zmin = vt - RA0 / gamma;
  else zmin = vt - sqrt(( RA0 + impactParameter - xx ) * ( RA0 - impactParameter + xx ) ) / gamma;
  zmax = -zmin;

  // x-axis in overlap region
  if ( fabs( T ) <= tc2 ) xmax = RA0;
  else xmax = xx;
  xmin = impactParameter - xmax;

  // y-axis in overlap region
  ymax = sqrt( RA0 * RA0 - impactParameter * impactParameter / 4.0 - gvt * gvt );
  ymin = -ymax;

  // the maximum of nA(s,z-vt)*nB(s-b,z+vt) in the overlap region
  max = densityA( impactParameter / 2.0, -velocity * T ) * densityA( impactParameter / 2.0, velocity * T );

  // sampling of position
  do
  {
    X = ( xmax - xmin ) * ran2() + xmin;
    Y = ( ymax - ymin ) * ran2() + ymin;
    Z = ( zmax - zmin ) * ran2() + zmin;

    bA = sqrt( X * X + Y * Y );
    bB = sqrt(( X - impactParameter ) * ( X - impactParameter ) + Y * Y );

    if (( bA > RA0 ) || ( bB > RA0 ) ) fds = -1.0;
    else
    {
      zA = sqrt(( RA0 + bA ) * ( RA0 - bA ) ) / gamma;
      zB = sqrt(( RA0 + bB ) * ( RA0 - bB ) ) / gamma;
      c1 = fabs( Z - velocity * T );
      c2 = fabs( Z + velocity * T );
      if (( c1 > zA ) || ( c2 > zB ) ) fds = -1.0;
      else fds = densityA( bA, c1 ) * densityA( bB, c2 );
    }

    if ( fds > max ) cout << "problem in sample_posit_one_partcl" << endl;

    fd = max * ran2();
  }
  while ( fds < fd );

  _particles[number].X = X;
  _particles[number].Y = Y;
  _particles[number].Z = Z;

  // sample if there are also soft partons at this position
  densityA_max = densityA( impactParameter / 2.0, velocity * T );
//   L_z = 1.0/gama*sqrt(pow(dA*log(gama*m-1.0)+RA,2.0)-pow(X,2.0)-pow(Y,2.0)); // only valid for central collision b(=impactParameter)=0 -> I think this formula is wrong> next line should be right!
  L_z = 2.0 / gamma * sqrt( pow( RA, 2.0 ) - pow( X, 2.0 ) - pow( Y, 2.0 ) ); // only valid for central collision b(=impactParameter)=0
  p_soft = 1.0 / sigma / densityA_max / L_z;

//   if(p_soft > 1.0)
//   {
//     cout << "error, p_soft>1 in init_pos() (but doen't matter; so it's clearly a soft event), p_soft=" << p_soft << endl;
//   }

//   if(impactParameter != 0.0) // b != 0 would cause an error since consideration above is only for b=0
//     soft = false;
//   else
  if ( ran2() < p_soft )
  {
    soft = true;
  }
  else
  {
    soft = false;
  }
}






bool pythiaInitialDistribution::samplingDataSetsExist() const
{
  bool exists = true;
  std::ifstream testFile;

  //--------------------------------  check file 1 --------------------------------
  testFile.clear();
  testFile.open( filename_samplingData_collisionTimes.c_str(), std::ios::in );
  if ( testFile.good() )  //check if quark file actually exists
  {
    testFile.close();
  }
  else
  {
    exists = false;
    cout << "Could not open: " << filename_samplingData_collisionTimes << endl;
  }
  //--------------------------------  check file 1 --------------------------------

  return exists;
}




// same as in miniJets
//-----------------------------------------------//
// time distribution, calculation runs only once.//
//-----------------------------------------------//
double pythiaInitialDistribution::generateTimeDistribution( int& _count )
{
  int nn, nc, ncut;
  double dt, tmp;
  double tgral, sd, chi2a;
  double *dist, *tt, *distin;
  double Tab;

  integrand_time ftime;
  ftime.setB( impactParameter );
  ftime.setWoodSaxonParameter( WoodSaxonParameter );

  std::fstream file_samplT( filename_samplingData_collisionTimes.c_str(), std::ios::out | std::ios::trunc );
  std::fstream file_tdist( "data/tdist.dat", std::ios::out | std::ios::trunc );
  file_samplT.precision( 10 );
  file_tdist.precision( 10 );
  string sep = "\t";

  double tmax = sqrt( 4 * pow( RA0, 2 ) - pow( impactParameter, 2 ) ) / ( 2 * velocity * gamma );

  nn = 100;
  dt = tmax / nn;
  nc = 2 * nn + 1;
  dist = new double[nc];
  tt = new double[nc];
  distin = new double[nc];

  for ( int i = 0; i < nc; i++ )
  {
    dist[i] = 0.0;
  }
  tt[nn] = 0.0;
  for ( int j = nn + 1; j < nc; j++ )
  {
    tt[j] = tt[j-1] + dt;
    tt[nc-1-j] = tt[nc-j] - dt;
  }

  ncut = 0;
  for ( int i = nn;i < nc;i++ )
  {
    ftime.setTime( tt[i] );
    vegas( 3, ftime, &tgral, &sd, &chi2a );
    dist[i] = dist[nc-1-i] = double( tgral );
    if ( tgral < 1.0e-7 )
    {
      i = nc;
    }
  }

  distin[0] = 0.0;
  for ( int j = 1;j < nc;j++ )
  {
    distin[j] = distin[j-1] + 0.5 * ( dist[j-1] + dist[j] ) * dt;
  }

  Tab = distin[nc-1] * 2.0 * velocity / 10.;// 1/10 due to the rescale 1/fm²->1/10mb

  ncut = nc;
  do
  {
    ncut -= 1;
    tmp = ( distin[ncut] - distin[ncut-1] ) / distin[nc-1];
  }
  while ( tmp < 1.0e-6 );
  _count = nc - 2 * ( nc - 1 - ncut );

  for ( int i = nc - 1 - ncut;i < ncut + 1;i++ )
  {
    file_samplT << tt[i] << sep << distin[i] / distin[ncut] << endl;
    file_tdist << tt[i] << sep << dist[i] / distin[ncut] << endl;
  }

  delete[] dist;
  delete[] tt;
  delete[] distin;

  file_samplT.close();
  file_tdist.close();

  cout << "++++  generated data set: " << filename_samplingData_collisionTimes << "  with " << _count << " entries." << endl;
  return Tab;
}





// // same as in miniJets
// void integrand_time::operator()( const int *ndim, const double xx[], const int *ncomp, double ff[] ) const
// {
//   double wgt;
//   ff[0] = this->operator()( xx, wgt );
// }
// 
// 
// // same as in miniJets
// double integrand_time::operator()( const double x[], double wgt ) const  //wgt is a dummy variable needed by VEGAS routine, x[] is {x,z,y}
// {
//   double gama, velocity, RA0;
//   double vt, gvt, c1, c2;
//   double V, xmax, xmin, xx;
//   double zAl, zAr, zBl, zBr, zmax, zmin, zz;
//   double yA2, yB2, ymax, yy;
//   double b1, b2;
// 
//   gama = woodSaxonParameter.gamma;
//   velocity = woodSaxonParameter.velocity;
//   RA0 = woodSaxonParameter.RA0;
// 
//   vt = velocity * time;
//   gvt = gama * vt;
// 
//   c1 = 4.0 * pow( RA0, 2.0 ) / ( pow( bImp, 2.0 ) + 4.0 * pow( gvt, 2.0 ) );
// 
//   // choose a x[fm] value
//   if ( 2.0*gvt > sqrt( bImp*( 2.0*RA0 - bImp ) ) )
//     xmax = bImp / 2.0 + gvt * sqrt( c1 - 1.0 );
//   else
//     xmax = RA0;
//   xmin = bImp - xmax;
// 
//   if ( fabs( xmax - xmin ) < 1.0e-8 )
//     return 0.0;
// 
//   V = xmax - xmin;
//   xx = V * double( x[1] ) + xmin;
// 
//   // choose a z[fm] value
//   c1 = sqrt(( RA0 + xx ) * ( RA0 - xx ) ) / gama;
//   c2 = sqrt(( RA0 + xx - bImp ) * ( RA0 - xx + bImp ) ) / gama;
//   zAl = vt - c1;
//   zAr = vt + c1;
//   zBl = -vt - c2;
//   zBr = -vt + c2;
//   if ( zAl > zBl )
//     zmin = zAl;
//   else
//     zmin = zBl;
//   if ( zAr < zBr )
//     zmax = zAr;
//   else
//     zmax = zBr;
// 
//   if ( fabs( zmax - zmin ) < 1.0e-8 )
//     return 0.0;
// 
//   V = V * ( zmax - zmin );
//   zz = ( zmax - zmin ) * double( x[2] ) + zmin;
// 
//   // choose a y[fm] value
//   yA2 = ( RA0 + gama * ( zz - vt ) ) * ( RA0 - gama * ( zz - vt ) ) - xx * xx;
//   yB2 = ( RA0 + gama * ( zz + vt ) ) * ( RA0 - gama * ( zz + vt ) ) - ( bImp - xx ) * ( bImp - xx );
// 
//   if (( fabs( yA2 ) < 1.0e-8 ) || ( fabs( yB2 ) < 1.0e-8 ) )
//     return 0.0;
// 
//   if ( yA2 <= yB2 )
//     ymax = sqrt( yA2 );
//   else
//     ymax = sqrt( yB2 );
// 
//   V = V * 2.0 * ymax;
//   yy = 2.0 * ymax * double( x[3] ) - ymax;
// 
//   b1 = sqrt( xx * xx + yy * yy );
//   b2 = sqrt(( xx - bImp ) * ( xx - bImp ) + yy * yy );
// 
//   return double( V*pythiaInitialDistribution::densityA( b1, zz - vt, woodSaxonParameter )*pythiaInitialDistribution::densityA( b2, zz + vt, woodSaxonParameter ) );
// }



// same as in miniJets
//--------------------//
//   interpolation    //
//--------------------//
void pythiaInitialDistribution::polint( const double xa[], const double ya[], const int n, const double x, double *y, double *dy ) const
{

  int ns = 1;
  double den, dif, dift, ho, hp, w;
  double *c, *d;

  dif = fabs( x - xa[1] );

  c = new double[n+1];
  d = new double[n+1];

  for ( int i = 1;i <= n;i++ )
  {
    if (( dift = fabs( x - xa[i] ) ) < dif )
    {
      ns = i;
      dif = dift;
    }
    c[i] = ya[i];
    d[i] = ya[i];
  }

  *y = ya[ns--];

  for ( int m = 1;m < n;m++ )
  {
    for ( int i = 1;i <= n - m;i++ )
    {
      ho = xa[i] - x;
      hp = xa[i+m] - x;
      w = c[i+1] - d[i];
      if (( den = ho - hp ) == 0.0 )
      {
        cout << "Error in routine polint" << endl;
        cout << x << endl;
        cout << xa[1] << "\t" << ya[1] << endl;
        cout << xa[2] << "\t" << ya[2] << endl;
        cout << xa[3] << "\t" << ya[3] << endl;
        cout << xa[4] << "\t" << ya[4] << endl;
      }
      den = w / den;
      d[i] = hp * den;
      c[i] = ho * den;
    }

    *y += ( *dy = ( 2 * ns < ( n - m ) ? c[ns+1] : d[ns--] ) );
  }

  delete[] c;
  delete[] d;
}



// same as in miniJets
double pythiaInitialDistribution::densityA( double b, double z ) const
{
  return gamma*n0A / ( exp(( sqrt( b*b + pow( double( gamma*z ), 2.0 ) ) - RA ) / dA ) + 1.0 );
}

// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;
