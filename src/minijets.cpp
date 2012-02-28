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
#include <iostream>
#include <fstream>
#include <list>

#include "configuration.h"
#include "minijets.h"
#include "psf.h"
#include "random.h"
#include "particle.h"
#include "FPT_compare.h"

using std::cout;
using std::endl;
using namespace ns_casc;


miniJets::miniJets( const config& _config, WoodSaxon& _WoodSaxonParameter, STORED_TABLE_USAGE _storedTableUsage ) :
    filename_samplingData_collisionTimes( "data/samplT.dat" ),
    filename_samplingData_PT( "data/samplPT.dat" ),
    filename_samplingData_PT_fine( "data/samplPT_fine.dat" ),
    filename_samplingData_maximumPT( "data/maxPT.dat" ),
    gamma( 1 ),     // set in computeWoodSaxonParameter()
    velocity( 0 ),  // set in computeWoodSaxonParameter()
    RA0( 0 ),       // set in computeWoodSaxonParameter()
    RA( 0 ),        // set in computeWoodSaxonParameter()
    dA( 0 ),        // set in computeWoodSaxonParameter()
    n0A( 0 ),       // set in computeWoodSaxonParameter()
    nEntries_collisionTimes( 0 ),
    nEntries_PT( 0 ),
    nEntries_PT_fine( 0 ),
    numberOfParticlesToGenerate( 0 ),
    numberOfTestparticles( _config.getTestparticles() ),
    A( _config.getA() ),
    Aatomic( _config.getAatomic() ),
    B( _config.getB() ),
    Batomic( _config.getBatomic() )
{
  impactParameter = _config.getImpactParameter();
  sqrtS_perNN = _config.getSqrtS();
  P0 = _config.getPtCutoff();
  computeWoodSaxonParameters( _config, _WoodSaxonParameter );
  WoodSaxonParameter = _WoodSaxonParameter;
  
  setDataFilesProperties( _config, _storedTableUsage );
}



void miniJets::setDataFilesProperties( const config& _config, STORED_TABLE_USAGE _storedTableUsage )
{
  bool correct_p0 = false;
  if (  ( FPT_COMP_E(_config.getPtCutoff(), 1.4 ) && sqrtS_perNN ==  200.0 ) ||
        ( FPT_COMP_E(_config.getPtCutoff(), 3.3 ) && sqrtS_perNN == 5500.0 ) ||
        ( FPT_COMP_E(_config.getPtCutoff(), 3.5 ) && sqrtS_perNN == 2760.0 )     )
  {
    correct_p0 = true;
  }

  if ( _storedTableUsage == useStoredTables && correct_p0 )
  {
    cout << "======= Using stored data sets for sampling of initial state =======" << endl;
    
    if(sqrtS_perNN == 5500.0 && impactParameter == 0.0) // LHC, 5.5 TeV
    {  
      // for p0=3.3 GeV
      nEntries_collisionTimes = 127;
      nEntries_PT = 108;
      nEntries_PT_fine = 29;
      
      numberOfParticlesToGenerate = 17316 * _config.getTestparticles(); 
      
      filename_samplingData_PT = "data/samplPT.dat";
      filename_samplingData_PT_fine = "data/samplPT_fine.dat";
      filename_samplingData_maximumPT = "data/maxPT.dat";
    }
    else if(sqrtS_perNN == 2760.0) // LHC, 2.76 TeV
    {  
      // for p0=3.5 GeV
      nEntries_collisionTimes = 129;
      nEntries_PT_fine = 24;
      
      filename_samplingData_PT = "data/samplPT.dat";
      filename_samplingData_PT_fine = "data/samplPT_fine.dat";
      filename_samplingData_maximumPT = "data/maxPT.dat";
      
      if ( FPT_COMP_E(_config.getImpactParameter(), 0.0 ) )
      {
        nEntries_PT = 122;
        numberOfParticlesToGenerate = 4832 * _config.getTestparticles(); 
      }
      else if ( FPT_COMP_E(_config.getImpactParameter(), 2.0 ) )
      {
        nEntries_PT = 107;
        numberOfParticlesToGenerate =  5576 * _config.getTestparticles(); 
      }
      else if ( FPT_COMP_E(_config.getImpactParameter(), 3.0 ) )
      {
        nEntries_PT = 96;
        numberOfParticlesToGenerate = 4978 * _config.getTestparticles();    
      }
      else if ( FPT_COMP_E(_config.getImpactParameter(), 4.5 ) )
      {
        nEntries_PT = 125;
        numberOfParticlesToGenerate = 3918 * _config.getTestparticles();      
      }
      else if ( FPT_COMP_E(_config.getImpactParameter(), 5.4 ) )
      {
        nEntries_PT = 112;
        numberOfParticlesToGenerate = 3264 * _config.getTestparticles();           
      }
      else if ( FPT_COMP_E(_config.getImpactParameter(), 6.1 ) )
      {
        nEntries_PT = 92;
        numberOfParticlesToGenerate = 2758 * _config.getTestparticles();       
      }
      else if ( FPT_COMP_E(_config.getImpactParameter(), 6.9 ) )
      {
        nEntries_PT = 91;
        numberOfParticlesToGenerate = 2202 * _config.getTestparticles(); 
      }
      else if ( FPT_COMP_E(_config.getImpactParameter(), 7.7 ) )
      {
        nEntries_PT = 96;
        numberOfParticlesToGenerate = 1684 * _config.getTestparticles(); 
      }
      else if ( FPT_COMP_E(_config.getImpactParameter(), 8.3 ) )
      {
        nEntries_PT = 79;
        numberOfParticlesToGenerate = 1344 * _config.getTestparticles(); 
      }
      else if ( FPT_COMP_E(_config.getImpactParameter(), 9.1 ) )
      {
        nEntries_PT = 96;
        numberOfParticlesToGenerate = 952 * _config.getTestparticles();  
      }
      else if ( FPT_COMP_E(_config.getImpactParameter(), 10.3 ) )
      {
        nEntries_PT = 101;
        numberOfParticlesToGenerate = 496 * _config.getTestparticles();        
      }
      else if ( FPT_COMP_E(_config.getImpactParameter(), 11.4 ) )
      {
        nEntries_PT = 96;
        numberOfParticlesToGenerate = 226 * _config.getTestparticles();          
      }      
      else if ( FPT_COMP_E(_config.getImpactParameter(), 12.3 ) )
      {
        nEntries_PT = 74;
        numberOfParticlesToGenerate = 102 * _config.getTestparticles();       
      }
      else if ( FPT_COMP_E(_config.getImpactParameter(), 13.2 ) )
      {
        nEntries_PT = 76;
        numberOfParticlesToGenerate = 38 * _config.getTestparticles();       
      }
      else
      {
        cout << "======= Generating data sets for sampling of initial state =======" << endl;
        generateSamplingDataSets();
        cout << "======= data sets generated =======" << endl;
        
        filename_samplingData_collisionTimes = "data/samplT.dat";
        filename_samplingData_PT = "data/samplPT.dat";
        filename_samplingData_PT_fine = "data/samplPT_fine.dat";
        filename_samplingData_maximumPT = "data/maxPT.dat";
      }
    }
    else if(sqrtS_perNN == 200.0) // RHIC 200 MeV
    {  
      nEntries_collisionTimes = 143;
      nEntries_PT_fine = 24;
      
      filename_samplingData_PT = "data/samplPT.dat";
      filename_samplingData_PT_fine = "data/samplPT_fine.dat";
      filename_samplingData_maximumPT = "data/maxPT.dat";
      
      if ( FPT_COMP_E(_config.getImpactParameter(), 0.0 ) )
      {
        nEntries_PT = 122;
        numberOfParticlesToGenerate = 4832 * _config.getTestparticles(); 
      }
      else if ( FPT_COMP_E(_config.getImpactParameter(), 2.0 ) )
      {
        nEntries_PT = 101;
        numberOfParticlesToGenerate = 4400 * _config.getTestparticles(); 
      }
      else if ( FPT_COMP_E(_config.getImpactParameter(), 2.8 ) )
      {
        nEntries_PT = 116;
        numberOfParticlesToGenerate = 4000 * _config.getTestparticles(); 
      }
      else if ( FPT_COMP_E(_config.getImpactParameter(), 3.4 ) )
      {
        nEntries_PT = 111;
        numberOfParticlesToGenerate = 3740 * _config.getTestparticles(); 
      }
      else if ( FPT_COMP_E(_config.getImpactParameter(), 4.0 ) )
      {
        nEntries_PT = 103;
        numberOfParticlesToGenerate = 3358 * _config.getTestparticles(); 
      }
      else if ( FPT_COMP_E(_config.getImpactParameter(), 4.5 ) )
      {
        nEntries_PT = 111;
        numberOfParticlesToGenerate = 3066 * _config.getTestparticles(); 
      }
      else if ( FPT_COMP_E(_config.getImpactParameter(), 5.0 ) )
      {
        nEntries_PT = 110;
        numberOfParticlesToGenerate = 2782 * _config.getTestparticles(); 
      }
      else if ( FPT_COMP_E(_config.getImpactParameter(), 5.6 ) )
      {
        nEntries_PT = 96;
        numberOfParticlesToGenerate = 2408 * _config.getTestparticles(); 
      }
      else if ( FPT_COMP_E(_config.getImpactParameter(), 6.3 ) )
      {
        nEntries_PT = 116;
        numberOfParticlesToGenerate = 2008 * _config.getTestparticles(); 
      }
      else if ( FPT_COMP_E(_config.getImpactParameter(), 7.0 ) )
      {
        nEntries_PT = 111;
        numberOfParticlesToGenerate = 1638 * _config.getTestparticles();
      }    
      else if ( FPT_COMP_E(_config.getImpactParameter(), 7.8 ) )
      {
        nEntries_PT = 99;
        numberOfParticlesToGenerate = 1230 * _config.getTestparticles(); 
      }
      else if ( FPT_COMP_E(_config.getImpactParameter(), 8.6 ) )
      {
        nEntries_PT = 106;
        numberOfParticlesToGenerate = 900 * _config.getTestparticles(); 
      }
      else if ( FPT_COMP_E(_config.getImpactParameter(), 9.6 ) )
      {
        nEntries_PT = 98;
        numberOfParticlesToGenerate = 540 * _config.getTestparticles(); 
      }
      else if ( FPT_COMP_E(_config.getImpactParameter(), 10.4 ) )
      {
        nEntries_PT = 120;
        numberOfParticlesToGenerate = 332 * _config.getTestparticles(); 
      }
      else if ( FPT_COMP_E(_config.getImpactParameter(), 11.0 ) )
      {
        nEntries_PT = 114;
        numberOfParticlesToGenerate = 216 * _config.getTestparticles(); 
      }
      else
      {
        cout << "======= Generating data sets for sampling of initial state =======" << endl;
        generateSamplingDataSets();
        cout << "======= data sets generated =======" << endl;
        
        filename_samplingData_collisionTimes = "data/samplT.dat";
        filename_samplingData_PT = "data/samplPT.dat";
        filename_samplingData_PT_fine = "data/samplPT_fine.dat";
        filename_samplingData_maximumPT = "data/maxPT.dat";
      }
    }
    else
    {
      cout << "======= Generating data sets for sampling of initial state =======" << endl;
      generateSamplingDataSets();
      cout << "======= data sets generated =======" << endl;
      
      filename_samplingData_collisionTimes = "data/samplT.dat";
      filename_samplingData_PT = "data/samplPT.dat";
      filename_samplingData_PT_fine = "data/samplPT_fine.dat";
      filename_samplingData_maximumPT = "data/maxPT.dat";
    }
  }
  else
  {
    cout << "======= Generating data sets for sampling of initial state (mini jets) =======" << endl;
    generateSamplingDataSets();
    cout << "======= data sets generated =======" << endl;
    
    filename_samplingData_collisionTimes = "data/samplT.dat";
    filename_samplingData_PT = "data/samplPT.dat";
    filename_samplingData_PT_fine = "data/samplPT_fine.dat";
    filename_samplingData_maximumPT = "data/maxPT.dat";
  }
  
  if( !samplingDataSetsExist() )
  {
    cout << "======= Generating data sets for sampling of initial state =======" << endl;
    generateSamplingDataSets();
    cout << "======= data sets generated =======" << endl;
    
    filename_samplingData_collisionTimes = "data/samplT.dat";
    filename_samplingData_PT = "data/samplPT.dat";
    filename_samplingData_PT_fine = "data/samplPT_fine.dat";
    filename_samplingData_maximumPT = "data/maxPT.dat";
  }
  
}





void miniJets::generateSamplingDataSets()
{
  double Tab = generateTimeDistribution( nEntries_collisionTimes );
  cout << "++++  Tab = " << Tab << "1/mb" << endl;

  double sigmaJet = generatePtDistribution( nEntries_PT, nEntries_PT_fine );
  cout << "++++  sigma_jet = " << sigmaJet << " mb" << endl;

  numberOfParticlesToGenerate = 2 * static_cast<int>( Tab * sigmaJet ) * numberOfTestparticles;
  cout << "++++ N = " << numberOfParticlesToGenerate / numberOfTestparticles << " * " << numberOfTestparticles
       << " testparticles" << endl;
}



void miniJets::populateParticleVector( std::vector< ParticleOffline >& _particles )
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


void miniJets::populateParticleVector( std::vector< ParticleOffline >& _particles, const int _numberOfParticlesToGenerate, const double _minimumPT )
{
  minimumPT = _minimumPT;
  numberOfParticlesToGenerate = _numberOfParticlesToGenerate;
  _particles.reserve( numberOfParticlesToGenerate * 1.5 );
  _particles.resize( numberOfParticlesToGenerate );
  sampleMomenta( _particles );
  samplePositions( _particles );
}



void miniJets::computeWoodSaxonParameters( const config& _config, WoodSaxon& _WoodSaxonParameter )
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
    throw eMiniJet_error( errMsg );
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



void miniJets::sampleMomenta( std::vector< ParticleOffline >& _particles )
{
  sample_PXY( _particles );
  sample_PZE( _particles );
  sample_FLAV( _particles );
}




void miniJets::samplePositions( std::vector< ParticleOffline >& _particles )
{
  sample_T( _particles );
  sample_XYZ( _particles );
}



void miniJets::sample_T( std::vector< ParticleOffline >& _particles ) const
{
  const int ord = 4;
  int nn;
  double *xx, *yy;
  double sample_y, t, dt;
  double xa[ord+1], ya[ord+1];

  std::fstream read( filename_samplingData_collisionTimes.c_str(), std::ios::in );

  xx = new double[nEntries_collisionTimes];
  yy = new double[nEntries_collisionTimes];
  for ( int i = 0; i < nEntries_collisionTimes; i++ )
  {
    read >> xx[i] >> yy[i];
  }

  for ( int index = 0; index < _particles.size(); index += 2 )
  {
    sample_y = ran2();

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

    for ( int i = 1; i <= ord; i++ )
    {
      ya[i] = yy[nn+i-1];
      xa[i] = xx[nn+i-1];
    }

    polint( ya, xa, ord, sample_y, &t, &dt );
    _particles[index].T = _particles[index+1].T = t;
  }

  delete[] xx;
  delete[] yy;
}



void miniJets::sample_XYZ( std::vector< ParticleOffline >& _particles ) const
{
  double X, Y, Z, T;
  double vt, gv, gvt, tc1, tc2, tc3, xx, bA, bB, zA, zB, fds, fd;
  double zmax, zmin, xmax, xmin, ymax, ymin, max;
  double c1, c2;

  for ( int index = 0; index < _particles.size(); index += 2 )
  {
    T = _particles[index].T;
    vt = velocity * fabs( T );
    gv = gamma * velocity;
    gvt = gv * fabs( T );

    c1 = sqrt(( RA0 + impactParameter ) * fabs( RA0 - impactParameter ) );
    tc1 = ( RA0 - c1 ) / ( 2 * gv );
    tc2 = sqrt( impactParameter * ( 2.0 * RA0 - impactParameter ) ) / ( 2 * gv );
    tc3 = ( RA0 + c1 ) / ( 2 * gv );

    xx = impactParameter / 2.0 + gvt * sqrt( pow( RA0, 2 ) / ( pow( impactParameter, 2 ) / 4.0  + pow( gvt, 2 ) ) - 1.0 );

    // z-axis in overlap region
    if (( RA0 > impactParameter ) && ( fabs( T ) >= tc1 ) && ( fabs( T ) <= tc3 ) )
    {
      zmin = vt - RA0 / gamma;
    }
    else
    {
      zmin = vt - sqrt(( RA0 + impactParameter - xx ) * ( RA0 - impactParameter + xx ) ) / gamma;
    }
    zmax = -zmin;

    // x-axis in overlap region
    if ( fabs( T ) <= tc2 )
    {
      xmax = RA0;
    }
    else
    {
      xmax = xx;
    }
    xmin = impactParameter - xmax;

    // y-axis in overlap region
    ymax = sqrt( pow( RA0, 2 ) - pow( impactParameter, 2 ) / 4.0 - pow( gvt, 2 ) );
    ymin = -ymax;

    // the maximum of nA(s,z-vt)*nB(s-b,z+vt) in the overlap region
    max = densityA( impactParameter / 2.0, -velocity * T ) * densityA( impactParameter / 2.0, velocity * T ) ;

    // sampling
    do
    {
      X = ( xmax - xmin ) * ran2() + xmin;
      Y = ( ymax - ymin ) * ran2() + ymin;
      Z = ( zmax - zmin ) * ran2() + zmin;

      bA = sqrt( pow( X, 2 ) + pow( Y, 2 ) );
      bB = sqrt( pow(( X - impactParameter ), 2 ) + pow( Y, 2 ) );

      if (( bA > RA0 ) || ( bB > RA0 ) )
      {
        fds = -1.0;
      }
      else
      {
        zA = sqrt(( RA0 + bA ) * ( RA0 - bA ) ) / gamma;
        zB = sqrt(( RA0 + bB ) * ( RA0 - bB ) ) / gamma;
        c1 = fabs( Z - velocity * T );
        c2 = fabs( Z + velocity * T );
        if (( c1 > zA ) || ( c2 > zB ) )
        {
          fds = -1.0;
        }
        else
        {
          fds = densityA( bA, c1 ) * densityA( bB, c2 );
        }
      }

      if ( fds > max )
        cout << "problem in sample_posit" << endl;

      fd = max * ran2();
    }
    while ( fds < fd );

    _particles[index].X = _particles[index+1].X = X;
    _particles[index].X = _particles[index].X - impactParameter / 2.0;    // shift into the correct coordinate system
    _particles[index+1].X = _particles[index+1].X - impactParameter / 2.0;    // shift into the correct coordinate system
    _particles[index].Y = _particles[index+1].Y = Y;
    _particles[index].Z = _particles[index+1].Z = Z;
  }
}



void miniJets::sample_PXY( std::vector< ParticleOffline >& _particles ) const
{
  const int ord = 4;
  int count, nn;
  double theta;
  double *xx, *yy, *xxFine, *yyFine;
  double sample_y, PT, dpt;
  double xa[ord+1], ya[ord+1];

  count = nEntries_PT;               //number of entries in samplPT
  xx = new double[count];
  yy = new double[count];

  std::fstream file_samplPT( filename_samplingData_PT.c_str(), std::ios::in );  //samplPT.dat gives PT and cross section for jet between P0 and PT
  for ( int i = 0; i < count; i++ )
  {
    file_samplPT >> xx[i] >> yy[i];
  }
  file_samplPT.close();

  for ( int j = 0; j < count; j++ )
  {
    yy[j] = yy[j] / yy[count-1];         //scales yy[j] to be in [0,1]; yy[count-1] is total cs for jet between P0 and Pmax
  }

  int countFine = nEntries_PT_fine;
  xxFine = new double[countFine];
  yyFine = new double[countFine];
  std::fstream file_samplPT_fine( filename_samplingData_PT_fine.c_str(), std::ios::in );
  for ( int i = 0; i < countFine; i++ )
  {
    file_samplPT_fine >> xxFine[i] >> yyFine[i];
  }
  file_samplPT_fine.close();


  for ( int index = 0; index < _particles.size(); index += 2 )
  {
    sample_y = ran2(); //a cross section is randomly chosen (within [0,1], i.e. relative to total c.s.)

    nn = 0;
    while ( yy[++nn] < sample_y );
    nn = nn - ord / 2;

    if ( nn < 0 )
    {
      nn = 0;
    }
    else if (( nn + ord - 1 > count - 1 ) || ( yy[nn+ord-1] == 1.0 ) )
    {
      while ( yy[++nn] < 1.0 );
      nn = nn - ord + 1;
    }

    for ( int j = 1; j <= ord; j++ )
    {
      ya[j] = yy[nn+j-1];
      xa[j] = xx[nn+j-1];
    }

    //interpolates between the values given in samplPT.dat and returns the PT corresponding to sample_y
    polint( ya, xa, ord, sample_y, &PT, &dpt );

    //---- for high-pt sample with higher resolution ----
    if ( PT > xxFine[0] )
    {
      sample_y = ran2();
      nn = 0;
      while ( yyFine[++nn] < sample_y );
      nn = nn - ord / 2;

      if ( nn < 0 )
      {
        nn = 0;
      }
      else if (( nn + ord - 1 > count - 1 ) || ( yyFine[nn+ord-1] == 1.0 ) )
      {
        while ( yyFine[++nn] < 1.0 );
        nn = nn - ord + 1;
      }
      for ( int j = 1; j <= ord; j++ )
      {
        ya[j] = yyFine[nn + j-1];
        xa[j] = xxFine[nn + j-1];
      }
      polint( ya, xa, ord, sample_y, &PT, &dpt );
    }
    //---------------------------------------------------

    theta = 2.0 * M_PI * ran2();
    _particles[index].PX = PT * cos( theta );
    _particles[index].PY = PT * sin( theta );
    _particles[index+1].PX = - _particles[index].PX;
    _particles[index+1].PY = - _particles[index].PY;
  }

  delete[] xx;
  delete[] yy;
  delete[] xxFine;
  delete[] yyFine;
}




void miniJets::sample_PZE( std::vector< ParticleOffline >& _particles ) const
{
  int count, nn;
  double *xx, *yy;
  double PT, Y1, Y2;
  double XT, Y1max, Y2max, Y2min, max, fg, fgs;

  std::fstream read( filename_samplingData_maximumPT.c_str(), std::ios::in );

  count = nEntries_PT;
  xx = new double[count];
  yy = new double[count];
  for ( int j = 0; j < count; j++ )
  {
    read >> xx[j] >> yy[j];
  }

  for ( int index = 0; index < _particles.size(); index += 2 )
  {
    nn = 0;
    PT = sqrt( pow( _particles[index].PX, 2 ) + pow( _particles[index].PY, 2 ) ) ;

    if ( PT > xx[count-1] )
    {
      std::string errMsg = "PT too large in sample_PZE";
      throw eMiniJet_error( errMsg );
    }

    while ( xx[++nn] <= PT );
    max = yy[nn - 1];

    XT = 2.0 * PT / sqrtS_perNN;
    Y1max = 1 / XT + sqrt( 1 / pow( XT, 2 ) - 1.0 );

    if ( Y1max <= 0 )
    {
      std::string errMsg = "problem1 in sample_Y12";
      throw eMiniJet_error( errMsg );
    }
    Y1max = log( Y1max );
    //-- set the rapidity interval --
    //   Yset=0.5;
    //   if(Y1max > Yset) Y1max=Yset;
    //-------------------------------

    //sampling via rejection method: comparison function is
    //f(y)=cs(y)=const=max, i.e. maximum value of cross section
    //for given PT when varying Y1 and Y2. Uniformly distributed
    //point is choosen: [Y,f(Y)] = [Ymin..Ymax , 0..max]. If
    //csjet(Y) < f(Y) the value for Y is accepted, otherwise it is rejected.

    do
    {
      Y1 = 2 * Y1max * ran2() - Y1max;
      Y2 = 2 * Y1max * ran2() - Y1max;

      Y2max = 2 / XT - exp( Y1 );
      Y2min = 2 / XT - exp( -Y1 );

      if (( Y2max <= 0 ) || ( Y2min <= 0 ) )
      {
        std::string errMsg = "problem1 in sample_Y12";
        throw eMiniJet_error( errMsg );
      }

      Y2max = log( Y2max );
      Y2min = -log( Y2min );

      if (( Y2 > Y2max ) || ( Y2 < Y2min ) )
      {
        fgs = -1.0;
      }
      else
      {
        fgs = csjet( sqrtS_perNN, PT, Y1, Y2 );
      }
      if ( fgs > max )
      {
        std::string errMsg = "problem1 in sample_Y12";
        throw eMiniJet_error( errMsg );
      }

      fg = max * ran2();
    }
    while ( fgs < fg );

    _particles[index].PZ = 0.5 * PT * ( exp( Y1 ) - exp( -Y1 ) );
    _particles[index+1].PZ = 0.5 * PT * ( exp( Y2 ) - exp( -Y2 ) );
    _particles[index].E = sqrt( pow( _particles[index].PZ, 2 ) + pow( PT, 2 ) );
    _particles[index+1].E = sqrt( pow( _particles[index+1].PZ, 2 ) + pow( PT, 2 ) );
    _particles[index].m = ParticleOffline::getMass( _particles[index].FLAVOR );
    _particles[index+1].m = ParticleOffline::getMass( _particles[index+1].FLAVOR );
  }

  delete[] xx;
  delete[] yy;
}




// elementary PQCD cross section
#define gggg(t,u) (9./2.*(3.-t*u-u/(t*t)-t/(u*u)))
#define ggqq(t,u) (1./6.*(t/u+u/t)-3./8.*(t*t+u*u))
#define qqgg(t,u) (64./9.*ggqq(t,u))
#define gqgq(t,u) (-4./9.*(1./u+u)+(1.+u*u)/(t*t))
#define qq12(t,u) (4./9.*(1.+u*u)/(t*t))
#define qqqq(t,u) (4./9.*((1.+u*u)/(t*t)+(1.+t*t)/(u*u))-8./27./t/u)
#define qaqs(t,u) (4./9.*((1.+u*u)/(t*t)+u*u+t*t)-8./27.*u*u/t)
#define qaqd(t,u) (4./9.*(t*t+u*u))

void miniJets::sample_FLAV( std::vector< ParticleOffline >& _particles ) const
{
  double XT, PT2, Y1, Y2, x1, x2, s, t, u;
  double F1[9], F2[9];
  double gg1, gg2, gg, gq, qqb, qq1, qq2, qq;
  double total, sf, tmp, c0, c1, c2, c3, c4, c5, d0, d1, d2, d3;
  double qqbp[5], gqp[9], qq2p[4];
  int ni, nn, n1, n2;
  unsigned int f1, f2;

  void xpDO( double Q2, double x1, double x2, int n1, int n2, double F1[9], double F2[9] );
  void xpGRV( double Q2, double x1, double x2, int n1, int n2, double F1[9], double F2[9] );

  for ( int index = 0; index < _particles.size(); index += 2 )
  {
    PT2 = pow( _particles[index].PX, 2 ) + pow( _particles[index].PY, 2 );
    XT = 2.0 * sqrt( PT2 ) / sqrtS_perNN;
    Y1 = 0.5 * log(( _particles[index].E + _particles[index].PZ ) / ( _particles[index].E-_particles[index].PZ ) );
    Y2 = 0.5 * log(( _particles[index+1].E + _particles[index+1].PZ ) / ( _particles[index+1].E-_particles[index+1].PZ ) );
    x1 = 0.5 * XT * ( exp( Y1 ) + exp( Y2 ) );
    x2 = 0.5 * XT * ( exp( -Y1 ) + exp( -Y2 ) );

    if ((( 1.0 - x1 ) < 0.0001 ) || (( 1.0 - x2 ) < 0.0001 ) )
    {
      std::string errMsg = "Problem in sample_FLAV()";
      throw eMiniJet_error( errMsg );
    }

    //The nucleons are sampled to being either proton(1) or neutron(-1).
    if ( ran2() < ( Aatomic / A ) )
    {
      n1 = 1;
    }
    else
    {
      n1 = -1;
    }
    if ( ran2() < ( Batomic / B ) )
    {
      n2 = 1;
    }
    else
    {
      n2 = -1;
    }

    // Parton distribution function
    xpGRV( PT2, x1, x2, n1, n2, F1, F2 );

    s = x1 * x2 * pow( sqrtS_perNN, 2 );
    t = -PT2 * ( 1. + exp( Y2 - Y1 ) );
    u = -PT2 * ( 1. + exp( Y1 - Y2 ) );
    t = t / s;
    u = u / s;

    // 0->g,1->u,2->ub,3->d,4->db,5->s,6->sb,7->c,8->cb
    gg1 = F1[0] * F2[0] * gggg( t, u );
    gg2 = ( F1[1] * F2[2] + F1[3] * F2[4] + F1[5] * F2[6] + F1[7] * F2[8] + F2[1] * F1[2] + F2[3] * F1[4] + F2[5] * F1[6] + F2[7] * F1[8] ) * qqgg( t, u );
    gg = gg1 + gg2;

    c0 = gqgq( t, u ) + gqgq( u, t );
    gq = 0.0;
    for ( int j = 1; j <= 8; j++ )
    {
      gqp[j] = ( F1[0] * F2[j] + F2[0] * F1[j] ) * c0;
      gq += gqp[j];
    }

    c1 = qaqs( t, u ) + qaqs( u, t );
    c2 = F1[1] * F2[2] + F1[3] * F2[4] + F1[5] * F2[6] + F1[7] * F2[8]
         + F2[1] * F1[2] + F2[3] * F1[4] + F2[5] * F1[6] + F2[7] * F1[8];
    c3 = qaqd( t, u ) + qaqd( u, t );
    c4 = F1[0] * F2[0] * ( ggqq( t, u ) + ggqq( u, t ) );

    qqb = 0.0;
    for ( int i = 1; i <= 4; i++ )
    {
      qqbp[i] = F1[2*i-1] * F2[2*i] + F2[2*i-1] * F1[2*i];
      qqbp[i] = qqbp[i] * c1 + ( c2 - qqbp[i] ) * c3 + c4;
      qqb += qqbp[i];
    }

    qq1 = ( F1[1] * F2[1] + F1[2] * F2[2] + F1[3] * F2[3] + F1[4] * F2[4] + F1[5] * F2[5]
            + F1[6] * F2[6] + F1[7] * F2[7] + F1[8] * F2[8] ) * qqqq( t, u );
    c5 = qq12( t, u ) + qq12( u, t );
    qq2p[1] = (( F1[1] + F1[2] ) * ( F2[3] + F2[4] + F2[5] + F2[6] + F2[7] + F2[8] )
               + ( F2[1] + F2[2] ) * ( F1[3] + F1[4] + F1[5] + F1[6] + F1[7] + F1[8] ) ) * c5;
    qq2p[2] = (( F1[3] + F1[4] ) * ( F2[5] + F2[6] + F2[7] + F2[8] )
               + ( F2[3] + F2[4] ) * ( F1[5] + F1[6] + F1[7] + F1[8] ) ) * c5;
    qq2p[3] = (( F1[5] + F1[6] ) * ( F2[7] + F2[8] )
               + ( F2[5] + F2[6] ) * ( F1[7] + F1[8] ) ) * c5;
    qq2 = qq2p[1] + qq2p[2] + qq2p[3];
    qq = qq1 + qq2;

    total = gg + gq + qqb + qq;
    sf = total * ran2();

    // gg,qqb->gg
    if ( sf <= gg )
    {
      f1 = f2 = 0;
    }
    // gq->gq
    else if ( sf <= ( gg + gq ) )
    {
      sf = sf - gg;
      tmp = 0.0;
      for ( int j = 1; j <= 8; j++ )
      {
        if ( gqp[j] > 1.0e-8 )
        {
          f1 = j;// Vorsichtsausnahme wegen
          f2 = 0;// "Rundungsfehler"
        }
        tmp += gqp[j];
        // q?=u,ub,d,db,s,sb,c,cb
        if ( sf <= tmp )
        {
          sf = sf - tmp + gqp[j];
          tmp = F1[0] * F2[j] * gqgq( t, u ) + F2[0] * F1[j] * gqgq( u, t );
          // 1->g, 2->q
          if ( sf <= tmp )
          {
            f1 = 0;
            f2 = j;
          }
          // 1->q, 2->g
          else
          {
            f1 = j;
            f2 = 0;
          }
          j = 8;
        }
      }
    }
    // gg,qqb,q'q'b->qqb
    else if ( sf <= ( gg + gq + qqb ) )
    {
      sf = sf - gg - gq;
      tmp = 0.0;
      for ( int i = 1; i <= 4; i++ )
      {
        if ( qqbp[i] > 1.0e-8 )
        {
          f1 = 2 * i;  // Vorsichtsausnahme wegen
          f2 = 2 * i - 1;// "Rundungsfehler"
        }
        tmp += qqbp[i];
        // qqb?=uub,ddb,ssb,ccb
        if ( sf <= tmp )
        {
          sf = sf - tmp + qqbp[i];
          d1 = F1[2*i-1] * F2[2*i];
          d2 = F2[2*i-1] * F1[2*i];
          tmp = d1 * qaqs( t, u ) + d2 * qaqs( u, t ) + 0.5 * (( c2 - d1 - d2 ) * c3 + c4 );
          // 1->q, 2->qb
          if ( sf <= tmp )
          {
            f1 = 2 * i - 1;
            f2 = f1 + 1;
          }
          // 1->qb, 2->q
          else
          {
            f1 = 2 * i;
            f2 = f1 - 1;
          }
          i = 4;
        }
      }
    }
    // qq'->qq'
    else
    {
      sf = sf - gg - gq - qqb;
      // q=q'
      if ( sf < qq1 )
      {
        tmp = 0.0;
        for ( int i = 1;i <= 8;i++ )
        {
          d0 = F1[i] * F2[i] * qqqq( t, u );
          if ( d0 > 1.0e-8 )
          {
            f1 = f2 = i;// Vorsichtsausnahme wegen "Rundungsfehler"
          }
          tmp += d0;
          // q?=u,ub,d,db,s,sb,c,cb
          if ( sf <= tmp )
          {
            f1 = f2 = i;
            i = 8;
          }
        }
      }
      // q!=q'
      else
      {
        sf = sf - qq1;
        tmp = 0.0;
        for ( int j = 1;j <= 3;j++ )
        {

          if ( qq2p[j] > 1.0e-8 )
          {
            for ( int n = 0;n <= 1;n++ )
            {
              nn = 2 * j - n;
              for ( int m = 8;m >= 2*j + 1;m-- )
              {
                d1 = F1[nn] * F2[m];
                d2 = F2[nn] * F1[m];
                d3 = ( d1 + d2 ) * c5;
                if ( d3 > 1.0e-8 )
                {
                  f1 = m; // Vorsichtsausnahme wegen
                  f2 = nn;// "Rundungsfehler"
                  m = 0;
                  n = 2;
                }
              }
            }
          }

          tmp += qq2p[j];
          // q?= u(ub),d(db),s(sb)
          if ( sf <= tmp )
          {
            sf = sf - tmp + qq2p[j];
            d1 = d2 = 0.0;
            for ( int i = 2 * j + 1;i <= 8;i++ )
            {
              d1 += F2[i];
              d2 += F1[i];
            }
            tmp = ( F1[2*j-1] * d1 + F2[2*j-1] * d2 ) * c5;
            // q=u,d or s
            if ( sf <= tmp )
              nn = 2 * j - 1;
            // q=ub,db or sb
            else
            {
              nn = 2 * j;
              sf = sf - tmp;
            }
            tmp = 0.0;
            for ( int l = 2 * j + 1;l <= 8;l++ )
            {
              d1 = F1[nn] * F2[l];
              d2 = F2[nn] * F1[l];
              d3 = ( d1 + d2 ) * c5;
              if ( d3 > 1.0e-8 )
              {
                f1 = l; // Vorsichtsausnahme wegen
                f2 = nn;// "Rundungsfehler"
              }
              tmp += d3;
              // q'?= d,db,s,sb,c,cb
              if ( sf <= tmp )
              {
                sf = sf - tmp + d3;
                // 1->q, 2->q'
                if ( sf <= ( d1*qq12( t, u ) + d2*qq12( u, t ) ) )
                {
                  f1 = nn;
                  f2 = l;
                }
                // 1->q', 2->q
                else
                {
                  f1 = l;
                  f2 = nn;
                }
                l = 8;
              }
            }
            j = 3;
          }
        }
      }
    }
    if (( f1 > 8 ) || ( f1 < 0 ) )
    {
      std::string errMsg = "problem in sample_FLAV";
      throw eMiniJet_error( errMsg );
    }
    if (( f2 > 8 ) || ( f2 < 0 ) )
    {
      std::string errMsg = "problem in sample_FLAV";
      throw eMiniJet_error( errMsg );
    }
    
    //!! temporary hack: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    double random = ran2();
    if( random < 0.6 )
    {
      f1 = 7;
      f2 = 8;
      _particles[index].N_EVENT_pp = index+1;
    _particles[index+1].N_EVENT_pp = index+1;
    }
//     else if( random < 0.9 )
// //     else
//     {
//       f1 = 9;
//       f2 = 10;
//     }
    else
    {
      f1 = 50;
      f2 = 50;
      _particles[index].N_EVENT_pp = index+1;
      _particles[index+1].N_EVENT_pp = index+2;
    }
//     else
//     {
//       f1 = 7;
//       f2 = 8;
//     }
    
    _particles[index].FLAVOR = static_cast<FLAVOR_TYPE>( f1 );
    _particles[index+1].FLAVOR = static_cast<FLAVOR_TYPE>( f2 );
    
    
    //!! temporary hack: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    _particles[index].m = ParticleOffline::getMass( _particles[index].FLAVOR );
    _particles[index+1].m = ParticleOffline::getMass( _particles[index+1].FLAVOR );
    
    _particles[index].E = sqrt( pow( _particles[index].E, 2 ) + pow( _particles[index].m, 2 ) );
    _particles[index+1].E = sqrt( pow( _particles[index+1].E, 2 ) + pow( _particles[index+1].m, 2 ) );
    
    
    
    
    
  }
}



bool miniJets::samplingDataSetsExist() const
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

  //--------------------------------  check file 2 --------------------------------
  testFile.clear();
  testFile.open( filename_samplingData_PT.c_str(), std::ios::in );
  if ( testFile.good() )  //check if quark file actually exists
  {
    testFile.close();
  }
  else
  {
    exists = false;
    cout << "Could not open: " << filename_samplingData_PT << endl;
  }
  //--------------------------------  check file 2 --------------------------------

  //--------------------------------  check file 3 --------------------------------
  testFile.clear();
  testFile.open( filename_samplingData_PT_fine.c_str(), std::ios::in );
  if ( testFile.good() )  //check if quark file actually exists
  {
    testFile.close();
  }
  else
  {
    exists = false;
    cout << "Could not open: " << filename_samplingData_PT_fine << endl;
  }
  //--------------------------------  check file 3 --------------------------------

  //--------------------------------  check file 4 --------------------------------
  testFile.clear();
  testFile.open( filename_samplingData_maximumPT.c_str(), std::ios::in );
  if ( testFile.good() )  //check if quark file actually exists
  {
    testFile.close();
  }
  else
  {
    exists = false;
    cout << "Could not open: " << filename_samplingData_maximumPT << endl;
  }
  //--------------------------------  check file 4 --------------------------------

  return exists;
}




//-----------------------------------------------//
// time distribution, calculation runs only once.//
//-----------------------------------------------//
double miniJets::generateTimeDistribution( int& _count )
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



//-----------------------------------------------------------//
// int_P0^PT{d(sigma)/d(PT)}, the calculation runs only once.//
//-----------------------------------------------------------//
double miniJets::generatePtDistribution( int& count, int& count_fine )
{
  double PT;
  double factor, end, dpt, old, c, s, y;
  double tgral, sd, chi2a, jetdist;
  double sigma_jet;

  std::list<double> fine;
  std::list<double> finePT;
  double fineStart = 0;

  std::ofstream file_samplPT( filename_samplingData_PT.c_str(), std::ios::out | std::ios::trunc );
  std::ofstream file_samplPT_fine( filename_samplingData_PT_fine.c_str(), std::ios::out | std::ios::trunc );
  std::ofstream file_ptdist( "data/ptdist.dat", std::ios::out | std::ios::trunc );
  std::ofstream file_maxPT( filename_samplingData_maximumPT.c_str(), std::ios::out | std::ios::trunc );
  std::ofstream file_maxYPT( "data/maxYPT.dat", std::ios::out | std::ios::trunc );
  file_samplPT.precision( 12 );
  file_ptdist.precision( 12 );
  file_maxPT.precision( 12 );
  file_maxYPT.precision( 12 );
  file_samplPT_fine.precision( 12 );
  string sep = "\t";

  integrand_distPT fdistPT;
  fdistPT.setSqrtS( sqrtS_perNN );

  factor = 1. / 2.5682;// transform dimension from 1/GeV^2 to mb
  factor = factor * K;

  end = 35.0;//large PT cut < sqrtS/2
  dpt = 0.01;

  old = infinity;
  s = 0.0;
  c = 0.0;
  count = 0;
  count_fine = 0;

  PT = P0;
  do
  {
    fdistPT.setPt( PT );
    vegas( 2, fdistPT, &tgral, &sd, &chi2a );  //integrate dsigma/(dPT^2*dY1*dY2) over Y1 and Y2
    jetdist = double( tgral ) * factor * 2.0 * PT;

    ++count;
    s += c * 0.5 * ( old + jetdist ) * dpt;        //cross section for a jet between P0 and the current PT

    if ( PT >= 11.0 )
    {
      if ( fine.empty() )
        fineStart = s / 2.0;
      fine.push_back( s / 2.0 - fineStart );
      finePT.push_back( PT );
      ++count_fine;
    }

    file_samplPT << PT << sep << s / 2.0 << endl;
    // file_samplPT << sd << sep << chi2a;
    file_ptdist << PT << sep << log( jetdist ) << endl;
    file_maxPT << PT << sep << max( PT, y ) << endl;  //max(PT,y) gives maximum value of dsigma for given PT
    file_maxYPT << PT << sep << y << endl;  //y corresponding to the max. value of dsigma for given PT, calc. in max(PT,y)

    if (( old - jetdist ) < 0.001 )
    {
      dpt = dpt * 10.;
      dpt = ( dpt < 1.0 ) ? dpt : 1.0;
    }
    old = jetdist;
    c = 1.0;

    PT += dpt;
  }
  while ( PT <= end );

  std::list<double>::const_iterator cIt;
  std::list<double>::const_iterator cItLabels = finePT.begin();
  for ( cIt = fine.begin(); cIt != fine.end(); ++cIt )
  {
    file_samplPT_fine << *cItLabels << sep << ( *cIt / fine.back() ) << endl;
    ++cItLabels;
  }

  sigma_jet = s / 2.0;
  return sigma_jet;
}





void integrand_distPT::operator()( const int *ndim, const double xx[], const int *ncomp, double ff[] ) const
{
  double wgt;
  ff[0] = this->operator()( xx, wgt );
}


double integrand_distPT::operator()( const double x[], double wgt ) const
{
  //Y1 and Y2 run from min[Yi] to max[Yi] (i=1,2) as x[1] and x[2] run from 0 to 1
  double XT, Y1max, Y2max, Y2min, Y1, Y2, V;

  XT = 2.0 * pt / sqrtS;
  Y1max = 1.0 / XT + sqrt( 1.0 / XT / XT - 1.0 );

  if ( Y1max <= 0. )
  {
    cout << "problem0 in fdistPT" << endl;
    return 0.0;
  }

  Y1max = log( Y1max );
  V = Y1max;
  Y1 = V * double( x[1] );

  Y2max = 2.0 / XT - exp( Y1 );
  Y2min = 2.0 / XT - exp( -Y1 );

  if (( Y2max <= 0. ) || ( Y2min <= 0. ) )
  {
    cout << "problem1 in fdistPT" << endl;
    return 0.0;
  }

  Y2max = log( Y2max );
  Y2min = -log( Y2min );

  Y2 = ( Y2max - Y2min ) * double( x[2] ) + Y2min;

  V = 2.0 * V * ( Y2max - Y2min );//factor 2.0 due to the symmetry in Y1

  return double( V*csjet( sqrtS, pt, Y1, Y2 ) );
}


//----------------------------------------------------------//
// to give the maximum value of G~d(sigma)/d(pt2)d(y1)d(y2) //
// for a fixed PT for preparing the sampling in y1-y2 plane //
//----------------------------------------------------------//
double miniJets::max( double PT, double& y )
{
  const int nm = 1000;
  double XT, Ymax, Y1, Y2, tmp, maxx;

  XT = 2.0 * PT / sqrtS_perNN;
  Ymax = 1.0 / XT + sqrt( 1.0 / XT / XT - 1.0 );
  if ( Ymax <= 0. )
  {
    cout << "problem in max" << endl;
    return 0.0;
  }
  Ymax = log( Ymax );

  maxx = csjet( sqrtS_perNN, PT, 0.0, 0.0 );
  y = 0.0;
  for ( int i = 0;i < nm;i++ )
  {
    Y1 = Ymax * ran2();
    Y2 = -Y1;
    tmp = csjet( sqrtS_perNN, PT, Y1, Y2 );
    if ( tmp > maxx )
    {
      maxx = tmp;
      y = Y1;
    }
  }
  return 1.03*maxx;// 3% uncertainty
}



//---------------------------------------------//
// the cross section d(sigma)/d(PT2)d(y1)d(y2) //
//---------------------------------------------//
double csjet( double sqrtS, double PT, double Y1, double Y2 )
{

  const int flav = 4;//u,d,s,c Quark

  double XT, PT2, x1, x2, s, t, u, alpha_s, factor;
  double gg1, gg2, gq, qq1, qq2, qq3, qq4, qq5, gg, qq;
  double F1[9], F2[9];

  XT = 2.0 * PT / sqrtS;
  PT2 = PT * PT;
  x1 = 0.5 * XT * ( exp( Y1 ) + exp( Y2 ) );
  x2 = 0.5 * XT * ( exp( -Y1 ) + exp( -Y2 ) );

  if ((( 1.0 - x1 ) < 0.0001 ) || (( 1.0 - x2 ) < 0.0001 ) ) return 0.0;
  // Parton distribution function
  //   xpDO(PT2,x1,x2,1,1,F1,F2);
  //   F1[7]=F1[8]=F2[7]=F2[8]=0.0;

  xpGRV( PT2, x1, x2, 1, 1, F1, F2 );// for a proton-proton-collision, PDFs according to Glück-Reya-Vogt parametrisation

  s = x1 * x2 * sqrtS * sqrtS;
  t = -PT2 * ( 1. + exp( Y2 - Y1 ) );
  u = -PT2 * ( 1. + exp( Y1 - Y2 ) );
  t = t / s;
  u = u / s;

  gg1 = F1[0] * F2[0] * gggg( t, u );

  gg2 = ( F1[1] * F2[2] + F1[3] * F2[4] + F1[5] * F2[6] + F1[7] * F2[8]
          + F2[1] * F1[2] + F2[3] * F1[4] + F2[5] * F1[6] + F2[7] * F1[8] ) * qqgg( t, u );

  gq = ( F1[0] * ( F2[1] + F2[2] + F2[3] + F2[4] + F2[5] + F2[6] + F2[7] + F2[8] ) +
         + F2[0] * ( F1[1] + F1[2] + F1[3] + F1[4] + F1[5] + F1[6] + F1[7] + F1[8] ) ) * ( gqgq( t, u ) + gqgq( u, t ) );

  qq1 = (( F1[1] + F1[2] ) * ( F2[3] + F2[4] + F2[5] + F2[6] + F2[7] + F2[8] )
         + ( F1[3] + F1[4] ) * ( F2[5] + F2[6] + F2[7] + F2[8] ) + ( F1[5] + F1[6] ) * ( F2[7] + F2[8] )
         + ( F2[1] + F2[2] ) * ( F1[3] + F1[4] + F1[5] + F1[6] + F1[7] + F1[8] )
         + ( F2[3] + F2[4] ) * ( F1[5] + F1[6] + F1[7] + F1[8] ) + ( F2[5] + F2[6] ) * ( F1[7] + F1[8] ) )
        * ( qq12( t, u ) + qq12( u, t ) );

  qq2 = ( F1[1] * F2[1] + F1[2] * F2[2] + F1[3] * F2[3] + F1[4] * F2[4] + F1[5] * F2[5] + F1[6] * F2[6]
          + F1[7] * F2[7] + F1[8] * F2[8] ) * qqqq( t, u );

  qq3 = ( F1[1] * F2[2] + F1[3] * F2[4] + F1[5] * F2[6] + F1[7] * F2[8]
          + F2[1] * F1[2] + F2[3] * F1[4] + F2[5] * F1[6] + F2[7] * F1[8] ) * ( qaqs( t, u ) + qaqs( u, t ) );

  qq4 = ( flav - 1 ) * ( F1[1] * F2[2] + F1[3] * F2[4] + F1[5] * F2[6] + F1[7] * F2[8]
                         + F2[1] * F1[2] + F2[3] * F1[4] + F2[5] * F1[6] + F2[7] * F1[8] ) * ( qaqd( t, u ) + qaqd( u, t ) );

  qq5 = flav * F1[0] * F2[0] * ( ggqq( t, u ) + ggqq( u, t ) );

  gg = gg1 + gg2;
  gq = gq;
  qq = qq1 + qq2 + qq3 + qq4 + qq5;

  alpha_s = 12.*M_PI / ( 33. - 2.*flav ) / log( PT2 / lambda2 );
  factor = M_PI * alpha_s * alpha_s / ( s * s );

  return factor*( gg + gq + qq );
}




//--------------------//
//   interpolation    //
//--------------------//
void miniJets::polint( const double xa[], const double ya[], const int n, const double x, double *y, double *dy ) const
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





double miniJets::densityA( double b, double z ) const
{
  return gamma*n0A / ( exp(( sqrt( b*b + pow( double( gamma*z ), 2.0 ) ) - RA ) / dA ) + 1.0 );
}
// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;
