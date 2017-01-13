//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL: svn+ssh://senzel@th.physik.uni-frankfurt.de/home/bamps/svn/full/offlineAnalysis/trunk/src/initialmodel_pythiashower.cpp $
//$LastChangedDate: $
//$LastChangedRevision: -1 $
//$LastChangedBy: $
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
#include "initialmodel_pythiashower.h"
#include "random.h"
#include "particle.h"
#include <binning.h>
#include "FPT_compare.h"

using namespace std;

extern "C" {
  void fixedshowerevent_( const double* px, const double* py, const double* pz, int* pythiaFlavor1, int* pythiaFlavor2, uint32_t* seed );
  void pythiashowerevent_( const double* ptMin, uint32_t* seed );
  void photonshowerevent_( const double* ptMin, uint32_t* seed );
  void hqshowerevent_( const double* ptMin, uint32_t* seed, const double* m_charm, const double* m_bottom, const int* heavy_jet_flavor );
  struct
  {
    double pa[100][7];
  } bamps_;
}


initialModel_PYTHIAShower::initialModel_PYTHIAShower( const config& _config, WoodSaxon& _WoodSaxonParameter, const SHOWER_TYPE _shower_type, const double _minimumPT, const FLAVOR_TYPE _initialPartonFlavor ):
  initialModelWS(_config),
  nEventsToGenerate( _config.getNaddedEvents() ),
  seed( _config.getSeed() ),
  filename_prefix( _config.getStandardOutputDirectoryName() + "/" + _config.getJobName() ),
  shower_type( _shower_type ),
  initialPartonFlavor( _initialPartonFlavor )
{
  impactParameter = _config.getImpactParameter();
  sqrtS_perNN = _config.getSqrtS();
  P0 = _minimumPT;
  if (!WoodSaxonParameter.Calculate( A, impactParameter, sqrtS_perNN))
  {
    std::string errMsg = "Impact parameter b too large. b > 2 R_A0";
    throw ePYTHIAShower_error( errMsg );
  }
  _WoodSaxonParameter = WoodSaxonParameter;

  double Tab;
  cout << "======= Generating data sets for sampling of initial state =======" << endl;

  generateTimeDistributionWS(Tab);
  cout << "++++  Tab = " << Tab << "1/mb" << endl;
  cout << "==================================================================" << endl;
}


void initialModel_PYTHIAShower::populateParticleVector( std::vector< Particle >& _particles )
{
  switch( shower_type )
  {
    case fixedShower:
    {
      vector<Particle> tempParticles;

      string filename = filename_prefix + "_" + "unshoweredParticles";
      fstream file( filename.c_str(), ios::out | ios::trunc );
      string sep = "\t";

      file << "#event" << sep << "px" << sep << "py" << sep << "pz" << sep << "E" << sep << "x"
           << sep << "y" << sep << "z" << sep << "t" << endl;

      for( int n = 0; n < nEventsToGenerate; n++ )
      {
        vector<Particle> particleShower;
        double px = P0;
        double py = 0.0;
        double pz = 0.0;
        
        FLAVOR_TYPE flavorA, flavorB;
        switch( initialPartonFlavor )
        {
          case 0:
            flavorA = gluon;
            flavorB = gluon;
            break;
          case 1:
            flavorA = up;
            flavorB = anti_up;
            break;
          case 2:
            flavorA = gluon;
            flavorB = up;
            break;
          default:
            throw ePYTHIAShower_error( "Unknown initial parton flavor. Unrecoverable error!");
        }
      
        particleShower = getFixedShowerEvent( px, py, pz, flavorA, flavorB );

        sample_TXYZ_singleParticle( particleShower[0] );
 
        for( int j = 0; j < particleShower.size(); j++ )
        {
          particleShower[j].N_EVENT_pp = n;
          particleShower[j].N_EVENT_AA = n;
          particleShower[j].Pos.X() = particleShower[0].Pos.X();
          particleShower[j].Pos.Y() = particleShower[0].Pos.Y();
          particleShower[j].Pos.Z() = particleShower[0].Pos.Z();
          particleShower[j].Pos.T() = particleShower[0].Pos.T();

          particleShower[j].init = true;

          tempParticles.push_back( particleShower[j] );
        }
        
        file << n << sep << P0 << sep << 0.0 << sep << 0.0 << sep << P0
          << sep << particleShower[0].Pos.X() << sep << particleShower[0].Pos.Y() << sep << particleShower[0].Pos.Z() << sep << particleShower[0].Pos.T()
          << sep << flavorA << endl;
        file << n << sep << -P0 << sep << 0.0 << sep << 0.0 << sep << P0
          << sep << particleShower[0].Pos.X() << sep << particleShower[0].Pos.Y() << sep << particleShower[0].Pos.Z() << sep << particleShower[0].Pos.T()
          << sep << flavorB << endl;

      }
      _particles.clear();
      _particles = tempParticles;
      break;
    }
    case pythiaShower:
    {
      binningLog gluonBinning( filename_prefix + "_unshowered_gluons.dat", 30.0, 55.0, 25 );
      binningLog quarkBinning( filename_prefix + "_unshowered_quarks.dat", 30.0, 55.0, 25 );
      string filename = filename_prefix + "_" + "unshoweredParticles";
      fstream file( filename.c_str(), ios::out | ios::trunc );
      string sep = "\t";

      file << "#event" << sep << "px" << sep << "py" << sep << "pz" << sep << "E" << sep << "x"
           << sep << "y" << sep << "z" << sep << "t" << endl;

      for ( int n = 0; n < nEventsToGenerate; n++ )
      {
        vector< Particle > tempParticles, initialPartons;
        getPythiaShowerEvent( tempParticles, initialPartons );
      
        sample_TXYZ_singleParticle( tempParticles[0] );
        for ( int i = 0; i < tempParticles.size(); i++ )
        {
          tempParticles[i].Pos.T() = tempParticles[0].Pos.T();
          tempParticles[i].Pos.X() = tempParticles[0].Pos.X();
          tempParticles[i].Pos.Y() = tempParticles[0].Pos.Y();
          tempParticles[i].Pos.Z() = tempParticles[0].Pos.Z();
          tempParticles[i].N_EVENT_pp = n;
          tempParticles[i].N_EVENT_AA = n;
          _particles.push_back( tempParticles[i] );
        }
      
        for( int i = 0; i < initialPartons.size(); i++ )
        {
          file << n << sep << initialPartons[i].Mom.Px() << sep << initialPartons[i].Mom.Py() << sep << initialPartons[i].Mom.Pz() << sep << initialPartons[i].Mom.E() 
          << sep << tempParticles[0].Pos.X() << sep << tempParticles[0].Pos.Y() << sep << tempParticles[0].Pos.Z() << sep << tempParticles[0].Pos.T()
          << sep << initialPartons[i].FLAVOR << endl;
          
          if ( abs( initialPartons[i].Mom.Pseudorapidity(0.0) ) < 0.8 )
          {
            if ( initialPartons[i].FLAVOR == gluon )
              gluonBinning.add( initialPartons[i].Mom.Pt() );
            else
              quarkBinning.add( initialPartons[i].Mom.Pt() );
          }
        }
      }
      file.close();
      gluonBinning.print();
      quarkBinning.print();
      break;
    }
    case fixedParton:
    {
      vector<Particle> tempParticles;

      string filename = filename_prefix + "_" + "unshoweredParticles";
      fstream file( filename.c_str(), ios::out | ios::trunc );
      string sep = "\t";

      file << "#event" << sep << "px" << sep << "py" << sep << "pz" << sep << "E" << sep << "x"
           << sep << "y" << sep << "z" << sep << "t" << endl;

      Particle initialPartonA, initialPartonB;
        
      sample_TXYZ_singleParticle( initialPartonA );
      initialPartonB = initialPartonA;
      initialPartonA.Mom.E() = initialPartonB.Mom.E() = P0;
      initialPartonA.Mom.Px() = P0;
      initialPartonB.Mom.Px() = -P0;
      initialPartonA.Mom.Py() = initialPartonB.Mom.Py() = 0.0;
      initialPartonA.Mom.Pz() = initialPartonB.Mom.Pz() = 0.0;
      
      FLAVOR_TYPE flavorA, flavorB;
      switch( initialPartonFlavor )
      {
        case 0:
          initialPartonA.FLAVOR = gluon;
          initialPartonB.FLAVOR = gluon;
          break;
        case 1:
          initialPartonA.FLAVOR = up;
          initialPartonB.FLAVOR = anti_up;
          break;
        case 2:
          initialPartonA.FLAVOR = gluon;
          initialPartonB.FLAVOR = up;
          break;
        default:
          throw ePYTHIAShower_error( "Unknown initial parton flavor. Unrecoverable error!");
      }
      
      for( int n = 0; n < nEventsToGenerate; n++ )
      {
        initialPartonA.N_EVENT_pp = initialPartonB.N_EVENT_pp = n;
        initialPartonA.N_EVENT_AA = initialPartonB.N_EVENT_AA = n;

        initialPartonA.init = initialPartonB.init = true;

        tempParticles.push_back( initialPartonA );
        tempParticles.push_back( initialPartonB );
        
        file << n << sep << P0 << sep << 0.0 << sep << 0.0 << sep << P0
          << sep << initialPartonA.Pos.X() << sep << initialPartonA.Pos.Y() << sep << initialPartonA.Pos.Z() << sep << initialPartonA.Pos.T()
          << sep << flavorA << endl;
        file << n << sep << -P0 << sep << 0.0 << sep << 0.0 << sep << P0
          << sep << initialPartonB.Pos.X() << sep << initialPartonB.Pos.Y() << sep << initialPartonB.Pos.Z() << sep << initialPartonB.Pos.T()
          << sep << flavorB << endl;
      }

      _particles.clear();
      _particles = tempParticles;
      break;
    }
    case photonShower:
    {
      string filename = filename_prefix + "_" + "unshoweredParticles";
      fstream file( filename.c_str(), ios::out | ios::trunc );
      string sep = "\t";

      file << "#event" << sep << "px" << sep << "py" << sep << "pz" << sep << "E" << sep << "x"
           << sep << "y" << sep << "z" << sep << "t" << endl;

      for ( int n = 0; n < nEventsToGenerate; n++ )
      {
        vector< Particle > tempParticles, initialPartons;
        getPhotonShowerEvent( tempParticles, initialPartons );
      
        sample_TXYZ_singleParticle( tempParticles[0] );
        for ( int i = 0; i < tempParticles.size(); i++ )
        {
          tempParticles[i].Pos.T() = tempParticles[0].Pos.T();
          tempParticles[i].Pos.X() = tempParticles[0].Pos.X();
          tempParticles[i].Pos.Y() = tempParticles[0].Pos.Y();
          tempParticles[i].Pos.Z() = tempParticles[0].Pos.Z();
          tempParticles[i].N_EVENT_pp = n;
          tempParticles[i].N_EVENT_AA = n;
          _particles.push_back( tempParticles[i] );
        }
      
        for( int i = 0; i < initialPartons.size(); i++ )
        {
          file << n << sep << initialPartons[i].Mom.Px() << sep << initialPartons[i].Mom.Py() << sep << initialPartons[i].Mom.Pz() << sep << initialPartons[i].Mom.E() 
          << sep << tempParticles[0].Pos.X() << sep << tempParticles[0].Pos.Y() << sep << tempParticles[0].Pos.Z() << sep << tempParticles[0].Pos.T()
          << sep << initialPartons[i].FLAVOR << endl;
        }
      }
      file.close();
      break;
    }
    case heavyQuarkShower:
    {
      string filename = filename_prefix + "_" + "unshoweredParticles";
      fstream file( filename.c_str(), ios::out | ios::trunc );
      string sep = "\t";

      file << "#event" << sep << "px" << sep << "py" << sep << "pz" << sep << "E" << sep << "x"
      << sep << "y" << sep << "z" << sep << "t" << endl;

      for ( int n = 0; n < nEventsToGenerate; n++ )
      {
        vector< Particle > tempParticles, initialPartons;
        getHQShowerEvent( tempParticles, initialPartons );

        sample_TXYZ_singleParticle( tempParticles[0] );
        for ( int i = 0; i < tempParticles.size(); i++ )
        {
          tempParticles[i].Pos.T() = tempParticles[0].Pos.T();
          tempParticles[i].Pos.X() = tempParticles[0].Pos.X();
          tempParticles[i].Pos.Y() = tempParticles[0].Pos.Y();
          tempParticles[i].Pos.Z() = tempParticles[0].Pos.Z();
          tempParticles[i].N_EVENT_pp = n;
          tempParticles[i].N_EVENT_AA = n;
          _particles.push_back( tempParticles[i] );
        }

        for( int i = 0; i < initialPartons.size(); i++ )
        {
          file << n << sep << initialPartons[i].Mom.Px() << sep << initialPartons[i].Mom.Py() << sep << initialPartons[i].Mom.Pz() << sep << initialPartons[i].Mom.E()
          << sep << tempParticles[0].Pos.X() << sep << tempParticles[0].Pos.Y() << sep << tempParticles[0].Pos.Z() << sep << tempParticles[0].Pos.T()
          << sep << initialPartons[i].FLAVOR << endl;
        }
      }
      file.close();
      break;
    }
    default:
      throw ePYTHIAShower_error( "Unknow initial shower type. Unrecoverable error!" );
  }
}

vector< Particle > initialModel_PYTHIAShower::getFixedShowerEvent(const double _px, const double _py, const double _pz, const FLAVOR_TYPE _flavorA, const FLAVOR_TYPE _flavorB)
{

// PYTHIA seed has max value 9E8
  while( seed > 900000000 )
    seed -= 900000000;

  vector<Particle> particlesToAdd;
  int flavor1, flavor2;

  switch( _flavorA )
  {
  case down:
    flavor1 = 1;
    break;
  case anti_down:
    flavor1 = -1;
    break;
  case up:
    flavor1 = 2;
    break;
  case anti_up:
    flavor1 = -2;
    break;
  case strange:
    flavor1 = 3;
    break;
  case anti_strange:
    flavor1 = -3;
    break;
  case gluon:
    flavor1 = 21;
    break;
  default:
    cout << "Unknown flavor type...." << endl;
  }

  switch( _flavorB )
  {
  case down:
    flavor2 = 1;
    break;
  case anti_down:
    flavor2 = -1;
    break;
  case up:
    flavor2 = 2;
    break;
  case anti_up:
    flavor2 = -2;
    break;
  case strange:
    flavor2 = 3;
    break;
  case anti_strange:
    flavor2 = -3;
    break;
  case gluon:
    flavor2 = 21;
    break;
  default:
    cout << "Unknown flavor type...." << endl;
  }

  int attempt = 0;
  do
  {
    particlesToAdd.clear();

    fixedshowerevent_( &_px, &_py, &_pz, &flavor1, &flavor2, &seed );

    int index = 0;
    while( bamps_.pa[index][0] != 0 )
    {
      Particle tempParticle;
      tempParticle.FLAVOR = mapToPYTHIAflavor( bamps_.pa[index][0] );
      tempParticle.Mom.Px() = bamps_.pa[index][2];
      tempParticle.Mom.Py() = bamps_.pa[index][3];
      tempParticle.Mom.Pz() = bamps_.pa[index][4];
      tempParticle.Mom.E() = sqrt( tempParticle.Mom.Px() * tempParticle.Mom.Px() + tempParticle.Mom.Py() * tempParticle.Mom.Py() + tempParticle.Mom.Pz() * tempParticle.Mom.Pz() );
      tempParticle.m = 0.0;
      particlesToAdd.push_back( tempParticle );
      index++;
    }

    attempt++;
  }
  while( particlesToAdd.size() == 0 );

  if( attempt > 10 )
    cout << attempt << " Attempts needed to get an allowed shower." << endl;

  double sumE = 0.0;
  double E1 = sqrt( _px * _px + _py * _py + _pz * _pz );
  double E2 = sqrt( _px * _px + _py * _py + _pz * _pz );
  for( int i = 0; i < particlesToAdd.size(); i++ )
  {
    sumE += particlesToAdd[i].Mom.E();
  }
  if( FPT_COMP_GE( abs( sumE - ( E1 + E2 ) ) / sumE, 0.05 ) )
  {
    stringstream errMsg;
    errMsg << "Total energy of shower particles differs from energy of shower-initiating partons more than 5%:\t"
           << sumE << "\t" << E1 + E2 << "Unrecoverable error!" << endl;
    throw ePYTHIAShower_error( errMsg.str() );
  }

  return particlesToAdd;
}

void initialModel_PYTHIAShower::getPythiaShowerEvent( vector< Particle >& _particles, vector< Particle >& _initialPartons )
{
  _particles.clear();
  vector< Particle > tempParticles;  
  
  // PYTHIA seed has max value 9E8
  while( seed > 900000000 )
    seed -= 900000000;

  int attempt = 0;
  do
  {
    tempParticles.clear();
    
    pythiashowerevent_( &P0, &seed );

    int index = 0;
    while( bamps_.pa[index][0] != 0.0 )
    {
      Particle tempParticle;
      tempParticle.FLAVOR = mapToPYTHIAflavor( bamps_.pa[index][0] );
      tempParticle.Mom.Px() = bamps_.pa[index][2];
      tempParticle.Mom.Py() = bamps_.pa[index][3];
      tempParticle.Mom.Pz() = bamps_.pa[index][4];
      tempParticle.Mom.E() = sqrt( tempParticle.Mom.Px() * tempParticle.Mom.Px() + tempParticle.Mom.Py() * tempParticle.Mom.Py() + tempParticle.Mom.Pz() * tempParticle.Mom.Pz() );
      tempParticle.m = 0.0;
      tempParticles.push_back( tempParticle );
      index++;
    }
    attempt++;
  }
  while( tempParticles.size() == 0 );

  if( attempt > 10 )
    throw ePYTHIAShower_error( "Too many Attempts needed to get an allowed PYTHIA event." );

  _initialPartons.clear();
  _initialPartons.push_back( tempParticles[6] );
  _initialPartons.push_back( tempParticles[7] );
  
  for ( int i = 8; i < tempParticles.size(); i++ )
  {
    if ( tempParticles[i].FLAVOR != allFlavors )
      _particles.push_back( tempParticles[i] );
  }
}

void initialModel_PYTHIAShower::getPhotonShowerEvent( vector< Particle >& _particles, vector< Particle >& _initialPartons )
{
  _particles.clear();
  vector< Particle > tempParticles;  
  
  // PYTHIA seed has max value 9E8
  while( seed > 900000000 )
    seed -= 900000000;

  int attempt = 0;
  do
  {
    tempParticles.clear();
    
    photonshowerevent_( &P0, &seed );

    int index = 0;
    while( bamps_.pa[index][0] != 0.0 )
    {
      Particle tempParticle;
      tempParticle.FLAVOR = mapToPYTHIAflavor( bamps_.pa[index][0] );
      tempParticle.Mom.Px() = bamps_.pa[index][2];
      tempParticle.Mom.Py() = bamps_.pa[index][3];
      tempParticle.Mom.Pz() = bamps_.pa[index][4];
      tempParticle.Mom.E() = sqrt( tempParticle.Mom.Px() * tempParticle.Mom.Px() + tempParticle.Mom.Py() * tempParticle.Mom.Py() + tempParticle.Mom.Pz() * tempParticle.Mom.Pz() );
      tempParticle.m = 0.0;
      tempParticles.push_back( tempParticle );
      index++;
    }
    attempt++;
  }
  while( tempParticles.size() == 0 );

  if( attempt > 10 )
    throw ePYTHIAShower_error( "Too many Attempts needed to get an allowed PYTHIA event." );

  _initialPartons.clear();
  _initialPartons.push_back( tempParticles[7] );
  _initialPartons.push_back( tempParticles[8] );
  
  for ( int i = 9; i < tempParticles.size(); i++ )
  {
    if ( tempParticles[i].FLAVOR != allFlavors )
      _particles.push_back( tempParticles[i] );
  }
}

void initialModel_PYTHIAShower::getHQShowerEvent( vector< Particle >& _particles, vector< Particle >& _initialPartons )
{
  _particles.clear();
  vector< Particle > tempParticles;

  // PYTHIA seed has max value 9E8
  while( seed > 900000000 )
    seed -= 900000000;

  int attempt = 0;
  do
  {
    tempParticles.clear();

    double m_charm = Particle::Mcharm;
    double m_bottom = Particle::Mbottom;

    int jetFlavor;
    switch( initialPartonFlavor )
    {
      case charm:
        jetFlavor = 4;
        break;
      case bottom:
        jetFlavor = 5;
        break;
      default:
        string errMsg = "Heavy quark shower selected but no valid initial parton flavor given.";
        throw eInitialModel_error( errMsg );
    }
    hqshowerevent_( &P0, &seed, &m_charm, &m_bottom, &jetFlavor );

    int index = 0;
    while( bamps_.pa[index][0] != 0.0 )
    {
      Particle tempParticle;
      tempParticle.FLAVOR = mapToPYTHIAflavor( bamps_.pa[index][0] );
      tempParticle.Mom.Px() = bamps_.pa[index][2];
      tempParticle.Mom.Py() = bamps_.pa[index][3];
      tempParticle.Mom.Pz() = bamps_.pa[index][4];
      tempParticle.m = bamps_.pa[index][6];

      if( tempParticle.FLAVOR <= 2 * Particle::max_N_light_flavor )
      {
        tempParticle.m = 0;
      }

      tempParticle.Mom.E() = sqrt( tempParticle.Mom.vec2() + pow( tempParticle.m, 2.0 ) );
      tempParticles.push_back( tempParticle );
      index++;
    }
    attempt++;
  }
  while( tempParticles.size() == 0 );

  if( attempt > 10 )
    throw ePYTHIAShower_error( "Too many Attempts needed to get an allowed PYTHIA event." );

  _initialPartons.clear();
  _initialPartons.push_back( tempParticles[6] );
  _initialPartons.push_back( tempParticles[7] );

  for ( int i = 8; i < tempParticles.size(); i++ )
  {
    if ( tempParticles[i].FLAVOR != allFlavors )
      _particles.push_back( tempParticles[i] );
  }
}
