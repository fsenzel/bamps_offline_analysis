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
  void pythiashowerevent_( const int* shower_type, const int* initialPartonFlavor, const double* sqrtS_perNN, const double* ptMin, const uint32_t* seed );
  struct
  {
    double pa[100][7];
  } bamps_;
}


initialModel_PYTHIAShower::initialModel_PYTHIAShower( const config& _config, WoodSaxon& _WoodSaxonParameter, const SHOWER_TYPE _shower_type, const double _ptCutOff ):
  initialModelWS( _config.getA(), _config.getAatomic(), _config.getB(), _config.getBatomic() ),
  nEventsToGenerate( _config.getNaddedEvents() ),
  seed( _config.getSeed() ),
  filename_prefix( _config.getStandardOutputDirectoryName() + "/" + _config.getJobName() ),
  shower_type( _shower_type ),
  ptCutOff( _ptCutOff ),
  initial_parton_flavor( _config.getInitialPartonFlavor() )
{
  impactParameter = _config.getImpactParameter();
  sqrtS_perNN = _config.getSqrtS();
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


initialModel_PYTHIAShower::initialModel_PYTHIAShower( WoodSaxon& _WoodSaxonParameter, const SHOWER_TYPE _shower_type,
                                                      const int _Nevents, const uint32_t _seed, const string _filename_prefix,
                                                      const double _A, const double _Aatomic, const double _B, const double _Batomic,
                                                      const double _sqrtS_perNN, const double _impactParameter,
                                                      const double _ptCutOff, const FLAVOR_TYPE _initialPartonFlavor ):
        initialModelWS( _A, _Aatomic, _B, _Batomic ),
        nEventsToGenerate( _Nevents ),
        seed( _seed ),
        filename_prefix( _filename_prefix ),
        shower_type( _shower_type ),
        ptCutOff( _ptCutOff ),
        initial_parton_flavor( _initialPartonFlavor )
{
  impactParameter = _impactParameter;
  sqrtS_perNN = _sqrtS_perNN;
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


void initialModel_PYTHIAShower::populateParticleVector( std::vector<Particle>& _particles )
{
  vector<Particle> dummy_initial_partons;
  populateParticleVector( _particles, dummy_initial_partons );
}

void initialModel_PYTHIAShower::populateParticleVector( std::vector< Particle >& _particles, std::vector< Particle >& _initial_partons )
{
  string filename = filename_prefix + "_" + "unshoweredParticles";
  fstream file( filename.c_str(), ios::out | ios::trunc );
  string sep = "\t";

  file << "#event" << sep << "px" << sep << "py" << sep << "pz" << sep << "E" << sep << "x"
       << sep << "y" << sep << "z" << sep << "t" << endl;

  for ( int n = 0; n < nEventsToGenerate; n++ )
  {
    vector< Particle > tempParticles, initialPartonPair;
    initialPartonPair.resize( 2 );

    if( shower_type == fixedShower || shower_type == fixedParton )
    {
      initialPartonPair[0].Mom = VectorEPxPyPz( ptCutOff, ptCutOff, 0.0, 0.0 );
      initialPartonPair[1].Mom = VectorEPxPyPz( ptCutOff, -ptCutOff, 0.0, 0.0 );

      switch ( initial_parton_flavor )
      {
        case 0:
          initialPartonPair[0].FLAVOR = gluon;
          initialPartonPair[1].FLAVOR = gluon;
          break;
        case 1:
          initialPartonPair[0].FLAVOR = up;
          initialPartonPair[1].FLAVOR = anti_up;
          break;
        case 2:
          initialPartonPair[0].FLAVOR = gluon;
          initialPartonPair[1].FLAVOR = up;
          break;
        default:
          throw ePYTHIAShower_error( "Unknown initial parton flavor. Unrecoverable error!" );
      }
    }

    if( shower_type != fixedParton )
    {
      getShowerEvent( tempParticles, initialPartonPair );
    }
    else if( shower_type == fixedParton )
    {
      tempParticles = initialPartonPair;
    }
    else
    {
      string errMsg = "Unknown shower type. Unrecoverable error!";
      throw ePYTHIAShower_error( errMsg );
    }

    sample_TXYZ_singleParticle( tempParticles[0] );
    for ( int i = 0; i < tempParticles.size(); i++ )
    {
      tempParticles[i].Pos = tempParticles[0].Pos;
      tempParticles[i].N_EVENT_pp = n;
      tempParticles[i].N_EVENT_AA = n;
      tempParticles[i].init = true;

      _particles.push_back( tempParticles[i] );
    }

    for( int i = 0; i < initialPartonPair.size(); i++ )
    {
      initialPartonPair[i].Pos = tempParticles[0].Pos;
      initialPartonPair[i].N_EVENT_pp = n;
      initialPartonPair[i].N_EVENT_AA = n;

      file << n << sep << initialPartonPair[i].Mom.Px() << sep << initialPartonPair[i].Mom.Py() << sep << initialPartonPair[i].Mom.Pz() << sep << initialPartonPair[i].Mom.E()
           << sep << initialPartonPair[i].Pos.X() << sep << initialPartonPair[i].Pos.Y() << sep << initialPartonPair[i].Pos.Z() << sep << initialPartonPair[i].Pos.T()
           << sep << initialPartonPair[i].FLAVOR << endl;

      _initial_partons.push_back( initialPartonPair[i] );
    }
  }
  file.close();
}

void initialModel_PYTHIAShower::getShowerEvent( vector< Particle >& _particles, vector< Particle >& _initial_parton_pair )
{
  _particles.clear();
  vector< Particle > tempParticles;

  // PYTHIA seed has max value 9E8
  while( seed > 900000000 )
    seed -= 900000000;

  int index_firstparton = 6;

//  set heavy quark masses for HQ shower
  if( ( initial_parton_flavor == charm || initial_parton_flavor == anti_charm ||
    initial_parton_flavor == bottom || initial_parton_flavor == anti_bottom ||
    shower_type == charmQuarkShower || shower_type == bottomQuarkShower ) && ( Particle::Mcharm != 1.3 || Particle::Mbottom != 4.6 ) )
  {
    string errMsg = "Heavy quark masses are not set to standard values (Mcharm=1.3 GeV, Mbottom=4.6 GeV) in shower routines. Unrecoverable error!";
    throw ePYTHIAShower_error( errMsg );
  }
  int shower_type_int = static_cast<int>( shower_type );
  int initial_parton_flavor_int = mapToBAMPSflavor( initial_parton_flavor );

  if( shower_type == fixedShower )
  {
    int flavorA = mapToBAMPSflavor( _initial_parton_pair[0].FLAVOR );
    int flavorB = mapToBAMPSflavor( _initial_parton_pair[1].FLAVOR );
    double px = _initial_parton_pair[0].Mom.Px();
    double py = _initial_parton_pair[0].Mom.Py();
    double pz = _initial_parton_pair[0].Mom.Pz();
    fixedshowerevent_( &px, &py, &pz, &flavorA, &flavorB, &seed );
    index_firstparton = -2;
  }
  else
  {
    pythiashowerevent_( &shower_type_int, &initial_parton_flavor_int, &sqrtS_perNN, &ptCutOff, &seed );
  }

  int index = 0;
  while( bamps_.pa[index][0] != 0.0 )
  {
    Particle tempParticle;
    tempParticle.FLAVOR = mapToPYTHIAflavor( bamps_.pa[index][0] );
    tempParticle.m = Particle::getMass( tempParticle.FLAVOR );
    tempParticle.Mom.Px() = bamps_.pa[index][2];
    tempParticle.Mom.Py() = bamps_.pa[index][3];
    tempParticle.Mom.Pz() = bamps_.pa[index][4];
    tempParticle.Mom.E() = sqrt( tempParticle.Mom.vec2() + pow( tempParticle.m, 2.0 ) );

    tempParticles.push_back( tempParticle );
    index++;
  }

  if( shower_type == fixedShower )
  {
    double sumE = 0.0;
    for( int i = 0; i < tempParticles.size(); i++ )
    {
      sumE += tempParticles[i].Mom.E();
    }
    if( FPT_COMP_GE( abs( sumE - ( _initial_parton_pair[0].Mom.E() + _initial_parton_pair[1].Mom.E() ) ) / sumE, 0.05 ) )
    {
      stringstream errMsg;
      errMsg << "Total energy of shower particles differs from energy of shower-initiating partons more than 5%:\t"
             << sumE << "\t" << _initial_parton_pair[0].Mom.E() + _initial_parton_pair[1].Mom.E() << "Unrecoverable error!" << endl;
      throw ePYTHIAShower_error( errMsg.str() );
    }
  }
  else
  {
    _initial_parton_pair.clear();
    _initial_parton_pair.push_back( tempParticles[index_firstparton] );
    _initial_parton_pair.push_back( tempParticles[index_firstparton + 1] );
  }

  for ( int i = index_firstparton + 2; i < tempParticles.size(); i++ )
  {
    if( tempParticles[i].FLAVOR != allFlavors && tempParticles[i].FLAVOR != photon )
      _particles.push_back( tempParticles[i] );
  }
}