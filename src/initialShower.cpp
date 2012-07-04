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
#include <sstream>
#include <list>

#include "configuration.h"
#include "initialShower.h"
#include "initialmodel_minijets.h"
#include "pdfinterface.h"
#include "random.h"
#include "particle.h"
#include "FPT_compare.h"

using std::cout;
using std::endl;
using namespace ns_casc;


extern "C" {
  void shower_(double *px, double *py, double *pz1, double *pz2, int *pythiaFlavor1, int *pythiaFlavor2, double *tauf, long int *seed);
  struct{
    double pa[100][6];
    } bamps_;
}

initialShower::initialShower( const config& _config, WoodSaxon& _WoodSaxonParameter, const double _minimumPT, const int _nToGenerate ) :
  nEventsToGenerate( _config.getNaddedEvents() ),
  insertionTime( _config.getInsertionTime() ),
  seed( _config.getSeed() )
{
}

initialShower::~initialShower()
{
}

void initialShower::showerParticles(vector< Particle > _particles)
{
  for ( int index = 0; index < _particles.size(); index+=2 )
  {
    _particles[index].N_EVENT_pp = _particles[index+1].N_EVENT_pp = static_cast<int>( index / 2 );
  }
  
  vector<Particle> tempParticles;
  for (int i=0; i < _particles.size(); i+=2)
  {
    vector<Particle> showerParticles;
    double px = _particles[i].PX;
    double py = _particles[i].PY;
    double pz1 = _particles[i].PZ;
    double pz2 = _particles[i+1].PZ;
    FLAVOR_TYPE flavor1 = _particles[i].FLAVOR;
    FLAVOR_TYPE flavor2 = _particles[i+1].FLAVOR;
    showerParticles = getShowerParticles(px,py,pz1,pz2,flavor1,flavor2);
    if (showerParticles.size() == 0)
      cout << "No showered particles..." << endl;
      
    for (int j=0; j < showerParticles.size(); j++)
    {
      showerParticles[j].N_EVENT_pp = _particles[i].N_EVENT_pp;
      showerParticles[j].X = _particles[i].X;
      showerParticles[j].Y = _particles[i].Y;
      showerParticles[j].Z = _particles[i].Z;
      showerParticles[j].T = _particles[i].T;
        
      showerParticles[j].init = true;
        
      tempParticles.push_back( showerParticles[j] );
    }
  }
  _particles.clear();
  _particles = tempParticles;
  cout << "#### " << _particles.size() << " showered particles added after showering" << endl;
}

vector<Particle> initialShower::getShowerParticles( const double _px, const double _py, const double _pz1, const double _pz2, const FLAVOR_TYPE _flavor1, const FLAVOR_TYPE _flavor2 )
{
  vector<Particle> particlesToAdd;
  int flavor1, flavor2;
  
  switch (_flavor1)
  {
    case down: flavor1 = 1; break;
    case anti_down: flavor1 = -1; break;
    case up: flavor1 = 2; break;
    case anti_up: flavor1 = -2; break;
    case strange: flavor1 = 3; break;
    case anti_strange: flavor1 = -3; break;
    case gluon:  flavor1 = 21; break;
    default: cout << "Unknown flavor type...." << endl;
  }
  switch (_flavor2)
  {
    case down: flavor2 = 1; break;
    case anti_down: flavor2 = -1; break;
    case up: flavor2 = 2; break;
    case anti_up: flavor2 = -2; break;
    case strange: flavor2 = 3; break;
    case anti_strange: flavor2 = -3; break;
    case gluon:  flavor2 = 21; break;
    default: cout << "Unknown flavor type...." << endl;
  }

  int attempt = 0;
  do
  {
    particlesToAdd = createShower( flavor1, flavor2, _px, _py, _pz1, _pz2 );
    attempt++;
  } while (particlesToAdd.size() == 0);

  if (attempt > 10)
    cout << attempt << " Attempts needed to get an allowed shower." << endl;
  
  
  double sumE = 0.0;
  double E1 = sqrt( _px*_px + _py*_py + _pz1*_pz1 );
  double E2 = sqrt( _px*_px + _py*_py + _pz2*_pz2 );
  for (int i = 0; i < particlesToAdd.size(); i++)
  {
     sumE += particlesToAdd[i].E;
  }
  if (FPT_COMP_GE(abs(sumE-(E1+E2))/sumE,0.01))
    cout << "Total energy of shower particles is unequal energy of shower-initiating partons:\t" 
         << sumE << "\t" << E1 + E2 << endl;
  
  return particlesToAdd;
}


vector< Particle > initialShower::createShower(int flavor1, int flavor2, double px, double py, double pz1, double pz2)
{
  vector<Particle> showerParticles;
  
  shower_(&px,&py,&pz1,&pz2,&flavor1,&flavor2,&insertionTime,&seed);
  
  int index = 0;
  while (bamps_.pa[index][0] != 0)
  {
    ParticleOffline tempParticle;
    switch (static_cast<int>(bamps_.pa[index][0]))
    {
      case 21: tempParticle.FLAVOR = gluon; break;
      case 2: tempParticle.FLAVOR = up; break;
      case 1: tempParticle.FLAVOR = down; break;
      case -2: tempParticle.FLAVOR = anti_up; break;
      case -1: tempParticle.FLAVOR = anti_down; break;
      case 3: tempParticle.FLAVOR = strange; break;
      case -3: tempParticle.FLAVOR = anti_strange; break;
      default: cout << "Unknown flavor type:\t" << bamps_.pa[index][0] << endl;
    }
    
    tempParticle.PX = bamps_.pa[index][1];
    tempParticle.PY = bamps_.pa[index][2];
    tempParticle.PZ = bamps_.pa[index][3];
    tempParticle.E = sqrt( tempParticle.PX * tempParticle.PX + tempParticle.PY * tempParticle.PY + tempParticle.PZ * tempParticle.PZ);
    tempParticle.m = 0.0;
    showerParticles.push_back( tempParticle );
    index++;
  }

  cout << "Created pythia shower with cut-off time tau_i = " << insertionTime;
  cout << "PYTHIA seed:\t" << seed << endl;
  
  return showerParticles;
}



// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;
