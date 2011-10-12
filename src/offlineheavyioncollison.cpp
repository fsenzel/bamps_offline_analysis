#include <vector>
#include <list>
#include <fstream>
#include <iostream>
#include <deque>
#include <math.h>
#include <limits.h>

#include "offlineheavyioncollison.h"
#include "additionalparticlesdistribution.h"
#include "configuration.h"
#include "coordinateBins.h"
#include "cellcontainer.h"
#include "scattering22.h"
#include "scattering23.h"
#include "prefactors23.h"
#include "random.h"
#include "binary_cross_sections.h"
#include "offlineoutput.h"

using namespace ns_casc;
using namespace std;


extern int number, numberAdded, IX, IY, IZ;


namespace
{
  const int cellcut = 4;
  double dx, dy, dv, transLen;
  double timenow, timenext;
  double simulationTime;

  vector<cellContainer> cells;
  vector<cellContainer> cellsAdded;
  vector<cellContainer> cellsCopy;
  vector<cellContainer> cellsAddedCopy;
  list<int> edgeCell;
  list<int> formGeom;
  list<int> edgeCellAdded;
  list<int> formGeomAdded;

  list<int> deadParticleList;
  coordinateEtaBins etaBins;

  int gG, gQ;

  long ncoll, ncoll32, ncoll23, ncoll22, ncolle;

  double randomShiftEta, randomShiftX, randomShiftY;
  bool dodo1;
  int nn_ana1;
  
  double p23_collected_gluon;
  int n23_collected_gluon;
  double p23_collected_quark;
  int n23_collected_quark;
  double lambdaJet_gluon;
  double lambdaJet_quark;
}


offlineHeavyIonCollision::offlineHeavyIonCollision( config* const _config, offlineOutputInterface* const _offlineInterface ) :
    theConfig( _config ), stoptime_last( 0 ), stoptime( 5.0 ), currentNumber( 0 ),
    rings( _config->getRingNumber(), _config->getCentralRingRadius(), _config->getDeltaR() ),
    testpartcl( _config->getTestparticles() ),
    offlineInterface( _offlineInterface )
{
}



offlineHeavyIonCollision::~offlineHeavyIonCollision()
{

}



void offlineHeavyIonCollision::init()
{
  additionalParticlesDistribution addedStuff( theConfig, miniJetsInitialState );
  addedStuff.populateParticleVector( addedParticles, WoodSaxonParameter );

  int Nbefore = addedParticles.size();
  for ( int j = 0; j < addedParticles.size(); j++ )
  {
    if ( addedParticles[j].FLAVOR > ( Particle::N_light_flavor * 2 ) )
    {
      while ( addedParticles.back().FLAVOR > ( Particle::N_light_flavor * 2 ) )
      {
        addedParticles.pop_back();
      }
      addedParticles[j] = addedParticles.back();
      addedParticles.pop_back();
    }
  }
  cout << "#### " << addedParticles.size() << " out of " << Nbefore << " added ( N_f = " << Particle::N_light_flavor << " )." << endl;
  numberAdded = addedParticles.size();

  addedStuff.prepareParticles( addedParticles );
}



void offlineHeavyIonCollision::mainFramework( analysis& aa )
{
  double dt_cascade = 0;
  double nexttime;
  int ncoll_backup = 0;
  int ncoll22_backup = 0;
  int ncoll23_backup = 0;
  int ncoll32_backup = 0;
  int ncolle_backup = 0;
  
  list<int> edgeCellCopy, edgeCellAddedCopy;

  gG = 2 * ( pow( Ncolor, 2 ) - 1 );
  gQ = 2 * Ncolor;

  dx = theConfig->get_dx();
  dy = theConfig->get_dy();
  transLen = theConfig->getTransLen();
  
  edgeCellAdded.clear();
  edgeCell.clear();

  bool endOfDataFiles = false;
  bool doAnalysisStep = false;
  bool doMovieStep = false;
  bool doMovieStepMedium = false;
  bool again = false;
  int ncell = IX * IY * IZ;
  cells.resize( ncell );
  cellsCopy.resize( ncell );
  cellsAdded.resize( ncell );
  cellsAddedCopy.resize( ncell );
  for ( int i = 0; i < ncell; i++ )
  {
    cells[i].setCoordinates( i, dx, IX, transLen, dy, IY, transLen );
    cellsAdded[i].setCoordinates( i, dx, IX, transLen, dy, IY, transLen );
    cells[i].rates.normalizeRates();
    cellsCopy[i].rates.normalizeRates();
    cellsAdded[i].rates.normalizeRates();
    cellsAddedCopy[i].rates.normalizeRates();
  }

  vector<double> tempVec( rings.size(), 0 );
  rateGluons.assign( IZ, tempVec );
  rateQuarks.assign( IZ, tempVec );
  rateAntiQuarks.assign( IZ, tempVec );
  rateGluons_prev.assign( IZ, tempVec );
  rateQuarks_prev.assign( IZ, tempVec );
  rateAntiQuarks_prev.assign( IZ, tempVec );

  stoptime = theConfig -> getRuntime();
  int nn_ana = 0;
  int nn_ana_movie = 0;
  int jumpMovieSteps = 0;
  double factor_dt = theConfig->getFactor_dt();
  cout << "scale time steps dt by facot " << factor_dt << endl;

  aa.initialOutput();
  if ( theConfig->movieOutputJets )
  {
    aa.movieOutput( 0, jumpMovieSteps );
  }
  if ( theConfig->movieOutputBackground )
  {
    aa.movieOutputMedium( 0, jumpMovieSteps );
  }
  aa.collectPtDataInitial();
  aa.collectYDataInitial();
  aa.collectEtDataInitial();
  simulationTime = theConfig->getTimefirst(); //fm/c

  while ( simulationTime >= aa.tstep[nn_ana] )
  {
    nn_ana++;
  }
  while ( simulationTime >= aa.tstep_movie[nn_ana_movie] )
  {
    nn_ana_movie++;
    jumpMovieSteps++;
  }
  
  // propagate added particles to current time
  double cc;
  for ( int i = 0;i < addedParticles.size(); i++ )
  {
    if ( addedParticles[i].T <= simulationTime )
    {
      cc = ( simulationTime - addedParticles[i].T ) / addedParticles[i].E;
      addedParticles[i].T = simulationTime;
      addedParticles[i].X = addedParticles[i].X + addedParticles[i].PX * cc;
      addedParticles[i].Y = addedParticles[i].Y + addedParticles[i].PY * cc;
      addedParticles[i].Z = addedParticles[i].Z + addedParticles[i].PZ * cc;

      addedParticles[i].X_traveled += sqrt( pow( addedParticles[i].PX * cc, 2.0 ) + pow( addedParticles[i].PY * cc, 2.0 ) + pow( addedParticles[i].PZ * cc, 2.0 ) );
      aa.addJetEvent_initial( i );
    }
  }

  int n_dt = 0;
  double dt_sum = 0.0;
  int n_again = 0;

  int n_up = 0;
  int n_down = 0;
  int n_strange = 0;
  int n_anti_up = 0;
  int n_anti_down = 0;
  int n_anti_strange = 0;
  int n_gluon = 0;
  int n_quarks_temp = 0, n_anti_quarks_temp = 0;
  

  do
  {
    n_up = 0;
    n_down = 0;
    n_strange = 0;
    n_anti_up = 0;
    n_anti_down = 0;
    n_anti_strange = 0;
    n_gluon = 0;
    
    // remove init tag from particles after their formation
    for ( int i = 0; i < addedParticles.size(); i++ )
    {
      if ( addedParticles[i].T <= simulationTime )
      {
        addedParticles[i].init = false;
      }
      
      double pt_temp = sqrt( pow( addedParticles[i].PX, 2 ) + pow( addedParticles[i].PY, 2 ) );
      if ( pt_temp > theConfig->getMinimumPT() )
      {
        switch( addedParticles[i].FLAVOR )
        {
          case gluon:
            ++n_gluon;
            break;
          case up:
            ++n_up;
            break;
          case down:
            ++n_down;
            break;
          case strange:
            ++n_strange;
            break;
          case anti_up:
            ++n_anti_up;
            break;
          case anti_down:
            ++n_anti_down;
            break;            
          case anti_strange:
            ++n_anti_strange;
            break;
          default:
            cout << "##quarks error" << endl;
            break;
        }
      } 
    }
    n_quarks_temp = n_up + n_down + n_strange;
    n_anti_quarks_temp = n_anti_up + n_anti_down + n_anti_strange;

    // evolution of the medium to present time
    double dt_cascade_from_data = evolveMedium( simulationTime, endOfDataFiles );
    
    // specify time step
    if ( theConfig->DtSpecified() )
    {
      dt = theConfig->getDt();
    }
    else if ( dt_cascade_from_data >= 0 )
    {
      dt_cascade = dt_cascade_from_data;
      dt = dt_cascade * factor_dt;  // use time step from original medium evolution but make it slightly smaller
    }
    
    if ( endOfDataFiles )
    {
      cout << "* End of data files reached. Aborting." << endl;
      cout << "# time = " << simulationTime << "    dt = " << dt << endl;
      break;
    }
    cout << "# time = " << simulationTime << "    dt = " << dt << endl;

    p23_collected_gluon = 0;
    n23_collected_gluon = 0;
    p23_collected_quark = 0;
    n23_collected_quark = 0;
    lambdaJet_gluon = 0;
    lambdaJet_quark = 0;
      
    if ( doMovieStepMedium && theConfig->movieOutputBackground )
    {
      aa.movieOutputMedium( nn_ana_movie - 1, jumpMovieSteps );
      doMovieStepMedium = false;
    }

    //------------ make copies ---------------
    particles_atTimeNowCopy = particles_atTimeNow;
    addedParticlesCopy = addedParticles;
    cellsCopy = cells;
    cellsAddedCopy = cellsAdded;
    edgeCellCopy = edgeCell;
    edgeCellAddedCopy = edgeCellAdded;
    //--------------------------------------

    ncoll_backup = ncoll;
    ncoll22_backup = ncoll22;
    ncoll23_backup = ncoll23;
    ncoll32_backup = ncoll32;
    ncolle_backup = ncolle;
    //--------------------------

    nexttime = simulationTime + dt;

    // ask if it is time for analysis
    if ( nexttime >= aa.tstep[nn_ana] )
    {
      nexttime = aa.tstep[nn_ana];
      dt = nexttime - simulationTime;
      doAnalysisStep = true;
      cout << "profile " << nexttime << endl;
    }
    
    // ask if it is time for movie output
    if ( nexttime >= aa.tstep_movie[nn_ana_movie] )
    {
      nexttime = aa.tstep_movie[nn_ana_movie];
      dt = nexttime - simulationTime;
      doMovieStep = true;
      cout << "** movie: " << nexttime << endl;
    }
    
    cell_ID( nexttime );

    // collide added particles with gluonic medium
    deadParticleList.clear();
    scattering( nexttime, again, aa );

    // if time step is too large -> collide again with smaller time step
    while ( again )
    {
      p23_collected_gluon = 0;
      n23_collected_gluon = 0;
      p23_collected_quark = 0;
      n23_collected_quark = 0;
      lambdaJet_gluon = 0;
      lambdaJet_quark = 0;
      
      doAnalysisStep = false;
      doMovieStep = false;
      n_again++;

      dt_cascade = dt / factor_dt; // use same dt for next time steps until cascade data defines new dt.

      particles_atTimeNow = particles_atTimeNowCopy;
      addedParticles = addedParticlesCopy;
      numberAdded = addedParticles.size();
      cells = cellsCopy;
      cellsAdded = cellsAddedCopy;
      edgeCell = edgeCellCopy;
      edgeCellAdded = edgeCellAddedCopy;
      
      ncoll = ncoll_backup;
      ncoll22 = ncoll22_backup;
      ncoll23 = ncoll23_backup;
      ncoll32 = ncoll32_backup;
      ncolle = ncolle_backup;

      nexttime = simulationTime + dt;

      cell_ID( nexttime );

      deadParticleList.clear();
      scattering( nexttime, again, aa );
    }

    scatterEdgeParticles( edgeCell, edgeCellAdded, nexttime );

    removeDeadParticles( aa );
    
    if ( doAnalysisStep )
    {
      aa.intermediateOutput( nn_ana );
      aa.collectPtData( nn_ana );
      aa.collectYData( nn_ana );
      aa.collectEtData( nn_ana );
      nn_ana++;
      doAnalysisStep = false;
    }
    
    if ( doMovieStep )
    {
      aa.mfpJetsOutput( nn_ana_movie, jumpMovieSteps );
      
      if ( theConfig->movieOutputJets )
      {
        aa.movieOutput( nn_ana_movie, jumpMovieSteps );
      }
      nn_ana_movie++;
      doMovieStep = false;
      doMovieStepMedium = true;
    }
    aa.printCentralDensities( simulationTime );

    // analyse timesteps
    dt_sum += dt;
    n_dt++;

    simulationTime = nexttime;
  }
  while ( simulationTime < stoptime && !endOfDataFiles );//fm/c

  aa.finalOutput( stoptime );
  aa.addJetEvents_final();

  cout << "number of jet collisions: " << ncoll << endl;

  cout << "average dt = " << dt_sum / n_dt << endl;
  cout << "last dt = " << dt << endl;
  cout << "number of agains: = " << n_again << endl;
}




double offlineHeavyIonCollision::evolveMedium( const double evolveToTime, bool& _endOfDataFiles )
{
  double dt_cascade = -1;
  bool stop = false;
  int iscat, jscat, kscat, dead;
  double time = 0;
  double timei = 0;
  double timej = 0;
  double pxi, pyi, pzi, pxj, pyj, pzj;
  FLAVOR_TYPE F1, F2, F3;
  double c;

  // give partcl from cascade the values for cell structurement which could have changed in collisions()
  for ( int k = 0; k < particles.size(); k++ )
  {
    particles[k].init = particles_atTimeNow[k].init;
    particles[k].edge = particles_atTimeNow[k].edge;
    particles[k].free = particles_atTimeNow[k].free;
  }

  if ( evolveToTime <= stoptime_last )
  {
    for ( int i = 0; i < theConfig->getN_init(); i++ )
    {
      particles[i].init = true;

      particles[i].T = particles_init[i].T;
      particles[i].X = particles_init[i].X;
      particles[i].Y = particles_init[i].Y;
      particles[i].Z = particles_init[i].Z;

      particles[i].PXold = particles[i].PX = particles_init[i].PX;
      particles[i].PYold = particles[i].PY = particles_init[i].PY;
      particles[i].PZold = particles[i].PZ = particles_init[i].PZ;
      particles[i].Eold = particles[i].E = particles_init[i].E;
    }

    number = theConfig->getN_init();
  }

  offlineEventType actiontype = event_dummy;

  while (( actiontype != event_endOfCascade ) && ( !stop ) )
  {
    try
    {
      boost::shared_ptr< offlineDataEventType > ptrEventType = offlineInterface->readOfflineDataFromArchive< offlineDataEventType >();
      actiontype = ptrEventType->event;
    }    
    catch ( boost::archive::archive_exception& err )
    {
      stop = true;
      _endOfDataFiles = true;
    }
    
    if ( actiontype == event_newTimestep )
    {
      boost::shared_ptr< offlineDataCellConfiguration > ptrCellStructure = offlineInterface->readOfflineDataFromArchive< offlineDataCellConfiguration >();
      etaBins = ptrCellStructure->etaBins;
      timenow = ptrCellStructure->timenow;
      timenext = ptrCellStructure->timenext;
      randomShiftX = ptrCellStructure->randomShiftX;
      randomShiftY = ptrCellStructure->randomShiftY;
      randomShiftEta = ptrCellStructure->randomShiftEta;
         
      dt_cascade = timenext - timenow;

      rateGluons_prev = rateGluons;
      rateQuarks_prev = rateQuarks;
      rateAntiQuarks_prev = rateAntiQuarks;
      for ( int i = 0; i < IZ; i++ )
      {
        rateGluons[i].assign( rings.size(), 0 );
        rateQuarks[i].assign( rings.size(), 0 );
        rateAntiQuarks[i].assign( rings.size(), 0 );
      }

      boost::shared_ptr< offlineDataInteractionRates > ptrRates = offlineInterface->readOfflineDataFromArchive< offlineDataInteractionRates >();
      rateGluons = ptrRates->gluonRates;
      rateQuarks = ptrRates->quarkRates;
      rateAntiQuarks = ptrRates->antiQuarkRates;
    }
    else if ( actiontype == event_interaction22 )
    {
      boost::shared_ptr< offlineDataInteraction22 > ptrInteraction22 = offlineInterface->readOfflineDataFromArchive< offlineDataInteraction22 >();
      time = ptrInteraction22->time;
      iscat = ptrInteraction22->iscat;
      jscat = ptrInteraction22->jscat;
      pxi = ptrInteraction22->pix;
      pyi = ptrInteraction22->piy;
      pzi = ptrInteraction22->piz;
      pxj = ptrInteraction22->pjx;
      pyj = ptrInteraction22->pjy;
      pzj = ptrInteraction22->pjz;
      F1 = ptrInteraction22->F1;
      F2 = ptrInteraction22->F2;
      
      if ( time <= ( evolveToTime + 1.0e-6 ) )
      {
        particles[iscat].init = particles[jscat].init = false;

        c = ( time - particles[iscat].T ) / particles[iscat].E;
        particles[iscat].T = time;
        particles[iscat].X = particles[iscat].X + particles[iscat].PX * c;
        particles[iscat].Y = particles[iscat].Y + particles[iscat].PY * c;
        particles[iscat].Z = particles[iscat].Z + particles[iscat].PZ * c;

        c = ( time - particles[jscat].T ) / particles[jscat].E;
        particles[jscat].T = time;
        particles[jscat].X = particles[jscat].X + particles[jscat].PX * c;
        particles[jscat].Y = particles[jscat].Y + particles[jscat].PY * c;
        particles[jscat].Z = particles[jscat].Z + particles[jscat].PZ * c;

        if ( time < particles[iscat].T )
        {
          cout << "back22_i" << endl;
          int zz;
          cin >> zz;
        }
        if ( time < particles[jscat].T )
        {
          cout << "back22_j" << endl;
          int zz;
          cin >> zz;
        }

        particles[iscat].FLAVOR = static_cast<FLAVOR_TYPE>( F1 );
        particles[iscat].PX = pxi;
        particles[iscat].PY = pyi;
        particles[iscat].PZ = pzi;
        particles[iscat].E = sqrt( pow( pxi, 2 ) + pow( pyi, 2 ) + pow( pzi, 2 ) );

        particles[jscat].FLAVOR = static_cast<FLAVOR_TYPE>( F2 );
        particles[jscat].PX = pxj;
        particles[jscat].PY = pyj;
        particles[jscat].PZ = pzj;
        particles[jscat].E = sqrt( pow( pxj, 2 ) + pow( pyj, 2 ) + pow( pzj, 2 ) );
      }
      else
      {
        stop = true;
        boost::shared_ptr< offlineDataEventType > tempPtr( new offlineDataEventType(event_interaction22) );
        offlineInterface->temporaryStoreData< offlineDataEventType >( tempPtr );
        offlineInterface->temporaryStoreData< offlineDataInteraction22 >( ptrInteraction22 );
      }
    }
    else if ( actiontype == event_interaction23 )
    {
      boost::shared_ptr< offlineDataInteraction23 > ptrInteraction23 = offlineInterface->readOfflineDataFromArchive< offlineDataInteraction23 >();
      time = ptrInteraction23->time;
      iscat = ptrInteraction23->iscat;
      jscat = ptrInteraction23->jscat;
      kscat = ptrInteraction23->newp;
      pxi = ptrInteraction23->pix;
      pyi = ptrInteraction23->piy;
      pzi = ptrInteraction23->piz;
      pxj = ptrInteraction23->pjx;
      pyj = ptrInteraction23->pjy;
      pzj = ptrInteraction23->pjz;
      F1 = static_cast<FLAVOR_TYPE>( ptrInteraction23->F1 );
      F2 = static_cast<FLAVOR_TYPE>( ptrInteraction23->F2 );
      F3 = static_cast<FLAVOR_TYPE>( ptrInteraction23->F3 );
      particles[kscat].X = ptrInteraction23->newx;
      particles[kscat].Y = ptrInteraction23->newy;
      particles[kscat].Z = ptrInteraction23->newz;
      particles[kscat].PX = ptrInteraction23->newpx;
      particles[kscat].PY = ptrInteraction23->newpy;
      particles[kscat].PZ = ptrInteraction23->newpz;

      if ( time <= evolveToTime + 1.0e-6 )
      {
        particles[kscat].T = time;
        particles[kscat].FLAVOR = static_cast<FLAVOR_TYPE>( F3 );
        particles[kscat].E = sqrt( pow( particles[kscat].PX, 2 ) + pow( particles[kscat].PY, 2 ) + pow( particles[kscat].PZ, 2 ) );
        particles[kscat].free = false;

        particles[iscat].init = particles[jscat].init = particles[kscat].init = false;

        c = ( time - particles[iscat].T ) / particles[iscat].E;
        particles[iscat].T = time;
        particles[iscat].X = particles[iscat].X + particles[iscat].PX * c;
        particles[iscat].Y = particles[iscat].Y + particles[iscat].PY * c;
        particles[iscat].Z = particles[iscat].Z + particles[iscat].PZ * c;

        c = ( time - particles[jscat].T ) / particles[jscat].E;
        particles[jscat].T = time;
        particles[jscat].X = particles[jscat].X + particles[jscat].PX * c;
        particles[jscat].Y = particles[jscat].Y + particles[jscat].PY * c;
        particles[jscat].Z = particles[jscat].Z + particles[jscat].PZ * c;

        if ( time < particles[iscat].T )
        {
          cout << "back23_i" << endl;
          int zz;
          cin >> zz;
        }
        if ( time < particles[jscat].T )
        {
          cout << "back23_j" << endl;
          int zz;
          cin >> zz;
        }

        particles[iscat].FLAVOR = static_cast<FLAVOR_TYPE>( F1 );
        particles[iscat].PX = pxi;
        particles[iscat].PY = pyi;
        particles[iscat].PZ = pzi;
        particles[iscat].E = sqrt( pxi * pxi + pyi * pyi + pzi * pzi );

        particles[jscat].FLAVOR = static_cast<FLAVOR_TYPE>( F2 );
        particles[jscat].PX = pxj;
        particles[jscat].PY = pyj;
        particles[jscat].PZ = pzj;
        particles[jscat].E = sqrt( pxj * pxj + pyj * pyj + pzj * pzj );

        number++;//production
      }
      else
      {
        stop = true;
        boost::shared_ptr< offlineDataEventType > tempPtr( new offlineDataEventType(event_interaction23) );
        offlineInterface->temporaryStoreData< offlineDataEventType >( tempPtr );
        offlineInterface->temporaryStoreData< offlineDataInteraction23 >( ptrInteraction23 );
      }
    }
    else if ( actiontype == event_interaction32 )
    {
      boost::shared_ptr< offlineDataInteraction32 > ptrInteraction32 = offlineInterface->readOfflineDataFromArchive< offlineDataInteraction32 >();
      time = ptrInteraction32->time;
      iscat = ptrInteraction32->iscat;
      jscat = ptrInteraction32->jscat;
      dead = ptrInteraction32->dead;
      pxi = ptrInteraction32->pix;
      pyi = ptrInteraction32->piy;
      pzi = ptrInteraction32->piz;
      pxj = ptrInteraction32->pjx;
      pyj = ptrInteraction32->pjy;
      pzj = ptrInteraction32->pjz;
      F1 = ptrInteraction32->F1;
      F2 = ptrInteraction32->F2;

      if ( time <= evolveToTime + 1.0e-6 )
      {
        particles[iscat].init = particles[jscat].init = false;

        c = ( time - particles[iscat].T ) / particles[iscat].E;
        particles[iscat].T = time;
        particles[iscat].X = particles[iscat].X + particles[iscat].PX * c;
        particles[iscat].Y = particles[iscat].Y + particles[iscat].PY * c;
        particles[iscat].Z = particles[iscat].Z + particles[iscat].PZ * c;

        c = ( time - particles[jscat].T ) / particles[jscat].E;
        particles[jscat].T = time;
        particles[jscat].X = particles[jscat].X + particles[jscat].PX * c;
        particles[jscat].Y = particles[jscat].Y + particles[jscat].PY * c;
        particles[jscat].Z = particles[jscat].Z + particles[jscat].PZ * c;

        if ( time < particles[iscat].T )
        {
          cout << "back32_i" << endl;
          int zz;
          cin >> zz;
        }
        if ( time < particles[jscat].T )
        {
          cout << "back32_j" << endl;
          int zz;
          cin >> zz;
        }

        particles[iscat].FLAVOR = static_cast<FLAVOR_TYPE>( F1 );
        particles[iscat].PX = pxi;
        particles[iscat].PY = pyi;
        particles[iscat].PZ = pzi;
        particles[iscat].E = sqrt( pxi * pxi + pyi * pyi + pzi * pzi );

        particles[jscat].FLAVOR = static_cast<FLAVOR_TYPE>( F2 );
        particles[jscat].PX = pxj;
        particles[jscat].PY = pyj;
        particles[jscat].PZ = pzj;
        particles[jscat].E = sqrt( pxj * pxj + pyj * pyj + pzj * pzj );
        number--;
      }
      else
      {
        stop = true;
        boost::shared_ptr< offlineDataEventType > tempPtr( new offlineDataEventType(event_interaction32) );
        offlineInterface->temporaryStoreData< offlineDataEventType >( tempPtr );
        offlineInterface->temporaryStoreData< offlineDataInteraction32 >( ptrInteraction32 );
      }
    }
    else if ( actiontype == event_interactionElastic )
    {
      boost::shared_ptr< offlineDataInteractionElastic > ptrInteractionElastic = offlineInterface->readOfflineDataFromArchive< offlineDataInteractionElastic >();
      timei = ptrInteractionElastic->ct_i;
      timej = ptrInteractionElastic->ct_j;
      iscat = ptrInteractionElastic->iscat;
      jscat = ptrInteractionElastic->jscat;
      pxi = ptrInteractionElastic->pix;
      pyi = ptrInteractionElastic->piy;
      pzi = ptrInteractionElastic->piz;
      pxj = ptrInteractionElastic->pjx;
      pyj = ptrInteractionElastic->pjy;
      pzj = ptrInteractionElastic->pjz;
      
      if (( timei <= evolveToTime + 1.0e-6 ) || ( timej <= evolveToTime + 1.0e-6 ) )
      {
        particles[iscat].init = particles[jscat].init = false;

        particles[iscat].PXold = particles[iscat].PX;
        particles[iscat].PYold = particles[iscat].PY;
        particles[iscat].PZold = particles[iscat].PZ;
        particles[iscat].Eold = particles[iscat].E;

        particles[jscat].PXold = particles[jscat].PX;
        particles[jscat].PYold = particles[jscat].PY;
        particles[jscat].PZold = particles[jscat].PZ;
        particles[jscat].Eold = particles[jscat].E;

        c = ( timei - particles[iscat].T ) / particles[iscat].E;
        particles[iscat].T = timei;
        particles[iscat].X = particles[iscat].X + particles[iscat].PX * c;
        particles[iscat].Y = particles[iscat].Y + particles[iscat].PY * c;
        particles[iscat].Z = particles[iscat].Z + particles[iscat].PZ * c;

        particles[iscat].PX = pxi;
        particles[iscat].PY = pyi;
        particles[iscat].PZ = pzi;
        particles[iscat].E = sqrt( pxi * pxi + pyi * pyi + pzi * pzi );

        c = ( timej - particles[jscat].T ) / particles[jscat].E;
        particles[jscat].T = timej;
        particles[jscat].X = particles[jscat].X + particles[jscat].PX * c;
        particles[jscat].Y = particles[jscat].Y + particles[jscat].PY * c;
        particles[jscat].Z = particles[jscat].Z + particles[jscat].PZ * c;

        particles[jscat].PX = pxj;
        particles[jscat].PY = pyj;
        particles[jscat].PZ = pzj;
        particles[jscat].E = sqrt( pxj * pxj + pyj * pyj + pzj * pzj );

        if ( timei < particles[iscat].T )
        {
          cout << "back22e_i" << endl;
          int zz;
          cin >> zz;
        }
        if ( timej < particles[jscat].T )
        {
          cout << "back22e_j" << endl;
          int zz;
          cin >> zz;
        }
      }
      else
      {
        stop = true;
        boost::shared_ptr< offlineDataEventType > tempPtr( new offlineDataEventType(event_interactionElastic) );
        offlineInterface->temporaryStoreData< offlineDataEventType >( tempPtr );
        offlineInterface->temporaryStoreData< offlineDataInteractionElastic >( ptrInteractionElastic );
      }
    }
    else if ( actiontype == event_particleIdSwap )
    {
      boost::shared_ptr< offlineDataParticleIdSwap > ptrSwap = offlineInterface->readOfflineDataFromArchive< offlineDataParticleIdSwap >();
      iscat = ptrSwap->removedParticleID;
      jscat = ptrSwap->replacingParticleID;

      particles[iscat] = particles[jscat];
    }
  }

  for ( int i = 0; i < number; i++ )
  {
    if ( particles[i].T < evolveToTime + 1.0e-8 )
    {
      particles[i].init = false;
    }
  }

  // duplicate partcl from cascade to partclAtTimeNow
  particles_atTimeNow = particles;


  // propagate particles_atTimeNow to current time
  for ( int i = 0; i < number; i++ )
  {
    if ( particles_atTimeNow[i].T < evolveToTime )
    {
      c = ( evolveToTime - particles_atTimeNow[i].T ) / particles_atTimeNow[i].E;
      particles_atTimeNow[i].T = evolveToTime;
      particles_atTimeNow[i].X = particles_atTimeNow[i].X + particles_atTimeNow[i].PX * c;
      particles_atTimeNow[i].Y = particles_atTimeNow[i].Y + particles_atTimeNow[i].PY * c;
      particles_atTimeNow[i].Z = particles_atTimeNow[i].Z + particles_atTimeNow[i].PZ * c;
    }
  }

  stoptime_last = evolveToTime;
  return dt_cascade;
}




void offlineHeavyIonCollision::cell_ID( double _time )
{
  int nx, ny, nz, cell_id;
  double halfsize, eta;
  int kk, len, sum, sum1;

  halfsize = transLen / 2.0;
  int centralEtaIndex = static_cast<int>( etaBins.size() ) / 2;

  int NinAct, Nfree1, NinFormInit, NinFormGeom;
  NinAct = Nfree1 = NinFormInit = NinFormGeom = 0;

  formGeom.clear();
  for ( int i = 0; i < cells.size(); i++ )
  {
    cells[i].particleList.clear();
  }

  for ( int i = 0; i < number; i++ )
  {
    if ( particles_atTimeNow[i].T < _time )
    {
      NinAct++;
      particles_atTimeNow[i].init = false;

      if ( fabs( particles_atTimeNow[i].X - randomShiftX - halfsize ) < 1.0e-6 )
      {
        nx = IX - 1;
      }
      else
      {
        nx = static_cast<int>((( particles_atTimeNow[i].X - randomShiftX ) / transLen + 0.5 ) * IX );
      }

      if ( fabs( particles_atTimeNow[i].Y - randomShiftY - halfsize ) < 1.0e-6 )
      {
        ny = IY - 1;
      }
      else
      {
        ny = static_cast<int>((( particles_atTimeNow[i].Y - randomShiftY ) / transLen + 0.5 ) * IY );
      }

      eta = 0.5 * log(( particles_atTimeNow[i].T + particles_atTimeNow[i].Z ) / ( particles_atTimeNow[i].T - particles_atTimeNow[i].Z ) );
      nz = centralEtaIndex;
      if ( eta >= 0 )
      {
        while ( eta >= etaBins[nz].right && nz < IZ )
        {
          ++nz;
        }
      }
      else
      {
        while ( eta < etaBins[nz].left && nz >= 0 )
        {
          --nz;
        }
      }

      if (( nz >= IZ ) || ( nz < 0 ) )
      {
        std::string errMsg = "error in cell_id, nz out of range";
        throw eHIC_error( errMsg );
      }

      if (( nz < etaBins.min_index() ) || ( nz > etaBins.max_index() ) )
      {
        cell_id = -2;// -2:edge
        particles_atTimeNow[i].rate = 0.0;//GeV
      }
      else if (( nx < 0 ) || ( nx >= IX ) || ( ny < 0 ) || ( ny >= IY ) )
      {
        cell_id = -100;// -100:free
        Nfree1++;

        particles_atTimeNow[i].free = true;
        particles_atTimeNow[i].edge = false;

        particles_atTimeNow[i].rate = 0.0;//GeV
      }
      else
      {
        cell_id = nx + IX * ny + IX * IY * nz;
      }

      particles_atTimeNow[i].cell_id = cell_id;
      if ( cell_id >= 0 )
      {
        cells[cell_id].particleList.push_back( i );
      }
      else if ( cell_id == -2 )
      {
        if ( !particles_atTimeNow[i].edge ) // new member
        {
          edgeCell.push_back( i );
        }
      }
    }
    else
    {
      if ( particles_atTimeNow[i].init )
      {
        NinFormInit++;
      }
      else
      {
        NinFormGeom++;
        formGeom.push_back( i );
      }

      particles_atTimeNow[i].cell_id = -1;
      particles_atTimeNow[i].rate = 0.0;//GeV
    }
  }


  //------------------- added particles ---------------------------------
  for ( int i = 0; i < cells.size(); i++ )
  {
    cellsAdded[i].particleList.clear();
  }
  formGeomAdded.clear();

  for ( int i = 0; i < numberAdded; i++ )
  {
    if ( addedParticles[i].T < _time )
    {
      addedParticles[i].init = false;

      if ( fabs( addedParticles[i].X - randomShiftX - halfsize ) < 1.0e-6 )
      {
        nx = IX - 1;
      }
      else
      {
        nx = static_cast<int>((( addedParticles[i].X - randomShiftX ) / transLen + 0.5 ) * IX );
      }

      if ( fabs( addedParticles[i].Y - randomShiftY - halfsize ) < 1.0e-6 )
      {
        ny = IY - 1;
      }
      else
      {
        ny = static_cast<int>((( addedParticles[i].Y - randomShiftY ) / transLen + 0.5 ) * IY );
      }

      eta = 0.5 * log(( addedParticles[i].T + addedParticles[i].Z ) / ( addedParticles[i].T - addedParticles[i].Z ) );
      nz = centralEtaIndex;
      if ( eta >= 0 )
      {
        while ( eta >= etaBins[nz].right && nz < IZ )
        {
          ++nz;
        }
      }
      else
      {
        while ( eta < etaBins[nz].left && nz >= 0 )
        {
          --nz;
        }
      }

      if (( nz >= IZ ) || ( nz < 0 ) )
      {
        std::string errMsg = "error in cell_id, nz out of range for added particles";
        throw eHIC_error( errMsg );
      }

      if (( nz < etaBins.min_index() ) || ( nz > etaBins.max_index() ) )
      {
        cell_id = -2;// -2:edge
        addedParticles[i].rate = 0.0;//GeV
      }
      else if (( nx < 0 ) || ( nx >= IX ) || ( ny < 0 ) || ( ny >= IY ) )
      {
        cell_id = -100;// -100:free

        addedParticles[i].free = true;
        addedParticles[i].edge = false;

        addedParticles[i].rate = 0.0;//GeV
      }
      else
      {
        cell_id = nx + IX * ny + IX * IY * nz;
      }

      addedParticles[i].cell_id = cell_id;
      if ( cell_id >= 0 )
      {
        cellsAdded[cell_id].particleList.push_back( i );
      }
      else if ( cell_id == -2 )
      {
        if ( !addedParticles[i].edge ) // new member
        {
          edgeCellAdded.push_back( i );
          addedParticles[i].edge = true;
          addedParticles[i].collisionTime = infinity;
          addedParticles[i].collisionPartner = -1;
          addedParticles[i].md2g = -1.0;
          addedParticles[i].md2q = -1.0;
        }
      }
    }
    else
    {
      if ( addedParticles[i].init )
      {
//         nCharmInit++;
      }
      else
      {
//         nCharmOther++;
        formGeomAdded.push_back( i );
      }

      addedParticles[i].cell_id = -1;
      addedParticles[i].rate = 0.0;//GeV
    }
  }

}






void offlineHeavyIonCollision::scattering( const double nexttime, bool& again, analysis& aa )
{
  double xt;
  double ee;
  
  int nGluons = 0;
  int nAllQuarks = 0;
  int nAllAntiQuarks = 0;
  int nGluonsAdded = 0;
  int nAllQuarksAdded = 0;
  int nAllAntiQuarksAdded = 0;
  
  int IXY, id, nc;
  
  double beta_rf[4];
  double dz, XT, pr, dvv, cc, xx, yy, zz, eta, vrr;
  bool free;
  list<int> formGeomCopy;
  vector<int> gluonList, allParticlesList;
  vector<int> gluonListAdded, allParticlesListAdded;
  
  again = false;
  
  IXY = IX * IY;
  
  formGeomCopy = formGeom;
  
  gG = 2 * ( pow( Ncolor, 2 ) - 1 );
  
  if ( !formGeomCopy.empty() )
  {
    list<int>::iterator iIt;
    int id = -1;
    for ( iIt = formGeomCopy.begin(); iIt != formGeomCopy.end(); )
    {
      id = *iIt;
      cc = ( timenow - particles[id].T ) / sqrt( pow( particles[id].PXold, 2 ) + pow( particles[id].PYold, 2 ) + pow( particles[id].PZold, 2 ) );
      zz = particles[id].Z + particles[id].PZold * cc;
      eta = 0.5 * log(( timenow + zz ) / ( timenow - zz ) );
      
      if (( eta < etaBins[etaBins.min_index()].left ) || ( eta > etaBins[etaBins.max_index()].right ) )
      {
        iIt = formGeomCopy.erase( iIt );
      }
      else
      {
        ++iIt;
      }
    }
  }
  
  
  int centralEtaIndex = static_cast<int>( etaBins.size() ) / 2;
  
  // go through cells
  for ( int etaSliceIndex = etaBins.min_index(); etaSliceIndex <= etaBins.max_index(); etaSliceIndex++ )
  {
    //---------- populate the ring structure for averages ----------
    dz = timenow * ( tanh( etaBins[etaSliceIndex].right ) - tanh( etaBins[etaSliceIndex].left ) );
    
    rings.clear();
    rings.setLongitudinalGeometry( etaBins[etaSliceIndex].left, etaBins[etaSliceIndex].right, timenow );
    
    for ( int j = IXY * etaSliceIndex; j < IXY * ( etaSliceIndex + 1 ); j++ )
    {
      list<int>::const_iterator iIt;
      for ( iIt = cells[j].particleList.begin(); iIt != cells[j].particleList.end(); iIt++ )
      {
        xt = sqrt( pow( particles_atTimeNow [( *iIt )].X, 2 )  + pow( particles_atTimeNow[( *iIt )].Y, 2 ) );
        
        rings.addParticle( particles_atTimeNow[( *iIt )] );
        //           rings.addRates( particles_atTimeNow[(*iIt)] );
      }
    }
    
    list<int>::iterator iIt;
    for ( iIt = formGeomCopy.begin(); iIt != formGeomCopy.end(); iIt++ )
    {
      int id = *iIt;
      
      ee = sqrt( pow( particles_atTimeNow[id].PXold, 2 ) + pow( particles_atTimeNow[id].PYold, 2 ) + pow( particles_atTimeNow[id].PZold, 2 ) );
      cc = ( timenow - particles_atTimeNow[id].T ) / ee;
      zz = particles_atTimeNow[id].Z + particles_atTimeNow[id].PZold * cc;
      eta = 0.5 * log(( timenow + zz ) / ( timenow - zz ) );
      
      if (( eta >= etaBins[etaSliceIndex].left ) && ( eta <= etaBins[etaSliceIndex].right ) )
      {
        rings.addParticleInFormGeom( particles_atTimeNow[id], timenow );
        iIt = formGeomCopy.erase( iIt );    // erase this particle such that it needs not be looped over for the next rapidity slab
      }
    }
    
    rings.prepareAverages( dz, testpartcl );
    if ( etaSliceIndex == centralEtaIndex )
    {
      aa.centralRingsCopyFromCascade = rings;
    }
    //---------- populate the ring structure for averages ----------
    
    // scatterings in cell or not?
    dv = dx * dy * dz;
    for ( int j = IXY * etaSliceIndex; j < IXY * ( etaSliceIndex + 1 ); j++ )
    {
      if ( !( cells[j].empty() || cellsAdded[j].empty() ) )
      {
        cells[j].resetStoredValues();
        cellsAdded[j].resetStoredValues();
        gluonList.clear();
        allParticlesList.clear();
        gluonList.reserve( cells[j].size() );
        allParticlesList.reserve( cells[j].size() );
        gluonListAdded.clear();
        allParticlesListAdded.clear();
        gluonListAdded.reserve( cellsAdded[j].size() );
        allParticlesListAdded.reserve( cellsAdded[j].size() );
        
        //------------------------- check whether all particles in this cell are free --------------------
        free = true;
        for ( iIt = cells[j].particleList.begin(); iIt != cells[j].particleList.end(); iIt++ )
        {
          int id = *iIt;
          xt = sqrt( pow( particles_atTimeNow[id].X, 2 ) + pow( particles_atTimeNow[id].Y, 2 ) );
          
          nc = rings.getIndexPure( xt );
          
          if ( nc < rings.size() )
          {
            if ( rings[nc].getEnergyDensity() < theConfig->getFreezeOutEnergyDensity() )
            {
              particles_atTimeNow[id].free = true;
            }
            else
            {
              particles_atTimeNow[id].free = false;
            }
          }
          else
          {
            particles_atTimeNow[id].free = true;
          }
          
          if ( !particles_atTimeNow[id].free )
          {
            free = false;
          }
        }
        
        for ( iIt = cellsAdded[j].particleList.begin(); iIt != cellsAdded[j].particleList.end(); iIt++ )
        {
          int id = *iIt;
          xt = sqrt( pow( addedParticles[id].X, 2 ) + pow( addedParticles[id].Y, 2 ) );
          
          nc = rings.getIndexPure( xt );
          
          if ( nc < rings.size() )
          {
            if ( rings[nc].getEnergyDensity() < theConfig->getFreezeOutEnergyDensity() )
            {
              addedParticles[id].free = true;
            }
            else
            {
              addedParticles[id].free = false;
            }
          }
          else
          {
            addedParticles[id].free = true;
          }
          
          if ( !addedParticles[id].free )
          {
            free = false;
          }
        }
        //------------------------- check whether all particles in this cell are free --------------------
        
        if ( free )
        {
          list<int>::const_iterator iIt;
          for ( iIt = cells[j].particleList.begin(); iIt != cells[j].particleList.end(); iIt++ )
          {
            id = *iIt;
            particles_atTimeNow[id].edge = false;
          }
          
          for ( iIt = cellsAdded[j].particleList.begin(); iIt != cellsAdded[j].particleList.end(); iIt++ )
          {
            id = *iIt;
            addedParticles[id].edge = false;
            
            cc = ( nexttime - addedParticles[id].T ) / addedParticles[id].E;
            addedParticles[id].T = nexttime;
            addedParticles[id].X = addedParticles[id].X + addedParticles[id].PX * cc;
            addedParticles[id].Y = addedParticles[id].Y + addedParticles[id].PY * cc;
            addedParticles[id].Z = addedParticles[id].Z + addedParticles[id].PZ * cc;
          }
          cells[j].rates.normalizeRates();
          cellsAdded[j].rates.normalizeRates();
        }
        else
        {
          //             if (( cells[j].size() + cellsAdded[j].size() ) >= cellcut )     // enough particles in cell -> scatter
          if ( true )
          {
            nGluons = 0;
            nGluonsAdded = 0;
            free = false;
            vector<int> nQuarks( Particle::N_light_flavor, 0 );
            vector<int> nAntiQuarks( Particle::N_light_flavor, 0 );
            vector<int> nQuarksAdded( Particle::N_light_flavor, 0 );
            vector<int> nAntiQuarksAdded( Particle::N_light_flavor, 0 );
            
            list<int>::const_iterator iIt;
            for ( iIt = cells[j].particleList.begin(); iIt != cells[j].particleList.end(); iIt++ )
            {
              id = *iIt;
              particles_atTimeNow[id].edge = false;
              particles_atTimeNow[id].free = false;
              allParticlesList.push_back( id );
              
              if ( particles_atTimeNow[id].FLAVOR == 0 )
              {
                nGluons++;
                gluonList.push_back( *iIt );
              }
              else
              {
                switch ( particles_atTimeNow[id].FLAVOR )
                {
                  case up:
                    ++nQuarks[0];
                    break;
                  case down:
                    ++nQuarks[1];
                    break;
                  case strange:
                    ++nQuarks[2];
                    break;
                  case anti_up:
                    ++nAntiQuarks[0];
                    break;
                  case anti_down:
                    ++nAntiQuarks[1];
                    break;
                  case anti_strange:
                    ++nAntiQuarks[2];
                    break;
                  default:
                    break;
                }
              }
              
              xt = sqrt( pow( particles_atTimeNow[id].X, 2 ) + pow( particles_atTimeNow[id].Y, 2 ) );
              nc = rings.getIndex( xt );
              
              particles_atTimeNow[id].md2g = rings[nc].getAveraged_md2g();
              particles_atTimeNow[id].md2q = rings[nc].getAveraged_md2q();
              if ( particles_atTimeNow[id].FLAVOR == gluon )
              {
                particles_atTimeNow[id].rate = rateGluons[etaSliceIndex][nc];
                particles_atTimeNow[id].ratev = rateGluons_prev[etaSliceIndex][nc];
              }
              else if ( particles_atTimeNow[id].FLAVOR == up || particles_atTimeNow[id].FLAVOR == down || particles_atTimeNow[id].FLAVOR == strange || particles_atTimeNow[id].FLAVOR == light_quark )
              {
                particles_atTimeNow[id].rate = rateQuarks[etaSliceIndex][nc];
                particles_atTimeNow[id].ratev = rateQuarks_prev[etaSliceIndex][nc];
              }
              else
              {
                particles_atTimeNow[id].rate = rateAntiQuarks[etaSliceIndex][nc];
                particles_atTimeNow[id].ratev = rateAntiQuarks_prev[etaSliceIndex][nc];
              }
            }
            for ( iIt = cellsAdded[j].particleList.begin(); iIt != cellsAdded[j].particleList.end(); iIt++ )
            {
              id = *iIt;
              addedParticles[id].edge = false;
              addedParticles[id].free = false;
              allParticlesListAdded.push_back( id );
              
              if ( addedParticles[id].FLAVOR == 0 )
              {
                nGluonsAdded++;
                gluonListAdded.push_back( *iIt );
              }
              else
              {
                switch ( addedParticles[id].FLAVOR )
                {
                  case up:
                    ++nQuarksAdded[0];
                    break;
                  case down:
                    ++nQuarksAdded[1];
                    break;
                  case strange:
                    ++nQuarksAdded[2];
                    break;
                  case anti_up:
                    ++nAntiQuarksAdded[0];
                    break;
                  case anti_down:
                    ++nAntiQuarksAdded[1];
                    break;
                  case anti_strange:
                    ++nAntiQuarksAdded[2];
                    break;
                  default:
                    break;
                }
              }
              
              xt = sqrt( pow( addedParticles[id].X, 2 ) + pow( addedParticles[id].Y, 2 ) );
              nc = rings.getIndex( xt );
              
              addedParticles[id].md2g = rings[nc].getAveraged_md2g();
              addedParticles[id].md2q = rings[nc].getAveraged_md2q();
              if ( addedParticles[id].FLAVOR == gluon )
              {
                addedParticles[id].rate = rateGluons[etaSliceIndex][nc];
                addedParticles[id].ratev = rateGluons_prev[etaSliceIndex][nc];
              }
              else if ( addedParticles[id].FLAVOR == up || addedParticles[id].FLAVOR == down || addedParticles[id].FLAVOR == strange || addedParticles[id].FLAVOR == light_quark )
              {
                addedParticles[id].rate = rateQuarks[etaSliceIndex][nc];
                addedParticles[id].ratev = rateQuarks_prev[etaSliceIndex][nc];
              }
              else
              {
                addedParticles[id].rate = rateAntiQuarks[etaSliceIndex][nc];
                addedParticles[id].ratev = rateAntiQuarks_prev[etaSliceIndex][nc];
              }
            }
            
            int n32 = 0;
            
            scatt32_withAddedParticles( cells[j], allParticlesList, gluonList, cellsAdded[j], allParticlesListAdded, gluonListAdded, n32, again, aa, nexttime );
            if ( again )
            {
              return;
            }
            
            double scaleFactor = 1;
            if ( nGluons > n32 )
            {
              scaleFactor = static_cast<double>( nGluons ) / static_cast<double>( nGluons - n32 );
            }
            
            analysisRingStructure tempRing( aa.rings.size(), aa.rings.getCentralRadius(), aa.rings.getDeltaR() );
            scatt2223_withAddedParticles( cells[j], allParticlesList, gluonList, cellsAdded[j], allParticlesListAdded, gluonListAdded, scaleFactor, again, aa, nexttime, tempRing );
            if ( etaSliceIndex == centralEtaIndex )
            {
              aa.rings += tempRing;
            }
            
            if ( again )
            {
              return;
            }
            
            if ( again )
            {
              return;
            }
          }
          else  //if(nmb < cellcut)
            {
              for ( iIt = cells[j].particleList.begin(); iIt != cells[j].particleList.end(); iIt++ )
              {
                int id = *iIt;
                
                particles_atTimeNow[id].free = false;
                
                xt = sqrt( pow( particles_atTimeNow[id].X, 2 ) + pow( particles_atTimeNow[id].Y, 2 ) );
                nc = rings.getIndex( xt );
                
                if ( !particles_atTimeNow[id].edge )
                {
                  //      edgecell.add(id);
                  particles_atTimeNow[id].edge = true;
                  particles_atTimeNow[id].md2g = rings[nc].getAveraged_md2g();
                  particles_atTimeNow[id].md2q = rings[nc].getAveraged_md2q();
                }
                
                if ( particles_atTimeNow[id].FLAVOR == gluon )
                {
                  particles_atTimeNow[id].ratev = rateGluons[etaSliceIndex][nc];
                }
                else if ( particles_atTimeNow[id].FLAVOR == up || particles_atTimeNow[id].FLAVOR == down || particles_atTimeNow[id].FLAVOR == strange || particles_atTimeNow[id].FLAVOR == light_quark )
                {
                  particles_atTimeNow[id].ratev = rateQuarks[etaSliceIndex][nc];
                }
                else
                {
                  particles_atTimeNow[id].ratev = rateAntiQuarks[etaSliceIndex][nc];
                }
                particles_atTimeNow[id].rate = 0.0;
              }
              
              for ( iIt = cellsAdded[j].particleList.begin(); iIt != cellsAdded[j].particleList.end(); iIt++ )
              {
                int id = *iIt;
                addedParticles[id].free = false;
                
                xt = sqrt( pow( addedParticles[id].X, 2 ) + pow( addedParticles[id].Y, 2 ) );
                nc = rings.getIndex( xt );
                
                if ( !addedParticles[id].edge )
                {
                  //      edgecell.add(id);
                  addedParticles[id].edge = true;
                  addedParticles[id].md2g = rings[nc].getAveraged_md2g();
                  addedParticles[id].md2q = rings[nc].getAveraged_md2q();
                }
                
                if ( addedParticles[id].FLAVOR == gluon )
                {
                  addedParticles[id].ratev = rateGluons[etaSliceIndex][nc];
                }
                else if ( addedParticles[id].FLAVOR == up || addedParticles[id].FLAVOR == down || addedParticles[id].FLAVOR == strange || addedParticles[id].FLAVOR == light_quark )
                {
                  addedParticles[id].ratev = rateQuarks[etaSliceIndex][nc];
                }
                else
                {
                  addedParticles[id].ratev = rateAntiQuarks[etaSliceIndex][nc];
                }
                addedParticles[id].rate = 0.0;
              }
            }
        }
        
      }
      else if ( !cellsAdded[j].empty() ) // belongs to:  if ( !( cells[j].empty() || cellsAdded[j].empty() ) )
        {
          for ( iIt = cellsAdded[j].particleList.begin(); iIt != cellsAdded[j].particleList.end(); iIt++ )
          {
            int iscat = *iIt;
            cc = ( nexttime - addedParticles[iscat].T ) / addedParticles[iscat].E;
            addedParticles[iscat].T = nexttime;
            addedParticles[iscat].X = addedParticles[iscat].X + addedParticles[iscat].PX * cc;
            addedParticles[iscat].Y = addedParticles[iscat].Y + addedParticles[iscat].PY * cc;
            addedParticles[iscat].Z = addedParticles[iscat].Z + addedParticles[iscat].PZ * cc;
          }
        }
    }
  }
  
  formGeomCopy.clear();
}





void offlineHeavyIonCollision::scatt2223_withAddedParticles( cellContainer& _cell, std::vector< int >& _allParticlesList, std::vector< int >& _gluonList,
    cellContainer& _cellAdded, std::vector< int >& _allParticlesListAdded, std::vector< int >& _gluonListAdded,
    const double scaleFactor, bool& again, analysis& aa, const double nexttime, analysisRingStructure& _analysisRings )
{
  const double epsilon = 1.0e-4;
  int iscat, jscat, typ;
  double P1[4], P2[4];
  double s, as, csgg, cs23, cs22, cs23t, cs22t, Vrel, lambda_scaled;
  double probab22, probab23, probab2322;
  double averagedRate;
  double xt;
  double betaDistEntry;
  double md2g, md2q;
  int ringIndex;
  int initialStateIndex = -1;
  FLAVOR_TYPE F1, F2;


  const int nTotal = _cell.particleList.size();
  int nGluons = 0;
  vector<int> nQuarks( Particle::N_light_flavor, 0 );
  vector<int> nAntiQuarks( Particle::N_light_flavor, 0 );

  list<int>::const_iterator iIt;
  list<int>::const_iterator jIt;

  scattering23 scatt23_object( &theI23 );
  scattering22 scatt22_object;
  
  double lambdaJet = 0;
  double pt_addedParticle = 0;

  //   const int nAllQuarks = std::accumulate( nQuarks.begin(), nQuarks.end(), 0 ) + std::accumulate( nAntiQuarks.begin(), nAntiQuarks.end(), 0 );
//   const int allPairs = binomial( nTotal, 2 );
//   const int consideredPairs = 25;

  for ( jIt = _cellAdded.particleList.begin(); jIt != _cellAdded.particleList.end(); jIt++ ) 
  {
    jscat = *jIt;
    pt_addedParticle = sqrt( pow( addedParticles[jscat].PX, 2.0 ) + pow( addedParticles[jscat].PY, 2.0 ) );
    
    if ( pt_addedParticle < 4.0 )
    {
      continue; // jump to next particle in the list
    }
    
    if ( pt_addedParticle > 8.0 )
    {
      lambdaJet = iterateMFP( _allParticlesList, _gluonList, jscat, dt, dv );
      
      xt = sqrt( pow( addedParticles[jscat].X, 2 )  + pow( addedParticles[jscat].Y, 2 ) );
      ringIndex = rings.getIndex( xt );
      
      if ( addedParticles[jscat].FLAVOR == gluon )
      {
        _analysisRings[ringIndex].lambdaGluon += lambdaJet;
        _analysisRings[ringIndex].collectedGluon++;
      }
      else
      {
        _analysisRings[ringIndex].lambdaQuark += lambdaJet;
        _analysisRings[ringIndex].collectedQuark++;
      }
      
//       lambdaJet = 0;
    }
    else
    {
      lambdaJet = -1;
    }

    for ( iIt =  _cell.particleList.begin(); iIt != _cell.particleList.end(); iIt++ )
    {
      iscat = *iIt;  

      F1 = particles_atTimeNow[iscat].FLAVOR;
      particles_atTimeNow[iscat].getMomentumArray( P1 );
      F2 = addedParticles[jscat].FLAVOR;
      addedParticles[jscat].getMomentumArray( P2 );

      s = pow(( P1[0] + P2[0] ), 2 ) - pow(( P1[1] + P2[1] ), 2 ) - pow(( P1[2] + P2[2] ), 2 ) - pow(( P1[3] + P2[3] ), 2 );

      xt = ( sqrt( pow( particles_atTimeNow[iscat].X, 2 )  + pow( particles_atTimeNow[iscat].Y, 2 ) ) +
             sqrt( pow( addedParticles[jscat].X, 2 )  + pow( addedParticles[jscat].Y, 2 ) ) ) / 2.0;

      ringIndex = rings.getIndex( xt );

      averagedRate = ( particles_atTimeNow[iscat].rate + addedParticles[jscat].rate +
                       particles_atTimeNow[iscat].ratev + addedParticles[jscat].ratev ) / ( 2.0 * 2.0 );


      if ( s < 1.1*lambda2 )
      {
        probab2322 = -1.0;
        cs22 = cs22t = cs23 = cs23t = 0.0; //1/GeV^2
      }
      else
      {
        Vrel = s / ( 2.0 * P1[0] * P2[0] );
        as = alpha_s( s );

        //factor as (alpha_s) is not included in definitions of partcl[jscat].md2g, therefore multiplied here
        md2g = as * ( particles_atTimeNow[iscat].md2g + addedParticles[jscat].md2g ) / 2.0;
        md2q = as * ( particles_atTimeNow[iscat].md2q + addedParticles[jscat].md2q ) / 2.0;
        
        // HACK for N_f = 0 background!!!!
//         md2g = as * ( particles_atTimeNow[iscat].md2g + addedParticles[jscat].md2g ) / 2.0 * 2;
//         md2q = md2g * 2.0 / 9.0;

        scatt22_object.setParameter( P1, P2, F1, F2, s, md2g, md2q );
        cs22 = scatt22_object.getXSection22( initialStateIndex );

        probab22 = pow( 0.197, 2.0 ) * cs22 * Vrel * dt / ( dv * testpartcl );

//         probab23 = 0;
        if ( lambdaJet > 0 )
        {
          lambda_scaled = lambdaJet * sqrt( s );
        }
        else
        {
          if ( averagedRate > epsilon )
          {
            lambda_scaled = sqrt( s ) / averagedRate;  //dimensionless
          }
          else
          {
            xsection_gg_gg csObj( s, md2g, md2q );
            csgg = csObj.totalCrossSection();
            lambda_scaled = ( dv * rings[ringIndex].getGamma() * testpartcl * sqrt( s ) ) / ( nTotal * csgg * pow( 0.197, 3.0 ) );
          }
        }
        
        double vx = rings[ringIndex].getAveraged_v_x();
        double vy = rings[ringIndex].getAveraged_v_y();
        double vz = rings[ringIndex].getAveraged_v_z();

        betaDistEntry = scatt23_object.setParameter( vx, vy, vz, P1, P2, F1, F2, sqrt( s ), md2g / s, lambda_scaled, _gluonList.size() );

        cs23 = 1 / s * Ncolor * pow( as, 3 ) * scatt23_object.getIntegral23();

        probab23 = pow( 0.197, 2.0 ) * cs23 * Vrel * dt / ( dv * testpartcl );
        if ( F1 == gluon )
        {
          probab22 *= scaleFactor;
          probab23 *= scaleFactor;
        }
        if ( F2 == gluon )
        {
          probab22 *= scaleFactor;
          probab23 *= scaleFactor;
        }

        if ( pt_addedParticle > 10.0 )
        {
          if ( ParticleOffline::mapToGenericFlavorType( F2 ) == gluon )
          {
            ++n23_collected_gluon;
            p23_collected_gluon += probab23;
            lambdaJet_gluon += lambda_scaled;
          }
          else if ( ParticleOffline::mapToGenericFlavorType( F2 ) == light_quark )
          {
            ++n23_collected_quark;
            p23_collected_quark += probab23;
            lambdaJet_quark += lambda_scaled;
          }
        }

        cs22t = 18.0 * M_PI * as * as / s * ( log( 1.0 + 0.25 / md2g ) + 1.0 / ( 1.0 + 0.25 / md2g ) - 1.0 );//1/GeV^2

        if ( cs22 > 0.0 )
        {
          ++_cellAdded.nCollected22;
          _cellAdded.md2g_scaled_22 += md2g / s;
          _cellAdded.md2q_scaled_22 += md2q / s;
        }
        if ( cs23 > 0.0 )
        {
          ++_cellAdded.nCollected23;
          _cellAdded.md2g_scaled_23 += md2g / s;
          _cellAdded.md2q_scaled_23 += md2q / s;
          _cellAdded.lambdaScaled += lambda_scaled;
        }


        probab2322 = probab22 + probab23;
      }

      ++_cellAdded.nCollectedAll2223;
      _cellAdded.sigma_22 += cs22;               //1/GeV^2
      _cellAdded.sigma_23 += cs23;               //1/GeV^2

      if ( probab2322 > 1.0 )
      {
//         cout << "P2322=" << probab2322 << ">1" << endl;
        again = true;
//         cout << "dt (old) = " << dt << endl;
        dt = 0.5 / ( probab2322 / dt );
//         cout << "dt (new) = " << dt << endl;
        return;
      }

      if ( ran2() < probab2322 )
      {
        double pt_jscat = sqrt( pow( P2[1], 2.0 ) + pow( P2[2], 2.0 ) );
        double pt_nmb;

        if ( ran2() * probab2322 < probab23 )
        {
          int jetEventIndex = -1;
          if ( pt_jscat > aa.getJetTracking_PT() )
          {
            jetEventIndex = aa.addJetEvent_in( iscat, -1, jscat, c2to3, cs23, _cell.index, lambda_scaled / sqrt( s ) );
          }

          int newIndex = scatt23_utility( scatt23_object, _cell, iscat, jscat, again, nexttime );

          pt_jscat = sqrt( pow( addedParticles[jscat].PX, 2.0 ) + pow( addedParticles[jscat].PY, 2.0 ) );
          pt_nmb = sqrt( pow( addedParticles[newIndex].PX, 2.0 ) + pow( addedParticles[newIndex].PY, 2.0 ) );
          if ( jetEventIndex != -1 || pt_jscat > aa.getJetTracking_PT() || pt_nmb > aa.getJetTracking_PT() )
          {
            aa.addJetEvent_out( jetEventIndex, jscat, iscat, newIndex, c2to3 );
          }
        }
        else
        {
          int jetEventIndex = -1;
          if ( pt_jscat > aa.getJetTracking_PT() )
          {
            jetEventIndex = aa.addJetEvent_in( iscat, -1, jscat, c2to2, cs22, _cell.index, lambda_scaled / sqrt( s ) );
          }

          scatt22_utility( scatt22_object, iscat, jscat, typ, nexttime );

          pt_jscat = sqrt( pow( addedParticles[jscat].PX, 2.0 ) + pow( addedParticles[jscat].PY, 2.0 ) );
          if ( jetEventIndex != -1 || pt_jscat > aa.getJetTracking_PT() )
          {
            aa.addJetEvent_out( jetEventIndex, jscat, iscat, -1, c2to2 );
          }
        }
      }
    }
  }

  double cc;
  for ( iIt =  _cellAdded.particleList.begin(); iIt != _cellAdded.particleList.end(); iIt++ )
  {
    iscat = *iIt;
    cc = ( nexttime - addedParticles[iscat].T ) / addedParticles[iscat].E;
    addedParticles[iscat].T = nexttime;
    addedParticles[iscat].X = addedParticles[iscat].X + addedParticles[iscat].PX * cc;
    addedParticles[iscat].Y = addedParticles[iscat].Y + addedParticles[iscat].PY * cc;
    addedParticles[iscat].Z = addedParticles[iscat].Z + addedParticles[iscat].PZ * cc;
  }
}




void offlineHeavyIonCollision::scatt32_withAddedParticles( cellContainer& _cell, std::vector< int >& _allParticlesList, std::vector< int >& _gluonList, cellContainer& _cellAdded, std::vector< int >& _allParticlesListAdded, std::vector< int >& _gluonListAdded, int& n32, bool& again, analysis& aa, const double nexttime )
{
  const double epsilon = 1.0e-4;
  int iscat, jscat, kscat;
  int m1, m2, m3;
  int ringIndex;
  int consideredTriplets;  //number of triplets to consider

  FLAVOR_TYPE F1, F2, F3;

  double P1[4], P2[4], P3[4];
  double averagedRate;
  double s, as, csgg, I32, probab32, e1, e3, cg, xt, ratio;
  double md2g, md2q;
  double lambda_scaled;
  double ran2out;
  double betaDistEntry;
  int order = -1;

  //n32=0;

  const int nTotal = _cell.particleList.size();
  int nTotalAdded = _allParticlesListAdded.size();
  const int nGluons = _gluonList.size();
  int nGluonsAdded = _gluonListAdded.size();
  const int nAllQuarks = nTotal - nGluons;
  int nAllQuarksAdded = nTotalAdded - nGluonsAdded;
  
  for ( int i = 0; i < _allParticlesListAdded.size(); i++ )
  {
    double pt_addedParticle =  sqrt( pow( addedParticles[_allParticlesListAdded[i]].PX, 2.0 ) + pow( addedParticles[_allParticlesListAdded[i]].PY, 2.0 ) ); 
    if ( pt_addedParticle < 4.0 )
    {
      -- nTotalAdded;
      if ( addedParticles[_allParticlesListAdded[i]].FLAVOR == gluon )
      {
        --nGluonsAdded;
      }
      else
      {
        --nAllQuarksAdded;
      }
    }
  }

  const int allTriplets = ( binomial( nTotal, 2 ) * nTotalAdded )  - ( binomial( nAllQuarks, 2 ) * nAllQuarksAdded );

  scattering32 scatt32_object;
  
  double lambdaJet = -1;
  double pt_addedParticle = 0;

  if ( allTriplets > 20 )
  {
    //----------- compute iterated MFP for one (the first) jet particle -------------
    // (in general there shouldn't be much more)
    for ( int jj = 0; jj < _allParticlesListAdded.size(); jj++ )
    {
      kscat = _allParticlesListAdded[jj];
      pt_addedParticle = sqrt( pow( addedParticles[kscat].PX, 2.0 ) + pow( addedParticles[kscat].PY, 2.0 ) ); 
    
      if ( pt_addedParticle > 8 )
      {
        lambdaJet = iterateMFP( _allParticlesList, _gluonList, kscat, dt, dv );      
        break;
      }
    }
    //--------------------------------------------------------------------------
    
    if ( nTotalAdded * 3 >= 20 )
    {
      consideredTriplets = nTotalAdded * 3;  // number of triplets to consider for 3 -> 2 processes
    }
    else
    {
      consideredTriplets = 20;  //20 is an empirical value
    }

    double scaleForSelectedTriplets = static_cast<double>( allTriplets ) / static_cast<double>( consideredTriplets );

    for ( int i = 0; i < consideredTriplets && ( _gluonList.size() + _gluonListAdded.size() ) > 0 && _allParticlesListAdded.size() > 0 ; i++ )
    {
      do
      {
        do
        {
          m1 = int ( _allParticlesList.size() * ran2() );
          if ( m1 == _allParticlesList.size() )
          {
            m1 = _allParticlesList.size() - 1;
          }
          iscat = _allParticlesList[m1];
        }
        while ( particles_atTimeNow[iscat].dead );
        F1 = particles_atTimeNow[iscat].FLAVOR;

        do
        {
          m2 = int ( _allParticlesList.size() * ran2() );
          if ( m2 == _allParticlesList.size() )
          {
            m2 = _allParticlesList.size() - 1;
          }
          jscat = _allParticlesList[m2];
        }
        while (( m2 == m1 ) || ( particles_atTimeNow[jscat].dead ) );
        F2 = particles_atTimeNow[jscat].FLAVOR;

        do
        {
          m3 = int ( _allParticlesListAdded.size() * ran2() );
          if ( m3 == _allParticlesListAdded.size() )
          {
            m3 = _allParticlesListAdded.size() - 1;
          }
          kscat = _allParticlesListAdded[m3];
        }
        while ( addedParticles[kscat].dead );
        F3 = addedParticles[kscat].FLAVOR;
      }
      while ( !( F1 == gluon || F2 == gluon || F3 == gluon ) );
      
      pt_addedParticle = sqrt( pow( addedParticles[kscat].PX, 2.0 ) + pow( addedParticles[kscat].PY, 2.0 ) ); 
      if ( pt_addedParticle < 4.0 )
      {
        continue; //go to next particle triplet 
      }

      particles_atTimeNow[iscat].getMomentumArray( P1 );
      particles_atTimeNow[jscat].getMomentumArray( P2 );
      addedParticles[kscat].getMomentumArray( P3 );

      s = pow(( P1[0] + P2[0] + P3[0] ), 2 ) - pow(( P1[1] + P2[1] + P3[1] ), 2 ) - pow(( P1[2] + P2[2] + P3[2] ), 2 ) - pow(( P1[3] + P2[3] + P3[3] ), 2 );

      averagedRate = ( particles_atTimeNow[iscat].rate + particles_atTimeNow[jscat].rate + addedParticles[kscat].rate +
                       particles_atTimeNow[iscat].ratev + particles_atTimeNow[jscat].ratev + addedParticles[kscat].ratev ) / ( 3.0 * 2.0 );

      if ( s < 1.1*lambda2 )
      {
        probab32 = -1.0;
      }
      else
      {
        as = alpha_s( s );

        //factor as (alpha_s) is not included in definitions of partcl[jscat].md2g, therefore multiplied here
        md2g = as * ( particles_atTimeNow[iscat].md2g + particles_atTimeNow[jscat].md2g + addedParticles[kscat].md2g ) / 3.0;
        md2q = as * ( particles_atTimeNow[iscat].md2q + particles_atTimeNow[jscat].md2q + addedParticles[kscat].md2q ) / 3.0;

        // HACK for N_f = 0 background!!!!
//         md2g = as * ( particles_atTimeNow[iscat].md2g + particles_atTimeNow[jscat].md2g + addedParticles[kscat].md2g ) / 3.0 * 2;
//         md2q = md2g * 2.0 / 9.0;
        
        xt = ( sqrt( pow( particles_atTimeNow[iscat].X, 2 )  + pow( particles_atTimeNow[iscat].Y, 2 ) ) +
               sqrt( pow( particles_atTimeNow[jscat].X, 2 )  + pow( particles_atTimeNow[jscat].Y, 2 ) ) +
               sqrt( pow( addedParticles[kscat].X, 2 )  + pow( addedParticles[kscat].Y, 2 ) ) ) / 3.0;

        ringIndex = rings.getIndex( xt );

        if ( lambdaJet > 0 )
        {
          lambda_scaled = lambdaJet * sqrt( s );
        }
        else
        {
          if ( averagedRate > epsilon )
          {
            lambda_scaled = sqrt( s ) / averagedRate;  //dimensionless
          }
          else
          {
            xsection_gg_gg csObj( s, md2g, md2q );
            csgg = csObj.totalCrossSection();
            lambda_scaled = ( dv * rings[ringIndex].getGamma() * testpartcl * sqrt( s ) ) / ( nTotal * csgg * pow( 0.197, 3.0 ) );
          }
        }

        // create scattering32 object for the given 3 particles
        double vx = rings[ringIndex].getAveraged_v_x();
        double vy = rings[ringIndex].getAveraged_v_y();
        double vz = rings[ringIndex].getAveraged_v_z();

        betaDistEntry = scatt32_object.setParameter( vx, vy, vz, P1, P2, P3, F1, F2, F3, sqrt( s ), md2g / s, lambda_scaled, _gluonList.size() );
        I32 = scatt32_object.getIntegral32_fast();                        // get the integral I32 for the given 3 particles

        probab32 = scaleForSelectedTriplets * pow( 0.197, 5.0 ) * 9.0 * M_PI * Ncolor * pow( as, 3.0 ) * dt * I32 / ( gG * pow( dv, 2 ) * pow( testpartcl, 2 ) * s * P1[0] * P2[0] * P3[0] );
      }

      if ( probab32 > 1.0 )
      {
        cout << "P32 = " << probab32 << " > 1" << endl;
        cout << "dt (old) = " << dt << endl;
        again = true;
        dt = 0.5 / ( probab32 / dt );
        cout << "dt (new) = " << dt << endl;
        return;
      }

      ran2out = ran2();
      if ( ran2out < probab32 )
      {
        double pt_iscat, pt_jscat;
        double pt_kscat = sqrt( pow( P3[1], 2.0 ) + pow( P3[2], 2.0 ) );

        int jetEventIndex = -1;
        if ( pt_kscat > aa.getJetTracking_PT() )
        {
          jetEventIndex = aa.addJetEvent_in( iscat, jscat, kscat, c3to2, I32, _cell.index, lambda_scaled / sqrt( s ) );
        }

        order = scatt32_utility( scatt32_object, _cellAdded.particleList, _allParticlesListAdded, _gluonListAdded, iscat, jscat, kscat, n32, ran2out / probab32, nexttime );

        if ( scatt32_object.getCollisionStatus() )
        {
          if ( order == 4 || order == 6 )
          {
            pt_jscat = sqrt( pow( particles_atTimeNow[jscat].PX, 2.0 ) + pow( particles_atTimeNow[jscat].PY, 2.0 ) );
            pt_kscat = sqrt( pow( addedParticles[kscat].PX, 2.0 ) + pow( addedParticles[kscat].PY, 2.0 ) );
            if ( jetEventIndex != -1 || pt_jscat > aa.getJetTracking_PT() || pt_kscat > aa.getJetTracking_PT() )
            {
              aa.addJetEvent_out( jetEventIndex, kscat, jscat, -1, c3to2 );
            }
          }
          else if ( order == 2 || order == 5 )
          {
            pt_iscat = sqrt( pow( particles_atTimeNow[iscat].PX, 2.0 ) + pow( particles_atTimeNow[iscat].PY, 2.0 ) );
            pt_kscat = sqrt( pow( addedParticles[kscat].PX, 2.0 ) + pow( addedParticles[kscat].PY, 2.0 ) );
            if ( jetEventIndex != -1 || pt_iscat > aa.getJetTracking_PT() || pt_kscat > aa.getJetTracking_PT() )
            {
              aa.addJetEvent_out( jetEventIndex, kscat, iscat, -1, c3to2 );
            }
          }
          else if ( order == 1 || order == 3 )
          {
            pt_iscat = sqrt( pow( particles_atTimeNow[iscat].PX, 2.0 ) + pow( particles_atTimeNow[iscat].PY, 2.0 ) );
            pt_jscat = sqrt( pow( particles_atTimeNow[jscat].PX, 2.0 ) + pow( particles_atTimeNow[jscat].PY, 2.0 ) );
            if ( jetEventIndex != -1 || pt_iscat > aa.getJetTracking_PT() || pt_jscat > aa.getJetTracking_PT() )
            {
              aa.addJetEvent_out( jetEventIndex, kscat, jscat, -1, c3to2 );
            }
          }
        }
        else
        {
          if ( jetEventIndex != -1 )
          {
            aa.removeJetEvent_in( jetEventIndex );
          }
        }
      }
    }
  }
  else
  {
    for ( int m3 = 0; m3 < _allParticlesListAdded.size(); m3++ )
    {
      kscat = _allParticlesListAdded[m3];
      pt_addedParticle = sqrt( pow( addedParticles[kscat].PX, 2.0 ) + pow( addedParticles[kscat].PY, 2.0 ) ); 
      if ( pt_addedParticle < 4.0 )
      {
        continue; //go to next particle triplet 
      }
      
      if ( pt_addedParticle > 8.0 )
      {
        lambdaJet = iterateMFP( _allParticlesList, _gluonList, kscat, dt, dv );
      }
      else
      {
        lambdaJet = -1;
      }
      
      for ( int m1 = 0; m1 < _allParticlesList.size() - 1; m1++ )
      {
        iscat = _allParticlesList[m1];
  
        for ( int m2 = m1 + 1; m2 < _allParticlesList.size(); m2++ )
        {
          jscat = _allParticlesList[m2];

          F1 = particles_atTimeNow[iscat].FLAVOR;
          F2 = particles_atTimeNow[jscat].FLAVOR;
          F3 = addedParticles[kscat].FLAVOR;

          // at least one of the particles must be a gluon
          // otherwise go to next step in the loop
          if ( !( F1 == gluon || F2 == gluon || F3 == gluon ) )
          {
            continue;
          }

          particles_atTimeNow[iscat].getMomentumArray( P1 );
          particles_atTimeNow[jscat].getMomentumArray( P2 );
          addedParticles[kscat].getMomentumArray( P3 );

          s = pow(( P1[0] + P2[0] + P3[0] ), 2 ) - pow(( P1[1] + P2[1] + P3[1] ), 2 ) - pow(( P1[2] + P2[2] + P3[2] ), 2 ) - pow(( P1[3] + P2[3] + P3[3] ), 2 );

          averagedRate = ( particles_atTimeNow[iscat].rate + particles_atTimeNow[jscat].rate + addedParticles[kscat].rate +
                           particles_atTimeNow[iscat].ratev + particles_atTimeNow[jscat].ratev + addedParticles[kscat].ratev ) / ( 3.0 * 2.0 );

          if ( s < 1.1*lambda2 )
          {
            probab32 = -1.0;
          }
          else
          {
            as = alpha_s( s );

            //factor as (alpha_s) is not included in definitions of partcl[jscat].md2g, therefore multiplied here
            md2g = as * ( particles_atTimeNow[iscat].md2g + particles_atTimeNow[jscat].md2g + addedParticles[kscat].md2g ) / 3.0;
            md2q = as * ( particles_atTimeNow[iscat].md2q + particles_atTimeNow[jscat].md2q + addedParticles[kscat].md2q ) / 3.0;

            // HACK for N_f = 0 background!!!!
//             md2g = as * ( particles_atTimeNow[iscat].md2g + particles_atTimeNow[jscat].md2g + addedParticles[kscat].md2g ) / 3.0 * 2;
//             md2q = md2g * 2.0 / 9.0;
            
            xt = ( sqrt( pow( particles_atTimeNow[iscat].X, 2 )  + pow( particles_atTimeNow[iscat].Y, 2 ) ) +
                   sqrt( pow( particles_atTimeNow[jscat].X, 2 )  + pow( particles_atTimeNow[jscat].Y, 2 ) ) +
                   sqrt( pow( addedParticles[kscat].X, 2 )  + pow( addedParticles[kscat].Y, 2 ) ) ) / 3.0;

            ringIndex = rings.getIndex( xt );

            if ( lambdaJet > 0 )
            {
              lambda_scaled = lambdaJet * sqrt( s );
            }
            else
            {
              if ( averagedRate > epsilon )
              {
                lambda_scaled = sqrt( s ) / averagedRate;  //dimensionless
              }
              else
              {
                xsection_gg_gg csObj( s, md2g, md2q );
                csgg = csObj.totalCrossSection();
                lambda_scaled = ( dv * rings[ringIndex].getGamma() * testpartcl * sqrt( s ) ) / ( nTotal * csgg * pow( 0.197, 3.0 ) );
              }
            }

            // create scattering32 object for the given 3 particles
            double vx = rings[ringIndex].getAveraged_v_x();
            double vy = rings[ringIndex].getAveraged_v_y();
            double vz = rings[ringIndex].getAveraged_v_z();

            betaDistEntry = scatt32_object.setParameter( vx, vy, vz, P1, P2, P3, F1, F2, F3, sqrt( s ), md2g / s, lambda_scaled, _gluonList.size() );
            I32 = scatt32_object.getIntegral32_fast();                        // get the integral I32 for the given 3 particles

            probab32 = pow( 0.197, 5.0 ) * 9.0 * M_PI * Ncolor * pow( as, 3.0 ) * dt * I32 / ( gG * pow( dv, 2 ) * pow( testpartcl, 2 ) * s * P1[0] * P2[0] * P3[0] );
          }

          if ( probab32 > 1.0 )
          {
            cout << "P32 = " << probab32 << " > 1" << endl;
            cout << "dt (old) = " << dt << endl;
            again = true;
            dt = 0.5 / ( probab32 / dt );
            cout << "dt (new) = " << dt << endl;
            return;
          }

          ran2out = ran2();
          if ( ran2out < probab32 )
          {
            double pt_iscat, pt_jscat;
            double pt_kscat = sqrt( pow( P3[1], 2.0 ) + pow( P3[2], 2.0 ) );

            int jetEventIndex = -1;
            if ( pt_kscat > aa.getJetTracking_PT() )
            {
              jetEventIndex = aa.addJetEvent_in( iscat, jscat, kscat, c3to2, I32, _cell.index, lambda_scaled / sqrt( s ) );
            }

            order = scatt32_utility( scatt32_object, _cellAdded.particleList, _allParticlesListAdded, _gluonListAdded, iscat, jscat, kscat, n32, ran2out / probab32, nexttime );

            if ( scatt32_object.getCollisionStatus() )
            {
              if ( order == 4 || order == 6 )
              {
                pt_jscat = sqrt( pow( particles_atTimeNow[jscat].PX, 2.0 ) + pow( particles_atTimeNow[jscat].PY, 2.0 ) );
                pt_kscat = sqrt( pow( addedParticles[kscat].PX, 2.0 ) + pow( addedParticles[kscat].PY, 2.0 ) );
                if ( jetEventIndex != -1 || pt_jscat > aa.getJetTracking_PT() || pt_kscat > aa.getJetTracking_PT() )
                {
                  aa.addJetEvent_out( jetEventIndex, kscat, jscat, -1, c3to2 );
                }
              }
              else if ( order == 2 || order == 5 )
              {
                pt_iscat = sqrt( pow( particles_atTimeNow[iscat].PX, 2.0 ) + pow( particles_atTimeNow[iscat].PY, 2.0 ) );
                pt_kscat = sqrt( pow( addedParticles[kscat].PX, 2.0 ) + pow( addedParticles[kscat].PY, 2.0 ) );
                if ( jetEventIndex != -1 || pt_iscat > aa.getJetTracking_PT() || pt_kscat > aa.getJetTracking_PT() )
                {
                  aa.addJetEvent_out( jetEventIndex, kscat, iscat, -1, c3to2 );
                }
              }
              else if ( order == 1 || order == 3 )
              {
                m3--;

                pt_iscat = sqrt( pow( particles_atTimeNow[iscat].PX, 2.0 ) + pow( particles_atTimeNow[iscat].PY, 2.0 ) );
                pt_jscat = sqrt( pow( particles_atTimeNow[jscat].PX, 2.0 ) + pow( particles_atTimeNow[jscat].PY, 2.0 ) );
                if ( jetEventIndex != -1 || pt_iscat > aa.getJetTracking_PT() || pt_jscat > aa.getJetTracking_PT() )
                {
                  aa.addJetEvent_out( jetEventIndex, kscat, jscat, -1, c3to2 );
                }
              }
            }
            else
            {
              if ( jetEventIndex != -1 )
              {
                aa.removeJetEvent_in( jetEventIndex );
              }
            }
          }
        }
      }
    }
  }
}





int offlineHeavyIonCollision::scatt23_utility( scattering23& scatt23_obj, cellContainer& _cell, int iscat, const int jscat, bool& again, const double nexttime )
{
  FLAVOR_TYPE F1, F2;
  int nc, IXY;
  double cx, cy, cz, dz;
  double P1[4], P2[4], P3[4], R1[4], R2[4];
  double Tmax, TT, cc, pt1, pt3, y, phi, pz1;
  int typ;

  int centralEtaIndex = static_cast<int>( etaBins.size() ) / 2;
  int etaIndex = centralEtaIndex;
  double eta = 0.5 * log(( particles_atTimeNow[iscat].T + particles_atTimeNow[iscat].Z ) / ( particles_atTimeNow[iscat].T - particles_atTimeNow[iscat].Z ) );
  if ( eta >= 0 )
  {
    while ( eta >= etaBins[etaIndex].right && etaIndex < IZ )
    {
      ++etaIndex;
    }
  }
  else
  {
    while ( eta < etaBins[etaIndex].left && etaIndex >= 0 )
    {
      --etaIndex;
    }
  }
  double leftZ = timenow * tanh( etaBins[etaIndex].left );
  double deltaZ = timenow * tanh( etaBins[etaIndex].right ) - leftZ;
  double leftY = _cell.corner.y_min;
  double leftX = _cell.corner.x_min;

  particles_atTimeNow[iscat].getMomentumArray( P1 );
  particles_atTimeNow[iscat].getCoordinateArray( R1 );
  F1 = particles_atTimeNow[iscat].FLAVOR;

  addedParticles[jscat].getMomentumArray( P2 );
  addedParticles[jscat].getCoordinateArray( R2 );
  F2 = addedParticles[jscat].FLAVOR;

  ncoll++;
  ncoll23++;
  addedParticles[jscat].coll_id = ncoll;

  if ( R1[0] > R2[0] )
    Tmax = R1[0];
  else
    Tmax = R2[0];
  TT = ( nexttime - Tmax ) * ran2() + Tmax;

  cc = ( TT - R2[0] ) / P2[0];
  addedParticles[jscat].T = TT;
  addedParticles[jscat].X = R2[1] + P2[1] * cc;
  addedParticles[jscat].Y = R2[2] + P2[2] * cc;
  addedParticles[jscat].Z = R2[3] + P2[3] * cc;

  addedParticles[jscat].X_traveled += sqrt( pow( P2[1] * cc, 2.0 ) + pow( P2[2] * cc, 2.0 ) + pow( P2[3] * cc, 2.0 ) );

  // at this point scattering take place, therefore set properties of last scattering point
  addedParticles[jscat].X_lastInt = addedParticles[jscat].X;
  addedParticles[jscat].Y_lastInt = addedParticles[jscat].Y;
  addedParticles[jscat].Z_lastInt = addedParticles[jscat].Z;
  addedParticles[jscat].T_lastInt = addedParticles[jscat].T;

  scatt23_obj.getMomenta23( pt1, pt3, y, phi, pz1, typ, F1, F2 );
  scatt23_obj.setNewMomenta23( P1, P2, P3, R1, R2, pt1, pt3, y, phi, pz1 );

  double pt_out1 = sqrt( pow( P1[1], 2 ) + pow( P1[2], 2 ) );
  double pt_out2 = sqrt( pow( P2[1], 2 ) + pow( P2[2], 2 ) );

  //<<---------------------------------------------
  // set new properties for added particle
  if ( pt_out1 > pt_out2 )
  {
    addedParticles[jscat].FLAVOR = F1;
    addedParticles[jscat].PX = P1[1];
    addedParticles[jscat].PY = P1[2];
    addedParticles[jscat].PZ = P1[3];
    addedParticles[jscat].E = sqrt( pow( P1[1], 2 ) + pow( P1[2], 2 ) + pow( P1[3], 2 ) );
  }
  else
  {
    addedParticles[jscat].FLAVOR = F2;
    addedParticles[jscat].PX = P2[1];
    addedParticles[jscat].PY = P2[2];
    addedParticles[jscat].PZ = P2[3];
    addedParticles[jscat].E = sqrt( pow( P2[1], 2 ) + pow( P2[2], 2 ) + pow( P2[3], 2 ) );
  }

  ParticleOffline tempParticle;
  tempParticle.FLAVOR = gluon;
  tempParticle.PX = P3[1];
  tempParticle.PY = P3[2];
  tempParticle.PZ = P3[3];
  tempParticle.E = sqrt( P3[1] * P3[1] + P3[2] * P3[2] + P3[3] * P3[3] );
  tempParticle.T = TT;
  tempParticle.X = leftX + dx * ran2();
  tempParticle.Y = leftY + dy * ran2();
  tempParticle.Z = leftZ + deltaZ * ran2();
  tempParticle.X_lastInt = tempParticle.X;
  tempParticle.Y_lastInt = tempParticle.Y;
  tempParticle.Z_lastInt = tempParticle.Z;
  tempParticle.T_lastInt = tempParticle.T;
  tempParticle.X_init = tempParticle.X;
  tempParticle.Y_init = tempParticle.X;
  tempParticle.Z_init = tempParticle.X;
  tempParticle.PX_init = tempParticle.PX;
  tempParticle.PY_init = tempParticle.PY;
  tempParticle.PZ_init = tempParticle.PZ;
  tempParticle.E_init = -tempParticle.E;   // negative to indicate creation via 2->3 process
  tempParticle.m = 0.0;
  tempParticle.md2g = ( particles[iscat].md2g + particles[jscat].md2g ) / 2.0;
  tempParticle.md2q = ( particles[iscat].md2q + particles[jscat].md2q ) / 2.0;
  tempParticle.edge = false;
  tempParticle.dead = false;
  tempParticle.free = false;
  tempParticle.cell_id = -1;
  tempParticle.coll_id = ncoll;
  tempParticle.collisionTime = infinity;
  tempParticle.collisionPartner = -1;
  tempParticle.rate = particles[iscat].rate;    //GeV
  tempParticle.ratev = particles[iscat].ratev;    //GeV
  tempParticle.unique_id = ParticleOffline::unique_id_counter_added;
  --ParticleOffline::unique_id_counter_added;


//   if ( sqrt( pow( tempParticle.PX, 2) + pow( tempParticle.PY, 2) ) > 3.0 )
//   {
    addedParticles.push_back( tempParticle );
    int newIndex = addedParticles.size() - 1;
    ++numberAdded;

    cc = ( nexttime - addedParticles[newIndex].T ) / addedParticles[newIndex].E;
    addedParticles[newIndex].T = nexttime;
    addedParticles[newIndex].X = addedParticles[newIndex].X + addedParticles[newIndex].PX * cc;
    addedParticles[newIndex].Y = addedParticles[newIndex].Y + addedParticles[newIndex].PY * cc;
    addedParticles[newIndex].Z = addedParticles[newIndex].Z + addedParticles[newIndex].PZ * cc;
//   }
  
  return newIndex;
}



void offlineHeavyIonCollision::scatt22_utility( scattering22& scatt22_obj, const int iscat, const int jscat, int& typ, const double nexttime )
{
  FLAVOR_TYPE F1, F2;
  double P1[4], P2[4], R1[4], R2[4];
  double as;
  double Tmax, TT, cc, PT2, t;

  particles_atTimeNow[iscat].getMomentumArray( P1 );
  particles_atTimeNow[iscat].getCoordinateArray( R1 );
  F1 = particles_atTimeNow[iscat].FLAVOR;

  addedParticles[jscat].getMomentumArray( P2 );
  addedParticles[jscat].getCoordinateArray( R2 );
  F2 = addedParticles[jscat].FLAVOR;


  ncoll++;
  ncoll22++;
  addedParticles[jscat].coll_id = ncoll;

  if ( R1[0] > R2[0] )
    Tmax = R1[0];
  else
    Tmax = R2[0];
  TT = ( nexttime - Tmax ) * ran2() + Tmax;


  cc = ( TT - R2[0] ) / P2[0];
  addedParticles[jscat].T = TT;
  addedParticles[jscat].X = R2[1] + P2[1] * cc;
  addedParticles[jscat].Y = R2[2] + P2[2] * cc;
  addedParticles[jscat].Z = R2[3] + P2[3] * cc;

  addedParticles[jscat].X_traveled += sqrt( pow( P2[1] * cc, 2.0 ) + pow( P2[2] * cc, 2.0 ) + pow( P2[3] * cc, 2.0 ) );

  // at this point scattering take place, therefore set properties of last scattering point
  addedParticles[jscat].X_lastInt = addedParticles[jscat].X;
  addedParticles[jscat].Y_lastInt = addedParticles[jscat].Y;
  addedParticles[jscat].Z_lastInt = addedParticles[jscat].Z;
  addedParticles[jscat].T_lastInt = addedParticles[jscat].T;

  scatt22_obj.getMomenta22( PT2, typ, F1, F2 );
  scatt22_obj.setNewMomenta22( P1, P2, R1, R2, PT2 );

  double pt_out1 = sqrt( pow( P1[1], 2 ) + pow( P1[2], 2 ) );
  double pt_out2 = sqrt( pow( P2[1], 2 ) + pow( P2[2], 2 ) );

  //<<---------------------------------------------
  // set new properties for added particle
  if ( pt_out1 > pt_out2 )
  {
    addedParticles[jscat].FLAVOR = F1;
    addedParticles[jscat].PX = P1[1];
    addedParticles[jscat].PY = P1[2];
    addedParticles[jscat].PZ = P1[3];
    addedParticles[jscat].E = sqrt( pow( P1[1], 2 ) + pow( P1[2], 2 ) + pow( P1[3], 2 ) );
  }
  else
  {
    addedParticles[jscat].FLAVOR = F2;
    addedParticles[jscat].PX = P2[1];
    addedParticles[jscat].PY = P2[2];
    addedParticles[jscat].PZ = P2[3];
    addedParticles[jscat].E = sqrt( pow( P2[1], 2 ) + pow( P2[2], 2 ) + pow( P2[3], 2 ) );
  }
}



int offlineHeavyIonCollision::scatt32_utility( scattering32& scatt32_obj, std::list< int >& _cellMembersAdded, std::vector< int >& _allParticlesListAdded, std::vector< int >& _gluonListAdded, const int iscat, const int jscat, const int kscat, int& n32, const double picked_ratio, const double nexttime )
{
  double P1[4], P2[4], P3[4], R1[4], R2[4], R3[4];
  double Tmax, TT, u, phi, cc;
  int order;
  FLAVOR_TYPE F1, F2, F3;
  int typ;

  particles_atTimeNow[iscat].getMomentumArray( P1 );
  particles_atTimeNow[iscat].getCoordinateArray( R1 );
  F1 = particles_atTimeNow[iscat].FLAVOR;

  particles_atTimeNow[jscat].getMomentumArray( P2 );
  particles_atTimeNow[jscat].getCoordinateArray( R2 );
  F2 = particles_atTimeNow[jscat].FLAVOR;

  addedParticles[kscat].getMomentumArray( P3 );
  addedParticles[kscat].getCoordinateArray( R3 );
  F3 = addedParticles[kscat].FLAVOR;

  ncoll++;
  ncoll32++;
  addedParticles[kscat].coll_id = ncoll;

  if ( R1[0] > R2[0] )
    Tmax = R1[0];
  else
    Tmax = R2[0];
  if ( R3[0] > Tmax )
    Tmax = R3[0];
  TT = ( nexttime - Tmax ) * ran2() + Tmax;

  cc = ( TT - R3[0] ) / P3[0];
  addedParticles[kscat].T = TT;
  addedParticles[kscat].X = R3[1] + P3[1] * cc;
  addedParticles[kscat].Y = R3[2] + P3[2] * cc;
  addedParticles[kscat].Z = R3[3] + P3[3] * cc;

  scatt32_obj.getMomenta32( u, phi, picked_ratio, typ, F1, F2 );

  if ( scatt32_obj.getCollisionStatus() )
  {
    scatt32_obj.setNewMomenta32( P1, P2, u, phi );
    order = scatt32_obj.getOrder();

    double pt_out1 = sqrt( pow( P1[1], 2 ) + pow( P1[2], 2 ) );
    double pt_out2 = sqrt( pow( P2[1], 2 ) + pow( P2[2], 2 ) );

    if (( order == 1 ) || ( order == 3 ) )   //123  or  213
    {
      //mark absorbed gluon for removal from global particle list
      deadParticleList.push_back( kscat );
      addedParticles[kscat].dead = true;

      //erase absorbed gluon from the local list of gluons in this cell
      removeElementFromVector<int>( _gluonListAdded, kscat );
      //erase absorbed gluon from the local list of all particles in this cell
      removeElementFromVector<int>( _allParticlesListAdded, kscat );

      //erase absorbed gluon from the list of all particles in this cell
      _cellMembersAdded.remove( kscat );
    }
    else if (( order == 2 ) || ( order == 5 ) )  //132  or  312
    {
      if ( pt_out1 > pt_out2 )
      {
        addedParticles[kscat].FLAVOR = F1;
        addedParticles[kscat].PX = P1[1];
        addedParticles[kscat].PY = P1[2];
        addedParticles[kscat].PZ = P1[3];
        addedParticles[kscat].E = sqrt( P1[1] * P1[1] + P1[2] * P1[2] + P1[3] * P1[3] );
      }
      else
      {
        addedParticles[kscat].FLAVOR = F2;
        addedParticles[kscat].PX = P2[1];
        addedParticles[kscat].PY = P2[2];
        addedParticles[kscat].PZ = P2[3];
        addedParticles[kscat].E = sqrt( P2[1] * P2[1] + P2[2] * P2[2] + P2[3] * P2[3] );
      }
    }
    else //if((order == 4) || (order == 6))  //231 or 321
    {
      if ( pt_out1 > pt_out2 )
      {
        addedParticles[kscat].FLAVOR = F1;
        addedParticles[kscat].PX = P1[1];
        addedParticles[kscat].PY = P1[2];
        addedParticles[kscat].PZ = P1[3];
        addedParticles[kscat].E = sqrt( P1[1] * P1[1] + P1[2] * P1[2] + P1[3] * P1[3] );
      }
      else
      {
        addedParticles[kscat].FLAVOR = F2;
        addedParticles[kscat].PX = P2[1];
        addedParticles[kscat].PY = P2[2];
        addedParticles[kscat].PZ = P2[3];
        addedParticles[kscat].E = sqrt( P2[1] * P2[1] + P2[2] * P2[2] + P2[3] * P2[3] );
      }
    }

  }

  return order;
}





int offlineHeavyIonCollision::binomial( const int N, const int k ) const
{
  if ( N < k )
  {
    return 0;
  }
  else if ( k == 3 )
  {
    return ( N * ( N - 1 ) * ( N - 2 ) / 6 );
  }
  else if ( k == 2 )
  {
    return ( N * ( N - 1 ) / 2 );
  }
  else if ( k <= 0 )
  {
    return 0;
  }
  else if ( N == k )
  {
    return 1;
  }
  else
  {
    int nominator = 1;
    for ( int i = 0; i < k; i++ )
    {
      nominator *= ( N - i );
    }

    int denominator = 1;
    for ( int i = 1; i <= k; i++ )
    {
      denominator *= i;
    }

    return ( nominator / denominator );
  }

  return 0;
}




void offlineHeavyIonCollision::scatterEdgeParticles( std::list< int >& _offlineParticleList, std::list< int >& _addedParticleList, const double nexttime )
{
  int index;
  double c;

  list<int>::iterator iIt;

  for ( iIt = _addedParticleList.begin(); iIt != _addedParticleList.end(); )
  {
    index = *iIt;

    if ( !addedParticles[index].edge )
    {
      iIt = _addedParticleList.erase( iIt );
    }
    else
    {
      ++iIt;
    }
  }
  
  for ( iIt = _offlineParticleList.begin(); iIt != _offlineParticleList.end(); )
  {
    index = *iIt;

    if ( !particles_atTimeNow[index].edge )
    {
      iIt = _offlineParticleList.erase( iIt );
    }
    else
    {
      ++iIt;
    }
  }
  

//   int particleIndexOfFirstCollision = updateGeometricCollisionTimesA( _addedParticleList );
//   doGeometricCollisions( _addedParticleList, particleIndexOfFirstCollision );

  for ( iIt = _addedParticleList.begin(); iIt != _addedParticleList.end(); iIt++ )
  {
    index = *iIt;
    if ( addedParticles[index].T < nexttime )
    {
      c = ( nexttime - addedParticles[index].T ) / addedParticles[index].E;
      addedParticles[index].T = nexttime;
      addedParticles[index].X = addedParticles[index].X + addedParticles[index].PX * c;
      addedParticles[index].Y = addedParticles[index].Y + addedParticles[index].PY * c;
      addedParticles[index].Z = addedParticles[index].Z + addedParticles[index].PZ * c;
    }
  }
}



/**
* This routine removes particlre stored in deadParticleList (global to this translation unit via unnamed namespace)
*/
void offlineHeavyIonCollision::removeDeadParticles( analysis& _aa )
{
  double pt_new;
  int lastIndex = -1;
  deadParticleList.sort();
  while ( addedParticles.size() != 0 && addedParticles.back().dead )
  {
    addedParticles.pop_back();
    --numberAdded;
    deadParticleList.pop_back();
  }

  list<int>::iterator it;
  for ( it = deadParticleList.begin(); it != deadParticleList.end(); it++ )
  {
    lastIndex = addedParticles.size() - 1;

    addedParticles[( *it )] = addedParticles.back();
    if ( addedParticles.back().edge )
    {
      edgeCellAdded.remove( lastIndex );
      edgeCellAdded.remove( *it ); // just to be sure it's not added twice
      edgeCellAdded.push_back( *it );

//       list<int>::iterator jIt;
//       for ( jIt = edgeCellAdded.begin(); jIt != edgeCellAdded.end() ; jIt++ )
//       {
//         if ( addedParticles[*jIt].collisionPartner == lastIndex )
//         {
//           addedParticles[*jIt].collisionPartner = *it;
//         }
//       }
    }
    addedParticles.pop_back();
    --numberAdded;

    while ( addedParticles.back().dead )
    {
      addedParticles.pop_back();
      --numberAdded;
      deadParticleList.pop_back();
    }
  }

}



//provides iterative calculation of MFP for high-pt particles, returned lambda and start value in GeV^-1
double offlineHeavyIonCollision::iterateMFP( std::vector< int >& _allParticlesList, std::vector< int >& _gluonList, const int jetID, const double dt, const double dv )
{
  double P1[4], P2[4], P3[4];
  int iscat, jscat, kscat;
  int n32 = 0, n22 = 0, n23 = 0;
  double lambda_scaled, s;
  double probab22 = 0, probab23 = 0, probab32 = 0;
  double cs22, cs23, I32, I23;
  double R22, R23, R32;
  double as, Vrel, md2g, md2q;
  double lambda, lambdaAvr;
  int iter = 0;
  deque<double> lambdaArray;
  double betaDistEntry;
  
  int initialStateIndex = -1;
  FLAVOR_TYPE F1, F2, F3;
  
  const double gG = 2 * ( Ncolor * Ncolor - 1 );
  const double epsilon = 0.01;
  const int nIterationsMax = 8;
  bool converged = false;
  
  
  //----------------------------- to which ring does the jet belong? ---------------------------
  int nc = rings.getIndex( addedParticles[jetID] );
  //--------------------------------------------------------------------------------------------
  
  double vx = rings[nc].getAveraged_v_x();
  double vy = rings[nc].getAveraged_v_y();
  double vz = rings[nc].getAveraged_v_z();
  
  
  F1 = addedParticles[jetID].FLAVOR;
  P1[0] = addedParticles[jetID].E;
  P1[1] = addedParticles[jetID].PX;
  P1[2] = addedParticles[jetID].PY;
  P1[3] = addedParticles[jetID].PZ;
  
  
  list<int>::const_iterator iIt;
  list<int>::const_iterator jIt;
  
  
  const int nTotal = _allParticlesList.size();
  const int nTotalAdded = 1;
  const int nGluons = _gluonList.size();
  int nGluonsAdded = 0;
  if ( addedParticles[jetID].FLAVOR == gluon )
  {
    nGluonsAdded = 1; 
  }
  const int nAllQuarks = nTotal - nGluons;
  const int nAllQuarksAdded = nTotalAdded - nGluonsAdded;
  
  const int allTriplets = ( binomial( nTotal, 2 ) * nTotalAdded )  - ( binomial( nAllQuarks, 2 ) * nAllQuarksAdded );
  
  
  //------------------------------ determine initial value of lambda ---------------------------
  double csgg = 0;
  const double small = 1.0e-4;
  double jetRate = ( addedParticles[jetID].rate + addedParticles[jetID].ratev ) / 2.0;
  if ( jetRate > small )
  {
    lambda = 1 / jetRate;
  }
  else
  {
    int selected = 0;
    do  // pick random particle (that is neihter dead nor the jet) from the _particleList
    {
      selected = static_cast<int>( ran2() * _allParticlesList.size() );
      if ( selected >= _allParticlesList.size() )
      {
        continue;
      }      
      iscat = _allParticlesList[selected];
    }
    while ( particles_atTimeNow[iscat].dead );
    
    P2[0] = particles_atTimeNow[iscat].E;
    P2[1] = particles_atTimeNow[iscat].PX;
    P2[2] = particles_atTimeNow[iscat].PY;
    P2[3] = particles_atTimeNow[iscat].PZ;
    
    as = alpha_s( s );
    md2g = as * ( addedParticles[jetID].md2g + particles_atTimeNow[iscat].md2g ) / 2.0;
    md2q = as * ( addedParticles[jetID].md2q + particles_atTimeNow[iscat].md2q ) / 2.0;

    // HACK for N_f = 0 background
//     md2g = as * ( addedParticles[jetID].md2g + particles_atTimeNow[iscat].md2g ) / 2.0 * 2;
//     md2q = md2g * 2.0 / 9.0;
    
    s = pow(( P1[0] + P2[0] ), 2 ) - pow(( P1[1] + P2[1] ), 2 ) - pow(( P1[2] + P2[2] ), 2 ) - pow(( P1[3] + P2[3] ), 2 );
    
    xsection_gg_gg csObj( s, md2g, md2q );
    csgg = csObj.totalCrossSection();
    lambda = ( dv * rings[nc].getGamma() * testpartcl   ) / ( pow( 0.197, 3.0 ) * _allParticlesList.size() * csgg );  
  }
  //--------------------------------------------------------------------------------------------
  
  scattering22 scatt22_object;
  scattering32 scatt32_object;
  scattering23 scatt23_object( &theI23 );
  
  do
  {    
    //------------------------ 2<->2 & 2->3-----------------------
    for ( int m1 = 0; m1 < _allParticlesList.size(); m1++ )
    {
      jscat = _allParticlesList[m1];
      if ( !particles_atTimeNow[jscat].dead )
      {
        F2 = particles_atTimeNow[jscat].FLAVOR;
        P2[0] = particles_atTimeNow[jscat].E;
        P2[1] = particles_atTimeNow[jscat].PX;
        P2[2] = particles_atTimeNow[jscat].PY;
        P2[3] = particles_atTimeNow[jscat].PZ;
        
        s = pow(( P1[0] + P2[0] ), 2 ) - pow(( P1[1] + P2[1] ), 2 ) - pow(( P1[2] + P2[2] ), 2 ) - pow(( P1[3] + P2[3] ), 2 );
        
        if ( s > 1.1*lambda2 )
        {
          n22++;
          Vrel = s / ( 2.0 * P1[0] * P2[0] );
          as = alpha_s( s );
          
          md2g = as * ( addedParticles[jetID].md2g + particles_atTimeNow[jscat].md2g ) / 2.0;
          md2q = as * ( addedParticles[jetID].md2q + particles_atTimeNow[jscat].md2q ) / 2.0;
          // HACK for N_f = 0 background
//           md2g = as * ( addedParticles[jetID].md2g + particles_atTimeNow[jscat].md2g ) / 2.0 * 2;          
//           md2q = md2g * 2.0 / 9.0;
          
          scatt22_object.setParameter( P1, P2, F1, F2, s, md2g, md2q );
          cs22 = scatt22_object.getXSection22();
          probab22 += pow( 0.197, 2.0 ) * cs22 * Vrel * dt / ( dv * testpartcl );
          
          n23++;
          lambda_scaled = lambda * sqrt( s );
          
          betaDistEntry = scatt23_object.setParameter( vx, vy, vz, P1, P2, F1, F2 , sqrt( s ), md2g / s, lambda_scaled, _gluonList.size() );
          I23 = scatt23_object.getIntegral23( initialStateIndex );
          cs23 = Ncolor * pow( as, 3 ) / s * I23;    //1/GeV^2
          probab23 += pow( 0.197, 2.0 ) * cs23 * Vrel * dt / ( dv * testpartcl );
        }
      }
    }
    R22 = probab22 / dt * rings[nc].getGamma();
    R23 = probab23 / dt * rings[nc].getGamma();
    //-------------------------------------------------------------
    
    const int consideredTriplets = 5;
    double scaleForSelectedTriplets = 1;
    if ( allTriplets > 10 )
    {      
      int m2, m3;
      scaleForSelectedTriplets = static_cast<double>( allTriplets ) / static_cast<double>( consideredTriplets );
      
      for ( int i = 0; i < consideredTriplets && _gluonList.size() > 0; i++ )
      {
        do
        {
          do
          {
            m2 = int ( _allParticlesList.size() * ran2() );
            if ( m2 == _allParticlesList.size() )
            {
              m2 = _allParticlesList.size() - 1;
            }
            iscat = _allParticlesList[m2];
          }
          while ( particles_atTimeNow[iscat].dead );
          F2 = particles_atTimeNow[iscat].FLAVOR;
          
          do
          {
            m3 = int ( _allParticlesList.size() * ran2() );
            if ( m3 == _allParticlesList.size() )
            {
              m3 = _allParticlesList.size() - 1;
            }
            jscat = _allParticlesList[m3];
          }
          while (( m3 == m2 ) || ( particles_atTimeNow[jscat].dead ) );
          F3 = particles_atTimeNow[jscat].FLAVOR;          
        }
        while ( !( F1 == gluon || F2 == gluon || F3 == gluon ) );
        
        P2[0] = particles_atTimeNow[iscat].E;
        P2[1] = particles_atTimeNow[iscat].PX;
        P2[2] = particles_atTimeNow[iscat].PY;
        P2[3] = particles_atTimeNow[iscat].PZ;
        
        P3[0] = particles_atTimeNow[jscat].E;
        P3[1] = particles_atTimeNow[jscat].PX;
        P3[2] = particles_atTimeNow[jscat].PY;
        P3[3] = particles_atTimeNow[jscat].PZ;
        
        s = pow(( P1[0] + P2[0] + P3[0] ), 2 ) - pow(( P1[1] + P2[1] + P3[1] ), 2 ) - pow(( P1[2] + P2[2] + P3[2] ), 2 ) - pow(( P1[3] + P2[3] + P3[3] ), 2 );
        
        if ( s > 1.1*lambda2 )
        {
          n32++;
          lambda_scaled = lambda * sqrt( s );
          
          as = alpha_s( s );
          
          md2g = as * ( particles_atTimeNow[iscat].md2g + particles_atTimeNow[jscat].md2g + addedParticles[jetID].md2g ) / 3.0;
          md2g = as * ( particles_atTimeNow[iscat].md2q + particles_atTimeNow[jscat].md2q + addedParticles[jetID].md2q ) / 3.0;
          // HACK for N_f = 0 background
//           md2g = as * ( particles_atTimeNow[iscat].md2g + particles_atTimeNow[jscat].md2g + addedParticles[jetID].md2g ) / 3.0 * 2;
//           md2q = md2g * 2.0 / 9.0;
          
          // create scattering32 object for the given 3 particles
          betaDistEntry = scatt32_object.setParameter( vx, vy, vz, P1, P2, P3, F1, F2, F3, sqrt( s ), md2g / s, lambda_scaled, _gluonList.size() );
          I32 = scatt32_object.getIntegral32_fast();                        // get the integral I32 for the given 3 particles
          
          probab32 += pow( 0.197, 5.0 ) * 9.0 * M_PI * Ncolor * pow( as, 3.0 ) * dt * I32 / ( gG * pow( dv, 2.0 ) * pow( testpartcl, 2.0 ) * s * P1[0] * P2[0] * P3[0] );
        }
      }
    }
    else
    {
      //---------------------------- 3->2 ---------------------------
      for ( int m1 = 0; m1 < _allParticlesList.size() - 1; m1++ )
      {
        iscat = _allParticlesList[m1];
        if ( !particles_atTimeNow[iscat].dead )
        {   
          for ( int m2 = m1 + 1; m2 < _allParticlesList.size(); m2++ )
          {
            jscat = _allParticlesList[m2];
            if ( !particles_atTimeNow[jscat].dead )
            {
              F2 = particles_atTimeNow[iscat].FLAVOR;
              P2[0] = particles_atTimeNow[iscat].E;
              P2[1] = particles_atTimeNow[iscat].PX;
              P2[2] = particles_atTimeNow[iscat].PY;
              P2[3] = particles_atTimeNow[iscat].PZ;
              
              F3 = particles_atTimeNow[jscat].FLAVOR;
              P3[0] = particles_atTimeNow[jscat].E;
              P3[1] = particles_atTimeNow[jscat].PX;
              P3[2] = particles_atTimeNow[jscat].PY;
              P3[3] = particles_atTimeNow[jscat].PZ;
              
              s = pow(( P1[0] + P2[0] + P3[0] ), 2 ) - pow(( P1[1] + P2[1] + P3[1] ), 2 ) - pow(( P1[2] + P2[2] + P3[2] ), 2 ) - pow(( P1[3] + P2[3] + P3[3] ), 2 );
              
              if ( s > 1.1*lambda2 )
              {
                n32++;
                lambda_scaled = lambda * sqrt( s );
                
                as = alpha_s( s );
                
                md2g = as * ( particles_atTimeNow[iscat].md2g + particles_atTimeNow[jscat].md2g + addedParticles[jetID].md2g ) / 3.0;
                md2g = as * ( particles_atTimeNow[iscat].md2q + particles_atTimeNow[jscat].md2q + addedParticles[jetID].md2q ) / 3.0;
                // HACK for N_f = 0 background
//                 md2g = as * ( particles_atTimeNow[iscat].md2g + particles_atTimeNow[jscat].md2g + addedParticles[jetID].md2g ) / 3.0 * 2;
//                 md2q = md2g * 2.0 / 9.0;
                
                // create scattering32 object for the given 3 particles
                betaDistEntry = scatt32_object.setParameter( vx, vy, vz, P1, P2, P3, F1, F2, F3, sqrt( s ), md2g / s, lambda_scaled, _gluonList.size() );
                I32 = scatt32_object.getIntegral32_fast();                        // get the integral I32 for the given 3 particles
                
                probab32 += pow( 0.197, 5.0 ) * 9.0 * M_PI * Ncolor * pow( as, 3.0 ) * dt * I32 / ( gG * pow( dv, 2.0 ) * pow( testpartcl, 2.0 ) * s * P1[0] * P2[0] * P3[0] );
              }
            }
          }
        }
      }
    }
    R32 = scaleForSelectedTriplets * probab32 / dt * rings[nc].getGamma();  // fm^-1
    //-------------------------------------------------------------
    
    lambda = 1 / ( R22 + R23 + R32 ); //fm
    lambda = lambda / 0.197;  //GeV^-1
    
    if ( lambdaArray.size() < 4 )
    {
      lambdaArray.push_back( lambda );
    }
    else
    {
      lambdaArray.push_back( lambda );
      lambdaArray.pop_front();
    }
    
    if ( lambdaArray.size() == 4 )
    {
      lambdaAvr = 0;
      for ( int m = 0; m < lambdaArray.size(); m++ )
      {
        lambdaAvr += lambdaArray[m];
      }
      lambdaAvr = lambdaAvr / lambdaArray.size();
      
      if (( fabs( lambdaArray[3] - lambdaAvr ) / lambdaAvr < epsilon ) &&
        ( fabs( lambdaArray[2] - lambdaAvr ) / lambdaAvr < epsilon ) &&
        ( fabs( lambdaArray[1] - lambdaAvr ) / lambdaAvr < epsilon ) &&
        ( fabs( lambdaArray[0] - lambdaAvr ) / lambdaAvr < epsilon ) )
      {
        converged = true;
      }
      else
      {
        converged = false;
      }
    }
    else
    {
      converged = false;
    }
    
    ++iter;
    probab22 = probab23 = probab32 = 0;
  }
  while ( !converged && iter < nIterationsMax );
  
  lambdaAvr = 0;
  for ( int m = 0; m < lambdaArray.size(); m++ )
  {
    lambdaAvr += lambdaArray[m];
  }
  lambdaAvr = lambdaAvr / lambdaArray.size();
  
  return lambdaAvr; //GeV^-1
}


// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;
