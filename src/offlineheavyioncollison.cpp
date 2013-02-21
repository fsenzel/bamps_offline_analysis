//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------
//---------------------------------------------


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
#include "coupling.h"
#include "scattering22.h"
#include "scattering23.h"
#include "prefactors23.h"
#include "random.h"
#include "binary_cross_sections.h"
#include "offlineoutput.h"
#include "FPT_compare.h"
#include "hadronization_hq.h"
#include "lorentz.h"

#include "configBAMPS.h"
#ifdef Pythia_FOUND
#include "mesonDecay.h"
#endif

using namespace ns_casc;
using namespace std;


extern int IX, IY, IZ;


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

namespace ns_heavy_quarks
{
  int jpsi_dissociation = 0;
  int jpsi_dissociation_from_temperature = 0;
  int jpsicreation = 0;
  int charmAnnihil = 0;
}




offlineHeavyIonCollision::offlineHeavyIonCollision( config* const _config, offlineOutputInterface* const _offlineInterface, analysis* const _analysis ) :
    theConfig( _config ), stoptime_last( 0 ), stoptime( 5.0 ), currentNumber( 0 ), numberEvolvingParticles( _config->getN_init() ),
    rings( _config->getRingNumber(), _config->getCentralRingRadius(), _config->getDeltaR() ),
    testpartcl( _config->getTestparticles() ),
    offlineInterface( _offlineInterface ),
    theI23( _config->doScattering_23() ), theAnalysis( _analysis )
{
  theI22.configure( theConfig->isCouplingRunning(), Particle::N_light_flavor, Particle::N_heavy_flavor, Particle::Mcharm, Particle::Mbottom );
}



offlineHeavyIonCollision::~offlineHeavyIonCollision()
{

}



void offlineHeavyIonCollision::initialize()
{
  if ( theConfig->getNumberOfParticlesToAdd() != 0 )
  {
    additionalParticlesDistribution addedStuff( theConfig, theConfig->getInitialStateType() );
    addedStuff.populateParticleVector( addedParticles, WoodSaxonParameter );
    
    int N_light_flav_added = theConfig->getNlightFlavorsAdded();
    int N_heavy_flav_added = theConfig->getNheavyFlavorsAdded();

    int Nbefore = addedParticles.size();
    for(int j = 0; j < addedParticles.size(); j++ )
    {
      if( ( addedParticles[j].FLAVOR > 2 * N_light_flav_added ) && 
          ( ( addedParticles[j].FLAVOR <= 2 * Particle::max_N_light_flavor ) || 
            ( addedParticles[j].FLAVOR > 2 * ( Particle::max_N_light_flavor + N_heavy_flav_added ) ) ) &&
          !( addedParticles[j].FLAVOR >= 50 && addedParticles[j].FLAVOR < 50 + Particle::N_psi_states )
        )
      {
        // delete last particle if also not active otherwise switch position with particle to be deleted
        while( ( addedParticles.back().FLAVOR > 2 * N_light_flav_added ) && 
              ( ( addedParticles.back().FLAVOR <= 2 * Particle::max_N_light_flavor ) || 
                ( addedParticles.back().FLAVOR > 2 * ( Particle::max_N_light_flavor + N_heavy_flav_added ) ) ) &&
              !( addedParticles.back().FLAVOR >= 50 && addedParticles.back().FLAVOR < 50 + Particle::N_psi_states )&& 
              ( j != addedParticles.size() - 1 ) ) // if particle j is the last particle in the particle list it is deleted here and the then last in the list below as well, which is not correct.
        {
          addedParticles.pop_back();
        }
        addedParticles[j] = addedParticles.back();
        addedParticles.pop_back();
      }
    }
    cout << "#### " << addedParticles.size() << " out of " << Nbefore << " particles kept for simulation ( N_f = N_f_light_quarks + N_f_heavy_quarks = " << N_light_flav_added << " + " <<   N_heavy_flav_added << " )." << endl;
    
    // List particle numbers for all flavors
    cout << "==========================" << endl;
    cout << "Added particles:" << endl;
    cout << "# of all=" << addedParticles.size() << endl;
    int fl_sum = 0;
    cout << "Flavor       Number    Number wo testparticles per Event" << endl;
    for( int i = 0; i <= 10; i++ )
    {
      for( int j = 0; j < addedParticles.size(); j++ )
        if(addedParticles[j].FLAVOR == i)
          fl_sum++;
      cout.width(5);
      cout << i;
      cout.width(12);
      cout << fl_sum;
      cout.width(16);
      cout << double( fl_sum ) / testpartcl / theConfig->getNaddedEvents() << endl;
      fl_sum = 0;
    }
    for( int j = 0; j < addedParticles.size(); j++ )
      if(addedParticles[j].FLAVOR == jpsi)
        fl_sum++;
    cout.width(5);
    cout << "Jpsi";
    cout.width(12);
    cout << fl_sum;
    cout.width(16);
    cout << double( fl_sum ) / testpartcl / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles() << endl;
    fl_sum = 0;
    double sum = 0.0;
    for( int j = 0; j < addedParticles.size(); j++ )
      sum += addedParticles[j].E;
    cout << "E(init) = " << sum << endl;
    cout << "==========================" << endl;

    addedStuff.prepareParticles( addedParticles );
  }
  else
  {
    addedParticles.clear();
    cout << "#### 0 particles added" << endl;    
  }
}



void offlineHeavyIonCollision::mainFramework()
{
  double dt_cascade = 0;
  double dt_backup = 0;
  double nexttime;
  int ncoll_backup = 0;
  int ncoll22_backup = 0;
  int ncoll23_backup = 0;
  int ncoll32_backup = 0;
  int ncolle_backup = 0;
  int charmAnnihil_backup = 0;
  int jpsicreation_backup = 0;
  int jpsi_dissociation_from_temperature_backup = 0;
  int jpsi_dissociation_backup = 0;
  
  list<int> edgeCellCopy, edgeCellAddedCopy;

  gG = 2 * ( pow( static_cast<double>( Ncolor ) , 2 ) - 1 );
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
  cout << "scale time steps dt by factor " << factor_dt << endl;
  
  if( theConfig->isHadronizationHQ() )
  {
    hadronization_hq ppHadronization_hq;
    ppHadronization_hq.heavyQuarkFragmentation();
  }
  if( theConfig->isMesonDecay() )
  {
#ifdef Pythia_FOUND
    mesonDecay ppMesonDecay( theConfig->getNumberElectronStat(), theConfig->isMuonsInsteadOfElectrons(), theConfig->isStudyNonPromptJpsiInsteadOfElectrons() );
    ppMesonDecay.decayToElectronsPythia();
#else
    string errMsg( "Could not perform decay of heavy mesons to electron because PYTHIA was not found." );
    throw eHIC_error( errMsg );
#endif
  }

  theAnalysis->initialOutput();
  if ( theConfig->doOutput_movieOutputJets() )
  {
    theAnalysis->movieOutput( 0, jumpMovieSteps );
  }
  if ( theConfig->doOutput_movieOutputBackground() )
  {
    theAnalysis->movieOutputMedium( 0, jumpMovieSteps );
  }
  theAnalysis->collectPtDataInitial();
  theAnalysis->collectYDataInitial();
  theAnalysis->collectEtDataInitial();
  simulationTime = theConfig->getTimefirst(); //fm/c

  while( simulationTime >= theAnalysis->tstep[nn_ana] )
  {
    nn_ana++;
  }
  while( simulationTime >= theAnalysis->tstep_movie[nn_ana_movie] )
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
      theAnalysis->addJetEvent_initial( i );
    }
  }

  int n_dt = 0;
  double dt_sum = 0.0;
  int n_again = 0;

  do
  {
    // remove init tag from particles after their formation
    for ( int i = 0; i < addedParticles.size(); i++ )
    {
      if ( addedParticles[i].T <= simulationTime )
      {
        addedParticles[i].init = false;
      }
    }

    // evolution of the medium to present time
    double dt_cascade_from_data = evolveMedium( simulationTime, endOfDataFiles );
    
    // specify time step
    if ( theConfig->useFixed_dt() )
    {
      dt = theConfig->getFixed_dt();
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
//     cout << "# time = " << simulationTime << "    dt = " << dt << endl;

    p23_collected_gluon = 0;
    n23_collected_gluon = 0;
    p23_collected_quark = 0;
    n23_collected_quark = 0;
    lambdaJet_gluon = 0;
    lambdaJet_quark = 0;
      
    if ( doMovieStepMedium && theConfig->doOutput_movieOutputBackground() )
    {
      theAnalysis->movieOutputMedium( nn_ana_movie - 1, jumpMovieSteps );
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
    charmAnnihil_backup = ns_heavy_quarks::charmAnnihil;
    jpsicreation_backup = ns_heavy_quarks::jpsicreation;
    jpsi_dissociation_from_temperature_backup = ns_heavy_quarks::jpsi_dissociation_from_temperature;
    jpsi_dissociation_backup = ns_heavy_quarks::jpsi_dissociation;
    //--------------------------

    nexttime = simulationTime + dt;

    if( nexttime >= theAnalysis->tstep_movie[nn_ana_movie] || nexttime >= theAnalysis->tstep[nn_ana] )
    {
      dt_backup = dt;
      
      if( nexttime >= theAnalysis->tstep[nn_ana] )  // ask if it is time for analysis
      {
        nexttime = theAnalysis->tstep[nn_ana];
        dt = nexttime - simulationTime;
        doAnalysisStep = true;
        cout << "profile " << nexttime << endl;
      }
      if( nexttime >= theAnalysis->tstep_movie[nn_ana_movie] )  // ask if it is time for movie output
      {
        nexttime = theAnalysis->tstep_movie[nn_ana_movie];
        dt = nexttime - simulationTime;
        doMovieStep = true;
        if ( theConfig->doOutput_movieOutputJets() || theConfig->doOutput_movieOutputBackground() )
        {
          cout << "** movie: " << nexttime << endl;
        }
      }
      
      if ( doAnalysisStep && doMovieStep )
      {
        if( !FPT_COMP_E( theAnalysis->tstep[nn_ana], theAnalysis->tstep_movie[nn_ana_movie] ) )
        {
          string errMsg( "time steps for movie output and general analysis output do not match" );
          throw eHIC_error( errMsg );
        }
      }      
    }
    
    cell_ID( nexttime );

    // collide added particles with gluonic medium
    deadParticleList.clear();
    scattering( nexttime, again );

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
      cells = cellsCopy;
      cellsAdded = cellsAddedCopy;
      edgeCell = edgeCellCopy;
      edgeCellAdded = edgeCellAddedCopy;
      
      ncoll = ncoll_backup;
      ncoll22 = ncoll22_backup;
      ncoll23 = ncoll23_backup;
      ncoll32 = ncoll32_backup;
      ncolle = ncolle_backup;
      ns_heavy_quarks::charmAnnihil = charmAnnihil_backup;
      ns_heavy_quarks::jpsicreation = jpsicreation_backup;
      ns_heavy_quarks::jpsi_dissociation_from_temperature = jpsi_dissociation_from_temperature_backup;
      ns_heavy_quarks::jpsi_dissociation = jpsi_dissociation_backup;

      nexttime = simulationTime + dt;

      cell_ID( nexttime );

      deadParticleList.clear();
      scattering( nexttime, again );
    }

    scatterEdgeParticles( edgeCell, edgeCellAdded, nexttime );

    removeDeadParticles();
    if( theConfig->doOutput_progressLog() )
    {
      theAnalysis->registerProgressInformationForOutput( simulationTime, dt, addedParticles.size(), particles_atTimeNow.size(), ncoll, ncoll22, ncoll23, ncoll32 );
    }
    
    if ( doAnalysisStep )
    {
      theAnalysis->intermediateOutput( nn_ana );
      theAnalysis->collectPtData( nn_ana );
      theAnalysis->collectYData( nn_ana );
      theAnalysis->collectEtData( nn_ana );
      nn_ana++;
      doAnalysisStep = false;
      dt = dt_backup;
    }
    
    if ( doMovieStep )
    {
      theAnalysis->mfpJetsOutput( nn_ana_movie, jumpMovieSteps );
      
      if ( theConfig->doOutput_movieOutputJets() )
      {
        theAnalysis->movieOutput( nn_ana_movie, jumpMovieSteps );
      }
      nn_ana_movie++;
      doMovieStep = false;
      doMovieStepMedium = true;
      dt = dt_backup;
    }
    theAnalysis->printCentralDensities( simulationTime );
    
    
//     // just error checking if masses and flavors are correct
//     for ( int j = 0; j < addedParticles.size(); j++ )
//     {
//       double E_check = sqrt( pow( addedParticles[j].PX, 2 ) + pow( addedParticles[j].PY, 2 ) + pow( addedParticles[j].PZ, 2 ) + pow( addedParticles[j].m, 2 ) );
//       if( !FPT_COMP_E( addedParticles[j].E, E_check ) )
//         cout << "Error! Particle " << j << " does not fulfill E^2 = p^2 + m^2." << endl;
//       
//       if( !FPT_COMP_E( addedParticles[j].m, Particle::getMass( addedParticles[j].FLAVOR ) ) )
//         cout << "Error! Particle " << j << " with flavor " << addedParticles[j].FLAVOR << " has wrong mass: " << addedParticles[j].m << "  " << Particle::getMass( addedParticles[j].FLAVOR )  << endl;
//       
//       if( ( addedParticles[j].FLAVOR > 2 * theConfig->getNlightFlavorsAdded() ) && 
//           ( ( addedParticles[j].FLAVOR <= 2 * 3 ) || 
//             ( addedParticles[j].FLAVOR > 2 * ( 3 + theConfig->getNheavyFlavorsAdded() ) ) ) &&
//           !( addedParticles[j].FLAVOR >= 50 && addedParticles[j].FLAVOR < 50 + Particle::N_psi_states )
//         )
//         cout << "Error! Particle " << j << " with flavor " << addedParticles[j].FLAVOR << " should not exist: Nf = " << theConfig->getNlightFlavorsAdded() << " + " <<   theConfig->getNheavyFlavorsAdded() <<  endl;
//     }

    // analyse timesteps
    dt_sum += dt;
    n_dt++;
    
    simulationTime = nexttime;
  }
  while ( simulationTime < stoptime && !endOfDataFiles );//fm/c

  if( theConfig->isHadronizationHQ() )
  {
    hadronization_hq theHadronization_hq;
    theHadronization_hq.heavyQuarkFragmentation();
  }
  if( theConfig->isMesonDecay() )
  {
#ifdef Pythia_FOUND
    mesonDecay theMesonDecay( theConfig->getNumberElectronStat(), theConfig->isMuonsInsteadOfElectrons(), theConfig->isStudyNonPromptJpsiInsteadOfElectrons() );
    theMesonDecay.decayToElectronsPythia();
#else
    string errMsg( "Could not perform decay of heavy mesons to electron because PYTHIA was not found." );
    throw eHIC_error( errMsg );
#endif
  }

  theAnalysis->finalOutput( stoptime );
  theAnalysis->addJetEvents_final();
  
  // List particle numbers for all flavors
  cout << "==========================" << endl;
  cout << "Added particles:" << endl;
  cout << "# of all=" << addedParticles.size() << endl;
  int fl_sum = 0;
  cout << "Flavor       Number    Number wo testparticles per Event" << endl;
  for( int i = 0; i <= 10; i++ )
  {
    for( int j = 0; j < addedParticles.size(); j++ )
      if(addedParticles[j].FLAVOR == i)
        fl_sum++;
    cout.width(5);
    cout << i;
    cout.width(12);
    cout << fl_sum;
    cout.width(16);
    cout << double( fl_sum ) / testpartcl / theConfig->getNaddedEvents() << endl;
    fl_sum = 0;
  }
  for( int j = 0; j < addedParticles.size(); j++ )
    if(addedParticles[j].FLAVOR == jpsi)
      fl_sum++;
  cout.width(5);
  cout << "Jpsi";
  cout.width(12);
  cout << fl_sum;
  cout.width(16);
  cout << double( fl_sum ) / testpartcl / theConfig->getNaddedEvents() / theConfig->getJpsiTestparticles() << endl;
  fl_sum = 0;
  double sum = 0.0;
  for( int j = 0; j < addedParticles.size(); j++ )
    sum += addedParticles[j].E;
  cout << "E(end) = " << sum << endl;
  cout << "==========================" << endl;

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
  for ( int k = 0; k < particles_atTimeNow.size(); k++ )
  {
    particlesEvolving[k].init = particles_atTimeNow[k].init;
    particlesEvolving[k].edge = particles_atTimeNow[k].edge;
    particlesEvolving[k].free = particles_atTimeNow[k].free;
  }

  if ( evolveToTime <= stoptime_last )
  {
    for ( int i = 0; i < theConfig->getN_init(); i++ )
    {
      particlesEvolving[i].init = true;

      particlesEvolving[i].T = particles_init[i].T;
      particlesEvolving[i].X = particles_init[i].X;
      particlesEvolving[i].Y = particles_init[i].Y;
      particlesEvolving[i].Z = particles_init[i].Z;

      particlesEvolving[i].PXold = particlesEvolving[i].PX = particles_init[i].PX;
      particlesEvolving[i].PYold = particlesEvolving[i].PY = particles_init[i].PY;
      particlesEvolving[i].PZold = particlesEvolving[i].PZ = particles_init[i].PZ;
      particlesEvolving[i].Eold = particlesEvolving[i].E = particles_init[i].E;
    }

    numberEvolvingParticles = theConfig->getN_init();
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
        particlesEvolving[iscat].init = particlesEvolving[jscat].init = false;

        c = ( time - particlesEvolving[iscat].T ) / particlesEvolving[iscat].E;
        particlesEvolving[iscat].T = time;
        particlesEvolving[iscat].X = particlesEvolving[iscat].X + particlesEvolving[iscat].PX * c;
        particlesEvolving[iscat].Y = particlesEvolving[iscat].Y + particlesEvolving[iscat].PY * c;
        particlesEvolving[iscat].Z = particlesEvolving[iscat].Z + particlesEvolving[iscat].PZ * c;

        c = ( time - particlesEvolving[jscat].T ) / particlesEvolving[jscat].E;
        particlesEvolving[jscat].T = time;
        particlesEvolving[jscat].X = particlesEvolving[jscat].X + particlesEvolving[jscat].PX * c;
        particlesEvolving[jscat].Y = particlesEvolving[jscat].Y + particlesEvolving[jscat].PY * c;
        particlesEvolving[jscat].Z = particlesEvolving[jscat].Z + particlesEvolving[jscat].PZ * c;

        if ( time < particlesEvolving[iscat].T )
        {
          cout << "back22_i" << endl;
          int zz;
          cin >> zz;
        }
        if ( time < particlesEvolving[jscat].T )
        {
          cout << "back22_j" << endl;
          int zz;
          cin >> zz;
        }

        particlesEvolving[iscat].FLAVOR = static_cast<FLAVOR_TYPE>( F1 );
        particlesEvolving[iscat].PX = pxi;
        particlesEvolving[iscat].PY = pyi;
        particlesEvolving[iscat].PZ = pzi;
        particlesEvolving[iscat].E = sqrt( pow( pxi, 2 ) + pow( pyi, 2 ) + pow( pzi, 2 ) );

        particlesEvolving[jscat].FLAVOR = static_cast<FLAVOR_TYPE>( F2 );
        particlesEvolving[jscat].PX = pxj;
        particlesEvolving[jscat].PY = pyj;
        particlesEvolving[jscat].PZ = pzj;
        particlesEvolving[jscat].E = sqrt( pow( pxj, 2 ) + pow( pyj, 2 ) + pow( pzj, 2 ) );
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
      particlesEvolving[kscat].X = ptrInteraction23->newx;
      particlesEvolving[kscat].Y = ptrInteraction23->newy;
      particlesEvolving[kscat].Z = ptrInteraction23->newz;
      particlesEvolving[kscat].PX = ptrInteraction23->newpx;
      particlesEvolving[kscat].PY = ptrInteraction23->newpy;
      particlesEvolving[kscat].PZ = ptrInteraction23->newpz;
      particlesEvolving[kscat].unique_id = Particle::unique_id_counter;
      Particle::unique_id_counter++;

      if ( time <= evolveToTime + 1.0e-6 )
      {
        particlesEvolving[kscat].T = time;
        particlesEvolving[kscat].FLAVOR = static_cast<FLAVOR_TYPE>( F3 );
        particlesEvolving[kscat].E = sqrt( pow( particlesEvolving[kscat].PX, 2 ) + pow( particlesEvolving[kscat].PY, 2 ) + pow( particlesEvolving[kscat].PZ, 2 ) );
        particlesEvolving[kscat].free = false;

        particlesEvolving[iscat].init = particlesEvolving[jscat].init = particlesEvolving[kscat].init = false;

        c = ( time - particlesEvolving[iscat].T ) / particlesEvolving[iscat].E;
        particlesEvolving[iscat].T = time;
        particlesEvolving[iscat].X = particlesEvolving[iscat].X + particlesEvolving[iscat].PX * c;
        particlesEvolving[iscat].Y = particlesEvolving[iscat].Y + particlesEvolving[iscat].PY * c;
        particlesEvolving[iscat].Z = particlesEvolving[iscat].Z + particlesEvolving[iscat].PZ * c;

        c = ( time - particlesEvolving[jscat].T ) / particlesEvolving[jscat].E;
        particlesEvolving[jscat].T = time;
        particlesEvolving[jscat].X = particlesEvolving[jscat].X + particlesEvolving[jscat].PX * c;
        particlesEvolving[jscat].Y = particlesEvolving[jscat].Y + particlesEvolving[jscat].PY * c;
        particlesEvolving[jscat].Z = particlesEvolving[jscat].Z + particlesEvolving[jscat].PZ * c;

        if ( time < particlesEvolving[iscat].T )
        {
          cout << "back23_i" << endl;
          int zz;
          cin >> zz;
        }
        if ( time < particlesEvolving[jscat].T )
        {
          cout << "back23_j" << endl;
          int zz;
          cin >> zz;
        }

        particlesEvolving[iscat].FLAVOR = static_cast<FLAVOR_TYPE>( F1 );
        particlesEvolving[iscat].PX = pxi;
        particlesEvolving[iscat].PY = pyi;
        particlesEvolving[iscat].PZ = pzi;
        particlesEvolving[iscat].E = sqrt( pxi * pxi + pyi * pyi + pzi * pzi );

        particlesEvolving[jscat].FLAVOR = static_cast<FLAVOR_TYPE>( F2 );
        particlesEvolving[jscat].PX = pxj;
        particlesEvolving[jscat].PY = pyj;
        particlesEvolving[jscat].PZ = pzj;
        particlesEvolving[jscat].E = sqrt( pxj * pxj + pyj * pyj + pzj * pzj );

        numberEvolvingParticles++;//production
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
        particlesEvolving[iscat].init = particlesEvolving[jscat].init = false;

        c = ( time - particlesEvolving[iscat].T ) / particlesEvolving[iscat].E;
        particlesEvolving[iscat].T = time;
        particlesEvolving[iscat].X = particlesEvolving[iscat].X + particlesEvolving[iscat].PX * c;
        particlesEvolving[iscat].Y = particlesEvolving[iscat].Y + particlesEvolving[iscat].PY * c;
        particlesEvolving[iscat].Z = particlesEvolving[iscat].Z + particlesEvolving[iscat].PZ * c;

        c = ( time - particlesEvolving[jscat].T ) / particlesEvolving[jscat].E;
        particlesEvolving[jscat].T = time;
        particlesEvolving[jscat].X = particlesEvolving[jscat].X + particlesEvolving[jscat].PX * c;
        particlesEvolving[jscat].Y = particlesEvolving[jscat].Y + particlesEvolving[jscat].PY * c;
        particlesEvolving[jscat].Z = particlesEvolving[jscat].Z + particlesEvolving[jscat].PZ * c;

        if ( time < particlesEvolving[iscat].T )
        {
          cout << "back32_i" << endl;
          int zz;
          cin >> zz;
        }
        if ( time < particlesEvolving[jscat].T )
        {
          cout << "back32_j" << endl;
          int zz;
          cin >> zz;
        }

        particlesEvolving[iscat].FLAVOR = static_cast<FLAVOR_TYPE>( F1 );
        particlesEvolving[iscat].PX = pxi;
        particlesEvolving[iscat].PY = pyi;
        particlesEvolving[iscat].PZ = pzi;
        particlesEvolving[iscat].E = sqrt( pxi * pxi + pyi * pyi + pzi * pzi );

        particlesEvolving[jscat].FLAVOR = static_cast<FLAVOR_TYPE>( F2 );
        particlesEvolving[jscat].PX = pxj;
        particlesEvolving[jscat].PY = pyj;
        particlesEvolving[jscat].PZ = pzj;
        particlesEvolving[jscat].E = sqrt( pxj * pxj + pyj * pyj + pzj * pzj );
        
        particlesEvolving[dead].dead = true;
        
        numberEvolvingParticles--;
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
        particlesEvolving[iscat].init = particlesEvolving[jscat].init = false;

        particlesEvolving[iscat].PXold = particlesEvolving[iscat].PX;
        particlesEvolving[iscat].PYold = particlesEvolving[iscat].PY;
        particlesEvolving[iscat].PZold = particlesEvolving[iscat].PZ;
        particlesEvolving[iscat].Eold = particlesEvolving[iscat].E;

        particlesEvolving[jscat].PXold = particlesEvolving[jscat].PX;
        particlesEvolving[jscat].PYold = particlesEvolving[jscat].PY;
        particlesEvolving[jscat].PZold = particlesEvolving[jscat].PZ;
        particlesEvolving[jscat].Eold = particlesEvolving[jscat].E;

        c = ( timei - particlesEvolving[iscat].T ) / particlesEvolving[iscat].E;
        particlesEvolving[iscat].T = timei;
        particlesEvolving[iscat].X = particlesEvolving[iscat].X + particlesEvolving[iscat].PX * c;
        particlesEvolving[iscat].Y = particlesEvolving[iscat].Y + particlesEvolving[iscat].PY * c;
        particlesEvolving[iscat].Z = particlesEvolving[iscat].Z + particlesEvolving[iscat].PZ * c;

        particlesEvolving[iscat].PX = pxi;
        particlesEvolving[iscat].PY = pyi;
        particlesEvolving[iscat].PZ = pzi;
        particlesEvolving[iscat].E = sqrt( pxi * pxi + pyi * pyi + pzi * pzi );

        c = ( timej - particlesEvolving[jscat].T ) / particlesEvolving[jscat].E;
        particlesEvolving[jscat].T = timej;
        particlesEvolving[jscat].X = particlesEvolving[jscat].X + particlesEvolving[jscat].PX * c;
        particlesEvolving[jscat].Y = particlesEvolving[jscat].Y + particlesEvolving[jscat].PY * c;
        particlesEvolving[jscat].Z = particlesEvolving[jscat].Z + particlesEvolving[jscat].PZ * c;

        particlesEvolving[jscat].PX = pxj;
        particlesEvolving[jscat].PY = pyj;
        particlesEvolving[jscat].PZ = pzj;
        particlesEvolving[jscat].E = sqrt( pxj * pxj + pyj * pyj + pzj * pzj );

        if ( timei < particlesEvolving[iscat].T )
        {
          cout << "back22e_i" << endl;
          int zz;
          cin >> zz;
        }
        if ( timej < particlesEvolving[jscat].T )
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

      particlesEvolving[iscat] = particlesEvolving[jscat];
    }
  }

  for ( int i = 0; i < numberEvolvingParticles; i++ )
  {
    if ( particlesEvolving[i].T < evolveToTime + 1.0e-8 )
    {
      particlesEvolving[i].init = false;
    }
  }

  // duplicate partcl from cascade to partclAtTimeNow
  particles_atTimeNow = particlesEvolving;
  
  // particles vector is larger and in size constant. But consider for partclAtTimenow only actual physical present particles (number is numberEvolvingParticles)
  particles_atTimeNow.resize( numberEvolvingParticles ); 

  // propagate particles_atTimeNow to current time
  for ( int i = 0; i < particles_atTimeNow.size(); i++ )
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

  for ( int i = 0; i < particles_atTimeNow.size(); i++ )
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

  for ( int i = 0; i < addedParticles.size(); i++ )
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
          addedParticles[i].temperature = -1.0;
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






void offlineHeavyIonCollision::scattering( const double nexttime, bool& again )
{
  double xt;
  double ee;
  
  int nGluons = 0;
//   int nAllLightQuarks = 0;
//   int nAllAntiLightQuarks = 0;
  int nGluonsAdded = 0;
//   int nAllLightQuarksAdded = 0;
//   int nAllAntiLightQuarksAdded = 0;
  
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

  gG = 2 * ( pow( static_cast<double>( Ncolor ) , 2 ) - 1 );
  
  if ( !formGeomCopy.empty() )
  {
    list<int>::iterator iIt;
    int id = -1;
    for ( iIt = formGeomCopy.begin(); iIt != formGeomCopy.end(); )
    {
      id = *iIt;
      cc = ( timenow - particles_atTimeNow[id].T ) / sqrt( pow( particles_atTimeNow[id].PXold, 2 ) + pow( particles_atTimeNow[id].PYold, 2 ) + pow( particles_atTimeNow[id].PZold, 2 ) + pow( particles_atTimeNow[id].m, 2 ) );
      zz = particles_atTimeNow[id].Z + particles_atTimeNow[id].PZold * cc;
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
    dz = simulationTime * ( tanh( etaBins[etaSliceIndex].right ) - tanh( etaBins[etaSliceIndex].left ) );
    
    rings.clear();
    rings.setLongitudinalGeometry( etaBins[etaSliceIndex].left, etaBins[etaSliceIndex].right, simulationTime );
    
    for ( int j = IXY * etaSliceIndex; j < IXY * ( etaSliceIndex + 1 ); j++ )
    {
      list<int>::const_iterator iIt;
      for ( iIt = cells[j].particleList.begin(); iIt != cells[j].particleList.end(); iIt++ )
      {
        xt = sqrt( pow( particles_atTimeNow [( *iIt )].X, 2 )  + pow( particles_atTimeNow[( *iIt )].Y, 2 ) );

        // calculation of energy density & co. is only written for massless particles, only add those -> for the medium properties heavy quarks do not contribute anyhow
        if ( !particles_atTimeNow[( *iIt )].dead && particles_atTimeNow[( *iIt )].FLAVOR <= 2 * Particle::max_N_light_flavor )
        {
          rings.addParticle( particles_atTimeNow[( *iIt )] );
//           rings.addRates( particles_atTimeNow[(*iIt)] );
        }
      }
    }
    
    list<int>::iterator iIt;
    for ( iIt = formGeomCopy.begin(); iIt != formGeomCopy.end(); iIt++ )
    {
      int id = *iIt;
      
      if ( particles_atTimeNow[id].dead )
      {
        continue; // jump to next particle in the list if this one is dead (should not happen)
      }
      
      ee = sqrt( pow( particles_atTimeNow[id].PXold, 2 ) + pow( particles_atTimeNow[id].PYold, 2 ) + pow( particles_atTimeNow[id].PZold, 2 ) + pow( particles_atTimeNow[id].m, 2 ) );
      cc = ( simulationTime - particles_atTimeNow[id].T ) / ee;
      zz = particles_atTimeNow[id].Z + particles_atTimeNow[id].PZold * cc;
      eta = 0.5 * log(( simulationTime + zz ) / ( simulationTime - zz ) );
      
      if (( eta >= etaBins[etaSliceIndex].left ) && ( eta <= etaBins[etaSliceIndex].right ) )
      {
        // calculation of energy density & co. is only written for massless particles, only add those -> for the medium properties heavy quarks do not contribute anyhow
        if( particles_atTimeNow[( *iIt )].FLAVOR <= 2 * Particle::max_N_light_flavor ) 
        {
          rings.addParticleInFormGeom( particles_atTimeNow[id], simulationTime );
        }
        iIt = formGeomCopy.erase( iIt );    // erase this particle such that it needs not be looped over for the next rapidity slab
      }
    }
    
    rings.prepareAverages( dz, testpartcl );
    if ( etaSliceIndex == centralEtaIndex )
    {
      theAnalysis->centralRingsCopyFromCascade = rings;
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
            vector<int> nLightQuarks( Particle::N_light_flavor , 0 );
            vector<int> nAntiLightQuarks( Particle::N_light_flavor, 0 );
            int nCharmQuarks = 0;
            int nAntiCharmQuarks = 0;
            int nBottomQuarks = 0;
            int nAntiBottomQuarks = 0;
            vector<int> nLightQuarksAdded( Particle::N_light_flavor , 0 );
            vector<int> nAntiLightQuarksAdded( Particle::N_light_flavor, 0 );
            int nCharmQuarksAdded = 0;
            int nAntiCharmQuarksAdded = 0;
            int nBottomQuarksAdded = 0;
            int nAntiBottomQuarksAdded = 0;
            
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
                    ++nLightQuarks[0];
                    break;
                  case down:
                    ++nLightQuarks[1];
                    break;
                  case strange:
                    ++nLightQuarks[2];
                    break;
                  case charm:
                    ++nCharmQuarks;
                    break;
                  case bottom:
                    ++nBottomQuarks;
                    break;
                  case anti_up:
                    ++nAntiLightQuarks[0];
                    break;
                  case anti_down:
                    ++nAntiLightQuarks[1];
                    break;
                  case anti_strange:
                    ++nAntiLightQuarks[2];
                    break;
                  case anti_charm:
                    ++nAntiCharmQuarks;
                    break;
                  case anti_bottom:
                    ++nAntiBottomQuarks;
                    break;
                  default:
                    break;
                }
              }
              
              xt = sqrt( pow( particles_atTimeNow[id].X, 2 ) + pow( particles_atTimeNow[id].Y, 2 ) );
              nc = rings.getIndex( xt );
              
              particles_atTimeNow[id].md2g = rings[nc].getAveraged_md2g();
              particles_atTimeNow[id].md2q = rings[nc].getAveraged_md2q();
              particles_atTimeNow[id].temperature = rings[nc].getEffectiveTemperature();
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
                    ++nLightQuarks[0];
                    break;
                  case down:
                    ++nLightQuarks[1];
                    break;
                  case strange:
                    ++nLightQuarks[2];
                    break;
                  case charm:
                    ++nCharmQuarks;
                    break;
                  case bottom:
                    ++nBottomQuarks;
                    break;
                  case anti_up:
                    ++nAntiLightQuarks[0];
                    break;
                  case anti_down:
                    ++nAntiLightQuarks[1];
                    break;
                  case anti_strange:
                    ++nAntiLightQuarks[2];
                    break;
                  case anti_charm:
                    ++nAntiCharmQuarks;
                    break;
                  case anti_bottom:
                    ++nAntiBottomQuarks;
                    break;
                  default:
                    break;
                }
              }
              
              xt = sqrt( pow( addedParticles[id].X, 2 ) + pow( addedParticles[id].Y, 2 ) );
              nc = rings.getIndex( xt );
              
              addedParticles[id].md2g = rings[nc].getAveraged_md2g();
              addedParticles[id].md2q = rings[nc].getAveraged_md2q();
              addedParticles[id].temperature = rings[nc].getEffectiveTemperature();
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
            
            if( theConfig->isScatt_offlineWithAddedParticles() && theConfig->doScattering_32() )
            {
              scatt32_offlineWithAddedParticles( cells[j], allParticlesList, gluonList, cellsAdded[j], allParticlesListAdded, gluonListAdded, n32, again, nexttime );
              if ( again )
              {
                return;
              }
            }
            
            
            double scaleFactor = 1;
            if ( nGluons > n32 )
            {
              scaleFactor = static_cast<double>( nGluons ) / static_cast<double>( nGluons - n32 );
            }
            
            if( theConfig->isScatt_amongAddedParticles() && theConfig->doScattering_22() )
            {
              scatt22_amongAddedParticles( cellsAdded[j], allParticlesListAdded, scaleFactor, again, nexttime );
              
              if ( again )
              {
                return;
              }
            }
            
            if( theConfig->isScatt_offlineWithAddedParticles() )
            {
              analysisRingStructure tempRing( theAnalysis->rings.size(), theAnalysis->rings.getCentralRadius(), theAnalysis->rings.getDeltaR() );
              scatt2223_offlineWithAddedParticles( cells[j], allParticlesList, gluonList, cellsAdded[j], allParticlesListAdded, gluonListAdded, scaleFactor, again, nexttime, tempRing );
              if ( etaSliceIndex == centralEtaIndex )
              {
                theAnalysis->rings += tempRing;
              }
              
              if ( again )
              {
                return;
              }
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
                  particles_atTimeNow[id].temperature = rings[nc].getEffectiveTemperature();
                }
                
                if ( particles_atTimeNow[id].FLAVOR == gluon )
                {
                  particles_atTimeNow[id].ratev = rateGluons[etaSliceIndex][nc];
                }
                else if ( particles_atTimeNow[id].FLAVOR == up || particles_atTimeNow[id].FLAVOR == down || particles_atTimeNow[id].FLAVOR == strange || particles_atTimeNow[id].FLAVOR == light_quark )
                {
                  particles_atTimeNow[id].ratev = rateQuarks[etaSliceIndex][nc];
                }
                else if ( particles_atTimeNow[id].FLAVOR == anti_up || particles_atTimeNow[id].FLAVOR == anti_down || particles_atTimeNow[id].FLAVOR == anti_strange || particles_atTimeNow[id].FLAVOR == anti_light_quark )
                {
                  particles_atTimeNow[id].ratev = rateAntiQuarks[etaSliceIndex][nc];
                }
                else
                {
                  particles_atTimeNow[id].ratev = 0.0;
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
                  addedParticles[id].temperature = rings[nc].getEffectiveTemperature();
                }
                
                if ( addedParticles[id].FLAVOR == gluon )
                {
                  addedParticles[id].ratev = rateGluons[etaSliceIndex][nc];
                }
                else if ( addedParticles[id].FLAVOR == up || addedParticles[id].FLAVOR == down || addedParticles[id].FLAVOR == strange || addedParticles[id].FLAVOR == light_quark )
                {
                  addedParticles[id].ratev = rateQuarks[etaSliceIndex][nc];
                }
                else  if ( addedParticles[id].FLAVOR == anti_up || addedParticles[id].FLAVOR == anti_down || addedParticles[id].FLAVOR == anti_strange || addedParticles[id].FLAVOR == anti_light_quark )
                {
                  addedParticles[id].ratev = rateAntiQuarks[etaSliceIndex][nc];
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
  
  // J/psi dissociation: if temperature in cell is higher than Td = 2 Tc, decay J/psi to two charm quarks
  if( Particle::N_psi_states > 0 )
    jpsi_dissociation_td( nexttime );
}





void offlineHeavyIonCollision::scatt2223_offlineWithAddedParticles( cellContainer& _cell, std::vector< int >& _allParticlesList, std::vector< int >& _gluonList,
    cellContainer& _cellAdded, std::vector< int >& _allParticlesListAdded, std::vector< int >& _gluonListAdded,
    const double scaleFactor, bool& again, const double nexttime, analysisRingStructure& _analysisRings )
{
  const double epsilon = 1.0e-4;
  int iscat, jscat, typ;
  double P1[4], P2[4];
  double s, as, csgg, cs23, cs22, Vrel, lambda_scaled;
  double M1, M2, P1P2;
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
//   vector<int> nQuarks( Particle::N_light_flavor, 0 );
//   vector<int> nAntiQuarks( Particle::N_light_flavor, 0 );

  scattering23 scatt23_object( &theI23 );
  scattering22 scatt22_object( &theI22 );
  
  double lambdaJet = 0;
  double pt_addedParticle = 0;

  //   const int nAllQuarks = std::accumulate( nQuarks.begin(), nQuarks.end(), 0 ) + std::accumulate( nAntiQuarks.begin(), nAntiQuarks.end(), 0 );
//   const int allPairs = binomial( nTotal, 2 );
//   const int consideredPairs = 25;

  for ( int j = 0; j < static_cast<int>( _allParticlesListAdded.size() ); j++ )
  {
    jscat = _allParticlesListAdded[j];

    pt_addedParticle = sqrt( pow( addedParticles[jscat].PX, 2.0 ) + pow( addedParticles[jscat].PY, 2.0 ) );
    
    if ( pt_addedParticle < theConfig->getMinimumPT() || addedParticles[jscat].dead )
    {
      continue; // jump to next particle in the list
    }
    
    if ( pt_addedParticle > 8.0 && addedParticles[jscat].FLAVOR <= 2 * Particle::max_N_light_flavor )
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

    for ( int i = 0; i < static_cast<int>( _allParticlesList.size() ); i++ )
    {
      iscat = _allParticlesList[i];
      
      if ( particles_atTimeNow[iscat].dead || addedParticles[jscat].dead ) // The last statement to check again whether jscat is dead is necessary since jscat could be deleted in the scattering (Jpsi dissociation which is not "real" since we use more testparticles for Jpsi in which case the Jpsi is deleted) before but the scatterings with all remaining jscat would be carried out.
      {
        continue; // jump to next particle from _cell.particleList
      }
      
      F1 = particles_atTimeNow[iscat].FLAVOR;
      particles_atTimeNow[iscat].getMomentumArray( P1 );
      M1 = particles_atTimeNow[iscat].m;
      F2 = addedParticles[jscat].FLAVOR;
      addedParticles[jscat].getMomentumArray( P2 );
      M2 = addedParticles[jscat].m;

      s = pow(( P1[0] + P2[0] ), 2 ) - pow(( P1[1] + P2[1] ), 2 ) - pow(( P1[2] + P2[2] ), 2 ) - pow(( P1[3] + P2[3] ), 2 );

      xt = ( sqrt( pow( particles_atTimeNow[iscat].X, 2 )  + pow( particles_atTimeNow[iscat].Y, 2 ) ) +
             sqrt( pow( addedParticles[jscat].X, 2 )  + pow( addedParticles[jscat].Y, 2 ) ) ) / 2.0;

      ringIndex = rings.getIndex( xt );

      averagedRate = ( particles_atTimeNow[iscat].rate + addedParticles[jscat].rate +
                       particles_atTimeNow[iscat].ratev + addedParticles[jscat].ratev ) / ( 2.0 * 2.0 );


      if ( s < 1.1*lambda2 )
      {
        probab2322 = -1.0;
        cs22 = cs23 = 0.0; //1/GeV^2
      }
      else
      {
        P1P2 = P1[0]*P2[0] - P1[1]*P2[1] - P1[2]*P2[2] - P1[3]*P2[3];
        Vrel = pow((pow(P1P2,2.0)-pow(M1,2.0)*pow(M2,2.0)),0.5)/(P1[0]*P2[0]); // general relative velocity
        as = coupling::get_constant_coupling();

        //factor as (alpha_s) is not included in definitions of partcl[jscat].md2g, therefore multiplied here
        md2g = as * ( particles_atTimeNow[iscat].md2g + addedParticles[jscat].md2g ) / 2.0;
        md2q = as * ( particles_atTimeNow[iscat].md2q + addedParticles[jscat].md2q ) / 2.0;
        
        // HACK for N_f = 0 background!!!!
//         md2g = as * ( particles_atTimeNow[iscat].md2g + addedParticles[jscat].md2g ) / 2.0 * 2;
//         md2q = md2g * 2.0 / 9.0;
        
        if( theConfig->doScattering_22() )
        {
          scatt22_object.setParameter( P1, P2, F1, F2, M1, M2, s, md2g / as , md2q / as,
                                      theConfig->getKggQQb(), theConfig->getKgQgQ(), theConfig->getKappa_gQgQ(), 
                                      theConfig->isConstantCrossSecGQ(),
                                      theConfig->getConstantCrossSecValueGQ(), theConfig->isIsotropicCrossSecGQ() ); // md2g, md2q are debye masses without the factor alpha_s which is multiplied in scattering22.cpp
          cs22 = scatt22_object.getXSection22( initialStateIndex );

          probab22 = pow( 0.197, 2.0 ) * cs22 * Vrel * dt / ( dv * testpartcl );
        }
        else
        {
          cs22 = probab22 = 0.0;
        }
        

        if( theConfig->doScattering_23() )
        {
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

          betaDistEntry = scatt23_object.setParameter( vx, vy, vz, P1, P2, F1, F2, sqrt( s ), md2g / s, lambda_scaled, as, _gluonList.size() );

          cs23 = scatt23_object.getXSection23();

          probab23 = pow( 0.197, 2.0 ) * cs23 * Vrel * dt / ( dv * testpartcl );
        }
        else
        {  
          probab23 = 0;
          cs23 = 0;
          lambda_scaled = 0;
        }
        
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
          if( pt_jscat > theAnalysis->getJetTracking_PT() )
          {
            jetEventIndex = theAnalysis->addJetEvent_in( iscat, -1, jscat, c2to3, cs23, _cell.index, lambda_scaled / sqrt( s ) );
          }

          int newIndex = scatt23_offlineWithAddedParticles_utility( scatt23_object, _cell, iscat, jscat, again, nexttime );

          pt_jscat = sqrt( pow( addedParticles[jscat].PX, 2.0 ) + pow( addedParticles[jscat].PY, 2.0 ) );
          pt_nmb = sqrt( pow( addedParticles[newIndex].PX, 2.0 ) + pow( addedParticles[newIndex].PY, 2.0 ) );
          if( jetEventIndex != -1 || pt_jscat > theAnalysis->getJetTracking_PT() || pt_nmb > theAnalysis->getJetTracking_PT() )
          {
            theAnalysis->addJetEvent_out( jetEventIndex, jscat, iscat, newIndex, c2to3 );
          }
        }
        else
        {
          int jetEventIndex = -1;
          if( pt_jscat > theAnalysis->getJetTracking_PT() )
          {
            jetEventIndex = theAnalysis->addJetEvent_in( iscat, -1, jscat, c2to2, cs22, _cell.index, lambda_scaled / sqrt( s ) );
          }

          scatt22_offlineWithAddedParticles_utility( scatt22_object, _cellAdded.particleList, _allParticlesListAdded, iscat, jscat, typ, nexttime );

          pt_jscat = sqrt( pow( addedParticles[jscat].PX, 2.0 ) + pow( addedParticles[jscat].PY, 2.0 ) );
          if( jetEventIndex != -1 || pt_jscat > theAnalysis->getJetTracking_PT() )
          {
            theAnalysis->addJetEvent_out( jetEventIndex, jscat, iscat, -1, c2to2 );
          }
        }
      }
    }
  }

  double cc;
  list<int>::const_iterator iIt;
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




void offlineHeavyIonCollision::scatt22_amongAddedParticles( cellContainer& _cellAdded, std::vector< int >& _allParticlesListAdded, const double scaleFactor, bool& again, const double nexttime )
{
  int iscat, jscat, typ;
  double P1[4], P2[4];
  double s, as, cs22, Vrel;
  double M1, M2, P1P2;
  double probab22;
  double md2g, md2q;
  int initialStateIndex = -1;
  FLAVOR_TYPE F1, F2;
  double temperature;


  list<int>::const_iterator iIt;
  list<int>::const_iterator jIt;

  scattering22 scatt22_object( &theI22 );
  
  for ( int i = 0; i < static_cast<int>( _allParticlesListAdded.size() ) - 1; i++ )
  {
    iscat = _allParticlesListAdded[i];

    if ( !( addedParticles[iscat].FLAVOR == charm || addedParticles[iscat].FLAVOR == anti_charm) || addedParticles[iscat].dead ) // currently scattering among added particles is only used for c+cbar -> Jpsi + g
    {
      continue; // jump to next particle in the list
    }

    for ( int j = i + 1; j < _allParticlesListAdded.size(); j++ )
    {
      jscat = _allParticlesListAdded[j];

      if ( !( addedParticles[jscat].FLAVOR == charm || addedParticles[jscat].FLAVOR == anti_charm) || addedParticles[jscat].dead || addedParticles[iscat].dead  ) // currently scattering among added particles is only used for c+cbar -> Jpsi + g. The last statement to check again whether iscat is dead is necessary since iscat could be deleted in the scattering before but the scatterings with all remaining jscat would be carried out.
      {
        continue; // jump to next particle from _cell.particleList
      }
      
      F1 = addedParticles[iscat].FLAVOR;
      addedParticles[iscat].getMomentumArray( P1 );
      M1 = addedParticles[iscat].m;
      F2 = addedParticles[jscat].FLAVOR;
      addedParticles[jscat].getMomentumArray( P2 );
      M2 = addedParticles[jscat].m;

      s = pow(( P1[0] + P2[0] ), 2 ) - pow(( P1[1] + P2[1] ), 2 ) - pow(( P1[2] + P2[2] ), 2 ) - pow(( P1[3] + P2[3] ), 2 );

      if ( s < 1.1*lambda2 )
      {
        probab22 = -1.0;
        cs22 = 0.0; //1/GeV^2
      }
      else
      {
        P1P2 = P1[0]*P2[0] - P1[1]*P2[1] - P1[2]*P2[2] - P1[3]*P2[3];
        Vrel = pow((pow(P1P2,2.0)-pow(M1,2.0)*pow(M2,2.0)),0.5)/(P1[0]*P2[0]); // general relative velocity
        as = coupling::get_constant_coupling();

        //factor as (alpha_s) is not included in definitions of partcl[jscat].md2g, therefore multiplied here
        md2g = as * ( addedParticles[iscat].md2g + addedParticles[jscat].md2g ) / 2.0;
        md2q = as * ( addedParticles[iscat].md2q + addedParticles[jscat].md2q ) / 2.0;
        
        // determine largest temp of both charm quarks
        if( addedParticles[iscat].temperature > addedParticles[jscat].temperature )
          temperature = addedParticles[iscat].temperature;
        else
          temperature = addedParticles[jscat].temperature;

        scatt22_object.setParameter( P1, P2, F1, F2, M1, M2, s, md2g / as , md2q / as,
                                     theConfig->getKggQQb(), theConfig->getKgQgQ(), theConfig->getKappa_gQgQ(), theConfig->isConstantCrossSecGQ(),
                                     theConfig->getConstantCrossSecValueGQ(), theConfig->isIsotropicCrossSecGQ(),
                                     temperature, theConfig->getTdJpsi(), theConfig->isConstantCrossSecJpsi(), theConfig->getConstantCrossSecValueJpsi()
                                     ); // md2g, md2q are debye masses without the factor alpha_s which is multiplied in scattering22.cpp
      
        cs22 = scatt22_object.getXSection22( initialStateIndex );

        // multipy with Jpsi testparticle number, produce Jpsi in every choosen collision, but annihilate ccbar only in every 1/N_test_Jpsi collision to get the correct rate
        probab22 = pow( 0.197, 2.0 ) * cs22 * Vrel * dt / ( dv * testpartcl * theConfig->getNaddedEvents() ) * theConfig->getJpsiTestparticles();

        if ( F1 == gluon )
        {
          probab22 *= scaleFactor;
        }
        if ( F2 == gluon )
        {
          probab22 *= scaleFactor;
        }
      }

      if ( probab22 > 1.0 )
      {
//         cout << "P2322=" << probab2322 << ">1" << endl;
        again = true;
//         cout << "dt (old) = " << dt << endl;
        dt = 0.5 / ( probab22 / dt );
//         cout << "dt (new) = " << dt << endl;
        return;
      }

      if ( ran2() < probab22 )
      {
        scatt22_amongAddedParticles_utility( scatt22_object, _cellAdded.particleList, _allParticlesListAdded, iscat, jscat, typ, nexttime );
      }
    }
  }
}



void offlineHeavyIonCollision::scatt32_offlineWithAddedParticles( cellContainer& _cell, std::vector< int >& _allParticlesList, std::vector< int >& _gluonList, cellContainer& _cellAdded, std::vector< int >& _allParticlesListAdded, std::vector< int >& _gluonListAdded, int& n32, bool& again, const double nexttime )
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
    if ( pt_addedParticle < theConfig->getMinimumPT() )
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
      if ( pt_addedParticle < theConfig->getMinimumPT() )
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
        as = coupling::get_constant_coupling();

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

        betaDistEntry = scatt32_object.setParameter( vx, vy, vz, P1, P2, P3, F1, F2, F3, sqrt( s ), md2g / s, lambda_scaled, as, _gluonList.size() );
        I32 = scatt32_object.getIntegral32_withPrefactors();                        // get the integral I32 for the given 3 particles

        probab32 = scaleForSelectedTriplets * I32 * dt / ( pow( dv, 2 ) * pow( static_cast<double>( testpartcl ), 2 ) );
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
        if( pt_kscat > theAnalysis->getJetTracking_PT() )
        {
          jetEventIndex = theAnalysis->addJetEvent_in( iscat, jscat, kscat, c3to2, I32, _cell.index, lambda_scaled / sqrt( s ) );
        }

        order = scatt32_offlineWithAddedParticles_utility( scatt32_object, _cellAdded.particleList, _allParticlesListAdded, _gluonListAdded, iscat, jscat, kscat, n32, nexttime );

        if ( order == 4 || order == 6 )
        {
          pt_jscat = sqrt( pow( particles_atTimeNow[jscat].PX, 2.0 ) + pow( particles_atTimeNow[jscat].PY, 2.0 ) );
          pt_kscat = sqrt( pow( addedParticles[kscat].PX, 2.0 ) + pow( addedParticles[kscat].PY, 2.0 ) );
          if( jetEventIndex != -1 || pt_jscat > theAnalysis->getJetTracking_PT() || pt_kscat > theAnalysis->getJetTracking_PT() )
          {
            theAnalysis->addJetEvent_out( jetEventIndex, kscat, jscat, -1, c3to2 );
          }
        }
        else if ( order == 2 || order == 5 )
        {
          pt_iscat = sqrt( pow( particles_atTimeNow[iscat].PX, 2.0 ) + pow( particles_atTimeNow[iscat].PY, 2.0 ) );
          pt_kscat = sqrt( pow( addedParticles[kscat].PX, 2.0 ) + pow( addedParticles[kscat].PY, 2.0 ) );
          if( jetEventIndex != -1 || pt_iscat > theAnalysis->getJetTracking_PT() || pt_kscat > theAnalysis->getJetTracking_PT() )
          {
            theAnalysis->addJetEvent_out( jetEventIndex, kscat, iscat, -1, c3to2 );
          }
        }
        else if ( order == 1 || order == 3 )
        {
          pt_iscat = sqrt( pow( particles_atTimeNow[iscat].PX, 2.0 ) + pow( particles_atTimeNow[iscat].PY, 2.0 ) );
          pt_jscat = sqrt( pow( particles_atTimeNow[jscat].PX, 2.0 ) + pow( particles_atTimeNow[jscat].PY, 2.0 ) );
          if( jetEventIndex != -1 || pt_iscat > theAnalysis->getJetTracking_PT() || pt_jscat > theAnalysis->getJetTracking_PT() )
          {
            theAnalysis->addJetEvent_out( jetEventIndex, kscat, jscat, -1, c3to2 );
          }
        }
        else
        {
          if ( jetEventIndex != -1 )
          {
            theAnalysis->removeJetEvent_in( jetEventIndex );
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
      if ( pt_addedParticle < theConfig->getMinimumPT() )
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
      
      for ( int m1 = 0; m1 < static_cast<int>( _allParticlesList.size() ) - 1; m1++ )
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
            as = coupling::get_constant_coupling();

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

            betaDistEntry = scatt32_object.setParameter( vx, vy, vz, P1, P2, P3, F1, F2, F3, sqrt( s ), md2g / s, lambda_scaled, as, _gluonList.size() );
            I32 = scatt32_object.getIntegral32_withPrefactors();                        // get the integral I32 for the given 3 particles

            probab32 = I32 * dt / ( pow( dv, 2 ) * pow( static_cast<double>( testpartcl ), 2 ) ); 
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
            if( pt_kscat > theAnalysis->getJetTracking_PT() )
            {
              jetEventIndex = theAnalysis->addJetEvent_in( iscat, jscat, kscat, c3to2, I32, _cell.index, lambda_scaled / sqrt( s ) );
            }

            order = scatt32_offlineWithAddedParticles_utility( scatt32_object, _cellAdded.particleList, _allParticlesListAdded, _gluonListAdded, iscat, jscat, kscat, n32, nexttime );

            if ( order == 4 || order == 6 )
            {
              pt_jscat = sqrt( pow( particles_atTimeNow[jscat].PX, 2.0 ) + pow( particles_atTimeNow[jscat].PY, 2.0 ) );
              pt_kscat = sqrt( pow( addedParticles[kscat].PX, 2.0 ) + pow( addedParticles[kscat].PY, 2.0 ) );
              if( jetEventIndex != -1 || pt_jscat > theAnalysis->getJetTracking_PT() || pt_kscat > theAnalysis->getJetTracking_PT() )
              {
                theAnalysis->addJetEvent_out( jetEventIndex, kscat, jscat, -1, c3to2 );
              }
            }
            else if ( order == 2 || order == 5 )
            {
              pt_iscat = sqrt( pow( particles_atTimeNow[iscat].PX, 2.0 ) + pow( particles_atTimeNow[iscat].PY, 2.0 ) );
              pt_kscat = sqrt( pow( addedParticles[kscat].PX, 2.0 ) + pow( addedParticles[kscat].PY, 2.0 ) );
              if( jetEventIndex != -1 || pt_iscat > theAnalysis->getJetTracking_PT() || pt_kscat > theAnalysis->getJetTracking_PT() )
              {
                theAnalysis->addJetEvent_out( jetEventIndex, kscat, iscat, -1, c3to2 );
              }
            }
            else if ( order == 1 || order == 3 )
            {
              pt_iscat = sqrt( pow( particles_atTimeNow[iscat].PX, 2.0 ) + pow( particles_atTimeNow[iscat].PY, 2.0 ) );
              pt_jscat = sqrt( pow( particles_atTimeNow[jscat].PX, 2.0 ) + pow( particles_atTimeNow[jscat].PY, 2.0 ) );
              if( jetEventIndex != -1 || pt_iscat > theAnalysis->getJetTracking_PT() || pt_jscat > theAnalysis->getJetTracking_PT() )
              {
                theAnalysis->addJetEvent_out( jetEventIndex, kscat, jscat, -1, c3to2 );
              }
              
              //                 hack for exiting both inner loops and continuing with next addedParticle 
              m1 = static_cast<int>( _allParticlesList.size() ) - 1; 
              m2 = _allParticlesList.size();        
              m3--; // m3 is decreased because in the next loop step it is again increased by 1 
            }
            else
            {
              if ( jetEventIndex != -1 )
              {
                theAnalysis->removeJetEvent_in( jetEventIndex );
              }
            }
          }
        }
      }
    }
  }
}





int offlineHeavyIonCollision::scatt23_offlineWithAddedParticles_utility( scattering23& scatt23_obj, cellContainer& _cell, int iscat, const int jscat, bool& again, const double nexttime )
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

  if( theConfig->doOutput_scatteredMediumParticles() )
  {
    scatteredMediumParticles.push_back( particles_atTimeNow[iscat] );
    scatteredMediumParticles.back().N_EVENT_pp = addedParticles[jscat].N_EVENT_pp;
  }

  particles_atTimeNow[iscat].getMomentumArray( P1 );
  particles_atTimeNow[iscat].getCoordinateArray( R1 );
  F1 = particles_atTimeNow[iscat].FLAVOR;

  addedParticles[jscat].getMomentumArray( P2 );
  addedParticles[jscat].getCoordinateArray( R2 );
  F2 = addedParticles[jscat].FLAVOR;
  
  // these routines are not written for heavy flavors
  if( F1 > 2*Particle::max_N_light_flavor || F2 > 2*Particle::max_N_light_flavor )
  {
    std::string errMsg = "Heavy flavor in scatt23_offlineWithAddedParticles_utility.";
    throw eHIC_error( errMsg );
  }

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
  // consider outgoing particle with highest pt if it not a tagged jet (charm, bottom, jpsi, etc)
  if ( pt_out1 > pt_out2 && addedParticles[jscat].FLAVOR <= 2 * Particle::max_N_light_flavor )
  {
    addedParticles[jscat].FLAVOR = F1;
    addedParticles[jscat].PX = P1[1];
    addedParticles[jscat].PY = P1[2];
    addedParticles[jscat].PZ = P1[3];
    addedParticles[jscat].E = sqrt( pow( P1[1], 2 ) + pow( P1[2], 2 ) + pow( P1[3], 2 ) );
    ParticleOffline tempParticle = particles_atTimeNow[iscat];
    tempParticle.FLAVOR = F2;
    tempParticle.PX = P2[1];
    tempParticle.PY = P2[2];
    tempParticle.PZ = P2[3];
    tempParticle.E = sqrt( P2[1] * P2[1] + P2[2] * P2[2] + P2[3] * P2[3] );
    tempParticle.unique_id = particles_atTimeNow[iscat].unique_id;
    tempParticle.N_EVENT_pp = addedParticles[jscat].N_EVENT_pp;
    tempParticle.N_EVENT_AA = addedParticles[jscat].N_EVENT_AA;
    if( theConfig->isScatt_furtherOfflineParticles() )
      addedParticles.push_back( tempParticle );
  }
  else
  {
    addedParticles[jscat].FLAVOR = F2;
    addedParticles[jscat].PX = P2[1];
    addedParticles[jscat].PY = P2[2];
    addedParticles[jscat].PZ = P2[3];
    addedParticles[jscat].E = sqrt( pow( P2[1], 2 ) + pow( P2[2], 2 ) + pow( P2[3], 2 ) );
    ParticleOffline tempParticle = particles_atTimeNow[iscat];
    tempParticle.FLAVOR = F1;
    tempParticle.PX = P1[1];
    tempParticle.PY = P1[2];
    tempParticle.PZ = P1[3];
    tempParticle.E = sqrt( P1[1] * P1[1] + P1[2] * P1[2] + P1[3] * P1[3] );
    tempParticle.unique_id = particles_atTimeNow[iscat].unique_id;
    tempParticle.N_EVENT_pp = addedParticles[jscat].N_EVENT_pp;
    tempParticle.N_EVENT_AA = addedParticles[jscat].N_EVENT_AA;
    if( theConfig->isScatt_furtherOfflineParticles() )
      addedParticles.push_back( tempParticle );
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
  tempParticle.Y_init = tempParticle.Y;
  tempParticle.Z_init = tempParticle.Z;
  tempParticle.PX_init = tempParticle.PX;
  tempParticle.PY_init = tempParticle.PY;
  tempParticle.PZ_init = tempParticle.PZ;
  tempParticle.E_init = -tempParticle.E;   // negative to indicate creation via 2->3 process
  tempParticle.m = 0.0;
  tempParticle.md2g = ( particles_atTimeNow[iscat].md2g + addedParticles[jscat].md2g ) / 2.0;
  tempParticle.md2q = ( particles_atTimeNow[iscat].md2q + addedParticles[jscat].md2q ) / 2.0;
  tempParticle.temperature = ( particles_atTimeNow[iscat].temperature + addedParticles[jscat].temperature ) / 2.0;
  tempParticle.edge = false;
  tempParticle.dead = false;
  tempParticle.free = false;
  tempParticle.cell_id = -1;
  tempParticle.coll_id = ncoll;
  tempParticle.collisionTime = infinity;
  tempParticle.collisionPartner = -1;
  tempParticle.rate = particles_atTimeNow[iscat].rate;    //GeV
  tempParticle.ratev = particles_atTimeNow[iscat].ratev;    //GeV
  tempParticle.initially_produced = false;
  tempParticle.init = false;
  tempParticle.X_traveled = 0.0;
  tempParticle.T_creation = tempParticle.T;
  tempParticle.unique_id = ParticleOffline::unique_id_counter_added;
  --ParticleOffline::unique_id_counter_added;
  tempParticle.N_EVENT_pp = addedParticles[jscat].N_EVENT_pp;
  tempParticle.N_EVENT_AA = addedParticles[jscat].N_EVENT_AA;

//   if ( sqrt( pow( tempParticle.PX, 2) + pow( tempParticle.PY, 2) ) > 3.0 )
//   {
    addedParticles.push_back( tempParticle );
    int newIndex = addedParticles.size() - 1;

    cc = ( nexttime - addedParticles[newIndex].T ) / addedParticles[newIndex].E;
    addedParticles[newIndex].T = nexttime;
    addedParticles[newIndex].X = addedParticles[newIndex].X + addedParticles[newIndex].PX * cc;
    addedParticles[newIndex].Y = addedParticles[newIndex].Y + addedParticles[newIndex].PY * cc;
    addedParticles[newIndex].Z = addedParticles[newIndex].Z + addedParticles[newIndex].PZ * cc;
//   }
  
  return newIndex;
}



void offlineHeavyIonCollision::scatt22_offlineWithAddedParticles_utility( scattering22& scatt22_obj, std::list< int >& _cellMembersAdded, std::vector< int >& _allParticlesListAdded, const int iscat, const int jscat, int& typ, const double nexttime )
{
  FLAVOR_TYPE F1, F2;
  double P1[4], P2[4], R1[4], R2[4];
  double as, Tmax, TT, cc, t;
  double M1, M2;
  double t_hat;

  if( theConfig->doOutput_scatteredMediumParticles() )
  {
    scatteredMediumParticles.push_back( particles_atTimeNow[iscat] );
    scatteredMediumParticles.back().N_EVENT_pp = addedParticles[jscat].N_EVENT_pp;
  }

  particles_atTimeNow[iscat].getMomentumArray( P1 );
  particles_atTimeNow[iscat].getCoordinateArray( R1 );
  F1 = particles_atTimeNow[iscat].FLAVOR;
  M1 = particles_atTimeNow[iscat].m;

  addedParticles[jscat].getMomentumArray( P2 );
  addedParticles[jscat].getCoordinateArray( R2 );
  F2 = addedParticles[jscat].FLAVOR;
  M2 = addedParticles[jscat].m;


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

  // determine type of scattering, momentum transfer t_hat, new flavor and new masses
  scatt22_obj.getMomentaAndMasses22( F1, F2, M1, M2, t_hat, typ );
  // translate momemtum transfer t_hat into actual momenta of outgoing particles
  scatt22_obj.setNewMomenta22( P1, P2, R1, R2, t_hat );

  double pt_out1 = sqrt( pow( P1[1], 2 ) + pow( P1[2], 2 ) );
  double pt_out2 = sqrt( pow( P2[1], 2 ) + pow( P2[2], 2 ) );

  //<<---------------------------------------------
  // set new properties for added particle
  // consider outgoing particle with highest pt if it not a tagged jet (charm, bottom, jpsi, etc)
  if ( pt_out1 > pt_out2 && addedParticles[jscat].FLAVOR <= 2 * Particle::max_N_light_flavor )
  {
    addedParticles[jscat].FLAVOR = F1;
    addedParticles[jscat].m = M1;
    addedParticles[jscat].PX = P1[1];
    addedParticles[jscat].PY = P1[2];
    addedParticles[jscat].PZ = P1[3];
    addedParticles[jscat].E = P1[0];
    ParticleOffline tempParticle = particles_atTimeNow[iscat];
    tempParticle.FLAVOR = F2;
    tempParticle.PX = P2[1];
    tempParticle.PY = P2[2];
    tempParticle.PZ = P2[3];
    tempParticle.E = sqrt( P2[1] * P2[1] + P2[2] * P2[2] + P2[3] * P2[3] );
    tempParticle.unique_id = particles_atTimeNow[iscat].unique_id;
    tempParticle.N_EVENT_pp = addedParticles[jscat].N_EVENT_pp;
    tempParticle.N_EVENT_AA = addedParticles[jscat].N_EVENT_AA;
    if( theConfig->isScatt_furtherOfflineParticles() )
      addedParticles.push_back( tempParticle );
  }
  else
  {
    addedParticles[jscat].FLAVOR = F2;
    addedParticles[jscat].m = M2;
    addedParticles[jscat].PX = P2[1];
    addedParticles[jscat].PY = P2[2];
    addedParticles[jscat].PZ = P2[3];
    addedParticles[jscat].E = P2[0];
    ParticleOffline tempParticle = particles_atTimeNow[iscat];
    tempParticle.FLAVOR = F1;
    tempParticle.PX = P1[1];
    tempParticle.PY = P1[2];
    tempParticle.PZ = P1[3];
    tempParticle.E = sqrt( P1[1] * P1[1] + P1[2] * P1[2] + P1[3] * P1[3] );
    tempParticle.unique_id = particles_atTimeNow[iscat].unique_id;
    tempParticle.N_EVENT_pp = addedParticles[jscat].N_EVENT_pp;
    tempParticle.N_EVENT_AA = addedParticles[jscat].N_EVENT_AA;
    if( theConfig->isScatt_furtherOfflineParticles() )
      addedParticles.push_back( tempParticle );
  }
  
  if ( typ == 2240 ) // J/psi + g -> c + cb
  {
    ns_heavy_quarks::jpsi_dissociation++;
    
    // figure out if this is a real event, that is, the ccbar is created from the Jpsi (in every 1/N_test_Jpsi collision). The Jpsi is deleted in any case.
    if( ran2() < 1.0 / theConfig->getJpsiTestparticles() )
    {
      // add second charm quark
      ParticleOffline tempParticle;
      tempParticle.FLAVOR = F1;
      tempParticle.m = Particle::getMass( F1 );
      tempParticle.PX = P1[1];
      tempParticle.PY = P1[2];
      tempParticle.PZ = P1[3];
      tempParticle.E = P1[0];
      tempParticle.T = addedParticles[jscat].T;
      tempParticle.X = addedParticles[jscat].X + dx * ( ran2() - 0.5 );
      tempParticle.Y = addedParticles[jscat].Y + dy * ( ran2() - 0.5 );
      tempParticle.Z = addedParticles[jscat].Z + dx * ( ran2() - 0.5 );
      tempParticle.X_lastInt = tempParticle.X;
      tempParticle.Y_lastInt = tempParticle.Y;
      tempParticle.Z_lastInt = tempParticle.Z;
      tempParticle.T_lastInt = tempParticle.T;
      tempParticle.X_init = tempParticle.X;
      tempParticle.Y_init = tempParticle.Y;
      tempParticle.Z_init = tempParticle.Z;
      tempParticle.PX_init = tempParticle.PX;
      tempParticle.PY_init = tempParticle.PY;
      tempParticle.PZ_init = tempParticle.PZ;
      tempParticle.E_init = -tempParticle.E;   // negative to indicate creation via 2->3 process
      tempParticle.md2g = addedParticles[jscat].md2g;
      tempParticle.md2q = addedParticles[jscat].md2q;
      tempParticle.temperature = addedParticles[jscat].temperature;
      tempParticle.edge = false;
      tempParticle.dead = false;
      tempParticle.free = false;
      tempParticle.init = false;
      tempParticle.initially_produced = false;
      tempParticle.T_creation = tempParticle.T;
      tempParticle.X_traveled = 0.0;
      tempParticle.N_EVENT_pp = addedParticles[jscat].N_EVENT_pp;
      tempParticle.N_EVENT_AA = addedParticles[jscat].N_EVENT_AA;
      tempParticle.cell_id = -1;
      tempParticle.coll_id = -1;
      tempParticle.collisionTime = infinity;
      tempParticle.collisionPartner = -1;
      tempParticle.rate = addedParticles[jscat].rate;    //GeV
      tempParticle.ratev = addedParticles[jscat].ratev;    //GeV
      tempParticle.jpsi_dissociation_number = ns_heavy_quarks::jpsi_dissociation;
      if( addedParticles[jscat].initially_produced )
        tempParticle.N_EVENT_pp = addedParticles[jscat].N_EVENT_pp;
      else
        tempParticle.N_EVENT_pp = addedParticles[jscat].N_EVENT_Cbar;
      tempParticle.unique_id = ParticleOffline::unique_id_counter_added;
      --ParticleOffline::unique_id_counter_added;

      addedParticles.push_back( tempParticle );

      int newIndex = addedParticles.size() - 1;

      cc = ( nexttime - addedParticles[newIndex].T ) / addedParticles[newIndex].E;
      addedParticles[newIndex].T = nexttime;
      addedParticles[newIndex].X = addedParticles[newIndex].X + addedParticles[newIndex].PX * cc;
      addedParticles[newIndex].Y = addedParticles[newIndex].Y + addedParticles[newIndex].PY * cc;
      addedParticles[newIndex].Z = addedParticles[newIndex].Z + addedParticles[newIndex].PZ * cc;
      
      addedParticles[jscat].initially_produced = false;
      addedParticles[jscat].jpsi_dissociation_number = ns_heavy_quarks::jpsi_dissociation;
    }
    else
    {
      // mark Jpsi for removal from global particle list
      deadParticleList.push_back( jscat );
      addedParticles[jscat].dead = true;
      //erase Jpsi from the local list of all particles in this cell
      removeElementFromVector<int>( _allParticlesListAdded, jscat );
      //erase Jpsi from the list of all particles in this cell
      _cellMembersAdded.remove( jscat );
    }
  }
}



void offlineHeavyIonCollision::scatt22_amongAddedParticles_utility( scattering22& scatt22_obj, std::list< int >& _cellMembersAdded, std::vector< int >& _allParticlesListAdded, const int iscat, const int jscat, int& typ, const double nexttime )
{
  FLAVOR_TYPE F1, F2;
  double P1[4], P2[4], R1[4], R2[4];
  double as, Tmax, TT, cc, t;
  double M1, M2;
  double t_hat;
  
  ParticleOffline temp_particle_iscat = addedParticles[iscat];
  ParticleOffline temp_particle_jscat = addedParticles[jscat];

  addedParticles[iscat].getMomentumArray( P1 );
  addedParticles[iscat].getCoordinateArray( R1 );
  F1 = addedParticles[iscat].FLAVOR;
  M1 = addedParticles[iscat].m;

  addedParticles[jscat].getMomentumArray( P2 );
  addedParticles[jscat].getCoordinateArray( R2 );
  F2 = addedParticles[jscat].FLAVOR;
  M2 = addedParticles[jscat].m;

  if ( R1[0] > R2[0] )
    Tmax = R1[0];
  else
    Tmax = R2[0];
  TT = ( nexttime - Tmax ) * ran2() + Tmax;

  cc = ( TT - R1[0] ) / P1[0];
  addedParticles[iscat].T = TT;
  addedParticles[iscat].X = R1[1] + P1[1] * cc;
  addedParticles[iscat].Y = R1[2] + P1[2] * cc;
  addedParticles[iscat].Z = R1[3] + P1[3] * cc;

  addedParticles[jscat].X_traveled += sqrt( pow( P1[1] * cc, 2.0 ) + pow( P1[2] * cc, 2.0 ) + pow( P1[3] * cc, 2.0 ) );

  cc = ( TT - R2[0] ) / P2[0];
  addedParticles[jscat].T = TT;
  addedParticles[jscat].X = R2[1] + P2[1] * cc;
  addedParticles[jscat].Y = R2[2] + P2[2] * cc;
  addedParticles[jscat].Z = R2[3] + P2[3] * cc;

  addedParticles[jscat].X_traveled += sqrt( pow( P2[1] * cc, 2.0 ) + pow( P2[2] * cc, 2.0 ) + pow( P2[3] * cc, 2.0 ) );

  // at this point scattering take place, therefore set properties of last scattering point
  addedParticles[iscat].X_lastInt = addedParticles[iscat].X;
  addedParticles[iscat].Y_lastInt = addedParticles[iscat].Y;
  addedParticles[iscat].Z_lastInt = addedParticles[iscat].Z;
  addedParticles[iscat].T_lastInt = addedParticles[iscat].T;
  
  addedParticles[jscat].X_lastInt = addedParticles[jscat].X;
  addedParticles[jscat].Y_lastInt = addedParticles[jscat].Y;
  addedParticles[jscat].Z_lastInt = addedParticles[jscat].Z;
  addedParticles[jscat].T_lastInt = addedParticles[jscat].T;

  // determine type of scattering, momentum transfer t_hat, new flavor and new masses
  scatt22_obj.getMomentaAndMasses22( F1, F2, M1, M2, t_hat, typ );
  // translate momemtum transfer t_hat into actual momenta of outgoing particles
  scatt22_obj.setNewMomenta22( P1, P2, R1, R2, t_hat );
  
  // figure out if this is a real event, that is, the ccbar is annihilated (in every 1/N_test_Jpsi collision). The Jpsi is created in any case
  if( ran2() < 1.0 / theConfig->getJpsiTestparticles() )
  {
    if ( typ == 2212 || typ == 2213 ) // ccbar->gg or ccbar->qqbar
    {
      ns_heavy_quarks::charmAnnihil++;
      
      //mark charm quark for removal from global particle list
      deadParticleList.push_back( iscat );
      addedParticles[iscat].dead = true;
      //erase charm quark from the local list of all particles in this cell
      removeElementFromVector<int>( _allParticlesListAdded, iscat );
      //erase charm quark from from the list of all particles in this cell
      _cellMembersAdded.remove( iscat );
      
      //mark second charm quark for removal from global particle list
      deadParticleList.push_back( jscat );
      addedParticles[jscat].dead = true;
      //erase charm quark from the local list of all particles in this cell
      removeElementFromVector<int>( _allParticlesListAdded, jscat );
      //erase charm quark from from the list of all particles in this cell
      _cellMembersAdded.remove( jscat );
    }
    else if ( typ == 2241 ) // c + cbar -> J/psi + g
    {
      ns_heavy_quarks::jpsicreation++;
      // set new properties of outgoing particle -> J/psi (is always particle 2)
      addedParticles[iscat].FLAVOR = F2;
      addedParticles[iscat].m = M2;
      addedParticles[iscat].PX = P2[1];
      addedParticles[iscat].PY = P2[2];
      addedParticles[iscat].PZ = P2[3];
      addedParticles[iscat].E = P2[0];
      addedParticles[iscat].initially_produced = false;
      addedParticles[iscat].N_EVENT_Cbar = addedParticles[jscat].N_EVENT_pp;
      addedParticles[iscat].jpsi_dissociation_number = -1; // delete the old property of the charm quark if it was produced via jpsi dissociation before.

      //mark second charm quark for removal from global particle list
      deadParticleList.push_back( jscat );
      addedParticles[jscat].dead = true;
      //erase charm quark from the local list of all particles in this cell
      removeElementFromVector<int>( _allParticlesListAdded, jscat );
      //erase charm quark from from the list of all particles in this cell
      _cellMembersAdded.remove( jscat );
    }
    else
      cout << "error in scatt22AmongAddedPartcl_utility()" << endl;
  }
  else
  {
    if ( typ == 2212 || typ == 2213 ) // ccbar->gg or ccbar->qqbar
    {
      // count this as collision, but do not annihilate charm
      ns_heavy_quarks::charmAnnihil++;
    }
    else if ( typ == 2241 ) // c + cbar -> J/psi + g
    {
      ns_heavy_quarks::jpsicreation++;
      // set new properties of outgoing particle -> J/psi (is always particle 2)
      // add second charm quark
      ParticleOffline tempParticle;
      tempParticle.FLAVOR = F2;
      tempParticle.m = M2;
      tempParticle.PX = P2[1];
      tempParticle.PY = P2[2];
      tempParticle.PZ = P2[3];
      tempParticle.E = P2[0];
      tempParticle.T = addedParticles[iscat].T;
      tempParticle.X = addedParticles[iscat].X;
      tempParticle.Y = addedParticles[iscat].Y;
      tempParticle.Z = addedParticles[iscat].Z;
      tempParticle.X_lastInt = tempParticle.X;
      tempParticle.Y_lastInt = tempParticle.Y;
      tempParticle.Z_lastInt = tempParticle.Z;
      tempParticle.T_lastInt = tempParticle.T;
      tempParticle.X_init = tempParticle.X;
      tempParticle.Y_init = tempParticle.Y;
      tempParticle.Z_init = tempParticle.Z;
      tempParticle.PX_init = tempParticle.PX;
      tempParticle.PY_init = tempParticle.PY;
      tempParticle.PZ_init = tempParticle.PZ;
      tempParticle.E_init = -tempParticle.E;   // negative to indicate creation via 2->3 process
      tempParticle.md2g = addedParticles[iscat].md2g;
      tempParticle.md2q = addedParticles[iscat].md2q;
      tempParticle.temperature = addedParticles[iscat].temperature;
      tempParticle.edge = false;
      tempParticle.dead = false;
      tempParticle.free = false;
      tempParticle.init = false;
      tempParticle.initially_produced = false;
      tempParticle.T_creation = tempParticle.T;
      tempParticle.X_traveled = 0.0;
      tempParticle.N_EVENT_pp = addedParticles[iscat].N_EVENT_pp;
      tempParticle.N_EVENT_AA = addedParticles[iscat].N_EVENT_AA;
      tempParticle.cell_id = -1;
      tempParticle.coll_id = -1;
      tempParticle.collisionTime = infinity;
      tempParticle.collisionPartner = -1;
      tempParticle.rate = addedParticles[iscat].rate;    //GeV
      tempParticle.ratev = addedParticles[iscat].ratev;    //GeV
      tempParticle.jpsi_dissociation_number = -1; // delete the old property of the charm quark if it was produced via jpsi dissociation before.
      tempParticle.N_EVENT_Cbar = addedParticles[iscat].N_EVENT_pp;
      tempParticle.unique_id = ParticleOffline::unique_id_counter_added;
      --ParticleOffline::unique_id_counter_added;

      addedParticles.push_back( tempParticle );
      
      int newIndex = addedParticles.size() - 1;

      cc = ( nexttime - addedParticles[newIndex].T ) / addedParticles[newIndex].E;
      addedParticles[newIndex].T = nexttime;
      addedParticles[newIndex].X = addedParticles[newIndex].X + addedParticles[newIndex].PX * cc;
      addedParticles[newIndex].Y = addedParticles[newIndex].Y + addedParticles[newIndex].PY * cc;
      addedParticles[newIndex].Z = addedParticles[newIndex].Z + addedParticles[newIndex].PZ * cc;
    }
    else
      cout << "error in scatt22AmongAddedPartcl_utility()" << endl;
    
    // reset charm quarks
    addedParticles[iscat] = temp_particle_iscat;
    addedParticles[jscat] = temp_particle_jscat;
  }
}


int offlineHeavyIonCollision::scatt32_offlineWithAddedParticles_utility( scattering32& scatt32_obj, std::list< int >& _cellMembersAdded, std::vector< int >& _allParticlesListAdded, std::vector< int >& _gluonListAdded, const int iscat, const int jscat, const int kscat, int& n32, const double nexttime )
{
  double P1[4], P2[4], P3[4], R1[4], R2[4], R3[4];
  double Tmax, TT, u, phi, cc;
  int order;
  FLAVOR_TYPE F1, F2, F3;
  int typ;

  if( theConfig->doOutput_scatteredMediumParticles() )
  {
    scatteredMediumParticles.push_back( particles_atTimeNow[iscat] );
    scatteredMediumParticles.back().N_EVENT_pp = addedParticles[kscat].N_EVENT_pp;
    scatteredMediumParticles.push_back( particles_atTimeNow[jscat] );
    scatteredMediumParticles.back().N_EVENT_pp = addedParticles[kscat].N_EVENT_pp;
  }

  particles_atTimeNow[iscat].getMomentumArray( P1 );
  particles_atTimeNow[iscat].getCoordinateArray( R1 );
  F1 = particles_atTimeNow[iscat].FLAVOR;

  particles_atTimeNow[jscat].getMomentumArray( P2 );
  particles_atTimeNow[jscat].getCoordinateArray( R2 );
  F2 = particles_atTimeNow[jscat].FLAVOR;

  addedParticles[kscat].getMomentumArray( P3 );
  addedParticles[kscat].getCoordinateArray( R3 );
  F3 = addedParticles[kscat].FLAVOR;
  
  // these routines are not written for heavy flavors
  if( F1 > 2*Particle::max_N_light_flavor || F2 > 2*Particle::max_N_light_flavor || F3 > 2*Particle::max_N_light_flavor )
  {
    std::string errMsg = "Heavy flavor in scatt32_offlineWithAddedParticles_utility.";
    throw eHIC_error( errMsg );
  }

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

  scatt32_obj.getMomenta32( u, phi, typ, F1, F2 );

  scatt32_obj.setNewMomenta32( P1, P2, u, phi );
  order = scatt32_obj.getOrder();

  double pt_out1 = sqrt( pow( P1[1], 2 ) + pow( P1[2], 2 ) );
  double pt_out2 = sqrt( pow( P2[1], 2 ) + pow( P2[2], 2 ) );

  if (( order == 1 ) || ( order == 3 ) )   //123  or  213
  {
//       TODO:// add absorbed gluon to annihilated particle list
//       if ( theConfig->doOutput_scatteredMediumParticles() )
//       {
//         mediumParticles.push_back( addedParticles[kscat] );
//         mediumParticles.back().N_EVENT_pp = addedParticles[kscat].N_EVENT_pp;
//       }

    //mark absorbed gluon for removal from global particle list
    deadParticleList.push_back( kscat );
    addedParticles[kscat].dead = true;

    //erase absorbed gluon from the local list of gluons in this cell
    removeElementFromVector<int>( _gluonListAdded, kscat );
    //erase absorbed gluon from the local list of all particles in this cell
    removeElementFromVector<int>( _allParticlesListAdded, kscat );

    //erase absorbed gluon from the list of all particles in this cell
    _cellMembersAdded.remove( kscat );
    ParticleOffline tempParticle = particles_atTimeNow[iscat];
    tempParticle.FLAVOR = F1;
    tempParticle.PX = P1[1];
    tempParticle.PY = P1[2];
    tempParticle.PZ = P1[3];
    tempParticle.E = sqrt( P1[1] * P1[1] + P1[2] * P1[2] + P1[3] * P1[3] );
    tempParticle.unique_id = particles_atTimeNow[iscat].unique_id;
    tempParticle.N_EVENT_pp = addedParticles[kscat].N_EVENT_pp;
    tempParticle.N_EVENT_AA = addedParticles[kscat].N_EVENT_AA;
    if( theConfig->isScatt_furtherOfflineParticles() )
      addedParticles.push_back( tempParticle );
    ParticleOffline tempParticle2 = particles_atTimeNow[jscat];
    tempParticle2.FLAVOR = F2;
    tempParticle2.PX = P2[1];
    tempParticle2.PY = P2[2];
    tempParticle2.PZ = P2[3];
    tempParticle2.E = sqrt( P2[1] * P2[1] + P2[2] * P2[2] + P2[3] * P2[3] );
    tempParticle2.unique_id = particles_atTimeNow[jscat].unique_id;
    tempParticle2.N_EVENT_pp = addedParticles[kscat].N_EVENT_pp;
    tempParticle2.N_EVENT_AA = addedParticles[kscat].N_EVENT_AA;
    if( theConfig->isScatt_furtherOfflineParticles() )
      addedParticles.push_back( tempParticle2 );
  }
  else if (( order == 2 ) || ( order == 5 ) )  //132  or  312
  {
    if ( pt_out1 > pt_out2 )
    {

//         if jet particle changes to a gluon, it is added to gluon list in cell
      if ( (addedParticles[kscat].FLAVOR != gluon) && (F1 == gluon))
        _gluonListAdded.push_back( kscat );

      addedParticles[kscat].FLAVOR = F1;
      addedParticles[kscat].PX = P1[1];
      addedParticles[kscat].PY = P1[2];
      addedParticles[kscat].PZ = P1[3];
      addedParticles[kscat].E = sqrt( P1[1] * P1[1] + P1[2] * P1[2] + P1[3] * P1[3] );
      ParticleOffline tempParticle = particles_atTimeNow[iscat];
      tempParticle.FLAVOR = F2;
      tempParticle.PX = P2[1];
      tempParticle.PY = P2[2];
      tempParticle.PZ = P2[3];
      tempParticle.E = sqrt( P2[1] * P2[1] + P2[2] * P2[2] + P2[3] * P2[3] );
      tempParticle.unique_id = particles_atTimeNow[iscat].unique_id;
      tempParticle.N_EVENT_pp = addedParticles[kscat].N_EVENT_pp;
      tempParticle.N_EVENT_AA = addedParticles[kscat].N_EVENT_AA;
      if( theConfig->isScatt_furtherOfflineParticles() )
        addedParticles.push_back( tempParticle );
    }
    else
    {
      if( ( addedParticles[kscat].FLAVOR != gluon ) && ( F2 == gluon ) )
        _gluonListAdded.push_back( kscat ); //  if jet particle flavor changes to a gluon, it is added to gluon list in cell

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

      if ( (addedParticles[kscat].FLAVOR != gluon) && (F1 == gluon))
        _gluonListAdded.push_back( kscat ); //  if jet particle flavor changes to a gluon, it is added to gluon list in cell
      
      addedParticles[kscat].FLAVOR = F1;
      addedParticles[kscat].PX = P1[1];
      addedParticles[kscat].PY = P1[2];
      addedParticles[kscat].PZ = P1[3];
      addedParticles[kscat].E = sqrt( P1[1] * P1[1] + P1[2] * P1[2] + P1[3] * P1[3] );
      ParticleOffline tempParticle = particles_atTimeNow[jscat];
      tempParticle.FLAVOR = F2;
      tempParticle.PX = P2[1];
      tempParticle.PY = P2[2];
      tempParticle.PZ = P2[3];
      tempParticle.E = sqrt( P2[1] * P2[1] + P2[2] * P2[2] + P2[3] * P2[3] );
      tempParticle.unique_id = particles_atTimeNow[jscat].unique_id;
      tempParticle.N_EVENT_pp = addedParticles[kscat].N_EVENT_pp;
      tempParticle.N_EVENT_AA = addedParticles[kscat].N_EVENT_AA;
      if( theConfig->isScatt_furtherOfflineParticles() )
        addedParticles.push_back( tempParticle );
    }
    else
    {
      if ( (addedParticles[kscat].FLAVOR != gluon) && (F2 == gluon))
        _gluonListAdded.push_back( kscat ); //  if jet particle flavor changes to a gluon, it is added to gluon list in cell
      
      addedParticles[kscat].FLAVOR = F2;
      addedParticles[kscat].PX = P2[1];
      addedParticles[kscat].PY = P2[2];
      addedParticles[kscat].PZ = P2[3];
      addedParticles[kscat].E = sqrt( P2[1] * P2[1] + P2[2] * P2[2] + P2[3] * P2[3] );
      ParticleOffline tempParticle = particles_atTimeNow[jscat];
      tempParticle.FLAVOR = F1;
      tempParticle.PX = P1[1];
      tempParticle.PY = P1[2];
      tempParticle.PZ = P1[3];
      tempParticle.E = sqrt( P1[1] * P1[1] + P1[2] * P1[2] + P1[3] * P1[3] );
      tempParticle.unique_id = particles_atTimeNow[jscat].unique_id;
      tempParticle.N_EVENT_pp = addedParticles[kscat].N_EVENT_pp;
      tempParticle.N_EVENT_AA = addedParticles[kscat].N_EVENT_AA;
      if( theConfig->isScatt_furtherOfflineParticles() )
        addedParticles.push_back( tempParticle );
    }
  }

  return order;
}




void offlineHeavyIonCollision::jpsi_dissociation_td( const double time )
{
  double *P1, *P2, *P1cm, *P2cm, *P_jpsi, *P_jpsi_cm, *beta_vec, *inv_beta;
  double pp, phi, costheta, sintheta;
  
  P1 = new double[4];
  P2 = new double[4];
  P1cm = new double[4];
  P2cm = new double[4];
  P_jpsi = new double[4];
  P_jpsi_cm = new double[4];
  beta_vec = new double[4];
  inv_beta = new double[4];
  
  for ( int j = 0;j < addedParticles.size();j++ )
  {
    if ( addedParticles[j].FLAVOR == jpsi && addedParticles[j].temperature >= theConfig->getTdJpsi() && FPT_COMP_LE( addedParticles[j].T, time ) && !addedParticles[j].dead ) // Jpsi in cell with temperature higher than Td
    {
      ns_heavy_quarks::jpsi_dissociation++;
      ns_heavy_quarks::jpsi_dissociation_from_temperature++;
      
      // decay Jpsi to charm and anti-charm quarks

      // figure out if this is a real event, that is, a ccbar is created from the Jpsi (in every 1/N_test_Jpsi collision). The Jpsi is deleted in any case.
      if( ran2() < 1.0 / theConfig->getJpsiTestparticles() )
      {
        P_jpsi[0] = addedParticles[j].E;
        P_jpsi[1] = addedParticles[j].PX;
        P_jpsi[2] = addedParticles[j].PY;
        P_jpsi[3] = addedParticles[j].PZ;
        for(int i=1; i<=3; i++)
        {
          beta_vec[i] = ( P_jpsi[i] ) / P_jpsi[0];
          inv_beta[i] = -beta_vec[i];
        }
        
        // boost to cms frame
        lorentz(beta_vec,P_jpsi,P_jpsi_cm);
        
        // conserve energy
        P1cm[0] = P_jpsi_cm[0] / 2.0;
        // momentum
        pp = sqrt( pow(P1cm[0],2.0) - pow( Particle::getMass( charm ) , 2.0 ) );
        if(pp < 0.0)
        {
          cout << "error in jpsi_dissociation_temperature(): " << pp << "  " << P1cm[0] << endl;
          pp = 0.0;
          P1cm[0] = Particle::getMass( charm );
        }
        
        // sample momentum isotropic for first charm quark
        phi = 2.0 * M_PI * ran2();
        costheta = 2.0 * ran2() - 1.0;
        sintheta= sqrt( 1.0 - costheta*costheta );
        P1cm[1] = pp * sintheta * cos(phi);
        P1cm[2] = pp * sintheta * sin(phi);
        P1cm[3] = pp * costheta;
        
        // second charm quark
        P2cm[0] = P1cm[0];
        for(int i=1; i<=3; i++)
          P2cm[i] = -P1cm[i];
        
        // boost both charm quarks back
        lorentz(inv_beta,P1cm,P1);
        lorentz(inv_beta,P2cm,P2);
        
        // update momenta of particles
        addedParticles[j].FLAVOR = charm;
        addedParticles[j].m = Particle::getMass( charm );
        addedParticles[j].PX = P2[1];
        addedParticles[j].PY = P2[2];
        addedParticles[j].PZ = P2[3];
        addedParticles[j].E = P2[0];
        
        
        
        // add second charm quark
        ParticleOffline tempParticle;
        tempParticle.FLAVOR = anti_charm;
        tempParticle.m = Particle::getMass( charm );
        tempParticle.PX = P1[1];
        tempParticle.PY = P1[2];
        tempParticle.PZ = P1[3];
        tempParticle.E = P1[0];
        tempParticle.T = addedParticles[j].T;
        tempParticle.X = addedParticles[j].X + dx * ( ran2() - 0.5 );
        tempParticle.Y = addedParticles[j].Y + dy * ( ran2() - 0.5 );
        tempParticle.Z = addedParticles[j].Z + dx * ( ran2() - 0.5 );
        tempParticle.X_lastInt = tempParticle.X;
        tempParticle.Y_lastInt = tempParticle.Y;
        tempParticle.Z_lastInt = tempParticle.Z;
        tempParticle.T_lastInt = tempParticle.T;
        tempParticle.X_init = tempParticle.X;
        tempParticle.Y_init = tempParticle.Y;
        tempParticle.Z_init = tempParticle.Z;
        tempParticle.PX_init = tempParticle.PX;
        tempParticle.PY_init = tempParticle.PY;
        tempParticle.PZ_init = tempParticle.PZ;
        tempParticle.E_init = -tempParticle.E;   // negative to indicate creation via 2->3 process
        tempParticle.md2g = addedParticles[j].md2g;
        tempParticle.md2q = addedParticles[j].md2q;
        tempParticle.temperature = addedParticles[j].temperature;
        tempParticle.edge = false;
        tempParticle.dead = false;
        tempParticle.free = false;
        tempParticle.init = false;
        tempParticle.initially_produced = false;
        tempParticle.T_creation = tempParticle.T;
        tempParticle.X_traveled = 0.0;
        tempParticle.N_EVENT_pp = addedParticles[j].N_EVENT_pp;
        tempParticle.N_EVENT_AA = addedParticles[j].N_EVENT_AA;
        tempParticle.cell_id = -1;
        tempParticle.coll_id = -1;
        tempParticle.collisionTime = infinity;
        tempParticle.collisionPartner = -1;
        tempParticle.rate = addedParticles[j].rate;    //GeV
        tempParticle.ratev = addedParticles[j].ratev;    //GeV
        tempParticle.jpsi_dissociation_number = ns_heavy_quarks::jpsi_dissociation;
        if( addedParticles[j].initially_produced )
          tempParticle.N_EVENT_pp = addedParticles[j].N_EVENT_pp;
        else
          tempParticle.N_EVENT_pp = addedParticles[j].N_EVENT_Cbar;
        tempParticle.unique_id = ParticleOffline::unique_id_counter_added;
        --ParticleOffline::unique_id_counter_added;

        addedParticles.push_back( tempParticle );
        
        addedParticles[j].initially_produced = false;
        addedParticles[j].jpsi_dissociation_number = ns_heavy_quarks::jpsi_dissociation;

      }
      else
      {
        // mark Jpsi for removal from global particle list
        deadParticleList.push_back( j );
        addedParticles[j].dead = true;
      }
    }
  }
  
  delete[] P1; P1 = NULL;
  delete[] P2; P2 = NULL;
  delete[] P1cm; P1cm = NULL;
  delete[] P2cm; P2cm = NULL;
  delete[] P_jpsi; P_jpsi = NULL;
  delete[] P_jpsi_cm; P_jpsi_cm = NULL;
  delete[] beta_vec; beta_vec = NULL;
  delete[] inv_beta; inv_beta = NULL;
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
void offlineHeavyIonCollision::removeDeadParticles( )
{
  double pt_new;
  int lastIndex = -1;
  deadParticleList.sort();
  while ( addedParticles.size() != 0 && addedParticles.back().dead )
  {
    addedParticles.pop_back();
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

    while ( addedParticles.back().dead )
    {
      addedParticles.pop_back();
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
  double cs22, cs23, I32;
  double R22, R23, R32;
  double as, Vrel, md2g, md2q;
  double lambda, lambdaAvr;
  int iter = 0;
  deque<double> lambdaArray;
  double betaDistEntry;
  double M1, M2, P1P2; 
  
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
  M1 = addedParticles[jetID].m;
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
    
    as = coupling::get_constant_coupling();
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
  
  scattering22 scatt22_object( &theI22 );
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
        M2 = particles_atTimeNow[jscat].m;
        P2[0] = particles_atTimeNow[jscat].E;
        P2[1] = particles_atTimeNow[jscat].PX;
        P2[2] = particles_atTimeNow[jscat].PY;
        P2[3] = particles_atTimeNow[jscat].PZ;
        
        s = pow(( P1[0] + P2[0] ), 2 ) - pow(( P1[1] + P2[1] ), 2 ) - pow(( P1[2] + P2[2] ), 2 ) - pow(( P1[3] + P2[3] ), 2 );
        
        if ( s > 1.1*lambda2 )
        {
          n22++;
          P1P2 = P1[0]*P2[0] - P1[1]*P2[1] - P1[2]*P2[2] - P1[3]*P2[3];
          Vrel = pow((pow(P1P2,2.0)-pow(M1,2.0)*pow(M2,2.0)),0.5)/(P1[0]*P2[0]); // general relative velocity
          as = coupling::get_constant_coupling();
          
          md2g = as * ( addedParticles[jetID].md2g + particles_atTimeNow[jscat].md2g ) / 2.0;
          md2q = as * ( addedParticles[jetID].md2q + particles_atTimeNow[jscat].md2q ) / 2.0;
          // HACK for N_f = 0 background
//           md2g = as * ( addedParticles[jetID].md2g + particles_atTimeNow[jscat].md2g ) / 2.0 * 2;          
//           md2q = md2g * 2.0 / 9.0;
          
          if( theConfig->doScattering_22() )
          {
            scatt22_object.setParameter( P1, P2, F1, F2, M1, M2, s, md2g / as , md2q / as,
                                          theConfig->getKggQQb(), theConfig->getKgQgQ(), theConfig->getKappa_gQgQ(), 
                                          theConfig->isConstantCrossSecGQ(),
                                          theConfig->getConstantCrossSecValueGQ(), theConfig->isIsotropicCrossSecGQ() ); // md2g, md2q are debye masses without the factor alpha_s which is multiplied in scattering22.cpp
            cs22 = scatt22_object.getXSection22();
            probab22 += pow( 0.197, 2.0 ) * cs22 * Vrel * dt / ( dv * testpartcl );
          }
          else
          {
            cs22 = 0.0;
            probab22 = 0.0;
          }
          
          if( theConfig->doScattering_23() )
          {
            n23++;
            lambda_scaled = lambda * sqrt( s );
            
            betaDistEntry = scatt23_object.setParameter( vx, vy, vz, P1, P2, F1, F2 , sqrt( s ), md2g / s, lambda_scaled, as, _gluonList.size() );
            cs23 = scatt23_object.getXSection23( initialStateIndex );   //1/GeV^2
            probab23 += pow( 0.197, 2.0 ) * cs23 * Vrel * dt / ( dv * testpartcl );
          }
          else
          {
            probab23 = 0;
            cs23 = 0;
            lambda_scaled = 0;
          }
        }
      }
    }
    R22 = probab22 / dt * rings[nc].getGamma();
    R23 = probab23 / dt * rings[nc].getGamma();
    //-------------------------------------------------------------
    
    if( theConfig->doScattering_32() )
    {
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
            
            as = coupling::get_constant_coupling();
            
            md2g = as * ( particles_atTimeNow[iscat].md2g + particles_atTimeNow[jscat].md2g + addedParticles[jetID].md2g ) / 3.0;
            md2g = as * ( particles_atTimeNow[iscat].md2q + particles_atTimeNow[jscat].md2q + addedParticles[jetID].md2q ) / 3.0;
            // HACK for N_f = 0 background
  //           md2g = as * ( particles_atTimeNow[iscat].md2g + particles_atTimeNow[jscat].md2g + addedParticles[jetID].md2g ) / 3.0 * 2;
  //           md2q = md2g * 2.0 / 9.0;
            
            // create scattering32 object for the given 3 particles
            betaDistEntry = scatt32_object.setParameter( vx, vy, vz, P1, P2, P3, F1, F2, F3, sqrt( s ), md2g / s, lambda_scaled, as, _gluonList.size() );
            I32 = scatt32_object.getIntegral32_withPrefactors();                        // get the integral I32 for the given 3 particles
            
            probab32 += I32 * dt / ( pow( dv, 2.0 ) * pow( static_cast<double>( testpartcl ) , 2.0 ) );
          }
        }
      }
      else
      {
        //---------------------------- 3->2 ---------------------------
        for ( int m1 = 0; m1 < static_cast<int>( _allParticlesList.size() ) - 1; m1++ )
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
                  
                  as = coupling::get_constant_coupling();
                  
                  md2g = as * ( particles_atTimeNow[iscat].md2g + particles_atTimeNow[jscat].md2g + addedParticles[jetID].md2g ) / 3.0;
                  md2g = as * ( particles_atTimeNow[iscat].md2q + particles_atTimeNow[jscat].md2q + addedParticles[jetID].md2q ) / 3.0;
                  // HACK for N_f = 0 background
  //                 md2g = as * ( particles_atTimeNow[iscat].md2g + particles_atTimeNow[jscat].md2g + addedParticles[jetID].md2g ) / 3.0 * 2;
  //                 md2q = md2g * 2.0 / 9.0;
                  
                  // create scattering32 object for the given 3 particles
                  betaDistEntry = scatt32_object.setParameter( vx, vy, vz, P1, P2, P3, F1, F2, F3, sqrt( s ), md2g / s, lambda_scaled, as, _gluonList.size() );
                  I32 = scatt32_object.getIntegral32_withPrefactors();                        // get the integral I32 for the given 3 particles
                  
                  probab32 += I32 * dt / ( pow( dv, 2.0 ) * pow( static_cast<double>( testpartcl ), 2.0 ) ); 
                }
              }
            }
          }
        }
      }
      R32 = scaleForSelectedTriplets * probab32 / dt * rings[nc].getGamma();  // fm^-1
    }
    else
    {
      R32 = 0.0;
    }
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
