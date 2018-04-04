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
#include <cmath>
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

  long ncoll, ncoll32, ncoll23, ncoll22, ncolle, nemiss12;

  double randomShiftEta, randomShiftX, randomShiftY;
  //  bool dodo1;
  //  int nn_ana1;
  
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
    theI23_massless( false ), theI23_charm_m1( false ), theI23_charm_m2( false ), theI23_bottom_m1( false ), theI23_bottom_m2( false ), // do not load data files right at construction, but after configure() has been called below
    offlineInterface( _offlineInterface ),
    theMFP( _config ),
    theAnalysis( _analysis )
{  
  // if( theConfig->doScattering_12() )
  //   theI12.configure();

  // load 2->2 cross section interpolation data
  if( theConfig->doScattering_22() )
    theI22.configure( theConfig->isCouplingRunning(), Particle::N_light_flavor, Particle::N_heavy_flavor, Particle::Mcharm, Particle::Mbottom, theConfig->getMaxRunningCoupling(), theConfig->getfixedCouplingValue() );
  
  // load 2->3 cross section interpolation data
  if( theConfig->doScattering_23() )
  {
    if( theConfig->getNlightFlavorsAdded() >= 0 )
    {
      theI23_massless.configure( theConfig->I23onlineIntegrationIsSet(), 1, 0.0, theConfig->getKappa23LightPartons(), theConfig->get23GluonFormationTimeTyp(), theConfig->getMatrixElement23(), theConfig->isMd2CounterTermInI23(), theConfig->get23FudgeFactorLpm(), theConfig->getInterpolation23Mode(), theConfig->isMatrixElement23_22qt() );
    }
    if( Particle::N_heavy_flavor > 0 )
    {
      theI23_charm_m1.configure( theConfig->I23onlineIntegrationIsSet(), 1, Particle::Mcharm, theConfig->getKappa23HeavyQuarks(), theConfig->get23GluonFormationTimeTyp(), theConfig->getMatrixElement23(), theConfig->isMd2CounterTermInI23(), theConfig->get23FudgeFactorLpm(), theConfig->getInterpolation23Mode(), theConfig->isMatrixElement23_22qt() );
      theI23_charm_m2.configure( theConfig->I23onlineIntegrationIsSet(), 2, Particle::Mcharm, theConfig->getKappa23HeavyQuarks(), theConfig->get23GluonFormationTimeTyp(), theConfig->getMatrixElement23(), theConfig->isMd2CounterTermInI23(), theConfig->get23FudgeFactorLpm(), theConfig->getInterpolation23Mode(), theConfig->isMatrixElement23_22qt() );
    }
    if( Particle::N_heavy_flavor > 1 )
    {
      theI23_bottom_m1.configure( theConfig->I23onlineIntegrationIsSet(), 1, Particle::Mbottom, theConfig->getKappa23HeavyQuarks(), theConfig->get23GluonFormationTimeTyp(), theConfig->getMatrixElement23(), theConfig->isMd2CounterTermInI23(), theConfig->get23FudgeFactorLpm(), theConfig->getInterpolation23Mode(), theConfig->isMatrixElement23_22qt() );
      theI23_bottom_m2.configure( theConfig->I23onlineIntegrationIsSet(), 2, Particle::Mbottom, theConfig->getKappa23HeavyQuarks(), theConfig->get23GluonFormationTimeTyp(), theConfig->getMatrixElement23(), theConfig->isMd2CounterTermInI23(), theConfig->get23FudgeFactorLpm(), theConfig->getInterpolation23Mode(), theConfig->isMatrixElement23_22qt() );
      
    }
    
    if( theConfig->getJetMfpComputationType() == computeMfpInterpolation || theConfig->getJetMfpComputationType() == thermalMfpGluon )
      theMFP.loadData(); 
  }
  
  nGet23Errors = 0;
  nGet32Errors = 0;
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
    for(unsigned int j = 0; j < addedParticles.size(); j++ )
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
      for( unsigned int j = 0; j < addedParticles.size(); j++ )
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
    for( unsigned int j = 0; j < addedParticles.size(); j++ )
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
    for( unsigned int j = 0; j < addedParticles.size(); j++ )
      sum += addedParticles[j].Mom.E();
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
  int nemiss12_backup = 0;
  int charmAnnihil_backup = 0;
  int jpsicreation_backup = 0;
  int jpsi_dissociation_from_temperature_backup = 0;
  int jpsi_dissociation_backup = 0;
  
  list<int> edgeCellCopy, edgeCellAddedCopy;

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
    addedParticlesCopy = ppHadronization_hq.heavyQuarkFragmentation( addedParticles );
  }
  if( theConfig->isMesonDecay() )
  {
#ifdef Pythia_FOUND
    mesonDecay ppMesonDecay( theConfig->getNumberElectronStat(), theConfig->isMuonsInsteadOfElectrons(), theConfig->isStudyNonPromptJpsiInsteadOfElectrons() );
    addedPartcl_electron = ppMesonDecay.decayToElectronsPythia( addedParticlesCopy );
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
  for ( unsigned int i = 0;i < addedParticles.size(); i++ )
  {
    if ( addedParticles[i].Pos.T() <= simulationTime )
    {
      addedParticles[i].Propagate( simulationTime, addedParticles[i].X_traveled );
      theAnalysis->addJetEvent_initial( i );
    }
  }

  int n_dt = 0;
  double dt_sum = 0.0;
  int n_again = 0;

  do
  {
    // remove init tag from particles after their formation
    for ( unsigned int i = 0; i < addedParticles.size(); i++ )
    {
      if ( addedParticles[i].Pos.T() <= simulationTime )
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
    nemiss12_backup = nemiss12;
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
      nemiss12 = nemiss12_backup;
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
      theAnalysis->registerProgressInformationForOutput( simulationTime, dt, addedParticles.size(), particles_atTimeNow.size(), ncoll, ncoll22, ncoll23, ncoll32, nemiss12 );
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
    addedParticlesCopy = theHadronization_hq.heavyQuarkFragmentation( addedParticles );
  }
  if( theConfig->isMesonDecay() )
  {
#ifdef Pythia_FOUND
    mesonDecay theMesonDecay( theConfig->getNumberElectronStat(), theConfig->isMuonsInsteadOfElectrons(), theConfig->isStudyNonPromptJpsiInsteadOfElectrons() );
    addedPartcl_electron = theMesonDecay.decayToElectronsPythia( addedParticlesCopy );
#else
    string errMsg( "Could not perform decay of heavy mesons to electron because PYTHIA was not found." );
    throw eHIC_error( errMsg );
#endif
  }

  theAnalysis->finalOutput( stoptime );
  theAnalysis->addJetEvents_final();
  
  
  cout << "number of errors in get32(...) = " << nGet32Errors << endl;
  cout << "number of errors in get23(...) = " << nGet23Errors << endl;
  
  // List particle numbers for all flavors
  cout << "==========================" << endl;
  cout << "Added particles:" << endl;
  cout << "# of all=" << addedParticles.size() << endl;
  int fl_sum = 0;
  cout << "Flavor       Number    Number wo testparticles per Event" << endl;
  for( int i = 0; i <= 10; i++ )
  {
    for( unsigned int j = 0; j < addedParticles.size(); j++ )
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
  for( unsigned int j = 0; j < addedParticles.size(); j++ )
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
  for( unsigned int j = 0; j < addedParticles.size(); j++ )
    sum += addedParticles[j].Mom.E();
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
  VectorEPxPyPz Mom1, Mom2;
  FLAVOR_TYPE F1, F2, F3;

  // give partcl from cascade the values for cell structurement which could have changed in collisions()
  for ( unsigned int k = 0; k < particles_atTimeNow.size(); k++ )
  {
    particlesEvolving[k].init = particles_atTimeNow[k].init;
    particlesEvolving[k].edge = particles_atTimeNow[k].edge;
    particlesEvolving[k].free = particles_atTimeNow[k].free;
    particlesEvolving[k].isAlreadyInAddedParticles = particles_atTimeNow[k].isAlreadyInAddedParticles;
  }

  if ( evolveToTime <= stoptime_last )
  {
    for ( int i = 0; i < theConfig->getN_init(); i++ )
    {
      particlesEvolving[i].init = true;
      particlesEvolving[i].Pos = particles_init[i].Pos;
      particlesEvolving[i].Old = particlesEvolving[i].Mom = particles_init[i].Mom;
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
      Mom1 = ptrInteraction22->Mom1;
      Mom2 = ptrInteraction22->Mom2;
      F1 = ptrInteraction22->F1;
      F2 = ptrInteraction22->F2;
      
      if ( time <= ( evolveToTime + 1.0e-6 ) )
      {
        particlesEvolving[iscat].init = particlesEvolving[jscat].init = false;

        particlesEvolving[iscat].Propagate( time );
        particlesEvolving[jscat].Propagate( time );

        particlesEvolving[iscat].FLAVOR = static_cast<FLAVOR_TYPE>( F1 );
        particlesEvolving[iscat].Mom = Mom1;
        particlesEvolving[iscat].Mom.E() = sqrt( Mom1.vec2() ); // necessary ???

        particlesEvolving[jscat].FLAVOR = static_cast<FLAVOR_TYPE>( F2 );
        particlesEvolving[jscat].Mom = Mom2;
        particlesEvolving[jscat].Mom.E() = sqrt( Mom2.vec2() ); // necessary ???

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
      Mom1 = ptrInteraction23->Mom1;
      Mom2 = ptrInteraction23->Mom2;
      F1 = static_cast<FLAVOR_TYPE>( ptrInteraction23->F1 );
      F2 = static_cast<FLAVOR_TYPE>( ptrInteraction23->F2 );
      F3 = static_cast<FLAVOR_TYPE>( ptrInteraction23->F3 );
      particlesEvolving[kscat].Pos = ptrInteraction23->Pos3;
      particlesEvolving[kscat].Mom = ptrInteraction23->Mom3;
      particlesEvolving[kscat].unique_id = Particle::unique_id_counter;
      Particle::unique_id_counter++;

      if ( time <= evolveToTime + 1.0e-6 )
      {
        particlesEvolving[kscat].Pos.T() = time;
        particlesEvolving[kscat].FLAVOR = static_cast<FLAVOR_TYPE>( F3 );
        particlesEvolving[kscat].Mom.E() = sqrt( particlesEvolving[kscat].Mom.vec2() );
        particlesEvolving[kscat].free = false;

        particlesEvolving[iscat].init = particlesEvolving[jscat].init = particlesEvolving[kscat].init = false;

        particlesEvolving[iscat].Propagate( time );
        particlesEvolving[jscat].Propagate( time );

        particlesEvolving[iscat].FLAVOR = static_cast<FLAVOR_TYPE>( F1 );
        particlesEvolving[iscat].Mom = Mom1;
        particlesEvolving[iscat].Mom.E() = sqrt( Mom1.vec2() ); // necessary ???

        particlesEvolving[jscat].FLAVOR = static_cast<FLAVOR_TYPE>( F2 );
        particlesEvolving[jscat].Mom = Mom2;
        particlesEvolving[jscat].Mom.E() = sqrt( Mom2.vec2() ); // necessary ???

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
      Mom1 = ptrInteraction32->Mom1;
      Mom2 = ptrInteraction32->Mom2;
      F1 = ptrInteraction32->F1;
      F2 = ptrInteraction32->F2;

      if ( time <= evolveToTime + 1.0e-6 )
      {
        particlesEvolving[iscat].init = particlesEvolving[jscat].init = false;

        particlesEvolving[iscat].Propagate( time );
        particlesEvolving[jscat].Propagate( time );

        particlesEvolving[iscat].FLAVOR = static_cast<FLAVOR_TYPE>( F1 );
        particlesEvolving[iscat].Mom = Mom1;
        particlesEvolving[iscat].Mom.E() = sqrt( Mom1.vec2() ); // necessary ???

        particlesEvolving[jscat].FLAVOR = static_cast<FLAVOR_TYPE>( F2 );
        particlesEvolving[jscat].Mom = Mom2;
        particlesEvolving[jscat].Mom.E() = sqrt( Mom2.vec2() ); // necessary ???
        
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
      Mom1 = ptrInteractionElastic->Mom1;
      Mom2 = ptrInteractionElastic->Mom2;
      
      if (( timei <= evolveToTime + 1.0e-6 ) || ( timej <= evolveToTime + 1.0e-6 ) )
      {
        particlesEvolving[iscat].init = particlesEvolving[jscat].init = false;

        particlesEvolving[iscat].Old = particlesEvolving[iscat].Mom;
        particlesEvolving[jscat].Old = particlesEvolving[jscat].Mom;

        particlesEvolving[iscat].Propagate( timei );
        particlesEvolving[iscat].Mom = Mom1;
        particlesEvolving[iscat].Mom.E() = sqrt( Mom1.vec2() ); // necessary ???

        particlesEvolving[jscat].Propagate( timej );
        particlesEvolving[jscat].Mom = Mom2;
        particlesEvolving[jscat].Mom.E() = sqrt( Mom2.vec2() ); // necessary ???

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
    if ( particlesEvolving[i].Pos.T() < evolveToTime + 1.0e-8 )
    {
      particlesEvolving[i].init = false;
    }
  }

  // duplicate partcl from cascade to partclAtTimeNow
  particles_atTimeNow = particlesEvolving;
  
  // particles vector is larger and in size constant. But consider for partclAtTimenow only actual physical present particles (number is numberEvolvingParticles)
  particles_atTimeNow.resize( numberEvolvingParticles ); 

  // propagate particles_atTimeNow to current time
  for ( unsigned int i = 0; i < particles_atTimeNow.size(); i++ )
  {
    if ( particles_atTimeNow[i].Pos.T() < evolveToTime )
    {
      particles_atTimeNow[i].Propagate( evolveToTime );
    }
  }

  stoptime_last = evolveToTime;
  return dt_cascade;
}




void offlineHeavyIonCollision::cell_ID( double _time )
{
  int nx, ny, nz, cell_id;
  double halfsize, eta;

  halfsize = transLen / 2.0;
  int centralEtaIndex = static_cast<int>( etaBins.size() ) / 2;

  int NinAct, Nfree1, NinFormInit, NinFormGeom;
  NinAct = Nfree1 = NinFormInit = NinFormGeom = 0;

  formGeom.clear();
  for ( unsigned int i = 0; i < cells.size(); i++ )
  {
    cells[i].particleList.clear();
  }

  for ( unsigned int i = 0; i < particles_atTimeNow.size(); i++ )
  {
    if ( particles_atTimeNow[i].Pos.T() < _time )
    {
      NinAct++;
      particles_atTimeNow[i].init = false;

      if ( fabs( particles_atTimeNow[i].Pos.X() - randomShiftX - halfsize ) < 1.0e-6 )
      {
        nx = IX - 1;
      }
      else
      {
        nx = static_cast<int>((( particles_atTimeNow[i].Pos.X() - randomShiftX ) / transLen + 0.5 ) * IX );
      }

      if ( fabs( particles_atTimeNow[i].Pos.Y() - randomShiftY - halfsize ) < 1.0e-6 )
      {
        ny = IY - 1;
      }
      else
      {
        ny = static_cast<int>((( particles_atTimeNow[i].Pos.Y() - randomShiftY ) / transLen + 0.5 ) * IY );
      }

      eta = particles_atTimeNow[i].Pos.Rapidity();
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
  for ( unsigned int i = 0; i < cells.size(); i++ )
  {
    cellsAdded[i].particleList.clear();
  }
  formGeomAdded.clear();

  for ( unsigned int i = 0; i < addedParticles.size(); i++ )
  {
    if ( addedParticles[i].Pos.T() < _time )
    {
      addedParticles[i].init = false;

      if ( fabs( addedParticles[i].Pos.X() - randomShiftX - halfsize ) < 1.0e-6 )
      {
        nx = IX - 1;
      }
      else
      {
        nx = static_cast<int>((( addedParticles[i].Pos.X() - randomShiftX ) / transLen + 0.5 ) * IX );
      }

      if ( fabs( addedParticles[i].Pos.Y() - randomShiftY - halfsize ) < 1.0e-6 )
      {
        ny = IY - 1;
      }
      else
      {
        ny = static_cast<int>((( addedParticles[i].Pos.Y() - randomShiftY ) / transLen + 0.5 ) * IY );
      }

      eta = addedParticles[i].Pos.Rapidity();
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
  
  int nGluons = 0;
  int nGluonsAdded = 0;
  
  int IXY, id, nc;
  
  double dz, cc, zz, eta;
  bool free;
  list<int> formGeomCopy;
  vector<int> gluonList, allParticlesList;
  vector<int> gluonListAdded, allParticlesListAdded;
  
  again = false;
  
  IXY = IX * IY;
  
  formGeomCopy = formGeom;

  if ( !formGeomCopy.empty() )
  {
    list<int>::iterator iIt;
    int id = -1;
    for ( iIt = formGeomCopy.begin(); iIt != formGeomCopy.end(); )
    {
      id = *iIt;
      cc = ( timenow - particles_atTimeNow[id].Pos.T() ) / sqrt( particles_atTimeNow[id].Old.vec2() + pow( particles_atTimeNow[id].m, 2 ) );
      zz = particles_atTimeNow[id].Pos.Z() + particles_atTimeNow[id].Old.Pz() * cc;
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
        xt = particles_atTimeNow[( *iIt )].Pos.Perp();

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
      
      cc = ( timenow - particles_atTimeNow[id].Pos.T() ) / sqrt( particles_atTimeNow[id].Old.vec2() + pow( particles_atTimeNow[id].m, 2 ) );
      zz = particles_atTimeNow[id].Pos.Z() + particles_atTimeNow[id].Old.Pz() * cc;
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
          xt = particles_atTimeNow[id].Pos.Perp();
          
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
          xt = addedParticles[id].Pos.Perp();
          
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
            addedParticles[id].Propagate( nexttime );
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
            
        //------------------------- calculate flow of this cell --------------------
            // VectorEPxPyPz p_cell = VectorEPxPyPz( 0.0, 0.0, 0.0, 0.0 );
            // list<int>::const_iterator iIt;
            // for ( iIt = cells[j].particleList.begin(); iIt != cells[j].particleList.end(); iIt++ )
            // {
            //   id = *iIt;
            //   p_cell += particles_atTimeNow[id].Mom;
            // }
            // const VectorEPxPyPz v_cell = (cells[j].particleList.size() >= 5)? p_cell.NormalizeToE() : VectorEPxPyPz( 1.0, 0.0, 0.0, 0.0 );

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
              
              xt = particles_atTimeNow[id].Pos.Perp();
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
              
              xt = addedParticles[id].Pos.Perp();
              nc = rings.getIndex( xt );
              
              addedParticles[id].md2g = rings[nc].getAveraged_md2g();
              addedParticles[id].md2q = rings[nc].getAveraged_md2q();
              addedParticles[id].temperature = rings[nc].getEffectiveTemperature();
              addedParticles[id].temperatureAMY = rings[nc].getEffectiveTemperatureAMY();
              addedParticles[id].v_cell = rings[nc].getAveraged_v();
              
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
            
            // calculate actual mean free path of added particles within this timestep in order to use it for both 3->2 and 2->3 processes
            if ( theConfig->getJetMfpComputationType() == computeMfpLastTimestep && ( theConfig->doScattering_23() || theConfig->doScattering_32() ) )
            {
              for ( iIt = cellsAdded[j].particleList.begin(); iIt != cellsAdded[j].particleList.end(); iIt++ )
              {
                int id = *iIt;
                
                if ( addedParticles[id].Mom.Perp() > theConfig->getMinimumPT() )
                {
                  addedParticles[id].lambda_added_old = addedParticles[id].lambda_added;
                  
                  double xt = addedParticles[id].Pos.Perp();
                  int ringIndex = rings.getIndex( xt );
                  
                  if(  addedParticles[id].lambda_added_old <= 0.0 || addedParticles[id].rate_added <= 0.0 )
                    addedParticles[id].lambda_added = iterate_mfp_bisection( allParticlesList, gluonList, id, dt, dv, addedParticles[id].lambda_added_old ); // fm
                  else
                  {
                    // velocity of cell
                    VectorXYZ v_cell = rings[ringIndex].getAveraged_v();
                
                    // velocity of added particle in lab frame
                    VectorTXYZ v_jet = addedParticles[id].Mom * (1.0/addedParticles[id].Mom.E());
          //           double velocity_lab = sqrt( 1.0 - pow( addedParticles[id].m ,2.0) / pow( addedParticles[id].E ,2.0) );
                    
                    // Compute velocity of added particle in rest frame of fluid. The general expression for adding two velocities is not symmetric in v1 and v2. Therefore compute both cases and take average.
                    double velocity_rest_1 = addVelocities( v_cell.X(), v_cell.Y(), v_cell.Z(), v_jet.X(), v_jet.Y(), v_jet.Z() );
                    double velocity_rest_2 = addVelocities( v_jet.X(), v_jet.Y(), v_jet.Z(), v_cell.X(), v_cell.Y(), v_cell.Z() );
                    double velocity_rest = ( velocity_rest_1 + velocity_rest_2 ) / 2.0;
                    
                    // mean free path in rest frame. Consequently the velocity in the rest frame is needed.
                    addedParticles[id].lambda_added = velocity_rest / addedParticles[id].rate_added; // fm
                  }
                }
              }
            }

            int n32 = 0;

            if( theConfig->doScattering_12() )
            {
              emission12_offlineParticles( allParticlesListAdded, again, nexttime );
              if( again && false )
              {
                return;
              }
            }
            
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
                
                xt = particles_atTimeNow[id].Pos.Perp();
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
                
                xt = addedParticles[id].Pos.Perp();
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
            addedParticles[iscat].Propagate( nexttime );
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
  double s, csgg, cs23, cs22, Vrel, lambda_scaled;
  double M1, M2;
  double probab22, probab23, probab2322;
  double averagedRate;
  double rate_added_sum; // 1/fm
  double xt;
  double betaDistEntry;
  double md2g_wo_as, md2q_wo_as;
  int ringIndex;
  int initialStateIndex = -1;
  FLAVOR_TYPE F1, F2;


  const int nTotal = _cell.particleList.size();
//   vector<int> nQuarks( Particle::N_light_flavor, 0 );
//   vector<int> nAntiQuarks( Particle::N_light_flavor, 0 );

  scattering23 scatt23_object( &theI23_massless, &theI23_charm_m1, &theI23_charm_m2, &theI23_bottom_m1, &theI23_bottom_m2 );
  scattering22 scatt22_object( &theI22 );
  
  double lambda = 0; // fm
  double pt_addedParticle = 0;

  //   const int nAllQuarks = std::accumulate( nQuarks.begin(), nQuarks.end(), 0 ) + std::accumulate( nAntiQuarks.begin(), nAntiQuarks.end(), 0 );
//   const int allPairs = binomial( nTotal, 2 );
//   const int consideredPairs = 25;

  for ( int j = 0; j < static_cast<int>( _allParticlesListAdded.size() ); j++ )
  {
    jscat = _allParticlesListAdded[j];
    
    rate_added_sum = 0;

    pt_addedParticle = addedParticles[jscat].Mom.Perp();
    const double tau = sqrt( addedParticles[jscat].Pos.T() * addedParticles[jscat].Pos.T() - addedParticles[jscat].Pos.Z() * addedParticles[jscat].Pos.Z() );

    if( pt_addedParticle < theConfig->getMinimumPT() || addedParticles[jscat].dead || tau < theConfig->getTauMin() )
    {
      continue; // jump to next particle in the list
    }
    
    if( theConfig->doScattering_23() )
    {
      // calculate their actual mean free path
      if( nTotal > 0 ) // if there are actually medium particles in the cell. Also: think about whether calculating a mean free path with just a few other particle in a cell makes sense or if a larger cut off must be implemented
      {
        xt = addedParticles[jscat].Pos.Perp();
        ringIndex = rings.getIndex( xt );
        
        switch(theConfig->getJetMfpComputationType())
        {
          case 0: // computeMfpLastTimestep
            if( addedParticles[jscat].lambda_added > 0.0 && addedParticles[jscat].lambda_added_old > 0.0 )
              lambda = ( addedParticles[jscat].lambda_added + addedParticles[jscat].lambda_added_old ) / 2.0 ; // mean free path of addedParticles in fm
            else if( addedParticles[jscat].lambda_added > 0.0 )
              lambda = addedParticles[jscat].lambda_added; // mean free path of addedParticles in fm
            else
            {
              cout << "error in scattOfflinePartclWithAddedPartcl: lambda negative: lambda = " << addedParticles[jscat].lambda_added << "  " << addedParticles[jscat].lambda_added_old << endl;
              cout << addedParticles[jscat].rate_added << "  " << iterate_mfp_bisection( _allParticlesList, _gluonList, jscat, dt, dv, addedParticles[jscat].lambda_added_old ) << endl;
              std::string errMsg = "Error in scattOfflinePartclWithAddedPartcl: lambda negative.";
              throw eHIC_error( errMsg );
            }
            break;
          case 1: // computeMfpIteration
            lambda = iterate_mfp_bisection( _allParticlesList, _gluonList, jscat, dt, dv, addedParticles[jscat].lambda_added_old ); // fm
//             lambda = iterateMFP( _allParticlesList, _gluonList, jscat, dt, dv ); //fm
            addedParticles[jscat].lambda_added = lambda;
            break;
          case 2: // computeMfpInterpolation
            lambda = theMFP.getMeanFreePath( addedParticles[jscat].Mom.E(), addedParticles[jscat].FLAVOR, rings[ringIndex].getEffectiveTemperature(), rings[ringIndex].getGluonDensity(), rings[ringIndex].getQuarkDensity(), fm );
            break;
          case 3: // fixedMfp
            lambda = theConfig->getFixedMfpAdded();
            break;
          case 4: // thermalMfpGluon
            lambda = theMFP.getMeanFreePath( 3.0 * rings[ringIndex].getEffectiveTemperature(), gluon, rings[ringIndex].getEffectiveTemperature(), rings[ringIndex].getGluonDensity(), rings[ringIndex].getQuarkDensity(), fm );
            break;
          default: 
            std::string errMsg = "Error in scattOfflinePartclWithAddedPartcl: wrong mfp determination type.";
            throw eHIC_error( errMsg );
        }
      }
      else
      {
        std::string errMsg = "Error in scattOfflinePartclWithAddedPartcl: not enough particles in cells.";
        throw eHIC_error( errMsg );
      }
    }

    for ( int i = 0; i < static_cast<int>( _allParticlesList.size() ); i++ )
    {
      iscat = _allParticlesList[i];
      
      if ( particles_atTimeNow[iscat].dead || addedParticles[jscat].dead ) // The last statement to check again whether jscat is dead is necessary since jscat could be deleted in the scattering (Jpsi dissociation which is not "real" since we use more testparticles for Jpsi in which case the Jpsi is deleted) before but the scatterings with all remaining jscat would be carried out.
      {
        continue; // jump to next particle from _cell.particleList
      }
      
      F1 = particles_atTimeNow[iscat].FLAVOR;
      M1 = particles_atTimeNow[iscat].m;

      F2 = addedParticles[jscat].FLAVOR;
      M2 = addedParticles[jscat].m;

      s = (particles_atTimeNow[iscat].Mom + addedParticles[jscat].Mom).M2();

      xt = ( particles_atTimeNow[iscat].Pos.Perp() + addedParticles[jscat].Pos.Perp() ) / 2;
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
        Vrel = VelRel(particles_atTimeNow[iscat].Mom, addedParticles[jscat].Mom, M1, M2);

        //factor as (alpha_s) is not included in definitions of partcl[jscat].md2g, it will be multiplied in the scattering routines
        md2g_wo_as = ( particles_atTimeNow[iscat].md2g + addedParticles[jscat].md2g ) / 2.0;
        md2q_wo_as = ( particles_atTimeNow[iscat].md2q + addedParticles[jscat].md2q ) / 2.0;
        
        if( theConfig->doScattering_22() )
        {
          scatt22_object.setParameter( particles_atTimeNow[iscat].Mom, addedParticles[jscat].Mom,
                                       F1, F2, M1, M2, s, md2g_wo_as , md2q_wo_as,
                                      theConfig->getKggQQb(), theConfig->getKgQgQ(), theConfig->getKappa_gQgQ(), 
                                      theConfig->isConstantCrossSecGQ(),
                                      theConfig->getConstantCrossSecValueGQ(), theConfig->isIsotropicCrossSecGQ(), theConfig->getKfactor_light() ); // md2g, md2q are debye masses without the factor alpha_s which is multiplied in scattering22.cpp
          
          switch( theConfig->getCrossSectionMethod() )
          {
            case csMethod_pQCD:
              cs22 = scatt22_object.getXSection22( initialStateIndex );
              break;
            case csMethod_constCS:
              cs22 = theConfig->getInputCrossSectionValue() / pow(0.197,2) / 10.0;//1/GeV^2
              break;
            default:
              string errMsg = "Unknown cross-section type in scatt2223_offlineWithAddedParticles... Unrecoverable error!";
              throw eHIC_error( errMsg );
          }
          
          probab22 = pow( 0.197, 2.0 ) * cs22 * Vrel * dt / ( dv * testpartcl );
        }
        else
        {
          cs22 = probab22 = 0.0;
        }
        

        if( theConfig->doScattering_23() )
        {
          if ( lambda > 0 )
          {
            lambda_scaled = lambda * sqrt( s ) / 0.197; // lambda in fm, sqrt(s) in GeV, lambda_scaled dimensionless
          }
          else
          {
            if ( averagedRate > epsilon )
            {
              lambda_scaled = sqrt( s ) / averagedRate;  //dimensionless
            }
            else
            {
              xsection_gg_gg csObj( s, md2g_wo_as, md2q_wo_as, &theI22, theConfig->getKfactor_light() );
              csgg = csObj.totalCrossSection();
              lambda_scaled = ( dv * rings[ringIndex].getGamma() * testpartcl * sqrt( s ) ) / ( nTotal * csgg * pow( 0.197, 3.0 ) );
            }
          }
          
          betaDistEntry = scatt23_object.setParameter( rings[ringIndex].getAveraged_v(), 
                                                       particles_atTimeNow[iscat].Mom, addedParticles[jscat].Mom,
                        F1, F2, M1, M2, sqrt( s ), 
                        md2g_wo_as / s, lambda_scaled, 
                        theConfig->getK23LightPartons(), theConfig->getK23HeavyQuarks(),
                        theConfig->getKappa23LightPartons(), theConfig->getKappa23HeavyQuarks(),
                        theConfig->I23onlineIntegrationIsSet(),
                        theConfig->get23GluonFormationTimeTyp(), theConfig->getMatrixElement23(), theConfig->isMd2CounterTermInI23(), theConfig->get23FudgeFactorLpm(), _gluonList.size() );   
          
          cs23 = scatt23_object.getXSection23( initialStateIndex ); //1/GeV^2

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
          else if ( ParticleOffline::mapToGenericFlavorType( F2 ) == light_quark || ParticleOffline::mapToGenericFlavorType( F2 ) == anti_light_quark )
          {
            ++n23_collected_quark;
            p23_collected_quark += probab23;
            lambdaJet_quark += lambda_scaled;
          }
        }

        if ( cs22 > 0.0 )
        {
          ++_cellAdded.nCollected22;
          _cellAdded.md2g_scaled_22 += md2g_wo_as / s; // check!!!
          _cellAdded.md2q_scaled_22 += md2q_wo_as / s; // check!!!
        }
        if ( cs23 > 0.0 )
        {
          ++_cellAdded.nCollected23;
          _cellAdded.md2g_scaled_23 += md2g_wo_as / s; // check!!!
          _cellAdded.md2q_scaled_23 += md2q_wo_as / s; // check!!!
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
      
      if ( probab22 > 0.0 )
      {
        rate_added_sum += probab22 / dt * rings[ringIndex].getGamma(); // fm^-1
      }
      
      if ( probab23 > 0.0 )
      {
        rate_added_sum += probab23 / dt * rings[ringIndex].getGamma(); // fm^-1
      }
      
      if ( ran2() < probab2322 )
      {
        double pt_jscat = addedParticles[jscat].Mom.Perp();
        double pt_nmb;

        if ( ran2() * probab2322 < probab23 )
        {
          int jetEventIndex = -1;
          if( pt_jscat > theAnalysis->getJetTracking_PT() )
          {
            jetEventIndex = theAnalysis->addJetEvent_in( iscat, -1, jscat, c2to3, cs23, _cell.index, lambda_scaled / sqrt( s ) );
          }

          int newIndex = scatt23_offlineWithAddedParticles_utility( scatt23_object, _cell, iscat, jscat, again, nexttime );

          if( newIndex >= 0 )
          {
            pt_jscat = addedParticles[jscat].Mom.Perp();
            pt_nmb =addedParticles[newIndex].Mom.Perp();
            if( jetEventIndex != -1 || pt_jscat > theAnalysis->getJetTracking_PT() || pt_nmb > theAnalysis->getJetTracking_PT() )
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

          pt_jscat = addedParticles[jscat].Mom.Perp();
          if( jetEventIndex != -1 || pt_jscat > theAnalysis->getJetTracking_PT() )
          {
            theAnalysis->addJetEvent_out( jetEventIndex, jscat, iscat, -1, c2to2 );
          }
        }
      }
    }
    
    // add 3->2 rate of this timestep calculated in scatt32_offlineWithAddedParticles and reset it afterwards
    rate_added_sum += addedParticles[jscat].rate_added_32;
    addedParticles[jscat].rate_added_32 = 0.0;
    
    const double epsilon_rate_added = 1E-3; // 0.001 fm^-1 which corresponds to a lambda = 1000 fm
    // set new rate
    // If the cross sections are very small or 0, R22 = R23 = 0 or very small which causes lambda to be very large. Therefore, epsilon serves as a cut-off.  However, since the cross section is small probably no scattering would take place anyhow. Return a negative value which forces that the mean free path will be calculated iteratively.
    if( rate_added_sum < epsilon_rate_added ) 
      addedParticles[jscat].rate_added = -1.0;
    else
      addedParticles[jscat].rate_added = rate_added_sum; // fm^-1
  }

  list<int>::const_iterator iIt;
  for ( iIt =  _cellAdded.particleList.begin(); iIt != _cellAdded.particleList.end(); iIt++ )
  {
    iscat = *iIt;
    addedParticles[iscat].Propagate( nexttime, addedParticles[iscat].X_traveled );
  }
}

void offlineHeavyIonCollision::emission12_offlineParticles( std::vector< int >& _allParticlesListAdded, bool& again, const double nexttime )
{
  int jscat;
  double xt;
  int ringIndex;
  FLAVOR_TYPE F1, F2;

  emission12 emission12_object( &theI12 );

  double pt_addedParticle = 0;

  for ( int j = 0; j < static_cast<int>( _allParticlesListAdded.size() ); j++ )
  {
    jscat = _allParticlesListAdded[j];
    
    pt_addedParticle = addedParticles[jscat].Mom.Perp();
    const double tau = sqrt( addedParticles[jscat].Pos.T() * addedParticles[jscat].Pos.T() - addedParticles[jscat].Pos.Z() * addedParticles[jscat].Pos.Z() );

    if( pt_addedParticle < theConfig->getMinimumPT() || addedParticles[jscat].dead || tau < theConfig->getTauMin() )
    {
      continue; // jump to next particle in the list
    }
    
    F1 = addedParticles[jscat].FLAVOR;
    xt = ( addedParticles[jscat].Pos.Perp() ) / 2;
    ringIndex = rings.getIndex( xt );

    const double T = addedParticles[jscat].temperatureAMY;
    const VectorTXYZ v_cell = addedParticles[jscat].v_cell;
        
    emission12_object.setParameter( addedParticles[jscat].Mom, F1, T, v_cell );
        
    const double rate12 = emission12_object.getIntegratedRate(); // rate in lab-frame
    
    const double probab12 = rate12 * dt;

    // int N12 = static_cast<int>( probab12 );
         
    // if( ran2() < (probab12 - static_cast<double>( N12 ) ) )
    //   N12++;

    // for( int i12 = 0; i12 < N12; i12++ )
    // {
    if ( probab12 > 1.0 )
    {
      cout << "P12=" << probab12 << ">1" << endl;
      again = true;
      cout << "dt (old) = " << dt << endl;
      dt = 0.5 / ( probab12 / dt );
      cout << "dt (new) = " << dt << endl;
      return;
    }
          
    if ( ran2() < probab12 )
    {
      // double pt_jscat = addedParticles[jscat].Mom.Perp();
      // double pt_nmb;

      // int jetEventIndex = -1;
      // if( pt_jscat > theAnalysis->getJetTracking_PT() )
      // {
      //   jetEventIndex = theAnalysis->addJetEvent_in( iscat, -1, jscat, c2to3, cs23, _cell.index, lambda_scaled / sqrt( s ) );
      // }

      int newIndex = emission12_offlineParticles_utility( emission12_object, jscat, again, nexttime );
    }
  }
}



void offlineHeavyIonCollision::scatt22_amongAddedParticles( cellContainer& _cellAdded, std::vector< int >& _allParticlesListAdded, const double scaleFactor, bool& again, const double nexttime )
{
  int iscat, jscat, typ;
  double s, cs22, Vrel;
  double M1, M2;
  double probab22;
  double md2g_wo_as, md2q_wo_as;
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

    for ( unsigned int j = i + 1; j < _allParticlesListAdded.size(); j++ )
    {
      jscat = _allParticlesListAdded[j];

      if ( !( addedParticles[jscat].FLAVOR == charm || addedParticles[jscat].FLAVOR == anti_charm) || addedParticles[jscat].dead || addedParticles[iscat].dead  ) // currently scattering among added particles is only used for c+cbar -> Jpsi + g. The last statement to check again whether iscat is dead is necessary since iscat could be deleted in the scattering before but the scatterings with all remaining jscat would be carried out.
      {
        continue; // jump to next particle from _cell.particleList
      }
      
      F1 = addedParticles[iscat].FLAVOR;
      M1 = addedParticles[iscat].m;
      F2 = addedParticles[jscat].FLAVOR;
      M2 = addedParticles[jscat].m;

      s = (addedParticles[iscat].Mom + addedParticles[jscat].Mom).M2();

      if ( s < 1.1*lambda2 )
      {
        probab22 = -1.0;
        cs22 = 0.0; //1/GeV^2
      }
      else
      {
        Vrel = VelRel( addedParticles[iscat].Mom, addedParticles[jscat].Mom, M1, M2 );

        //factor as (alpha_s) is not included in definitions of partcl[jscat].md2g, it will be multiplied in the scattering routines
        md2g_wo_as = ( addedParticles[iscat].md2g + addedParticles[jscat].md2g ) / 2.0;
        md2q_wo_as = ( addedParticles[iscat].md2q + addedParticles[jscat].md2q ) / 2.0;
        
        // determine largest temp of both charm quarks
        if( addedParticles[iscat].temperature > addedParticles[jscat].temperature )
          temperature = addedParticles[iscat].temperature;
        else
          temperature = addedParticles[jscat].temperature;

        scatt22_object.setParameter( addedParticles[iscat].Mom, addedParticles[jscat].Mom,
                                     F1, F2, M1, M2, s, md2g_wo_as , md2q_wo_as,
                                     theConfig->getKggQQb(), theConfig->getKgQgQ(), theConfig->getKappa_gQgQ(), theConfig->isConstantCrossSecGQ(),
                                     theConfig->getConstantCrossSecValueGQ(), theConfig->isIsotropicCrossSecGQ(), theConfig->getKfactor_light(),
                                     temperature, theConfig->getTdJpsi(), theConfig->isConstantCrossSecJpsi(), theConfig->getConstantCrossSecValueJpsi() ); // md2g, md2q are debye masses without the factor alpha_s which is multiplied in scattering22.cpp
      
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

  double averagedRate;
  double s, as, csgg, I32, probab32, xt;
  double md2g, md2q;
  double lambda_scaled;
  double ran2out;
  double betaDistEntry;
  int absorbedGluon = -1;
  double rate_added_sum = 0.0;
  
  //n32=0;

  const int nTotal = _cell.particleList.size();
  int nTotalAdded = _allParticlesListAdded.size();
  const int nGluons = _gluonList.size();
  int nGluonsAdded = _gluonListAdded.size();
  const int nAllQuarks = nTotal - nGluons;
  int nAllQuarksAdded = nTotalAdded - nGluonsAdded;
  
  for ( unsigned int i = 0; i < _allParticlesListAdded.size(); i++ )
  {
    double pt_addedParticle = addedParticles[_allParticlesListAdded[i]].Mom.Perp();
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
  
  double lambda = -1; // fm
  double pt_addedParticle = 0;

  double scaleForSelectedTriplets = 1.0;
  if ( allTriplets > 20 )
  {
    if ( nTotalAdded * 3 >= 20 )
    {
      consideredTriplets = nTotalAdded * 3;  // number of triplets to consider for 3 -> 2 processes
    }
    else
    {
      consideredTriplets = 20;  //20 is an empirical value
    }

    scaleForSelectedTriplets = static_cast<double>( allTriplets ) / static_cast<double>( consideredTriplets );

    for ( int i = 0; i < consideredTriplets && ( _gluonList.size() + _gluonListAdded.size() ) > 0 && _allParticlesListAdded.size() > 0 ; i++ )
    {
      rate_added_sum = 0.0;
      
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
      
      pt_addedParticle = addedParticles[kscat].Mom.Perp();
      if ( pt_addedParticle < theConfig->getMinimumPT() )
      {
        continue; //go to next particle triplet 
      }

      // compute actual mean free path
      xt = addedParticles[kscat].Pos.Perp();
      ringIndex = rings.getIndex( xt );
      
      switch(theConfig->getJetMfpComputationType())
      {
        case 0: // computeMfpLastTimestep
          
          if( addedParticles[kscat].lambda_added > 0.0 && addedParticles[kscat].lambda_added_old > 0.0 )
            lambda = ( addedParticles[kscat].lambda_added + addedParticles[kscat].lambda_added_old ) / 2.0 ; // mean free path of addedParticles in fm
          else if( addedParticles[kscat].lambda_added > 0.0 )
            lambda = addedParticles[kscat].lambda_added; // mean free path of addedParticles in fm
          else
          {
            cout << "error in scattOfflinePartclWithAddedPartcl: lambda negative: lambda = " << addedParticles[kscat].lambda_added << "  " << addedParticles[kscat].lambda_added_old << endl;
            cout << addedParticles[kscat].rate_added << "  " << iterate_mfp_bisection( _allParticlesList, _gluonList, kscat, dt, dv, addedParticles[kscat].lambda_added_old ) << endl;
            std::string errMsg = "Error in scattOfflinePartclWithAddedPartcl: lambda negative.";
            throw eHIC_error( errMsg );
          }
          break;
        case 1: // computeMfpIteration
          lambda = iterate_mfp_bisection( _allParticlesList, _gluonList, kscat, dt, dv, addedParticles[kscat].lambda_added_old ); // fm
//             lambda = iterateMFP( _allParticlesList, _gluonList, kscat, dt, dv ); //fm
          addedParticles[kscat].lambda_added = lambda;
          break;
        case 2: // computeMfpInterpolation
          lambda = theMFP.getMeanFreePath( addedParticles[kscat].Mom.E(), addedParticles[kscat].FLAVOR, rings[ringIndex].getEffectiveTemperature(), rings[ringIndex].getGluonDensity(), rings[ringIndex].getQuarkDensity(), fm );
          break;
        case 3: // fixedMfp
          lambda = theConfig->getFixedMfpAdded();
          break;
        case 4: // thermalMfpGluon
          lambda = theMFP.getMeanFreePath( 3.0 * rings[ringIndex].getEffectiveTemperature(), gluon, rings[ringIndex].getEffectiveTemperature(), rings[ringIndex].getGluonDensity(), rings[ringIndex].getQuarkDensity(), fm );
          break;
        default: 
          std::string errMsg = "Error in scattOfflinePartclWithAddedPartcl: wrong mfp determination type.";
          throw eHIC_error( errMsg );
      }

      s = ( particles_atTimeNow[iscat].Mom + particles_atTimeNow[jscat].Mom + addedParticles[kscat].Mom ).M2();

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
        
        xt = ( particles_atTimeNow[iscat].Pos.Perp() +
               particles_atTimeNow[jscat].Pos.Perp() +
               addedParticles[kscat].Pos.Perp() ) / 3.0;
        ringIndex = rings.getIndex( xt );

        if ( lambda > 0 )
        {
          lambda_scaled = lambda * sqrt( s ) / 0.197; // lambda in fm, sqrt(s) in GeV, lambda_scaled dimensionless
        }
        else
        {
          if ( averagedRate > epsilon )
          {
            lambda_scaled = sqrt( s ) / averagedRate;  //dimensionless
          }
          else
          {
            xsection_gg_gg csObj( s, md2g, md2q, &theI22, theConfig->getKfactor_light() );
            csgg = csObj.totalCrossSection();
            lambda_scaled = ( dv * rings[ringIndex].getGamma() * testpartcl * sqrt( s ) ) / ( nTotal * csgg * pow( 0.197, 3.0 ) );
          }
        }

        // create scattering32 object for the given 3 particles

        betaDistEntry = scatt32_object.setParameter( rings[ringIndex].getAveraged_v(),
                                                     particles_atTimeNow[iscat].Mom,
                                                     particles_atTimeNow[jscat].Mom,
                                                     addedParticles[kscat].Mom,
                                                     F1, F2, F3, sqrt( s ), md2g / s, lambda_scaled, as, theConfig->isMd2CounterTermInI23(), theConfig->isMatrixElement23_22qt(), theConfig->get23FudgeFactorLpm(), _gluonList.size() );  // create scattering32 object for the given 3 particles
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
      
      if ( probab32 > 0.0 )
      {
        rate_added_sum += scaleForSelectedTriplets * probab32 / dt * rings[ringIndex].getGamma();  // fm^-1
      }


      ran2out = ran2();
      if ( ran2out < probab32 )
      {
        double pt_iscat, pt_jscat;
        double pt_kscat = addedParticles[kscat].Mom.Perp();

        int jetEventIndex = -1;
        if( pt_kscat > theAnalysis->getJetTracking_PT() )
        {
          jetEventIndex = theAnalysis->addJetEvent_in( iscat, jscat, kscat, c3to2, I32, _cell.index, lambda_scaled / sqrt( s ) );
        }

        absorbedGluon = scatt32_offlineWithAddedParticles_utility( scatt32_object, _cellAdded.particleList, _allParticlesListAdded, _gluonListAdded, iscat, jscat, kscat, n32, nexttime );

        if ( absorbedGluon == 1 )
        {
          pt_jscat = particles_atTimeNow[jscat].Mom.Perp();
          pt_kscat = addedParticles[kscat].Mom.Perp();
          if ( jetEventIndex != -1 || pt_jscat > theAnalysis->getJetTracking_PT() || pt_kscat > theAnalysis->getJetTracking_PT() )
          {
            theAnalysis->addJetEvent_out( jetEventIndex, kscat, jscat, -1, c3to2 );
          }
        }
        else if ( absorbedGluon == 2 )
        {
          pt_iscat = particles_atTimeNow[iscat].Mom.Perp();
          pt_kscat = addedParticles[kscat].Mom.Perp();
          if ( jetEventIndex != -1 || pt_iscat > theAnalysis->getJetTracking_PT() || pt_kscat > theAnalysis->getJetTracking_PT() )
          {
            theAnalysis->addJetEvent_out( jetEventIndex, kscat, iscat, -1, c3to2 );
          }
        }
        else if ( absorbedGluon == 3 )
        {
          pt_iscat = particles_atTimeNow[iscat].Mom.Perp();
          pt_jscat = particles_atTimeNow[jscat].Mom.Perp();
          if ( jetEventIndex != -1 || pt_iscat > theAnalysis->getJetTracking_PT() || pt_jscat > theAnalysis->getJetTracking_PT() )
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
    
      addedParticles[kscat].rate_added_32 += rate_added_sum;
    }
  }
  else
  {
    for ( int m3 = 0; m3 < _allParticlesListAdded.size(); m3++ )
    {
      rate_added_sum = 0.0;
      
      kscat = _allParticlesListAdded[m3];
      pt_addedParticle = addedParticles[kscat].Mom.Perp();
      if ( pt_addedParticle < theConfig->getMinimumPT() )
      {
        continue; //go to next particle triplet 
      }
      
      // compute actual mean free path
      xt = addedParticles[kscat].Pos.Perp();
      ringIndex = rings.getIndex( xt );
      
      switch(theConfig->getJetMfpComputationType())
      {
        case 0: // computeMfpLastTimestep
          
          if( addedParticles[kscat].lambda_added > 0.0 && addedParticles[kscat].lambda_added_old > 0.0 )
            lambda = ( addedParticles[kscat].lambda_added + addedParticles[kscat].lambda_added_old ) / 2.0 ; // mean free path of addedParticles in fm
          else if( addedParticles[kscat].lambda_added > 0.0 )
            lambda = addedParticles[kscat].lambda_added; // mean free path of addedParticles in fm
          else
          {
            cout << "error in scattOfflinePartclWithAddedPartcl: lambda negative: lambda = " << addedParticles[kscat].lambda_added << "  " << addedParticles[kscat].lambda_added_old << endl;
            cout << addedParticles[kscat].rate_added << "  " << iterate_mfp_bisection( _allParticlesList, _gluonList, kscat, dt, dv, addedParticles[kscat].lambda_added_old ) << endl;
            std::string errMsg = "Error in scattOfflinePartclWithAddedPartcl: lambda negative.";
            throw eHIC_error( errMsg );
          }
          break;
        case 1: // computeMfpIteration
          lambda = iterate_mfp_bisection( _allParticlesList, _gluonList, kscat, dt, dv, addedParticles[kscat].lambda_added_old ); // fm
//             lambda = iterateMFP( _allParticlesList, _gluonList, kscat, dt, dv ); //fm
          addedParticles[kscat].lambda_added = lambda;
          break;
        case 2: // computeMfpInterpolation
          lambda = theMFP.getMeanFreePath( addedParticles[kscat].Mom.E(), addedParticles[kscat].FLAVOR, rings[ringIndex].getEffectiveTemperature(), rings[ringIndex].getGluonDensity(), rings[ringIndex].getQuarkDensity(), fm );
          break;
        case 3: // fixedMfp
          lambda = theConfig->getFixedMfpAdded();
          break;
        case 4: // thermalMfpGluon
          lambda = theMFP.getMeanFreePath( 3.0 * rings[ringIndex].getEffectiveTemperature(), gluon, rings[ringIndex].getEffectiveTemperature(), rings[ringIndex].getGluonDensity(), rings[ringIndex].getQuarkDensity(), fm );
          break;
        default: 
          std::string errMsg = "Error in scattOfflinePartclWithAddedPartcl: wrong mfp determination type.";
          throw eHIC_error( errMsg );
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

          s = (particles_atTimeNow[iscat].Mom + particles_atTimeNow[jscat].Mom + addedParticles[kscat].Mom).M2();


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
            

            xt = ( particles_atTimeNow[iscat].Pos.Perp() +
                   particles_atTimeNow[jscat].Pos.Perp() +
                   addedParticles[kscat].Pos.Perp() ) / 3.0;
            ringIndex = rings.getIndex( xt );

            if ( lambda > 0 )
            {
              lambda_scaled = lambda * sqrt( s ) / 0.197; // lambda in fm, sqrt(s) in GeV, lambda_scaled dimensionless
            }
            else
            {
              if ( averagedRate > epsilon )
              {
                lambda_scaled = sqrt( s ) / averagedRate;  //dimensionless
              }
              else
              {
                xsection_gg_gg csObj( s, md2g, md2q, &theI22, theConfig->getKfactor_light() );
                csgg = csObj.totalCrossSection();
                lambda_scaled = ( dv * rings[ringIndex].getGamma() * testpartcl * sqrt( s ) ) / ( nTotal * csgg * pow( 0.197, 3.0 ) );
              }
            }

            // create scattering32 object for the given 3 particles
            betaDistEntry = scatt32_object.setParameter( rings[ringIndex].getAveraged_v(),
                                                         particles_atTimeNow[iscat].Mom,
                                                         particles_atTimeNow[jscat].Mom,
                                                         addedParticles[kscat].Mom,
                                                         F1, F2, F3, sqrt( s ), md2g / s, lambda_scaled, as, theConfig->isMd2CounterTermInI23(), theConfig->isMatrixElement23_22qt(), theConfig->get23FudgeFactorLpm(), _gluonList.size() );  // create scattering32 object for the given 3 particles
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

          if ( probab32 > 0.0 )
          {
            rate_added_sum += scaleForSelectedTriplets * probab32 / dt * rings[ringIndex].getGamma();  // fm^-1
          }
          
          ran2out = ran2();
          if ( ran2out < probab32 )
          {
            double pt_iscat, pt_jscat;
            double pt_kscat = addedParticles[kscat].Mom.Perp();

            int jetEventIndex = -1;
            if( pt_kscat > theAnalysis->getJetTracking_PT() )
            {
              jetEventIndex = theAnalysis->addJetEvent_in( iscat, jscat, kscat, c3to2, I32, _cell.index, lambda_scaled / sqrt( s ) );
            }

            absorbedGluon = scatt32_offlineWithAddedParticles_utility( scatt32_object, _cellAdded.particleList, _allParticlesListAdded, _gluonListAdded, iscat, jscat, kscat, n32, nexttime );

            if ( absorbedGluon == 1 )
            {
              pt_jscat = particles_atTimeNow[jscat].Mom.Perp();
              pt_kscat = addedParticles[kscat].Mom.Perp();
              if ( jetEventIndex != -1 || pt_jscat > theAnalysis->getJetTracking_PT() || pt_kscat > theAnalysis->getJetTracking_PT() )
              {
                theAnalysis->addJetEvent_out( jetEventIndex, kscat, jscat, -1, c3to2 );
              }
            }
            else if ( absorbedGluon == 2 )
            {
              pt_iscat = particles_atTimeNow[iscat].Mom.Perp();
              pt_kscat = addedParticles[kscat].Mom.Perp();
              if ( jetEventIndex != -1 || pt_iscat > theAnalysis->getJetTracking_PT() || pt_kscat > theAnalysis->getJetTracking_PT() )
              {
                theAnalysis->addJetEvent_out( jetEventIndex, kscat, iscat, -1, c3to2 );
              }
            }
            else if ( absorbedGluon == 3 )
            {
              pt_iscat = particles_atTimeNow[iscat].Mom.Perp();
              pt_jscat = particles_atTimeNow[jscat].Mom.Perp();
              if ( jetEventIndex != -1 || pt_iscat > theAnalysis->getJetTracking_PT() || pt_jscat > theAnalysis->getJetTracking_PT() )
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
    
      addedParticles[kscat].rate_added_32 = rate_added_sum;
    }
  }
}





int offlineHeavyIonCollision::scatt23_offlineWithAddedParticles_utility( scattering23& scatt23_obj, cellContainer& _cell, int iscat, const int jscat, bool& again, const double nexttime )
{
  FLAVOR_TYPE F1, F2;
  double Tmax, TT, pt1, pt3, y, phi, pz1;
  int typ;

  int centralEtaIndex = static_cast<int>( etaBins.size() ) / 2;
  int etaIndex = centralEtaIndex;
  double eta = particles_atTimeNow[iscat].Pos.Rapidity();

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

  if( theConfig->doOutput_scatteredMediumParticles() && !particles_atTimeNow[iscat].isAlreadyInAddedParticles[addedParticles[jscat].N_EVENT_pp] )
  {
    scatteredMediumParticles.push_back( particles_atTimeNow[iscat] );
    scatteredMediumParticles.back().N_EVENT_pp = addedParticles[jscat].N_EVENT_pp;
  }

  F1 = particles_atTimeNow[iscat].FLAVOR;
  F2 = addedParticles[jscat].FLAVOR;

  ncoll++;
  ncoll23++;
  addedParticles[jscat].coll_id = ncoll;

  Tmax = std::max( particles_atTimeNow[iscat].Pos.T(), addedParticles[jscat].Pos.T());
  TT = ( nexttime - Tmax ) * ran2() + Tmax;

  addedParticles[jscat].Propagate( TT, addedParticles[jscat].X_traveled );

  // at this point scattering take place, therefore set properties of
  // last scattering point
  addedParticles[jscat].lastInt = addedParticles[jscat].Pos;

  nGet23Errors += scatt23_obj.getMomenta23( pt1, pt3, y, phi, pz1, typ, F1, F2 );
  VectorEPxPyPz P1new, P2new, P3new;
  scatt23_obj.setNewMomenta23( P1new, P2new, P3new, 
                               particles_atTimeNow[iscat].Pos, addedParticles[jscat].Pos,
                               pt1, pt3, y, phi, pz1 );

  double pt_out1 = P1new.Perp();
  double pt_out2 = P2new.Perp();

  //<<---------------------------------------------
  // set new properties for added particle
  // consider outgoing particle with highest pt if it not a tagged jet (charm, bottom, jpsi, etc)
  if ( pt_out1 > pt_out2 && !theConfig->isJetTagged() )
  {
    addedParticles[jscat].FLAVOR = F1;
    addedParticles[jscat].Mom = P1new;
    ParticleOffline tempParticle = particles_atTimeNow[iscat];
    tempParticle.FLAVOR = F2;
    tempParticle.Mom = P2new;
    tempParticle.unique_id = particles_atTimeNow[iscat].unique_id; // necessary ???
    tempParticle.N_EVENT_pp = addedParticles[jscat].N_EVENT_pp;
    tempParticle.N_EVENT_AA = addedParticles[jscat].N_EVENT_AA;
    if( theConfig->isScatt_furtherOfflineParticles() && !particles_atTimeNow[iscat].isAlreadyInAddedParticles[addedParticles[jscat].N_EVENT_pp] )
    {
      addedParticles.push_back( tempParticle );
      particles_atTimeNow[iscat].isAlreadyInAddedParticles[addedParticles[jscat].N_EVENT_pp] = true;
    }
  }
  else
  {
    addedParticles[jscat].FLAVOR = F2;
    addedParticles[jscat].Mom = P2new;
    ParticleOffline tempParticle = particles_atTimeNow[iscat];
    tempParticle.FLAVOR = F1;
    tempParticle.Mom = P1new;
    tempParticle.unique_id = particles_atTimeNow[iscat].unique_id; // necessary ???
    tempParticle.N_EVENT_pp = addedParticles[jscat].N_EVENT_pp;
    tempParticle.N_EVENT_AA = addedParticles[jscat].N_EVENT_AA;
    if( theConfig->isScatt_furtherOfflineParticles() && !particles_atTimeNow[iscat].isAlreadyInAddedParticles[addedParticles[jscat].N_EVENT_pp] )
    {
      addedParticles.push_back( tempParticle );
      particles_atTimeNow[iscat].isAlreadyInAddedParticles[addedParticles[jscat].N_EVENT_pp] = true;
    }
  }

  int newIndex = -1;
  if( theConfig->getNlightFlavorsAdded() >= 0 )
  {
    ParticleOffline tempParticle;
    tempParticle.FLAVOR = gluon;
    tempParticle.Mom = P3new;
    tempParticle.Pos = VectorTXYZ(TT,leftX + dx * ran2(),leftY + dy * ran2(),leftZ + deltaZ * ran2());
    tempParticle.PosInit = tempParticle.lastInt = tempParticle.Pos;
    tempParticle.MomInit = tempParticle.Mom;
    tempParticle.MomInit.E() = -tempParticle.MomInit.E(); // negative to indicate creation via 2->3 process

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
    tempParticle.T_creation = tempParticle.Pos.T();
    tempParticle.unique_id = ParticleOffline::unique_id_counter_added;
    --ParticleOffline::unique_id_counter_added;
    tempParticle.N_EVENT_pp = addedParticles[jscat].N_EVENT_pp;
    tempParticle.N_EVENT_AA = addedParticles[jscat].N_EVENT_AA;

  //   if ( sqrt( pow( tempParticle.PX, 2) + pow( tempParticle.PY, 2) ) > 3.0 )
  //   {
      addedParticles.push_back( tempParticle );
      newIndex = addedParticles.size() - 1;

      addedParticles[newIndex].Propagate( nexttime, addedParticles[newIndex].X_traveled );
//   }
  }
  
  return newIndex;
}



void offlineHeavyIonCollision::scatt22_offlineWithAddedParticles_utility( scattering22& scatt22_obj, std::list< int >& _cellMembersAdded, std::vector< int >& _allParticlesListAdded, const int iscat, const int jscat, int& typ, const double nexttime )
{
  FLAVOR_TYPE F1, F2;
  double Tmax, TT;
  double M1, M2;
  double t_hat;
  bool identical_particle_position_flipped = false; // is only important if jet is tagged and if both particles are identical
  bool qqbar_position_flipped = false; // is only important if jet is tagged and if outgoing particle is quark and anti-quark. Then use not always the second particle, which is the anti-quark, but randomize.

  if( theConfig->doOutput_scatteredMediumParticles() && !particles_atTimeNow[iscat].isAlreadyInAddedParticles[addedParticles[jscat].N_EVENT_pp] )
  {
    scatteredMediumParticles.push_back( particles_atTimeNow[iscat] );
    scatteredMediumParticles.back().N_EVENT_pp = addedParticles[jscat].N_EVENT_pp;
  }

  F1 = particles_atTimeNow[iscat].FLAVOR;
  M1 = particles_atTimeNow[iscat].m;

  F2 = addedParticles[jscat].FLAVOR;
  M2 = addedParticles[jscat].m;


  ncoll++;
  ncoll22++;
  addedParticles[jscat].coll_id = ncoll;

  Tmax = std::max( particles_atTimeNow[iscat].Pos.T(), addedParticles[jscat].Pos.T() );
  TT = ( nexttime - Tmax ) * ran2() + Tmax;
  addedParticles[jscat].Propagate( TT, addedParticles[jscat].X_traveled );

  // at this point scattering take place, therefore set properties of
  // last scattering point
  addedParticles[jscat].lastInt = addedParticles[jscat].Pos;

  if ( theConfig->isIsotropicCrossSection() == false )
  {
    // determine type of scattering, momentum transfer t_hat, new flavor and new masses for pQCD cross-sections
    scatt22_obj.getMomentaAndMasses22( F1, F2, M1, M2, t_hat, typ );
  }
  else
  {
    // determine type of scattering, momentum transfer t_hat, new flavor and new masses for isotropic processes
    scatt22_obj.getMomentaAndMasses22_isotropic( F1, F2, M1, M2, t_hat );
  }
  
  // translate momemtum transfer t_hat into actual momenta of outgoing
  // particles
  VectorEPxPyPz P1new, P2new;
  scatt22_obj.setNewMomenta22( P1new, P2new, 
                               particles_atTimeNow[iscat].Pos, particles_atTimeNow[jscat].Pos,
                               t_hat );
  
  // if tagged jet and identical particles, and large momentum transfer, the u channel is active which flips both particles. Since this is only an effect for identical particles and we cannot distinguish them anyhow, we take the particle which is going in the same direction as the incoming added particle. So for large t_hat we have to choose the other outgoing particle.
  if( theConfig->isJetTagged() && F1 == F2 )
  {
    double s = ( P1new + P2new ).M2();
    if( fabs( t_hat ) > s / 2.0 )
      identical_particle_position_flipped = true;
  }
  
  // if tagged jet and the outgoing particles are newly produced quark and anti-quarks (from gg->qqbar or qqbar->q'qbar' with typ 222 and 225, respectively) do not consider always the second particle, which is the anti-quark, but randomize
  if( theConfig->isJetTagged() && ( typ == 222 || typ == 225 ) )
  {
    if( ran2() < 0.5 )
      qqbar_position_flipped = true;
  }

  double pt_out1 = P1new.Perp();
  double pt_out2 = P2new.Perp();

  //<<---------------------------------------------
  // set new properties for added particle
  // consider outgoing particle with highest pt if it not a tagged jet (charm, bottom, jpsi, etc)
  if ( ( pt_out1 > pt_out2 && !theConfig->isJetTagged() ) || ( theConfig->isJetTagged() && ( identical_particle_position_flipped || qqbar_position_flipped ) ) )
  {
    addedParticles[jscat].FLAVOR = F1;
    addedParticles[jscat].m = M1;
    addedParticles[jscat].Mom = P1new;
    ParticleOffline tempParticle = particles_atTimeNow[iscat];
    tempParticle.FLAVOR = F2;
    tempParticle.Mom = P2new;
    tempParticle.unique_id = particles_atTimeNow[iscat].unique_id; // necessary ?
    tempParticle.N_EVENT_pp = addedParticles[jscat].N_EVENT_pp;
    tempParticle.N_EVENT_AA = addedParticles[jscat].N_EVENT_AA;
    if( theConfig->isScatt_furtherOfflineParticles() && !particles_atTimeNow[iscat].isAlreadyInAddedParticles[addedParticles[jscat].N_EVENT_pp] )
    {
      addedParticles.push_back( tempParticle );
      particles_atTimeNow[iscat].isAlreadyInAddedParticles[addedParticles[jscat].N_EVENT_pp] = true;
    }
  }
  else
  {
    addedParticles[jscat].FLAVOR = F2;
    addedParticles[jscat].m = M2;
    addedParticles[jscat].Mom = P2new;
    ParticleOffline tempParticle = particles_atTimeNow[iscat];
    tempParticle.FLAVOR = F1;
    tempParticle.Mom = P1new;
    tempParticle.unique_id = particles_atTimeNow[iscat].unique_id; // necessary ?
    tempParticle.N_EVENT_pp = addedParticles[jscat].N_EVENT_pp;
    tempParticle.N_EVENT_AA = addedParticles[jscat].N_EVENT_AA;
    if( theConfig->isScatt_furtherOfflineParticles() && !particles_atTimeNow[iscat].isAlreadyInAddedParticles[addedParticles[jscat].N_EVENT_pp] )
    {
      addedParticles.push_back( tempParticle );
      particles_atTimeNow[iscat].isAlreadyInAddedParticles[addedParticles[jscat].N_EVENT_pp] = true;
    }
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
      tempParticle.Mom = P1new;

      
      tempParticle.Pos = addedParticles[jscat].Pos + VectorTXYZ(0,dx * ( ran2() - 0.5 ),dy * ( ran2() - 0.5 ),dx * ( ran2() - 0.5 ));
      tempParticle.PosInit = tempParticle.lastInt = tempParticle.Pos;
      tempParticle.MomInit = tempParticle.Mom;
      tempParticle.MomInit.E() = -tempParticle.MomInit.E(); // negative to indicate creation via 2->3 process
      tempParticle.md2g = addedParticles[jscat].md2g;
      tempParticle.md2q = addedParticles[jscat].md2q;
      tempParticle.temperature = addedParticles[jscat].temperature;
      tempParticle.edge = false;
      tempParticle.dead = false;
      tempParticle.free = false;
      tempParticle.init = false;
      tempParticle.initially_produced = false;
      tempParticle.T_creation = tempParticle.Pos.T();
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

      addedParticles[newIndex].Propagate( nexttime, addedParticles[newIndex].X_traveled );
      
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


int offlineHeavyIonCollision::emission12_offlineParticles_utility( emission12& emission12_object, const int jscat, bool& again, const double nexttime )
{
  nemiss12++;

  FLAVOR_TYPE F1, F2;
  int newIndex = -1;
  
  F1 = addedParticles[jscat].FLAVOR;
  F2 = allFlavors;

  const double t = addedParticles[jscat].Pos.T();
  const double tt = ( nexttime - t ) * ran2() + t;

  addedParticles[jscat].Propagate( tt, addedParticles[jscat].X_traveled );

  // at this point scattering take place, therefore set properties of
  // last scattering point
  addedParticles[jscat].lastInt = addedParticles[jscat].Pos;

  int process_type;
  const double sampled_k = emission12_object.getMomenta12( process_type );
        
  VectorEPxPyPz P1new, P2new;

  P2new = emission12_object.setNewMomenta12( process_type, sampled_k, P1new, F1, F2 );
 
  addedParticles[jscat].FLAVOR = F1;
  addedParticles[jscat].Mom = P1new;

  if( F2 != allFlavors )
  {
    ParticleOffline tempParticle;
    tempParticle.FLAVOR = F2;
    tempParticle.Mom = P2new;
    tempParticle.Pos = addedParticles[jscat].Pos;
    tempParticle.PosInit = tempParticle.lastInt = tempParticle.Pos;
    tempParticle.MomInit = tempParticle.Mom;
    tempParticle.MomInit.E() = -tempParticle.MomInit.E(); // negative to indicate creation via 2->3 process

    tempParticle.m = 0.0;
    tempParticle.md2g = addedParticles[jscat].md2g;
    tempParticle.md2q = addedParticles[jscat].md2q;
    tempParticle.temperature = addedParticles[jscat].temperature;
    tempParticle.edge = false;
    tempParticle.dead = false;
    tempParticle.free = false;
    tempParticle.cell_id = -1;
    tempParticle.coll_id = ncoll;
    tempParticle.collisionTime = infinity;
    tempParticle.collisionPartner = -1;
    tempParticle.initially_produced = false;
    tempParticle.init = false;
    tempParticle.X_traveled = 0.0;
    tempParticle.T_creation = tempParticle.Pos.T();
    tempParticle.unique_id = ParticleOffline::unique_id_counter_added;
    --ParticleOffline::unique_id_counter_added;
    tempParticle.N_EVENT_pp = addedParticles[jscat].N_EVENT_pp;
    tempParticle.N_EVENT_AA = addedParticles[jscat].N_EVENT_AA;

    addedParticles.push_back( tempParticle );
    newIndex = addedParticles.size() - 1;

    addedParticles[newIndex].Propagate( nexttime, addedParticles[newIndex].X_traveled );          
  }

  return newIndex;
}

void offlineHeavyIonCollision::scatt22_amongAddedParticles_utility( scattering22& scatt22_obj, std::list< int >& _cellMembersAdded, std::vector< int >& _allParticlesListAdded, const int iscat, const int jscat, int& typ, const double nexttime )
{
  FLAVOR_TYPE F1, F2;
  double Tmax, TT;
  double M1, M2;
  double t_hat;
  
  ParticleOffline temp_particle_iscat = addedParticles[iscat];
  ParticleOffline temp_particle_jscat = addedParticles[jscat];

  F1 = addedParticles[iscat].FLAVOR;
  M1 = addedParticles[iscat].m;

  F2 = addedParticles[jscat].FLAVOR;
  M2 = addedParticles[jscat].m;

  Tmax = std::max( addedParticles[iscat].Pos.T(), addedParticles[jscat].Pos.T() );
  TT = ( nexttime - Tmax ) * ran2() + Tmax;
  addedParticles[iscat].Propagate( TT, addedParticles[iscat].X_traveled );
  addedParticles[jscat].Propagate( TT, addedParticles[jscat].X_traveled );

  // at this point scattering take place, therefore set properties of last scattering point
  addedParticles[iscat].lastInt = addedParticles[iscat].Pos;
  addedParticles[jscat].lastInt = addedParticles[jscat].Pos;

  // determine type of scattering, momentum transfer t_hat, new flavor and new masses
  scatt22_obj.getMomentaAndMasses22( F1, F2, M1, M2, t_hat, typ );
  // translate momemtum transfer t_hat into actual momenta of outgoing particles
  VectorEPxPyPz P1new, P2new;
  scatt22_obj.setNewMomenta22( P1new, P2new, 
                               addedParticles[iscat].Pos, addedParticles[jscat].Pos,
                               t_hat );
  
  // figure out if this is a real event, that is, the ccbar is annihilated (in every 1/N_test_Jpsi collision). The Jpsi is created in any case
  if( ran2() * theConfig->getJpsiTestparticles() < 1.0 )
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
      addedParticles[iscat].Mom = P2new;
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
      tempParticle.Mom = P2new;
      tempParticle.PosInit = tempParticle.lastInt = tempParticle.Pos = addedParticles[iscat].Pos;
      tempParticle.MomInit = tempParticle.Mom;
      tempParticle.MomInit.E() = -tempParticle.MomInit.E(); // negative to indicate creation via 2->3 process
      tempParticle.md2g = addedParticles[iscat].md2g;
      tempParticle.md2q = addedParticles[iscat].md2q;
      tempParticle.temperature = addedParticles[iscat].temperature;
      tempParticle.edge = false;
      tempParticle.dead = false;
      tempParticle.free = false;
      tempParticle.init = false;
      tempParticle.initially_produced = false;
      tempParticle.T_creation = tempParticle.Pos.T();
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

      addedParticles[newIndex].Propagate( nexttime, addedParticles[newIndex].X_traveled );
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
  double Tmax, TT, u, phi;
  FLAVOR_TYPE F1, F2, F3;
  int typ;

  if( theConfig->doOutput_scatteredMediumParticles() && !particles_atTimeNow[iscat].isAlreadyInAddedParticles[addedParticles[kscat].N_EVENT_pp] )
  {
    scatteredMediumParticles.push_back( particles_atTimeNow[iscat] );
    scatteredMediumParticles.back().N_EVENT_pp = addedParticles[kscat].N_EVENT_pp;
  }
  
  if( theConfig->doOutput_scatteredMediumParticles() && !particles_atTimeNow[jscat].isAlreadyInAddedParticles[addedParticles[kscat].N_EVENT_pp] )
  {
    scatteredMediumParticles.push_back( particles_atTimeNow[jscat] );
    scatteredMediumParticles.back().N_EVENT_pp = addedParticles[kscat].N_EVENT_pp;
  }

  F1 = particles_atTimeNow[iscat].FLAVOR;
  F2 = particles_atTimeNow[jscat].FLAVOR;
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

  Tmax = std::max( std::max( particles_atTimeNow[iscat].Pos.T(), 
                             particles_atTimeNow[jscat].Pos.T() ),
                   addedParticles[kscat].Pos.T() );
  TT = ( nexttime - Tmax ) * ran2() + Tmax;

  addedParticles[kscat].Propagate( TT, addedParticles[kscat].X_traveled );
  
  int absorbedGluon = scatt32_obj.getMomenta32( u, phi, typ, F1, F2 );

  VectorEPxPyPz P1new, P2new;
  scatt32_obj.setNewMomenta32( P1new, P2new, u, phi );
// Attention: Has to be modified after implementing heavy-quark 3->2 processes.
  P1new.E() = sqrt( P1new.vec2() );
  P2new.E() = sqrt( P2new.vec2() );

  double pt_out1 = P1new.Perp();
  double pt_out2 = P2new.Perp();

  if ( absorbedGluon == 3 )
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
    ParticleOffline tempParticle = particles_atTimeNow[iscat];
    tempParticle.FLAVOR = F1;
    tempParticle.Mom = P1new;
    tempParticle.unique_id = particles_atTimeNow[iscat].unique_id; // necessary ???
    tempParticle.N_EVENT_pp = addedParticles[kscat].N_EVENT_pp;
    tempParticle.N_EVENT_AA = addedParticles[kscat].N_EVENT_AA;
    if( theConfig->isScatt_furtherOfflineParticles() && !particles_atTimeNow[iscat].isAlreadyInAddedParticles[addedParticles[kscat].N_EVENT_pp] )
    {
      addedParticles.push_back( tempParticle );
      particles_atTimeNow[iscat].isAlreadyInAddedParticles[addedParticles[kscat].N_EVENT_pp] = true;
    }
    
    ParticleOffline tempParticle2 = particles_atTimeNow[jscat];
    tempParticle2.FLAVOR = F2;
    tempParticle2.Mom = P2new;
    tempParticle2.unique_id = particles_atTimeNow[jscat].unique_id; // necessary ???
    tempParticle2.N_EVENT_pp = addedParticles[kscat].N_EVENT_pp;
    tempParticle2.N_EVENT_AA = addedParticles[kscat].N_EVENT_AA;
    if( theConfig->isScatt_furtherOfflineParticles() && !particles_atTimeNow[jscat].isAlreadyInAddedParticles[addedParticles[kscat].N_EVENT_pp] )
    {
      addedParticles.push_back( tempParticle );
      particles_atTimeNow[jscat].isAlreadyInAddedParticles[addedParticles[kscat].N_EVENT_pp] = true;
    }
  }
  else if ( absorbedGluon == 2 ) 
  {
    if ( !theConfig->isJetTagged() && pt_out1 > pt_out2 )
    {

//         if jet particle changes to a gluon, it is added to gluon list in cell
      if ( (addedParticles[kscat].FLAVOR != gluon) && (F1 == gluon))
        _gluonListAdded.push_back( kscat );

      addedParticles[kscat].FLAVOR = F1;
      addedParticles[kscat].Mom = P1new;

      if( theConfig->isScatt_furtherOfflineParticles() )
      {
        string errMsg = "Recoiled particles should be evovlved further while not using tagged jets. Unrecoverable error!";
        throw eHIC_error( errMsg );
      }
    }
    else
    {
      if( ( addedParticles[kscat].FLAVOR != gluon ) && ( F2 == gluon ) )
        _gluonListAdded.push_back( kscat ); //  if jet particle flavor changes to a gluon, it is added to gluon list in cell

      addedParticles[kscat].FLAVOR = F2;
      addedParticles[kscat].Mom = P2new;
      ParticleOffline tempParticle = particles_atTimeNow[iscat];
      tempParticle.FLAVOR = F1;
      tempParticle.Mom = P1new;
      tempParticle.unique_id = particles_atTimeNow[iscat].unique_id; // necessary ???
      tempParticle.N_EVENT_pp = addedParticles[kscat].N_EVENT_pp;
      tempParticle.N_EVENT_AA = addedParticles[kscat].N_EVENT_AA;
      if( theConfig->isScatt_furtherOfflineParticles() && !particles_atTimeNow[iscat].isAlreadyInAddedParticles[addedParticles[kscat].N_EVENT_pp] )
      {
        addedParticles.push_back( tempParticle );
        particles_atTimeNow[iscat].isAlreadyInAddedParticles[addedParticles[kscat].N_EVENT_pp] = true;
      }
    }
  }
  else // ( absorbedGluon == 1 ) 
  {
    if ( !theConfig->isJetTagged() && pt_out1 > pt_out2 )
    {

      if ( (addedParticles[kscat].FLAVOR != gluon) && (F1 == gluon))
        _gluonListAdded.push_back( kscat ); //  if jet particle flavor changes to a gluon, it is added to gluon list in cell
        
      addedParticles[kscat].FLAVOR = F1;
      addedParticles[kscat].Mom = P1new;

      if( theConfig->isScatt_furtherOfflineParticles() )
      {
        string errMsg = "Recoiled particles should be evovlved further while not using tagged jets. Unrecoverable error!";
        throw eHIC_error( errMsg );
      }
    }
    else
    {
      if ( (addedParticles[kscat].FLAVOR != gluon) && (F2 == gluon))
        _gluonListAdded.push_back( kscat ); //  if jet particle flavor changes to a gluon, it is added to gluon list in cell
      
      addedParticles[kscat].FLAVOR = F2;
      addedParticles[kscat].Mom = P2new;
      ParticleOffline tempParticle = particles_atTimeNow[jscat];
      tempParticle.FLAVOR = F1;
      tempParticle.Mom = P1new;
      tempParticle.unique_id = particles_atTimeNow[jscat].unique_id; // necessary
      tempParticle.N_EVENT_pp = addedParticles[kscat].N_EVENT_pp;
      tempParticle.N_EVENT_AA = addedParticles[kscat].N_EVENT_AA;
      if( theConfig->isScatt_furtherOfflineParticles() && !particles_atTimeNow[jscat].isAlreadyInAddedParticles[addedParticles[kscat].N_EVENT_pp] )
      {
        addedParticles.push_back( tempParticle );
        particles_atTimeNow[jscat].isAlreadyInAddedParticles[addedParticles[kscat].N_EVENT_pp] = true;
      }
    }
  }

  return absorbedGluon;
}




void offlineHeavyIonCollision::jpsi_dissociation_td( const double time )
{
  VectorEPxPyPz P_jpsi, P_jpsi_cm, P1cm, P2cm, P1, P2;
  lorentz LL;
  double pp, phi, costheta, sintheta;
  
  for ( int j = 0;j < addedParticles.size();j++ )
  {
    if ( addedParticles[j].FLAVOR == jpsi && addedParticles[j].temperature >= theConfig->getTdJpsi() && FPT_COMP_LE( addedParticles[j].Pos.T(), time ) && !addedParticles[j].dead ) // Jpsi in cell with temperature higher than Td
    {
      ns_heavy_quarks::jpsi_dissociation++;
      ns_heavy_quarks::jpsi_dissociation_from_temperature++;
      
      // decay Jpsi to charm and anti-charm quarks

      // figure out if this is a real event, that is, a ccbar is created from the Jpsi (in every 1/N_test_Jpsi collision). The Jpsi is deleted in any case.
      if( ran2()*theConfig->getJpsiTestparticles() < 1.0 )
      {
        LL.setBetaCM( addedParticles[j].Mom );
        
        // boost to cms frame
        LL.boost( P_jpsi, P_jpsi_cm);
        
        // conserve energy
        P1cm.E() = P_jpsi_cm.E() / 2.0;
        // momentum
        if ( P1cm.E() < pow( Particle::getMass( charm ), 2.0 ) )
        {
          cout << "error in jpsi_dissociation_temperature(): " << P1cm.E() << endl;
          pp = 0.0;
          P1cm.E() = Particle::getMass( charm );
        }
        else
        {
          pp = sqrt( pow( P1cm.E(), 2.0 ) - pow( Particle::getMass( charm ) , 2.0 ) );
        }
        
        // sample momentum isotropic for first charm quark
        phi = 2.0 * M_PI * ran2();
        costheta = 2.0 * ran2() - 1.0;
        sintheta= sqrt( 1.0 - costheta*costheta );
        P1cm.Px() = pp * sintheta * cos(phi);
        P1cm.Py() = pp * sintheta * sin(phi);
        P1cm.Pz() = pp * costheta;
        
        // second charm quark
        P2cm = P1cm;
        P2cm.Minus3();
        
        // boost both charm quarks back
        LL.boostInv(P1cm,P2cm, P1,P2);
        
        // update momenta of particles
        addedParticles[j].FLAVOR = charm;
        addedParticles[j].m = Particle::getMass( charm );
        addedParticles[j].Mom = P2;
        
        // add second charm quark
        ParticleOffline tempParticle;
        tempParticle.FLAVOR = anti_charm;
        tempParticle.m = Particle::getMass( charm );
        tempParticle.Mom = P1;
        tempParticle.Pos = addedParticles[j].Pos + VectorTXYZ(0,dx * ( ran2() - 0.5 ),dy * ( ran2() - 0.5 ),dx * ( ran2() - 0.5 ));
        tempParticle.PosInit = tempParticle.lastInt = tempParticle.Pos;
        tempParticle.MomInit = tempParticle.Mom;
        tempParticle.MomInit.E() = -tempParticle.MomInit.E(); // negative to indicate creation via 2->3 process
        tempParticle.md2g = addedParticles[j].md2g;
        tempParticle.md2q = addedParticles[j].md2q;
        tempParticle.temperature = addedParticles[j].temperature;
        tempParticle.edge = false;
        tempParticle.dead = false;
        tempParticle.free = false;
        tempParticle.init = false;
        tempParticle.initially_produced = false;
        tempParticle.T_creation = tempParticle.Pos.T();
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
    if ( addedParticles[index].Pos.T() < nexttime )
    {
      addedParticles[index].Propagate( nexttime );
    }
  }
}



/**
* This routine removes particle stored in deadParticleList (global to this translation unit via unnamed namespace)
*/
void offlineHeavyIonCollision::removeDeadParticles( )
{
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



//provides iterative calculation of MFP for high-pt particles, returned lambda and start value in fm
double offlineHeavyIonCollision::iterateMFP( std::vector< int >& _allParticlesList, std::vector< int >& _gluonList, const int jetID, const double dt, const double dv )
{
  VectorEPxPyPz P1, P2, P3;
  int iscat, jscat;
  int n32 = 0, n22 = 0, n23 = 0;
  double lambda_scaled, s;
  double probab22 = 0, probab23 = 0, probab32 = 0;
  double cs22, cs23, I32;
  double R22, R23, R32;
  double Vrel, md2g_wo_as, md2q_wo_as;
  double lambda, lambdaAvr; // fm
  int iter = 0;
  deque<double> lambdaArray;
  double betaDistEntry;
  double M1, M2;
  double velocity_jet_restframe;
  
  int initialStateIndex = -1;
  FLAVOR_TYPE F1, F2, F3;
  
  const double epsilon = 0.01;
  const int nIterationsMax = 8;
  bool converged = false;
  
  
  //----------------------------- to which ring does the jet belong? ---------------------------
  int nc = rings.getIndex( addedParticles[jetID] );
  //--------------------------------------------------------------------------------------------
  
  F1 = addedParticles[jetID].FLAVOR;
  M1 = addedParticles[jetID].m;

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
    lambda = 1 / jetRate * 0.197; // fm
  }
  else
  {
    do  // ensure found cross section is nonzero
    {
      do  // pick random particle (that is neihter dead nor the jet) from the _particleList
      {
	int selected = static_cast<int>( ran2() * _allParticlesList.size() );
	if ( selected >= _allParticlesList.size() )
	{
	  continue;
	}      
	iscat = _allParticlesList[selected];
      }
      while ( particles_atTimeNow[iscat].dead );
      
      md2g_wo_as = ( addedParticles[jetID].md2g + particles_atTimeNow[iscat].md2g ) / 2.0;
      md2q_wo_as = ( addedParticles[jetID].md2q + particles_atTimeNow[iscat].md2q ) / 2.0;
      
      s = (addedParticles[jetID].Mom + particles_atTimeNow[iscat].Mom).M2();
      
      xsection_gg_gg csObj( s, md2g_wo_as, md2q_wo_as, &theI22, theConfig->getKfactor_light() );
      csgg = csObj.totalCrossSection();
    } while (csgg < 1.0e-7);

    lambda = ( dv * rings[nc].getGamma() * testpartcl   ) / ( pow( 0.197, 3.0 ) * _allParticlesList.size() * csgg ) * 0.197; // fm  

  }
  //--------------------------------------------------------------------------------------------
  
  scattering22 scatt22_object( &theI22 );
  scattering32 scatt32_object;
  scattering23 scatt23_object( &theI23_massless, &theI23_charm_m1, &theI23_charm_m2, &theI23_bottom_m1, &theI23_bottom_m2 );
  
  // the rate for 2->2 does not depend on lambda and therefore does not change during the iteration. So it is more efficient to calculate it here:
  //------------------------ 2<->2-----------------------
  if( theConfig->doScattering_22() )
  {
    for ( int m1 = 0; m1 < _allParticlesList.size(); m1++ )
    {
      jscat = _allParticlesList[m1];
      if ( !particles_atTimeNow[jscat].dead )
      {
        F2 = particles_atTimeNow[jscat].FLAVOR;
        M2 = particles_atTimeNow[jscat].m;

        s = (addedParticles[jetID].Mom + particles_atTimeNow[jscat].Mom).M2();
        
        if ( s > 1.1*lambda2 )
        {
          n22++;
          Vrel = VelRel( addedParticles[jetID].Mom, particles_atTimeNow[jscat].Mom, M1,M2 );
          
          md2g_wo_as = ( addedParticles[jetID].md2g + particles_atTimeNow[jscat].md2g ) / 2.0;
          md2q_wo_as = ( addedParticles[jetID].md2q + particles_atTimeNow[jscat].md2q ) / 2.0;
          
          scatt22_object.setParameter( addedParticles[jetID].Mom, particles_atTimeNow[jscat].Mom,
                                       F1, F2, M1, M2, s, md2g_wo_as , md2q_wo_as,
                                        theConfig->getKggQQb(), theConfig->getKgQgQ(), theConfig->getKappa_gQgQ(), 
                                        theConfig->isConstantCrossSecGQ(),
                                        theConfig->getConstantCrossSecValueGQ(), theConfig->isIsotropicCrossSecGQ(),
                                        theConfig->getKfactor_light() ); // md2g, md2q are debye masses without the factor alpha_s which is multiplied in scattering22.cpp
          cs22 = scatt22_object.getXSection22();
          probab22 += pow( 0.197, 2.0 ) * cs22 * Vrel * dt / ( dv * testpartcl );
        }
      }
    }
  }
  R22 = probab22 / dt * rings[nc].getGamma();
  
  
  do
  {    
    //------------------------ 2->3 -----------------------
    if( theConfig->doScattering_23() )
    {
      for ( int m1 = 0; m1 < _allParticlesList.size(); m1++ )
      {
        jscat = _allParticlesList[m1];
        if ( !particles_atTimeNow[jscat].dead )
        {
          F2 = particles_atTimeNow[jscat].FLAVOR;
          M2 = particles_atTimeNow[jscat].m;

          s = (addedParticles[jetID].Mom + particles_atTimeNow[jscat].Mom).M2();
          
          if ( s > 1.1*lambda2 )
          {
            Vrel = VelRel( addedParticles[jetID].Mom, particles_atTimeNow[jscat].Mom, M1,M2 );
            
            md2g_wo_as = ( addedParticles[jetID].md2g + particles_atTimeNow[jscat].md2g ) / 2.0;
            md2q_wo_as = ( addedParticles[jetID].md2q + particles_atTimeNow[jscat].md2q ) / 2.0;

            n23++;
            lambda_scaled = ( lambda / 0.197 * sqrt( s ) ); // lambda in fm, sqrt(s) in GeV, lambda_scaled dimensionless
            
            betaDistEntry = scatt23_object.setParameter( rings[nc].getAveraged_v(),
                                                         addedParticles[jetID].Mom,
                                                         particles_atTimeNow[jscat].Mom,
                                                         F1, F2, M1, M2, sqrt( s ), 
                        md2g_wo_as / s, lambda_scaled, 
                        theConfig->getK23LightPartons(), theConfig->getK23HeavyQuarks(),
                        theConfig->getKappa23LightPartons(), theConfig->getKappa23HeavyQuarks(),
                        theConfig->I23onlineIntegrationIsSet(),
                        theConfig->get23GluonFormationTimeTyp(), theConfig->getMatrixElement23(), theConfig->isMd2CounterTermInI23(), theConfig->get23FudgeFactorLpm(), _gluonList.size(), theConfig->isMatrixElement23_22qt() );    
            cs23 = scatt23_object.getXSection23( initialStateIndex ); //1/GeV^2

            probab23 += pow( 0.197, 2.0 ) * cs23 * Vrel * dt / ( dv * testpartcl );
            
          }
        }
      }
    }
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
          
          s = ( addedParticles[jetID].Mom + particles_atTimeNow[iscat].Mom + particles_atTimeNow[jscat].Mom ).M2();
          
          if ( s > 1.1*lambda2 )
          {
            n32++;
            lambda_scaled = ( lambda / 0.197 * sqrt( s ) ); // lambda in fm, sqrt(s) in GeV, lambda_scaled dimensionless
            
            md2g_wo_as = ( particles_atTimeNow[iscat].md2g + particles_atTimeNow[jscat].md2g + addedParticles[jetID].md2g ) / 3.0;
            md2q_wo_as = ( particles_atTimeNow[iscat].md2q + particles_atTimeNow[jscat].md2q + addedParticles[jetID].md2q ) / 3.0;
            
            // create scattering32 object for the given 3 particles
            betaDistEntry = scatt32_object.setParameter( rings[nc].getAveraged_v(),
                                                         addedParticles[jetID].Mom,
                                                         particles_atTimeNow[iscat].Mom,
                                                         particles_atTimeNow[jscat].Mom,
                                                         F1, F2, F3, sqrt( s ), md2g_wo_as * coupling::get_constant_coupling() / s, lambda_scaled, coupling::get_constant_coupling(), false, false, theConfig->get23FudgeFactorLpm(), _gluonList.size() );
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
                F3 = particles_atTimeNow[jscat].FLAVOR;

                // at least one of the particles must be a gluon
                // otherwise go to next step in the loop
                if ( !( F1 == gluon || F2 == gluon || F3 == gluon ) )
                {
                  continue;
                }
                s = ( addedParticles[jetID].Mom + particles_atTimeNow[iscat].Mom + particles_atTimeNow[jscat].Mom ).M2();
                
                if ( s > 1.1*lambda2 )
                {
                  n32++;
                  lambda_scaled = ( lambda / 0.197 * sqrt( s ) ); // lambda in fm, sqrt(s) in GeV, lambda_scaled dimensionless
                  
                  md2g_wo_as = ( particles_atTimeNow[iscat].md2g + particles_atTimeNow[jscat].md2g + addedParticles[jetID].md2g ) / 3.0;
                  md2g_wo_as = ( particles_atTimeNow[iscat].md2q + particles_atTimeNow[jscat].md2q + addedParticles[jetID].md2q ) / 3.0;
                  // HACK for N_f = 0 background
  //                 md2g = as * ( particles_atTimeNow[iscat].md2g + particles_atTimeNow[jscat].md2g + addedParticles[jetID].md2g ) / 3.0 * 2;
  //                 md2q = md2g * 2.0 / 9.0;
                  
                  // create scattering32 object for the given 3 particles
                  betaDistEntry = scatt32_object.setParameter( rings[nc].getAveraged_v(),
                                                               addedParticles[jetID].Mom,
                                                               particles_atTimeNow[iscat].Mom,
                                                               particles_atTimeNow[jscat].Mom,
                                                               F1, F2, F3, sqrt( s ), md2g_wo_as * coupling::get_constant_coupling() / s, lambda_scaled, coupling::get_constant_coupling(), false, false, theConfig->get23FudgeFactorLpm(), _gluonList.size() );
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
    
    // velocity of cell
    VectorXYZ v_cell = rings[nc].getAveraged_v();
              
    // velocity of added particle in lab frame
    VectorTXYZ v_jet = addedParticles[jetID].Mom * (1.0 / addedParticles[jetID].Mom.E());
//     double velocity_lab = sqrt( 1.0 - pow( addedParticles[jetID].m ,2.0) / pow( addedParticles[jetID].E ,2.0) );
    
    // Compute velocity of added particle in rest frame of fluid. The general expression for adding two velocities is not symmetric in v1 and v2. Therefore compute both cases and take average.
    double velocity_rest_1 = addVelocities( v_cell.X(), v_cell.Y(), v_cell.Z(), v_jet.X(), v_jet.Y(), v_jet.Z() );
    double velocity_rest_2 = addVelocities( v_jet.X(), v_jet.Y(), v_jet.Z(), v_cell.X(), v_cell.Y(), v_cell.Z() );
    velocity_jet_restframe = ( velocity_rest_1 + velocity_rest_2 ) / 2.0;

    const double epsilon_rate = 1E-3; // 0.001 fm^-1 which corresponds to a lambda = 1000 fm
    // get new lambda which one obtains with the previous employed lambda for the cross section
    // If the cross sections are very small or 0, R22 = R23 = R32 = 0 or very small which causes lambda to be very large. Therefore, epsilon serves as a cut-off.  However, since the cross section is small probably no scattering would take place anyhow. So it does not matter which value is returned...
    if( (R22 + R23 + R32) < epsilon_rate ) 
      lambda = velocity_jet_restframe / epsilon_rate; //fm
    else
      lambda = velocity_jet_restframe / (R22 + R23 + R32); //fm
    
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
    probab23 = probab32 = 0;
  }
  while ( !converged && iter < nIterationsMax );
  
  lambdaAvr = 0;
  for ( int m = 0; m < lambdaArray.size(); m++ )
  {
    lambdaAvr += lambdaArray[m];
  }
  lambdaAvr = lambdaAvr / lambdaArray.size();
  
  return lambdaAvr; // fm
}




//provides iterative calculation of MFP for added particles, returns lambda in fm
// lambda_old in fm
double offlineHeavyIonCollision::iterate_mfp_bisection( std::vector< int >& _allParticlesList, std::vector< int >& _gluonList, const int jetID, const double dt, const double dv, const double lambda_old )
{
  int iscat, jscat, kscat;
  int n32 = 0, n22 = 0, n23 = 0;
  double lambda_scaled, s;
  double probab22 = 0, probab23 = 0, probab32 = 0;
  double cs22, cs23, I32;
  double R22, R23, R32;
  double as, Vrel, md2g_wo_as, md2q_wo_as;
  double betaDistEntry;
  double M1, M2;
  double velocity_jet_restframe;
  
  int initialStateIndex = -1;
  FLAVOR_TYPE F1, F2, F3;

  double lambda,lambdaRes; // fm
  bool lambda_converged, converged, lambda_lambdaRes_converged, lambdaRes_converged;
  double lambdaAvr, lambdaResAvr; // fm
  deque<double> lambdaArray, lambdaResArray;
  
  const double epsilon = 0.01;
  const int nIterationsMax = 30;
  
  double lambdaMax, lambdaMin;
  const double lambda_range_variation = theConfig->getMfpAddedRangeVariation(); // in %
  const double offset = 2.0 * ran2() - 1.0;
  if( lambda_old > 0.0 )
  {
    lambdaMax = lambda_old * ( 1.0 + lambda_range_variation / 100.0 + lambda_range_variation / 1000.0 * offset );
    lambdaMin = lambda_old * ( 1.0 - lambda_range_variation / 100.0 + lambda_range_variation / 1000.0 * offset );
    if( lambdaMin < 0.0 )
      lambdaMin = 0.0;
    if( lambdaMax > 2.0 )
      lambdaMax = 2.0;
  }
  else
  {
    lambdaMax = 2.0;  //!!
    lambdaMin = 0.0;
  }
  lambda = (lambdaMax - lambdaMin) / 2.0 + lambdaMin; // fm

  int nIterations = 0;

  
  
  //----------------------------- to which ring does the jet belong? ---------------------------
  int nc = rings.getIndex( addedParticles[jetID] );
  //--------------------------------------------------------------------------------------------
  
  F1 = addedParticles[jetID].FLAVOR;
  M1 = addedParticles[jetID].m;
  
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

  scattering22 scatt22_object( &theI22 );
  scattering32 scatt32_object;
  scattering23 scatt23_object( &theI23_massless, &theI23_charm_m1, &theI23_charm_m2, &theI23_bottom_m1, &theI23_bottom_m2 );
  
  // the rate for 2->2 does not depend on lambda and therefore does not change during the iteration. So it is more efficient to calculate it here:
  //------------------------ 2<->2-----------------------
  if( theConfig->doScattering_22() )
  {
    for ( int m1 = 0; m1 < _allParticlesList.size(); m1++ )
    {
      jscat = _allParticlesList[m1];
      if ( !particles_atTimeNow[jscat].dead )
      {
        F2 = particles_atTimeNow[jscat].FLAVOR;
        M2 = particles_atTimeNow[jscat].m;

        s = ( addedParticles[jetID].Mom + particles_atTimeNow[jscat].Mom ).M2();
        
        if ( s > 1.1*lambda2 )
        {
          n22++;
          Vrel = VelRel( addedParticles[jetID].Mom, particles_atTimeNow[jscat].Mom, M1,M2 );
          
          md2g_wo_as = ( addedParticles[jetID].md2g + particles_atTimeNow[jscat].md2g ) / 2.0;
          md2q_wo_as = ( addedParticles[jetID].md2q + particles_atTimeNow[jscat].md2q ) / 2.0;
          
          scatt22_object.setParameter( addedParticles[jetID].Mom, particles_atTimeNow[jscat].Mom,
                                       F1, F2, M1, M2, s, md2g_wo_as , md2q_wo_as,
                                        theConfig->getKggQQb(), theConfig->getKgQgQ(), theConfig->getKappa_gQgQ(), 
                                        theConfig->isConstantCrossSecGQ(),
                                        theConfig->getConstantCrossSecValueGQ(), theConfig->isIsotropicCrossSecGQ() ); // md2g_wo_as, md2q_wo_as are debye masses without the factor alpha_s which is multiplied in scattering22.cpp
          cs22 = scatt22_object.getXSection22();
          probab22 += pow( 0.197, 2.0 ) * cs22 * Vrel * dt / ( dv * testpartcl );
        }
      }
    }
  }
  R22 = probab22 / dt * rings[nc].getGamma();
  
  
  do
  {    
    //------------------------ 2->3 -----------------------
    if( theConfig->doScattering_23() )
    {
      for ( int m1 = 0; m1 < _allParticlesList.size(); m1++ )
      {
        jscat = _allParticlesList[m1];
        if ( !particles_atTimeNow[jscat].dead )
        {
          F2 = particles_atTimeNow[jscat].FLAVOR;
          M2 = particles_atTimeNow[jscat].m;

          s = (addedParticles[jetID].Mom + particles_atTimeNow[jscat].Mom).M2();
          
          if ( s > 1.1*lambda2 )
          {
            Vrel = VelRel( addedParticles[jetID].Mom, particles_atTimeNow[jscat].Mom, M1,M2 );
          
            md2g_wo_as = ( addedParticles[jetID].md2g + particles_atTimeNow[jscat].md2g ) / 2.0;
            md2q_wo_as = ( addedParticles[jetID].md2q + particles_atTimeNow[jscat].md2q ) / 2.0;

            n23++;
            lambda_scaled = ( lambda / 0.197 * sqrt( s ) ); // lambda in fm, sqrt(s) in GeV, lambda_scaled dimensionless
            
            betaDistEntry = scatt23_object.setParameter( rings[nc].getAveraged_v(),
                                                         addedParticles[jetID].Mom,
                                                         particles_atTimeNow[jscat].Mom,
                                                         F1, F2, M1, M2, sqrt( s ), 
                        md2g_wo_as / s, lambda_scaled, 
                        theConfig->getK23LightPartons(), theConfig->getK23HeavyQuarks(),
                        theConfig->getKappa23LightPartons(), theConfig->getKappa23HeavyQuarks(),
                        theConfig->I23onlineIntegrationIsSet(),
                        theConfig->get23GluonFormationTimeTyp(), theConfig->getMatrixElement23(), theConfig->isMd2CounterTermInI23(), theConfig->get23FudgeFactorLpm(), _gluonList.size(), theConfig->isMatrixElement23_22qt() );    
            cs23 = scatt23_object.getXSection23( initialStateIndex ); //1/GeV^2

            probab23 += pow( 0.197, 2.0 ) * cs23 * Vrel * dt / ( dv * testpartcl );
            
          }
        }
      }
    }
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
          
          s = ( addedParticles[jetID].Mom + particles_atTimeNow[iscat].Mom + particles_atTimeNow[jscat].Mom ).M2();
          
          if ( s > 1.1*lambda2 )
          {
            n32++;
            lambda_scaled = ( lambda / 0.197 * sqrt( s ) ); // lambda in fm, sqrt(s) in GeV, lambda_scaled dimensionless
            
            md2g_wo_as = ( particles_atTimeNow[iscat].md2g + particles_atTimeNow[jscat].md2g + addedParticles[jetID].md2g ) / 3.0;
            md2q_wo_as = ( particles_atTimeNow[iscat].md2q + particles_atTimeNow[jscat].md2q + addedParticles[jetID].md2q ) / 3.0;
            
            // create scattering32 object for the given 3 particles
            betaDistEntry = scatt32_object.setParameter( rings[nc].getAveraged_v(),
                                                         addedParticles[jetID].Mom,
                                                         particles_atTimeNow[iscat].Mom,
                                                         particles_atTimeNow[jscat].Mom,
                                                         F1, F2, F3, sqrt( s ), md2g_wo_as * coupling::get_constant_coupling() / s, lambda_scaled, coupling::get_constant_coupling(), false, false, theConfig->get23FudgeFactorLpm(), _gluonList.size() );
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
                F3 = particles_atTimeNow[jscat].FLAVOR;

                // at least one of the particles must be a gluon
                // otherwise go to next step in the loop
                if ( !( F1 == gluon || F2 == gluon || F3 == gluon ) )
                {
                  continue;
                }

                s = ( addedParticles[jetID].Mom + particles_atTimeNow[iscat].Mom + particles_atTimeNow[jscat].Mom ).M2();
                
                if ( s > 1.1*lambda2 )
                {
                  n32++;
                  lambda_scaled = ( lambda / 0.197 * sqrt( s ) ); // lambda in fm, sqrt(s) in GeV, lambda_scaled dimensionless
                  
                  md2g_wo_as = ( particles_atTimeNow[iscat].md2g + particles_atTimeNow[jscat].md2g + addedParticles[jetID].md2g ) / 3.0;
                  md2g_wo_as = ( particles_atTimeNow[iscat].md2q + particles_atTimeNow[jscat].md2q + addedParticles[jetID].md2q ) / 3.0;
                  // HACK for N_f = 0 background
  //                 md2g = as * ( particles_atTimeNow[iscat].md2g + particles_atTimeNow[jscat].md2g + addedParticles[jetID].md2g ) / 3.0 * 2;
  //                 md2q = md2g * 2.0 / 9.0;
                  
                  // create scattering32 object for the given 3 particles
                  betaDistEntry = scatt32_object.setParameter( rings[nc].getAveraged_v(),
                                                               addedParticles[jetID].Mom,
                                                               particles_atTimeNow[iscat].Mom,
                                                               particles_atTimeNow[jscat].Mom,
                                                               F1, F2, F3, sqrt( s ), md2g_wo_as * coupling::get_constant_coupling() / s, lambda_scaled, coupling::get_constant_coupling(), false, false, theConfig->get23FudgeFactorLpm(), _gluonList.size() );
                  I32 = scatt32_object.getIntegral32_withPrefactors();                        // get the integral I32 for the given 3 particles
                  
                  probab32 += I32 * dt / ( pow( dv, 2.0 ) * pow( static_cast<double>( testpartcl ) , 2.0 ) );
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

    // velocity of cell
    VectorXYZ v_cell = rings[nc].getAveraged_v();
              
    // velocity of added particle in lab frame
    VectorTXYZ v_jet = addedParticles[jetID].Mom * (1.0 / addedParticles[jetID].Mom.E());
//     double velocity_lab = sqrt( 1.0 - pow( addedParticles[jetID].m ,2.0) / pow( addedParticles[jetID].E ,2.0) );
    
    // Compute velocity of added particle in rest frame of fluid. The general expression for adding two velocities is not symmetric in v1 and v2. Therefore compute both cases and take average.
    double velocity_rest_1 = addVelocities( v_cell.X(), v_cell.Y(), v_cell.Z(), v_jet.X(), v_jet.Y(), v_jet.Z() );
    double velocity_rest_2 = addVelocities( v_jet.X(), v_jet.Y(), v_jet.Z(), v_cell.X(), v_cell.Y(), v_cell.Z() );
    velocity_jet_restframe = ( velocity_rest_1 + velocity_rest_2 ) / 2.0;
    
    const double epsilon_rate = 1E-3; // 0.001 fm^-1 which corresponds to a lambda = 1000 fm
    // get new lambda which one obtains with the previous employed lambda for the cross section
    // If the cross sections are very small or 0, R22 = R23 = R32 = 0 or very small which causes lambda to be very large. Therefore, epsilon serves as a cut-off.  However, since the cross section is small probably no scattering would take place anyhow. So it does not matter which value is returned...
    if( (R22 + R23 + R32) < epsilon_rate ) 
      lambdaRes = velocity_jet_restframe / epsilon_rate; //fm
    else
      lambdaRes = velocity_jet_restframe / (R22 + R23 + R32); //fm

    // error checking  
    if( std::isnan(lambdaRes) || std::isinf(lambdaRes) || FPT_COMP_E( lambdaRes, 0.0 ) )
    {
      if( std::isinf(lambdaRes) && (R22 + R23 + R32) != 0.0 )
        cout << "lambda = " << lambdaRes  << "  " <<  sqrt( s )  << "  " <<  (R22 + R23 + R32) << endl;
      
      if( std::isnan(lambdaRes) )
        cout << "lambda = " << lambdaRes  << "  " <<  sqrt( s ) << endl;
      
      if( lambdaRes == 0.0 )
        cout << "lambda = " <<  lambdaRes  << "  " <<  sqrt( s ) << endl;
      
      lambdaRes = 0.0001;
      return lambdaRes;
    }
    
    if( lambdaRes > 10000.0 ) // fm
    {
      cout << "lambda is very large = " << lambdaRes  << "  " <<  sqrt( s )  << "  " <<  (R22 + R23 + R32) << endl;

      cout << "lambda = " << lambda << "  lambdaRes = " << lambdaRes << endl;
      cout <<  "  " << velocity_jet_restframe << "  " << R22 <<"  " << R23 <<"  " << R32 <<"  " << endl;
      cout << nTotal << endl;
      //      cout << vx << "  " << vy << "  " << vz << "  " << P1[0] << "  " << P2[0] << "  " << F1 << "  " << F2 << "  " <<  M1 << "  " <<  M2 << "  " <<  sqrt( s ) << "  " <<  md2g_wo_as / s << "  " <<  lambda_scaled << endl;

      lambdaRes = 10000.0;
      return lambdaRes;
    }
    
    // constrain for bisection: we know that the function lambdaRes( lambda ) is monotonly falling
    if( lambdaRes < lambda )
      lambdaMax = lambda;
    else
      lambdaMin = lambda;
    
    
    if (lambdaArray.size() < 4)
    {
      lambdaArray.push_back( lambda );
      lambdaResArray.push_back( lambdaRes );
    }
    else
    {
      lambdaArray.push_back( lambda );
      lambdaResArray.push_back( lambdaRes );
      
      lambdaArray.pop_front();
      lambdaResArray.pop_front();
    }

    if ( lambdaArray.size() == 4 )
    {
      lambdaAvr = 0;
      for ( int m = ( int( lambdaArray.size() ) - 1 ); m >= ( int( lambdaArray.size() ) - 4 ); m-- )
        lambdaAvr += lambdaArray[m];
      lambdaAvr = lambdaAvr / 4;
      
      lambdaResAvr = 0;
      for ( int m = ( int( lambdaResArray.size() ) - 1 ); m >= ( int( lambdaResArray.size() ) - 4 ); m-- )
        lambdaResAvr += lambdaResArray[m];
      lambdaResAvr = lambdaResAvr / 4;
      
      lambda_lambdaRes_converged = fabs( lambdaAvr - lambdaResAvr ) / lambdaAvr < epsilon;
      
      lambda_converged = ( ( fabs( lambdaArray[lambdaArray.size()-1] - lambdaAvr ) / lambdaAvr < epsilon ) &&
                    ( fabs( lambdaArray[lambdaArray.size()-2] - lambdaAvr ) / lambdaAvr < epsilon ) &&
                    ( fabs( lambdaArray[lambdaArray.size()-3] - lambdaAvr ) / lambdaAvr < epsilon ) &&
                    ( fabs( lambdaArray[lambdaArray.size()-4] - lambdaAvr ) / lambdaAvr < epsilon ) );
                    
      lambdaRes_converged = ( ( fabs( lambdaResArray[lambdaResArray.size()-1] - lambdaResAvr ) / lambdaResAvr < epsilon ) &&
                    ( fabs( lambdaResArray[lambdaResArray.size()-2] - lambdaResAvr ) / lambdaResAvr < epsilon ) &&
                    ( fabs( lambdaResArray[lambdaResArray.size()-3] - lambdaResAvr ) / lambdaResAvr < epsilon ) &&
                    ( fabs( lambdaResArray[lambdaResArray.size()-4] - lambdaResAvr ) / lambdaResAvr < epsilon ) );
                    
      converged = lambda_lambdaRes_converged && lambda_converged;
                    
      if( lambda_converged && !converged )
      {
        if( lambdaResAvr > lambdaMax )
        {
          lambdaMax = lambdaResAvr;
//           cout << "rechange max = " << lambdaMax << endl;
          lambdaArray.clear();
          lambdaResArray.clear();
        }
        else if( lambdaResAvr < lambdaMin )
        {
          lambdaMin = lambdaResAvr;
//           cout << "rechange min = " << lambdaMin << endl;
          lambdaArray.clear();
          lambdaResArray.clear();
        }
      }
    }
    else
      converged = false;

    
    lambda = (lambdaMax - lambdaMin) / 2.0 + lambdaMin;
    
    nIterations++;
    probab23 = probab32 = 0;
  }
  while ( !converged && nIterations < nIterationsMax );

  
  lambda = ( lambdaAvr + lambdaResAvr ) / 2.0;
  
  if(!converged)
  {
    double fraction = lambdaResAvr / lambdaAvr * 100.0; // in %
    if( ( fraction < 10.0 || fraction > 1000.0 ) && nTotal > 3 ) // %
//     if( ( fraction < 80.0 || fraction > 120.0 ) && nTotal > 0 ) // %
    {
      cout << "not converged: lambda = " << lambda << " fm   lambdaResAvr = " << fraction << " % of lambdaAvr   nmb in cell: " << nTotal << endl;
//       cout << lambdaAvr << "\t";
//       for( int i = 1; ( i <= 4  && i <= lambdaArray.size() ); i++ )
//         cout << lambdaArray[lambdaArray.size()-i] << "\t";
//       cout << endl;
//       cout << lambdaResAvr << "\t";
//       for( int i = 1; ( i <= 4  && i <= lambdaResArray.size() ); i++ )
//         cout << lambdaResArray[lambdaResArray.size()-i] << "\t";
//       cout << endl;
    }

    if( fabs(fraction) > 1E4 ) // %
    {
      cout << "ultra large lambda, not converged: lambda = " << lambda << " fm   difference = " << fraction << " %   nmb in cell: " << nTotal << endl;
      cout << lambdaAvr << "\t";
      for( int i = 1; ( i <= 4  && i <= lambdaArray.size() ); i++ )
        cout << lambdaArray[lambdaArray.size()-i] << "\t";
      cout << endl;
      cout << lambdaResAvr << "\t";
      for( int i = 1; ( i <= 4  && i <= lambdaResArray.size() ); i++ )
        cout << lambdaResArray[lambdaResArray.size()-i] << "\t";
      cout << endl;
    }
    else
    {
//       n_lambda23NotConverged++;
//       diff_lambda23NotConverged += fabs(fraction);
    }
  }
  
  if( lambda < 0.0001 )
  {
    double fraction = lambdaResAvr / lambdaAvr * 100.0; // in %
    cout << "very small lambda: lambda = " << lambda << " fm   difference = " << fraction << " %   nmb in cell: " << nTotal << "  converged: " << converged << endl;
    cout << lambdaAvr << "\t";
    for( int i = 1; ( i <= 4  && i <= lambdaArray.size() ); i++ )
      cout << lambdaArray[lambdaArray.size()-i] << "\t";
    cout << endl;
    cout << lambdaResAvr << "\t";
    for( int i = 1; ( i <= 4  && i <= lambdaResArray.size() ); i++ )
      cout << lambdaResArray[lambdaResArray.size()-i] << "\t";
    cout << endl;
    
    lambda = 0.0001;
  }

  return lambda; // fm
}


// add two velocity relativistically
// returns the absolut value of the added velocity
// the general expression for this is not symmetric in v1 and v2
double offlineHeavyIonCollision::addVelocities( const double vx_1, const double vy_1, const double vz_1, const double vx_2, const double vy_2, const double vz_2 )
{
  double v1[4], v2[4], v2_parallel[4], v2_perp[4], result[4];
  double scalar_v1_v2, v1_squared, result_scalar;
  
  v1[1] = vx_1;
  v1[2] = vy_1;
  v1[3] = vz_1;
  
  v2[1] = vx_2;
  v2[2] = vy_2;
  v2[3] = vz_2;
  
  scalar_v1_v2 = 0;
  for( int i = 1; i <= 3; i++ )
    scalar_v1_v2 += v1[i] * v2[i];
  
  v1_squared = 0;
  for( int i = 1; i <= 3; i++ )
    v1_squared += v1[i] * v1[i];
  
  if( v1_squared == 0.0 )
    return sqrt( vx_2 * vx_2 + vy_2 * vy_2 + vz_2 * vz_2 );
  
  if( v1_squared >= 1.01 ) // velocity larger than 1
    cout << "error in addVelocities: velocity larger than 1: " << v1_squared << endl;
  else if( v1_squared >= 1.0 ) // avoid rounding errors, set velocity to 1
    v1_squared = 1.0;
  
  
  for( int i = 1; i <= 3; i++ )
    v2_parallel[i] = scalar_v1_v2 / v1_squared * v1[i];
  
  for( int i = 1; i <= 3; i++ )
    v2_perp[i] = v2[i] - v2_parallel[i];
  
  for( int i = 1; i <= 3; i++ )
    result[i] = ( v1[i] + v2_parallel[i] + sqrt( 1.0 - v1_squared ) * v2_perp[i] ) / ( 1.0 + scalar_v1_v2 );
  
  result_scalar = 0;
  for( int i = 1; i <= 3; i++ )
    result_scalar += result[i] * result[i];
  
  return sqrt( result_scalar );
}


// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;
