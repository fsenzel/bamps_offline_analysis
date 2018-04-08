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
  //WARNING Cellcut changed from 4 to 2
  const int cellcut = 2;
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

  long ncoll, ncoll32, ncoll23, ncoll22, ncolle;

  double randomShiftEta, randomShiftX, randomShiftY;
  //  bool dodo1;
  //  int nn_ana1;
  
  double p23_collected_gluon;
  int n23_collected_gluon;
  double p23_collected_quark;
  int n23_collected_quark;
  double lambdaJet_gluon;
  double lambdaJet_quark;
  
  double lambdaT400,countLambdaatT400;
  
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
    theI23_massless( false ), theI23_charm_m1( false ), theI23_charm_m2( false ), theI23_bottom_m1( false ), theI23_bottom_m2( false ), theI23_photons( false ), // do not load data files right at construction, but after configure() has been called below
    offlineInterface( _offlineInterface ),
    theMFP( _config ),
    theAnalysis( _analysis )
{  
  // load 2->2 cross section interpolation data
  // Scattering23_photons has mean free path routine (scatt22ForRates) which uses the I22.
  if( theConfig->doScattering_22() || theConfig->doScattering_22_photons() || theConfig->doScattering_23_photons())
    theI22.configure( theConfig->isCouplingRunning(), Particle::N_light_flavor, Particle::N_heavy_flavor, Particle::Mcharm, Particle::Mbottom, theConfig->getMaxRunningCoupling(), theConfig->getfixedCouplingValue(), theConfig->doScattering_22_photons(), theConfig->getDebyeModePhotons(), theConfig->getVertexModePhotons() );
   
  //Loeschen,alt:
  //if( theConfig->doScattering_22() )
  //  theI22.configure( theConfig->isCouplingRunning(), Particle::N_light_flavor, Particle::N_heavy_flavor, Particle::Mcharm, Particle::Mbottom, theConfig->getMaxRunningCoupling(), theConfig->getfixedCouplingValue() );

  // load 2->3 cross section interpolation data
  if( theConfig->doScattering_23() )
  {
    if( theConfig->getNlightFlavorsAdded() >= 0 )
    {
      theI23_massless.configure( theConfig->I23onlineIntegrationIsSet(),false, 1, 0.0, theConfig->getKappa23LightPartons(), theConfig->get23GluonFormationTimeTyp(), theConfig->getMatrixElement23(), theConfig->isMd2CounterTermInI23(), theConfig->get23FudgeFactorLpm(), theConfig->getInterpolation23Mode(), theConfig->isMatrixElement23_22qt());
    }
    if( Particle::N_heavy_flavor > 0 )
    {
      theI23_charm_m1.configure( theConfig->I23onlineIntegrationIsSet(),false, 1, Particle::Mcharm, theConfig->getKappa23HeavyQuarks(), theConfig->get23GluonFormationTimeTyp(), theConfig->getMatrixElement23(), theConfig->isMd2CounterTermInI23(), theConfig->get23FudgeFactorLpm(), theConfig->getInterpolation23Mode(), theConfig->isMatrixElement23_22qt() );
      theI23_charm_m2.configure( theConfig->I23onlineIntegrationIsSet(),false, 2, Particle::Mcharm, theConfig->getKappa23HeavyQuarks(), theConfig->get23GluonFormationTimeTyp(), theConfig->getMatrixElement23(), theConfig->isMd2CounterTermInI23(), theConfig->get23FudgeFactorLpm(), theConfig->getInterpolation23Mode(), theConfig->isMatrixElement23_22qt() );
    }
    if( Particle::N_heavy_flavor > 1 )
    {
      theI23_bottom_m1.configure( theConfig->I23onlineIntegrationIsSet(),false, 1, Particle::Mbottom, theConfig->getKappa23HeavyQuarks(), theConfig->get23GluonFormationTimeTyp(), theConfig->getMatrixElement23(), theConfig->isMd2CounterTermInI23(), theConfig->get23FudgeFactorLpm(), theConfig->getInterpolation23Mode(), theConfig->isMatrixElement23_22qt() );
      theI23_bottom_m2.configure( theConfig->I23onlineIntegrationIsSet(),false, 2, Particle::Mbottom, theConfig->getKappa23HeavyQuarks(), theConfig->get23GluonFormationTimeTyp(), theConfig->getMatrixElement23(), theConfig->isMd2CounterTermInI23(), theConfig->get23FudgeFactorLpm(), theConfig->getInterpolation23Mode(), theConfig->isMatrixElement23_22qt() );
      
    }
    
    if( theConfig->getJetMfpComputationType() == computeMfpInterpolation || theConfig->getJetMfpComputationType() == thermalMfpGluon )
      theMFP.loadData(); 
  }

  if(theConfig->doScattering_23_photons() && !theConfig->I23onlineIntegrationPhotonsIsSet())
  {
    cout << "Do the 23 photonproduction from Tables." << endl;
    theI23_photons.configurePhotons(theConfig->I23onlineIntegrationIsSet(), theConfig->I23onlineIntegrationPhotonsIsSet(), 1, 0.0, 0, "", "", 0, theConfig->get23FudgeFactorLpm(), theConfig->getInterpolation23Mode(), 0, normal_photons, theConfig->getDebyeModePhotons(),theConfig->getVertexModePhotons(),theConfig->getLPMModePhotons());
    theMFP.loadData();    
  }else if(theConfig->doScattering_23_photons())
  {
    theMFP.loadData(); 
  }
  
  if(theConfig->doScattering_AMY23_photons())
  {
    AMY.configure(0.05);
  }  
  //if( theConfig->doScattering_23_photons() )
  //{
  //  theI23_photons.configure(false,theConfig->I23onlineIntegrationPhotonsIsSet(), 1, 0.0, 0, "", "", 0, theConfig->get23FudgeFactorLpm(), theConfig->getInterpolation23Mode(), 0, normal_photons);
  //  theMFP.loadData(); 
  //}
  
  nGet23Errors = 0;
  nGet32Errors = 0;
  totalPhotonNumber = 0;
  totalDileptonNumber = 0;
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
  
  totalPhotonNumber = 0;
  totalDileptonNumber = 0;
  theConfig->v2average_debug = 0;
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
  bool onlyMediumEvolution=theConfig->isOnlyMediumEvolution();
  
  //TEST
//   lambdaT400=0;
//   countLambdaatT400=0;
  //
  
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
  
  //TEST
  //theAnalysis->printCellV2Distribution(nexttime, simpleTimestepcount);
  //simpleTimestepcount++;
  
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
  
  if(!theConfig->useInitialFormationTimesForPhotonproduction())
  {
    cout << "No formation times for photon scattering" << endl;
  }

  cout << "Maximum Run Time set to " << stoptime << endl;
  if(onlyMediumEvolution)
  {
    cout << "ONLY MEDIUM EVOLUTION" << endl;
  }
  
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

    //
    //cout << "noninteractingParticles.size = " << noninteractingParticles.size() << endl;
    
    if(onlyMediumEvolution)
    {
      // time for next analysis or movie
      double time_ana = theAnalysis->tstep[nn_ana];
      double time_movie = theAnalysis->tstep_movie[nn_ana_movie];
      // determine which one is earlier/smaller
      if( time_movie < time_ana && theConfig->doOutput_movieOutputBackground() )
        simulationTime = time_movie;
      else
        simulationTime = time_ana;

      if(simulationTime > stoptime)
      {        
        break;
      }
      
      // evolution of the medium to this time
      double dt_cascade_from_data = evolveMedium( simulationTime, endOfDataFiles );
      
    }
    else
    {
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
      
    if ( doMovieStepMedium && theConfig->doOutput_movieOutputBackground() )
    {
      theAnalysis->movieOutputMedium( nn_ana_movie - 1, jumpMovieSteps );
      doMovieStepMedium = false;
    }
  
  
    if (!onlyMediumEvolution)
    {
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
      scattering( nexttime, again);
      
//       {
//         double v2avg=0.;
//         int v2avgN=0;
//         double v2,pt,eta;
//         for(int i=0; i< particles_atTimeNow.size(); i++)
//         {
//           pt = particles_atTimeNow[i].Mom.Pt();
//           // for most particle species we are interested in the pseudorapidity (for massless particle it does not matter anyhow)
//           eta = particles_atTimeNow[i].Mom.Rapidity();
//           v2 = ( pow( particles_atTimeNow[i].Mom.Px(), 2.0 ) - pow( particles_atTimeNow[i].Mom.Py(), 2.0 ) ) / pow( pt, 2.0 );
//           if(pt<0.1 && -0.5 < eta && 0.5 > eta )
//           {
//             v2avg+=v2;
//             v2avgN++;
//           }
//         }
//         cout << "Parton v2= " << 100*v2avg/v2avgN << " % N = " << v2avgN  << endl; 
//       }
//       {
//         double v2avg=0.;
//         int v2avgN=0;
//         double v2,pt,eta;
//         
//         for(int i=0; i< noninteractingParticles.size(); i++)
//         {
//           pt = noninteractingParticles[i].Mom.Pt();
//           // for most particle species we are interested in the pseudorapidity (for massless particle it does not matter anyhow)
//           eta = noninteractingParticles[i].Mom.Rapidity();
//           v2 = ( pow( noninteractingParticles[i].Mom.Px(), 2.0 ) - pow( noninteractingParticles[i].Mom.Py(), 2.0 ) ) / pow( pt, 2.0 );
//           if(pt<0.1 && -0.5 < eta && 0.5 > eta )
//           {
//             v2avg+=v2;
//             v2avgN++;
//           }
//         }
//         if(v2avgN>0)
//         cout << "Photon v2= " << 100*v2avg/v2avgN << " %  of N=" << v2avgN << endl;       
//       }
      
      
      //theAnalysis->printCellV2Distribution(nexttime, simpleTimestepcount);
      //simpleTimestepcount++; 
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
        //theAnalysis->printCellV2Distribution(nexttime, simpleTimestepcount);
        //simpleTimestepcount++;
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
        

      // analyse timesteps
      dt_sum += dt;
      n_dt++;
      
      simulationTime = nexttime;
      
    }
    else
    {
      
      /*cell_ID( simulationTime );

      // collide added particles with gluonic medium
      deadParticleList.clear();
    
      scattering( simulationTime, again );*/
      {
        double v2avg=0.;
        int v2avgN=0;
        double v2,pt,eta;
        for(int i=0; i< particles_atTimeNow.size(); i++)
        {
          pt = particles_atTimeNow[i].Mom.Pt();
          // for most particle species we are interested in the pseudorapidity (for massless particle it does not matter anyhow)
          eta = particles_atTimeNow[i].Mom.Rapidity();
          v2 = ( pow( particles_atTimeNow[i].Mom.Px(), 2.0 ) - pow( particles_atTimeNow[i].Mom.Py(), 2.0 ) ) / pow( pt, 2.0 );
          if(pt<0.1 && -0.5 < eta && 0.5 > eta )
          {
            v2avg+=v2;
            v2avgN++;
          }
        }
        cout << "Parton v2= " << 100*v2avg/v2avgN << " % N = " << v2avgN  << endl; 
      }
      
      
      if ( FPT_COMP_E( simulationTime, theAnalysis->tstep[nn_ana] ) ) // ask if it is time for analysis
      {
        doAnalysisStep = true;
        cout << "Analyse: " << simulationTime << endl;
      }
      if ( doAnalysisStep )
      {
        theAnalysis->intermediateOutput( nn_ana );
        theAnalysis->collectPtData( nn_ana );
        theAnalysis->collectYData( nn_ana );
        //theAnalysis->collectEtData( nn_ana );
        nn_ana++;
        doAnalysisStep = false;
      }
    }
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

  cout << "Do final output: " << endl;
  theAnalysis->finalOutput( stoptime );
  //theAnalysis->addJetEvents_final();
  
  //TEST
  //cout << "lambda/fm = " << lambdaT400/countLambdaatT400  << endl;  
  //
  
  
  if(!onlyMediumEvolution)
  {
    cout << "Photon spectra printed" << endl;
    theAnalysis->photonSpectrumOutput();
    
    cout << "number of errors in get32(...) = " << nGet32Errors << endl;
    cout << "number of errors in get23(...) = " << nGet23Errors << endl;
    
    cout << "==========================" << endl;
    cout << "Total Number of Photons = " << totalPhotonNumber << endl;
    cout << "Total Number of Dilepton Pairs = " << totalDileptonNumber << endl;
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
    actiontype = event_dummy;
    try
    {
      //Kai's idea:
      //actiontype = offlineInterface->readOfflineDataFromArchive< offlineDataEventType >()->event;
      
      boost::shared_ptr< offlineDataEventType > ptrEventType = offlineInterface->readOfflineDataFromArchive< offlineDataEventType >();
      actiontype = ptrEventType->event;
    }    
    catch ( boost::archive::archive_exception& err )
    {
      cout << " try catch for actiontype = offlineInterface->readOfflineDataFromArchive< offlineDataEventType >()->event;" << endl;
      stop = true;
      _endOfDataFiles = true;
      actiontype = event_dummy;
    } catch(...) 
    {
      // this executes if f() throws std::string or int or any other unrelated type
      cout << " unknown try catch for actiontype = offlineInterface->readOfflineDataFromArchive< offlineDataEventType >()->event;" << endl;
      stop = true;
      _endOfDataFiles = true;
      actiontype = event_dummy;
    }
    

    if ( actiontype == event_newTimestep )
    {
      try
      {
        boost::shared_ptr< offlineDataCellConfiguration > ptrCellStructure = offlineInterface->readOfflineDataFromArchive< offlineDataCellConfiguration >();
        etaBins = ptrCellStructure->etaBins;
        timenow = ptrCellStructure->timenow;
        timenext = ptrCellStructure->timenext;
        randomShiftX = ptrCellStructure->randomShiftX;
        randomShiftY = ptrCellStructure->randomShiftY;
        randomShiftEta = ptrCellStructure->randomShiftEta;
          
        dt_cascade = timenext - timenow;      
      }
      catch( boost::archive::archive_exception& err )
      {
        cout << "try catch for actiontype == event_newTimestep" << endl;
      }
      rateGluons_prev = rateGluons;
      rateQuarks_prev = rateQuarks;
      rateAntiQuarks_prev = rateAntiQuarks;
      for ( int i = 0; i < IZ; i++ )
      {
        rateGluons[i].assign( rings.size(), 0 );
        rateQuarks[i].assign( rings.size(), 0 );
        rateAntiQuarks[i].assign( rings.size(), 0 );
      }
      try
      {
        boost::shared_ptr< offlineDataInteractionRates > ptrRates = offlineInterface->readOfflineDataFromArchive< offlineDataInteractionRates >();
        rateGluons = ptrRates->gluonRates;
        rateQuarks = ptrRates->quarkRates;
        rateAntiQuarks = ptrRates->antiQuarkRates;              
      }
      catch( boost::archive::archive_exception& err )
      {
        cout << "try catch for boost::shared_ptr< offlineDataInteractionRates > ptrRates = offlineInterface->readOfflineDataFromArchive< offlineDataInteractionRates >();" << endl;
      }
    }
    else if ( actiontype == event_interaction22 )
    {
      try
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
      catch( boost::archive::archive_exception& err )
      {
        cout << "try catch for boost::shared_ptr< offlineDataInteraction22 > ptrInteraction22 = offlineInterface->readOfflineDataFromArchive< offlineDataInteraction22 >();" << endl;
      }  
    }
    else if ( actiontype == event_interaction23 )
    {
      try
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
      catch( boost::archive::archive_exception& err )
      {
        cout << "try catch for boost::shared_ptr< offlineDataInteraction23 > ptrInteraction23 = offlineInterface->readOfflineDataFromArchive< offlineDataInteraction23 >();" << endl;
      }
    }
    else if ( actiontype == event_interaction32 )
    {
      try
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
      catch( boost::archive::archive_exception& err )
      {
        cout << "try catch for boost::shared_ptr< offlineDataInteraction32 > ptrInteraction32 = offlineInterface->readOfflineDataFromArchive< offlineDataInteraction32 >();" << endl;
      }
    }
    else if ( actiontype == event_interactionElastic )
    {
      try
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
      catch( boost::archive::archive_exception& err )
      {
        cout << "try catch for boost::shared_ptr< offlineDataInteractionElastic > ptrInteractionElastic = offlineInterface->readOfflineDataFromArchive< offlineDataInteractionElastic >();" << endl;
      }
    }
    else if ( actiontype == event_particleIdSwap )
    {
      try
      {
        boost::shared_ptr< offlineDataParticleIdSwap > ptrSwap = offlineInterface->readOfflineDataFromArchive< offlineDataParticleIdSwap >();
        iscat = ptrSwap->removedParticleID;
        jscat = ptrSwap->replacingParticleID;
      }
      catch( boost::archive::archive_exception& err )
      {
        cout << "try catch for boost::shared_ptr< offlineDataParticleIdSwap > ptrSwap = offlineInterface->readOfflineDataFromArchive< offlineDataParticleIdSwap >();" << endl;
      }
      particlesEvolving[iscat] = particlesEvolving[jscat];
    }
  }


  for ( int i = 0; i < numberEvolvingParticles; i++ )
  {
    if ( particlesEvolving[i].Pos.T() < evolveToTime + 1.0e-8 )
    {
      particlesEvolving[i].init = false;
    }
    if( (theConfig->doScattering_AMY23_photons() || theConfig->doScattering_22_photons() || theConfig->doScattering_23_photons() ) && !(theConfig->useInitialFormationTimesForPhotonproduction()) )
    {
      particlesEvolving[i].init = false;
      particlesEvolving[i].Pos.T() = evolveToTime - 1.0e-6 ;
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


// void offlineHeavyIonCollision::computeVelocityOfCells(int cell_index)
// {
//   list<int>::iterator iIt;
//   VectorEPxPyPz sumofMomentaIncell;
//   VectorXYZ boostbetaVector;
//   for ( iIt = cells[cell_index].particleList.begin(); iIt != cells[cell_index].particleList.end(); iIt++ )
//   {
//     int id = *iIt;
//     sumofMomentaIncell += particles_atTimeNow[id].Mom;
//   }
//   cells[cell_index].boostLRFvelocity = sumofMomentaIncell.NormalizeToE();
// }


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
    if ( particles_atTimeNow[i].Pos.T() < _time || !(theConfig->useInitialFormationTimesForPhotonproduction()) )
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
  vector<int> gluonList, allParticlesList, allParticlesListWithNeighbors;
  vector<int> gluonListAdded, allParticlesListAdded;
  
  again = false;
  
  IXY = IX * IY;
  
  //used for MFP analysis:
  bool computeAveragedMfp = theConfig->getMfpCellAveraging();
  int  numberOfCellsAveragedForMfp = 5;
  double lambdaSpecific1=0.;
  double lambdaSpecific2=0.;
  double lambdaSpecific3=0.;
  long int count_for_average_lambda1 = 0;
  long int count_for_average_lambda2 = 0;
  long int count_for_average_lambda3 = 0;
  long int averageNumberOfParticles  = 0;
  long int count_for_average_NumberOfParticles = 0;
  averageQuarkNumber = 0.;
  countForAverageaverageQuarkNumber=0;
  
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
  
  //v2 of partons in this routine:    
  double v2avg=0.;
  int v2avgN=0;
  
  double Tepsilonthreep=0.;
  double TAMY=0.;
  int numberTotalForTemp=0;
  
  
  
  
  // go through cells
  for ( int etaSliceIndex = etaBins.min_index(); etaSliceIndex <= etaBins.max_index(); etaSliceIndex++ )
  {
    //cout << "Eta slice: " << etaSliceIndex << endl;
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
    //For each eta-slice, the cells have the same volume.
    dv = dx * dy * dz;
    //Loop through cells
    for ( int j = IXY * etaSliceIndex; j < IXY * ( etaSliceIndex + 1 ); j++ )
    { 
      
      /*if ( !( cells[j].empty()) )
      { 
        cells[j].resetStoredValues();
        cout << cells[j].size() << endl;
        allParticlesListWithNeighbors.clear(); 
        allParticlesListWithNeighbors.reserve( cells[j].size() );
      }*/

      
      int EdgecellOfList = 0;
      if(computeAveragedMfp==true)
      {      
        allParticlesListWithNeighbors.clear();  
        
        if( (j== IXY * (etaSliceIndex+1) - 1 ) && (!( cells[j].empty())  || !( cells[j-1].empty())  ||  !( cells[j-2].empty())  ||  !( cells[j-3].empty())  ||  !( cells[j-4].empty())))
        { 
          numberOfCellsAveragedForMfp = 5;
          EdgecellOfList = 2;
          allParticlesListWithNeighbors.reserve( round(  (cells[j-4].size()+cells[j-3].size()+cells[j-2].size() + cells[j-1].size()+cells[j].size())*1.1 ) );
          
        }else
        if(  (j== IXY *  etaSliceIndex  )  && ( !( cells[j].empty()) || !( cells[j+1].empty()) || !( cells[j+2].empty()) ||  !( cells[j+3].empty())  ||  !( cells[j+4].empty())))
        {
          numberOfCellsAveragedForMfp = 5;
          EdgecellOfList = 1;
          allParticlesListWithNeighbors.reserve( round(   (cells[j].size() + cells[j+1].size()+cells[j+2].size()+cells[j+3].size()+cells[j+4].size())*1.1 ) );
          
        }else  
        if( (j>IXY * etaSliceIndex) && (j< (IXY * ( etaSliceIndex + 1 ) - 1) ) && (  !( cells[j].empty())  ||  !( cells[j-1].empty())  ||  !( cells[j+1].empty())  ||  !( cells[j-2].empty()) ||  !( cells[j+2].empty()) )      )
        {  
          numberOfCellsAveragedForMfp = 5;
          EdgecellOfList = 0;   
          allParticlesListWithNeighbors.reserve( round(  (cells[j-2].size()+cells[j-1].size() + cells[j].size()+cells[j+1].size()+cells[j+2].size())*1.1  ) );
        }else
        {  
          numberOfCellsAveragedForMfp = 1;         
          allParticlesListWithNeighbors.reserve(cells[j].size());
        }
        
      }
      
      //WARNING! For the use without added particles.
      // if ( !( cells[j].empty() || cellsAdded[j].empty() ) )
      if ( !( cells[j].empty()) )
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
        int numberInThisCell=0;
        for ( iIt = cells[j].particleList.begin(); iIt != cells[j].particleList.end(); iIt++ )
        {
          
          
          
          int id = *iIt;
          xt = particles_atTimeNow[id].Pos.Perp();
          numberInThisCell++;
         
          
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
        //cout << "Number of Particles in This Cell = " << numberInThisCell << endl;
        
        for ( iIt = cellsAdded[j].particleList.begin(); iIt != cellsAdded[j].particleList.end(); iIt++ )
        {
          cout << "Added?" << endl;
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
          //cout << "All free" << endl;
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
          //WARNING Cellcut? included 11.4.
          if (( cells[j].size() ) >= cellcut )     // enough particles in cell -> scatter
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
            
            
            if(computeAveragedMfp)
            {
              if(EdgecellOfList==0)
              {
              list<int>::const_iterator iIt;
                for ( iIt = cells[j].particleList.begin(); iIt != cells[j].particleList.end(); iIt++ )
                {
                  id = *iIt;
                  allParticlesListWithNeighbors.push_back( id );
                }       
                for ( iIt = cells[j-1].particleList.begin(); iIt != cells[j-1].particleList.end(); iIt++ )
                {
                  id = *iIt;
                  allParticlesListWithNeighbors.push_back( id );
                }  
                for ( iIt = cells[j+1].particleList.begin(); iIt != cells[j+1].particleList.end(); iIt++ )
                {
                  id = *iIt;
                  allParticlesListWithNeighbors.push_back( id );
                }  
                for ( iIt = cells[j-2].particleList.begin(); iIt != cells[j-2].particleList.end(); iIt++ )
                {
                  id = *iIt;
                  allParticlesListWithNeighbors.push_back( id );
                }  
                for ( iIt = cells[j+2].particleList.begin(); iIt != cells[j+2].particleList.end(); iIt++ )
                {
                  id = *iIt;
                  allParticlesListWithNeighbors.push_back( id );
                }  
              }
              if(EdgecellOfList==1)
              {
              list<int>::const_iterator iIt;
                for ( iIt = cells[j].particleList.begin(); iIt != cells[j].particleList.end(); iIt++ )
                {
                  id = *iIt;
                  allParticlesListWithNeighbors.push_back( id );
                }       
                for ( iIt = cells[j+1].particleList.begin(); iIt != cells[j+1].particleList.end(); iIt++ )
                {
                  id = *iIt;
                  allParticlesListWithNeighbors.push_back( id );
                }  
                for ( iIt = cells[j+2].particleList.begin(); iIt != cells[j+2].particleList.end(); iIt++ )
                {
                  id = *iIt;
                  allParticlesListWithNeighbors.push_back( id );
                }  
                for ( iIt = cells[j+3].particleList.begin(); iIt != cells[j+3].particleList.end(); iIt++ )
                {
                  id = *iIt;
                  allParticlesListWithNeighbors.push_back( id );
                }  
                for ( iIt = cells[j+4].particleList.begin(); iIt != cells[j+4].particleList.end(); iIt++ )
                {
                  id = *iIt;
                  allParticlesListWithNeighbors.push_back( id );
                }  
              }           
              if(EdgecellOfList==2)
              {
              list<int>::const_iterator iIt;
                for ( iIt = cells[j].particleList.begin(); iIt != cells[j].particleList.end(); iIt++ )
                {
                  id = *iIt;
                  allParticlesListWithNeighbors.push_back( id );
                }       
                for ( iIt = cells[j-1].particleList.begin(); iIt != cells[j-1].particleList.end(); iIt++ )
                {
                  id = *iIt;
                  allParticlesListWithNeighbors.push_back( id );
                }  
                for ( iIt = cells[j-2].particleList.begin(); iIt != cells[j-2].particleList.end(); iIt++ )
                {
                  id = *iIt;
                  allParticlesListWithNeighbors.push_back( id );
                }  
                for ( iIt = cells[j-3].particleList.begin(); iIt != cells[j-2].particleList.end(); iIt++ )
                {
                  id = *iIt;
                  allParticlesListWithNeighbors.push_back( id );
                }  
                for ( iIt = cells[j-4].particleList.begin(); iIt != cells[j-2].particleList.end(); iIt++ )
                {
                  id = *iIt;
                  allParticlesListWithNeighbors.push_back( id );
                }  
              }                
            }
            
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
              
              xt = particles_atTimeNow[id].Pos.Perp();
              nc = rings.getIndex( xt );
              
              particles_atTimeNow[id].md2g = rings[nc].getAveraged_md2g();
              particles_atTimeNow[id].md2q = rings[nc].getAveraged_md2q();
              particles_atTimeNow[id].temperature = rings[nc].getEffectiveTemperature();
              particles_atTimeNow[id].temperatureAMY = rings[nc].getEffectiveTemperatureFirstOverZerothMoment();
              
              
              if(etaSliceIndex == centralEtaIndex)
              {
                            
                //TEST
                Tepsilonthreep+=rings[nc].getEffectiveTemperature();
                TAMY+=rings[nc].getEffectiveTemperatureFirstOverZerothMoment();
                numberTotalForTemp+=1.;
                //
                
                //cout<< "Gamma: "<< rings[nc].getGamma() << "\t" << "sqrt(Debye-Mass^2*pi/8/(Nc+Nf)): "<< sqrt(particles_atTimeNow[id].md2g*M_PI/8.0/(Particle::N_light_flavor + Ncolor))<<"   vs Teff= " << particles_atTimeNow[id].temperature<< "\t" << particles_atTimeNow[id].temperature / sqrt(particles_atTimeNow[id].md2g*M_PI/8.0/(Particle::N_light_flavor + Ncolor))<< endl;
              }

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
            
            
            unsigned int tempPhotons = totalPhotonNumber;
            if( theConfig->isScatt_amongBackgroundParticles() )
            {
                if( theConfig->doScattering_22_photons())// && nexttime >0.2
                { 
                  scatt22_amongBackgroundParticles_photons( cells[j], allParticlesList, scaleFactor, again, nexttime );
//                   for ( int i = 0; i < static_cast<int>( allParticlesList.size() ) - 1; i++ )
//                   {
//                     int iscat = allParticlesList[i]; 
//                     double v2,pt,eta;  
//                     pt = particles_atTimeNow[iscat].Mom.Pt();
//                     // for most particle species we are interested in the pseudorapidity (for massless particle it does not matter anyhow)
//                     eta = particles_atTimeNow[iscat].Mom.Rapidity();
//                     v2 = ( pow( particles_atTimeNow[iscat].Mom.Px(), 2.0 ) - pow( particles_atTimeNow[iscat].Mom.Py(), 2.0 ) ) / pow( pt, 2.0 );
//                     if(pt<0.1 && -0.5 < eta && 0.5 > eta )
//                     {
//                       v2avg+=v2;
//                       v2avgN++;
//                     }
//                   }
                }
                
                if( theConfig->doScattering_22_dileptons())
                {
                  scatt22_amongBackgroundParticles_dileptons( cells[j], allParticlesList, scaleFactor, again, nexttime );
                }        
                
                if( theConfig->doScattering_AMY23_photons() )
                {
                  scatt23_amongBackgroundParticles_AMYphotons(cells[j], allParticlesList, scaleFactor, again, nexttime);
                }
                
                if( theConfig->doScattering_23_photons())
                {
                  cells[j].rates.clearSpecific();
                  
                  if (theConfig->getOutputScheme()==203 )                  
                  {     
                    if(etaSliceIndex == centralEtaIndex)
                    {
                      if(computeAveragedMfp)
                      {  
                        averageNumberOfParticles += allParticlesListWithNeighbors.size();
                      }else
                      {
                        averageNumberOfParticles += allParticlesList.size();
                      }
                      count_for_average_NumberOfParticles++;
                  
                      if(computeAveragedMfp)
                      {
                        scatt22ForRates(cells[j],allParticlesListWithNeighbors, theI22, numberOfCellsAveragedForMfp);
                      }else
                      {
                        scatt22ForRates(cells[j],allParticlesList, theI22, 1);                    
                      }
        
                      if(FPT_COMP_GZ(getLambdaFromRates( 1,1, cells[j].rates)))
                      {
                        lambdaSpecific1 += getLambdaFromRates( 1,1, cells[j].rates);
                        count_for_average_lambda1++;
                      }
                      if(FPT_COMP_GZ(getLambdaFromRates( 1,2, cells[j].rates)))
                      {                 
                        lambdaSpecific2 += getLambdaFromRates( 1,2, cells[j].rates);
                        count_for_average_lambda2++;
                      }
                      if(FPT_COMP_GZ(getLambdaFromRates( 1,3, cells[j].rates)))
                      {                
                        lambdaSpecific3 += getLambdaFromRates( 1,3, cells[j].rates);
                        count_for_average_lambda3++;    
                      }   
                    }
                  }else
                  {
                    if(computeAveragedMfp)
                    {
                      scatt22ForRates(cells[j],allParticlesListWithNeighbors, theI22, numberOfCellsAveragedForMfp);
                    }else
                    {
                      scatt22ForRates(cells[j],allParticlesList, theI22, 1);                    
                    }
                    /*
                    if(FPT_COMP_GZ(getLambdaFromRates( 1,1, cells[j].rates)))
                    {
                      lambdaSpecific1 += getLambdaFromRates( 1,1, cells[j].rates);
                      count_for_average_lambda1++;
                    }
                    if(FPT_COMP_GZ(getLambdaFromRates( 1,2, cells[j].rates)))
                    {                 
                      lambdaSpecific2 += getLambdaFromRates( 1,2, cells[j].rates);
                      count_for_average_lambda2++;
                    }
                    if(FPT_COMP_GZ(getLambdaFromRates( 1,3, cells[j].rates)))
                    {                
                      lambdaSpecific3 += getLambdaFromRates( 1,3, cells[j].rates);
                      count_for_average_lambda3++;    
                    } 
                    cout << "Scatter " << lambdaSpecific1/count_for_average_lambda1 << "\t" << lambdaSpecific2/count_for_average_lambda2 << "\t" << lambdaSpecific3/count_for_average_lambda3 << endl;
                    */
                    scatt23_amongBackgroundParticles_photons( cells[j], allParticlesList, scaleFactor, again, nexttime );
                  }
                }             
                

                tempPhotons = totalPhotonNumber - tempPhotons;
                
                theAnalysis->cellV2Distribution(  computeBackgroundv2OfCell(allParticlesList), tempPhotons  ); 
                             
              if ( again )
              {
                //cout << "AGAIN!" << endl;
                return;
              }
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

  
  
  
  //cout << "Scattering Parton v2= " << 100*v2avg/v2avgN << " % N = " << v2avgN  << endl; 

  cout << "T1= " << Tepsilonthreep/numberTotalForTemp << " T2= " << TAMY/numberTotalForTemp  << endl;
  
  
  
  
  formGeomCopy.clear();
  
  // J/psi dissociation: if temperature in cell is higher than Td = 2 Tc, decay J/psi to two charm quarks
  if( Particle::N_psi_states > 0 )
    jpsi_dissociation_td( nexttime );
     
  if(theConfig->getOutputScheme() == 203)
  {
    string full_filename = theConfig->getStandardOutputDirectoryName() + "/" +theConfig->getJobName() + "_MFPaveraged" +  ".dat";
    fstream outfile( full_filename.c_str(), ios::out | ios::app);
    outfile.precision( 8 );
    int avQNumber = 0; 
    if(countForAverageaverageQuarkNumber>0)
    {
      avQNumber = averageQuarkNumber/countForAverageaverageQuarkNumber;
    }
    if(count_for_average_NumberOfParticles>0)
    {
      cout    << nexttime << "\t" << lambdaSpecific1/count_for_average_lambda1 << "\t" << lambdaSpecific2/count_for_average_lambda2<< "\t" <<  lambdaSpecific3/count_for_average_lambda3 << "\t" << averageNumberOfParticles/count_for_average_NumberOfParticles << "\t" << avQNumber << endl;
      outfile << nexttime << "\t" << lambdaSpecific1/count_for_average_lambda1 << "\t" << lambdaSpecific2/count_for_average_lambda2<< "\t" <<  lambdaSpecific3/count_for_average_lambda3 << "\t" << averageNumberOfParticles/count_for_average_NumberOfParticles << "\t" << avQNumber << endl;
      outfile.close();
    }
  }
  
//TEST  
  cout << "Total # of photons " << totalPhotonNumber << endl;
  cout << "Total # of dileptons " << totalDileptonNumber << endl;

  
  
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

  scattering23 scatt23_object( &theI23_massless, &theI23_charm_m1, &theI23_charm_m2, &theI23_bottom_m1, &theI23_bottom_m2 , &theI23_photons);
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
    
    if ( pt_addedParticle < theConfig->getMinimumPT() || addedParticles[jscat].dead )
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

double offlineHeavyIonCollision::computeBackgroundv2OfCell( std::vector< int >& allParticlesList )
{
  int iscat;
  double v2=0.0;
  for ( int i = 0; i <  allParticlesList.size(); i++ )
  {
      iscat = allParticlesList[i];     
      if(fabs(particles_atTimeNow[iscat].Mom.Rapidity())<0.35 && particles_atTimeNow[iscat].Mom.Pt()>1.5 )
      {
        v2+=(pow(particles_atTimeNow[iscat].Mom.Px(),2.0)-pow(particles_atTimeNow[iscat].Mom.Py(),2.0))/pow(particles_atTimeNow[iscat].Mom.Pt(),2.0);
      }
  }
  return v2/allParticlesList.size();   
}



void offlineHeavyIonCollision::scatt22_amongBackgroundParticles_photons( cellContainer& _cells, std::vector< int >& _allParticlesList, const double scaleFactor, bool& again, const double nexttime )
{
  int iscat, jscat;
  unsigned int m1, m2;

  scattering22 scatt22_object( &theI22 );
  
  const int N = _allParticlesList.size();
  //this can be adjusted!
  const int consideredPairs = N;     //number of pairs to be considered (on the average)
  // The cross section needs to be scaled since only a certain number "consideredPairs" of all pairs
  // is taken for simulation. The scale factor depends on the inital state (gg->X has a different number of pairs than gq->X etc.)
  double scaleForSelectedPairs = 0;
  const int allPairs = binomial ( N, 2 );
  
  if ( allPairs >= 5000 )
  {   //only a certain number consideredPairs of all pairs will be considered  
    for ( int i = 0; i < consideredPairs; i++ )
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
      do
      {
        m2 = int ( _allParticlesList.size() * ran2() );
        if ( m2 == _allParticlesList.size() ) 
        {
          m2 = _allParticlesList.size() - 1;
        }
        jscat = _allParticlesList[m2]; 
      }
      while ( jscat == iscat || particles_atTimeNow[jscat].dead );
      
      if(theConfig->getRestrictParentPTForPhotons())
      {
        if(parentParticlesAllowed(iscat,jscat))
        {
          scaleForSelectedPairs =  static_cast<double> ( allPairs ) / static_cast<double> ( consideredPairs );
          scatt22_amongBackgroundParticles_photons_utility_1(scatt22_object, iscat, jscat, nexttime,  scaleForSelectedPairs, again);    
        }else
        {
          continue;
        }
      }else
      {
        scaleForSelectedPairs =  static_cast<double> ( allPairs ) / static_cast<double> ( consideredPairs );
        scatt22_amongBackgroundParticles_photons_utility_1(scatt22_object, iscat, jscat, nexttime,  scaleForSelectedPairs, again);         
      }
    }                                                    
  }
  else
  {
    for ( int i = 0; i < static_cast<int>( _allParticlesList.size() ) - 1; i++ )
    {
      iscat = _allParticlesList[i];    
      for ( unsigned int j = i+1; j < _allParticlesList.size(); j++ )
      {              
        jscat = _allParticlesList[j]; 
        scaleForSelectedPairs=1.0;  
        if(particles_atTimeNow[iscat].dead || particles_atTimeNow[jscat].dead)
        {
          continue;
        }else
        {
          if(theConfig->getRestrictParentPTForPhotons())
          {
            if(parentParticlesAllowed(iscat,jscat))
            {
              scatt22_amongBackgroundParticles_photons_utility_1(scatt22_object, iscat, jscat, nexttime,  scaleForSelectedPairs, again);    
            }else
            {
              continue;
            }
          }else
          {
            scatt22_amongBackgroundParticles_photons_utility_1(scatt22_object, iscat, jscat, nexttime,  scaleForSelectedPairs, again);         
          }
        } 
      }
    }
  }
}


/**
* Performs the 2->2 dilepton production for all particles in a given computational cell
*
* @param[in] cellMembers list with all particles in the cell
* @param[in] dt size of time step
* @param[in] time current time
* @param[in] factor scale factor = ngluon / (ngluon - n32) needed to correct rates in #scatt23 and #scatt22 after #scatt32
* @param[in] aa Reference to the analysis object, used for output operations and analyis.
*/
void offlineHeavyIonCollision::scatt22_amongBackgroundParticles_dileptons( cellContainer& _cells, std::vector< int >& _allParticlesList, const double scaleFactor, bool& again, const double nexttime )
{  
  int iscat, jscat;
  unsigned int m1, m2;
  
  scattering22 scatt22_object( &theI22 );
  
  const int N = _allParticlesList.size();
  //this can be adjusted!
  const int consideredPairs = N;     //number of pairs to be considered (on the average)
  // The cross section needs to be scaled since only a certain number "consideredPairs" of all pairs
  // is taken for simulation. The scale factor depends on the inital state (gg->X has a different number of pairs than gq->X etc.)
  double scaleForSelectedPairs = 0;
  const int allPairs = binomial ( N, 2 ); 

  if ( allPairs >= 50 )
  {   //only a certain number consideredPairs of all pairs will be considered  
    for ( int i = 0; i < consideredPairs; i++ )
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
      do
      {
        m2 = int ( _allParticlesList.size() * ran2() );
        if ( m2 == _allParticlesList.size() ) 
        {
          m2 = _allParticlesList.size() - 1;
        }
        jscat = _allParticlesList[m2]; 
      }
      while ( jscat == iscat || particles_atTimeNow[jscat].dead );
      
      if(theConfig->getRestrictParentPTForPhotons())
      {
        if(parentParticlesAllowed(iscat,jscat))
        {
          scaleForSelectedPairs =  static_cast<double> ( allPairs ) / static_cast<double> ( consideredPairs );
          scatt22_amongBackgroundParticles_dileptons_utility_1(scatt22_object, iscat, jscat, nexttime,  scaleForSelectedPairs, again);    
        }else
        {
          continue;
        }
      }else
      {
        scaleForSelectedPairs =  static_cast<double> ( allPairs ) / static_cast<double> ( consideredPairs );
        scatt22_amongBackgroundParticles_dileptons_utility_1(scatt22_object, iscat, jscat, nexttime,  scaleForSelectedPairs, again);         
      }
    }                                                    
  }
  else
  {
    for ( int i = 0; i < static_cast<int>( _allParticlesList.size() ) - 1; i++ )
    {
      iscat = _allParticlesList[i];    
      for ( unsigned int j = i+1; j < _allParticlesList.size(); j++ )
      {              
        jscat = _allParticlesList[j]; 
        scaleForSelectedPairs=1.0;  
        if(particles_atTimeNow[iscat].dead || particles_atTimeNow[jscat].dead)
        {
          continue;
        }else
        {
          if(theConfig->getRestrictParentPTForPhotons())
          {
            if(parentParticlesAllowed(iscat,jscat))
            {
              scatt22_amongBackgroundParticles_dileptons_utility_1(scatt22_object, iscat, jscat, nexttime,  scaleForSelectedPairs, again);    
            }else
            {
              continue;
            }
          }else
          {
            scatt22_amongBackgroundParticles_dileptons_utility_1(scatt22_object, iscat, jscat, nexttime,  scaleForSelectedPairs, again);         
          }
        } 
      }
    }
  }
}

/**
 * Performs the inelastic photonproduction.
 * 
 * @param[in] _cells
 * @param[in] _allParticlesList
 * @param[in] scaleFactor
 * @param[in] again
 * @param[in] nexttime
 * 
 * 
 */
void offlineHeavyIonCollision::scatt23_amongBackgroundParticles_photons( cellContainer& _cells, std::vector< int >& _allParticlesList, const double scaleFactor, bool& again, const double nexttime )
{
  int iscat, jscat, typ,ringIndex;
  double s, cs22, Vrel,temp1,temp2,temp3;
  double M1, M2;
  VectorEPxPyPz P1, P2;
  double md2g_wo_as, md2q_wo_as,s_cutoff_for_pqcd, lambda,lambda_scaled,xt;
  int initialStateIndex = -1;
  FLAVOR_TYPE F1, F2;
  double temperature;
  double cs23Total,cs23Photons,probab23,betaDistEntry;
  unsigned int m1, m2;
  list<int>::const_iterator iIt;
  list<int>::const_iterator jIt;
  
  scattering23 scatt23_object ( &theI23_massless, &theI23_charm_m1, &theI23_charm_m2, &theI23_bottom_m1, &theI23_bottom_m2, &theI23_photons );

  const int N = _allParticlesList.size();
  //this can be adjusted!
  const int consideredPairs = N;     //number of pairs to be considered (on the average)
  // The cross section needs to be scaled since only a certain number "consideredPairs" of all pairs
  // is taken for simulation. The scale factor depends on the inital state (gg->X has a different number of pairs than gq->X etc.)
  double scaleForSelectedPairs = 0;
  const int allPairs = binomial ( N, 2 );
  
  
  if ( allPairs >= 50 )
  {   //only a certain number consideredPairs of all pairs will be considered  
    for ( int i = 0; i < consideredPairs; i++ )
    {
      m1 = int ( _allParticlesList.size() * ran2() );
      if ( m1 == _allParticlesList.size() ) 
      {
        m1 = _allParticlesList.size() - 1;
      }
      iscat = _allParticlesList[m1];

      do
      {
        m2 = int ( _allParticlesList.size() * ran2() );
        if ( m2 == _allParticlesList.size() ) 
        {
          m2 = _allParticlesList.size() - 1;
        }
        jscat = _allParticlesList[m2]; 
      }
      while ( jscat == iscat );
      

      if(theConfig->getRestrictParentPTForPhotons())
      {
        if(parentParticlesAllowed(iscat,jscat))
        {
          scaleForSelectedPairs =  static_cast<double> ( allPairs ) / static_cast<double> ( consideredPairs );
          scatt23_amongBackgroundParticles_photons_utility_1( _cells, scatt23_object, iscat, jscat, nexttime, scaleForSelectedPairs, again); 
        }else
        {
          continue;
        }
      }else
      {
        scaleForSelectedPairs =  static_cast<double> ( allPairs ) / static_cast<double> ( consideredPairs );
        scatt23_amongBackgroundParticles_photons_utility_1( _cells, scatt23_object, iscat, jscat, nexttime, scaleForSelectedPairs, again);        
      }
    }                                                    
  }
  else
  {
    for ( int i = 0; i < static_cast<int>( _allParticlesList.size() ) - 1; i++ )
    {
      iscat = _allParticlesList[i];    
      for ( unsigned int j = i+1; j < _allParticlesList.size(); j++ )
      {              
        jscat = _allParticlesList[j]; 
        scaleForSelectedPairs=1.0;        
        if(theConfig->getRestrictParentPTForPhotons())
        {
          if(parentParticlesAllowed(iscat,jscat))
          {            
            scatt23_amongBackgroundParticles_photons_utility_1( _cells, scatt23_object, iscat, jscat, nexttime, scaleForSelectedPairs, again); 
          }else
          {
            continue;
          }
        }else
        {          
          scatt23_amongBackgroundParticles_photons_utility_1( _cells, scatt23_object, iscat, jscat, nexttime, scaleForSelectedPairs, again);        
        }
      }
    }
  }
}




/**
 * Performs the inelastic photonproduction using the AMY class.
 * 
 * @param[in] _cells
 * @param[in] _allParticlesList
 * @param[in] scaleFactor
 * @param[in] again
 * @param[in] nexttime
 * 
 * 
 */
void offlineHeavyIonCollision::scatt23_amongBackgroundParticles_AMYphotons( cellContainer& _cells, std::vector< int >& _allParticlesList, const double scaleFactor, bool& again, const double nexttime )
{
  int iscat;
  FLAVOR_TYPE F1;
  double thisEOverT,temperature,TotRateGeV,probab,photonEnergyAMY,photonEnergyAMYLAB;
  VectorEPxPyPz newPhoton,newPhotonLRF;
  ParticleOffline temp_particle_produced_photon1;
  
  if(_allParticlesList.size()>2)
  {
    VectorEPxPyPz sumofMomentaIncell;
    VectorXYZ boostbetaVector;
    for ( int i = 0; i < static_cast<int>( _allParticlesList.size() ) ; i++ )
    {
      iscat = _allParticlesList[i]; 
      sumofMomentaIncell += particles_atTimeNow[iscat].Mom;
    }
    
    //HACK
    //boostbetaVector =  sumofMomentaIncell.NormalizeToE();
    boostbetaVector.SetTXYZ(0,0,0,0);
    
    
    lorentz LorentzBoost;
    LorentzBoost.setBeta( boostbetaVector );
    VectorEPxPyPz momentumLRF;
    double timestepBoosted = dt * LorentzBoost.gammaVal();
    
    for ( int i = 0; i < static_cast<int>( _allParticlesList.size() ) ; i++ )
    {
      if(std::isnan(LorentzBoost.gammaVal()) || std::isinf(LorentzBoost.gammaVal()))
      {
        continue;       
      }
      iscat = _allParticlesList[i];   
      F1 = particles_atTimeNow[iscat].FLAVOR;
      temperature = particles_atTimeNow[iscat].temperatureAMY; //GeV
      if(std::isnan(temperature) || std::isinf(temperature) || FPT_COMP_E( temperature, 0.0 ))
      {
        continue;       
      }
      
      LorentzBoost.boost(particles_atTimeNow[iscat].Mom,momentumLRF);
      thisEOverT =  momentumLRF.E()/temperature;
      if( AMY.isInMomentumRange(thisEOverT,1.))
      {
        TotRateGeV = temperature * AMY.GetTotalRateOverT(thisEOverT, ParticlePrototype::getElectricCharge(ParticlePrototype::getParticleAntiparticleType( F1 )) )/testpartcl;
        probab = timestepBoosted/0.197 * TotRateGeV *theConfig->getkFactorEMprocesses23();

        int countHowManyPhotonsDefinitively=0;
        while (FPT_COMP_G(probab,1.0))
        {
          probab-=1.0;
          countHowManyPhotonsDefinitively++;
        }
        for(int k=1;k<=countHowManyPhotonsDefinitively;k++)
        {
          photonEnergyAMY=AMY.sampleMomentumK_METRO(thisEOverT)*temperature;
          newPhotonLRF = momentumLRF;
          newPhotonLRF.NormalizeToE();
          newPhotonLRF = newPhotonLRF * photonEnergyAMY;
          LorentzBoost.boostInv(newPhotonLRF,newPhoton);
          
          totalPhotonNumber++;
          theAnalysis->PtDistributionPhotons(newPhoton.Pt(), newPhoton.Rapidity(), newPhoton.E());
          temp_particle_produced_photon1.Mom = newPhoton;
          temp_particle_produced_photon1.m = 0.0;
          temp_particle_produced_photon1.initially_produced = false;
          temp_particle_produced_photon1.FLAVOR = photon;
          temp_particle_produced_photon1.production_time = nexttime;
          noninteractingParticles.push_back(temp_particle_produced_photon1);
        
        }
        //get probabilistic photon if probab < 1
        if (FPT_COMP_L(probab,1.0))
        {
          if(ran2()<probab)
          {
            photonEnergyAMY=AMY.sampleMomentumK_METRO(thisEOverT)*temperature;
            newPhotonLRF = momentumLRF;
            newPhotonLRF.NormalizeToE();
            newPhotonLRF = newPhotonLRF * photonEnergyAMY;
            LorentzBoost.boostInv(newPhotonLRF,newPhoton);

            totalPhotonNumber++;
            theAnalysis->PtDistributionPhotons(newPhoton.Pt(), newPhoton.Rapidity(), newPhoton.E());
            temp_particle_produced_photon1.Mom = newPhoton;
            temp_particle_produced_photon1.m = 0.0;
            temp_particle_produced_photon1.initially_produced = false;
            temp_particle_produced_photon1.FLAVOR = photon;
            temp_particle_produced_photon1.production_time = nexttime;
            noninteractingParticles.push_back(temp_particle_produced_photon1);

          }
        }
      }
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

          scatt22_object.setParameter( particles_atTimeNow[iscat].Mom, particles_atTimeNow[jscat].Mom,
                                      F1, F2, M1, M2, s, md2g_wo_as , md2q_wo_as,
                                      theConfig->getKggQQb(), theConfig->getKgQgQ(), theConfig->getKappa_gQgQ(), theConfig->isConstantCrossSecGQ(),
                                      theConfig->getConstantCrossSecValueGQ(), theConfig->isIsotropicCrossSecGQ(), theConfig->getKfactor_light(), theConfig->getkFactorEMProcesses22(),
                                      theConfig->getInfraredCutOffEMProcesses(), theConfig->getKappa22Photons(),theConfig->getDebyeModePhotons(),theConfig->getVertexModePhotons(),
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

  P1new(0) = sqrt( P1new.vec2() );
  P2new(0) = sqrt( P2new.vec2() );
  P3new(0) = sqrt( P3new.vec2() );

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

      addedParticles[newIndex].Propagate( nexttime );
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

bool offlineHeavyIonCollision::parentParticlesAllowed(const int iscat,const int jscat)
{
  double maxPT,minAllowed,maxAllowed;
  if( theConfig->getRestrictParentPTForPhotons()==false )
  {
    std::string errMsg = "No restriction  for parent PT and still tries to restrict.";
    throw eHIC_error( errMsg );
  }
  minAllowed = theConfig->getMinAllowedParentPT();
  maxAllowed = theConfig->getMaxAllowedParentPT();
  maxPT = std::max(particles_atTimeNow[iscat].Mom.Pt(),particles_atTimeNow[jscat].Mom.Pt());
  if(maxPT > minAllowed && maxPT < maxAllowed)
  {
    return true;
  }else
  {
    return false;
  }
}



void offlineHeavyIonCollision::scatt22_amongBackgroundParticles_photons_utility_1( scattering22& scatt22_obj,  const int iscat, const int jscat, const double nexttime, double scaleForSelectedPairs, bool & again )
{
  FLAVOR_TYPE F1, F2;
  double Tmax, TT;
  double M1, M2;
  double t_hat, s_cutoff_for_pqcd,s,xt ;
  double Vrel,md2g_wo_as,md2q_wo_as,cs22,probab22,md2_wo_as_gluon_use,md2_wo_as_quark_use;
  int ringIndex;
  double temperature = 0.0;
  
  ParticleOffline temp_particle_iscat = particles_atTimeNow[iscat];
  ParticleOffline temp_particle_jscat = particles_atTimeNow[jscat];
  

  F1 = particles_atTimeNow[iscat].FLAVOR;
  M1 = particles_atTimeNow[iscat].m;  
  
  F2 = particles_atTimeNow[jscat].FLAVOR;
  M2 = particles_atTimeNow[jscat].m;

  unsigned int _F1 = std::min( static_cast<unsigned int>( F1), static_cast<unsigned int>( F2 ) );
  unsigned int _F2 = std::max( static_cast<unsigned int>( F1), static_cast<unsigned int>( F2 ) );
        
  s = (particles_atTimeNow[iscat].Mom + particles_atTimeNow[jscat].Mom).M2();
  if ( theConfig->getCrossSectionMethod() == 0 )
  {
    if ( ( theConfig->doScattering_22_photons() ) && ( FPT_COMP_NZ ( theConfig->getInfraredCutOffEMProcesses() ) ) )
    {
      s_cutoff_for_pqcd = 2.0 * theConfig->getInfraredCutOffEMProcesses();
      //cout << s_cutoff_for_pqcd << endl;
    }
    else   // Debye-screened 22 photonproduction
    {
      //Same lower cutoff as used for tables.
      //WARNING
      s_cutoff_for_pqcd =  0.0001;//0.01;//( 1.1 * ns_casc::lambda2 ) ;
      //For very low pT effects, this cutoff must be set to a small number. This number however, must be inline with the value of the matrix element at the smallest
      //s which could happen. In the routine get_mandelstam_t() the matrix element is checked to be greater than 10e-15. For values of s smaller than 
      //10e-5 it is about 10e-15. So 10e-4 as lower s cutoff should be fine.
    }

  }
  else // const cross section or some such.
  {
    s_cutoff_for_pqcd = 0.0;
  } 
  if( FPT_COMP_G(s, s_cutoff_for_pqcd) )// if s > s_cutoff_for_pqcd...
  {
    Vrel = VelRel( particles_atTimeNow[iscat].Mom, particles_atTimeNow[jscat].Mom, M1, M2 );

    //factor as (alpha_s) is not included in definitions of partcl[jscat].md2g, it will be multiplied in the scattering routines
    md2g_wo_as = ( particles_atTimeNow[iscat].md2g + particles_atTimeNow[jscat].md2g ) / 2.0;
    md2q_wo_as = ( particles_atTimeNow[iscat].md2q + particles_atTimeNow[jscat].md2q ) / 2.0;
    
    int initialStateIndex = -1;
    
    //Compute LRF temperature, encoded in md2g
    xt = particles_atTimeNow[iscat].Pos.Perp();
    ringIndex = rings.getIndex( xt );
    double effectiveTemperatureFromRings = rings[ringIndex].getEffectiveTemperature();
    temperature = effectiveTemperatureFromRings;//Doesn't do anything here.
    
    double Nf=3.0;
    double LRF_md2g_wo_as = ( 8.0/M_PI*pow(effectiveTemperatureFromRings,2.0)*(Ncolor+Nf) );
    double LRF_md2q_wo_as = 1.0/9.0 * LRF_md2g_wo_as ;
    
    //WARNING: Decide, which Debye mass should be used. Either LRF_md2g_wo_as or md2g_wo_as
    if (theConfig->getDebyeModePhotons()==HTLDebyepQCDrunningCouplingScale2piT)
    {
      md2_wo_as_gluon_use  = md2g_wo_as;
      md2_wo_as_quark_use = md2q_wo_as;           
    }else if (theConfig->getDebyeModePhotons()==LatticeDebye || theConfig->getDebyeModePhotons()==HTLDebyepQCDrunningCouplingScaleMD2)  
    {
      md2_wo_as_gluon_use  = LRF_md2g_wo_as;
      md2_wo_as_quark_use  = LRF_md2q_wo_as;           
    }else 
    {
      md2_wo_as_gluon_use  = md2g_wo_as;
      md2_wo_as_quark_use = md2q_wo_as;    
    }
    
    
    
    //if (temperature > 1.2)cout << temperature << "\t"<<LRF_md2q_wo_as << endl;
    scatt22_obj.setParameter( particles_atTimeNow[iscat].Mom, particles_atTimeNow[jscat].Mom,
                                F1, F2, M1, M2, s, md2_wo_as_gluon_use , md2_wo_as_quark_use,
                                theConfig->getKggQQb(), theConfig->getKgQgQ(), theConfig->getKappa_gQgQ(), theConfig->isConstantCrossSecGQ(),
                                theConfig->getConstantCrossSecValueGQ(), theConfig->isIsotropicCrossSecGQ(), theConfig->getKfactor_light(), theConfig->getkFactorEMProcesses22(),
                                theConfig->getInfraredCutOffEMProcesses(),theConfig->getKappa22Photons(),theConfig->getDebyeModePhotons(),theConfig->getVertexModePhotons(),
                                temperature, theConfig->getTdJpsi(), theConfig->isConstantCrossSecJpsi(), theConfig->getConstantCrossSecValueJpsi() ); // md2g, md2q are debye masses without the factor alpha_s which is multiplied in scattering22.cpp
  
    switch ( theConfig->getCrossSectionMethod() )
    {
      case 0:  // PQCD cross sections
        cs22 = scatt22_obj.getXSection22onlyPhotons ( initialStateIndex );
        break;
      case 1:  // const cross section
        cs22 = theConfig->getInputCrossSectionValue() * 0.1 / ( pow ( 0.197, 2.0 ) );
        break;
      default
          :
        string errMsg = "Error in method of cross section";
        throw eHIC_error ( errMsg );
    }
    // 1/ GeV^2
    //100 mb = 10 fm^2 = 100 * 2.58 1/GeV^2
    
    probab22 = pow( 0.197, 2.0 ) * cs22 * Vrel * dt * scaleForSelectedPairs  / ( dv * testpartcl );
     
    /*TEST:
    double Avv2 =  1.0/2.0*(( pow( particles_atTimeNow[iscat].Mom.Px(), 2.0 ) - pow( particles_atTimeNow[iscat].Mom.Py(), 2.0 ) ) / pow( particles_atTimeNow[iscat].Mom.Pt(), 2.0 )+( pow( particles_atTimeNow[jscat].Mom.Px(), 2.0 ) - pow( particles_atTimeNow[jscat].Mom.Py(), 2.0 ) ) / pow( particles_atTimeNow[jscat].Mom.Pt(), 2.0 ) );
    if( FPT_COMP_GZ(Avv2) )
    {
      
      //cout << "v2 = " << Avv2 << "\t s = " << s << endl;
      //cout << "cs = " << cs22 << " 1/GeV^2" << " = " << cs22/2.58 << " mb" << endl;
      averageS_Bigv2 += s ;
      numberA++;
    }
    if( FPT_COMP_LZ(Avv2) )
      {
      //cout << "v2 = " << Avv2 << "\t s = " << s << endl;
      //cout << "cs = " << cs22 << " 1/GeV^2" << " = " << cs22/2.58 << " mb" << endl;
      averageS_smallv2 += s ;
      numberB++;
    }  
    cout << "Big: " << averageS_Bigv2/numberA << "\t Small: " << averageS_smallv2/numberB << endl;
    */
    
    //cout << "Probab = " << probab22 << "\tdv = " << dv << "\ttestpartcl = " << testpartcl << "\tdt = " << dt << "\tVrel= " << Vrel << "\tNaddedEvents" << theConfig->getNaddedEvents() << endl;
    
  

    if ( FPT_COMP_G(probab22, 2.0) )
    {
      cout << "P22 photons=" << probab22 << ">1" << endl;
      again = true;
      dt = 0.5 * dt;
    }else if ( FPT_COMP_G(probab22, 1.0) )
    {
      cout << "P22 photons=" << probab22 << ">1" << endl;
      again = true;
      cout << "dt (old) = " << dt << endl;
      dt = 0.5 / ( probab22 / dt );
      cout << "dt (new) = " << dt << endl;
      return;
    }
  
    if ( ran2() < probab22 )
    {
      //TEST:
      /*if( FPT_COMP_GZ(1.0/2.0*(( pow( particles_atTimeNow[iscat].Mom.Px(), 2.0 ) - pow( particles_atTimeNow[iscat].Mom.Py(), 2.0 ) ) / pow( particles_atTimeNow[iscat].Mom.Pt(), 2.0 )+(( pow( particles_atTimeNow[jscat].Mom.Px(), 2.0 ) - pow( particles_atTimeNow[jscat].Mom.Py(), 2.0 ) ) / pow( particles_atTimeNow[jscat].Mom.Pt(), 2.0 ) ))) )
      {
        theConfig->countPositiveV2++;
      }else
      {
        theConfig->countNegativeV2++;
      }*/

      scatt22_amongBackgroundParticles_photons_utility_2( scatt22_obj, iscat, jscat, nexttime );
    }  
  }else
  {
    probab22 = -1.0;
    cs22 = 0.0; //1/GeV^2 
  }
}

void offlineHeavyIonCollision::scatt22_amongBackgroundParticles_photons_utility_2( scattering22& scatt22_obj,  const int iscat, const int jscat, const double nexttime)
{
  FLAVOR_TYPE F1, F2;
  double Tmax, TT;
  double M1, M2;
  double t_hat;
  ELASTIC_MODE_PHOTONS mechanism;
  
  ParticleOffline temp_particle_iscat = particles_atTimeNow[iscat];
  ParticleOffline temp_particle_jscat = particles_atTimeNow[jscat];
  

  F1 = particles_atTimeNow[iscat].FLAVOR;
  M1 = particles_atTimeNow[iscat].m;  
  
  F2 = particles_atTimeNow[jscat].FLAVOR;
  M2 = particles_atTimeNow[jscat].m;

  unsigned int F1_ = std::min(static_cast<unsigned int>(F1), static_cast<unsigned int>(F2));
  unsigned int F2_ = std::max(static_cast<unsigned int>(F1), static_cast<unsigned int>(F2));
  
  if (FPT_COMP_Z(F1_ * F2_) && (F1 != F2_))
  {
    //Only Compton Process here:
    mechanism = OnlyCompton;
  }
  if ((F2_ - F1_) == 1 && (F2_ % 2) == 0)
  {
    //Only Annihilation Process here:
    mechanism = OnlyAnnihilation;
  }  
  
  
  
  Tmax = std::max( particles_atTimeNow[iscat].Pos.T(), particles_atTimeNow[jscat].Pos.T() );
  TT = ( nexttime - Tmax ) * ran2() + Tmax;
  temp_particle_iscat.Propagate( TT, particles_atTimeNow[iscat].X_traveled );
  temp_particle_jscat.Propagate( TT, particles_atTimeNow[jscat].X_traveled );

  // at this point scattering take place, therefore set properties of last scattering point
  temp_particle_iscat.lastInt = particles_atTimeNow[iscat].Pos;
  temp_particle_jscat.lastInt = particles_atTimeNow[jscat].Pos;

  // determine type of scattering, momentum transfer t_hat, new flavor and new masses
  int typ;
  
  scatt22_obj.getMomentaAndMasses22onlyPhotons( F1, F2, M1, M2, t_hat, typ ); 

  // translate momentum transfer t_hat into actual momenta of outgoing particles
  VectorEPxPyPz P1new = temp_particle_iscat.Mom; //will be overwritten anyway
  VectorEPxPyPz P2new = temp_particle_jscat.Mom; //will be overwritten anyway

  VectorTXYZ Pos1new = temp_particle_iscat.Pos;
  VectorTXYZ Pos2new = temp_particle_jscat.Pos;
  
  //HACK 
  //t_hat=0.;
  
  
  
  
  scatt22_obj.setNewMomenta22( P1new, P2new, Pos1new, Pos2new, t_hat );
  
  //TEST:
  /*
  if((particles_atTimeNow[iscat].Mom.Pt() > particles_atTimeNow[jscat].Mom.Pt())&& (F1 != gluon))
  {
    P1new = particles_atTimeNow[iscat].Mom;
    F1 = photon;
    F2 = up;
  }else if ((particles_atTimeNow[jscat].Mom.Pt() > particles_atTimeNow[iscat].Mom.Pt())&& (F2 != gluon))
  {
    P2new = particles_atTimeNow[jscat].Mom;
    F2 = photon;
    F1 = up; 
  }else
  {
    F1 = up;
    F2 = up;
  }*/
  //v2nachher = 1.0/2.0 * (       ( pow(P1new.Px(), 2.0 ) - pow( P1new.Py(), 2.0 ) ) / pow( P1new.Pt(), 2.0 )     +   ( pow(P2new.Px(), 2.0 ) - pow( P2new.Py(), 2.0 ) ) / pow( P2new.Pt(), 2.0 )      );
  
  /*
  cout << " After: i  Px:\t" << P1new.Px() << "\t Py:\t" << P1new.Py() << "\t Pz:\t" << P1new.Pz() << endl;
  cout << " After: j  Px:\t" << P2new.Px() << "\t Py:\t" << P2new.Py() << "\t Pz:\t" << P2new.Pz() << endl;
  
  cout << "v2 after: " <<  1.0/2.0 * (       ( pow(P1new.Px(), 2.0 ) - pow( P1new.Py(), 2.0 ) ) / pow( P1new.Pt(), 2.0 )     +   ( pow(P2new.Px(), 2.0 ) - pow( P2new.Py(), 2.0 ) ) / pow( P2new.Pt(), 2.0 )      ) << endl;
  */
  
  /*if(FPT_COMP_L(v2vorher,v2nachher))
  {
    theConfig->v2_bigger++;
  }else
  {
    theConfig->v2_smaller++;
  }*/

  //Set the collision tags
  /*if( ((F1 != 0) && (F2 != 0)))
  {
    particles_atTimeNow[jscat].collision_tag = true;
    particles_atTimeNow[iscat].collision_tag = true;
  }*/
 
  
  //TEST no scattering, just magic photon production! Only quarks
  /*if(F1!=gluon)
  {
    P1new = temp_particle_iscat.Mom;
    Pos1new = temp_particle_iscat.Pos;  
    F1 = photon;
    particles_atTimeNow[iscat].collision_tag = true;
  }
  if(F2!=gluon)
  {
    P2new = temp_particle_jscat.Mom;
    Pos2new = temp_particle_jscat.Pos;
    F2 = photon;
    particles_atTimeNow[jscat].collision_tag = true;
  }*/
  //cout << "Control: " << endl;
  //cout << temp_particle_jscat.Mom.Px() << '\t' << temp_particle_iscat.Mom.Px() << '\t' << P1new.Px() << '\t' << P2new.Px() << endl; 
      
  if(F1==photon )
  { 
    //cout << "Photon with Energy E=" << P2new.E() << " GeV produced!" << "PT= " << P2new.Pt() << " - Perp= " << P2new.Perp()<< endl;
    totalPhotonNumber++;
    theAnalysis->PtDistributionPhotons(P1new.Pt(), P1new.Rapidity(), P1new.E());  
    ParticleOffline temp_particle_produced_photon1;
    temp_particle_produced_photon1.FLAVOR = photon;
    temp_particle_produced_photon1.m = 0.0;
    temp_particle_produced_photon1.Mom = P1new;
    temp_particle_produced_photon1.Pos = Pos1new;
    temp_particle_produced_photon1.initially_produced = false;   
    temp_particle_produced_photon1.production_time = nexttime;
    temp_particle_produced_photon1.production_mechanism = mechanism;
    noninteractingParticles.push_back(temp_particle_produced_photon1);
  }   
  if(F2==photon )
  {
    //cout << "Photon with Energy E=" << P2new.E() << " GeV produced!" << "PT= " << P2new.Pt() << " - Perp= " << P2new.Perp()<< endl;
    totalPhotonNumber++;
    theAnalysis->PtDistributionPhotons(P2new.Pt(), P2new.Rapidity(), P2new.E());     
    ParticleOffline temp_particle_produced_photon2;
    temp_particle_produced_photon2.FLAVOR = photon;
    temp_particle_produced_photon2.m = 0.0;
    temp_particle_produced_photon2.Mom = P2new;
    temp_particle_produced_photon2.Pos = Pos2new;
    temp_particle_produced_photon2.initially_produced = false;    
    temp_particle_produced_photon2.production_time = nexttime;
    temp_particle_produced_photon2.production_mechanism = mechanism;
    noninteractingParticles.push_back(temp_particle_produced_photon2);
  }
}

void offlineHeavyIonCollision::scatt22_amongBackgroundParticles_dileptons_utility_1( scattering22& scatt22_obj,  const int iscat, const int jscat, const double nexttime, double scaleForSelectedPairs, bool & again )
{
  FLAVOR_TYPE F1, F2;
  double Tmax, TT;
  double M1, M2;
  double t_hat, s_cutoff_for_pqcd,s,xt ;
  double Vrel,md2g_wo_as,md2q_wo_as,cs22,probab22,md2_wo_as_gluon_use,md2_wo_as_quark_use;
  int ringIndex;
  double temperature = 0.0;
  
  ParticleOffline temp_particle_iscat = particles_atTimeNow[iscat];
  ParticleOffline temp_particle_jscat = particles_atTimeNow[jscat];
  

  F1 = particles_atTimeNow[iscat].FLAVOR;
  M1 = particles_atTimeNow[iscat].m;  
  
  F2 = particles_atTimeNow[jscat].FLAVOR;
  M2 = particles_atTimeNow[jscat].m;

  unsigned int _F1 = std::min( static_cast<unsigned int>( F1), static_cast<unsigned int>( F2 ) );
  unsigned int _F2 = std::max( static_cast<unsigned int>( F1), static_cast<unsigned int>( F2 ) );
        
  s = (particles_atTimeNow[iscat].Mom + particles_atTimeNow[jscat].Mom).M2();

  s_cutoff_for_pqcd = 0.0;
   
  if( FPT_COMP_G(s, s_cutoff_for_pqcd) )// if s > s_cutoff_for_pqcd...
  {
    Vrel = VelRel( particles_atTimeNow[iscat].Mom, particles_atTimeNow[jscat].Mom, M1, M2 );

    //factor as (alpha_s) is not included in definitions of partcl[jscat].md2g, it will be multiplied in the scattering routines
    md2g_wo_as = ( particles_atTimeNow[iscat].md2g + particles_atTimeNow[jscat].md2g ) / 2.0;
    md2q_wo_as = ( particles_atTimeNow[iscat].md2q + particles_atTimeNow[jscat].md2q ) / 2.0;
    
    int initialStateIndex = -1;
    
    //Compute LRF temperature, encoded in md2g
    xt = particles_atTimeNow[iscat].Pos.Perp();
    ringIndex = rings.getIndex( xt );
    double effectiveTemperatureFromRings = rings[ringIndex].getEffectiveTemperature();
    temperature = effectiveTemperatureFromRings;//Doesn't do anything here.
    
    double Nf=3.0;
    double LRF_md2g_wo_as = ( 8.0/M_PI*pow(effectiveTemperatureFromRings,2.0)*(Ncolor+Nf) );
    double LRF_md2q_wo_as = 1.0/9.0 * LRF_md2g_wo_as ;
    
    //WARNING: Decide, which Debye mass should be used. Either LRF_md2g_wo_as or md2g_wo_as
    if (theConfig->getDebyeModePhotons()==HTLDebyepQCDrunningCouplingScale2piT)
    {
      md2_wo_as_gluon_use  = md2g_wo_as;
      md2_wo_as_quark_use = md2q_wo_as;           
    }else if (theConfig->getDebyeModePhotons()==LatticeDebye || theConfig->getDebyeModePhotons()==HTLDebyepQCDrunningCouplingScaleMD2)  
    {
      md2_wo_as_gluon_use  = LRF_md2g_wo_as;
      md2_wo_as_quark_use  = LRF_md2q_wo_as;           
    }else 
    {
      md2_wo_as_gluon_use  = md2g_wo_as;
      md2_wo_as_quark_use = md2q_wo_as;    
    }
    

    scatt22_obj.setParameter( particles_atTimeNow[iscat].Mom, particles_atTimeNow[jscat].Mom,
                                F1, F2, M1, M2, s, md2_wo_as_gluon_use , md2_wo_as_quark_use,
                                theConfig->getKggQQb(), theConfig->getKgQgQ(), theConfig->getKappa_gQgQ(), theConfig->isConstantCrossSecGQ(),
                                theConfig->getConstantCrossSecValueGQ(), theConfig->isIsotropicCrossSecGQ(), theConfig->getKfactor_light(), theConfig->getkFactorEMProcesses22(),
                                theConfig->getInfraredCutOffEMProcesses(),theConfig->getKappa22Photons(),theConfig->getDebyeModePhotons(),theConfig->getVertexModePhotons(),
                                temperature, theConfig->getTdJpsi(), theConfig->isConstantCrossSecJpsi(), theConfig->getConstantCrossSecValueJpsi() ); // md2g, md2q are debye masses without the factor alpha_s which is multiplied in scattering22.cpp
  
    switch ( theConfig->getCrossSectionMethod() )
    {
      case 0:  // PQCD cross sections
        cs22 = scatt22_obj.getXSection22onlyDileptons ( initialStateIndex );
        break;
      case 1:  // const cross section
        cs22 = theConfig->getInputCrossSectionValue() * 0.1 / ( pow ( 0.197, 2.0 ) );
        break;
      default
          :
        string errMsg = "Error in method of cross section";
        throw eHIC_error ( errMsg );
    }
    // 1/ GeV^2
    //100 mb = 10 fm^2 = 100 * 2.58 1/GeV^2
    
    probab22 = pow( 0.197, 2.0 ) * cs22 * Vrel * dt * scaleForSelectedPairs  / ( dv * testpartcl );
     

    if ( FPT_COMP_G(probab22, 2.0) )
    {
      cout << "P23 dileptons=" << probab22 << ">1" << endl;
      again = true;
      dt = 0.5 * dt;
    }else if ( FPT_COMP_G(probab22, 1.0) )
    {
      cout << "P23 dileptons=" << probab22 << ">1" << endl;
      again = true;
      cout << "dt (old) = " << dt << endl;
      dt = 0.5 / ( probab22 / dt );
      cout << "dt (new) = " << dt << endl;
      return;
    }
  
    if ( ran2() < probab22 )
    {
      scatt22_amongBackgroundParticles_dileptons_utility_2( scatt22_obj, iscat, jscat, nexttime );
    }  
  }else
  {
    probab22 = -1.0;
    cs22 = 0.0; //1/GeV^2 
  }
}

void offlineHeavyIonCollision::scatt22_amongBackgroundParticles_dileptons_utility_2( scattering22& scatt22_obj,  const int iscat, const int jscat, const double nexttime)
{
  FLAVOR_TYPE F1, F2;
  double Tmax, TT;
  double M1, M2;
  double t_hat;
  ELASTIC_MODE_PHOTONS mechanism;
  
  ParticleOffline temp_particle_iscat = particles_atTimeNow[iscat];
  ParticleOffline temp_particle_jscat = particles_atTimeNow[jscat];
  

  F1 = particles_atTimeNow[iscat].FLAVOR;
  M1 = particles_atTimeNow[iscat].m;  
  
  F2 = particles_atTimeNow[jscat].FLAVOR;
  M2 = particles_atTimeNow[jscat].m;

  unsigned int F1_ = std::min(static_cast<unsigned int>(F1), static_cast<unsigned int>(F2));
  unsigned int F2_ = std::max(static_cast<unsigned int>(F1), static_cast<unsigned int>(F2));
  
  
  Tmax = std::max( particles_atTimeNow[iscat].Pos.T(), particles_atTimeNow[jscat].Pos.T() );
  TT = ( nexttime - Tmax ) * ran2() + Tmax;
  temp_particle_iscat.Propagate( TT, particles_atTimeNow[iscat].X_traveled );
  temp_particle_jscat.Propagate( TT, particles_atTimeNow[jscat].X_traveled );

  // at this point scattering take place, therefore set properties of last scattering point
  temp_particle_iscat.lastInt = particles_atTimeNow[iscat].Pos;
  temp_particle_jscat.lastInt = particles_atTimeNow[jscat].Pos;

  // determine type of scattering, momentum transfer t_hat, new flavor and new masses
  int typ;
  switch ( theConfig->getCrossSectionMethod() )
  {
    case 0:  // PQCD cross sections
      scatt22_obj.getMomentaAndMasses22onlyDileptons ( F1, F2, M1, M2, t_hat, typ );
      break;
    default
        :
      string errMsg = "Error in method of cross section for dileptons.";
      throw eHIC_error ( errMsg );
  }
  
  if ( F1 == lepton && F2 == lepton  )
  {
    // translate momentum transfer t_hat into actual momenta of outgoing particles
    VectorEPxPyPz P1new = temp_particle_iscat.Mom; //will be overwritten anyway
    VectorEPxPyPz P2new = temp_particle_jscat.Mom; //will be overwritten anyway
    VectorTXYZ Pos1new = temp_particle_iscat.Pos;
    VectorTXYZ Pos2new = temp_particle_jscat.Pos;
    
    scatt22_obj.setNewMomenta22( P1new, P2new, Pos1new, Pos2new, t_hat );
    
    //Get invariant Mass and Momentum
    ParticleOffline temp_particle_produced_dilepton_pair;
    temp_particle_produced_dilepton_pair.FLAVOR = dilepton;
    temp_particle_produced_dilepton_pair.m = sqrt(  ( P1new + P2new ).M2()  );
    temp_particle_produced_dilepton_pair.Mom = P1new + P2new;
    temp_particle_produced_dilepton_pair.dead = false;
    temp_particle_produced_dilepton_pair.cell_id = -1;
    temp_particle_produced_dilepton_pair.unique_id = Particle::unique_id_counter;
    temp_particle_produced_dilepton_pair.production_time = nexttime;
    temp_particle_produced_dilepton_pair.dilepton_pt_min = std::min(P1new.Pt(),P2new.Pt());
    temp_particle_produced_dilepton_pair.dilepton_pt_max = std::max(P1new.Pt(),P2new.Pt());
    temp_particle_produced_dilepton_pair.dilepton_y_min  = std::min(P1new.Rapidity(),P2new.Rapidity());
    temp_particle_produced_dilepton_pair.dilepton_y_max  = std::max(P1new.Rapidity(),P2new.Rapidity());    
    ++Particle::unique_id_counter;
    dileptons.push_back ( temp_particle_produced_dilepton_pair );
    
    totalDileptonNumber++;
  }   
}

int offlineHeavyIonCollision::getSpecificScatteringType(int _F1, int _F2)
{
  if( _F1 > _F2 )
  {
    cout << "F1: " << _F1 << "   F2: " << _F2 << endl;
    std::string errMsg = "Please sort the flavors in box::getSpecificScatteringType. _F1 must be smaller or equal than _F2";
    throw eHIC_error( errMsg );
  }
  int scatteringType =0;
  
  if(_F1 + _F2 == 0)
  {
    scatteringType = -1;
  }else if(_F1*_F2 == 0)
  {
    scatteringType = -1;
  }else if ( _F1 == _F2 )  // qq -> qq, qbarqbar -> qbarqbar  (only for light quarks)
  {
    scatteringType = 0;
  }else if( (_F2 - _F1) == 1 &&  (_F2 % 2) == 0 ) //qqbar->qqbar
  {
    scatteringType = 1;
  }else if( _F2 != _F1 )  //qqbar' -> qqbar'= qq' -> qq'
  {
    scatteringType = 2;
  }else
  {
    scatteringType = -1;
  }
  
  return scatteringType;
}

double offlineHeavyIonCollision::getLambdaFromRates(int _F1,int _F2, ratesManager& rates)
{
  //0.0001; 0.00001;
  double TuneMFPArbitrary = 1.0;
  int scatteringType = getSpecificScatteringType(_F1,_F2);
  double lambda;
  if(scatteringType > -1)
  {
    lambda = rates.getLambdaSpecific(scatteringType,fm)*TuneMFPArbitrary;
    //cout << lambda << endl;
    
    return lambda;
  }
  else
  {
    return 0;
  }
}

/**
 * Handles the Sampling of the inelastic photons and saves them in the photon vector.
 * @param[in] scatt23_obj
 * @param[in] iscat
 * @param[in] jscat
 * @param[in] nexttime
 * @param[in] scaleForSelectedPairs, 
 * @param[out] again
 * 
 */
void offlineHeavyIonCollision::scatt23_amongBackgroundParticles_photons_utility_1(cellContainer& _cells, scattering23& scatt23_obj, const int iscat, const int jscat, const double nexttime, double scaleForSelectedPairs, bool & again )
{
  FLAVOR_TYPE F1, F2;
  double M1,M2,s,probab22,cs22,Vrel,md2g_wo_as,md2q_wo_as,xt,lambda,lambda_scaled,cs23Photons,probab23,s_cutoff_for_pqcd,betaDistEntry,cs23Total;
  double md2_wo_as_gluon_use,md2_wo_as_quark_use;
  int ringIndex;
  VectorEPxPyPz P1, P2;
  
  F1 = particles_atTimeNow[iscat].FLAVOR;
  M1 = particles_atTimeNow[iscat].m;
  P1 = particles_atTimeNow[iscat].Mom;

  F2 = particles_atTimeNow[jscat].FLAVOR;
  M2 = particles_atTimeNow[jscat].m;
  P2 = particles_atTimeNow[jscat].Mom;
  
  unsigned int _F1 = std::min( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );
  unsigned int _F2 = std::max( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );
  

  //This s-cutoff can be adjusted!
  s_cutoff_for_pqcd = 1.1*lambda2; //1.1*LambdaQCD^2=1.1*0.2*0.2=0.044
  s = (particles_atTimeNow[iscat].Mom + particles_atTimeNow[jscat].Mom).M2();
  
  if( FPT_COMP_G(s, s_cutoff_for_pqcd) )// if s > s_cutoff_for_pqcd...
  {

    Vrel = VelRel( particles_atTimeNow[iscat].Mom, particles_atTimeNow[jscat].Mom, M1, M2 );

    //factor as (alpha_s) is not included in definitions of partcl[jscat].md2g, it will be multiplied in the scattering routines
    md2g_wo_as = ( particles_atTimeNow[iscat].md2g + particles_atTimeNow[jscat].md2g ) / 2.0;
    md2q_wo_as = ( particles_atTimeNow[iscat].md2q + particles_atTimeNow[jscat].md2q ) / 2.0;
    
    int initialStateIndex = -1;
    
    //Compute MFP
    lambda = getLambdaFromRates( _F1,_F2, _cells.rates);      
    lambda_scaled = lambda * sqrt( s ) / 0.197;// lambda in fm, sqrt(s) in GeV, lambda_scaled dimensionless

    xt = particles_atTimeNow[iscat].Pos.Perp();
    ringIndex = rings.getIndex( xt );
    double effectiveTemperatureFromRings = rings[ringIndex].getEffectiveTemperature();
    
    if(effectiveTemperatureFromRings>2.5)
    {
      effectiveTemperatureFromRings = 2.5;
    }
    //[lambda]=fm
    /*
    lambda = theMFP.getMeanFreePath( particles_atTimeNow[iscat].Mom.E(), particles_atTimeNow[iscat].FLAVOR, effectiveTemperatureFromRings, rings[ringIndex].getGluonDensity(), rings[ringIndex].getQuarkDensity(), fm );
    lambda = _cells.SpecificMFP[0];
    */
    
    //TEST:       
    //lambda = 10.5;//fm
    //cout << "MFP [fm] = " << lambda << endl;
    //lambda_scaled = lambda * sqrt( s ) / 0.197; // lambda in fm, sqrt(s) in GeV, lambda_scaled dimensionless
          
    //if(lambda_scaled>100.0)
    //{
    //  cout << sqrt( s ) << endl;
    //  cout << "Big lambda: " << lambda_scaled << "\t" << lambda << endl;
    //}
    /*if(lambda_scaled>1.0)
    {
      cout << lambda_scaled << endl;
    }*/  
    
          
//     if(lambda>0 && rings[ringIndex].getEffectiveTemperature()<0.41 && rings[ringIndex].getEffectiveTemperature()>0.39)
//     {
//       lambdaT400+=lambda;
//       countLambdaatT400+=1;
//     }
    
    /*TEST
    cout << "1/ns = " << 1.0/(rings[ringIndex].getGluonDensity()*1.0) << endl;
    cout << "E = " <<  particles_atTimeNow[iscat].Mom.E() << endl;
    cout << "T= : " << rings[ringIndex].getEffectiveTemperature() << endl;
    cout << rings[ringIndex].getGluonDensity() << endl;
    cout << rings[ringIndex].getQuarkDensity() << endl;
    cout <<"sqrt(s)="<< sqrt(s)<< " - New LPM ktCO: " << theConfig->get23FudgeFactorLpm()/lambda_scaled << "\t old(400MeV,X=0.3): " << theConfig->get23ktCutOffPhotons()/sqrt(s) << endl;
    */
    
    //TEST
    //cout << "theConfig->get23ktCutOffPhotons() " << theConfig->get23ktCutOffPhotons() << endl;
    //cout << "theConfig->getLPMModePhotons()    " << theConfig->getLPMModePhotons()    << endl;
    
    
    
    if( lambda > 0 )
    {          
      double Nf=3.0;
      double scaled_LRF_md2g_wo_as = ( 8.0/M_PI*pow(effectiveTemperatureFromRings,2.0)*(Ncolor+Nf) ) / s;
      double scaled_LRF_md2q_wo_as = 1.0/9.0 * scaled_LRF_md2g_wo_as ;
      
      //WARNING: Decide, which Debye mass should be used. Either scaled_LRF_md2g_wo_as or md2g_wo_as/s
      
      //WARNING: Decide, which Debye mass should be used. Either LRF_md2g_wo_as or md2g_wo_as
      if (theConfig->getDebyeModePhotons()==HTLDebyepQCDrunningCouplingScale2piT)
      {
        md2_wo_as_gluon_use  = md2g_wo_as/s;
        md2_wo_as_quark_use = md2q_wo_as/s;      
      }else if (theConfig->getDebyeModePhotons()==LatticeDebye || theConfig->getDebyeModePhotons()==HTLDebyepQCDrunningCouplingScaleMD2)  
      {
        md2_wo_as_gluon_use  = scaled_LRF_md2g_wo_as;
        md2_wo_as_quark_use  = scaled_LRF_md2g_wo_as;           
      }else 
      {
        md2_wo_as_gluon_use  = md2g_wo_as/s;
        md2_wo_as_quark_use = md2q_wo_as/s;  
      }     
 

      
      
      betaDistEntry = scatt23_obj.setParameter ( VectorXYZ ( 0, 0, 0 ), P1, P2, F1, F2, M1, M2, sqrt ( s ),
                            md2_wo_as_gluon_use, lambda_scaled,
                            theConfig->getK23LightPartons(), theConfig->getK23HeavyQuarks(),
                            theConfig->getKappa23LightPartons(), theConfig->getKappa23HeavyQuarks(),
                            theConfig->I23onlineIntegrationIsSet(),
                            theConfig->get23GluonFormationTimeTyp(), theConfig->getMatrixElement23(), theConfig->isMd2CounterTermInI23(),
                            theConfig->get23FudgeFactorLpm(), -1, theConfig->isMatrixElement23_22qt(),theConfig->I23onlineIntegrationPhotonsIsSet(),theConfig->get23ktCutOffPhotons(), theConfig->getkFactorEMprocesses23(), md2_wo_as_quark_use, theConfig->getDebyeModePhotons(),theConfig->getVertexModePhotons(),theConfig->getLPMModePhotons() );
      
      cs23Photons = scatt23_obj.getTotalCrossSectionForPhotons ( initialStateIndex ); //1/GeV^2
      cs23Total =  cs23Photons;
      // 1/ GeV^2
    }else
    {
      //TODO: correct error handling for zero mfp.
      //cout << lambda << "\t" << effectiveTemperatureFromRings << endl;
      cs23Total = 0;
    }       
    probab23 = pow( 0.197, 2.0 ) * cs23Total * Vrel * dt * scaleForSelectedPairs / ( dv * testpartcl );

    
    if ( FPT_COMP_G(probab23, 3.0) )
    {
      cout << "P23 photons=" << probab23 << ">3" << endl;
      again = true;
      dt = 0.5 * dt;
    }else if ( FPT_COMP_G(probab23, 1.0) )
    {
      //cout << "P23 photons=" << probab23 << ">1" << endl;
      again = true;
      //cout << "dt (old) = " << dt << endl;
      dt = 0.5 / ( probab23 / dt );
      //cout << "dt (new) = " << dt << endl;
      return;
    }
    
    if ( ran2() < probab23 )
    {       
        scatt23_amongBackgroundParticles_photons_utility_2( scatt23_obj, iscat, jscat, nexttime );
    }  
  }
  else
  {
    probab23 = -1.0;
    cs22 = 0.0; //1/GeV^2
  }
}



void offlineHeavyIonCollision::NumbersInCell( std::vector< int >& ThisCell, double & NumberGluonsInCell, double & NumberUpsInCell, double & NumberAntiupsInCell, double & NumberDownsInCell, double & NumberAntisdownsInCell, double &NumberStrangesInCell, double & NumberAntiStrangesInCell )
{ 
  int wscat;

  NumberGluonsInCell =0 ;
  NumberUpsInCell =0 ;
  NumberAntiupsInCell =0 ;
  NumberDownsInCell =0 ;
  NumberAntisdownsInCell =0 ;
  NumberStrangesInCell =0 ;
  NumberAntiStrangesInCell =0 ;

  for ( int i = 0; i < static_cast<int>( ThisCell.size() ) - 1; i++ )
  {
    wscat =  ThisCell[i];    
    if (particles_atTimeNow[wscat].FLAVOR == gluon){NumberGluonsInCell++;}
    if (particles_atTimeNow[wscat].FLAVOR == up){NumberUpsInCell++;}
    if (particles_atTimeNow[wscat].FLAVOR == anti_up){NumberAntiupsInCell++;}
    if (particles_atTimeNow[wscat].FLAVOR == down){NumberDownsInCell++;}
    if (particles_atTimeNow[wscat].FLAVOR == anti_down){NumberAntisdownsInCell++;}
    if (particles_atTimeNow[wscat].FLAVOR == strange){NumberStrangesInCell++;}
    if (particles_atTimeNow[wscat].FLAVOR == anti_strange){NumberAntiStrangesInCell++;}
  }      
}


void offlineHeavyIonCollision::scatt22ForRates(  cellContainer& _cells,std::vector< int >& _allParticlesList,  const interpolation22& theI22, const int _NumberOfCellsAveraged)
{
  int iscat, jscat;
   
  list<int>::const_iterator iIt;
  list<int>::const_iterator jIt;
 
  scattering22 scatt22_object( &theI22 );
 
  _cells.rates.clearSpecific();
  
  for ( int i = 0; i < static_cast<int>( _allParticlesList.size() ) - 1; i++ )
  {
    iscat = _allParticlesList[i];    
    for ( unsigned int j = i+1; j < _allParticlesList.size(); j++ )
    {              
      jscat = _allParticlesList[j]; 
      if(particles_atTimeNow[iscat].dead || particles_atTimeNow[jscat].dead)
      {
        continue;
      }else
      { 
        scatt22ForRatesUtility(scatt22_object, iscat, jscat, theI22, _cells,_allParticlesList, _NumberOfCellsAveraged);
      }
    }
  }
}

void offlineHeavyIonCollision::scatt22ForRatesUtility(scattering22& scatt22_obj, const int iscat, const int jscat, const interpolation22& theI22,  cellContainer& _cells,std::vector< int >& _allParticlesList, const int NumberOfCellsAveraged)
{
  FLAVOR_TYPE F1, F2;
  double M1, M2;
  double s,s_cutoff_for_pqcd;
  double Vrel,md2g_wo_as,md2q_wo_as,cs22,probab22,md2_wo_as_gluon_use,md2_wo_as_quark_use;
  double t_hat,xt ;
  int ringIndex;
  double temperature = 0.0;
  
  ParticleOffline temp_particle_iscat = particles_atTimeNow[iscat];
  ParticleOffline temp_particle_jscat = particles_atTimeNow[jscat];

  double ng,nu,nub,nd,ndb,ns,nsb;
  NumbersInCell(_allParticlesList,ng,nu,nub,nd,ndb,ns,nsb);  
  
  double AverageQuarkNumber = (nu+nd+ns+nub+ndb+nsb)/6.;
  double SumOfQuarkNumber = nu+nd+ns+nub+ndb+nsb;
  averageQuarkNumber += SumOfQuarkNumber;
  countForAverageaverageQuarkNumber++;
  //cout << ng << "\t" << nu << "\t" << nub << "\t" << nd << "\t" << ndb << "\t" << ns << "\t" << nsb << "\t" << endl;
  
  
  F1 = particles_atTimeNow[iscat].FLAVOR;
  M1 = particles_atTimeNow[iscat].m;  
  
  F2 = particles_atTimeNow[jscat].FLAVOR;
  M2 = particles_atTimeNow[jscat].m;
  
  unsigned int _F1 = std::min( static_cast<unsigned int>( F1), static_cast<unsigned int>( F2 ) );
  unsigned int _F2 = std::max( static_cast<unsigned int>( F1), static_cast<unsigned int>( F2 ) );

  s = (particles_atTimeNow[iscat].Mom + particles_atTimeNow[jscat].Mom).M2();
  Vrel = VelRel( particles_atTimeNow[iscat].Mom, particles_atTimeNow[jscat].Mom, M1, M2 );
  s_cutoff_for_pqcd = 0.1;
  
  if ( FPT_COMP_G ( s, s_cutoff_for_pqcd ) )
  {
    //factor as (alpha_s) is not included in definitions of partcl[jscat].md2g, it will be multiplied in the scattering routines
    md2g_wo_as = ( particles_atTimeNow[iscat].md2g + particles_atTimeNow[jscat].md2g ) / 2.0;
    md2q_wo_as = ( particles_atTimeNow[iscat].md2q + particles_atTimeNow[jscat].md2q ) / 2.0;
  
    int initialStateIndex = -1;
    //Compute LRF temperature, encoded in md2g
    xt = particles_atTimeNow[iscat].Pos.Perp();
    ringIndex = rings.getIndex( xt );
    double effectiveTemperatureFromRings = rings[ringIndex].getEffectiveTemperature();
    temperature = effectiveTemperatureFromRings;//Doesn't do anything here.
    
    double Nf=3.0;
    double LRF_md2g_wo_as = ( 8.0/M_PI*pow(effectiveTemperatureFromRings,2.0)*(Ncolor+Nf) );
    double LRF_md2q_wo_as = 1.0/9.0 * LRF_md2g_wo_as ;
    
    //WARNING: Decide, which Debye mass should be used. Either LRF_md2g_wo_as or md2g_wo_as
    if (theConfig->getDebyeModePhotons()==HTLDebyepQCDrunningCouplingScale2piT)
    {
      md2_wo_as_gluon_use  = md2g_wo_as;
      md2_wo_as_quark_use = md2q_wo_as;           
    }else if (theConfig->getDebyeModePhotons()==LatticeDebye || theConfig->getDebyeModePhotons()==HTLDebyepQCDrunningCouplingScaleMD2)  
    {
      md2_wo_as_gluon_use  = LRF_md2g_wo_as;
      md2_wo_as_quark_use  = LRF_md2q_wo_as;           
    }else 
    {
      md2_wo_as_gluon_use  = md2g_wo_as;
      md2_wo_as_quark_use = md2q_wo_as;    
    }
    
    scatt22_obj.setParameter ( particles_atTimeNow[iscat].Mom, particles_atTimeNow[jscat].Mom, F1, F2, M1, M2, s,md2_wo_as_gluon_use , md2_wo_as_quark_use,
                                  theConfig->getKggQQb(), theConfig->getKgQgQ(), theConfig->getKappa_gQgQ(), theConfig->isConstantCrossSecGQ(),
                                  theConfig->getConstantCrossSecValueGQ(), theConfig->isIsotropicCrossSecGQ(), theConfig->getKfactor_light(), theConfig->getkFactorEMProcesses22(),
                                  theConfig->getInfraredCutOffEMProcesses(),theConfig->getKappa22Photons(),theConfig->getDebyeModePhotons(),theConfig->getVertexModePhotons(), 0.0, theConfig->getTdJpsi(), theConfig->isConstantCrossSecJpsi(), theConfig->getConstantCrossSecValueJpsi()
                                ); // md2g_wo_as, md2q_wo_as are debye masses without the factor alpha_s which is multiplied in scattering22.cpp

    //****************************************************************************************************

    int scatteringType = getSpecificScatteringType(_F1,_F2);
    cs22=0.0;
    
    if(scatteringType == 0) // qq -> qq, qbarqbar -> qbarqbar  (only for light quarks)
    {
      cs22 = scatt22_obj.getXSection22Specific ( 227 );
      probab22 = pow( 0.197, 2.0 ) * cs22 * Vrel * dt  / ( NumberOfCellsAveraged * dv * testpartcl );
      if(_F1 != _F2)
      {
        std::string errMsg = "Wrong scattering type in scatt22ForRatesUtility.";
        throw eHIC_error( errMsg );
      }
      double specificQuarkNumberFactor=0;
      switch ( _F1 )
      {
        case 1:  // up quark
          specificQuarkNumberFactor = nu-1;
          break;
        case 2:  // up-bar quark
          specificQuarkNumberFactor = nub-1;
          break;
        case 3:  // down quark
          specificQuarkNumberFactor = nd-1;
          break;
        case 4:  // up quark
          specificQuarkNumberFactor = ndb-1;
          break;
        case 5:  // up quark
          specificQuarkNumberFactor = ns-1;
          break;
        case 6:  // up quark
          specificQuarkNumberFactor = nsb-1;
          break;
        /* Mapped according to:
        * 0 = g (gluon)
        * 1 = u (up)
        * 2 = ub (anti-up)
        * 3 = d (down)
        * 4 = db (anti-down)
        * 5 = s (strange)
        * 6 = sb (anti-strange)
        */        
      }  
      double FlavorAverage = 1./6.;
      _cells.rates.addSpecific(scatteringType,FlavorAverage*probab22/dt*2.0/(specificQuarkNumberFactor) , pow( 0.197, 2.0 ) * cs22 * Vrel/ ( dv * testpartcl ) );  
    }
    else if(scatteringType == 1)//qqbar->qqbar
    {      
      cs22 = scatt22_obj.getXSection22Specific ( 224 );
      probab22 = pow( 0.197, 2.0 ) * cs22 * Vrel * dt   / ( dv * testpartcl );
      double specificQuarkNumberFactor=0;
      switch ( _F1 )
      {
        case 1:  // up quark
          specificQuarkNumberFactor = (nu+nub)/2.0;
          break;
        case 2:  // up-bar quark
          specificQuarkNumberFactor = (nu+nub)/2.0;
          break;
        case 3:  // down quark
          specificQuarkNumberFactor = (nd+ndb)/2.0;
          break;
        case 4:  // up quark
          specificQuarkNumberFactor = (nd+ndb)/2.0;
          break;
        case 5:  // up quark
          specificQuarkNumberFactor = (ns+nsb)/2.0;
          break;
        case 6:  // up quark
          specificQuarkNumberFactor = (ns+nsb)/2.0;
          break;
        /* Mapped according to:
        * 0 = g (gluon)
        * 1 = u (up)
        * 2 = ub (anti-up)
        * 3 = d (down)
        * 4 = db (anti-down)
        * 5 = s (strange)
        * 6 = sb (anti-strange)
        */        
      } 
      double FlavorAverage = 1./3.;
      _cells.rates.addSpecific(scatteringType,FlavorAverage*probab22/dt/(specificQuarkNumberFactor) , pow( 0.197, 2.0 ) * cs22 * Vrel/ ( dv * testpartcl ) );  
    }
    //This process would not be included in the AMY scattering
    /*else if (scatteringType == 2)//qqbar' -> qqbar'= qq' -> qq'
    {
      cs22 = scatt22_obj.getXSection22Specific ( 228 );
      probab22 = pow( 0.197, 2.0 ) * cs22 * Vrel * dt   / ( dv * testpartcl );
      double specificQuarkNumberFactor=0;
      switch ( _F1 )
      {
        case 1:  // up quark
          specificQuarkNumberFactor = (nu+nub)/2.0;
          break;
        case 2:  // up-bar quark
          specificQuarkNumberFactor = (nu+nub)/2.0;
          break;
        case 3:  // down quark
          specificQuarkNumberFactor = (nd+ndb)/2.0;
          break;
        case 4:  // up quark
          specificQuarkNumberFactor = (nd+ndb)/2.0;
          break;
        case 5:  // up quark
          specificQuarkNumberFactor = (ns+nsb)/2.0;
          break;
        case 6:  // up quark
          specificQuarkNumberFactor = (ns+nsb)/2.0;
          break;
        // Mapped according to:
        // 0 = g (gluon)
        // 1 = u (up)
        // 2 = ub (anti-up)
        // 3 = d (down)
        // 4 = db (anti-down)
        // 5 = s (strange)
        // 6 = sb (anti-strange)
        //        
      }        
      double FlavorAverage = 1./6.;
      _cells.rates.addSpecific(scatteringType,FlavorAverage*probab22/dt/(specificQuarkNumberFactor) , pow( 0.197, 2.0 ) * cs22 * Vrel/ ( dv * testpartcl ) );    

    }*/ 
  }
}


/**
 * Handles the Sampling of the inelastic photons and saves them in the photon vector.
 * @param[in] scatt23_obj
 * @param[in] iscat
 * @param[in] jscat
 * @param[in] nexttime
 * 
 * scattering23& scatt23_obj, std::list< int >& _cellMembersAdded, std::vector< int >& _allParticlesList, const int iscat, const int jscat, int& typ, const double nexttime
 */
void offlineHeavyIonCollision::scatt23_amongBackgroundParticles_photons_utility_2( scattering23& scatt23_obj, const int iscat, const int jscat, const double nexttime )
{
  FLAVOR_TYPE F1, F2;
  double Tmax, TT;
  double M1, M2;
  double t_hat;
  
  VectorEPxPyPz P1new, P2new, P3new;
  double pt1, pt3, y, phi, pz1;
  
  ParticleOffline temp_particle_iscat = particles_atTimeNow[iscat];
  ParticleOffline temp_particle_jscat = particles_atTimeNow[jscat];
  
  F1 = temp_particle_iscat.FLAVOR;
  M1 = temp_particle_iscat.m;  

  F2 = temp_particle_jscat.FLAVOR;
  M2 = temp_particle_jscat.m;
                    
  int get23errors;
  get23errors += scatt23_obj.getPhotonMomenta23_metropolis( pt1, pt3, y, phi, pz1 );
  
  //TEST
  scatt23_obj.setNewMomenta23( P1new, P2new, P3new,temp_particle_iscat.Pos, temp_particle_jscat.Pos,pt1, pt3, y, phi, pz1 );
  //scatt23_obj.setNewMomenta23onlyCMBoost( P1new, P2new, P3new,temp_particle_iscat.Pos, temp_particle_jscat.Pos,pt1, pt3, y, phi, pz1 );
  
  ParticleOffline temp_particle_produced_photon1;
  totalPhotonNumber++;
  theAnalysis->PtDistributionPhotons(P3new.Pt(), P3new.Rapidity(), P3new.E());
  temp_particle_produced_photon1.Mom = P3new;
  temp_particle_produced_photon1.m = 0.0;
  temp_particle_produced_photon1.initially_produced = false;
  temp_particle_produced_photon1.FLAVOR = photon;
  temp_particle_produced_photon1.production_time = nexttime;
  noninteractingParticles.push_back(temp_particle_produced_photon1);
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
  scattering23 scatt23_object( &theI23_massless, &theI23_charm_m1, &theI23_charm_m2, &theI23_bottom_m1, &theI23_bottom_m2, &theI23_photons );
  
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
  scattering23 scatt23_object( &theI23_massless, &theI23_charm_m1, &theI23_charm_m2, &theI23_bottom_m1, &theI23_bottom_m2, &theI23_photons );
  
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
