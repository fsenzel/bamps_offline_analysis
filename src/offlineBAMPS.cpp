//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------
//---------------------------------------------


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <stdint.h>
#include <stdexcept>
#include "offlineBAMPS.h"
#include "configuration.h"
#include "additionalparticlesdistribution.h"
#include "analysis.h"
#include "binning.h"
#include "random.h"
#include "offlineheavyioncollison.h"
#include "offlineoutput.h"


using namespace std;
using namespace ns_casc;

int Ntest;


int main( int argc, char *argv[] )
{
  try
  {
    int ng;
    double stoptime, time, nexttime;
    bool readCellStructure = false;
    bool again;
    bool dodo = false;
    int nn_ana = 0;
    
    //--------------------------------------------------------------
    // create and initialize the main objects needed for
    // configuration,execution and analysis of the simulation
    config theConfig( argc, argv );
    
    offlineOutputInterface offlineInterface( theConfig.getPathdirOfflineDataChar() );
    offlineInterface.setAdditionalFilenameTag( theConfig.getOriginalName() );    
    
    theConfig.readAndPrepareInitialSettings( &offlineInterface );

    analysis theAnalysis( &theConfig );
    offlineHeavyIonCollision theHIC( &theConfig, &offlineInterface );
    //--------------------------------------------------------------

    //--------------------------------------------------------------
    // initialize the random number generator
    // defined globally in random.h and random.cpp
    uint32_t seed;

    if ( theConfig.getSeed() != 0 )
    {
      seed = theConfig.getSeed();
      ran2.setSeed( seed );
    }
    else
    {
      seed = ran2.setSeed();
    }
    theAnalysis.setSeed( seed );
    cout << "seed: " << seed << endl;
    //--------------------------------------------------------------


    //--------cascade-------------------------------------
    cout << "=============start===============" << endl;

    theHIC.initialize();
    theHIC.mainFramework( theAnalysis );

    cout << "==============end================" << endl;
    //--------end of cascade------------------------------
  }
  /**
  *
  * handle exceptions
  *
  */
  catch ( std::exception& e )
  {
    // output of the error information provided by the exception class
    cout << e.what() << endl;
    cout << "Program terminated." << endl;
    return EXIT_FAILURE;
  }
  catch ( ... )
  {
    cout << "Unhandled exception. Program terminated." << endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;
