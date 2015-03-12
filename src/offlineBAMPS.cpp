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
#include "revision.h"


using namespace std;
using namespace ns_casc;

int Ntest;


int main( int argc, char *argv[] )
{
  cout << "# Version: SVN " << SVN_REVISION << endl;
  
  try
  {
    //--------------------------------------------------------------
    // create and initialize the main objects needed for
    // configuration,execution and analysis of the simulation
    config theConfig;
    theConfig.readAndProcessProgramOptions( argc, argv );
    
    offlineOutputInterface offlineInterface( theConfig.getPathdirOfflineDataChar() );
    offlineInterface.setAdditionalFilenameTag( theConfig.getOriginalName() );    
    
    theConfig.readAndPrepareInitialSettings( &offlineInterface );

    analysis theAnalysis( &theConfig );
    offlineHeavyIonCollision theHIC( &theConfig, &offlineInterface, &theAnalysis );
    //--------------------------------------------------------------

    //--------------------------------------------------------------
    // initialize the random number generator
    // defined globally in random.h and random.cpp

    uint32_t seed = theConfig.getSeed();
    if (seed == 0) seed = ran2.findSeed();
    ran2.setSeed( seed );
    theAnalysis.setSeed( seed );
    theConfig.setSeed( seed );
    cout << "seed: " << seed << endl;
    //--------------------------------------------------------------


    //--------cascade-------------------------------------
    cout << "=============start===============" << endl;
    
    
    if( theConfig.isOnlyMediumEvolution() )
    {
      theHIC.onlyMediumEvolution( theAnalysis );
    }
    else
    {
      theHIC.initialize();
      theHIC.mainFramework();
    }

    cout << "==============end================" << endl;
    //--------end of cascade------------------------------
  }
  /**
  *
  * handle exceptions
  *
  */
  catch (int e) // provided to handle program termination after displaying usage and help messages
  {
    cout << "Signal caught: " << e << endl;
    if ( e == HELP_MESSAGE_REQUESTED )
    {
      return EXIT_SUCCESS;
    }
    else
    {
      return EXIT_FAILURE;
    }    
  }  
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
