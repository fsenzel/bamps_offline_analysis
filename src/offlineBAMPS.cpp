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



using namespace std;
using namespace ns_casc;

int Ntest;
extern int numberAdded, number;
extern int IX, IY, IZ;


extern double dt_cascade;

// ofstream printdens0("output/density.dat");
// ofstream printcol("output/collisions.dat");

// double md2g_min = 0.05 / gluon_alpha_s(0.1); // md2g = 0.05 corresponds to T = 100 MeV, / gluon_alpha_s(s) because it is compared to md2g/as
// double md2g_min = 0.0001 / gluon_alpha_s( 0.1 );  // md2g is tabled for cross section above this value, cannot deal with smaller one -> however, checked that a smaller value actually does not occur
// binning md2bins( 0.0, md2g_min*0.3, 200 );
// binning sbins( 0.0, 4.0, 500 );

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

    double factor_dt = 0.8;
    if ( argc == 3 )
    {
      factor_dt = atof( argv[2] );
    }
    cout << "factor for dt = " << factor_dt << endl;


    //--------------------------------------------------------------
    // create and initialize the main objects needed for
    // configuration,execution and analysis of the simulation
    config theConfig( argc, argv );
    theConfig.setting();

    analysis theAnalysis( &theConfig );
    offlineHeavyIonCollision theHIC( &theConfig );
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

    theHIC.init();
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
