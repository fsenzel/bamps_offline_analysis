#include <iostream>
#include <math.h>
#include <iomanip>
#include "hadronization_hq.h"
#include "particle.h"
#include "random.h"
#include "binning.h"
#include "configuration.h"

#include "config.h"

#include "Pythia.h"
using namespace Pythia8; 
using namespace ns_casc;
using namespace std;



void hadronization_hq::heavyQuarkFragmentation()
{
  double z;
  int flav;

  binning bins( "output/binsFragment.dat", 0.0, 1.0, 400 );
  for ( int i = 1; i <= number; i++ )
  {
    addedParticlesCopy[i]=addedParticles[i]; // copy all added particles, perform hadronization here (important to have a second instant for initial hadronization, simulating pp collisions)
     
    flav = addedParticlesCopy[i].FLAVOR;
    if ( flav >= 7 && flav <= 10 )
    {
      z = getFragmentationZ( flav );
      addedParticlesCopy[i].PX = addedParticlesCopy[i].PX * z;
      addedParticlesCopy[i].PY = addedParticlesCopy[i].PY * z;
      addedParticlesCopy[i].PZ = addedParticlesCopy[i].PZ * z;
      
      // D+ = c dbar, D0 = c ubar, D0bar = cbar u, D- = cbar d
      // B+ = bbar u, B0 = bbar d, B0bar = b dbar, B- = b ubar 
      // for flavor use PDG particle codes: 
      // 411  D+ ,   -411  D- ,   421  D0 ,  -421  D0bar
      // 521  B+ ,   -521  B- ,   511  B0 ,  -511  B0bar         
      const double fraction = 0.31; // fraction of c->D+ / (c->D0  +  c->D+), from PHENIX paper 1005.1627
      if( flav == 7 )
      {
        if( ran2() < fraction)
          addedParticlesCopy[i].FLAVOR = dmeson_plus;
        else
          addedParticlesCopy[i].FLAVOR = dmeson_zero;
      }
      else if( flav == 8 )
      {
        if( ran2() < fraction)
          addedParticlesCopy[i].FLAVOR = dmeson_minus;
        else
          addedParticlesCopy[i].FLAVOR = dmeson_zero_bar;
      }
      else if( flav == 9 )
      {
        if( ran2() < fraction)
          addedParticlesCopy[i].FLAVOR = bmeson_zero_bar;
        else
          addedParticlesCopy[i].FLAVOR = bmeson_minus;
      }
      else if( flav == 10 )
      {
        if( ran2() < fraction)
          addedParticlesCopy[i].FLAVOR = bmeson_zero;
        else
          addedParticlesCopy[i].FLAVOR = bmeson_plus;
      }
      
      addedParticlesCopy[i].m = ParticlePrototype::getMass( addedParticlesCopy[i].FLAVOR );        
        
      addedParticlesCopy[i].E = sqrt( pow( addedParticlesCopy[i].PX, 2.0 ) + pow( addedParticlesCopy[i].PY, 2.0 ) + pow( addedParticlesCopy[i].PZ, 2.0 ) + pow( addedParticlesCopy[i].m, 2.0 ) );

//       cout << addedParticlesCopy[i].PZ  << endl;

      bins.add( addedParticlesCopy[i].PZ );
    }
  }

//   bins.print();
}


// samples z with the metropolis algoritm according to fragmentation function getFragmentationFunction(z)
double hadronization_hq::getFragmentationZ(const int flav)
{
  double z;
  double z_new;
  double g = 0, g_new = 0;
  double ratio;
  double r;

  // the range in which the variable z needs to be sampled
  const double z_min = 0.0;
  const double z_max = 1.0;

  // select initial values of z
  do
  {
//     z = z_min + ran2() * ( z_max - z_min );
    z = pow( ran2() , 1.0 / 17.0 ); // inverse function of integrated comparison function x^16. Thus, samples z according to z^16
    g = getFragmentationFunction( z, flav );
  }
  while ( g == 0.0 );


  // number of steps in the Markov chain
  const int n_steps = 50;//50;

  // do n_steps steps
  // the current location is (t)
  // one steps consists of
  //  - proposing a new point (t') according to a certain proposal distribution
  //  - calculate the matrix element g(t') at this new point
  //  - if g(t') > g(t) accept the new point (t')
  //  - else accept the new point (t') with a probability g(t') / g(t) (last point is modified since function for sampling t is not symmetric
  //
  for ( int i = 0; i < n_steps; i++ )
  {
    do
    {
//       z_new = z_min + ran2() * ( z_max - z_min );      // propose new z using a uniform distribution over the entire range
      z_new = pow( ran2() , 1.0 / 17.0 );
    }
    while ( z_new < z_min || z_new > z_max );

    g_new = getFragmentationFunction( z_new, flav );            // calculate the matrix element at the proposed point

    ratio = g_new / g;                                    // ratio of g(z') / g(z)

    ratio = ratio * ( pow( z, 16.0 ) ) / ( pow( z_new, 16.0 ) ); // necessary if one does not use a symmetric propose fct. for t_new

    if ( ratio >= 1.0 || ran2() < ratio )     // accept if g(z') > g(z) or with probability "ratio"
    {
      z = z_new;
      g = g_new;
    }
  }

  return z;
}


// use Peterson fragmentation, cf. Peterson et al. Phys. Rev. D 27, 105-111 (1983)
double hadronization_hq::getFragmentationFunction( const double z, const int flav )
{
  double D, epsilon;
//   const double m_q = 0.02;
  const double epsilon_c = 0.05; //pow( m_q / Mcharm , 2.0 ); // standard values: epsilon_charm = 0.05, epsilon_bottom = 0.005
  const double epsilon_b = 0.005;
  if(flav == 7 || flav == 8)
    epsilon = epsilon_c;
  else if(flav == 9 || flav == 10)
    epsilon = epsilon_b;
  else
    cout << "error in hadronization_hq::getFragmentationFunction( const double z, const int flav )" << endl;

  D = 1.0 / ( z * pow( 1.0 - 1.0 / z - epsilon / ( 1.0 - z ) , 2.0 ) );

  return D;
}





mesonDecay::mesonDecay( const int number_arg, const int numberElectronStat_arg, const bool local_cluster_arg, const bool muonsInsteadOfElectrons_arg, const bool nonPromptJpsiInsteadOfElectrons_arg ) :
  number(number_arg), numberElectronStat(numberElectronStat_arg), local_cluster(local_cluster_arg), muonsInsteadOfElectrons(muonsInsteadOfElectrons_arg), nonPromptJpsiInsteadOfElectrons(nonPromptJpsiInsteadOfElectrons_arg)
{
  /**
  * Reserve memory for the Particle vector.
  */
  addedPartcl_electron.reserve( number * numberElectronStat );
  /**
  * Now the addedPartcl_electron vector is re-sized to hold the needed number of particles (this is NOT done when reserving memory!).
  * The particles are initialised with the standard constructor, actual attributes such as momentum etc. MUST be set later!
  */
  addedPartcl_electron.resize( number * numberElectronStat );
}

// void mesonDecay::decayToElectronsToyModel()
// {
//   double P[4];
//   int F;
// 
//   binning bins( "output/binsDecay.dat", 0.0, 1.0, 400 );
// 
//   for ( int i = 1; i <= number; i++ )
//   {
//     if ( addedParticlesCopy[i].FLAVOR == 700 || addedParticlesCopy[i].FLAVOR == 800 ) // d meson
//     {
//       P[0] = addedParticlesCopy[i].E;
//       P[1] = addedParticlesCopy[i].PX;
//       P[2] = addedParticlesCopy[i].PY;
//       P[3] = addedParticlesCopy[i].PZ;
//       F = addedParticlesCopy[i].FLAVOR;
//       
//       mesonToElectronToyModel( P, F );
//       
//       addedParticlesCopy[i].PX = P[1];
//       addedParticlesCopy[i].PY = P[2];
//       addedParticlesCopy[i].PZ = P[3];
//       addedParticlesCopy[i].FLAVOR = 1000; // electron
//       addedParticlesCopy[i].m = 0.000511; //electrons: 511 keV
//       addedParticlesCopy[i].E = sqrt( pow( addedParticlesCopy[i].PX, 2.0 ) + pow( addedParticlesCopy[i].PY, 2.0 ) + pow( addedParticlesCopy[i].PZ, 2.0 ) + pow( addedParticlesCopy[i].m, 2.0 ) );
// 
// //       cout << addedParticlesCopy[i].PZ  << endl;
// 
//       bins.add( addedParticlesCopy[i].PZ );
//     }
//   }
// 
//   bins.print();
// }


// use PYTHIA for Meson decay
void mesonDecay::decayToElectronsPythia()
{
  double p_tmp[4];
  int id, k_e;
  double mm,ee;

  binning bins( "output/binsDecay.dat", -3.0, 3.0, 400 );
  binning binsElec( "output/binsElectron.dat", 1000., 1002., 3 );
  binning binsTheta( "output/binsTheta.dat", 0., 3.15, 200 );

  // give correct xml directory to pythia object, only possible via constructor
  string xmlpath = PYTHIA_XML_DIR;
  
  // Generator; shorthand for event and particleData.                           
  Pythia pythia( xmlpath );
  Event& event      = pythia.event;
  ParticleData& pdt = pythia.particleData;

  // Key requirement: switch off ProcessLevel, and thereby also PartonLevel.
  pythia.readString("ProcessLevel:all = off");

  // Optionally switch off decays.
  //pythia.readString("HadronLevel:Decay = off");

  // Provide printout of initial information.        
  pythia.settings.listChanged();
 
  // Initialize.
  pythia.init();
  
  
  // switch off all decay channels for charmed mesons except decay to electrons
  if( nonPromptJpsiInsteadOfElectrons ) // only decay to muons
    setPythiaDecayChannelsNonPromptJpsi( pythia );
  else if( muonsInsteadOfElectrons ) // only decay to muons
    setPythiaDecayChannelsMuons( pythia );
  else // only decay to electrons
    setPythiaDecayChannelsElectrons( pythia );
  
//   pythia.particleData.listChanged();
//   pythia.particleData.list(443);
  

  for ( int i = 1; i <= number; i++ )
  {
    if ( ParticlePrototype::mapToGenericFlavorType( addedParticlesCopy[i].FLAVOR ) == dmeson_gen || ParticlePrototype::mapToGenericFlavorType( addedParticlesCopy[i].FLAVOR ) == bmeson_gen ) // d or b meson
    {
      // each d meson decays to numberElectronStat electrons
      for(int k = 1; k <= numberElectronStat; k++)
      {      
        k_e = ( i - 1 ) * numberElectronStat + k ;
        // do not copy all particle properties, just E, p, M, Flavor (see below))
//         addedPartcl_electron[k_e]=addedParticlesCopy[i]; // copy all added particles, perform decay here (important to have a second instant for initial electrons, simulating pp collisions)
        
        //Reset event record to allow for new event.
        event.reset();
        
        id = addedParticlesCopy[i].FLAVOR;

        // mass from Pythia
        mm = pdt.mass(id);
        //calculate energy with new mass
        ee = sqrt( pow( addedParticlesCopy[i].PX, 2.0 ) + pow( addedParticlesCopy[i].PY, 2.0 ) + pow( addedParticlesCopy[i].PZ, 2.0 ) + pow( mm, 2.0 ) );
        // Store the particle in the event record.
        // definition from event.h: int append(int id, int status, int color, int anticolor, double px, double py, double pz, double e, double m = 0.)
        event.append(  id, 1, 0,   0, addedParticlesCopy[i].PX, addedParticlesCopy[i].PY,  addedParticlesCopy[i].PZ, ee, mm); // add particle to event, status=1 ensures that no message about non vanishing total charge pops up

        // Generate events. Quit if failure.
        if (!pythia.next()) {
          cout << " Event generation aborted prematurely, owing to error!\n"; 
        }
    
//         if(k_e < 10)
//         {
//           // List first few events.
//           event.list();
//         }
        
        // for theta analysis
        p_tmp[1] = addedParticlesCopy[i].PX;
        p_tmp[2] = addedParticlesCopy[i].PY;
        p_tmp[3] = addedParticlesCopy[i].PZ;
	
        
        bool found_electron = false; // just for error checking if there are 2 electrons
        // Loop over all particles from this decay and search for the electron
        for (int j = 0; j < event.size(); ++j) 
        {
          // search for electron or positron, if the decay to muons is switched on, search for them instead. In the cascade however they get the ID of an electron to avoid writing the analysis routines new
          if( ( event[j].idAbs() == 11 && !muonsInsteadOfElectrons && !nonPromptJpsiInsteadOfElectrons ) || ( event[j].idAbs() == 13 && muonsInsteadOfElectrons && !nonPromptJpsiInsteadOfElectrons ) || ( event[j].idAbs() == 443 && !muonsInsteadOfElectrons && nonPromptJpsiInsteadOfElectrons ) ) 
          {
            if( event[j].isFinal() && !found_electron && ( event[j].mother1()==1 || event[j].mother2()==1 ) )
            {
              if( event[j].id() > 0) // electron
                addedPartcl_electron[k_e].FLAVOR = electron;
              else
                addedPartcl_electron[k_e].FLAVOR = positron;
              
              addedPartcl_electron[k_e].PX = event[j].px();
              addedPartcl_electron[k_e].PY = event[j].py();
              addedPartcl_electron[k_e].PZ = event[j].pz();
              
              addedPartcl_electron[k_e].m = event[j].m();
              addedPartcl_electron[k_e].E = sqrt( pow( addedPartcl_electron[k_e].PX, 2.0 ) + pow( addedPartcl_electron[k_e].PY, 2.0 ) + pow( addedPartcl_electron[k_e].PZ, 2.0 ) + pow( addedPartcl_electron[k_e].m, 2.0 ) );
              
              found_electron = true;
            }
            else if( event[j].isFinal() && found_electron && ( event[j].mother1()==1 || event[j].mother2()==1 ) )
            {
              cout << "error in meson decay to electron. Additional electron found for particle ID " << addedParticlesCopy[i].FLAVOR << "  (411&421 D meson, 511&521 B meson)" << endl;
              // List event
//               event.list();
            }
          }
        }
        
        if( !found_electron && !( ParticlePrototype::mapToGenericFlavorType( addedParticlesCopy[i].FLAVOR ) == dmeson_gen && nonPromptJpsiInsteadOfElectrons ) )
        {
          cout << "error in meson decay to electron. No electron found" << endl;
          // List event
          event.list();
        }
  //       else
  //       {
  //         bins.add( addedPartcl_electron[k_e].PZ );
  //         binsElec.add( addedPartcl_electron[k_e].FLAVOR );
  //         
  //         double costheta = ( p_tmp[1]*addedPartcl_electron[k_e].PX + p_tmp[2]*addedPartcl_electron[k_e].PY + p_tmp[3]*addedPartcl_electron[k_e].PZ ) 
  //         /  sqrt( pow( addedPartcl_electron[k_e].PX, 2.0 ) + pow( addedPartcl_electron[k_e].PY, 2.0 ) + pow( addedPartcl_electron[k_e].PZ, 2.0 ) ) 
  //         / sqrt( pow( p_tmp[1], 2.0 ) + pow( p_tmp[2], 2.0 ) + pow( p_tmp[3], 2.0 ) );
  //         binsTheta.add( acos( costheta ) );
  //       }


      }
    }
  }

//   bins.print();
//   binsElec.print();
//   binsTheta.print();

  // pythia changes output format. Set it back 
  cout.precision (6);
  cout.setf(ios::floatfield);


}


// use PYTHIA for Hadronisation and Meson decay
/*
  event.append(  4, 300, 101,   0, 0., 0.,  pp, ee, mm); 
  event.append(  -4, 300, 0,   101, 0., 0.,  0., mm, mm); */




