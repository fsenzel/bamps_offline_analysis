#include "mesonDecay.h"

using namespace Pythia8;
using namespace std;

mesonDecay::mesonDecay( const int numberElectronStat_arg, const bool muonsInsteadOfElectrons_arg, const bool nonPromptJpsiInsteadOfElectrons_arg ) :
  numberElectronStat(numberElectronStat_arg), muonsInsteadOfElectrons(muonsInsteadOfElectrons_arg), nonPromptJpsiInsteadOfElectrons(nonPromptJpsiInsteadOfElectrons_arg)
{
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
std::vector<ParticleOffline> mesonDecay::decayToElectronsPythia( std::vector<ParticleOffline> hadronsToDecay )
{
  vector<ParticleOffline> decayElectrons;
  
  /**
  * Reserve memory for the Particle vector.
  */
  decayElectrons.reserve( hadronsToDecay.size() * numberElectronStat );
  /**
  * Now the addedPartcl_electron vector is re-sized to hold the needed number of particles (this is NOT done when reserving memory!).
  * The particles are initialised with the standard constructor, actual attributes such as momentum etc. MUST be set later!
  */
  decayElectrons.resize( hadronsToDecay.size() * numberElectronStat );

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

  double p_tmp[4];

  for ( int i = 0; i < hadronsToDecay.size(); i++ )
  {
    if ( ParticlePrototype::mapToGenericFlavorType( hadronsToDecay[i].FLAVOR ) == dmeson_gen || ParticlePrototype::mapToGenericFlavorType( hadronsToDecay[i].FLAVOR ) == bmeson_gen ) // d or b meson
    {
      // each d meson decays to numberElectronStat electrons
      for( int k = 0; k < numberElectronStat; k++ )
      {      
        int k_e = ( i ) * numberElectronStat + k ;
        // do not copy all particle properties, just E, p, M, Flavor (see below))
//         addedPartcl_electron[k_e]=addedParticlesCopy[i]; // copy all added particles, perform decay here (important to have a second instant for initial electrons, simulating pp collisions)
        
        //Reset event record to allow for new event.
        event.reset();
        
        int id = hadronsToDecay[i].FLAVOR;

        // mass from Pythia
        double mm = pdt.mSel(id);
        //calculate energy with new mass
        double ee = sqrt( hadronsToDecay[i].Mom.Perp2() + mm*mm );
        // Store the particle in the event record.
        // definition from event.h: int append(int id, int status, int color, int anticolor, double px, double py, double pz, double e, double m = 0.)
        event.append(  id, 1, 0,   0, hadronsToDecay[i].Mom.Px(), hadronsToDecay[i].Mom.Py(),  hadronsToDecay[i].Mom.Pz(), ee, mm); // add particle to event, status=1 ensures that no message about non vanishing total charge pops up

        // Generate events. Quit if failure.
        if (!pythia.next()) 
        {
          std::string errMsg = "PYTHIA event generation aborted prematurely, owing to error!";
          throw eMesonDecay_error( errMsg );
        }
    
//         if(k_e < 10)
//         {
//           // List first few events.
//           event.list();
//         }
        
        // for theta analysis
        p_tmp[1] = hadronsToDecay[i].Mom.Px();
        p_tmp[2] = hadronsToDecay[i].Mom.Py();
        p_tmp[3] = hadronsToDecay[i].Mom.Pz();        
        
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
                decayElectrons[k_e].FLAVOR = electron;
              else
                decayElectrons[k_e].FLAVOR = positron;
              
              decayElectrons[k_e].m = event[j].m();
              decayElectrons[k_e].Mom = VectorEPxPyPz( event[i].e(),
                                                             event[i].px(),
                                                             event[i].py(),
                                                             event[i].pz() );
              decayElectrons[k_e].Mom.SetEbyM( event[j].m() );
              
              found_electron = true;
            }
            else if( event[j].isFinal() && found_electron && ( event[j].mother1()==1 || event[j].mother2()==1 ) )
            {
              cout << "error in meson decay to electron. Additional electron found for particle ID " << hadronsToDecay[i].FLAVOR << "  (411&421 D meson, 511&521 B meson)" << endl;
              // List event
//               event.list();
            }
          }
        }
        
        if( !found_electron && !( ParticlePrototype::mapToGenericFlavorType( hadronsToDecay[i].FLAVOR ) == dmeson_gen && nonPromptJpsiInsteadOfElectrons ) )
        {
          cout << "error in meson decay to electron. No electron found" << endl;
          // List event
          event.list();
        }
      }
    }
  }

  // pythia changes output format. Set it back 
  cout.precision (6);
  cout.setf(ios::floatfield);
}


// use PYTHIA for Hadronisation and Meson decay
/*
  event.append(  4, 300, 101,   0, 0., 0.,  pp, ee, mm); 
  event.append(  -4, 300, 0,   101, 0., 0.,  0., mm, mm); */




