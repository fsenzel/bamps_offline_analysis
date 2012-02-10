#ifndef HADRONIZATION_H
#define HADRONIZATION_H

#include "Pythia.h"
using namespace Pythia8; 

class hadronization
{
  public:

    hadronization( const int number_arg );
    void heavyQuarkFragmentation();

  private:

    int number;
    double getFragmentationZ(const int flav);
    double getFragmentationFunction( const double z, const int flav );

};


class mesonDecay
{
  public:

    mesonDecay( const int number_arg );
    void decayToElectronsToyModel();
    void decayToElectronsPythia();

  private:

    int number;
    void mesonToElectronToyModel(double P[4], int & F);
    void setPythiaDecayChannelsElectrons(Pythia& pythia);
    void setPythiaDecayChannelsMuons(Pythia& pythia);
    void setPythiaDecayChannelsNonPromptJpsi(Pythia& pythia);
};


#endif // HADRONIZATION_H
