#ifndef HADRONIZATION_HQ_H
#define HADRONIZATION_HQ_H

#include "Pythia.h"

class hadronization_hq
{
  public:
    hadronization_hq( const int number_arg ) : number(number_arg) { };
    void heavyQuarkFragmentation();

  private:
    int number;
    double getFragmentationZ(const int flav);
    double getFragmentationFunction( const double z, const int flav );
};


class mesonDecay
{
  public:
    mesonDecay( const int number_arg, const int numberElectronStat_arg, const bool local_cluster_arg = false, const bool muonsInsteadOfElectrons_arg = false, const bool nonPromptJpsiInsteadOfElectrons_arg = false );
    void decayToElectronsToyModel();
    void decayToElectronsPythia();

  private:
    int number;
    void mesonToElectronToyModel(double P[4], int & F);
    void setPythiaDecayChannelsElectrons(Pythia8::Pythia& pythia);
    void setPythiaDecayChannelsMuons(Pythia8::Pythia& pythia);
    void setPythiaDecayChannelsNonPromptJpsi(Pythia8::Pythia& pythia);
    
    bool muonsInsteadOfElectrons;
    bool nonPromptJpsiInsteadOfElectrons;
    bool local_cluster;
    int numberElectronStat;
};


#endif // HADRONIZATION_HQ_H
