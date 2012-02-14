#ifndef MESONDECAY_H
#define MESONDECAY_H

#include "Pythia.h"

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

#endif // MESONDECAY_H
