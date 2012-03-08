#ifndef MESONDECAY_H
#define MESONDECAY_H

#include "Pythia.h"

class mesonDecay
{
  public:
    mesonDecay( const int numberElectronStat_arg, const bool muonsInsteadOfElectrons_arg = false, const bool nonPromptJpsiInsteadOfElectrons_arg = false );
    void decayToElectronsToyModel();
    void decayToElectronsPythia();

  private:
    void mesonToElectronToyModel(double P[4], int & F);
    void setPythiaDecayChannelsElectrons(Pythia8::Pythia& pythia);
    void setPythiaDecayChannelsMuons(Pythia8::Pythia& pythia);
    void setPythiaDecayChannelsNonPromptJpsi(Pythia8::Pythia& pythia);
    
    bool muonsInsteadOfElectrons;
    bool nonPromptJpsiInsteadOfElectrons;
    int numberElectronStat;
};

#endif // MESONDECAY_H
