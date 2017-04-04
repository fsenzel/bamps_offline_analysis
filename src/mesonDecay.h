#ifndef MESONDECAY_H
#define MESONDECAY_H

#include <iostream>
#include <math.h>
#include <iomanip>
#include <stdexcept>
#include <string>

#include "Pythia8/Pythia.h"

#include "particleOffline.h"
#include "configBAMPS.h"

class mesonDecay
{
  public:
    mesonDecay( const int numberElectronStat_arg, const bool muonsInsteadOfElectrons_arg = false, const bool nonPromptJpsiInsteadOfElectrons_arg = false );
    void decayToElectronsToyModel();
    std::vector<ParticleOffline> decayToElectronsPythia( std::vector<ParticleOffline> hadronsToDecay );

  private:
    void mesonToElectronToyModel(double P[4], int & F);
    void setPythiaDecayChannelsElectrons(Pythia8::Pythia& pythia);
    void setPythiaDecayChannelsMuons(Pythia8::Pythia& pythia);
    void setPythiaDecayChannelsNonPromptJpsi(Pythia8::Pythia& pythia);
    
    bool muonsInsteadOfElectrons;
    bool nonPromptJpsiInsteadOfElectrons;
    int numberElectronStat;
};

/** @brief exception class for handling unexpected critical behaviour within meson decay routines  */
class eMesonDecay_error : public std::runtime_error
{
  public:
    explicit eMesonDecay_error(const std::string& what) : std::runtime_error(what) {};
    
    virtual ~eMesonDecay_error() throw() {};
};

#endif // MESONDECAY_H
