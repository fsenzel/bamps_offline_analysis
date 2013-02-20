#ifndef MESONDECAY_H
#define MESONDECAY_H

#include <stdexcept>
#include <string>
#include "Pythia.h"

class mesonDecay
{
  public:
    mesonDecay ( const int numberElectronStat_arg, const bool muonsInsteadOfElectrons_arg = false, const bool nonPromptJpsiInsteadOfElectrons_arg = false );
    void decayToElectronsToyModel();
    void decayToElectronsPythia();

  private:
    void mesonToElectronToyModel ( double P[4], int& F );
    void setPythiaDecayChannelsElectrons ( Pythia8::Pythia& pythia );
    void setPythiaDecayChannelsMuons ( Pythia8::Pythia& pythia );
    void setPythiaDecayChannelsNonPromptJpsi ( Pythia8::Pythia& pythia );

    bool muonsInsteadOfElectrons;
    bool nonPromptJpsiInsteadOfElectrons;
    int numberElectronStat;
};

/** @brief exception class for handling unexpected critical behaviour within meson decay routines  */
class eMesonDecay_error : public std::runtime_error
{
  public:
    explicit eMesonDecay_error ( const std::string& what ) : std::runtime_error ( what ) {};

    virtual ~eMesonDecay_error() throw() {};
};

#endif // MESONDECAY_H
