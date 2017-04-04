#ifndef HADRONIZATION_HQ_H
#define HADRONIZATION_HQ_H

#include <iostream>
#include <math.h>
#include "particleOffline.h"
#include "random.h"

class hadronization_hq
{
  public:
    std::vector<ParticleOffline> heavyQuarkFragmentation( std::vector<ParticleOffline> particlesToHadronize );

  private:
    double getFragmentationZ(const int flav);
    double getFragmentationFunction( const double z, const int flav );
};

#endif // HADRONIZATION_HQ_H
