#ifndef HADRONIZATION_HQ_H
#define HADRONIZATION_HQ_H

class hadronization_hq
{
  public:
    void heavyQuarkFragmentation();

  private:
    double getFragmentationZ(const int flav);
    double getFragmentationFunction( const double z, const int flav );
};

#endif // HADRONIZATION_HQ_H
