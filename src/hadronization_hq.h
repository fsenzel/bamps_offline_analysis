#ifndef HADRONIZATION_HQ_H
#define HADRONIZATION_HQ_H

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

#endif // HADRONIZATION_HQ_H
