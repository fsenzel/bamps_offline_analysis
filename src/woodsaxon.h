//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------
//---------------------------------------------


#ifndef WOODSAXON_H
#define WOODSAXON_H

class WoodSaxon
{
  public:
    WoodSaxon() : gamma( 0 ), velocity( 0 ), RA0( 0 ), RA( 0 ), dA( 0 ), n0A( 0 ) {};
    ~WoodSaxon() {};
    
    /** @brief nuclear density (Woods-Saxon distribution) static version */
    static double densityA(double b, double z, const WoodSaxon& _w);
    
    double gamma;
    double velocity;
    double RA0;
    double RA;
    double dA;
    double n0A;
};

#endif // WOODSAXON_H
