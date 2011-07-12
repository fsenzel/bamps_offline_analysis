//---------------------------------------------------------------------------------------
//provided by subversion
//---------------------------------------------------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------

/** @file
 * @brief declarations for the Particle class
 */

/** @author Oliver Fochler */


#ifndef PARTICLE_H
#define PARTICLE_H


/**
 * @brief enum type for the flavor of a particle
 *
 * With this, the names given in the enum list can be used instead of integers (though their use is perfectly legal still).
 *
 * Mapped according to:
 * 0 = g (gluon)
 * 1 = u (up)
 * 2 = ub (anti-up)
 * 3 = d (down)
 * 4 = db (anti-down)
 * 5 = s (strange)
 * 6 = sb (anti-strange)
 * 7 = c (charm)
 * 8 = cb (anti-charm)
 * 88 = quark (generic)
 * 99 = anti-quark (generic)
 * 111 = allFlavors (any flavor)
 */
enum FLAVOR_TYPE {gluon, up, anti_up, down, anti_down, strange, anti_strange, charm, anti_charm, quark = 88, anti_quark = 99, allFlavors = 111};


/**
 * @brief Provides properties of a particle needed in.the simulations.
 * @author Oliver Fochler
 *
 * This class encapsulates properties of particles that are needed for the simulation process, such as position and momentum variables.
 * Additionally it overloads the assignment operator and the copy constructor to ease usage.
 */
class Particle
{
public:
    Particle();
    ~Particle() {};
  
    double T_creation;
    double X_init, Y_init,Z_init, X_traveled;//fm
    double PX_init, PY_init,PZ_init, E_init;//GeV
    double X_lastInt, Y_lastInt, Z_lastInt, T_lastInt;
  

    /** @brief time-coordinate [fm/c] */
    double T;
    /** @brief X-coordinate [fm] */
    double X;
    /** @brief Y-coordinate [fm] */
    double Y;
    /** @brief Z-coordinate [fm] */
    double Z;

    /** @brief E, energy [GeV] */
    double E;
    /** @brief p_{x}, momentum in x direction [GeV/c] */
    double PX;
    /** @brief p_{x}, momentum in x direction [GeV/c] */
    double PY;
    /** @brief p_{x}, momentum in x direction [GeV/c] */
    double PZ;

    /** @brief space time rapidity \eta */
    double eta;

    /** @brief mass [GeV] */
    double m;

    /** @brief gluonic screening mass currently associated with the particle object, scaled by alpha_s as md2g / alpha_s [GeV^2] */
    double md2g;
    /** @brief quark screening mass currently associated with the particle object, scaled by alpha_s as md2g / alpha_s [GeV^2] */
    double md2q;

    /** @brief flavor of the particle mapped to int according to FLAVOR_TYPE */
    FLAVOR_TYPE FLAVOR;


    /** @brief the cell the particle belongs to */
    int cell_id;

    /** @brief switch that indicates whether this particle has been annihilated within this timestep */
    bool dead;

    /** @brief unique particle ID */
    int unique_id;
    /** @brief counter for unique particle IDs (static) */
    static long int unique_id_counter;
    /** @brief counter for unique particle IDs of added particles (static) */
    static long int unique_id_counter_added;

    /** @brief TODO */
    bool edge;

    /** @brief TODO */
    long coll_id;

    /** @brief Pythia event number */
    int N_EVENT;
    /** @brief Pythia hard or soft scattering */
    bool HARD; // true/1 if parton comes from hard scattering


    /** @brief collision ordering time [fm] (geometric collisions) */
    double collisionTime;
    /** @brief collision partner (geometric collisions) */
    int collisionPartner;


    /** @brief particle is projectile, special use for jet in box calculations */
    bool isProjectile;

    double Eold, PXold, PYold, PZold;
//     double as22,as23;
//     double rate23,rate32,rate22;//GeV
//     double rate23v,rate32v,rate22v;//GeV
    double rate, ratev;
//     double cs22,cs23,cs22t,cs23t;
//     double lambda_scaled;
//     double md2g_scaled_22;
//     double md2q_scaled_22;
//     double md2g_scaled_23;
    double md2q_scaled_23;
    bool free;
    bool init;
//     int step,tstep,taustep;//fm


    // /** @brief assignment operator */
    //Particle& operator=(Particle const &rhs);

    /** @brief writes the 4-momentum vector to an array */
    void getMomentumArray(double* const) const;
    /** @brief writes the space-time vector to an array */
    void getCoordinateArray(double* const) const;

    static double getMass( const FLAVOR_TYPE _flav ) { return 0; }
    static FLAVOR_TYPE mapToGenericFlavorType( const FLAVOR_TYPE _flav )
    {
      if ( _flav == gluon )
      {
        return gluon;
      }
      else if ( _flav == up || _flav == down || _flav == strange )
      {
        return quark;
      }
      else if ( _flav == anti_up || _flav == anti_down || _flav == anti_strange )
      {
        return anti_quark;
      }      
    }
    
    static int mapToPDGCodes( const FLAVOR_TYPE _flav )
    {
      switch ( _flav )
      {
        case gluon:
          return 21;
          break;
        case up:
          return 2;
          break;
        case anti_up:
          return -2;
          break;
        case down:
          return 1;
          break;
        case anti_down:
          return -1;
          break;
        case strange:
          return 3;
          break;
        case anti_strange:
          return -3;
          break;
        case charm:
          return 4;
          break;
        case anti_charm:
          return -4;
          break;
        default:
          return 0;
          break;
      }      
    }

private:

};


class scatteringInfoStorage22
{
public:
    scatteringInfoStorage22() : iscat( 0 ), jscat( 0 ), time( 0 ), pix( 0 ), piy( 0 ), piz( 0 ), pjx( 0 ), pjy( 0 ), pjz( 0 ) {};
    scatteringInfoStorage22( const int _iscat, const int _jscat, const double _time, const double _pix, const double _piy, const double _piz, const double _pjx, const double _pjy, const double _pjz ) :
            iscat( _iscat ), jscat( _jscat ), time( _time ), pix( _pix ), piy( _piy ), piz( _piz ), pjx( _pjx), pjy( _pjy ), pjz( _pjz ) {};
    ~scatteringInfoStorage22() {};

    int iscat,jscat;
    double time;
    double pix,piy,piz,pjx,pjy,pjz;
};



class scatteringInfoStorage23
{
public:
    scatteringInfoStorage23() : iscat( -1 ), jscat( -1 ), newp( -1 ), time( 0 ), pix( 0 ), piy( 0 ), piz( 0 ), pjx( 0 ), pjy( 0 ), pjz( 0 ),
            newx( 0 ), newy( 0 ), newz( 0 ), newpx( 0 ), newpy( 0 ), newpz( 0 ) {};
    scatteringInfoStorage23( const int _iscat, const int _jscat, const int _newp, const double _time, const double _pix, const double _piy, const double _piz,
            const double _pjx, const double _pjy, const double _pjz,
            const double _newx, const double _newy, const double _newz, const double _newpx, const double _newpy, const double _newpz ) :
            iscat( _iscat ), jscat( _jscat ), newp( _newp ), time( _time ), pix( _pix ), piy( _piy ), piz( _piz ), pjx( _pjx ), pjy( _pjy ), pjz( _pjz ),
            newx( _newx ), newy( _newy ), newz( _newz ), newpx( _newpx ), newpy( _newpy ), newpz( _newpz ) {};
    ~scatteringInfoStorage23() {};

    int iscat,jscat,newp;
    double time;
    double pix,piy,piz,pjx,pjy,pjz;
    double newx,newy,newz,newpx,newpy,newpz;
};



class scatteringInfoStorage32
{
public:
    scatteringInfoStorage32() : iscat( -1 ), jscat( -1 ), dead( -1 ), time( 0 ), pix( 0 ), piy( 0 ), piz( 0 ), pjx( 0 ), pjy( 0 ), pjz( 0 ) {};
    scatteringInfoStorage32( const int _iscat, const int _jscat, const int _dead, const double _time, const double _pix, const double _piy, const double _piz, const double _pjx, const double _pjy, const double _pjz ) :
            iscat( _iscat ), jscat( _jscat ), dead( _dead ), time( _time ), pix( _pix ), piy( _piy ), piz( _piz ), pjx( _pjx), pjy( _pjy ), pjz( _pjz ) {};
    ~scatteringInfoStorage32() {};

    int iscat,jscat,dead;
    double time;
    double pix,piy,piz,pjx,pjy,pjz;
};



class scatteringInfoStorageChange
{
public:
    scatteringInfoStorageChange() : dead( -1 ), newp( -1 ) {};
    scatteringInfoStorageChange( const int _dead, const int _newp ) : dead( _dead ), newp( _newp ) {};
    ~scatteringInfoStorageChange() {};

    int dead,newp;
};

#endif
