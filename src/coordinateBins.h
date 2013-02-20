//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------
//---------------------------------------------


#ifndef COORDINATE_BINS_H
#define COORDINATE_BINS_H

#include <stdexcept>
#include <vector>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/access.hpp>

#include "woodsaxon.h"

using std::vector;



class coordinateSubBin
{
  public:
    coordinateSubBin() : left ( 0 ), right ( 0 ), content ( 0 ) {};
    ~coordinateSubBin() {};

    void setLeftRight ( const double _l, const double _r )
    {
      left = _l;
      right = _r;
    }

    coordinateSubBin& operator++()
    {
      ++content;  //prefix
      return *this;
    }
    coordinateSubBin operator++ ( int unused )
    {
      coordinateSubBin temp = *this;  //postfix
      ++content;
      return temp;
    }

    void clear()
    {
      left = 0;
      right = 0;
      content = 0;
    }

    double left;
    double right;
    int content;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize ( Archive& ar, const unsigned int version )
    {
      ar& left;
      ar& right;
      ar& content;
    }
};



class coordinateBins
{
  public:
    coordinateBins() : _negative_indices ( false ), _min_real ( 0 ), _max_real ( 0 ), _min_index_limit ( 0 ),
      _max_index_limit ( 0 ),  _min_index_active ( 0 ), _max_index_active ( 0 ), _delta_x ( 0 )
    {
      bins.clear();
    }
    coordinateBins ( const int _size, const double _min, const double _max );
    ~coordinateBins() {};

    double min_real()
    {
      return _min_real;
    }
    double max_real()
    {
      return _max_real;
    }

    int min_index()
    {
      return _min_index_active;
    }
    int max_index()
    {
      return _max_index_active;
    }

    int min_index_limit()
    {
      return _min_index_limit;
    }
    int max_index_limit()
    {
      return _max_index_limit;
    }

    int size()
    {
      return bins.size();
    }
    void clear()
    {
      *this = coordinateBins ( 0, 0 , 0 );
    }
    void reshape ( const int _size, const double _min, const double _max );
    double get_dx()
    {
      return _delta_x;
    }

    const coordinateSubBin& operator[] ( const int _index ) const;

    int getIndex ( const double _x ) const;

    void increase ( const double _x )
    {
      ++bins[ getIndex ( _x ) ];
    }



  protected:
    vector< coordinateSubBin > bins;

    double _min_real;
    double _max_real;
    double _delta_x;

    int _min_index_limit;
    int _max_index_limit;

    int _min_index_active;
    int _max_index_active;

    bool _negative_indices;

  private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize ( Archive& ar, const unsigned int version )
    {
      ar& bins;
      ar& _min_real;
      ar& _max_real;
      ar& _delta_x;
      ar& _min_index_limit;
      ar& _max_index_limit;
      ar& _min_index_active;
      ar& _max_index_active;
      ar& _negative_indices;
    }
};



class coordinateEtaBins : public coordinateBins
{
  public:
    coordinateEtaBins() : coordinateBins(), NinEtaBin ( 0 ) {};
    coordinateEtaBins ( const int _size, const double _min, const double _max ) : coordinateBins ( _size, _min, _max ), NinEtaBin ( 0 ) {};
    ~coordinateEtaBins() {};

    void populateEtaBins ( coordinateBins& _dNdEta, const double _etaShift, double _timenow, double& _dt, const double _dx, const double _dEta_fine
                         );
    int constructEtaBins ( const int _NperCell, const double _b, const double _dx, const double _dy, WoodSaxon& _param, const int _nTest );
    int getCentralIndex() const;
    int getIndex ( const double _eta ) const;

    int getNinEtaBin() const
    {
      return NinEtaBin;
    }

  private:
    int NinEtaBin;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize ( Archive& ar, const unsigned int version )
    {
      ar& boost::serialization::base_object<coordinateBins> ( *this );
      ar& NinEtaBin;
    }
};



/** @brief exception class for handling unexpected critical behaviour within simulations of heavy ion collisions  */
class eCoordBins_error : public std::runtime_error
{
  public:
    explicit eCoordBins_error ( const std::string& what ) : std::runtime_error ( what ) {};

    virtual ~eCoordBins_error() throw() {};
};


#endif
