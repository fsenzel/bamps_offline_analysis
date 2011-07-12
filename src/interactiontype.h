//---------------------------------------------------------------------------------------
//provided by subversion
//---------------------------------------------------------------------------------------
//$HeadURL: file:///home/fochler/svnrepos/cascade_full/main/branches/restructure_branch/src/interactiontype.h $
//$LastChangedDate: 2010-07-13 00:25:52 +0200 (Tue, 13 Jul 2010) $
//$LastChangedRevision: 126 $
//$LastChangedBy: fochler $
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------



#ifndef INTERACTIONTYPE_H
#define INTERACTIONTYPE_H

#include <vector>
#include <map>
#include <string>
#include <stdexcept>
#include "particle.h"


/**
 * @brief enum for generic 22, 23 and 32 interactions
 *
 * Mapped according to:
 * 0 = 2 -> 2
 * 1 = 2 -> 3
 * 2 = 3 -> 2
 */
enum GENERIC_COLL_TYPE { c22, c23, c32 };  


/**
 * @brief Provide classification of scattering processes
 *
 * PLEASE NOTE: Many features and routines aren't used (yet). But may be helpful in the future.
 *
 *
 * Interactions types are mapped according to the following scheme:
 *
 * 9221 : g+g -> X
 * 9222 : g+q -> X
 * 9223 : g+qbar -> X
 * 9224 : q+qbar -> X
 * 9225 : q+q -> X
 * 9226 : qbar+qbar -> X
 * 9227 : q+q' -> X
 * 9228 : q+qbar' -> X
 * 9229 : qbar+qbar' -> X
 *
 * 221 : g+g -> g+g
 * 222 : g+g -> q+qbar
 * 223 : g+q -> g+q                2230 : g+qbar -> g+qbar
 * 224 : q+qbar -> q+qbar
 * 225 : q+qbar -> q'+qbar'
 * 226 : q+qbar -> g+g
 * 227 : q+q -> q+q                2270 : qbar+qbar -> qbar+qbar
 * 228 : q+q'-> q+q'               2280 : q+qbar'-> q+qbar'           22800 : qbar+qbar' -> qbar+qbar'
 *
 * 9231 : g+g -> X+g
 * 9232 : g+q -> X+g
 * 9233 : g+qbar -> X+g
 * 9234 : q+qbar -> X+g
 * 9235 : q+q -> X+g
 * 9236 : qbar+qbar -> X+g
 * 9237 : q+q' -> X+g
 * 9238 : q+qbar' -> X+g
 * 9239 : qbar+qbar' -> X+g
 *
 * 231 : g+g -> g+g+g
 * 232 : g+g -> q+qbar+g
 * 233 : g+q -> g+q+g              2330 : g+qbar -> g+qbar+g
 * 234 : q+qbar -> q+qbar+g
 * 235 : q+qbar -> q'+qbar'+g
 * 236 : q+qbar -> g+g+g
 * 237 : q+q -> q+q+g              2370 : qbar+qbar -> qbar+qbar+g
 * 238 : q+q'-> q+q'+g             2380 : q+qbar'-> q+qbar'+g           23800 : qbar+qbar' -> qbar+qbar'+g
 *
 * 9321 : g+g+g -> X
 * 9322 : g+q+g -> X
 * 9323 : g+qbar+g -> X
 * 9324 : q+qbar+g -> X
 * 9325 : q+q+g -> X
 * 9326 : qbar+qbar+g -> X
 * 9327 : q+q'+g -> X
 * 9328 : q+qbar'+g -> X
 * 9329 : qbar+qbar'+g -> X
 *
 * 321 : g+g+g -> g+g
 * 322 : g+g+g -> q+qbar
 * 323 : g+q+g -> g+q              3230 : g+qbar+g -> g+qbar
 * 324 : q+qbar+g -> q+qbar
 * 325 : q+qbar+g -> q'+qbar'
 * 326 : q+qbar+g -> g+g
 * 327 : q+q+g -> q+q              3270 : qbar+qbar+g -> qbar+qbar
 * 328 : q+q'+g-> q+q'             3280 : q+qbar'+g-> q+qbar'         32800 : qbar+qbar'+g -> qbar+qbar'
 *
 */
class interactionType
{
public:
  /** @brief Standard constructor */
  interactionType() : id( 0 ) {};
  /** @brief Standard destructor */
  ~interactionType() {};

  /** @brief id of the interaction type */
  int id;

  /** @brief Returns the number of particles of type _F involved in the initial state of interaction type id */
  double getInvolvedInitial( const FLAVOR_TYPE _F ) const
  {
    return getInvolvedInitial( _F, id );
  }
  /** @brief (static version) Returns the number of particles of type _F involved in the initial state of interaction type id */
  static double getInvolvedInitial( const FLAVOR_TYPE _F, const int _id );

  
  /** @brief Returns the number of particles of type _F involved in the initial state of inclusive interaction type id */
  double getInvolvedInitialInclusive( const FLAVOR_TYPE _F ) const
  {
    return getInvolvedInitialInclusive( _F, id );
  }
  /** @brief (static version) Returns the number of particles of type _F involved in the initial state of inclusive interaction type id */
  static double getInvolvedInitialInclusive( const FLAVOR_TYPE _F, const int _id );
  
  
  /** @brief Returns a vector of interaction type ids that involve particles with flavor _F */
  static std::vector<int> getInteractionTypes( const FLAVOR_TYPE _F );

  /** @brief Generic collision type (2->2 etc) from specific interaction type */
  static GENERIC_COLL_TYPE getCollType( const int _id );
  
  /** @brief Get the index as denoted in processIndices22 (23, 32) for the given process */
  int getIndexFromProcessType( const GENERIC_COLL_TYPE _genType ) const { return getIndexFromProcessType(id, _genType) ;}
  /** @brief (static version) Get the index as denoted in processIndices22 (23, 32) for the given process */
  static int getIndexFromProcessType( const int _id, const GENERIC_COLL_TYPE _genType );
  
  /** @brief (static version) Get the process id at a certain index for a given generic type (22, 23, 32) */
  static int getProcessTypeFromIndex( const int _index, GENERIC_COLL_TYPE _genType );
  /** @brief (static version) Get the inclusive process id at a certain index for a given generic type (22, 23, 32) */
  static int getInclusiveProcessTypeFromIndex( const int _index, GENERIC_COLL_TYPE _genType );
  
  
  /** @brief Get the id for an inclusive process given the initial flavor (2->X) */
  static int getInclusiveProcessType( const FLAVOR_TYPE _F1, const FLAVOR_TYPE _F2, const GENERIC_COLL_TYPE genericType);  
  /** @brief Get the id for an inclusive process given the initial flavor (3->2) */
  static int getInclusiveProcessType( const FLAVOR_TYPE _F1, const FLAVOR_TYPE _F2, const FLAVOR_TYPE _F3, const GENERIC_COLL_TYPE genericType );
  
  
  /** @brief Associates process id numbers for 2->2 processes with indices starting from 0 for use in accessing vector elements etc. */
  static std::map<int, int> processIndices22;
  /** @brief Associates process id numbers for 2->3 processes with indices starting from 0 for use in accessing vector elements etc. */
  static std::map<int, int> processIndices23;
 /** @brief Associates process id numbers for 3->2 processes with indices starting from 0 for use in accessing vector elements etc. */
  static std::map<int, int> processIndices32;
  
  /** @brief Associates indices starting from 0 (for use in accessing vector elements etc.) with id numbers for 2->2 processes  */
  static std::vector<int> indexProcesses22;
  /** @brief Associates indices starting from 0 (for use in accessing vector elements etc.) with id numbers for 2->3 processes  */
  static std::vector<int> indexProcesses23;
  /** @brief Associates indices starting from 0 (for use in accessing vector elements etc.) with id numbers for 3->2 processes  */
  static std::vector<int> indexProcesses32;
  
   /** @brief Associates process id numbers for inclusive (specified initial state, arbitrary final state) 2->2 processes with indices starting from 0 for use in accessing vector elements etc. */
  static std::map<int, int> processIndicesInclusive22;
  /** @brief Associates process id numbers for inclusive (specified initial state, arbitrary final state) 2->3 processes with indices starting from 0 for use in accessing vector elements etc. */
  static std::map<int, int> processIndicesInclusive23;
  /** @brief Associates process id numbers for inclusive (specified initial state, arbitrary final state) 3->2 processes with indices starting from 0 for use in accessing vector elements etc. */
  static std::map<int, int> processIndicesInclusive32;
  
  /** @brief Associates indices starting from 0 (for use in accessing vector elements etc.) with id numbers for inclusive 2->2 processes  */
  static std::vector<int> indexProcessesInclusive22;
  /** @brief Associates indices starting from 0 (for use in accessing vector elements etc.) with id numbers for inclusive 2->3 processes  */
  static std::vector<int> indexProcessesInclusive23;
  /** @brief Associates indices starting from 0 (for use in accessing vector elements etc.) with id numbers for inclusive 3->2 processes  */
  static std::vector<int> indexProcessesInclusive32;
  
  /** @brief Associates process id numbers for inclusive processes with the number of gluons in the initial state */
  static std::map<int, int> involvedInitialGluons;
  /** @brief Associates process id numbers for inclusive processes with the number of quarks in the initial state */
  static std::map<int, int> involvedInitialQuarks;
  /** @brief Associates process id numbers for inclusive processes with the number of anti-quarks in the initial state */
  static std::map<int, int> involvedInitialAntiQuarks;
  
private:
  
};



/** @brief Helper class for initializing std::map container
 *
 * Usage: std::map mymap = create_map<int, int >(1,2)(3,4)(5,6);
 */
template <typename T, typename U>
class create_map
{
private:
  std::map<T, U> m_map;
  
public:
  create_map( const T& key, const U& val )
  {
    m_map[key] = val;
  }

  create_map<T, U>& operator()( const T& key, const U& val )
  {
    m_map[key] = val;
    return *this;
  }

  operator std::map<T, U>()
  {
    return m_map;
  }
};



/** @brief Helper class for initializing std::vector container
 *
 * Usage: std::vector myvec = create_vector<int>(1)(2)(3);
 */
template <typename T>
class create_vector
{
private:
  std::vector<T> m_vec;
  
public:
  create_vector( const T& val )
  {
    m_vec.push_back( val );
  }

  create_vector<T>& operator()( const T& val )
  {
    m_vec.push_back( val );
    return *this;
  }

  operator std::vector<T>()
  {
    return m_vec;
  }
};





/** @brief exception class for handling unexpected behaviour when dealing with interaction types */
class eInteraction_type_error : public std::runtime_error
{
  public:
    explicit eInteraction_type_error(const std::string& what) : std::runtime_error(what) {};

    virtual ~eInteraction_type_error() throw() {};
};


#endif // INTERACTIONTYPE_H
// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;
