//---------------------------------------------------------------------------------------
//provided by subversion
//---------------------------------------------------------------------------------------
//$HeadURL$
//$LastChangedDate$
//$LastChangedRevision$
//$LastChangedBy$
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------


#include "interactiontype.h"

/**
 * Initialize map containers that associate interaction type id numbers with index numbers for use in vector access etc.
 * When modifications are necessary, also change indexProcessesXX (see below)!
 */
std::map<int, int> interactionType::processIndices22 = create_map<int, int >( 221, 0 )( 222, 1 )( 223, 2 )( 2230, 3 )( 224, 4 )( 225, 5 )( 226, 6 )( 227, 7 )( 2270, 8 )( 228, 9 )( 2280, 10 )( 22800, 11 );
std::map<int, int> interactionType::processIndices23 = create_map<int, int >( 231, 0 )( 232, 1 )( 233, 2 )( 2330, 3 )( 234, 4 )( 235, 5 )( 236, 6 )( 237, 7 )( 2370, 8 )( 238, 9 )( 2380, 10 )( 23800, 11 );
std::map<int, int> interactionType::processIndices32 = create_map<int, int >( 321, 0 )( 322, 1 )( 323, 2 )( 3230, 3 )( 324, 4 )( 325, 5 )( 326, 6 )( 327, 7 )( 3270, 8 )( 328, 9 )( 3280, 10 )( 32800, 11 );


/**
 * Initialize vector containers that associate index numbers (for use in vector access etc.) with interaction type id numbers
 * When modifications are necessary, also change processIndicesXX (see above)!
 */
std::vector<int> interactionType::indexProcesses22 = create_vector<int>( 221 )( 222 )( 223 )( 2230 )( 224 )( 225 )( 226 )( 227 )( 2270 )( 228 )( 2280 )( 22800 );
std::vector<int> interactionType::indexProcesses23 = create_vector<int>( 231 )( 232 )( 233 )( 2330 )( 234 )( 235 )( 236 )( 237 )( 2370 )( 238 )( 2380 )( 23800 );
std::vector<int> interactionType::indexProcesses32 = create_vector<int>( 321 )( 322 )( 323 )( 3230 )( 324 )( 325 )( 326 )( 327 )( 3270 )( 328 )( 3280 )( 32800 );


/**
 * Initialize map containers that associate interaction type id numbers for inclusive processes (specified initial state, arbitrary final state,
 * IDs starting with 9) with index numbers for use in vector access etc.
 * When modifications are necessary, also change indexProcessesInclusiveXX (see below)!
 */
std::map<int, int> interactionType::processIndicesInclusive22 = create_map<int, int >( 9221, 0 )( 9222, 1 )( 9223, 2 )( 9224, 3 )( 9225, 4 )( 9226, 5 )( 9227, 6 )( 9228, 7 )( 9229, 8 );
std::map<int, int> interactionType::processIndicesInclusive23 = create_map<int, int >( 9231, 0 )( 9232, 1 )( 9233, 2 )( 9234, 3 )( 9235, 4 )( 9236, 5 )( 9237, 6 )( 9238, 7 )( 9239, 8 );
std::map<int, int> interactionType::processIndicesInclusive32 = create_map<int, int >( 9321, 0 )( 9322, 1 )( 9323, 2 )( 9324, 3 )( 9325, 4 )( 9326, 5 )( 9327, 6 )( 9328, 7 )( 9329, 8 );


/**
 * Initialize vector containers that associate index numbers (for use in vector access etc.) with inclusive interaction type id numbers
 * When modifications are necessary, also change processIndicesInclusiveXX (see above)!
 */
std::vector<int> interactionType::indexProcessesInclusive22 = create_vector<int>( 9221 )( 9222 )( 9223 )( 9224 )( 9225 )( 9226 )( 9227 )( 9228 )( 9229 );
std::vector<int> interactionType::indexProcessesInclusive23 = create_vector<int>( 9231 )( 9232 )( 9233 )( 9234 )( 9235 )( 9236 )( 9237 )( 9238 )( 9239 );
std::vector<int> interactionType::indexProcessesInclusive32 = create_vector<int>( 9321 )( 9322 )( 9323 )( 9324 )( 9325 )( 9326 )( 9327 )( 9328 )( 9329 );


/**
* Initialize map containers that associate interaction type id numbers with the number of involved particles in the initial state.
* So far only IDs starting with 9, i.e. inclusive processes
*/
std::map<int, int> interactionType::involvedInitialGluons = create_map<int, int >( 9221, 2 )( 9222, 1 )( 9223, 1 )( 9224, 0 )( 9225, 0 )( 9226, 0 )( 9227, 0 )( 9228, 0 )( 9229, 0 )
    ( 9231, 2 )( 9232, 1 )( 9233, 1 )( 9234, 0 )( 9235, 0 )( 9236, 0 )( 9237, 0 )( 9238, 0 )( 9239, 0 )
    ( 9321, 3 )( 9322, 2 )( 9323, 2 )( 9324, 1 )( 9325, 1 )( 9326, 1 )( 9327, 1 )( 9328, 1 )( 9329, 1 );
std::map<int, int> interactionType::involvedInitialQuarks = create_map<int, int >( 9221, 0 )( 9222, 1 )( 9223, 0 )( 9224, 1 )( 9225, 2 )( 9226, 0 )( 9227, 2 )( 9228, 1 )( 9229, 0 )
    ( 9231, 0 )( 9232, 1 )( 9233, 0 )( 9234, 1 )( 9235, 2 )( 9236, 0 )( 9237, 2 )( 9238, 1 )( 9239, 0 )
    ( 9321, 0 )( 9322, 1 )( 9323, 0 )( 9324, 1 )( 9325, 2 )( 9326, 0 )( 9327, 2 )( 9328, 1 )( 9329, 0 );
std::map<int, int> interactionType::involvedInitialAntiQuarks = create_map<int, int >( 9221, 0 )( 9222, 0 )( 9223, 1 )( 9224, 1 )( 9225, 0 )( 9226, 2 )( 9227, 0 )( 9228, 1 )( 9229, 2 )
    ( 9231, 0 )( 9232, 0 )( 9233, 1 )( 9234, 1 )( 9235, 0 )( 9236, 2 )( 9237, 0 )( 9238, 1 )( 9239, 2 )
    ( 9321, 0 )( 9322, 0 )( 9323, 1 )( 9324, 1 )( 9325, 0 )( 9326, 2 )( 9327, 0 )( 9328, 1 )( 9329, 2 );






/**
 * Look up the process denoted by _id in the std::maps and return the appropriate index number.
 * The generic type (22, 23, 32) is also returned via a parameter.
 *
 * @param[in] _id process id to be looked up (231, 32800 etc.)
 * @param[in,out] _genType generic collision type the specific process belongs to
 * @return the index number as denoted in the std:maps
 */
int interactionType::getIndexFromProcessType( const int _id, const GENERIC_COLL_TYPE _genType )
{
  std::map<int, int>::const_iterator find_iter;

  switch ( _genType )
  {
  case c22:
    find_iter = processIndicesInclusive22.find( _id );
    if ( find_iter != processIndicesInclusive22.end() )
    {
      return find_iter->second;
    }
    else
    {
      find_iter = processIndices22.find( _id );
      if ( find_iter != processIndices22.end() )
      {
        return find_iter->second;
      }
      else
      {
        std::string errMsg = "Determination of process index for 2->2 process failed. Unrecoverable error.";
        throw eInteraction_type_error( errMsg );
      }
    }
    break;
  case c23:
    find_iter = processIndicesInclusive23.find( _id );
    if ( find_iter != processIndicesInclusive23.end() )
    {
      return find_iter->second;
    }
    else
    {
      find_iter = processIndices23.find( _id );
      if ( find_iter != processIndices23.end() )
      {
        return find_iter->second;
      }
      else
      {
        std::string errMsg = "Determination of process index for 2->3 process failed. Unrecoverable error.";
        throw eInteraction_type_error( errMsg );
      }
    }
    break;
  case c32:
    find_iter = processIndicesInclusive32.find( _id );
    if ( find_iter != processIndicesInclusive32.end() )
    {
      return find_iter->second;
    }
    else
    {
      processIndices32.find( _id );
      if ( find_iter != processIndices32.end() )
      {
        return find_iter->second;
      }
      else
      {
        std::string errMsg = "Determination of process index for 3->2 process failed. Unrecoverable error.";
        throw eInteraction_type_error( errMsg );
      }
    }
    break;
  default:
    std::string errMsg = "Determination of process index failed. Bad generic type. Unrecoverable error.";
    throw eInteraction_type_error( errMsg );
    break;
  }

}



/**
 * Look up the process denoted by _id in the std::maps and return the appropriate index number.
 *
 * @param[in] _index index that needs to be looked up
 * @param[in] _genType generic collision type the specific process belongs to
 * @return process id to be looked up (231, 32800 etc.)
 */
int interactionType::getProcessTypeFromIndex( const int _index, GENERIC_COLL_TYPE _genType )
{
  switch ( _genType )
  {
  case c22:
    if ( _index >= 0 && _index < indexProcesses22.size() )
    {
      return indexProcesses22[_index];
    }
    else
    {
      std::string errMsg = "Resolution of index for 2->2 process failed. Unrecoverable error.";
      throw eInteraction_type_error( errMsg );
    }
    break;
  case c23:
    if ( _index >= 0 && _index < indexProcesses23.size() )
    {
      return indexProcesses23[_index];
    }
    else
    {
      std::string errMsg = "Resolution of index for 2->3 process failed. Unrecoverable error.";
      throw eInteraction_type_error( errMsg );
    }
    break;
  case c32:
    if ( _index >= 0 && _index < indexProcesses32.size() )
    {
      return indexProcesses32[_index];
    }
    else
    {
      std::string errMsg = "Resolution of index for 3->2 process failed. Unrecoverable error.";
      throw eInteraction_type_error( errMsg );
    }
    break;
  default:
    std::string errMsg = "Resolution of index failed. Bad generic type. Unrecoverable error.";
    throw eInteraction_type_error( errMsg );
    break;
  }

}



/**
* Look up the inclusive process denoted by _id in the std::maps and return the appropriate index number.
*
* @param[in] _index index that needs to be looked up
* @param[in] _genType generic collision type the specific process belongs to
* @return process id to be looked up (9231, 9223 etc.)
*/
int interactionType::getInclusiveProcessTypeFromIndex( const int _index, GENERIC_COLL_TYPE _genType )
{
  switch ( _genType )
  {
  case c22:
    if ( _index >= 0 && _index < indexProcessesInclusive22.size() )
    {
      return indexProcessesInclusive22[_index];
    }
    else
    {
      std::string errMsg = "Resolution of index for 2->2 process failed. Unrecoverable error.";
      throw eInteraction_type_error( errMsg );
    }
    break;
  case c23:
    if ( _index >= 0 && _index < indexProcessesInclusive23.size() )
    {
      return indexProcessesInclusive23[_index];
    }
    else
    {
      std::string errMsg = "Resolution of index for 2->3 process failed. Unrecoverable error.";
      throw eInteraction_type_error( errMsg );
    }
    break;
  case c32:
    if ( _index >= 0 && _index < indexProcessesInclusive32.size() )
    {
      return indexProcessesInclusive32[_index];
    }
    else
    {
      std::string errMsg = "Resolution of index for 3->2 process failed. Unrecoverable error.";
      throw eInteraction_type_error( errMsg );
    }
    break;
  default:
    std::string errMsg = "Resolution of index failed. Bad generic type. Unrecoverable error.";
    throw eInteraction_type_error( errMsg );
    break;
  }

}




/**
 * Given the flavors of the incoming particles, return the process type ID associated with this inclusive process (for 2->X)
 *
 * @param[in] _F1 Flavor of incoming particle 1
 * @param[in] _F2 Flavor of incoming particle 2
 * @param[in] genericType generic collision type (c22, c23, c32) the specific process belongs to
 * @return inclusive process type
 */
int interactionType::getInclusiveProcessType( const FLAVOR_TYPE _F1, const FLAVOR_TYPE _F2, const GENERIC_COLL_TYPE genericType )
{
  int type = 0;

  // sort F1 and F2 such that comparisons below are easier
  // convert to unsigned int to allow for some math that speeds up the decisions
  unsigned int F1 = std::min( static_cast<unsigned int>( _F1 ), static_cast<unsigned int>( _F2 ) );
  unsigned int F2 = std::max( static_cast<unsigned int>( _F1 ), static_cast<unsigned int>( _F2 ) );

  if (( F1 + F2 ) == 0 ) // gg -> X
  {
    type = 9221;
  }
  else if (( F1 * F2 ) == 0 && ( F2 % 2 ) == 0 )  // gqbar -> X
  {
    type = 9223;
  }
  else if (( F1 * F2 ) == 0 ) // gq -> X
  {
    type = 9222;
  }
  else if (( F1 % 2 ) == 1 && ( F2 % 2 ) == 1 )  // qq -> X
  {
    if ( F1 == F2 )
    {
      type = 9225;
    }
    else
    {
      type = 9227;  //q -> q'
    }
  }
  else if (( F1 % 2 ) == 0 && ( F2 % 2 ) == 0 )  // qbarqbar -> X
  {
    if ( F1 == F2 )
    {
      type = 9226;
    }
    else
    {
      type = 9229;
    }
  }
  else if ((( F1 % 2 ) == 1 && ( F2 % 2 ) == 0 ) && ( (F2 - F1) == 1 ) )  // qqbar -> X
  {
    type = 9224;
  }
  else if ( (( F1 % 2 ) == 1 && ( F2 % 2 ) == 0 ) ||  (( F2 % 2 ) == 1 && ( F1 % 2 ) == 0 ) )  //qqbar' -> X
  {
    type = 9228;
  }
  else
  {
    std::string errMsg = "Inclusive process type could not be determined. Unrecoverable error.";
    throw eInteraction_type_error( errMsg );
  }

  //fix me: ugly!
  if ( genericType == c23 )
  {
    type += 10;  // 9221 becomes 9231 etc.
  }

  return type;
}



/**
* Given the flavors of the incoming particles, return the process type ID associated with this inclusive process (for 3->2)
*
* @param[in] _F1 Flavor of incoming particle 1
* @param[in] _F2 Flavor of incoming particle 2
* @param[in] _F3 Flavor of incoming particle 3
* @param[in] genericType generic collision type (c22, c23, c32) the specific process belongs to
* @return inclusive process type
*/
int interactionType::getInclusiveProcessType( const FLAVOR_TYPE _F1, const FLAVOR_TYPE _F2, const FLAVOR_TYPE _F3, const GENERIC_COLL_TYPE genericType )
{
  int type = 0;

  FLAVOR_TYPE __F1, __F2;
  if ( _F1 == gluon )
  {
    __F1 = _F2;
    __F2 = _F3;
  }
  else if ( _F2 == gluon )
  {
    __F1 = _F1;
    __F2 = _F3;
  }
  else if ( _F3 == gluon )
  {
    __F1 = _F1;
    __F2 = _F2;
  }
  else
  {
    std::string errMsg = "No gluon in 3->2 process. Unrecoverable error.";
    throw eInteraction_type_error( errMsg );
  }

  
  // sort F1 and F2 such that comparisons below are easier
  // convert to unsigned int to allow for some math that speeds up the decisions
  unsigned int F1 = std::min( static_cast<unsigned int>( __F1 ), static_cast<unsigned int>( __F2 ) );
  unsigned int F2 = std::max( static_cast<unsigned int>( __F1 ), static_cast<unsigned int>( __F2 ) );


  if (( F1 + F2 ) == 0 ) // ggg -> X
  {
    type = 9321;
  }
  else if (( F1 * F2 ) == 0 && ( F2 % 2 ) == 0 )  // gqbarg -> X
  {
    type = 9323;
  }
  else if (( F1 * F2 ) == 0 ) // gqg -> X
  {
    type = 9322;
  }
  else if (( F1 % 2 ) == 1 && ( F2 % 2 ) == 1 )  // qqg -> X
  {
    if ( F1 == F2 )
    {
      type = 9325;
    }
    else
    {
      type = 9327;
    }
  }
  else if (( F1 % 2 ) == 0 && ( F2 % 2 ) == 0 )  // qbarqbarg -> X
  {
    if ( F1 == F2 )
    {
      type = 9326;
    }
    else
    {
      type = 9329;
    }
  }
  else if ((( F1 % 2 ) == 1 && ( F2 % 2 ) == 0 ) && ( (F2 - F1) == 1 ) )  // qqbar+g -> X
  {
    type = 9324;
  }
  else if ( (( F1 % 2 ) == 1 && ( F2 % 2 ) == 0 ) ||  (( F2 % 2 ) == 1 && ( F1 % 2 ) == 0 ) )  //qqbar'+g -> X
  {
    type = 9328;
  }
  else
  {
    std::string errMsg = "Inclusive process type could not be determined (3->2). Unrecoverable error.";
    throw eInteraction_type_error( errMsg );
  }

  return type;
}




GENERIC_COLL_TYPE interactionType::getCollType( const int _id )
{
  int number = _id;
  std::vector<unsigned short int> digits;

  while ( number > 0 )
  {
    digits.push_back( number % 10 );
    number /= 10;
  }
  int _typ = digits[ digits.size() - 1 ] * 10 + digits[ digits.size() - 2 ];

  if ( _typ == 22 )
  {
    return c22;
  }
  else if ( _typ == 23 )
  {
    return c23;
  }
  else if ( _typ == 32 )
  {
    return c32;
  }
}


double interactionType::getInvolvedInitial( const FLAVOR_TYPE _F, const int _id )
{
  if ( _F == gluon )
  {
    switch ( _id )
    {
      // -------- 2 -> 2 --------
    case 221 :
      return 2;
      break;
    case 222 :
      return 2;
      break;
    case 223 :
      return 1;
      break;
    case 2230 :
      return 1;
      break;
      // -------- 2 -> 3 --------
    case 231 :
      return 2;
      break;
    case 232 :
      return 2;
      break;
    case 233 :
      return 1;
      break;
    case 2330 :
      return 1;
      break;
      // -------- 3 -> 2 --------
    case 321 :
      return 3;
      break;
    case 322 :
      return 3;
      break;
    case 323 :
      return 2;
      break;
    case 3230 :
      return 2;
      break;
    case 324 :
      return 1;
      break;
    case 325 :
      return 1;
      break;
    case 326:
      return 1;
      break;
    case 327 :
      return 1;
      break;
    case 3270 :
      return 1;
      break;
    case 328 :
      return 1;
      break;
    case 3280 :
      return 1;
      break;
    case 32800 :
      return 1;
      break;
    default :
      return 0;
      break;
    }
  }
  else if ( _F == quark || _F == up || _F == down || _F == charm )
  {
    switch ( _id )
    {
      // -------- 2 -> 2 --------
    case 223 :
      return 1;
      break;
    case 224 :
      return 1;
      break;
    case 225 :
      return 1;
      break;
    case 226 :
      return 1;
      break;
    case 227 :
      return 2;
      break;
    case 228 :
      return 2;
      break;
    case 2280:
      return 1;
      break;
      // -------- 2 -> 3 --------
    case 233 :
      return 1;
      break;
    case 234 :
      return 1;
      break;
    case 235 :
      return 1;
      break;
    case 236 :
      return 1;
      break;
    case 237 :
      return 2;
      break;
    case 238 :
      return 2;
      break;
    case 2380:
      return 1;
      break;
      // -------- 3 -> 2 --------
    case 323 :
      return 1;
      break;
    case 324 :
      return 1;
      break;
    case 325 :
      return 1;
      break;
    case 326 :
      return 1;
      break;
    case 327 :
      return 2;
      break;
    case 328 :
      return 2;
      break;
    case 3280:
      return 1;
      break;
    default :
      return 0;
      break;
    }
  }
  else if ( _F == anti_quark || _F == anti_up || _F == anti_down || _F == anti_charm )
  {
    switch ( _id )
    {
      // -------- 2 -> 2 --------
    case 2230 :
      return 1;
      break;
    case 224 :
      return 1;
      break;
    case 225 :
      return 1;
      break;
    case 226 :
      return 1;
      break;
    case 2270 :
      return 2;
      break;
    case 2280:
      return 1;
      break;
    case 22800:
      return 2;
      break;
      // -------- 2 -> 3 --------
    case 2330 :
      return 1;
      break;
    case 234 :
      return 1;
      break;
    case 235 :
      return 1;
      break;
    case 236 :
      return 1;
      break;
    case 2370 :
      return 2;
      break;
    case 2380:
      return 1;
      break;
    case 23800:
      return 2;
      break;
      // -------- 3 -> 2 --------
    case 3230 :
      return 1;
      break;
    case 324 :
      return 1;
      break;
    case 325 :
      return 1;
      break;
    case 326 :
      return 1;
      break;
    case 3270 :
      return 2;
      break;
    case 3280:
      return 1;
      break;
    case 32800:
      return 2;
      break;
    default :
      return 0;
      break;
    }
  }
  else
  {
    return -1;
  }

}



double interactionType::getInvolvedInitialInclusive( const FLAVOR_TYPE _F, const int _id )
{
  std::map<int, int>::const_iterator find_iter;

  switch ( _F )
  {
  case gluon:
    find_iter = involvedInitialGluons.find( _id );
    if ( find_iter != involvedInitialGluons.end() )
    {
      return find_iter->second;
    }
    else
    {
      std::string errMsg = "Determination of process index for 2->2 process failed. Unrecoverable error.";
      throw eInteraction_type_error( errMsg );
    }
    break;
  case quark:
    find_iter = involvedInitialQuarks.find( _id );
    if ( find_iter != involvedInitialQuarks.end() )
    {
      return find_iter->second;
    }
    else
    {
      std::string errMsg = "Determination of process index for 2->3 process failed. Unrecoverable error.";
      throw eInteraction_type_error( errMsg );
    }
    break;
  case anti_quark:
    find_iter = involvedInitialAntiQuarks.find( _id );
    if ( find_iter != involvedInitialAntiQuarks.end() )
    {
      return find_iter->second;
    }
    else
    {
      std::string errMsg = "Determination of process index for 3->2 process failed. Unrecoverable error.";
      throw eInteraction_type_error( errMsg );
    }
    break;
  default:
    std::string errMsg = "Determination of process index failed. Bad generic type. Unrecoverable error.";
    throw eInteraction_type_error( errMsg );
    break;
  }

}





std::vector<int> interactionType::getInteractionTypes( const FLAVOR_TYPE _F )
{
  std::vector<int> _types;
  _types.reserve( 20 );

  if ( _F == gluon )
  {
    int mytypes[] = { 221, 222, 223, 2230, 231, 232, 233, 2330, 321, 322, 323, 3230, 324, 325, 326, 327, 3270, 328, 3280, 32800 };
    _types.assign( mytypes, mytypes + sizeof( mytypes ) / sizeof( int ) );
  }
  else if ( _F == quark || _F == up || _F == down || _F == charm )
  {
    int mytypes[] = {  223, 224, 225, 226, 227, 228, 2280, 233, 234, 235, 236, 237, 238, 2380, 323, 324, 325, 326, 327, 328, 3280 };
    _types.assign( mytypes, mytypes + sizeof( mytypes ) / sizeof( int ) );
  }
  else if ( _F == anti_quark || _F == anti_up || _F == anti_down || _F == anti_charm )
  {
    int mytypes[] = {  2230, 224, 225, 226, 2270, 2280, 22800, 2330, 234, 235, 236, 2370, 2380, 23800, 3230, 324, 325, 326, 3270, 3280, 32800 };
    _types.assign( mytypes, mytypes + sizeof( mytypes ) / sizeof( int ) );
  }

  return _types;
}




// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;
