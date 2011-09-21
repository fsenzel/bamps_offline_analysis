//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/restructureOutput/src/offlineoutput.cpp $
//$LastChangedDate: 2011-08-25 16:15:11 +0200 (Thu, 25 Aug 2011) $
//$LastChangedRevision: 104 $
//$LastChangedBy: fochler $
//---------------------------------------------
//---------------------------------------------

#include <typeinfo>
#include <iostream>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "offlineoutput.h"


std::string  offlineDataGeneric::filenameIdentifier = "oflineDataGeneric";
std::string  offlineDataEventType::filenameIdentifier = "eventList";
std::string  offlineDataInteraction22::filenameIdentifier = "interaction22";
std::string  offlineDataInteraction23::filenameIdentifier = "interaction23";
std::string  offlineDataInteraction32::filenameIdentifier = "interaction32";
std::string  offlineDataInteractionElastic::filenameIdentifier = "interactionElastic";
std::string  offlineDataCellConfiguration::filenameIdentifier = "cellConfiguration";
std::string  offlineDataParticleIdSwap::filenameIdentifier = "particleIdSwap";
std::string  offlineDataParticleNumbers::filenameIdentifier = "particleNumbers";
std::string  offlineDataCollisionNumbers::filenameIdentifier = "collisionNumbers";
std::string  offlineDataInteractionRates::filenameIdentifier = "interactionRates";


offlineOutputInterface::~offlineOutputInterface()
{
  collectAndOutputEventList();
  
  archiveMap.clear(); // due to the usage of shared pointers this is sufficient to take care of all cleanup
}



void offlineOutputInterface::collectAndOutputEventList()
{
  for ( int i = 0; i < eventList.size(); i++ )
  {
    offlineDataEventType collectEventType( eventList[i] );
    this->submitOfflineDataForOutput( &collectEventType );
  }
}


void offlineOutputInterface::submitOfflineDataForOutput( const offlineDataGeneric*const _data )
{
  tPointerToFilestream stream;
  tPointerToArchive archive;
 
  //check whether an output archive for the given type of data has already been created, re-use if yes, create if no
  tArchiveMap::iterator it = archiveMap.find( type_index(typeid(*_data)) );
  if( it != archiveMap.end() )
  {
//     stream = it->second;
    archive = it->second;
  }
  else
  {
    stream.reset( new boost::filesystem::ofstream() );
    std::string filename = outputDirectory.string() + "/" + filenamePrefix + additionalFilenameTag + _data->getFilenameIdentifier() + filenameSuffix;
    
    #if BINARY_OFFLINE_OUTPUT > 0
      stream->open( filename.c_str(), std::ios::trunc | std::ios::binary );
    #else
      stream->open( filename.c_str(), std::ios::trunc );
    #endif
    
    outputMap.insert(std::make_pair( type_index(typeid(*_data)), stream));
    archive.reset( new tArchive( *stream ) );
    archiveMap.insert(std::make_pair( type_index(typeid(*_data)), archive));
  }
  
  // the actual output
  (*archive) & _data;
}


template<class T>
boost::shared_ptr<T> offlineOutputInterface::readOfflineDataFromArchive()
{
  tPointerToInputFilestream stream;
  tPointerToInputArchive archive;
  
//   check whether an output archive for the given type of data has already been created, re-use if yes, create if no
  tInputArchiveMap::iterator it = inputArchiveMap.find( type_index(typeid(T)) );
  if( it != inputArchiveMap.end() )
  {
    tInputFileStreamMap::iterator itFile = inputMap.find( type_index(typeid(T)) );
    if( itFile != inputMap.end() )
    {
      stream = itFile->second;    
    }
    archive = it->second;
  }
  else
  {
    stream.reset( new boost::filesystem::ifstream() );
    std::string filename = outputDirectory.string() + "/" + filenamePrefix + additionalFilenameTag + T::filenameIdentifier + filenameSuffix;
    
    #if BINARY_OFFLINE_OUTPUT > 0
      stream->open( filename.c_str(), std::ios::binary );
    #else
      stream->open( filename.c_str() );
    #endif
    
    inputMap.insert(std::make_pair( type_index(typeid(T)), stream));
    archive.reset( new tInputArchive( *stream ) );
    inputArchiveMap.insert(std::make_pair( type_index(typeid(T)), archive));
    
    lastStreamPositionMap.insert( std::make_pair( type_index(typeid(T)), 0 ));
  }
  
  // the actual reading
  tStreamPositionMap::iterator itLastPos = lastStreamPositionMap.find( type_index(typeid(T)) );
  if( itLastPos != lastStreamPositionMap.end() )
  {
    itLastPos->second = stream->tellg();    
  }
  offlineDataGeneric* _ptr = 0;
  (*archive) & _ptr;
  
  T* ptrToDerived = dynamic_cast<T*>(_ptr);
  if ( ptrToDerived == 0 )
  {
    std::string errMsg = "Bad cast. Attempted dynamic_cast from ";
    errMsg += type_index(typeid(offlineDataGeneric)).name();
    errMsg += " to ";
    errMsg += type_index(typeid(T)).name();
    throw eOfflineOutput_error( errMsg );
  }
  
  boost::shared_ptr<T> shrptr( ptrToDerived );
  return shrptr;
}
         
template boost::shared_ptr<offlineDataInteraction22> offlineOutputInterface::readOfflineDataFromArchive<offlineDataInteraction22>();
template boost::shared_ptr<offlineDataInteraction23> offlineOutputInterface::readOfflineDataFromArchive<offlineDataInteraction23>();
template boost::shared_ptr<offlineDataInteraction32> offlineOutputInterface::readOfflineDataFromArchive<offlineDataInteraction32>();
template boost::shared_ptr<offlineDataInteractionElastic> offlineOutputInterface::readOfflineDataFromArchive<offlineDataInteractionElastic>();
template boost::shared_ptr<offlineDataCellConfiguration> offlineOutputInterface::readOfflineDataFromArchive<offlineDataCellConfiguration>();
template boost::shared_ptr<offlineDataInteractionRates> offlineOutputInterface::readOfflineDataFromArchive<offlineDataInteractionRates>();
template boost::shared_ptr<offlineDataEventType> offlineOutputInterface::readOfflineDataFromArchive<offlineDataEventType>();
template boost::shared_ptr<offlineDataParticleIdSwap> offlineOutputInterface::readOfflineDataFromArchive<offlineDataParticleIdSwap>();



template<class T>
void offlineOutputInterface::undoLastReadOperation()
{
  tPointerToInputFilestream stream;
  
  std::streampos lastStreamPosition;
  tStreamPositionMap::iterator itLastPos = lastStreamPositionMap.find( type_index(typeid(T)) );
  if( itLastPos != lastStreamPositionMap.end() )
  {
    lastStreamPosition = itLastPos->second;    
  }
  else
  {
    std::string errMsg = "Attempting undo operation on stream that has not been initialized yet";
    throw eOfflineOutput_error( errMsg );
  }
  
//   check whether an output archive for the given type of data has already been created, re-use if yes, create if no
  tInputFileStreamMap::iterator it = inputMap.find( type_index(typeid(T)) );
  if( it != inputMap.end() )
  {
    stream = it->second;
    stream->seekg( lastStreamPosition );
  }
}

template void offlineOutputInterface::undoLastReadOperation<offlineDataInteraction22>();
template void offlineOutputInterface::undoLastReadOperation<offlineDataInteraction23>();
template void offlineOutputInterface::undoLastReadOperation<offlineDataInteraction32>();
template void offlineOutputInterface::undoLastReadOperation<offlineDataInteractionElastic>();
template void offlineOutputInterface::undoLastReadOperation<offlineDataInteractionRates>();
template void offlineOutputInterface::undoLastReadOperation<offlineDataCellConfiguration>();
template void offlineOutputInterface::undoLastReadOperation<offlineDataParticleIdSwap>();
template void offlineOutputInterface::undoLastReadOperation<offlineDataEventType>();



void offlineOutputInterface::registerOfflineDataForTemporaryStorage(const offlineEventType _eventType, tConstPointerToOfflineData _data)
{
  this->registerOfflineDataForTemporaryStorage( tConstPointerToOfflineData( new offlineDataEventType(_eventType) ) );
  this->registerOfflineDataForTemporaryStorage( _data );
}


void offlineOutputInterface::outputAndResetTemporaryStorage()
{
  for ( tTemporaryDataStorage::const_iterator it = temporaryStorage.begin(); it != temporaryStorage.end(); it++ )
  {
    this->submitOfflineDataForOutput( (*it).get() ); //submit a pointer to the object stored at the iterator it
  }
  
  temporaryStorage.clear();
}


BOOST_CLASS_EXPORT_IMPLEMENT( offlineDataInteraction22 )
BOOST_CLASS_EXPORT_IMPLEMENT( offlineDataInteraction23 )
BOOST_CLASS_EXPORT_IMPLEMENT( offlineDataInteraction32 )
BOOST_CLASS_EXPORT_IMPLEMENT( offlineDataInteractionElastic )
BOOST_CLASS_EXPORT_IMPLEMENT( offlineDataCellConfiguration )
BOOST_CLASS_EXPORT_IMPLEMENT( offlineDataParticleIdSwap )
BOOST_CLASS_EXPORT_IMPLEMENT( offlineDataParticleNumbers )
BOOST_CLASS_EXPORT_IMPLEMENT( offlineDataCollisionNumbers )
BOOST_CLASS_EXPORT_IMPLEMENT( offlineDataInteractionRates )
BOOST_CLASS_EXPORT_IMPLEMENT( offlineDataEventType )
