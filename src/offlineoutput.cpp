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
std::string  offlineDataSimulationParameters::filenameIdentifier = "simulationParameters";


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
BOOST_CLASS_EXPORT_IMPLEMENT( offlineDataSimulationParameters )
