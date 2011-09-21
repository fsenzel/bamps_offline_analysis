//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/restructureOutput/src/offlineoutput.h $
//$LastChangedDate: 2011-08-25 17:42:35 +0200 (Thu, 25 Aug 2011) $
//$LastChangedRevision: 107 $
//$LastChangedBy: fochler $
//---------------------------------------------
//---------------------------------------------


#ifndef OFFLINE_OUTPUT_H
#define OFFLINE_OUTPUT_H

// #define BINARY_OFFLINE_OUTPUT 0
#define BINARY_OFFLINE_OUTPUT 1

#define DEFAULT_EVENT_NUMBER_ESTIMATE 1000000
#define DEFAULT_TEMPORARY_EVENT_NUMBER_ESTIMATE 100

#include <iostream>
#include <stdexcept>
#include <string>
#include <map>
#include <vector>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/level.hpp>
#include <boost/smart_ptr.hpp>

// #define BOOST_FILESYSTEM_VERSION 3
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "coordinateBins.h"
#include "typeindex.h"
#include "particle.h"


enum offlineEventType { event_interaction22, event_interaction23, event_interaction32, event_interactionElastic, 
                        event_particleIdSwap, event_newTimestep, event_endOfCascade, event_dummy = 99 };

class offlineDataGeneric
{
  public:
    offlineDataGeneric() {};
    ~offlineDataGeneric() {};
    
    virtual std::string getFilenameIdentifier() const { return this->filenameIdentifier; }
    virtual size_t getSize() const = 0;
    
    static std::string filenameIdentifier;
    
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
    }
};


typedef boost::shared_ptr<offlineDataGeneric> tPointerToOfflineData;
typedef boost::shared_ptr<const offlineDataGeneric> tConstPointerToOfflineData;

class offlineOutputInterface
{
  public:
    offlineOutputInterface( const char* _outputDir, const bool _active = true ) : 
      outputDirectory( _outputDir ), 
      filenamePrefix( "offline_" ), 
      filenameSuffix( ".dat" ),
      additionalFilenameTag( "" )
      { 
        if ( _active )
        {
          checkAndCreateOutputDirectory( outputDirectory );
          eventList.reserve( DEFAULT_EVENT_NUMBER_ESTIMATE ); 
          temporaryStorage.reserve( DEFAULT_TEMPORARY_EVENT_NUMBER_ESTIMATE );
        }
      }
    offlineOutputInterface( const char* _outputDir, const int eventNumberEstimate, const bool _active = true ) : 
      outputDirectory( _outputDir ), 
      filenamePrefix( "offline_" ), 
      filenameSuffix( ".dat" ),
      additionalFilenameTag( "" )
      { 
        if ( _active )
        {
          checkAndCreateOutputDirectory( outputDirectory );
          eventList.reserve( eventNumberEstimate );
          temporaryStorage.reserve( DEFAULT_TEMPORARY_EVENT_NUMBER_ESTIMATE );
        }
      }
    ~offlineOutputInterface();
    
    void setAdditionalFilenameTag( std::string _tag ) { additionalFilenameTag = _tag + "_"; }
    
    void submitOfflineDataForOutput( const offlineDataGeneric* const _data );
    void submitOfflineDataForOutput( const offlineEventType _eventType )
    { 
      eventList.push_back( _eventType );
    } 
    
    void registerOfflineDataForTemporaryStorage( tConstPointerToOfflineData _data )
    {
      temporaryStorage.push_back( _data );
    }
    void registerOfflineDataForTemporaryStorage( const offlineEventType _eventType, tConstPointerToOfflineData _data );
    void resetTemporaryStorage() { temporaryStorage.clear(); }
    void outputAndResetTemporaryStorage();
    
    /** read the stuff */
    template<class T>
    boost::shared_ptr<T> readOfflineDataFromArchive();
    
    /** provide ability to undo read operation */
    template<class T>
    void undoLastReadOperation();
    
  private:
    void checkAndCreateOutputDirectory( boost::filesystem::path& _dir );
    
    typedef boost::shared_ptr<boost::filesystem::ofstream> tPointerToFilestream;
    typedef std::map< type_index, tPointerToFilestream > tFileStreamMap;
    tFileStreamMap outputMap;

    #if BINARY_OFFLINE_OUTPUT > 0 
      typedef boost::archive::binary_oarchive tArchive;
    #else
      typedef boost::archive::text_oarchive tArchive;
    #endif
    typedef boost::shared_ptr<tArchive> tPointerToArchive;
    typedef std::map< type_index, tPointerToArchive > tArchiveMap;
    tArchiveMap archiveMap;
    
    
    /** read the stuff */
    typedef boost::shared_ptr<boost::filesystem::ifstream> tPointerToInputFilestream;
    typedef std::map< type_index, tPointerToInputFilestream > tInputFileStreamMap;
    tInputFileStreamMap inputMap;

    #if BINARY_OFFLINE_OUTPUT > 0 
      typedef boost::archive::binary_iarchive tInputArchive;
    #else
      typedef boost::archive::text_iarchive tInputArchive;
    #endif
    typedef boost::shared_ptr<tInputArchive> tPointerToInputArchive;
    typedef std::map< type_index, tPointerToInputArchive > tInputArchiveMap;
    tInputArchiveMap inputArchiveMap;
    /** read the stuff */
    
    /** save the last read position */
    typedef std::map< type_index, std::streampos > tStreamPositionMap;
    tStreamPositionMap lastStreamPositionMap;
    /** save the last read position */
    
    typedef std::vector<tConstPointerToOfflineData> tTemporaryDataStorage;
    tTemporaryDataStorage temporaryStorage;
    
    std::vector<offlineEventType> eventList;
    void collectAndOutputEventList();    
    
    boost::filesystem::path outputDirectory;
    std::string filenamePrefix;
    std::string filenameSuffix;
    std::string additionalFilenameTag;
};


inline void offlineOutputInterface::checkAndCreateOutputDirectory(boost::filesystem::path& _dir)
{
  if ( boost::filesystem::exists( _dir ) )
  {
    if ( boost::filesystem::is_directory( _dir ) )
    {
      return;
    }
    else
    {
      boost::filesystem::path renamePath( _dir.string() + ".backup" );
      std::cout << "File with name " << _dir.string() << " blocks the creation of an output folder for offline reconstruction." << std::endl;
      std::cout << "It is renamed to " << renamePath.string() << std::endl;
      boost::filesystem::rename( _dir, renamePath );
      boost::filesystem::create_directory( _dir );       
    }
  }
  else
  {
    std::cout << "Creating output folder for offline reconstruction data: " << _dir.string() << std::endl;
    boost::filesystem::create_directory( _dir );    
  }
}



class offlineDataEventType : public offlineDataGeneric
{
  public:
    offlineDataEventType() : offlineDataGeneric(), event( event_interaction22 ) {};
    offlineDataEventType( const offlineEventType _event ) : offlineDataGeneric(), event(_event) {};
    ~offlineDataEventType() {};
    
    size_t getSize() const { return sizeof( offlineDataEventType ); }
    std::string getFilenameIdentifier() const { return this->filenameIdentifier; }
    static std::string filenameIdentifier;
    
    /** @brief event type */
    offlineEventType event;
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & boost::serialization::base_object<offlineDataGeneric>(*this);
      ar & event;
    }
};


class offlineDataInteraction22 : public offlineDataGeneric
{
  public:
    offlineDataInteraction22() : offlineDataGeneric(),
                              iscat(-1), jscat(-1), time(-1), pix(-999), piy(-999), piz(-999), pjx(-999), pjy(-999),
                              pjz(-999), F1(gluon), F2(gluon) {};
    offlineDataInteraction22( const int _iscat, const int _jscat, const double _time, const double _pix, const double _piy,
                              const double _piz, const double _pjx, const double _pjy, const double _pjz,
                              const FLAVOR_TYPE _F1, const FLAVOR_TYPE _F2 ) : offlineDataGeneric(),
                              iscat(_iscat), jscat(_jscat), time(_time), pix(_pix), piy(_piy), piz(_piz), pjx(_pjx), pjy(_pjy),
                              pjz(_pjz), F1(_F1), F2(_F2) {};
    ~offlineDataInteraction22() {};
    
    size_t getSize() const { return sizeof( offlineDataInteraction22 ); }
    std::string getFilenameIdentifier() const { return this->filenameIdentifier; }
    static std::string filenameIdentifier;
    
    /** @brief Current ID of particle 1 */
    int iscat;
    /** @brief Current ID of particle 2 */
    int jscat;
    /** @brief Current time of the collision */
    double time;
    /** @brief Momentum (p_x) of particle 1 AFTER the collison */
    double pix;
    /** @brief Momentum (p_y) of particle 1 AFTER the collison */
    double piy;
    /** @brief Momentum (p_z) of particle 1 AFTER the collison */
    double piz;
    /** @brief Momentum (p_x) of particle 2 AFTER the collison */
    double pjx;
    /** @brief Momentum (p_y) of particle 2 AFTER the collison */
    double pjy;
    /** @brief Momentum (p_z) of particle 2 AFTER the collison */
    double pjz;
    /** @brief Flavor of particle 1 AFTER the collison */
    FLAVOR_TYPE F1; 
    /** @brief Flavor of particle 1 AFTER the collison */
    FLAVOR_TYPE F2;
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & boost::serialization::base_object<offlineDataGeneric>(*this);
      ar & iscat;
      ar & jscat;
      ar & time;
      ar & pix;
      ar & piy;
      ar & piz;
      ar & pjx;
      ar & pjy;
      ar & pjz;
      ar & F1;
      ar & F2;
    }
};


class offlineDataInteraction23 : public offlineDataGeneric
{
  public:
    offlineDataInteraction23() : offlineDataGeneric(),
                              iscat(-1), jscat(-1), newp(-1), time(-1), pix(999), piy(999), piz(999),
                              pjx(999), pjy(999), pjz(999), newx(999), newy(999), newz(999),
                              newpx(999), newpy(999), newpz(999), F1(gluon), F2(gluon), F3(gluon) {};
    offlineDataInteraction23( const int _iscat, const int _jscat, const int _newp, const double _time, const double _pix, const double _piy,
                              const double _piz, const double _pjx, const double _pjy, const double _pjz,
                              const double _newx, const double _newy, const double _newz,
                              const double _newpx, const double _newpy, const double _newpz,
                              const FLAVOR_TYPE _F1, const FLAVOR_TYPE _F2, const FLAVOR_TYPE _F3 ) : offlineDataGeneric(),
                              iscat(_iscat), jscat(_jscat), newp(_newp), time(_time), pix(_pix), piy(_piy), piz(_piz),
                              pjx(_pjx), pjy(_pjy), pjz(_pjz), newx(_newx), newy(_newy), newz(_newz),
                              newpx(_newpx), newpy(_newpy), newpz(_newpz), F1(_F1), F2(_F2), F3(_F3) {};
    ~offlineDataInteraction23() {};
    
    size_t getSize() const { return sizeof( offlineDataInteraction23 ); }
    std::string getFilenameIdentifier() const { return this->filenameIdentifier; }
    static std::string filenameIdentifier;
    
    /** @brief Current ID of particle 1 */
    int iscat;
    /** @brief Current ID of particle 2 */
    int jscat;
    /** @brief ID of newly created particle */
    int newp;
    /** @brief Current time of the collision */
    double time;
    /** @brief Momentum (p_x) of particle 1 AFTER the collison */
    double pix;
    /** @brief Momentum (p_y) of particle 1 AFTER the collison */
    double piy;
    /** @brief Momentum (p_z) of particle 1 AFTER the collison */
    double piz;
    /** @brief Momentum (p_x) of particle 2 AFTER the collison */
    double pjx;
    /** @brief Momentum (p_y) of particle 2 AFTER the collison */
    double pjy;
    /** @brief Momentum (p_z) of particle 2 AFTER the collison */
    double pjz;
    /** @brief Position (x) of new particle (particle 3) */
    double newx;
    /** @brief Position (y) of new particle (particle 3) */
    double newy;
    /** @brief Position (z) of new particle (particle 3) */
    double newz;
    /** @brief Momentum (p_x) of new particle (particle 3) */
    double newpx;
    /** @brief Momentum (p_y) of new particle (particle 3) */
    double newpy;
    /** @brief Momentum (p_z) of new particle (particle 3) */
    double newpz;
    /** @brief Flavor of particle 1 AFTER the collison */
    int F1; 
    /** @brief Flavor of particle 1 AFTER the collison */
    int F2;
    /** @brief Flavor of particle 3 AFTER the collison */
    int F3;
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & boost::serialization::base_object<offlineDataGeneric>(*this);
      ar & iscat;
      ar & jscat;
      ar & newp;
      ar & time;
      ar & pix;
      ar & piy;
      ar & piz;
      ar & pjx;
      ar & pjy;
      ar & pjz;
      ar & newx;
      ar & newy;
      ar & newz;
      ar & newpx;
      ar & newpy;
      ar & newpz;
      ar & F1;
      ar & F2;
      ar & F3;
    }
};


class offlineDataInteraction32 : public offlineDataGeneric
{
  public:
    offlineDataInteraction32() : offlineDataGeneric(),
                              iscat(-1), jscat(-1), dead(-1), time(-1), pix(999), piy(999), piz(999), pjx(999), pjy(999),
                              pjz(999), F1(gluon), F2(gluon) {};
    offlineDataInteraction32( const int _iscat, const int _jscat, const int _dead, const double _time, const double _pix, const double _piy,
                              const double _piz, const double _pjx, const double _pjy, const double _pjz,
                              const FLAVOR_TYPE _F1, const FLAVOR_TYPE _F2 ) : offlineDataGeneric(),
                              iscat(_iscat), jscat(_jscat), dead(_dead), time(_time), pix(_pix), piy(_piy), piz(_piz), pjx(_pjx), pjy(_pjy),
                              pjz(_pjz), F1(_F1), F2(_F2) {};
    ~offlineDataInteraction32() {};
    
    size_t getSize() const { return sizeof( offlineDataInteraction32 ); }
    std::string getFilenameIdentifier() const { return this->filenameIdentifier; }
    static std::string filenameIdentifier;
    
    /** @brief Current ID of particle 1 */
    int iscat;
    /** @brief Current ID of particle 2 */
    int jscat;
    /** @brief ID of absorbed ("dead") particle */
    int dead;
    /** @brief Current time of the collision */
    double time;
    /** @brief Momentum (p_x) of particle 1 AFTER the collison */
    double pix;
    /** @brief Momentum (p_y) of particle 1 AFTER the collison */
    double piy;
    /** @brief Momentum (p_z) of particle 1 AFTER the collison */
    double piz;
    /** @brief Momentum (p_x) of particle 2 AFTER the collison */
    double pjx;
    /** @brief Momentum (p_y) of particle 2 AFTER the collison */
    double pjy;
    /** @brief Momentum (p_z) of particle 2 AFTER the collison */
    double pjz;
    /** @brief Flavor of particle 1 AFTER the collison */
    int F1; 
    /** @brief Flavor of particle 1 AFTER the collison */
    int F2;
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & boost::serialization::base_object<offlineDataGeneric>(*this);
      ar & iscat;
      ar & jscat;
      ar & dead;
      ar & time;
      ar & pix;
      ar & piy;
      ar & piz;
      ar & pjx;
      ar & pjy;
      ar & pjz;
      ar & F1;
      ar & F2;
    }
};


class offlineDataInteractionElastic : public offlineDataGeneric
{
  public:
    offlineDataInteractionElastic( ) : offlineDataGeneric(),
                                   iscat(-1), jscat(-1), ct_i(-1), ct_j(-1), pix(999), piy(999), piz(999), 
                                   pjx(999), pjy(999), pjz(999) {};
    offlineDataInteractionElastic( const int _iscat, const int _jscat, const double _ct_i, const double _ct_j, 
                                   const double _pix, const double _piy, const double _piz, const double _pjx, 
                                   const double _pjy, const double _pjz ) : offlineDataGeneric(),
                                   iscat(_iscat), jscat(_jscat), ct_i(_ct_i), ct_j(_ct_j), pix(_pix), piy(_piy), piz(_piz), 
                                   pjx(_pjx), pjy(_pjy), pjz(_pjz) {};
    ~offlineDataInteractionElastic() {};
    
    size_t getSize() const { return sizeof( offlineDataInteractionElastic ); }
    std::string getFilenameIdentifier() const { return this->filenameIdentifier; }
    static std::string filenameIdentifier;
    
    /** @brief Current ID of particle 1 */
    int iscat;
    /** @brief Current ID of particle 2 */
    int jscat;
    /** @brief Geometric collision time of particle 1 */
    double ct_i;
    /** @brief Geometric collision time of particle 2 */
    double ct_j;
    /** @brief Momentum (p_x) of particle 1 AFTER the collison */
    double pix;
    /** @brief Momentum (p_y) of particle 1 AFTER the collison */
    double piy;
    /** @brief Momentum (p_z) of particle 1 AFTER the collison */
    double piz;
    /** @brief Momentum (p_x) of particle 2 AFTER the collison */
    double pjx;
    /** @brief Momentum (p_y) of particle 2 AFTER the collison */
    double pjy;
    /** @brief Momentum (p_z) of particle 2 AFTER the collison */
    double pjz;
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & boost::serialization::base_object<offlineDataGeneric>(*this);
      ar & iscat;
      ar & jscat;
      ar & ct_i;
      ar & ct_j;
      ar & pix;
      ar & piy;
      ar & piz;
      ar & pjx;
      ar & pjy;
      ar & pjz;
    }
};


class offlineDataCellConfiguration : public offlineDataGeneric
{
  public:
    offlineDataCellConfiguration( ) : offlineDataGeneric(), timenow(-1), timenext(-1),
                                  randomShiftX(0), randomShiftY(0), randomShiftEta(0) {};
    offlineDataCellConfiguration( const double _timenow, const double _timenext, 
                                  const double _randomshiftX, const double _randomShiftY, const double _randomShiftEta,
                                  const coordinateEtaBins& _etaBins ) : timenow(_timenow), timenext(_timenext),
                                  randomShiftX(_randomshiftX), randomShiftY(_randomShiftY), randomShiftEta(_randomShiftEta),
                                  etaBins(_etaBins), offlineDataGeneric() {};
    ~offlineDataCellConfiguration() {};
    
    size_t getSize() const { return sizeof( offlineDataCellConfiguration ); }
    std::string getFilenameIdentifier() const { return this->filenameIdentifier; }
    static std::string filenameIdentifier;
    
    /** @brief Current simulation time */
    double timenow;
    /** @brief The next time step */
    double timenext;
    /** @brief Random shift of the cell grid in x-direction */
    double randomShiftX;
    /** @brief Random shift of the cell grid in y-direction */
    double randomShiftY;
    /** @brief Random shift of the cell grid in eta-direction (= z-direction) */
    double randomShiftEta;
    /** @brief The current cell configuration in eta-direction */
    coordinateEtaBins etaBins;    
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & boost::serialization::base_object<offlineDataGeneric>(*this);
      ar & timenow;
      ar & timenext;
      ar & randomShiftX;
      ar & randomShiftY;
      ar & randomShiftEta;
      ar & etaBins;
    }
};


class offlineDataParticleIdSwap : public offlineDataGeneric
{
  public:
    offlineDataParticleIdSwap() : removedParticleID(-1),
      replacingParticleID(-1), offlineDataGeneric() {};
    offlineDataParticleIdSwap( const int _removedParticleID, const int _replacingParticleID ) : removedParticleID(_removedParticleID),
      replacingParticleID(_replacingParticleID), offlineDataGeneric() {};
    ~offlineDataParticleIdSwap() {};
    
    size_t getSize() const { return sizeof( offlineDataParticleIdSwap ); }
    std::string getFilenameIdentifier() const { return this->filenameIdentifier; }
    static std::string filenameIdentifier;
    
    /** @brief ID of the removed ("dead") particle */
    int removedParticleID;
    /** @brief ID of the particle that takes the place of the removed particle */
    int replacingParticleID;
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & boost::serialization::base_object<offlineDataGeneric>(*this);
      ar & removedParticleID;
      ar & replacingParticleID;
    }
};


class offlineDataParticleNumbers : public offlineDataGeneric
{
  public:
    offlineDataParticleNumbers( ) : time(-1), dt(-1), nTotalParticles(0),
                                nParticlesInFormationIntialProduction(0),
                                nParticlesInFormationGeometricCollisions(0),
                                nParticlesActive(0), nParticlesInActiveCells(0),
                                edgeCellSizes(0), nFreeParticles( 0 ), offlineDataGeneric() {};
    offlineDataParticleNumbers( const double _time, const double _dt, const int _nTotalParticles,
                                const int _nParticlesInFormationIntialProduction, 
                                const int _nParticlesInFormationGeometricCollisions, const int _nParticlesActive,
                                const int _nParticlesInActiveCells, const std::vector<int>& _edgeCellSizes,
                                const int _nFreeParticles ) : time(_time), dt(_dt), nTotalParticles(_nTotalParticles),
                                nParticlesInFormationIntialProduction(_nParticlesInFormationIntialProduction),
                                nParticlesInFormationGeometricCollisions(_nParticlesInFormationGeometricCollisions),
                                nParticlesActive(_nParticlesActive), nParticlesInActiveCells(_nParticlesInActiveCells),
                                edgeCellSizes(_edgeCellSizes), nFreeParticles( _nFreeParticles ), offlineDataGeneric() {};
    ~offlineDataParticleNumbers() {};
    
    size_t getSize() const { return sizeof( offlineDataParticleNumbers ); }
    std::string getFilenameIdentifier() const { return this->filenameIdentifier; }
    static std::string filenameIdentifier;
    
    /** @brief Current simulation time */
    double time;
    /** @brief Current time step \Delta t */
    double dt;
    /** @brief Current size of the particle vector */
    int nTotalParticles;
//    /** @brief Number of particles that are in formation */
//     int nParticlesInFormation;
    /** @brief Number of particles that are in formation from initial production */
    int nParticlesInFormationIntialProduction;
    /** @brief Number of particles that are in formation due to geometric collisions */
    int nParticlesInFormationGeometricCollisions;
    /** @brief The number of active particles, i.e. particles for which T < timenext */
    int nParticlesActive;
    /** @brief Number of particles in "active" cells, i.e. in cells whose particle content is larger than ::cellcut */
    int nParticlesInActiveCells;
    /** @brief The sizes of the edge cells (= number of particles in these cells) */
    std::vector<int> edgeCellSizes;
    /** @brief Total number of particles that are free because they are a) outside the transversal grid or b) in free cells */
    int nFreeParticles;
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & boost::serialization::base_object<offlineDataGeneric>(*this);
      ar & time;
      ar & dt;
      ar & nTotalParticles;
//       ar & nParticlesInFormation;
      ar & nParticlesInFormationIntialProduction;
      ar & nParticlesInFormationGeometricCollisions;
      ar & nParticlesActive;
      ar & nParticlesInActiveCells;
      ar & edgeCellSizes;
      ar & nFreeParticles;
    }
};


class offlineDataCollisionNumbers : public offlineDataGeneric
{
  public:
    offlineDataCollisionNumbers( ) : time( -1 ),
                                 nCollisions22( 0 ), nCollisions23( 0 ), nCollisions32( 0 ),
                                 nCollisionsElastic( 0 ), offlineDataGeneric() {};
    offlineDataCollisionNumbers( const double _time, const int _nCollisions22, const int _nCollisions23, 
                                 const int _nCollisions32, const int _nCollisionsElastic ) : time( _time ),
                                 nCollisions22( _nCollisions22 ), nCollisions23( _nCollisions23 ), nCollisions32( _nCollisions32 ),
                                 nCollisionsElastic( _nCollisionsElastic ), offlineDataGeneric() {};
    ~offlineDataCollisionNumbers() {};
    
    size_t getSize() const { return sizeof( offlineDataCollisionNumbers ); }
    std::string getFilenameIdentifier() const { return this->filenameIdentifier; }
    static std::string filenameIdentifier;
    
    /** @brief Current simulation time */
    double time;
    /** @brief Number of 2->2 interactions (divided by the test particle number) */
    int nCollisions22;
    /** @brief Number of 2->3 interactions (divided by the test particle number) */
    int nCollisions23;
    /** @brief Number of 3->2 interactions (divided by the test particle number) */
    int nCollisions32;
    /** @brief Number of geometric collisions (divided by the test particle number) */
    int nCollisionsElastic;
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & boost::serialization::base_object<offlineDataGeneric>(*this);
      ar & time;
      ar & nCollisions22;
      ar & nCollisions23;
      ar & nCollisions32;
      ar & nCollisionsElastic;
    }
};


class offlineDataInteractionRates : public offlineDataGeneric
{
  public:
    offlineDataInteractionRates( ) : offlineDataGeneric() {};
    offlineDataInteractionRates( std::vector< vector<double> >& _gluonRates,
                                 std::vector< vector<double> >& _quarkRates,
                                 std::vector< vector<double> >& _antiQuarkRates ) : gluonRates(_gluonRates),
                                 quarkRates(_quarkRates), antiQuarkRates(_antiQuarkRates), offlineDataGeneric() {};
    ~offlineDataInteractionRates() {};

    size_t getSize() const { return sizeof( offlineDataInteractionRates ); }
    std::string getFilenameIdentifier() const { return this->filenameIdentifier; }
    static std::string filenameIdentifier;
                                 
    /** @brief Stores the rates of gluons for output needed for offline reconstruction */
    std::vector< std::vector<double> > gluonRates;
    /** @brief Stores the rates of light quarks for output needed for offline reconstruction */
    std::vector< std::vector<double> > quarkRates;
    /** @brief Stores the rates of light anti-quarks for output needed for offline reconstruction */
    std::vector< std::vector<double> > antiQuarkRates;
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & boost::serialization::base_object<offlineDataGeneric>(*this);
      ar & gluonRates;
      ar & quarkRates;
      ar & antiQuarkRates;
    }
};


// "Export" the classes derived from offlineDataGeneric in case they need to be archived via boost:serialization
// using pointers to the base class
BOOST_CLASS_EXPORT_KEY( offlineDataInteraction22 )
BOOST_CLASS_EXPORT_KEY( offlineDataInteraction23 )
BOOST_CLASS_EXPORT_KEY( offlineDataInteraction32 )
BOOST_CLASS_EXPORT_KEY( offlineDataInteractionElastic )
BOOST_CLASS_EXPORT_KEY( offlineDataCellConfiguration )
BOOST_CLASS_EXPORT_KEY( offlineDataParticleIdSwap )
BOOST_CLASS_EXPORT_KEY( offlineDataParticleNumbers )
BOOST_CLASS_EXPORT_KEY( offlineDataCollisionNumbers )
BOOST_CLASS_EXPORT_KEY( offlineDataInteractionRates )
BOOST_CLASS_EXPORT_KEY( offlineDataEventType )

BOOST_CLASS_IMPLEMENTATION( offlineDataGeneric, boost::serialization::object_serializable)
BOOST_CLASS_IMPLEMENTATION( offlineDataInteraction22, boost::serialization::object_serializable)
BOOST_CLASS_IMPLEMENTATION( offlineDataInteraction23, boost::serialization::object_serializable)
BOOST_CLASS_IMPLEMENTATION( offlineDataInteraction32, boost::serialization::object_serializable)
BOOST_CLASS_IMPLEMENTATION( offlineDataInteractionElastic, boost::serialization::object_serializable)
BOOST_CLASS_IMPLEMENTATION( offlineDataCellConfiguration, boost::serialization::object_serializable)
BOOST_CLASS_IMPLEMENTATION( offlineDataParticleIdSwap, boost::serialization::object_serializable)
BOOST_CLASS_IMPLEMENTATION( offlineDataParticleNumbers, boost::serialization::object_serializable)
BOOST_CLASS_IMPLEMENTATION( offlineDataCollisionNumbers, boost::serialization::object_serializable)
BOOST_CLASS_IMPLEMENTATION( offlineDataInteractionRates, boost::serialization::object_serializable)
BOOST_CLASS_IMPLEMENTATION( offlineDataEventType, boost::serialization::object_serializable)

BOOST_CLASS_TRACKING(offlineDataGeneric, boost::serialization::track_never)
BOOST_CLASS_TRACKING(offlineDataInteraction22, boost::serialization::track_never)
BOOST_CLASS_TRACKING(offlineDataInteraction23, boost::serialization::track_never)
BOOST_CLASS_TRACKING(offlineDataInteraction32, boost::serialization::track_never)
BOOST_CLASS_TRACKING(offlineDataInteractionElastic, boost::serialization::track_never)
BOOST_CLASS_TRACKING(offlineDataCellConfiguration, boost::serialization::track_never)
BOOST_CLASS_TRACKING(offlineDataParticleIdSwap, boost::serialization::track_never)
BOOST_CLASS_TRACKING(offlineDataParticleNumbers, boost::serialization::track_never)
BOOST_CLASS_TRACKING(offlineDataCollisionNumbers, boost::serialization::track_never)
BOOST_CLASS_TRACKING(offlineDataInteractionRates, boost::serialization::track_never)
BOOST_CLASS_TRACKING(offlineDataEventType, boost::serialization::track_never)


/** @brief exception class for handling unexpected critical behaviour within simulations of heavy ion collisions  */
class eOfflineOutput_error : public std::runtime_error
{
  public:
    explicit eOfflineOutput_error(const std::string& what) : std::runtime_error(what) {};
    
    virtual ~eOfflineOutput_error() throw() {};
};


#endif