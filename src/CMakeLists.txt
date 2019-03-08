# name of main executable:
SET( MAIN offlineBAMPS )

# files containing stuff specific to the targets in this directory:
SET( SOURCES_lib
  additionalparticlesdistribution.cpp
  analysis.cpp
  configuration.cpp
  fixedshower_event.f
  hadronization_hq.cpp
  initialmodel_jpsi.cpp
  initialmodel_mcatnlo.cpp
  initialmodel_minijets.cpp
  initialmodel_particlesFromFile.cpp
  initialmodel_pythia.cpp
  initialmodel_pythiashower.cpp
  interpolation_iniJpsi.cpp
  interpolation_nJpsi.cpp
  offlineheavyioncollison.cpp
  offlineoutput.cpp
  particleOffline.cpp
  particleInFormTime.cpp
  pythia6425.f
  pythiashower_event.f90
)

IF( Pythia_FOUND )
  SET( SOURCES_lib
    ${SOURCES_lib}
    mesonDecayData.cpp
    mesonDecay.cpp
  )
ENDIF( )

# files which are equal in 'full' and 'full/offline':
SET( SOURCES_lib_common
  cellcontainer.cpp
  coordinateBins.cpp
  initialmodel.cpp
  initialmodel_cgc.cpp
  offlineoutput.cpp
  particle.cpp
  pdfinterface.cpp
  ringcontainer.cpp
  ringstructure.cpp
  woodsaxon.cpp
)

# files containing the main program:
SET( SOURCES_main
  offlineBAMPS.cpp
)


######### Include the header files of the BAMPS library #########
IF ( LINK_TO_BAMPS_LIB AND BAMPS_FOUND )
  INCLUDE_DIRECTORIES(${BAMPS_INCLUDE_DIRS})
ELSE ()
  INCLUDE_DIRECTORIES(${BAMPS_LIBRARY_PATH})
ENDIF ()

######### Create revision.h at build time #########

# include the output directory, where the revision.h file is generated
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

# a custom (empty) target that is always built
ADD_CUSTOM_TARGET(revisiontag ALL)

# a command that is run when target "revisiontag" is executed, creates revision.h using cmake script
IF ( NOT DISABLE_BUILD_EXTRAS )
  ADD_CUSTOM_COMMAND(TARGET revisiontag
    COMMAND ${CMAKE_COMMAND} -DMAIN_DIR=${CMAKE_SOURCE_DIR} -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR} -DLIB_SOURCE_DIR=${BAMPS_LIBRARY_PATH} -P ${CMAKE_MODULE_PATH}/getsvn.cmake)
ELSE ()
  ADD_CUSTOM_COMMAND(TARGET revisiontag
    COMMAND ${CMAKE_COMMAND} -DMAIN_DIR=${CMAKE_SOURCE_DIR} -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR} -DDISABLE_SVN_EXTRAS=TRUE -P ${CMAKE_MODULE_PATH}/getsvn.cmake)
ENDIF ()

# revision.h is a generated file
SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/revision.h PROPERTIES GENERATED TRUE HEADER_FILE_ONLY TRUE)

######### create a configBAMPS.h #########
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/configBAMPS.h.in ${CMAKE_CURRENT_BINARY_DIR}/configBAMPS.h)

######### Build the executable #########
ADD_LIBRARY(offline ${SOURCES_lib})
ADD_LIBRARY(fullcommon ${SOURCES_lib_common})
ADD_EXECUTABLE(${MAIN} ${SOURCES_main})
TARGET_LINK_LIBRARIES(${MAIN} offline fullcommon)

### activate interprocedural optimization if the compiler supports it
IF ( CMAKE_BUILD_TYPE STREQUAL "Release" )
  SET_PROPERTY(TARGET ${MAIN} PROPERTY INTERPROCEDURAL_OPTIMIZATION 1)  # -ipo option for intel compilers
ENDIF ()

IF ( LINK_TO_BAMPS_LIB AND BAMPS_FOUND )
  TARGET_LINK_LIBRARIES(${MAIN} ${BAMPS_LIBRARIES})    ## Link to static BAMPS library ##
ELSE ()
  TARGET_LINK_LIBRARIES(${MAIN} bamps)     ## Link to static BAMPS library ##
  ADD_DEPENDENCIES(${MAIN} bamps)          ## Add dependency such that libbamps.a is always build first ##
ENDIF ()

TARGET_LINK_LIBRARIES(${MAIN} ${Boost_LIBRARIES})
TARGET_LINK_LIBRARIES(${MAIN} ${GSL_LIBRARIES})
TARGET_LINK_LIBRARIES(${MAIN} ${ROOT_LIBRARIES})

IF ( LHAPDF_FOUND AND Pythia_FOUND )
  TARGET_LINK_LIBRARIES(${MAIN} ${LHAPDF_LIBRARIES} ${Pythia_LIBRARIES})
ELSEIF ( Pythia_FOUND )
  TARGET_LINK_LIBRARIES(${MAIN} ${Pythia_LIBRARIES} ${Pythia_LHAPDFDummy_LIBRARY})
ELSEIF ( LHAPDF_FOUND )
  TARGET_LINK_LIBRARIES(${MAIN} ${LHAPDF_LIBRARIES})
ENDIF ( LHAPDF_FOUND AND Pythia_FOUND )
  
#install(TARGETS ${MAIN} ${INSTALL_TARGETS_DEFAULT_ARGS})

#ADD_DEPENDENCIES(offline revisiontag)
#ADD_DEPENDENCIES(fullcommon revisiontag)
ADD_DEPENDENCIES(${MAIN} revisiontag)
