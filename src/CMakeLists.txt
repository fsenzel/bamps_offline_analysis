SET(project_SRCS
additionalparticlesdistribution.cpp
analysis.cpp
binning2.cpp
cellcontainer.cpp
cgcInitialDistribution.cpp
configuration.cpp
coordinateBins.cpp
initialdistribution.cpp
initialstatemodel.cpp
integrand_time.cpp
minijets.cpp
offlineBAMPS.cpp
offlineheavyioncollison.cpp 
offlineoutput.cpp
particle.cpp
psf.cpp
pythiaInitialDistribution.cpp
ringcontainer.cpp
ringstructure.cpp
woodsaxon.cpp
)


######### Include the header files of the BAMPS library #########
IF( LINK_TO_BAMPS_LIB AND BAMPS_FOUND )
  INCLUDE_DIRECTORIES(${BAMPS_INCLUDE_DIRS})
ELSE( LINK_TO_BAMPS_LIB AND BAMPS_FOUND )
  INCLUDE_DIRECTORIES(${BAMPS_LIBRARY_PATH})
ENDIF( LINK_TO_BAMPS_LIB AND BAMPS_FOUND ) 

######### Create revision.h at build time #########

# include the output directory, where the svnversion.h file is generated
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

# a custom (empty) target that is always built
ADD_CUSTOM_TARGET(revisiontag ALL)

# a command that is run when target "revisiontag" is executed, creates revision.h using cmake script
IF ( NOT DISABLE_BUILD_EXTRAS )
  ADD_CUSTOM_COMMAND(TARGET revisiontag
    COMMAND ${CMAKE_COMMAND} -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR} -P ${CMAKE_CURRENT_SOURCE_DIR}/getsvn.cmake)
ELSE ( NOT DISABLE_BUILD_EXTRAS )
  ADD_CUSTOM_COMMAND(TARGET revisiontag
    COMMAND ${CMAKE_COMMAND} -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR} -DDISABLE_SVN_EXTRAS=TRUE -P ${CMAKE_CURRENT_SOURCE_DIR}/getsvn.cmake)
ENDIF ( NOT DISABLE_BUILD_EXTRAS )

# revision.h is a generated file
SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/revision.h PROPERTIES GENERATED TRUE HEADER_FILE_ONLY TRUE)



######### Build the executable #########
ADD_EXECUTABLE(offlineBAMPS ${project_SRCS} )

### activate interprocedural optimization if the compiler supports it
IF (CMAKE_BUILD_TYPE STREQUAL "Release")
  SET_PROPERTY(TARGET offlineBAMPS PROPERTY INTERPROCEDURAL_OPTIMIZATION 1)  # -ipo option for intel compilers
ENDIF()

IF( LINK_TO_BAMPS_LIB AND BAMPS_FOUND )
  TARGET_LINK_LIBRARIES(offlineBAMPS ${BAMPS_LIBRARIES})    ## Link to static BAMPS library ##
ELSE( LINK_TO_BAMPS_LIB AND BAMPS_FOUND )
  TARGET_LINK_LIBRARIES(offlineBAMPS bamps)     ## Link to static BAMPS library ##
  ADD_DEPENDENCIES(offlineBAMPS bamps)          ## Add dependency such that libbamps.a is always build first ##
ENDIF( LINK_TO_BAMPS_LIB AND BAMPS_FOUND ) 

TARGET_LINK_LIBRARIES(offlineBAMPS ${Boost_LIBRARIES} ${GSL_LIBRARIES})
#TARGET_LINK_LIBRARIES(offlineBAMPS ${Boost_LIBRARIES} ${GSL_LIBRARIES} ${Cuba_LIBRARIES})

#install(TARGETS offlineBAMPS ${INSTALL_TARGETS_DEFAULT_ARGS})

ADD_DEPENDENCIES(offlineBAMPS revisiontag)
