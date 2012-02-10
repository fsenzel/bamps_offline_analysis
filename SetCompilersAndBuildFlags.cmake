######### This sets custom build flags and includes options for various compilers #########

SET(BAMPS_BUILD_TYPE "" CACHE STRING "Possible build types are: GCC, Intel, IntelCSC, IntelLoeweCSC, AMD, AMDCSC, AMDLoeweCSC" )
MARK_AS_ADVANCED( BAMPS_BUILD_TYPE )

IF (BAMPS_BUILD_TYPE STREQUAL "Intel")
  MESSAGE(STATUS "** Compiling with Intel settings **")
  SET(ENV{CXX} "icpc")
  IF (CMAKE_BUILD_TYPE STREQUAL "Release")
    SET(ENV{CXXFLAGS} "-O3 -no-prec-div -w -fp-model fast=2 -xHost")
  ENDIF()
  ##
ELSEIF (BAMPS_BUILD_TYPE STREQUAL "IntelParallel")
  MESSAGE(STATUS "** Compiling with Intel parallel settings **")
  SET(ENV{CXX} "icpc")
  IF (CMAKE_BUILD_TYPE STREQUAL "Release")
    SET(ENV{CXXFLAGS} "-parallel -O3 -no-prec-div -w -fp-model fast=2 -xHost")
  ENDIF()  
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "IntelLoeweCSC")
  MESSAGE(STATUS "** Compiling with Intel settings for LOEWE CSC**")
  SET(ENV{CXX} "icpc")
  IF (CMAKE_BUILD_TYPE STREQUAL "Release")
    SET(ENV{CXXFLAGS} "-O3 -simd -msse3 -no-prec-div -w -fp-relaxed -fp-model fast=2 -xHost")
  ENDIF()
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "IntelCSC")
  MESSAGE(STATUS "** Compiling with IntelFast settings for CSC**")
  SET(ENV{CXX} "icpc")
  IF (CMAKE_BUILD_TYPE STREQUAL "Release")
    SET(ENV{CXXFLAGS} "-O3 -no-prec-div -fp-relaxed -w -fp-model fast=2 -xHost")
  ENDIF()
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "IntelProfGen")
  MESSAGE(STATUS "** Compiling with Intel generating profile information**")
  SET(ENV{CXX} "icpc")
  IF (CMAKE_BUILD_TYPE STREQUAL "Release")
    SET(ENV{CXXFLAGS} "-O2 -prof-gen -w")
  ELSEIF(CMAKE_BUILD_TYPE STREQUAL "Debug")
    SET(CMAKE_CXX_FLAGS_DEBUG "-g -prof-gen")
  ENDIF()
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "IntelProfUse")
  MESSAGE(STATUS "** Compiling with Intel using profile information**")
  SET(ENV{CXX} "icpc")
  IF (CMAKE_BUILD_TYPE STREQUAL "Release")
    SET(ENV{CXXFLAGS} "-prof-use -prof-dir${CMAKE_CURRENT_BINARY_DIR} -O3 -no-prec-div -w -fp-model fast=2 -xHost")
  ELSEIF(CMAKE_BUILD_TYPE STREQUAL "Debug")
    SET(CMAKE_CXX_FLAGS_DEBUG "-prof-use -prof-dir${CMAKE_CURRENT_BINARY_DIR} -g")
  ENDIF()
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "IntelProfUseParallel")
  MESSAGE(STATUS "** Compiling with Intel using profile information**")
  SET(ENV{CXX} "icpc")
  IF (CMAKE_BUILD_TYPE STREQUAL "Release")
    SET(ENV{CXXFLAGS} "-prof-use -prof-dir${CMAKE_CURRENT_BINARY_DIR} -parallel -O3 -no-prec-div -w -fp-model fast=2 -xHost")
  ENDIF()
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "AMDLoeweCSC")
  MESSAGE(STATUS "** Compiling with AMD settings for LOEWE CSC**")
  SET(ENV{CXX} "openCC")
  IF (CMAKE_BUILD_TYPE STREQUAL "Release")
    SET(ENV{CXXFLAGS} "-Ofast -DNDEBUG")
  ENDIF()
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "AMDCSC")
  MESSAGE(STATUS "** Compiling with AMD settings for LOEWE CSC**")
  SET(ENV{CXX} "openCC")
  IF (CMAKE_BUILD_TYPE STREQUAL "Release")
    SET(ENV{CXXFLAGS} "-Ofast -DNDEBUG")
  ENDIF()
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "AMD")
  MESSAGE(STATUS "** Compiling with AMD settings**")
  SET(ENV{CXX} "openCC")
  IF (CMAKE_BUILD_TYPE STREQUAL "Release")
    SET(ENV{CXXFLAGS} "-Ofast -DNDEBUG")
  ENDIF()
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "GCC")
  MESSAGE(STATUS "** Compiling with GCC (g++) settings**")
  SET(ENV{CXX} "g++")
  ##
ELSE()
  MESSAGE(STATUS "** Compiling with default settings**")
  MESSAGE(STATUS "** Use -DBAMPS_BUILD_TYPE to specify predefined sets of build options. Possible Values are: GCC, Intel, IntelCSC,
          IntelLoeweCSC, AMD, AMDCSC, AMDLoeweCSC **")
ENDIF()

IF( NOT CMAKE_BUILD_TYPE )
  SET(CMAKE_BUILD_TYPE Release)
ENDIF( NOT CMAKE_BUILD_TYPE )