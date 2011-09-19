######### This sets custom build flags and includes options for various compilers #########

SET(BAMPS_BUILD_TYPE "" CACHE STRING "Possible build types are: GCC, Intel, IntelCSC, IntelLoeweCSC, AMD, AMDCSC, AMDLoeweCSC" )
MARK_AS_ADVANCED( BAMPS_BUILD_TYPE )

IF (BAMPS_BUILD_TYPE STREQUAL "Intel")
  MESSAGE(STATUS "** Compiling with Intel settings **")
  SET(CMAKE_CXX_COMPILER "icpc")
  SET(CMAKE_CXX_FLAGS_RELEASE "-ipo -O3 -no-prec-div -fp-relaxed -fp-model fast=2 -xHost")
  SET(CMAKE_CXX_FLAGS_DEBUG "-g")
  ##
ELSEIF (BAMPS_BUILD_TYPE STREQUAL "IntelParallel")
  MESSAGE(STATUS "** Compiling with Intel settings **")
  SET(CMAKE_CXX_COMPILER "icpc")
  SET(CMAKE_CXX_FLAGS_RELEASE "-ipo -parallel -O3 -no-prec-div -fp-relaxed -fp-model fast=2 -xHost")
  SET(CMAKE_CXX_FLAGS_DEBUG "-g")
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "IntelLoeweCSC")
  MESSAGE(STATUS "** Compiling with Intel settings for LOEWE CSC**")
  SET(CMAKE_CXX_COMPILER "icpc")
  SET(CMAKE_CXX_FLAGS_RELEASE "-ipo -O3 -simd -msse3 -no-prec-div -fp-relaxed -fp-model fast=2 -xHost")
  SET(CMAKE_CXX_FLAGS_DEBUG "-g")
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "IntelCSC")
  MESSAGE(STATUS "** Compiling with IntelFast settings for CSC**")
  SET(CMAKE_CXX_COMPILER "icpc")
  SET(CMAKE_CXX_FLAGS_RELEASE "-ipo -O3 -no-prec-div -fp-relaxed -fp-model fast=2 -xHost")
  SET(CMAKE_CXX_FLAGS_DEBUG "-g")
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "IntelProfGen")
  MESSAGE(STATUS "** Compiling with Intel generating profile information**")
  SET(CMAKE_CXX_COMPILER "icpc")
  SET(CMAKE_CXX_FLAGS_RELEASE "-O2 -prof-gen")
  SET(CMAKE_CXX_FLAGS_DEBUG "-g -prof-gen")
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "IntelProfUse")
  MESSAGE(STATUS "** Compiling with Intel using profile information**")
  SET(CMAKE_CXX_COMPILER "icpc")
  SET(CMAKE_CXX_FLAGS_RELEASE "-prof-use -prof-dir./ -ipo -O3 -no-prec-div -fp-relaxed -fp-model fast=2 -xHost")
  SET(CMAKE_CXX_FLAGS_DEBUG "-prof-use -prof-dir./ -g")
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "IntelProfUseParallel")
  MESSAGE(STATUS "** Compiling with Intel using profile information**")
  SET(CMAKE_CXX_COMPILER "icpc")
  SET(CMAKE_CXX_FLAGS_RELEASE "-prof-use -parallel -ipo -O3 -no-prec-div -fp-relaxed -fp-model fast=2 -xHost")
  SET(CMAKE_CXX_FLAGS_DEBUG "-g")
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "AMDLoeweCSC")
  MESSAGE(STATUS "** Compiling with AMD settings for LOEWE CSC**")
  SET(CMAKE_CXX_COMPILER "openCC")
  SET(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG")
  SET(CMAKE_CXX_FLAGS_DEBUG "-g")
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "AMDCSC")
  MESSAGE(STATUS "** Compiling with AMD settings for LOEWE CSC**")
  SET(CMAKE_CXX_COMPILER "openCC")
  SET(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG")
  SET(CMAKE_CXX_FLAGS_DEBUG "-g")
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "AMD")
  MESSAGE(STATUS "** Compiling with AMD settings**")
  SET(CMAKE_CXX_COMPILER "openCC")
  SET(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG")
  SET(CMAKE_CXX_FLAGS_DEBUG "-g")
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "GCC")
  MESSAGE(STATUS "** Compiling with GCC (g++) settings**")
  SET(CMAKE_CXX_COMPILER "g++")
  SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
  SET(CMAKE_CXX_FLAGS_DEBUG "-g")
  ##
ELSE()
  MESSAGE(STATUS "** Compiling with default settings**")
  MESSAGE(STATUS "** Use -DBAMPS_BUILD_TYPE to specify predefined sets of build options. Possible Values are: GCC, Intel, IntelCSC,
          IntelLoeweCSC, AMD, AMDCSC, AMDLoeweCSC **")
ENDIF()

IF( NOT CMAKE_BUILD_TYPE )
  SET(CMAKE_BUILD_TYPE Release)
ENDIF( NOT CMAKE_BUILD_TYPE )