######### This sets custom build flags and includes options for various compilers #########

### This file has to be included and called before any calls to 'PROJECT()' or
### 'ENABLE_LANGUAGE()' in CMakeLists.txt, since we set the compiler via the 
### environment variable. This can be done only at the *initial* configuration.
###
### We also set some CMAKE_CXX_<TYPE> variables. They have to be initialized
### with the 'CACHE' keyword, since they would be overwritten by
### the first call of Modules/CMakeCXXInformation.cmake. 
### (The given descriptions are copied from this file.)
### Values given at the prompt are respected.
###
### Changing the value of CMAKE_BUILD_TYPE at the prompt is allowed.

SET(BAMPS_BUILD_TYPE "GCC" CACHE STRING "Possible build types are: GCC, Intel, IntelCSC, IntelLoeweCSC, AMD, AMDCSC, AMDLoeweCSC" )
MARK_AS_ADVANCED( BAMPS_BUILD_TYPE )

IF( NOT CMAKE_BUILD_TYPE )
  SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel" FORCE)
ENDIF( NOT CMAKE_BUILD_TYPE )


#########
## We may define own CMAKE_BUILD_TYPE:

## CMAKE_BUILD_TYPE == DebugProf :

SET(CMAKE_CXX_FLAGS_DEBUGPROF "-O3 -g -pg -Wall -Wunused-variable" CACHE STRING 
  "Flags used by the compiler during profiling builds.")
SET(CMAKE_C_FLAGS_DEBUGPROF "-O3 -g -pg -Wall -Wunused-variable" CACHE STRING 
  "Flags used by the compiler during profiling builds.")

#########

IF (BAMPS_BUILD_TYPE STREQUAL "Intel")
  MESSAGE(STATUS "** Compiling with Intel settings **")
  SET(ENV{CXX} "icpc")
  SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -no-prec-div -w -fp-model fast=2 -xHost" CACHE STRING 
    "Flags used by the compiler during release builds (/MD /Ob1 /Oi /Ot /Oy /Gs will produce slightly less optimized but smaller files).")
  SET(CMAKE_CXX_FLAGS_DEBUG "-g" CACHE STRING 
    "Flags used by the compiler during debug builds.") 
  ##
ELSEIF (BAMPS_BUILD_TYPE STREQUAL "IntelParallel")
  MESSAGE(STATUS "** Compiling with Intel parallel settings **")
  SET(ENV{CXX} "icpc")
  SET(CMAKE_CXX_FLAGS_RELEASE "-parallel -O3 -no-prec-div -w -fp-model fast=2 -xHost" CACHE STRING 
    "Flags used by the compiler during release builds (/MD /Ob1 /Oi /Ot /Oy /Gs will produce slightly less optimized but smaller files).")
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "IntelLoeweCSC")
  MESSAGE(STATUS "** Compiling with Intel settings for LOEWE CSC**")
  SET(ENV{CXX} "icpc")
  SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -simd -msse3 -no-prec-div -w -fp-relaxed -fp-model fast=2 -xHost" CACHE STRING 
    "Flags used by the compiler during release builds (/MD /Ob1 /Oi /Ot /Oy /Gs will produce slightly less optimized but smaller files).")
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "IntelCSC")
  MESSAGE(STATUS "** Compiling with IntelFast settings for CSC**")
  SET(ENV{CXX} "icpc")
  SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -no-prec-div -fp-relaxed -w -fp-model fast=2 -xHost" CACHE STRING 
    "Flags used by the compiler during release builds (/MD /Ob1 /Oi /Ot /Oy /Gs will produce slightly less optimized but smaller files).")
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "IntelProfGen")
  MESSAGE(STATUS "** Compiling with Intel generating profile information**")
  SET(ENV{CXX} "icpc")
  SET(CMAKE_CXX_FLAGS_RELEASE "-O2 -prof-gen -w" CACHE STRING 
    "Flags used by the compiler during release builds (/MD /Ob1 /Oi /Ot /Oy /Gs will produce slightly less optimized but smaller files).")
  SET(CMAKE_CXX_FLAGS_DEBUG "-g -prof-gen" CACHE STRING 
    "Flags used by the compiler during debug builds.")
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "IntelProfUse")
  MESSAGE(STATUS "** Compiling with Intel using profile information**")
  SET(ENV{CXX} "icpc")
  SET(CMAKE_CXX_FLAGS_RELEASE "-prof-use -prof-dir${CMAKE_CURRENT_BINARY_DIR} -O3 -no-prec-div -w -fp-model fast=2 -xHost" CACHE STRING 
    ".Flags used by the compiler during release builds (/MD /Ob1 /Oi /Ot /Oy /Gs will produce slightly less optimized but smaller files)...")
  SET(CMAKE_CXX_FLAGS_DEBUG "-prof-use -prof-dir${CMAKE_CURRENT_BINARY_DIR} -g" CACHE STRING 
    "Flags used by the compiler during debug builds.")
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "IntelProfUseParallel")
  MESSAGE(STATUS "** Compiling with Intel using profile information**")
  SET(ENV{CXX} "icpc")
  SET(CMAKE_CXX_FLAGS_RELEASE "-prof-use -prof-dir${CMAKE_CURRENT_BINARY_DIR} -parallel -O3 -no-prec-div -w -fp-model fast=2 -xHost" CACHE STRING 
    "Flags used by the compiler during release builds (/MD /Ob1 /Oi /Ot /Oy /Gs will produce slightly less optimized but smaller files).")
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "AMDLoeweCSC")
  MESSAGE(STATUS "** Compiling with AMD settings for LOEWE CSC**")
  SET(ENV{CXX} "openCC")
  SET(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG" CACHE STRING 
    "Flags used by the compiler during release builds (/MD /Ob1 /Oi /Ot /Oy /Gs will produce slightly less optimized but smaller files).")
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "AMDCSC")
  MESSAGE(STATUS "** Compiling with AMD settings for LOEWE CSC**")
  SET(ENV{CXX} "openCC")
  SET(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG" CACHE STRING 
    "Flags used by the compiler during release builds (/MD /Ob1 /Oi /Ot /Oy /Gs will produce slightly less optimized but smaller files).")
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "AMD")
  MESSAGE(STATUS "** Compiling with AMD settings**")
  SET(ENV{CXX} "openCC")
  SET(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG" CACHE STRING 
    "Flags used by the compiler during release builds (/MD /Ob1 /Oi /Ot /Oy /Gs will produce slightly less optimized but smaller files).")
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "GCC")
  MESSAGE(STATUS "** Compiling with GCC (g++) settings**")
  SET(ENV{CXX} "g++")
  ##
ELSE()
  MESSAGE(STATUS "** Compiling with default settings**")
  MESSAGE(STATUS "** Use -DBAMPS_BUILD_TYPE to specify predefined sets of build options. Possible Values are: GCC, Intel, IntelCSC, IntelLoeweCSC, AMD, AMDCSC, AMDLoeweCSC **")
ENDIF()

#########

SET(VECTOR_IMPL_TYPE "SSE" CACHE STRING "Possible vector implementation types are: Scalar, SSE" )
MARK_AS_ADVANCED( VECTOR_IMPL_TYPE )

IF (VECTOR_IMPL_TYPE STREQUAL "SSE")
  MESSAGE(STATUS "** Vector implementation: SSE")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2 -msse3" CACHE STRING 
    "Flags used by the compiler during all build types.") 
ELSEIF (VECTOR_IMPL_TYPE STREQUAL "Scalar")
  MESSAGE(STATUS "** Vector implementation: Scalar")	
ELSE()
  MESSAGE(STATUS "** Vector implementation: undefined ==> Scalar !!!")
  SET(VECTOR_IMPL_TYPE "Scalar")
ENDIF()

#########

# enable C++11 features:
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x" CACHE STRING "..." FORCE)

#########