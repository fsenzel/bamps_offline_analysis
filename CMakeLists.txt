CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
get_filename_component(LIBSRC_DIR "${CMAKE_SOURCE_DIR}/lib" REALPATH)
SET(CMAKE_MODULE_PATH "${LIBSRC_DIR}/../cmake_modules/")

######### This sets custom build flags and includes options for various compilers #########
INCLUDE(SetCompilersAndBuildFlags)

######### Set project name  #########
PROJECT(offlineBAMPS)

######### Enabling Fortran for PYTHIA shower  #########
ENABLE_LANGUAGE(Fortran)

######### Activate packaging, use as "make package_source"  #########
INCLUDE(CMakeCustomPackageOptions.cmake)
INCLUDE(CPack)

######### Setup host properties (Vector, RDTSCP, ...) #########
INCLUDE(HostProperties)
SetVectorImplType()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_FLAGS_MARCH}")

######### check that C++11 code compiles #########
include(CheckCXXSourceCompiles)
CHECK_CXX_SOURCE_COMPILES("
#include <initializer_list>
struct A {};
struct B : public A { using A::A; };
template<typename T = void> int f() { return 1; }
int main() {
  auto l = []() { return f(); };
  return l() - 1;
}
" CXX11_COMPILES)

if(NOT CXX11_COMPILES)
#   message(FATAL_ERROR "Your compiler did not manage to compile a simple C++11 program. Please get a newer C++ compiler.")
   message(WARNING "Your compiler did not manage to compile a simple C++11 program. Please get a newer C++ compiler."
   "It may be okay for this code. Keep fingers crossed that it manages the C++11 features used.")
endif()

# ######### check for -fdiagnostics-color compiler option #########
# include(CheckCXXCompilerFlag)
# CHECK_CXX_COMPILER_FLAG("-fdiagnostics-color=always" HAS_DIAGNOSTICSCOLOR)
# if(HAS_DIAGNOSTICSCOLOR)
# #  message(STATUS "Current compiler supports -fdiagnostics-color=always option.")
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")
# elseif()
# #   message(STATUS "Current compiler does not support -fdiagnostics-color=always option.")
# endif()

######### Find Boost libraries #########
IF(NOT BOOST_ROOT AND NOT DEFINED ENV{BOOST_ROOT})
  SET(BOOST_ROOT "$ENV{HOME}/usr")
ENDIF()
MESSAGE(STATUS "** Boost Hint: ${BOOST_ROOT}")
FIND_PACKAGE(Boost 1.40 COMPONENTS serialization system filesystem program_options regex REQUIRED)
MESSAGE(STATUS "** Boost Include: ${Boost_INCLUDE_DIR}")
MESSAGE(STATUS "** Boost Libraries: ${Boost_LIBRARY_DIRS}")

######### Find Cuba library #########
FIND_PACKAGE(Cuba)

######### Find LHAPDF library #########
FIND_PACKAGE(LHAPDF)

######### Find Pythia8 library #########
FIND_PACKAGE(Pythia)

########## Find GSL library #########
FIND_PACKAGE(GSL REQUIRED)

######### This option controls whether optional goodies (subversion revison number and doxygen generation) should be disabled #########
## Default: FALSE
## Change on the command line via -DDISABLE_BUILD_EXTRAS=TRUE
SET ( DISABLE_BUILD_EXTRAS FALSE CACHE BOOL "Whether to disable some goodies [svn rev number and doxygen targets]" )

######### Find svn and retrieve revision number #########
IF( NOT DISABLE_BUILD_EXTRAS )
  FIND_PACKAGE(Subversion)
  IF(Subversion_FOUND)
    Subversion_WC_INFO(${CMAKE_CURRENT_SOURCE_DIR} ER)
    SET(SUBVERSION_REVISION ${ER_WC_REVISION})
    SET(CPACK_PACKAGE_VERSION_PATCH "${SUBVERSION_REVISION}")
  ENDIF(Subversion_FOUND)
ENDIF( NOT DISABLE_BUILD_EXTRAS )

######### This option controls whether linking to an installed version of the BAMPS lib should be tried  #########
## instead of building the BAMPS library directly. 
## Default: FALSE
## Change on the command line via -DLINK_TO_BAMPS_LIB=TRUE
SET ( LINK_TO_BAMPS_LIB FALSE CACHE BOOL "Whether to link to an installed version of the BAMPS lib instead of building it" )

######### Try to find the BAMPS library in case direct linking is requested  #########
IF ( LINK_TO_BAMPS_LIB )
  FIND_PACKAGE( BAMPS )
  IF ( NOT BAMPS_FOUND )
    MESSAGE( STATUS "!! Using direct build of BAMPS library as fallback" )
  ENDIF( NOT BAMPS_FOUND )
ENDIF( LINK_TO_BAMPS_LIB )

######### This option controls whether OpenMP should be used #########
## Default: FALSE
## Change on the command line via -DWITH_OPENMP=TRUE
SET ( WITH_OPENMP FALSE CACHE BOOL "Whether to use OpenMP or not" )

######### Find OpenMP and enable/disable it according flag #########

IF( WITH_OPENMP )
  FIND_PACKAGE( OpenMP )
  IF( OPENMP_FOUND )
#    MESSAGE(STATUS "** OpenMP found and will be used.")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -D_GLIBCXX_PARALLEL")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  ENDIF()
ELSE()
#    MESSAGE(STATUS "** OpenMP found, but will *not* be used.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
ENDIF()

######### Handle doxygen files via CMake, adds custom target, use as "make doxygen" #########
IF( NOT DISABLE_BUILD_EXTRAS )
  INCLUDE(ConfigureDoxygenTarget)
ENDIF( NOT DISABLE_BUILD_EXTRAS )
  

######### Set include and link directories #########
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${Cuba_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIRS})
IF ( LHAPDF_FOUND )
  INCLUDE_DIRECTORIES(${LHAPDF_INCLUDE_DIRS})
ENDIF ( LHAPDF_FOUND )
IF ( Pythia_FOUND )
  INCLUDE_DIRECTORIES(${Pythia_INCLUDE_DIRS})
ENDIF ( Pythia_FOUND )

######### Set the explicit path to the BAMPS library source directory #########
SET( BAMPS_LIBRARY_PATH ${CMAKE_SOURCE_DIR}/lib/ CACHE PATH "directory in which the BAMPS library source files reside" )

######### Set an option to disable the generation of install targets for the BAMPS library #########
OPTION( GENERATE_BAMPS_LIB_INSTALL_TARGETS "Whether to generate install targets for the BAMPS library or not" OFF )


######### Add source directories  #########
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(tests)

## Only process lib/src/ subirectory in case LINK_TO_BAMPS_LIB is false OR the bamps library has not been found ##
IF ( ( NOT LINK_TO_BAMPS_LIB ) OR ( NOT BAMPS_FOUND ) )
##
  ADD_SUBDIRECTORY(${BAMPS_LIBRARY_PATH})
##
ENDIF ( ( NOT LINK_TO_BAMPS_LIB ) OR ( NOT BAMPS_FOUND ) )
