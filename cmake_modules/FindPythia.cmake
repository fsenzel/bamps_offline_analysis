# Try to find a Pythia 8.x installation
#
# Oliver Fochler
#
# This package defines
#  Pythia_FOUND - Pythia has been found
#  Pythia_INCLUDE_DIRS - The directory in which the Pythia headers reside
#  Pythia_LIBRARIES - The Pythia library / libraries
#  Pythia_LHAPDFDummy_LIBRARY - The LHAPDF dummy library, use if the real LHAPDF is not installed / found

# find the path of the include files, based on two headers
FIND_PATH( Pythia_INCLUDE_DIR 
    NAMES Pythia.h
    PATHS $ENV{HOME}/usr/include      # suggest a user based include tree
    PATH_SUFFIXES PYTHIA Pythia pythia PYTHIA8 Pythia8 pythia8  # suggest some path suffixes in which the headers could be located
)

# find the library
# explicitlly try to find the static library first, in case this fails, try the more generic version (i.e. search for libpythia8.so)
FIND_LIBRARY( Pythia_LIBRARY
    NAMES libpythia8.a
    PATHS $ENV{HOME}/usr/lib          # suggest a user based include tree
    PATH_SUFFIXES PYTHIA pythia Pythia PYTHIA8 Pythia8 pythia8  # suggest some path suffixes in which the headers could be located
)
IF( NOT Pythia_LIBRARY )
  FIND_LIBRARY( Pythia_LIBRARY
    NAMES pythia8
    PATHS $ENV{HOME}/usr/lib          # suggest a user based include tree
    PATH_SUFFIXES PYTHIA pythia Pythia PYTHIA8 Pythia8 pythia8  # suggest some path suffixes in which the headers could be located
  )
ENDIF()


# find the LHAPDF dummy library
# in case LHAPDF is NOT installed on the system, Pyhtia provides a dummy library it needs to be linked against
FIND_LIBRARY( Pythia_LHAPDFDummy_LIBRARY
  NAMES liblhapdfdummy.a
  PATHS $ENV{HOME}/usr/lib          # suggest a user based include tree
  PATH_SUFFIXES PYTHIA pythia Pythia PYTHIA8 Pythia8 pythia8  # suggest some path suffixes in which the headers could be located
)
IF( NOT Pythia_LHAPDFDummy_LIBRARY )
  FIND_LIBRARY( Pythia_LIBRARY
    NAMES lhapdfdummy
    PATHS $ENV{HOME}/usr/lib          # suggest a user based include tree
    PATH_SUFFIXES PYTHIA pythia Pythia PYTHIA8 Pythia8 pythia8  # suggest some path suffixes in which the headers could be located
  )
ENDIF()


# find the Pythia xmldoc directory that is needed at runtime
FIND_PATH( Pythia_xmldoc_PATH
  NAMES PartonDistributions.xml ParticleProperties.xmlParticleDecays.xml 
  HINTS $ENV{HOME}/usr/share ENV PYTHIA8DATA
  PATH_SUFFIXES PYTHIA pythia Pythia PYTHIA/xmldoc Pythia/xmldoc pythia/xmldoc
)


# adhere to the standard nomenclature for find_package routines and set some variables
SET( Pythia_INCLUDE_DIRS ${Pythia_INCLUDE_DIR} )
SET( Pythia_LIBRARIES ${Pythia_LIBRARY} )

# handle the QUIETLY and REQUIRED arguments and set Pythia_FOUND to TRUE
# if all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS( Pythia DEFAULT_MSG Pythia_LIBRARY Pythia_INCLUDE_DIR Pythia_xmldoc_PATH)
SET ( Pythia_FOUND ${PYTHIA_FOUND} CACHE INTERNAL "Provide Pythia_FOUND in addition to PYTHIA_FOUND" FORCE )


# display some status information
IF( Pythia_FOUND )
 MESSAGE( STATUS "** Pythia 8 library: ${Pythia_LIBRARIES}" )
 MESSAGE( STATUS "** Pythia 8 include: ${Pythia_INCLUDE_DIRS}" )
 MESSAGE( STATUS "** Pythia 8 xmldoc: ${Pythia_xmldoc_PATH}" )
ENDIF( Pythia_FOUND )

# the variables will only show up in the GUI in the "advanced" view
MARK_AS_ADVANCED(Pythia_INCLUDE_DIR Pythia_LIBRARY Pythia_LHAPDFDummy_LIBRARY )