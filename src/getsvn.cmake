# the FindSubversion.cmake module is part of the standard distribution
INCLUDE(FindSubversion)

# extract working copy information for SOURCE_DIR into MY_XXX variables
IF ( NOT DISABLE_SVN_EXTRAS )
  Subversion_WC_INFO(${SOURCE_DIR} MY)
ELSE ( NOT DISABLE_SVN_EXTRAS )
  SET( MY_WC_REVISION "\"no svn revision info available\"" )
ENDIF( NOT DISABLE_SVN_EXTRAS )

# write a file with the SVNVERSION define
FILE(WRITE revision.h.temp "#ifndef REVISION_H
#define REVISION_H

#define SVN_REVISION ${MY_WC_REVISION}

#endif")

# copy the file to the final header only if the version changes
# reduces needless rebuilds
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy_if_different revision.h.temp revision.h)