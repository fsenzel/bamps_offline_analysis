######### handle generation of doxygen files via CMake #########
## see: http://snikt.net/index.php/2010/04/01/howto-use-cmake-with-cc-projects


######### check whether doxygen is installed #########
FIND_PACKAGE(Doxygen)
IF (DOXYGEN_FOUND STREQUAL "NO")
  MESSAGE(FATAL_ERROR "!! Doxygen not found.")
ELSE()
  MESSAGE(STATUS "** Found doxygen")
ENDIF (DOXYGEN_FOUND STREQUAL "NO")

######### prepare doxygen configuration file #########
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

######### add doxygen as target #########
ADD_CUSTOM_TARGET(doxygen ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

######### cleanup $build/doc on "make clean" #########
SET_PROPERTY(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES doc)

#########  add doxygen as dependency to doc-target ######### 
GET_TARGET_PROPERTY(DOC_TARGET doc TYPE)
IF(NOT DOC_TARGET)
  ADD_CUSTOM_TARGET(doc)
ENDIF()
ADD_DEPENDENCIES(doc doxygen)