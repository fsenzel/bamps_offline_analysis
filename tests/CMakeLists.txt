
SET( BAMPS_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/src/" "${CMAKE_CURRENT_BINARY_DIR}/../src/" )
SET( BAMPS_LIBRARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/../src/" )
SET( BAMPS_LIBRARIES -lbamps )

INCLUDE_DIRECTORIES( ${BAMPS_INCLUDE_DIRS} )
IF ( ( NOT LINK_TO_BAMPS_LIB ) OR ( NOT BAMPS_FOUND ) )
  INCLUDE_DIRECTORIES( ${BAMPS_LIBRARY_PATH} )
ENDIF ( ( NOT LINK_TO_BAMPS_LIB ) OR ( NOT BAMPS_FOUND ) )

macro(bamps_add_test _name)
   set(name ${_name})
   add_executable(${_name} ${_name}.cpp)
   target_link_libraries(${_name} offline fullcommon bamps
      ${ROOT_LIBRARIES} 
      ${Boost_LIBRARIES} )
   IF ( Cuba_FOUND )
     TARGET_LINK_LIBRARIES( ${_name} ${Cuba_LIBRARIES} )
   ENDIF ( Cuba_FOUND )

   TARGET_LINK_LIBRARIES( ${_name} ${Boost_LIBRARIES} ${GSL_LIBRARIES} ${ROOT_LIBRARIES})

   IF( LHAPDF_FOUND AND Pythia_FOUND )
     TARGET_LINK_LIBRARIES(${_name} ${LHAPDF_LIBRARIES} ${Pythia_LIBRARIES})
   ELSEIF( Pythia_FOUND )
     TARGET_LINK_LIBRARIES(${_name} ${Pythia_LIBRARIES} ${Pythia_LHAPDFDummy_LIBRARY})
   ELSEIF( LHAPDF_FOUND )
     TARGET_LINK_LIBRARIES(${_name} ${LHAPDF_LIBRARIES})
   ENDIF( LHAPDF_FOUND AND Pythia_FOUND )

   add_test(${_name} "${CMAKE_CURRENT_BINARY_DIR}/${_name}")
endmacro(bamps_add_test)

ENABLE_TESTING()

bamps_add_test(21_pythiashower)

######### Copy additional files #########

function(copy_files)
  foreach (filename ${ARGV})
#    MESSAGE(STATUS "copying: ${filename}")
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/${filename} ${CMAKE_CURRENT_BINARY_DIR}/${filename})
  endforeach (filename)
endfunction(copy_files)

# copy_files( testNf3_32configurations.dat )
