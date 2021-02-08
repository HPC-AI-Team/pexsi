#   FindParMETIS.cmake
#
#   Finds the ParMETIS library.
#
#   This module will define the following variables:
#   
#     ParMETIS_FOUND        - System has found ParMETIS installation
#     ParMETIS_INCLUDE_DIR  - Location of ParMETIS headers
#     ParMETIS_LIBRARIES    - ParMETIS libraries
#     ParMETIS_USES_ILP64   - Whether ParMETIS was compiled with ILP64

#   This module can handle the following COMPONENTS
#
#     ilp64 - 64-bit index integers
#
#   This module will export the following targets if ParMETIS_FOUND
#
#     ParMETIS::ParMETIS
#
#
#
#
#   Proper usage:
#
#     project( TEST_FIND_ParMETIS C )
#     find_package( ParMETIS )
#
#     if( ParMETIS_FOUND )
#       add_executable( test test.cxx )
#       target_link_libraries( test ParMETIS::ParMETIS )
#     endif()
#
#
#
#
#   This module will use the following variables to change
#   default behaviour if set
#
#     ParMETIS_PREFIX
#     ParMETIS_INCLUDE_DIR
#     ParMETIS_LIBRARY_DIR
#     ParMETIS_LIBRARIES

#==================================================================
#   Copyright (c) 2018 The Regents of the University of California,
#   through Lawrence Berkeley National Laboratory.  
#
#   Author: David Williams-Young
#   
#   This file is part of cmake-modules. All rights reserved.
#   
#   Redistribution and use in source and binary forms, with or without
#   modification, are permitted provided that the following conditions are met:
#   
#   (1) Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#   (2) Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#   (3) Neither the name of the University of California, Lawrence Berkeley
#   National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
#   be used to endorse or promote products derived from this software without
#   specific prior written permission.
#   
#   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#   ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
#   ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
#   ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#   
#   You are under no obligation whatsoever to provide any bug fixes, patches, or
#   upgrades to the features, functionality or performance of the source code
#   ("Enhancements") to anyone; however, if you choose to make your Enhancements
#   available either publicly, or directly to Lawrence Berkeley National
#   Laboratory, without imposing a separate written license agreement for such
#   Enhancements, then you hereby grant the following license: a non-exclusive,
#   royalty-free perpetual license to install, use, modify, prepare derivative
#   works, incorporate into other computer software, distribute, and sublicense
#   such enhancements or derivative works thereof, in binary and source code form.
#
#==================================================================

cmake_minimum_required( VERSION 3.11 ) # Require CMake 3.11+

include(CMakeFindDependencyMacro)

include( ${CMAKE_CURRENT_LIST_DIR}/util/CommonFunctions.cmake )
fill_out_prefix( ParMETIS )
# DEPENDENCIES

# Make sure C is enabled
get_property( ParMETIS_languages GLOBAL PROPERTY ENABLED_LANGUAGES )
if( NOT "C" IN_LIST ParMETIS_languages )
  message( FATAL_ERROR "C Language Must Be Enabled for ParMETIS Linkage" )
endif()
unset( ParMETIS_languages )


# MPI
if( NOT TARGET MPI::MPI_C )
  find_dependency( MPI )
endif()








# Try to find the header
find_path( ParMETIS_INCLUDE_DIR 
  NAMES parmetis.h
  HINTS ${ParMETIS_PREFIX}
  PATHS ${ParMETIS_INCLUDE_DIR}
  PATH_SUFFIXES include
  DOC "Location of ParMETIS header"
)

# Try to find libraries if not already set
if( NOT ParMETIS_LIBRARIES )

  # METIS
  copy_meta_data( ParMETIS METIS )
  find_dependency( METIS )

  find_library( ParMETIS_LIBRARIES
    NAMES parmetis
    HINTS ${ParMETIS_PREFIX}
    PATHS ${ParMETIS_LIBRARY_DIR}
    PATH_SUFFIXES lib lib64 lib32
    DOC "ParMETIS Libraries"
  )

  if( ParMETIS_LIBRARIES )
    list( APPEND ParMETIS_LIBRARIES METIS::METIS )
  endif()

else()

  foreach( _lib ${ParMETIS_LIBRARIES} )
    if( _lib MATCHES "METIS" )
      find_dependency(METIS)
    endif()
  endforeach()

endif()

# MPI Always a dependency
if( ParMETIS_LIBRARIES AND NOT MPI::MPI_C IN_LIST ParMETIS_LIBRARIES )
  list( APPEND ParMETIS_LIBRARIES MPI::MPI_C )
  set( ParMETIS_LIBRARIES "${ParMETIS_LIBRARIES}" CACHE INTERNAL "ParMETIS Libraries" FORCE )
endif()

# TODO Check if Link OK



# Check version
if( EXISTS ${ParMETIS_INCLUDE_DIR}/parmetis.h )

  set( version_pattern 
  "^#define[\t ]+ParMETIS_(MAJOR|MINOR|SUBMINOR)_VERSION[\t ]+([0-9\\.]+)$"
  )
  file( STRINGS ${ParMETIS_INCLUDE_DIR}/parmetis.h parmetis_version
        REGEX ${version_pattern} )
  
  foreach( match ${parmetis_version} )
  
    if(ParMETIS_VERSION_STRING)
      set(ParMETIS_VERSION_STRING "${ParMETIS_VERSION_STRING}.")
    endif()
  
    string(REGEX REPLACE ${version_pattern} 
      "${ParMETIS_VERSION_STRING}\\2" 
      ParMETIS_VERSION_STRING ${match}
    )
  
    set(ParMETIS_VERSION_${CMAKE_MATCH_1} ${CMAKE_MATCH_2})
  
  endforeach()
  
  unset( parmetis_version )
  unset( version_pattern )

endif()

# Check ILP64 (Inherits from METIS)
if( METIS_FOUND )

  set( ParMETIS_USES_ILP64 ${METIS_USES_ILP64} )

endif()

# Handle components
if( ParMETIS_USES_ILP64 )
  set( ParMETIS_ilp64_FOUND TRUE )
endif()




# Determine if we've found ParMETIS
mark_as_advanced( ParMETIS_FOUND ParMETIS_INCLUDE_DIR ParMETIS_LIBRARIES )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args( ParMETIS
  REQUIRED_VARS ParMETIS_LIBRARIES ParMETIS_INCLUDE_DIR
  VERSION_VAR ParMETIS_VERSION_STRING
  HANDLE_COMPONENTS
)

set( ParMETIS_LIBRARIES "${ParMETIS_LIBRARIES}" CACHE STRING "ParMETIS LIBRARIES" FORCE )
set( ParMETIS_INCLUDE_DIR "${ParMETIS_INCLUDE_DIR}" CACHE STRING "ParMETIS INCLUDE_DIR" FORCE )

# Export target
if( ParMETIS_FOUND AND NOT TARGET ParMETIS::ParMETIS )

  add_library( ParMETIS::ParMETIS INTERFACE IMPORTED )
  set_target_properties( ParMETIS::ParMETIS PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${ParMETIS_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES      "${ParMETIS_LIBRARIES}" 
  )

endif()
