#   FindSuperLU_DIST.cmake
#
#   Finds the SuperLU_DIST library.
#
#   This module will define the following variables:
#   
#     SuperLU_DIST_FOUND        - System has found SuperLU_DIST installation
#     SuperLU_DIST_INCLUDE_DIR  - Location of SuperLU_DIST headers
#     SuperLU_DIST_LIBRARIES    - SuperLU_DIST libraries
#
#   This module will export the following targets if SuperLU_FOUND
#
#     SuperLU::SuperLU_DIST
#
#
#
#
#   Proper usage:
#
#     project( TEST_FIND_SuperLUDIST C )
#     find_package( SuperLU_DIST )
#
#     if( SuperLU_DIST_FOUND )
#       add_executable( test test.cxx )
#       target_link_libraries( test SuperLU::SuperLU_DIST )
#     endif()
#
#
#
#
#   This module will use the following variables to change
#   default behaviour if set
#
#     SuperLU_DIST_PREFIX
#     SuperLU_DIST_INCLUDE_DIR
#     SuperLU_DIST_LIBRARY_DIR
#     SuperLU_DIST_LIBRARIES
#
#
#   This module also calls FindLinAlg.cmake if no LinAlg::BLAS
#   TARGET is found. If this behaviour is not desired, ensure 
#   that there is a proper definition of LinAlg::BLAS prior
#   to invokation by either calling find_package( LinAlg ) 
#   or creating a user defined target which properly links to
#   blas
#
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

include( CMakePushCheckState )
include( CheckLibraryExists )
include( CheckSymbolExists )
include( FindPackageHandleStandardArgs )

include( ${CMAKE_CURRENT_LIST_DIR}/util/CommonFunctions.cmake )
fill_out_prefix( SuperLU_DIST )



# Dependencies
include(CMakeFindDependencyMacro)
# Inherits MPI from ParMETIS
find_dependency( BLAS )
find_dependency( ParMETIS )



# Try to find headers
find_path( SuperLU_DIST_INCLUDE_DIR
  NAMES superlu_defs.h
  HINTS ${SuperLU_DIST_PREFIX}
  PATHS ${SuperLU_DIST_INCLUDE_DIR}
  PATH_SUFFIXES include
  DOC "Local of SuperLU_DIST Headers"
)



# Try to find libraries if not already set
if( NOT SuperLU_DIST_LIBRARIES )

  find_library( SuperLU_DIST_LIBRARIES
    NAMES superlu_dist
    HINTS ${SuperLU_DIST_PREFIX}
    PATHS ${SuperLU_DIST_LIBRARY_DIR}
    PATH_SUFFIXES lib lib64 lib32
    DOC "SuperLU_DIST Libraries"
  )

else()
  foreach( _lib ${SuperLU_DIST_LIBRARIES} )
    if( _lib MATCHES "ParMETIS" )
      find_dependency(ParMETIS)
    if( _lib MATCHES "BLAS" )
      find_dependency(BLAS)
    endif()
    endif()
  endforeach()
endif()

if( SuperLU_DIST_LIBRARIES )
  list( APPEND SuperLU_DIST_LIBRARIES BLAS::BLAS ParMETIS::ParMETIS )
endif()

# Check version
if( EXISTS ${SuperLU_DIST_INCLUDE_DIR}/superlu_defs.h )
  set( version_pattern 
  "^#define[\t ]+SuperLU_DIST_(MAJOR|MINOR|PATCH)_VERSION[\t ]+([0-9\\.]+)$"
  )
  file( STRINGS ${SuperLU_DIST_INCLUDE_DIR}/superlu_defs.h SuperLU_DIST_version
        REGEX ${version_pattern} )
  
  foreach( match ${SuperLU_DIST_version} )
  
    if(SuperLU_DIST_VERSION_STRING)
      set(SuperLU_DIST_VERSION_STRING "${SuperLU_DIST_VERSION_STRING}.")
    endif()
  
    string(REGEX REPLACE ${version_pattern} 
      "${SuperLU_DIST_VERSION_STRING}\\2" 
      SuperLU_DIST_VERSION_STRING ${match}
    )
  
    set(SuperLU_DIST_VERSION_${CMAKE_MATCH_1} ${CMAKE_MATCH_2})
  
  endforeach()
  
  unset( SuperLU_DIST_version )
  unset( version_pattern )
endif()



# Determine if we've found SuperLU_DIST
mark_as_advanced( SuperLU_DIST_FOUND SuperLU_DIST_INCLUDE_DIR SuperLU_DIST_LIBRARIES )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args( SuperLU_DIST
  REQUIRED_VARS SuperLU_DIST_LIBRARIES SuperLU_DIST_INCLUDE_DIR
  VERSION_VAR SuperLU_DIST_VERSION_STRING
)

set( SuperLU_DIST_LIBRARIES "${SuperLU_DIST_LIBRARIES}" CACHE STRING "SuperLU_DIST LIBRARIES" FORCE )
set( SuperLU_DIST_INCLUDE_DIR "${SuperLU_DIST_INCLUDE_DIR}" CACHE STRING "SuperLU_DIST INCLUDE_DIR" FORCE )

# Export target
if( SuperLU_DIST_FOUND AND NOT TARGET SuperLU::SuperLU_DIST )

  add_library( SuperLU::SuperLU_DIST INTERFACE IMPORTED )
  set_target_properties( SuperLU::SuperLU_DIST PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${SuperLU_DIST_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES      "${SuperLU_DIST_LIBRARIES}" 
  )

endif()
