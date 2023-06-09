#   Copyright (c) 2018 The Regents of the University of California,
#   through Lawrence Berkeley National Laboratory.  
#
#   Author: David Williams-Young
#   
#   This file is part of PEXSI. All rights reserved.
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
cmake_minimum_required( VERSION 3.17 ) # Require CMake 3.17+

# Set up project definition + version information
project( PEXSI VERSION 2.0.0 LANGUAGES CXX C Fortran )


# PEXSI Options
option( PEXSI_ENABLE_PROFILE    "Enable Performance Profiling"        OFF )
#option( PEXSI_ENABLE_SYMPACK    "Enable interface to symPACK"         OFF )
option( PEXSI_ENABLE_OPENMP     "Enable OpenMP Bindings"              OFF )
option( PEXSI_ENABLE_FORTRAN    "Enable Fortran Bindings"             ON  )


# Append local cmake directory to find CMAKE Modules
if( CMAKE_MODULE_PATH )
  list( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
else()
  set( CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
endif()
list( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")

# Global config variables
set( PEXSI_PATCH_PATH ${PROJECT_SOURCE_DIR}/external/patch )

# Default to Release if not subproject
if( NOT CMAKE_BUILD_TYPE AND (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME))

  message( STATUS "Using Default CMAKE_BUILD_TYPE=Release" )
  set( CMAKE_BUILD_TYPE "Release" CACHE STRING "CMake Build Type" FORCE )

endif()


#include( PEXSICompileFlags ) # Compile Flags
#include( PEXSIBasicDepends ) # Basic Dependencies
#include( PEXSITPLDepends   ) # TPL Dependencies


# libPEXSI
add_subdirectory( src )

if( PEXSI_ENABLE_FORTRAN )
  add_subdirectory( fortran )
endif()

# examples
add_subdirectory( examples )
