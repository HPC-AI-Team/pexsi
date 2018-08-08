# - Try to find PTSCOTCH
# Once done this will define
#
#  PTSCOTCH_FOUND        - system has PTSCOTCH
#  PTSCOTCH_INCLUDE_DIRS - include directories for PTSCOTCH
#  PTSCOTCH_LIBRARIES    - libraries for PTSCOTCH
#
# Variables used by this module. They can change the default behaviour and
# need to be set before calling find_package:
#
#  PTSCOTCH_DIR          - Prefix directory of the PTSCOTCH installation
#  PTSCOTCH_INCLUDE_DIR  - Include directory of the PTSCOTCH installation
#                          (set only if different from ${PTSCOTCH_DIR}/include)
#  PTSCOTCH_LIB_DIR      - Library directory of the PTSCOTCH installation
#                          (set only if different from ${PTSCOTCH_DIR}/lib)
#  PTSCOTCH_LIB_SUFFIX   - Also search for non-standard library names with the
#                          given suffix appended

#=============================================================================
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
#=============================================================================

if(NOT PTSCOTCH_INCLUDE_DIR)
  find_path(PTSCOTCH_INCLUDE_DIR ptscotch.h
    HINTS ${PTSCOTCH_INCLUDE_DIR} ENV PTSCOTCH_INCLUDE_DIR ${PTSCOTCH_DIR} ENV PTSCOTCH_DIR
    PATH_SUFFIXES include
    DOC "Directory where the PTSCOTCH header files are located"
  )
endif()

if(NOT SCOTCH_INCLUDE_DIR)
  find_path(SCOTCH_INCLUDE_DIR scotch.h
    HINTS ${SCOTCH_INCLUDE_DIR} ENV SCOTCH_INCLUDE_DIR ${SCOTCH_DIR} ENV SCOTCH_DIR
    PATH_SUFFIXES include
    DOC "Directory where the SCOTCH header files are located"
  )
endif()

if(PTSCOTCH_LIBRARIES)
  set(PTSCOTCH_LIBRARY ${PTSCOTCH_LIBRARIES})
endif()
if(NOT PTSCOTCH_LIBRARY)
  find_library(PTSCOTCH_LIBRARY
    NAMES ptscotch ptscotch${PTSCOTCH_LIB_SUFFIX}
    HINTS ${PTSCOTCH_LIB_DIR} ENV PTSCOTCH_LIB_DIR ${PTSCOTCH_DIR} ENV PTSCOTCH_DIR
    PATH_SUFFIXES lib
    DOC "Directory where the PTSCOTCH library is located"
  )
  find_library(PTSCOTCHPARMETIS_LIBRARY
    NAMES ptscotchparmetis ptscotchparmetis${PTSCOTCH_LIB_SUFFIX}
    HINTS ${PTSCOTCH_LIB_DIR} ENV PTSCOTCH_LIB_DIR ${PTSCOTCH_DIR} ENV PTSCOTCH_DIR
    PATH_SUFFIXES lib
    DOC "Directory where the PTSCOTCH library is located"
  )
  find_library(PTSCOTCHERR_LIBRARY
    NAMES ptscotcherr ptscotcherr${PTSCOTCH_LIB_SUFFIX}
    HINTS ${PTSCOTCH_LIB_DIR} ENV PTSCOTCH_LIB_DIR ${PTSCOTCH_DIR} ENV PTSCOTCH_DIR
    PATH_SUFFIXES lib
    DOC "Directory where the PTSCOTCH library is located"
  )

  set( PTSCOTCH_LIBRARY ${PTSCOTCHPARMETIS_LIBRARY} ${PTSCOTCH_LIBRARY} ${PTSCOTCHERR_LIBRARY} )
endif()

if(SCOTCH_LIBRARIES)
  set(SCOTCH_LIBRARY ${SCOTCH_LIBRARIES})
endif()
if(NOT SCOTCH_LIBRARY)
  find_library(SCOTCH_LIBRARY
    NAMES scotch scotch${PTSCOTCH_LIB_SUFFIX}
    HINTS ${PTSCOTCH_LIB_DIR} ENV PTSCOTCH_LIB_DIR ${PTSCOTCH_DIR} ENV PTSCOTCH_DIR
    PATH_SUFFIXES lib
    DOC "Directory where the SCOTCH library is located"
  )
  find_library(SCOTCHERR_LIBRARY
    NAMES scotcherr scotcherr${PTSCOTCH_LIB_SUFFIX}
    HINTS ${PTSCOTCH_LIB_DIR} ENV PTSCOTCH_LIB_DIR ${PTSCOTCH_DIR} ENV PTSCOTCH_DIR
    PATH_SUFFIXES lib
    DOC "Directory where the SCOTCH library is located"
  )

  set( SCOTCH_LIBRARY ${SCOTCH_LIBRARY} ${SCOTCHERR_LIBRARY} )
endif()


# Standard package handling
include(FindPackageHandleStandardArgs)
if(CMAKE_VERSION VERSION_GREATER 2.8.2)
  find_package_handle_standard_args(PTSCOTCH
    REQUIRED_VARS PTSCOTCH_LIBRARY PTSCOTCH_INCLUDE_DIR 
    VERSION_VAR PTSCOTCH_VERSION_STRING)
else()
  find_package_handle_standard_args(PTSCOTCH
    REQUIRED_VARS PTSCOTCH_LIBRARY PTSCOTCH_INCLUDE_DIR)
endif()

if(PTSCOTCH_FOUND)
  set(PTSCOTCH_LIBRARIES ${PTSCOTCH_LIBRARY} ${SCOTCH_LIBRARY})
  set(PTSCOTCH_INCLUDE_DIRS ${PTSCOTCH_INCLUDE_DIR} ${SCOTCH_INCLUDE_DIR})
else()
  unset(SCOTCH_LIBRARY CACHE)
  unset(SCOTCH_INCLUDE_DIR CACHE)
endif()

mark_as_advanced(PTSCOTCH_INCLUDE_DIR SCOTCH_INCLUDE_DIR
  PTSCOTCH_LIBRARY SCOTCH_LIBRARY)
