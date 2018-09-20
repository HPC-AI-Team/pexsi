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


add_library( PEXSI::parmetis INTERFACE IMPORTED )


# Try to find METIS / ParMETIS
find_package( ParMETIS QUIET )

if( PARMETIS_FOUND )
  # If we found METIS/ParMETIS, set vars

  target_link_libraries( PEXSI::parmetis INTERFACE ParMETIS::parmetis )

else()
  # If not, create a target

  include(ExternalProject)

  message( STATUS "Opting to build METIS/ParMETIS" )

  set( PARMETIS_PREFIX        ${PROJECT_BINARY_DIR}/external/parmetis )
  set( PARMETIS_INCLUDE_DIRS  ${PARMETIS_PREFIX}/include              )
  set( PARMETIS_LIBRARY_DIRS  ${PARMETIS_PREFIX}/lib                  )
  set( PARMETIS_LIBRARIES     ${PARMETIS_LIBRARY_DIRS}/libparmetis.a;${PARMETIS_LIBRARY_DIRS}/libmetis.a )

  file( MAKE_DIRECTORY ${PARMETIS_INCLUDE_DIRS} ) # Passify TARGET import

  ExternalProject_Add(parmetis
    PREFIX ${PARMETIS_PREFIX}
    URL http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz
    PATCH_COMMAND patch < ${PEXSI_PATCH_PATH}/install_metis.patch
    CMAKE_ARGS
      -D GKLIB_PATH=${PARMETIS_PREFIX}/src/parmetis/metis/GKlib
      -D METIS_PATH=${PARMETIS_PREFIX}/src/parmetis/metis
      -D CMAKE_INSTALL_PREFIX=${PARMETIS_PREFIX}
      -D CMAKE_C_COMPILER=${MPI_C_COMPILER}
      -D CMAKE_CXX_COMPILER=${MPI_CXX_COMPILER}
      -D METIS_INSTALL=ON
  )
  
  add_dependencies( PEXSI::parmetis parmetis )
  set_target_properties( PEXSI::parmetis PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${PARMETIS_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES      "${PARMETIS_LIBRARIES}"
  )
  

endif()
