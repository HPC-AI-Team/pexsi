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


# Compile

set( PEXSI_SRC_C   pdsymbfact.c pzsymbfact.c get_perm_c_parmetis.c )
set( PEXSI_SRC_CXX interface.cpp ppexsi.cpp pole.cpp TreeBcast.cpp
                   superlu_dist_internal_complex.cpp 
                   superlu_dist_internal_real.cpp 
                   mpi_interf.cpp lapack.cpp blas.cpp utility.cpp 
                   global.cpp timer.cpp getPole.cpp )


if( PEXSI_ENABLE_FORTRAN )
  message( STATUS "PEXSI Will Enable FORTRAN Bindings" )
  set( PEXSI_SRC_F90 f_interface.f90 )
endif()

add_library( pexsi STATIC ${PEXSI_SRC_C} ${PEXSI_SRC_CXX} ${PEXSI_SRC_F90})
add_library( PEXSI::pexsi ALIAS pexsi )

# Set up Fortran module directory
if( PEXSI_ENABLE_FORTRAN )
  set_property(
    TARGET pexsi
    PROPERTY Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/include
  )
endif()



# TARGET properties

target_link_libraries( pexsi
  PUBLIC
    PEXSI::compile_definitions
    PEXSI::parallel_cxx
    PEXSI::parallel_c
    PEXSI::tpl_depends
)

if( PEXSI_ENABLE_FORTRAN )
  target_link_libraries( pexsi PUBLIC PEXSI::parallel_fortran )
endif()

target_include_directories( pexsi  
  PUBLIC
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
)



# INSTALL rules

include( GNUInstallDirs )
set( INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/PEXSI )

# Targets
install(TARGETS pexsi
    EXPORT pexsi-targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

set_target_properties( pexsi PROPERTIES EXPORT_NAME PEXSI )

# Install Headers
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} )

if( PEXSI_ENABLE_FORTRAN )
  install(DIRECTORY ${PROJECT_BINARY_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} )
endif()

# Export target to script
install(EXPORT pexsi-targets
  FILE
    PEXSITargets.cmake
  NAMESPACE
    PEXSI::
  DESTINATION
    ${INSTALL_CONFIGDIR}
)


#Create a ConfigVersion.cmake file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/PEXSIConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

# Setup PEXSIConfig.cmake
configure_package_config_file(${PROJECT_SOURCE_DIR}/cmake/PEXSIConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/PEXSIConfig.cmake
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

#Install the config, configversion and custom find modules
install(DIRECTORY 
  ${PROJECT_SOURCE_DIR}/cmake/ 
  DESTINATION ${INSTALL_CONFIGDIR} 
  FILES_MATCHING PATTERN "*.cmake"
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/PEXSIConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/PEXSIConfigVersion.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)
