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


# Turn off system search if TPL linear algebra has
# been specified.
if( TPL_BLAS_LIBRARIES OR TPL_LAPACK_LIBRARIES )
  set( PEXSI_USE_SYSTEM_LINALG OFF )
endif()



if( PEXSI_USE_SYSTEM_LINALG )

  # Use default CMake system to find BLAS/LAPACK
  find_package( BLAS    REQUIRED )
  find_package( LAPACK  REQUIRED )

else()

  # Use user-specified TPL_BLAS_LIBRARIES
  if( TPL_BLAS_LIBRARIES )
    set( BLAS_LIBRARIES ${TPL_BLAS_LIBRARIES} )
  endif()

  # Use user-specified TPL_LAPACK_LIBRARIES
  if( TPL_LAPACK_LIBRARIES )
    set( LAPACK_LIBRARIES ${TPL_LAPACK_LIBRARIES} )
  endif()

endif()

# Link to the libraries if specified
if( BLAS_LIBRARIES )
  message( STATUS "Using BLAS_LIBRARIES = ${BLAS_LIBRARIES}" )
  list( APPEND PEXSI_EXT_LINK ${BLAS_LIBRARIES} )
endif()

if( LAPACK_LIBRARIES )
  message( STATUS "Using LAPACK_LIBRARIES = ${LAPACK_LIBRARIES}" )
  list( APPEND PEXSI_EXT_LINK ${LAPACK_LIBRARIES} )
endif()


file( WRITE ${CMAKE_CURRENT_BINARY_DIR}/linalg_underscore.c 
  "void dgemm_(); int main() { dgemm_(); return 0; }" )
file( WRITE ${CMAKE_CURRENT_BINARY_DIR}/linalg_no_underscore.c 
  "void dgemm(); int main() { dgemm(); return 0; }" )




message( STATUS "Perfoming Test BLAS_USE_UNDERSCORE" )

try_compile(
  BLAS_USE_UNDERSCORE
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}/linalg_underscore.c
  LINK_LIBRARIES ${PEXSI_EXT_LINK}
)

if( BLAS_USE_UNDERSCORE )
  message( STATUS "Perfoming Test BLAS_USE_UNDERSCORE -- Yes" )
  add_definitions( "-DAdd_" )
else()
  message( STATUS "Perfoming Test BLAS_USE_UNDERSCORE -- No" )
endif()




message( STATUS "Perfoming Test BLAS_NO_UNDERSCORE" )

try_compile(
  BLAS_NO_UNDERSCORE
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}/linalg_no_underscore.c
  LINK_LIBRARIES ${PEXSI_EXT_LINK}
)

if( BLAS_NO_UNDERSCORE )
  message( STATUS "Perfoming Test BLAS_NO_UNDERSCORE -- Yes" )
else()
  message( STATUS "Perfoming Test BLAS_NO_UNDERSCORE -- No" )
endif()


if( (NOT BLAS_USE_UNDERSCORE) AND (NOT BLAS_NO_UNDERSCORE) )
  message( FATAL_ERROR "No Suitible BLAS library found" )
endif()




file( REMOVE ${CMAKE_CURRENT_BINARY_DIR}/linalg_underscore.c    )
file( REMOVE ${CMAKE_CURRENT_BINARY_DIR}/linalg_no_underscore.c )
