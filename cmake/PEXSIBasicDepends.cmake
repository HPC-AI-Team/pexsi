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


# Handle MPI
find_package( MPI REQUIRED )

include_directories( ${MPI_CXX_INCLUDE_PATH}    )
include_directories( ${MPI_C_INCLUDE_PATH}      )
include_directories( ${MPI_Fortan_INCLUDE_PATH} )

list( APPEND PEXSI_EXT_LINK ${MPI_CXX_LIBRARIES}     )
list( APPEND PEXSI_EXT_LINK ${MPI_C_LIBRARIES}       )
list( APPEND PEXSI_EXT_LINK ${MPI_Fortran_LIBRARIES} )


if( MPI_CXX_COMPILE_FLAGS )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
endif( MPI_CXX_COMPILE_FLAGS )

if( MPI_C_COMPILE_FLAGS )
  set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS}")
endif( MPI_C_COMPILE_FLAGS )

if( MPI_Fortran_COMPILE_FLAGS )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS}")
endif( MPI_Fortran_COMPILE_FLAGS )






# Handle OpenMP
if( PEXSI_ENABLE_OPENMP )

  find_package( OpenMP REQUIRED )

  add_definitions( "-D OPENMP" ) # FIXME: Redundant with _OPENMP

  set( CMAKE_CXX_FLAGS     "${CMAKE_CXX_FLAGS}     ${OpenMP_CXX_FLAGS}"     )
  set( CMAKE_C_FLAGS       "${CMAKE_C_FLAGS}       ${OpenMP_C_FLAGS}"       )
  set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}" )

endif( PEXSI_ENABLE_OPENMP )
