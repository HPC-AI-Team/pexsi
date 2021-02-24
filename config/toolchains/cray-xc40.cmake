# Typical Build Environment for Cray systems with Cray compiler wrappers (GNU,Intel,PGI)
set( CMAKE_C_COMPILER       cc   )
set( CMAKE_CXX_COMPILER     CC   )
set( CMAKE_Fortran_COMPILER ftn  )

set( SuperLU_DIST_PREFIX    "/path/to/superlu/install"  )
set( ParMETIS_PREFIX        "/path/to/parmetis/install" )

# BLAS and LAPACK are auto deduced
