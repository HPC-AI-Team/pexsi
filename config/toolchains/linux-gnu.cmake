# Typical Build Environment for GNU+Linux systems with GNU compilers
set( CMAKE_C_COMPILER       mpicc    )
set( CMAKE_CXX_COMPILER     mpicxx   )
set( CMAKE_Fortran_COMPILER mpif90   )

set( SuperLU_DIST_PREFIX    "/home/linlin/Software/SuperLU_DIST_install/v7.2.0" )
set( ParMETIS_PREFIX        "/home/linlin/Software/parmetis-4.0.3_install" )
set( BLAS_PREFIX            "/home/linlin/Software/OpenBLAS_install" )
set( LAPACK_PREFIX          "${BLAS_PREFIX}" ) # OpenBLAS has a LAPACK linker
