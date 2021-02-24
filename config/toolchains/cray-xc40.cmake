# Typical Build Environment for Cray systems with Cray compiler wrappers (GNU,Intel,PGI)
set( CMAKE_C_COMPILER       cc   )
set( CMAKE_CXX_COMPILER     CC   )
set( CMAKE_Fortran_COMPILER ftn  )

set( SuperLU_DIST_PREFIX    "/global/project/projectdirs/m1027/shared_libraries_cori_intel/SuperLU_DIST_6.4.0/install"  )
set( ParMETIS_PREFIX        "/project/projectdirs/m1027/shared_libraries_cori_intel/parmetis_v0.4.3/install" )

# BLAS and LAPACK are auto deduced
