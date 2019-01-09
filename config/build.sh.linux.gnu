CC=mpicc
CXX=mpicxx
FC=mpif90
PEXSI_INSTALL_DIR=/home/lin/Projects/pexsi/install
BLAS_LIB=$HOME/Software/OpenBLAS/build_release/lib/libopenblas.a
LAPACK_LIB=$HOME/Software/lapack-3.5.0/liblapack.a
DSUPERLU_DIR=$HOME/Software/SuperLU_DIST_6.1.0/release
PARMETIS_DIR=$HOME/Software/parmetis-4.0.3/build_release
METIS_DIR=$HOME/Software/metis-5.1.0/build_release
GFORTRAN_LIB=/usr/lib/gcc/x86_64-linux-gnu/4.8/libgfortran.a


cmake \
  -DCMAKE_INSTALL_PREFIX=$PEXSI_INSTALL_DIR \
  -DCMAKE_C_COMPILER=$CC \
  -DCMAKE_CXX_COMPILER=$CXX \
  -DCMAKE_Fortran_COMPILER=$FC \
  -Dlinalg_BLAS_LIBRARIES="$BLAS_LIB;-fopenmp;-lm;$GFORTRAN_LIB" \
  -Dlinalg_LAPACK_LIBRARIES=$LAPACK_LIB \
  -Dsuperlu_dist_INCLUDE_DIR=$DSUPERLU_DIR/include \
  -Dsuperlu_dist_LIBRARIES=$DSUPERLU_DIR/lib/libsuperlu_dist.a \
  -Dmetis_INCLUDE_DIR=$METIS_DIR/include \
  -Dmetis_LIBRARIES=$METIS_DIR/lib/libmetis.a \
  -Dparmetis_INCLUDE_DIR=$PARMETIS_DIR/include \
  -Dparmetis_LIBRARIES=$PARMETIS_DIR/lib/libparmetis.a \
  ..
