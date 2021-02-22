#!/bin/bash
# Example file to build parmetis. copy this file to the home directory of
# superlu_dist and make necessary modifications.
# refer to README.md under the home directory for more information.

INSTALL_DIR=/home/linlin/Software/SuperLU_DIST_install/v6.4.0
SuperLU_DIR=/home/linlin/Software/SuperLU_DIST_6.4.0
BLAS_LIB=/usr/lib/x86_64-linux-gnu/libblas.a
LAPACK_LIB=/usr/lib/x86_64-linux-gnu/liblapack.a
PARMETIS_DIR=$HOME/Software/parmetis-4.0.3_install
METIS_DIR=$HOME/Software/parmetis-4.0.3_install
GFORTRAN_LIB=/usr/lib/gcc/x86_64-linux-gnu/7/libgfortran.a
CC=mpicc
CXX=mpicxx

BDIR=build_cpu
rm -rf ${BDIR}
cmake -H. -B${BDIR} \
  -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
  -DTPL_PARMETIS_INCLUDE_DIRS="${PARMETIS_DIR}/include;${METIS_DIR}/include" \
  -DTPL_PARMETIS_LIBRARIES="${PARMETIS_DIR}/lib/libparmetis.a;${METIS_DIR}/lib/libmetis.a" \
  -DCMAKE_C_FLAGS="-std=c99 -g -DPRNTlevel=0 -DDEBUGlevel=0" \
  -DCMAKE_C_COMPILER=$CC \
  -DCMAKE_CXX_COMPILER=$CXX \
  -DCMAKE_CXX_FLAGS="-std=c++11" \
  -DTPL_ENABLE_BLASLIB=OFF \
  -DTPL_BLAS_LIBRARIES="-lblas" \
  -DBUILD_SHARED_LIBS=OFF \
  -Denable_openmp=FALSE 
cd ${BDIR}
make -j 12
make install
cd ..
