#!/bin/sh

# Create External Directory
depdir="$PWD/external"
patchdir="$depdir/patch"

CC=gcc
CXX=g++

MPIC=$HOME/local/mpich/3.2.1/gcc/8.1.1/bin/mpicc
MPICXX=$HOME/local/mpich/3.2.1/gcc/8.1.1/bin/mpicxx
MPIFC=$HOME/local/mpich/3.2.1/gcc/8.1.1/bin/mpifort

BLAS_LIB=/home/dbwy/local/openblas/git/gcc/8.1.1/lib/libopenblas.a
PARMETIS_LIB="$depdir/lib/libparmetis.a;$depdir/lib/libmetis.a"
PARMETIS_INC=$depdir/include




# METIS / ParMETIS

metisdir=$depdir
if [ ! -f $depdir/include/metis.h ] 
then

  echo "Building METIS/ParMETIS"

  # DL ParMETIS / METIS
  cd $metisdir
  wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.3.tar.gz
  tar xvf parmetis-4.0.3.tar.gz
  
  # APPLY Patch
  cd parmetis-4.0.3
  cp $patchdir/install_metis.patch .
  patch < install_metis.patch
  
  # Run CMAKE
  mkdir build_loc && cd build_loc
  cmake \
  	-DCMAKE_VERBOSE_MAKEFILE=1 \
  	-DGKLIB_PATH=$PWD/../metis/GKlib \
  	-DMETIS_PATH=$PWD/../metis \
  	-DCMAKE_INSTALL_PREFIX=$depdir \
  	-DCMAKE_C_COMPILER=$MPIC \
  	-DCMAKE_CXX_COMPILER=$MPICXX \
  	-DMETIS_INSTALL=ON \
  	..
  
  # BUILD
  make -j4 
  make install

else
  
  echo "Found METIS/ParMETIS Installation"

fi


# PT-SCOTCH
scotchdir=$depdir

if [ ! -f $depdir/include/ptscotch.h ]
then

  echo "Building SCOTCH / PT-SCOTCH"
  # DL SCOTCH / PT-SCOTCH
  cd $scotchdir
  wget https://gforge.inria.fr/frs/download.php/31831/scotch_6.0.0.tar.gz
  tar xvf scotch_6.0.0.tar.gz
  
  
  # Setup Makefile
  ptscotchinc=Makefile.inc.x86-64_pc_linux2
  cd scotch_6.0.0/src
  cp Make.inc/$ptscotchinc Makefile.inc
  
  # Specific to generic linux2 makefile
  sedCC="$CC"
  sedCC="${sedCC//\//\\/}"
  sedMPIC="$MPIC"
  sedMPIC="${sedMPIC//\//\\/}"
  
  sed -i -e "s/CCS		= gcc/CCS		= $sedCC/g" Makefile.inc
  sed -i -e "s/CCD		= gcc/CCD		= $sedCC/g" Makefile.inc
  sed -i -e "s/CCP		= mpicc/CCP		= $sedMPIC/g" Makefile.inc
  sed -i -e "s/LDFLAGS		= -lz -lm -lrt/LDFLAGS		= -lz -lm -lrt -lpthread/g" Makefile.inc

  sed -i -e "s/CFLAGS		= -O3 -DCOMMON_FILE_COMPRESS_GZ -DCOMMON_PTHREAD -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -DSCOTCH_PTHREAD -Drestrict=__restrict -DIDXSIZE64/CFLAGS		= -O3 -DCOMMON_FILE_COMPRESS_GZ -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -Drestrict=__restrict -DIDXSIZE64/g" Makefile.inc
  
  # Make / Install
  make ptscotch prefix=$depdir install
else

  echo "Found a SCOTCH / PT-SCOTCH"

fi






# SymPACK
# sympackdir=$depdir
# 
# if [ ! -f $depdir/include/sympack.hpp ]
# then
# 
#   echo "Building SymPACK"
#   
#   cd $sympackdir
#   git clone https://github.com/symPACK/symPACK.git
#   cd symPACK
#   
#   # Patch
#   cp $patchdir/gasnet_url.patch .
#   patch -p0 < gasnet_url.patch
#   
#   # Configure / Point SymPACK to SCOTCH/PT-SCOTCH
#   mkdir build && cd build
#   cmake \
#     -DCMAKE_BUILD_TYPE=Release \
#     -DCMAKE_INSTALL_PREFIX=$depdir \
#     -DSCOTCH_DIR=$depdir \
#     -DENABLE_SCOTCH=ON \
#     -DCMAKE_C_COMPILER=$MPIC \
#     -DCMAKE_CXX_COMPILER=$MPICXX \
#     -DCMAKE_Fortran_COMPILER=$MPIFC \
#     ..
#   
#   make -j4
#   make install
# 
# else
# 
#   echo "Found SymPACK Installation"
# 
# fi



# SuperLU-DIST
superludir=$depdir

if [ ! -f $depdir/include/superlu_defs.h ]
then

  echo "Building SuperLU_DIST"

  cd $superludir
  wget http://crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_dist_5.1.3.tar.gz
  tar xvf superlu_dist_5.1.3.tar.gz
  cd SuperLU_DIST_5.1.3

  # Apply Patches
  cp $patchdir/superlu* .
  patch -p0 < superlu_dmemory_print.patch
  patch -p0 < superlu_zmemory_print.patch
  
  mkdir build && cd build
  cmake \
    -DTPL_PARMETIS_INCLUDE_DIRS=$PARMETIS_INC \
    -DTPL_PARMETIS_LIBRARIES=$PARMETIS_LIB \
    -DTPL_BLAS_LIBRARIES=$BLAS_LIB \
    -Denable_blaslib=OFF \
    -DCMAKE_C_FLAGS="-std=c99 -O2" \
    -DCMAKE_C_COMPILER=$MPIC \
    -DCMAKE_Fortran_COMPILER=$MPIFC \
    -DCMAKE_INSTALL_PREFIX=$depdir \
    ..
  
  make -j4
  make install
else

  echo "Found SuperLI_DIST Installation"

fi
