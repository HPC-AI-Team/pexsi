#!/bin/bash
# Example file to build parmetis. copy this file to the home directory of
# parmetis and make necessary modifications.
# refer to BUILD.txt and Install.txt under the home directory for more
# information.

INSTALL_DIR=/home/linlin/Software/parmetis-4.0.3_install
CC=mpicc
CXX=mpicxx

# Build METIS. This is important. Otherwise ParMETIS would not build.
cd metis
make config prefix=${INSTALL_DIR} 
make -j 12
make install
cd ..
echo "Finished building metis"

# Build ParMETIS
make config prefix=${INSTALL_DIR} cc=${CC} cxx=${CXX}
make -j 12
make install
echo "Finished building parmetis"
