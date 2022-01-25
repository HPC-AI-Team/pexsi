#!/usr/bin/bash
COMPILE_MODE     = release
USE_PROFILE      = 0
#PAR_ND_LIBRARY   = ptscotch
#SEQ_ND_LIBRARY   = scotch
PAR_ND_LIBRARY   = parmetis
SEQ_ND_LIBRARY   = metis
USE_SYMPACK      = 0
USE_OPENMP       = 0
USE_COREDUMPER   = 0

# Different compiling and linking options.
SUFFIX       = linux_release_v2.0

# Compiler and tools
################################################################
CC           = mpicc 
CXX          = mpic++
FC           = mpif90
LOADER       = mpic++

AR           = ar
ARFLAGS      = rvcu
# For System V based machine without ranlib, like Cray and SGI,
# use touch instead.
#RANLIB      = touch
RANLIB       = ranlib

CP           = cp
RM           = rm
RMFLAGS      = -f
################################################################

ifeq (${USE_OPENMP}, 1)
  OPENMP_DEF   = -DOPENMP
	OPENMP_FLAG  = -fopenmp
endif


# PEXSI directory
PEXSI_DIR       = $(HOME)/Projects/pexsi
PEXSI_BUILD_DIR = $(PEXSI_DIR)/build_make

# Required libraries directories
SUPERLU_DIST_DIR  = $(HOME)/Software/SuperLU_DIST_install/v7.2.0
METIS_DIR     = $(HOME)/Software/parmetis-4.0.3_install
PARMETIS_DIR  = $(HOME)/Software/parmetis-4.0.3_install
PTSCOTCH_DIR  = $(HOME)/Software/scotch_6.0.0
SCOTCH_DIR = ${PTSCOTCH_DIR}
LAPACK_DIR    = /usr/lib/x86_64-linux-gnu
BLAS_DIR      = /usr/lib/x86_64-linux-gnu

ifeq (${USE_COREDUMPER}, 1)
  COREDUMPER_DIR = /home/lin/Software/coredumper-1.2.1/build
  COREDUMPER_INCLUDE = -I${COREDUMPER_DIR}/include
  COREDUMPER_LIB   = ${COREDUMPER_DIR}/lib/libcoredumper.a
endif


# Includes
PEXSI_INCLUDE    = -I${PEXSI_DIR}/include 
SUPERLU_DIST_INCLUDE = -I${SUPERLU_DIST_DIR}/include
PARMETIS_INCLUDE = -I${PARMETIS_DIR}/include
INCLUDES         = ${PEXSI_INCLUDE} ${SUPERLU_DIST_INCLUDE} ${PARMETIS_INCLUDE} ${COREDUMPER_INCLUDE}

# Libraries
CPP_LIB          = -lstdc++ -lmpi -lmpi_cxx
GFORTRAN_LIB     = -lgfortran
LAPACK_LIB       = -llapack
BLAS_LIB         = -lblas
SUPERLU_DIST_LIB     = ${SUPERLU_DIST_DIR}/lib/libsuperlu_dist.a
PEXSI_LIB        = ${PEXSI_DIR}/src/libpexsi_${SUFFIX}.a

# Graph partitioning libraries
METIS_LIB        = -L${METIS_DIR}/lib -lmetis
PARMETIS_LIB     = -L${PARMETIS_DIR}/lib -lparmetis 
SCOTCH_LIB       = -L${PTSCOTCH_DIR}/lib -lscotchmetis -lscotch -lscotcherr
PTSCOTCH_LIB     = -L${PTSCOTCH_DIR}/lib -lptscotchparmetis -lptscotch -lptscotcherr -lscotch

# Different compiling and linking options.
ifeq (${COMPILE_MODE}, release)
  ifeq (${USE_COREDUMPER}, 1)
    COMPILE_DEF    = -DRELEASE -DCOREDUMPER
  else
    COMPILE_DEF    = -DRELEASE
  endif
  COMPILE_FLAG   = -O3 -w -g  
endif
ifeq (${COMPILE_MODE}, debug)
  ifeq (${USE_COREDUMPER}, 1)
    COMPILE_DEF    = -DDEBUG=1 -DCOREDUMPER
  else
    COMPILE_DEF    = -DDEBUG=1
  endif
  COMPILE_FLAG   = -O2 -w -g
endif

ifeq (${PAR_ND_LIBRARY}, ptscotch)
  PAR_ND_LIB = ${PTSCOTCH_LIB}
else
  PAR_ND_LIB = ${PARMETIS_LIB}
endif 

ifeq (${SEQ_ND_LIBRARY}, scotch)
  SEQ_ND_LIB = ${SCOTCH_LIB}
else
  SEQ_ND_LIB = ${METIS_LIB}
endif 

ifeq (${USE_PROFILE}, 1)
  PROFILE_FLAG  = -DPROFILE
endif


LIBS  = ${PEXSI_LIB} ${SUPERLU_DIST_LIB} ${PAR_ND_LIB} ${SEQ_ND_LIB} ${LAPACK_LIB} ${BLAS_LIB} ${COREDUMPER_LIB} ${GFORTRAN_LIB} 
COMPILE_DEF  += -DAdd_ #-D_MIRROR_RIGHT_
CPPFLAG = -std=c++11

ifeq (${USE_SYMPACK}, 1)
  #symPACK related definitions
  SYMPACK_DIR = /home/lin/Projects/sympack/release
  include ${SYMPACK_DIR}/include/sympack.mak
  CPPFLAG += ${SYMPACK_INCLUDE} 
  LIBS+= ${SYMPACK_LIB} ${LAPACK_LIB} ${BLAS_LIB} ${GFORTRAN_LIB}
  COMPILE_DEF  += -DWITH_SYMPACK
  CPPFLAG += -std=c++11
endif


CFLAGS       = ${COMPILE_FLAG} ${OPENMP_FLAG} ${PROFILE_FLAG} ${INCLUDES}  -std=c99
FFLAGS       = ${COMPILE_FLAG} ${OPENMP_FLAG} ${PROFILE_FLAG} ${INCLUDES} 
CXXFLAGS     = ${COMPILE_FLAG} ${OPENMP_FLAG} ${CPPFLAG} ${PROFILE_FLAG} ${INCLUDES}  
CCDEFS       = ${COMPILE_DEF} ${OPENMP_DEF}
CPPDEFS      = ${COMPILE_DEF} ${OPENMP_DEF}
LOADOPTS     = ${PROFILE_FLAG} ${OPENMP_FLAG} ${LIBS} -Wl,--allow-multiple-definition
FLOADOPTS    = ${PROFILE_FLAG} ${OPENMP_FLAG} ${LIBS} ${CPP_LIB} -Wl,--allow-multiple-definition

# Generate auto-dependencies 
%.d: %.c
	@set -e; rm -f $@; \
	$(CC) -M $(CCDEFS) $(CFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@;\
	rm -f $@.$$$$

%.d: %.cpp
	@set -e; rm -f $@; \
	$(CXX) -M $(CPPDEFS) $(CXXFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@;\
	rm -f $@.$$$$
