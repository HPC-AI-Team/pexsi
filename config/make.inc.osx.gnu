#/usr/bin/bash
COMPILE_MODE     = debug

# Different compiling and linking options.
SUFFIX       = osx_v0.5.5

#compiler, options and tools
################################################################
CC           = mpicc 
CXX          = mpic++
FC           = mpif90
LOADER       = mpic++

DEBUG_PRINT_LEVEL    = 1
DEBUG_COMPILE_FLAG   = -O0 -w -g
RELEASE_COMPILE_FLAG = -O3 -w 

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


#pexsi directory
PEXSI_DIR     = $(HOME)/Work/postdoc-lbl/pexsi

#required libraries directories
DSUPERLU_DIR  = $(HOME)/software/release/SuperLU_DIST_3.3

#useful on some systems (OSX using homebrew compilers)
GFORTRAN_DIR  = /usr/local/Cellar/gfortran/4.8.2/gfortran/lib

# Graph partitioning libraries
#Metis or scotch for sequential Nested Disection
#SEQ_ND_DIR  = $(HOME)/software/release/parmetis-4.0.3/build/Darwin-x86_64
#SEQ_ND_LIB     = -L${SEQ_ND_DIR}/libmetis -lmetis
SEQ_ND_DIR  = $(HOME)/software/release/scotch_6.0.0
SEQ_ND_LIB     = -L${SEQ_ND_DIR}/lib -lscotchmetis -lscotch -lscotcherr
#ParMetis or PTScotch for parallel Nested Dissection
#PAR_ND_DIR  = $(HOME)/software/release/parmetis-4.0.3/build/Darwin-x86_64
#PAR_ND_LIB     = -L${PAR_ND_DIR}/libparmetis -lparmetis 
PAR_ND_DIR  = $(HOME)/software/release/scotch_6.0.0
PAR_ND_LIB     = -L${PAR_ND_DIR}/lib -lptscotchparmetis -lptscotch -lptscotcherr -lscotch

# Includes
PEXSI_INCLUDE    = -I${PEXSI_DIR}/include 
DSUPERLU_INCLUDE = -I${DSUPERLU_DIR}/SRC
INCLUDES         = ${PEXSI_INCLUDE} ${DSUPERLU_INCLUDE} ${NGCHOL_INCLUDE} 

# Libraries
LAPACK_LIB       = -llapack
BLAS_LIB         = -lblas
DSUPERLU_LIB     = ${DSUPERLU_DIR}/lib/libsuperlu_dist_3.3.a
PEXSI_LIB        = ${PEXSI_DIR}/lib/libpexsi_${SUFFIX}.a
GFORTRAN_LIB     = -L${GFORTRAN_DIR} -lgfortran


################ End of configuration #####################


# Different compiling and linking options.
ifeq (${COMPILE_MODE}, release)
  COMPILE_DEF    = -DRELEASE
  COMPILE_FLAG   = $(RELEASE_COMPILE_FLAG)
endif
ifeq (${COMPILE_MODE}, debug)
  COMPILE_DEF    = -DDEBUG=$(DEBUG_PRINT_LEVEL)
  COMPILE_FLAG   = $(DEBUG_COMPILE_FLAG)
endif

LIBS  = ${PEXSI_LIB} ${NGCHOL_LIB} ${DSUPERLU_LIB} ${PAR_ND_LIB} ${SEQ_ND_LIB} ${LAPACK_LIB} ${BLAS_LIB} ${GFORTRAN_LIB}

COMPILE_DEF += -DAdd_ -std=c++11

CFLAGS       = ${COMPILE_FLAG} ${INCLUDES}
FFLAGS       = ${COMPILE_FLAG} ${INCLUDES}
CXXFLAGS     = ${COMPILE_FLAG} ${INCLUDES} 
CCDEFS       = ${COMPILE_DEF} 
CPPDEFS      = ${COMPILE_DEF} 
LOADOPTS     = ${LIBS}

.PHONY: ${PEXSI_LIB}
