#/usr/bin/bash
COMPILE_MODE     = release
USE_COMPLEX      = 0
USE_PROFILE      = 0
PAR_ND_LIBRARY   = ptscotch
SEQ_ND_LIBRARY   = scotch




# Different compiling and linking options.
SUFFIX       = osx_v0.5.5

#compiler and tools
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


#pexsi directory
PEXSI_DIR     = $(HOME)/Work/postdoc-lbl/pexsi

#required libraries directories
GFORTRAN_DIR  = /usr/local/lib
DSUPERLU_DIR  = $(HOME)/software/release/SuperLU_DIST_3.3
METIS_DIR     = $(HOME)/software/parmetis-4.0.3/build/Darwin-x86_64
PARMETIS_DIR  = $(HOME)/software/parmetis-4.0.3/build/Darwin-x86_64
SCOTCH_DIR    = $(HOME)/software/release/scotch_6.0.0



# Includes
PEXSI_INCLUDE    = -I${PEXSI_DIR}/include 
DSUPERLU_INCLUDE = -I${DSUPERLU_DIR}/SRC
INCLUDES         = ${PEXSI_INCLUDE} ${DSUPERLU_INCLUDE} 

# Libraries
GFORTRAN_LIB     = -lgfortran
LAPACK_LIB       = -llapack
BLAS_LIB         = -lblas
DSUPERLU_LIB     = ${DSUPERLU_DIR}/lib/libsuperlu_dist_3.3.a
PEXSI_LIB        = ${PEXSI_DIR}/src/libpexsi_${SUFFIX}.a

# Graph partitioning libraries
METIS_LIB        = -L${METIS_DIR}/libmetis -lmetis
PARMETIS_LIB     = -L${PARMETIS_DIR}/libparmetis -lparmetis 
SCOTCH_LIB       = -L${SCOTCH_DIR}/lib -lscotchmetis -lscotch -lscotcherr
PTSCOTCH_LIB     = -L${SCOTCH_DIR}/lib -lptscotchparmetis -lptscotch -lptscotcherr -lscotch


# Different compiling and linking options.
ifeq (${COMPILE_MODE}, release)
  COMPILE_DEF    = -DRELEASE
  COMPILE_FLAG   = -O3 -w 
endif
ifeq (${COMPILE_MODE}, debug)
  COMPILE_DEF    = -DDEBUG=1
  COMPILE_FLAG   = -O2 -w -g
endif

ifeq(${PAR_ND_LIBRARY}, ptscotch)
  PAR_ND_LIB = ${PTSCOTCH_LIB}
else
  PAR_ND_LIB = ${PARMETIS_LIB}
endif 

ifeq(${SEQ_ND_LIBRARY}, scotch)
  SEQ_ND_LIB = ${SCOTCH_LIB}
else
  SEQ_ND_LIB = ${METIS_LIB}
endif 

LIBS  = ${PEXSI_LIB} ${DSUPERLU_LIB} ${PAR_ND_LIB} ${SEQ_ND_LIB} ${LAPACK_LIB} ${BLAS_LIB} ${GFORTRAN_LIB}


COMPILE_DEFS += -DAdd_

#FLOADOPTS    = ${LIBS} -L/usr/local/lib  -lstdc++ 

CFLAGS       = ${COMPILE_FLAG} ${PROFILE_FLAG} ${INCLUDES}
FFLAGS       = ${COMPILE_FLAG} ${PROFILE_FLAG} ${INCLUDES}
CXXFLAGS     = ${COMPILE_FLAG} ${PROFILE_FLAG} ${INCLUDES} 
CCDEFS       = ${COMPILE_DEF} 
CPPDEFS      = ${COMPILE_DEF} 
LOADOPTS     = ${PROFILE_FLAG} ${LIBS}

