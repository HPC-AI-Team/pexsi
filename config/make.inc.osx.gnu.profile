#/usr/bin/bash
PROFILE        =   0
PROF_MPI       =   0
USE_TAU        =   0
USE_AUTO_TAU   =   0


# Different compiling and linking options.
#MODE           = debug
MODE	         = release
PLATFORM       = osx_v0.5.5
SUFFIX         = $(MODE)_${PLATFORM}

ifdef USE_TAU
  ifeq ($(USE_TAU),1)
    SUFFIX=$(MODE)_tau_${PLATFORM}
  endif

  ifeq ($(PROFILE),1)
    SUFFIX=$(MODE)_tau_${PLATFORM}
  endif

  ifeq ($(PROF_MPI),1)
    SUFFIX=$(MODE)_tau_${PLATFORM}
  endif
endif


PEXSI_DIR     = $(HOME)/Projects/pexsi
DSUPERLU_DIR  = $(HOME)/Documents/Software/SuperLU_DIST_3.3
METIS_DIR     = $(HOME)/Software/metis-5.1.0/build_release
PARMETIS_DIR  = $(HOME)/Software/parmetis-4.0.2/build/Darwin-x86_64
SCOTCH_DIR    = $(HOME)/Software/scotch_6.0.0

# includes
PEXSI_INCLUDE    = -I${PEXSI_DIR}/include 
DSUPERLU_INCLUDE = -I${DSUPERLU_DIR}/SRC
INCLUDES         = ${PEXSI_INCLUDE} ${DSUPERLU_INCLUDE} 

ifeq ($(USE_TAU),1)
  INCLUDES += -I${TAUROOTDIR}/include 
endif

# Libraries
GFORTRAN_LIB     = /usr/local/lib/libgfortran.dylib
LAPACK_LIB       = -llapack
BLAS_LIB         = -lblas
METIS_LIB        = -L${METIS_DIR}/lib -lmetis
PARMETIS_LIB     = -L${PARMETIS_DIR}/libparmetis -lparmetis 
SCOTCH_LIB       = -L${SCOTCH_DIR}/lib -lptscotchparmetis -lptscotch -lptscotcherr -lscotch
DSUPERLU_LIB     = ${DSUPERLU_DIR}/build_print/libsuperlu_dist_3.3.a
PEXSI_LIB        = ${PEXSI_DIR}/src/libpexsi_${SUFFIX}.a


ifdef USE_TAU
  ifeq ($(USE_TAU),1)
TAU_LIB          = -L/opt/cray/mpt/5.5.2/gni/mpich2-pgi/119/lib -L/usr/common/acts/TAU/tau-2.22/craycnl/lib -lTauMpi-papi-mpi-pdt-pgi -lpthread -lrt -lmpichcxx -lmpich -lrt -L/usr/common/acts/TAU/tau-2.22/craycnl/lib -ltau-papi-mpi-pdt-pgi -R/opt/cray/papi/4.3.0.1/perf_events/no-cuda/lib -L/opt/cray/papi/4.3.0.1/perf_events/no-cuda/lib -lpapi -L/usr/common/acts/TAU/tau-2.22/craycnl/binutils-2.20/lib -L/usr/common/acts/TAU/tau-2.22/craycnl/binutils-2.20/lib64 -lbfd -liberty -lz -Wl,--export-dynamic -lrt -L/opt/pgi/12.5.0/linux86-64/12.5/bin/../lib -lstd -lC -lstdc++ -L/usr/common/acts/TAU/tau-2.22/craycnl/lib/static-papi-mpi-pdt-pgi
  endif
endif


LIBS_PARMETIS    = ${PEXSI_LIB} ${DSUPERLU_LIB} ${PARMETIS_LIB} ${METIS_LIB} ${TAU_LIB} ${IPM}  ${LAPACK_LIB} ${BLAS_LIB} ${GFORTRAN_LIB} 
LIBS_PTSCOTCH    = ${PEXSI_LIB} ${DSUPERLU_LIB} ${SCOTCH_LIB}  ${METIS_LIB} ${TAU_LIB} ${IPM} ${LAPACK_LIB} ${BLAS_LIB} ${GFORTRAN_LIB}
LIBS             = ${LIBS_PTSCOTCH}

CC           = mpicc 
CXX          = mpic++
FC           = mpif90
LOADER       = mpic++

ifdef USE_TAU
  ifeq ($(USE_TAU),1)
#  ifeq ($(USE_AUTO_TAU),1)
CC           = tau_cc.sh
CXX          = tau_cxx.sh
FC           = tau_f90.sh
LOADER       = tau_cxx.sh
#  endif
  endif
endif





AR           = ar 
ARFLAGS      = rvcu
# For System V based machine without ranlib, like Cray and SGI,
# use touch instead.
#RANLIB      = touch
RANLIB       = ranlib

CP           = cp
RM           = rm
RMFLAGS      = -f

# Different compiling and linking options.
ifeq ($(MODE), debug)
	COMMONDEFS   = -DDEBUG=1 -DAdd_  #-DUSE_REDUCE_L -DUSE_BCAST_UL -DPRINT_COMMUNICATOR_STAT #-DUSE_BCAST_UL #-DCOMPARE_LUPDATE #-DUSE_MPI_COLLECTIVES -DBLOCK_REDUCE
  CFLAGS       = -O0 -g -w ${INCLUDES} -DAdd_
  FFLAGS       = -O0 -g -w ${INCLUDES}
  CXXFLAGS     = -O0 -g -w ${INCLUDES} -DAdd_ #-DSANITY_CHECK -DSANITY_PRECISION=1e-5 #-DSELINV_TIMING -DSELINV_MEMORY -DNO_PARMETIS_FIX 
	CCDEFS       = ${COMMONDEFS}
	CPPDEFS      = ${COMMONDEFS}
  LOADOPTS     = ${LIBS}
  LOADOPTS_PARMETIS     = ${LIBS_PARMETIS}
  LOADOPTS_PTSCOTCH     = ${LIBS_PTSCOTCH}
  FLOADOPTS    = ${LIBS} -L/usr/local/lib  -lstdc++
endif

ifeq ($(MODE), release)
	COMMONDEFS   = -DRELEASE -g -DDEBUG=0 -DAdd_ #-DUSE_REDUCE_L -DUSE_BCAST_UL -DPRINT_COMMUNICATOR_STAT #-DUSE_BCAST_UL #-DCOMPARE_LUPDATE #-DUSE_MPI_COLLECTIVES 
  CFLAGS       = -O3  -w ${INCLUDES}  #-fopenmp
  FFLAGS       = -O3 -w ${INCLUDES}
  CXXFLAGS     = -O3 -w ${INCLUDES}  #-DNO_PARMETIS_FIX  -DSANITY_CHECK -DSANITY_PRECISION=1e-5 #-fopenmp
	CCDEFS       = ${COMMONDEFS}
	CPPDEFS      = ${COMMONDEFS}
  LOADOPTS     = ${LIBS}
  LOADOPTS_PARMETIS     = ${LIBS_PARMETIS}
  LOADOPTS_PTSCOTCH     = ${LIBS_PTSCOTCH}
  FLOADOPTS    = ${LIBS} -L/usr/local/lib  -lstdc++ 
endif



