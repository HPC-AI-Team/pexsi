include ../make.inc

SRCS_C   = pdsymbfact.c pzsymbfact.c get_perm_c_parmetis.c
SRCS_CPP = interface.cpp ppexsi.cpp pole.cpp TreeBcast.cpp\
					 superlu_dist_internal_complex.cpp superlu_dist_internal_real.cpp \
					 mpi_interf.cpp lapack.cpp blas.cpp utility.cpp global.cpp timer.cpp \
					 # getPole.cpp

SRCS_CPP_SLAVE = slave_hello_world.cpp \
				 slave_JB_Loop_z.cpp slave_JB_Loop_Index_Compress_z.cpp slave_JB_Loop_Index_Compress_dma_z.cpp \
				 slave_JB_Loop_d.cpp slave_JB_Loop_Index_Compress_d.cpp slave_JB_Loop_Index_Compress_dma_d.cpp
				 
					

SRCS_F90 = f_interface.f90


SRCS_C_OBJ = ${SRCS_C:.c=.o} 
SRCS_CPP_OBJ = ${SRCS_CPP:.cpp=.o}
SRCS_CPP_SLAVE_OBJ = ${SRCS_CPP_SLAVE:.cpp=.o}
SRCS_F90_OBJ = ${SRCS_F90:.f90=.o}

DEPS     = ${SRCS_CPP:.cpp=.d} ${SRCS_C:.c=.d} 
OBJS     = ${SRCS_C_OBJ} ${SRCS_CPP_OBJ} ${SRCS_CPP_SLAVE_OBJ} ${SRCS_F90_OBJ} getPole.o

all: ${PEXSI_LIB} 

# Compilation replacement rules

${SRCS_C_OBJ} : %.o: %.c
	${CC} -mhost -msimd -c ${CFLAGS} ${CCDEFS} $< 
${SRCS_CPP_OBJ} : %.o: %.cpp
	${CXX} -mhost -msimd -c $< ${CXXFLAGS} ${CPPDEFS} 
${SRCS_CPP_SLAVE_OBJ} : %.o: %.cpp
	${CXX} -mslave -msimd -c $< ${CXXFLAGS} ${CPPDEFS} 
%.o: %.f
	${FC} -c ${FFLAGS} $<
%.o: %.F
	${FC} -c ${FFLAGS} $<
%.o: %.f90
	${FC} -c ${FFLAGS} $<
%.o: %.F90
	${FC} -c ${FFLAGS} $<


-include ${DEPS}

#libpexsi.a: ${OBJS}
${PEXSI_LIB}: ${OBJS}
	$(AR) $(ARFLAGS) $(PEXSI_LIB) ${OBJS} 
	$(RANLIB) $(PEXSI_LIB)

#	$(RANLIB) $(shell basename $@ .a)_${SUFFIX}.a
#ifeq (${MODE},release)
#ifneq ($(USE_TAU),1)
#	${CP} $(shell basename $@ .a)_${SUFFIX}.a ${PEXSI_LIB}
#endif
#endif

clean:
	${RM} -f ${OBJS} ${DEPS} *.d.* *.mod

cleanall:
	${RM} -f ${OBJS} ${DEPS} ${PEXSI_DIR}/lib/libpexsi.a *.d.* *.mod
