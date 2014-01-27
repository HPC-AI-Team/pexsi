include ../make.inc

SRCS_C   = pdsymbfact.c pzsymbfact.c get_perm_c_parmetis.c
SRCS_CPP = interface.cpp ppexsi.cpp ppexsi_real.cpp pole.cpp \
					 pselinv.cpp superlu_dist_interf.cpp \
					 mpi_interf.cpp lapack.cpp blas.cpp utility.cpp global.cpp timer.cpp 

OBJS     = ${SRCS_CPP:.cpp=.o} ${SRCS_C:.c=.o} ${SRCS_F:.f=.o}
DEPS     = ${SRCS_CPP:.cpp=.d} ${SRCS_C:.c=.d} ${SRCS_F:.f=.d}

# Compilation replacement rules

%.o: %.c
	${CC} -c ${CFLAGS} ${CCDEFS} $< 
%.o: %.cpp
	${CXX} -c $< ${CXXFLAGS} ${CPPDEFS} 
%.o: %.f
	${FC} -c ${FFLAGS} $<
%.o: %.F
	${FC} -c ${FFLAGS} $<


all: libpexsi.a 

-include ${DEPS}

libpexsi.a: ${OBJS}
	$(AR) $(ARFLAGS) $(PEXSI_LIB) ${OBJS}
	$(RANLIB) $(shell basename $@ .a)_${SUFFIX}.a
#ifeq (${MODE},release)
#ifneq ($(USE_TAU),1)
#	${CP} $(shell basename $@ .a)_${SUFFIX}.a ${PEXSI_LIB}
#endif
#endif

clean:
	${RM} -f ${OBJS} ${DEPS} *.d.*

cleanall:
	${RM} -f ${OBJS} ${DEPS} libpexsi.a *.d.*
