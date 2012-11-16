include ../make.inc.${HOSTID}

SRCS_C   = pzsymbfact.c sp_ienv.c
SRCS_CPP = pexsi.cpp pole.cpp superlu_dist_interf.cpp selinv_interf.cpp \
			     lapack.cpp blas.cpp utility.cpp global.cpp 

OBJS     = ${SRCS_C:.c=.o} ${SRCS_CPP:.cpp=.o} 
DEPS     = ${SRCS_C:.c=.d} ${SRCS_CPP:.cpp=.d} 

# Compilation replacement rules

%.o: %.c
	${CC} -c ${CFLAGS} ${CCDEFS} $< 
%.o: %.cpp
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} $< 
%.o: %.f
	${FC} -c ${FFLAGS} $<
%.o: %.F
	${FC} -c ${FFLAGS} $<


all: libpexsi.a

-include ${DEPS}

libpexsi.a: ${OBJS}
	$(AR) $(ARFLAGS) $@ ${OBJS}
	$(RANLIB) $@

cleanall:
	${RM} -f ${OBJS} ${DEPS} libpexsi.a *.d.*
