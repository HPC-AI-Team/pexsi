include ../make.inc.${HOSTID}

SRCS = pexsi.cpp getpolef.cpp selinv_interf.cpp \
			 lapack.cpp blas.cpp utility.cpp global.cpp
#SRCS =   pexsi.cpp 

OBJS = ${SRCS:.cpp=.o}
DEPS = ${SRCS:.cpp=.d}

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
