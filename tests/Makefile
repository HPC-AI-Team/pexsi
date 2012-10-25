# NOTE: This Makefile does NOT support auto-dependency for the .h files.
# If the header files are changed, do "make clean" first.

include ../make.inc.${HOSTID}

SRCS = ex1.cpp ex2.cpp ex3.cpp 

OBJS = ${SRCS:.cpp=.o}
DEPS = ${SRCS:.cpp=.d}
EXES = ${SRCS:.cpp=}	

# Compilation replacement rules

%.o: %.c
	${CC} -c ${CFLAGS} ${CCDEFS} $< 
%.o: %.cpp
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} $< 
%.o: %.f
	${FC} -c ${FFLAGS} $<
%.o: %.F
	${FC} -c ${FFLAGS} $<



ex4: ex4.o ${PEXSI_LIB} ${PLUSELINV_LIB}
	($(LOADER) -o $@ ex4.o $(LOADOPTS) )

ex4.o: ex4.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} -D_HAS_DOUBLECOMPLEX_ $< 

ex3: ex3.o ${PEXSI_LIB} ${PLUSELINV_LIB}
	($(LOADER) -o $@ ex3.o $(LOADOPTS) )

ex2: ex2.o ${PEXSI_LIB} 
	($(LOADER) -o $@ ex2.o $(LOADOPTS) )

ex1: ex1.o ${PEXSI_LIB} 
	($(LOADER) -o $@ ex1.o $(LOADOPTS) )

-include ${DEPS}

${PEXSI_LIB}:
	(cd ${PEXSI_DIR}/src; make all)

cleanall:
	rm -f ${EXES} ${OBJS} ${DEPS} *.d.o *.o
