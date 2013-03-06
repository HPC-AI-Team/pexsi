# NOTE: This Makefile does NOT support auto-dependency for the .h files.
# If the header files are changed, do "make clean" first.

include ../make.inc

SRCS  = run_pexsi.cpp run_ppexsi.cpp

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

all: run_pexsi run_ppexsi

run_pexsi: run_pexsi.o ${PEXSI_LIB} ${PLUSELINV_LIB}
	($(LOADER) -o $@ run_pexsi.o  $(LOADOPTS) )

run_pexsi.o: run_pexsi.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} $< 

run_ppexsi: run_ppexsi.o ${PEXSI_LIB} ${PLUSELINV_LIB}
	($(LOADER) -o $@ run_ppexsi.o  $(LOADOPTS) -mkl $(IPM))

run_ppexsi.o: run_ppexsi.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} $< 

-include ${DEPS}

${PEXSI_LIB}:
	(cd ${PEXSI_DIR}/src; make all)

cleanall:
	rm -f ${EXES} ${OBJS} ${DEPS} *.d.* *.o
