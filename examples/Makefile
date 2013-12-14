# NOTE: This Makefile does NOT support auto-dependency for the .h files.
# If the header files are changed, do "make clean" first.

include ../make.inc

SRCS_CPP  = run_pexsi.cpp run_ppexsi.cpp run_pselinv.cpp 
SRCS_C    = driver_pselinv.c driver_ksdft.c

OBJS = ${SRCS_CPP:.cpp=.o} ${SRCS_C:.c=.o}
DEPS = ${SRCS_CPP:.cpp=.d} ${SRCS_C:.c=.d}
EXES = ${SRCS_CPP:.cpp=} ${SRCS_C:.c=}

# Compilation replacement rules

%.o: %.c
	${CC} -c ${CFLAGS} ${CCDEFS} $< 
%.o: %.cpp
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} $< 
%.o: %.f
	${FC} -c ${FFLAGS} $<
%.o: %.F
	${FC} -c ${FFLAGS} $<

all: run_pselinv run_pexsi run_ppexsi f_ppexsi driver_pselinv

driver_ksdft: driver_ksdft.o ${PEXSI_LIB} 
	($(LOADER) -o $@ driver_ksdft.o  $(LOADOPTS_PTSCOTCH) )

driver_pselinv: driver_pselinv.o ${PEXSI_LIB} 
	($(LOADER) -o $@ driver_pselinv.o  $(LOADOPTS_PTSCOTCH) )

run_pselinv: run_pselinv.o ${PEXSI_LIB} run_pselinv_ptscotch run_pselinv_parmetis
	($(LOADER) -o $@ run_pselinv.o  $(LOADOPTS_PTSCOTCH) )

run_pselinv_parmetis: run_pselinv.o ${PEXSI_LIB} 
	($(LOADER) -o $@_${SUFFIX} run_pselinv.o  $(LOADOPTS_PARMETIS) )

run_pselinv_ptscotch: run_pselinv.o   ${PEXSI_LIB} 
	($(LOADER) -o $@_${SUFFIX} run_pselinv.o  $(LOADOPTS_PTSCOTCH) )


run_pselinv.o: run_pselinv.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} $< 

f_ppexsi: f_ppexsi.o ${PEXSI_LIB} 
	  ($(FC) -o $@_${SUFFIX} f_ppexsi.o  $(FLOADOPTS) )

f_ppexsi.o: f_ppexsi.f90
	  ${FC} -c ${FFLAGS} $<

run_pexsi: run_pexsi.o ${PEXSI_LIB} 
	($(LOADER) -o $@ run_pexsi.o  $(LOADOPTS) )

run_pexsi.o: run_pexsi.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} $< 

run_ppexsi: run_ppexsi.o ${PEXSI_LIB} 
	($(LOADER) -o $@ run_ppexsi.o  $(LOADOPTS) )

run_ppexsi.o: run_ppexsi.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} $< 

-include ${DEPS}

${PEXSI_LIB}:
	(cd ${PEXSI_DIR}/src; make all)

cleanall:
	rm -f ${EXES} ${OBJS} ${DEPS} *.d.* *.o
