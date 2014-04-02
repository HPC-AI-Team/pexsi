# NOTE: This Makefile does NOT support auto-dependency for the .h files.
# If the header files are changed, do "make clean" first.

include ../make.inc

SRCS_CPP  = run_inertia.cpp run_ppexsi.cpp run_pselinv.cpp 
SRCS_C    = driver_pselinv_real.c driver_pselinv_complex.c driver_ksdft.c 
SRCS_F90  = f_driver_pselinv_real.f90 f_interface.f90

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

all: driver_pselinv_complex driver_ksdft 

f_driver_pselinv_real: f_interface.o f_driver_pselinv_real.o ${PEXSI_LIB} 
	($(FC) -o $@_${SUFFIX} f_driver_pselinv_real.o  $(FLOADOPTS) )

f_driver_pselinv_real.o: f_driver_pselinv_real.f90
	  ${FC} -c ${FFLAGS} $<

f_interface.o: f_interface.f90
	  ${FC} -c ${FFLAGS} $<

driver_ksdft: driver_ksdft.o ${PEXSI_LIB} 
	($(LOADER) -o $@_${SUFFIX} driver_ksdft.o  $(LOADOPTS) )

driver_pselinv_real: driver_pselinv_real.o ${PEXSI_LIB} 
	($(LOADER) -o $@_${SUFFIX} driver_pselinv_real.o  $(LOADOPTS) )

driver_pselinv_complex: driver_pselinv_complex.o ${PEXSI_LIB} 
	($(LOADER) -o $@_${SUFFIX} driver_pselinv_complex.o  $(LOADOPTS) )

run_inertia: run_inertia.o ${PEXSI_LIB} 
	($(LOADER) -o $@_${SUFFIX} run_inertia.o  $(LOADOPTS) )

run_pselinv: run_pselinv.o ${PEXSI_LIB} 
	($(LOADER) -o $@ run_pselinv.o  $(LOADOPTS) )

run_inertia.o: run_inertia.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} $< 

run_pselinv.o: run_pselinv.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} $< 

f_ppexsi: f_ppexsi.o ${PEXSI_LIB} 
	  ($(FC) -o $@_${SUFFIX} f_ppexsi.o  $(FLOADOPTS) )

f_ppexsi.o: f_ppexsi.f90
	  ${FC} -c ${FFLAGS} $<

run_ppexsi: run_ppexsi.o ${PEXSI_LIB} 
	($(LOADER) -o $@ run_ppexsi.o  $(LOADOPTS) )

run_ppexsi.o: run_ppexsi.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} $< 

-include ${DEPS}

${PEXSI_LIB}:
	(cd ${PEXSI_DIR}/src; make all)

cleanall:
	rm -f ${EXES} ${OBJS} ${DEPS} *.d.* *.o *.mod
