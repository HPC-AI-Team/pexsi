# NOTE: This Makefile does NOT support auto-dependency for the .h files.
# If the header files are changed, do "make clean" first.

include ../make.inc

#SRCS = ex1.cpp ex2.cpp ex3.cpp ex4.cpp ex5.cpp ex6.cpp ex7.cpp \
			 ex8.cpp ex9.cpp ex10.cpp ex11.cpp ex12.cpp ex13.cpp ex14.cpp ex15.cpp
SRCS   = ex23.cpp
SRCS_C = ex18.c
SRCS_F = ex22.f90 

OBJS = ${SRCS:.cpp=.o} ${SRCS_C:.c=.o} ${SRCS_F:.f90=.o}
DEPS = ${SRCS:.cpp=.d} ${SRCS_C:.c=.d} ${SRCS_F:.f90=.d}
EXES = ${SRCS:.cpp=}	 ${SRCS_C:.c=}   ${SRCS_F:.f90=}

# Compilation replacement rules

%.o: %.c
	${CC} -c ${CFLAGS} ${CCDEFS} $< 
%.o: %.cpp
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} $< 
%.o: %.f90
	${FC} -c ${FFLAGS} $<
%.o: %.F90
	${FC} -c ${FFLAGS} $<

ex25: ex25.o ${PEXSI_LIB} 
	($(LOADER) -o $@ ex25.o  $(LOADOPTS) )

ex24: ex24.o ${PEXSI_LIB} 
	($(LOADER) -o $@ ex24.o  $(LOADOPTS) )

ex23: ex23.o ${PEXSI_LIB} 
	($(LOADER) -o $@ ex23.o  $(LOADOPTS) )

ex22: ex22.o ${PEXSI_LIB} 
	($(FC) -o $@ ex22.o  $(FLOADOPTS) )

ex21: ex21.o ${PEXSI_LIB} 
	($(LOADER) -o $@ ex21.o  $(LOADOPTS) )

ex20: ex20.o ${PEXSI_LIB} 
	($(FC) -o $@ ex20.o  $(FLOADOPTS) )

ex19: ex19.o ${PEXSI_LIB} 
	($(FC) -o $@ ex19.o  $(FLOADOPTS) )

ex18: ex18.o ${PEXSI_LIB} 
	($(LOADER) -o $@ ex18.o  $(LOADOPTS) )

ex17: ex17.o ${PEXSI_LIB} 
	($(LOADER) -o $@ ex17.o  $(LOADOPTS) )

ex17.o: ex17.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} $< 

ex16: ex16.o ${PEXSI_LIB} 
	($(LOADER) -o $@ ex16.o  $(LOADOPTS) )

ex16.o: ex16.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} $< 


ex15: ex15.o ${PEXSI_LIB} 
	($(LOADER) -o $@ ex15.o  $(LOADOPTS) $(IPM) )

ex15.o: ex15.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} $< 

ex14: ex14.o ${PEXSI_LIB} ${PLUSELINV_LIB}
	($(LOADER) -o $@ ex14.o  $(LOADOPTS) )

ex14.o: ex14.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} $< 


ex13: ex13.o ${PEXSI_LIB} ${PLUSELINV_LIB}
	($(LOADER) -o $@ ex13.o  $(LOADOPTS) )

ex13.o: ex13.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} $< 


ex12: ex12.o ${PEXSI_LIB} ${PLUSELINV_LIB}
	($(LOADER) -o $@ ex12.o  $(LOADOPTS) )

ex12.o: ex12.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} -D_HAS_DOUBLECOMPLEX_ $< 

ex11: ex11.o ${PEXSI_LIB} ${PLUSELINV_LIB}
	($(LOADER) -o $@ ex11.o  $(LOADOPTS) )

ex11.o: ex11.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} -D_HAS_DOUBLECOMPLEX_ $< 

ex10: ex10.o ${PEXSI_LIB} ${PLUSELINV_LIB}
	($(LOADER) -o $@ ex10.o  $(LOADOPTS) )

ex10.o: ex10.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS}  $< 

ex9: ex9.o ${PEXSI_LIB} ${PLUSELINV_LIB}
	($(LOADER) -o $@ ex9.o  $(LOADOPTS) )

ex9.o: ex9.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} -D_HAS_DOUBLECOMPLEX_ $< 

ex8: ex8.o ${PEXSI_LIB} ${PLUSELINV_LIB}
	($(LOADER) -o $@ ex8.o  $(LOADOPTS) )

ex8.o: ex8.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} -D_HAS_DOUBLECOMPLEX_ $< 

ex7: ex7.o ${PEXSI_LIB} ${PLUSELINV_LIB}
	($(LOADER) -o $@ ex7.o $(LOADOPTS) )

ex7.o: ex7.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} -D_HAS_DOUBLECOMPLEX_ $< 

ex6: ex6.o ${PEXSI_LIB} ${PLUSELINV_LIB}
	($(LOADER) -o $@ ex6.o $(LOADOPTS) )

ex6.o: ex6.cpp 
	${CXX} -c ${CXXFLAGS} ${CPPDEFS} -D_HAS_DOUBLECOMPLEX_ $< 

ex5: ex5.o ${PEXSI_LIB} ${PLUSELINV_LIB}
	($(LOADER) -o $@ ex5.o $(LOADOPTS) )

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
