# NOTE: This Makefile does NOT support auto-dependency for the .h files.
# If the header files are changed, do "make clean" first.

include ../make.inc.${HOSTID}

SRCS = ex1.cpp ex2.cpp ex3.cpp ex4.cpp ex5.cpp ex6.cpp ex7.cpp ex8.cpp\
			 ex9.cpp ex10.cpp ex11.cpp ex12.cpp ex13.cpp ex14.cpp ex15.cpp\
			 ex16.cpp

OBJS = ${SRCS:.cpp=.o}
DEPS = ${SRCS:.cpp=.d}
EXES = ${SRCS:.cpp=}	

ex3: ex3.o ${PEXSI_LIB} 
	($(LOADER) -o $@ ex3.o $(LOADOPTS) )

ex2: ex2.o ${PEXSI_LIB} 
	($(LOADER) -o $@ ex2.o $(LOADOPTS) )

ex1: ex1.o ${PEXSI_LIB} 
	($(LOADER) -o $@ ex1.o $(LOADOPTS) )

-include ${DEPS}

${PEXSI_LIB}:
	(cd ${PEXSI_DIR}/src; make all)

cleanall:
	rm -f ${EXES} ${OBJS} ${DEPS} *.d.o
