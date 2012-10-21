include Makefile.in
include ${MAKEINCFILE}

SRCS_C   = zcreate_matrix.c dcreate_matrix.c sp_ienv.c dreadhb.c dreadtriple.c
SRCS_CPP = runselinv.cpp pluselinv.cpp read_and_dist_rua.cpp
OBJS     = ${SRCS_C:.c=.o} ${SRCS_CPP:.cpp=.o} 
DEPS     = ${SRCS_C:.c=.d} ${SRCS_CPP:.cpp=.d} 


all: libpluselinv.a

-include ${DEPS}

libpluselinv.a: ${OBJS}
	$(AR) $(ARFLAGS) $@ ${OBJS}
	$(RANLIB) $@

runselinv: runselinv.o libpluselinv.a
	$(LOADER) $(LOADOPTS) runselinv.o -o $@ 


#DEXM	= pddrive.o dcreate_matrix.o sp_ienv.o dreadhb.o dreadtriple.o 

#SELINV	= dcreate_matrix.o sp_ienv.o dreadhb.o dreadtriple.o pluselinv.o read_and_dist_rua.o

#lusolve: lusolve.o $(SELINV) $(DSUPERLULIB)
#	$(LOADER) $(LOADOPTS) lusolve.o $(SELINV) $(LIBS) $(IPM)  -lm  -o $@ 


#c5by5: c5by5.o $(DSUPERLULIB)
#	$(LOADER) $(LOADOPTS) c5by5.o $(LIBS) -lm -o $@

#rua2bin: rua2bin.o dreadhb_serial.o
#	$(LOADER) $(LOADOPTS) rua2bin.o dreadhb_serial.o $(LIBS) -lm -o $@

#pddrive: $(DEXM) $(DSUPERLULIB)
#	$(LOADER) $(LOADOPTS) $(DEXM) $(LIBS) -lm -o $@


cleanall:	
	${RM} -f ${OBJS} ${DEPS} libpluselinv.a *.d.* runselinv.o
