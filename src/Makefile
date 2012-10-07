include ../make.inc.${HOSTID}

SRCS = getpolef.o selinv_interf.o pepsi.o 

OBJS = ${SRCS:.cpp=.o}
DEPS = ${SRCS:.cpp=.d}

all: libpexsi.a

-include ${DEPS}

libpexsi.a: ${OBJS}
	$(AR) $(ARFLAGS) $@ ${OBJS}
	$(RANLIB) $@

cleanall:
	${RM} -f ${OBJS} ${DEPS} libpexsi.a *.d.*
