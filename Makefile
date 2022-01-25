include make.inc

lib: pexsi_lib

examples: pexsi_examples fortran_examples

all: lib examples

pexsi_lib:
	(cd src && ${MAKE})
	(cp src/f_ppexsi_interface.mod ${PEXSI_DIR}/include)

pexsi_examples: pexsi_lib
	cd examples && ${MAKE}

fortran_examples: pexsi_lib
	cd fortran && ${MAKE}

install: pexsi_lib
	(if [ ! -d "${PEXSI_BUILD_DIR}" ]; then mkdir ${PEXSI_BUILD_DIR}; fi)
	(if [ ! -d "${PEXSI_BUILD_DIR}/include" ]; then mkdir ${PEXSI_BUILD_DIR}/include; fi)
	(if [ ! -d "${PEXSI_BUILD_DIR}/lib" ]; then mkdir ${PEXSI_BUILD_DIR}/lib; fi)
	(cp src/libpexsi_${SUFFIX}.a ${PEXSI_BUILD_DIR}/lib)
	(cp include/c_pexsi_interface.h include/ppexsi.hpp ${PEXSI_BUILD_DIR}/include)
	(cp src/f_ppexsi_interface.mod ${PEXSI_BUILD_DIR}/include)

clean:
	cd src && ${MAKE} clean

cleanall:
	cd src && ${MAKE} cleanall
	cd examples && ${MAKE} cleanall
	cd fortran && ${MAKE} cleanall
