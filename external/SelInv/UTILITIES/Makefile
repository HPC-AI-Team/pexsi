include ../Makefile.in
include ../${MAKEINCFILE}

all: dconvmat.x cconvmat.x

dconvmat.x: dconvmat.o mmio.o
	$(FC) $(FFLAGS) $^ -o $@

dconvmat.o: dconvmat.f90
	$(FC) $(FFLAGS) $< -c

cconvmat.x: cconvmat.o mmio.o
	$(FC) $(FFLAGS) $^ -o $@

cconvmat.o: cconvmat.f90
	$(FC) $(FFLAGS) $< -c

mmio.o: mmio.f
	$(FC) $(FFLAGS) $< -c

clean: 
	rm -rf *.o *.x 
