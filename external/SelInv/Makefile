include Makefile.in
include ${MAKEINCFILE}

default:
	(cd Lib ; make ) 
	(cd UTILITIES; make ) 

lib:
	( cd Lib ; make )

clean:
	(cd Lib ; make clean )
	(cd UTILITIES; make clean )
	(rm libselinv.a)
