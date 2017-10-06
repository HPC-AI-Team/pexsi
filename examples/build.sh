#! /bin/bash
cd ../src; make cleanall; make USE_SYMPACK=$2 -j; cd ../examples; rm $1; rm $1.o; make  USE_SYMPACK=$2 $1 
