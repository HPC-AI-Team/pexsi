#! /bin/bash
cd ../src; make cleanall; make -j; cd ../examples; make cleanall; make $1 
