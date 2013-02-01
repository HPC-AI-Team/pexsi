/// @file ex18.c
/// @brief Test the interface between C and PEXSI.
/// @author Lin Lin
/// @date 2013-01-31
#include "c_pexsi_interface.h"

int
main ( int argc, char *argv[] )
{
	MPI_Init(&argc, &argv);
	
	int a = 2;
	DummyInterface( MPI_COMM_WORLD, a );
	MPI_Finalize();
  	
	return 0;
}				// ----------  end of function main  ---------- 
