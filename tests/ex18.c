/// @file ex18.c
/// @brief Test the interface between C and PEXSI.
/// @author Lin Lin
/// @date 2013-01-31
#include "c_pexsi_interface.h"

int
main ( int argc, char *argv[] )
{
	int a = 2;

	DummyInterface( &a );
  	
	return 0;
}				// ----------  end of function main  ---------- 
