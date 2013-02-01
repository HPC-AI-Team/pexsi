/// @file interface.cpp
/// @brief Interface subroutines of PEXSI that can be called by
/// C/FORTRAN
/// @author Lin Lin
/// @date 2013-01-31
#include "c_pexsi_interface.h"
#include "ppexsi.hpp"

using namespace PEXSI;

// *********************************************************************
// C interface
// *********************************************************************
void DummyInterface( int* a ){
	std::cout << "Dummy inteface is working and is outputing an integer " 
		<< *a << std::endl;
	return;
}	

// *********************************************************************
// FORTRAN interface
// *********************************************************************
extern "C" void FORTRAN(f_dummy_interface)( int* a ){
	DummyInterface( a );
	return;
}
