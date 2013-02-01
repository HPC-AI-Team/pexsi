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
extern "C" void DummyInterface( MPI_Comm comm, int a ){
	int mpirank, mpisize;
	MPI_Comm_rank( comm, &mpirank );
	MPI_Comm_size( comm, &mpisize );
	if( mpirank == 0 ){
		std::cout << "Comm rank = " << mpirank << std::endl;
		std::cout << "Comm size = " << mpisize << std::endl;
		std::cout << "Dummy inteface is working and is outputing an integer " 
			<< a << std::endl;
	}
	return;
}	

// *********************************************************************
// FORTRAN interface
// *********************************************************************
/// @brief Internal subroutine to convert FORTRAN communicator to C
extern "C" MPI_Comm f2c_comm(int *Fcomm)
{
	return MPI_Comm_f2c((MPI_Fint)(*Fcomm));
}


extern "C" void FORTRAN(f_dummy_interface)( int* Fcomm, int* a ){
	DummyInterface( f2c_comm(Fcomm), *a );
	return;
}
