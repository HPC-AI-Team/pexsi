#include  "mpi_interf.hpp"

namespace PEXSI{

namespace mpi{

// *********************************************************************
// Allgather
// *********************************************************************

void
Allgatherv ( 
		std::vector<Int>& localVec, 
		std::vector<Int>& allVec,
		MPI_Comm          comm )
{
#ifndef _RELEASE_
	PushCallStack("Allgatherv");
#endif
	Int mpirank, mpisize;
	MPI_Comm_rank( comm, &mpirank );
	MPI_Comm_size( comm, &mpisize );

	Int localSize = localVec.size();
	std::vector<Int>  localSizeVec( mpisize );
	std::vector<Int>  localSizeDispls( mpisize );
	MPI_Allgather( &localSize, 1, MPI_INT, &localSizeVec[0], 1, MPI_INT, comm );
	localSizeDispls[0] = 0;
	for( Int ip = 1; ip < mpisize; ip++ ){
    localSizeDispls[ip] = localSizeDispls[ip-1] + localSizeVec[ip-1];
	}
	Int totalSize = localSizeDispls[mpisize-1] + localSizeVec[mpisize-1];

	allVec.clear();
	allVec.resize( totalSize );

	MPI_Allgatherv( &localVec[0], localSize, MPI_INT, &allVec[0], 
		 &localSizeVec[0], &localSizeDispls[0], MPI_INT, comm	);

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
}		// -----  end of function Allgatherv  ----- 

} // namespace mpi

} // namespace PEXSI
