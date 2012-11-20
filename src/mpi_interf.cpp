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
	PushCallStack("mpi::Allgatherv");
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


// *********************************************************************
// Send / Recv
// *********************************************************************
void 
Send( const std::stringstream& sstm, Int dest, Int tagSize, Int tagContent, 
		MPI_Comm comm ){
#ifndef _RELEASE_
	PushCallStack("mpi::Send");
#endif
	const std::string sstr = sstm.str();
	Int sizeStm = sstr.length();
	MPI_Send( &sizeStm, 1, MPI_INT,  dest, tagSize, comm );
	MPI_Send( (void*)sstr.c_str(), sizeStm, MPI_BYTE, dest, tagContent, comm );
#ifndef _RELEASE_
	PopCallStack();
#endif
	return; 
} // -----  end of function Send ----- 


void
Recv ( std::stringstream& sstm, Int src, Int tagSize, Int tagContent, 
		MPI_Comm comm, MPI_Status& statSize, MPI_Status& statContent )
{
#ifndef _RELEASE_
	PushCallStack("mpi::Recv");
#endif
	std::vector<char> str;
	Int sizeStm;
	MPI_Recv( &sizeStm, 1, MPI_INT, src, tagSize, comm, &statSize );
	str.resize( sizeStm );
	MPI_Recv( (void*) &str[0], sizeStm, MPI_BYTE, src, tagContent, comm, &statContent );
	sstm.write( &str[0], sizeStm );
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
}		// -----  end of function Recv  ----- 

void
Recv ( std::stringstream& sstm, Int src, Int tagSize, Int tagContent, 
		MPI_Comm comm )
{
#ifndef _RELEASE_
	PushCallStack("mpi::Recv (MPI_STATUS_IGNORE)");
#endif
	std::vector<char> str;
	Int sizeStm;
	MPI_Recv( &sizeStm, 1, MPI_INT, src, tagSize, comm, MPI_STATUS_IGNORE );
	str.resize( sizeStm );
	MPI_Recv( (void*) &str[0], sizeStm, MPI_BYTE, src, tagContent, comm, MPI_STATUS_IGNORE );
	sstm.write( &str[0], sizeStm );
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
}		// -----  end of function Recv  ----- 


void
Isend ( const std::stringstream& sstm, Int dest, Int tagSize, Int tagContent, 
		MPI_Comm comm, MPI_Request& reqSize, MPI_Request& reqContent )
{
#ifndef _RELEASE_
	PushCallStack("Isend");
#endif
	const std::string sstr = sstm.str();
	Int sizeStm = sstr.length();
	MPI_Isend( &sizeStm, 1, MPI_INT,  dest, tagSize, comm, &reqSize );
	MPI_Isend( (void*)sstr.c_str(), sizeStm, MPI_BYTE, dest, tagContent, comm, &reqContent );

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
}		// -----  end of function Isend  ----- 


void
Irecv ( std::stringstream& sstm, Int src, Int tagSize, Int tagContent, 
		MPI_Comm comm, MPI_Request& reqSize, MPI_Request& reqContent )
{
#ifndef _RELEASE_
	PushCallStack("Irecv");
#endif
	std::vector<char> str;
	Int sizeStm;
	MPI_Irecv( &sizeStm, 1, MPI_INT, src, tagSize, comm, &reqSize );
	str.resize( sizeStm );
	MPI_Irecv( (void*) &str[0], sizeStm, MPI_BYTE, src, tagContent, comm, &reqContent );
	sstm.write( &str[0], sizeStm );
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
}		// -----  end of function Irecv  ----- 


// *********************************************************************
// Waitall
// *********************************************************************

void
Waitall ( std::vector<MPI_Request>& reqs, std::vector<MPI_Status>& stats )
{
#ifndef _RELEASE_
	PushCallStack("mpi::Waitall");
#endif
  if( reqs.size() != stats.size() ){
    throw std::runtime_error( "MPI_Request does not have the same as as MPI_Status." );
	}
	for( Int i = 0; i < reqs.size(); i++ ){
		MPI_Wait( &reqs[i], &stats[i] );
	}
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
}		// -----  end of function Waitall  ----- 

void
Waitall ( std::vector<MPI_Request>& reqs )
{
#ifndef _RELEASE_
	PushCallStack("mpi::Waitall (MPI_STATUS_IGNORE)");
#endif
	for( Int i = 0; i < reqs.size(); i++ ){
		MPI_Wait( &reqs[i], MPI_STATUS_IGNORE );
	}
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
}		// -----  end of function Waitall  ----- 


} // namespace mpi

} // namespace PEXSI
