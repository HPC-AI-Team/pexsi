#include  "mpi_interf.hpp"

namespace PEXSI{

namespace mpi{

// *********************************************************************
// Gather
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
	const std::string& sstr = sstm.str();
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
	std::vector<char> sstr;
	Int sizeStm;
	MPI_Recv( &sizeStm, 1, MPI_INT, src, tagSize, comm, &statSize );
	sstr.resize( sizeStm );
	MPI_Recv( (void*) &sstr[0], sizeStm, MPI_BYTE, src, tagContent, comm, &statContent );
	sstm.write( &sstr[0], sizeStm );
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



// *********************************************************************
// Wait
// *********************************************************************


void
Wait	( MPI_Request& req  )
{
#ifndef _RELEASE_
	PushCallStack("mpi::Wait");
#endif
  MPI_Wait( &req, MPI_STATUS_IGNORE );
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method Wait  ----- 

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


// *********************************************************************
// Reduce
// *********************************************************************


void
Reduce ( Real* sendbuf, Real* recvbuf, Int count, MPI_Op op, Int root, MPI_Comm comm )
{
#ifndef _RELEASE_
	PushCallStack("mpi::Reduce");
#endif
	MPI_Reduce( sendbuf,  recvbuf, count, MPI_DOUBLE, op, root, comm );
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
}		// -----  end of function Reduce  ----- 

void
Reduce ( Complex* sendbuf, Complex* recvbuf, Int count, MPI_Op op, Int root, MPI_Comm comm )
{
#ifndef _RELEASE_
	PushCallStack("mpi::Reduce");
#endif
	MPI_Reduce( (Real*)sendbuf,  (Real*)recvbuf, 2 * count, MPI_DOUBLE, op, root, comm );
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
}		// -----  end of function Reduce  ----- 


void
Allreduce ( Real* sendbuf, Real* recvbuf, Int count, MPI_Op op, MPI_Comm comm )
{
#ifndef _RELEASE_
	PushCallStack("mpi::Allreduce");
#endif
	MPI_Allreduce( sendbuf,  recvbuf, count, MPI_DOUBLE, 
			op, comm );
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
}		// -----  end of function Allreduce  ----- 


void
Allreduce ( Complex* sendbuf, Complex* recvbuf, Int count, MPI_Op op, MPI_Comm comm )
{
#ifndef _RELEASE_
	PushCallStack("mpi::Allreduce");
#endif
	MPI_Allreduce( (Real*)sendbuf, (Real*) recvbuf, 2*count, MPI_DOUBLE, 
			op, comm );
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
}		// -----  end of function Allreduce  ----- 

} // namespace mpi

} // namespace PEXSI
