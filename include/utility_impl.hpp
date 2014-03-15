using namespace std;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::cerr;

namespace PEXSI{

  template<typename T>
void
GetDiagonal ( const DistSparseMatrix<T>& A, 
		NumVec<T>& diag )
{
#ifndef _RELEASE_
	PushCallStack("GetDiagonal");
#endif
	Int mpirank, mpisize;
	MPI_Comm_rank( A.comm, &mpirank );
	MPI_Comm_size( A.comm, &mpisize );

  NumVec<T>	 diagLocal( A.size );
	SetValue( diagLocal, ZERO<T>() );
	diag.Resize( A.size );
	SetValue( diag, ZERO<T>() );

	Int numColFirst = A.size / mpisize;
	Int firstCol    = mpirank * numColFirst;
	Int numColLocal = A.colptrLocal.m() - 1;

#if ( _DEBUGlevel_ >= 1 )
	statusOFS << "numColFirst = " << numColFirst << std::endl;
	statusOFS << "A.nzvalLocal.size = " << A.nzvalLocal.m() << std::endl;
	statusOFS << "A.nnzLocal = " << A.nnzLocal << std::endl;
#endif

	// Note that the indices in DistSparseMatrix follows the FORTRAN convention
  for( Int j = 0; j < numColLocal; j++ ){
		Int jcol = j + firstCol + 1;
		Int numRow = A.colptrLocal(j+1) - A.colptrLocal(j);
		const Int* rowPtr = &A.rowindLocal( A.colptrLocal(j) - 1 );
		// NOTE: The rows in DistSparseMatrix are not necessarily ordered.
		// So lower_bound cannot be used here for fast searching. find has to be used. 
		const Int* ptr = find( rowPtr, rowPtr + numRow, jcol ); 
		if( ptr == rowPtr + numRow ){
			std::ostringstream msg;
			msg << "Serious problem. Did not find the row corresponding to the column." << std::endl
				<< "This happens when j = " << j << ", jcol = " << jcol << ", and the row indices are " << std::endl
				<< IntNumVec( numRow, false, const_cast<Int*>(rowPtr) ) << std::endl;
			throw std::logic_error( msg.str().c_str() );
		}
		Int diagIdx = ptr - A.rowindLocal.Data();
    diagLocal( jcol - 1 ) = A.nzvalLocal( diagIdx );
	}

	mpi::Allreduce( &diagLocal[0], &diag[0], A.size, MPI_SUM, A.comm );

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
}		// -----  end of function GetDiagonal  ----- 


}
