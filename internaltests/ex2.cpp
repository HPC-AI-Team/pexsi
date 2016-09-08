// *********************************************************************
// Solve the H x = lambda S x with ScaLAPACK.  Temporarily stopped.
// *********************************************************************

#include "pexsi.hpp"
#include  "lapack.hpp"

using namespace PEXSI;
using namespace std;

int main(int argc, char **argv) 
{
	MPI_Init(&argc, &argv);
	int mpirank, mpisize;
	MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
	MPI_Comm_size( MPI_COMM_WORLD, &mpisize );
	
	try{
		stringstream  ss;
		ss << "logPEXSI";
		statusOFS.open( ss.str().c_str() );

		// *********************************************************************
		// Read input matrix
		// *********************************************************************

		// Test code
 
		SparseMatrix<Real>  HMat, SMat;

		if(1)
		{
			ReadSparseMatrix( "H.ccs", HMat );
			ReadSparseMatrix( "S.ccs", SMat );

			// Make sure that the sparsity of H and S matches.
			if( HMat.size != SMat.size ||
					HMat.nnz  != SMat.nnz ){
					std::ostringstream msg;
					msg 
						<< "The dimensions colptr for H and S do not match" << std::endl
						<< "H.colptr.size = " << HMat.size << std::endl
						<< "H.colptr.nnz  = " << HMat.nnz  << std::endl
						<< "S.colptr.size = " << SMat.size << std::endl
						<< "S.colptr.nnz  = " << SMat.nnz  << std::endl;
				throw std::logic_error( msg.str().c_str() );
			}

      for( int j = 0; j < HMat.colptr.m(); j++ ){
				if( HMat.colptr(j) != SMat.colptr(j) ){
					std::ostringstream msg;
					msg 
						<< "Colptr of H and S do not match:" << std::endl
						<< "H.colptr(" << j << ") = " << HMat.colptr(j) << std::endl
						<< "S.colptr(" << j << ") = " << SMat.colptr(j) << std::endl;
					throw std::logic_error( msg.str().c_str() );	
				}
			}
		}

		// *********************************************************************
		// Solve
		// *********************************************************************

		cout << SMat.size << endl;

		cout << SMat.nzval << endl;

		lapack::Cholesky( 'L', SMat.size, SMat.nzval.Data(), SMat.size );
		statusOFS << "Cholesky finished" << endl;

		lapack::Hegst( 1, 'L', HMat.size, HMat.nzval.Data(), HMat.size,
				SMat.nzval.Data(), SMat.size );

		statusOFS << "Hegst finished" << endl;

		statusOFS.close();
	}
	catch( std::exception& e )
	{
		std::cerr << " caught exception with message: "
			<< e.what() << std::endl;
		DumpCallStack();
	}
	
	MPI_Finalize();

	return 0;
}
