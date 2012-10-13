#include "pexsi.hpp"

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

		PEXSIData pexsiData;


		// *********************************************************************
		// Input parameter
		// *********************************************************************
		pexsiData.gap              = 0.0;
		pexsiData.temperature      = 300;
		pexsiData.deltaE           = 15.0;
		pexsiData.numPole          = 80;
		pexsiData.permOrder        = -1;
		pexsiData.numElectronTolerance = 1e-4;
		pexsiData.muMaxIter        = 30;
		pexsiData.permOrder        = -1;
		pexsiData.poleTolerance    = 1e-4;

		// WaterPT
//		pexsiData.mu0              = -0.5;
//		pexsiData.numElectronExact = 1600.0;

		// DNA
		pexsiData.mu0                = 0.00;
		pexsiData.numElectronExact   = 2442.0;

		// *********************************************************************
		// Read input matrix
		// *********************************************************************

		// Test code
		if(0){
			Real timeSta, timeEnd;
			GetTime( timeSta );
			ReadSparseMatrix( "H.ccs", pexsiData.HMat );
			GetTime( timeEnd );
			cout << pexsiData.HMat.size << endl;
			cout << pexsiData.HMat.nnz  << endl;
			cout << pexsiData.HMat.colptr.m() << endl;
			cout << "Time for reading H is " << timeEnd - timeSta << endl;
		}
 

		if(1)
		{
			ReadSparseMatrix( "H.ccs", pexsiData.HMat );
			ReadSparseMatrix( "S.ccs", pexsiData.SMat );

			// Make sure that the sparsity of H and S matches.
			if( pexsiData.HMat.size != pexsiData.SMat.size ||
					pexsiData.HMat.nnz  != pexsiData.SMat.nnz ){
					std::ostringstream msg;
					msg 
						<< "The dimensions colptr for H and S do not match" << std::endl
						<< "H.colptr.size = " << pexsiData.HMat.size << std::endl
						<< "H.colptr.nnz  = " << pexsiData.HMat.nnz  << std::endl
						<< "S.colptr.size = " << pexsiData.SMat.size << std::endl
						<< "S.colptr.nnz  = " << pexsiData.SMat.nnz  << std::endl;
				throw std::logic_error( msg.str().c_str() );
			}

      for( int j = 0; j < pexsiData.HMat.colptr.m(); j++ ){
				if( pexsiData.HMat.colptr(j) != pexsiData.SMat.colptr(j) ){
					std::ostringstream msg;
					msg 
						<< "Colptr of H and S do not match:" << std::endl
						<< "H.colptr(" << j << ") = " << pexsiData.HMat.colptr(j) << std::endl
						<< "S.colptr(" << j << ") = " << pexsiData.SMat.colptr(j) << std::endl;
					throw std::logic_error( msg.str().c_str() );	
				}
			}
		}


		// *********************************************************************
		// Solve
		// *********************************************************************
		pexsiData.Setup( );
		Print( statusOFS, "PEXSI setup finished." );

		pexsiData.Solve( );
		Print( statusOFS, "PEXSI solve finished." );

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
