#include "pexsi.hpp"

using namespace PEXSI;
using namespace std;

int main(int argc, char **argv) 
{
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
		pexsiData.deltaE           = 1.0;
		pexsiData.numPole          = 10;
		pexsiData.permOrder        = -1;
		pexsiData.mu0              = -1.0;
		pexsiData.numElectronExact = 6;
		pexsiData.numElectronTolerance = 1e-4;
		pexsiData.muMaxIter        = 10;
		pexsiData.permOrder        = -1;


		// *********************************************************************
		// Read input matrix
		// *********************************************************************
		ReadSparseMatrix( "H.ccs", pexsiData.HMat );
    cout << pexsiData.HMat.size << endl;
    cout << pexsiData.HMat.nnz  << endl;
 

//		{
//			ReadSparseMatrix( "H.matrix", pexsiData.HMat );
//
//			ReadSparseMatrix( "S.matrix", pexsiData.SMat );
//
//			// Make sure that the sparsity of H and S matches.
//			if( pexsiData.HMat.size != pexsiData.SMat.size ||
//					pexsiData.HMat.nnz  != pexsiData.SMat.nnz ){
//					std::ostringstream msg;
//					msg 
//						<< "The dimensions colptr for H and S do not match" << std::endl
//						<< "H.colptr.size = " pexsiData.HMat.size << std::endl
//						<< "H.colptr.nnz  = " pexsiData.HMat.nnz  << std::endl
//						<< "S.colptr.size = " pexsiData.SMat.size << std::endl
//						<< "S.colptr.nnz  = " pexsiData.SMat.nnz  << std::endl;
//				throw std::logic_error( msg.str().c_str() );
//			}
//
//      for( int j = 0; j < pexsiData.HMat.colptr.m(); j++ ){
//				if( pexsiData.HMat.colptr(j) != pexsiData.SMat.colptr(j) ){
//					std::ostringstream msg;
//					msg 
//						<< "Colptr of H and S do not match:" << std::endl
//						<< "H.colptr(" << j << ") = " << pexsiData.HMat.colptr(j) << std::endl
//						<< "S.colptr(" << j << ") = " << pexsiData.SMat.colptr(j) << std::endl;
//					throw std::logic_error( msg.str().c_str() );	
//				}
//			}
//		}


		// *********************************************************************
		// Solve
		// *********************************************************************
//		pexsiData.Setup( );
//		Print( statusOFS, "PEXSI setup finished." );
//
//		pexsiData.Solve( );
//		Print( statusOFS, "PEXSI solve finished." );

		statusOFS.close();
	}
	catch( std::exception& e )
	{
		std::cerr << " caught exception with message: "
			<< e.what() << std::endl;
		DumpCallStack();
	}

	return 0;
}
