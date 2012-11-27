/// @file ex15.cpp
/// @brief Test for the PPEXSI module using SuperLU and PSelInv.
/// @author Lin Lin
/// @version 0.1
/// @date 2012-11-26
#include "ppexsi.hpp"

using namespace PEXSI;
using namespace std;

void Usage(){
  std::cout 
		<< "ex15 -mu0 [mu0] -numel [numel] -deltaE [deltaE] -H [Hfile] -S [Sfile] -npPerPole [npPole]" << std::endl
		<< "mu0:    Initial guess for chemical potential" << std::endl
		<< "numel:  Exact number of electrons (spin-restricted)" << std::endl
		<< "deltaE: guess for the width of the spectrum of H-mu S" << std::endl
		<< "H: Hamiltonian matrix (csc format, both lower triangular and upper triangular)" << std::endl
		<< "S: Overlap     matrix (csc format, both lower triangular and upper triangular)" << std::endl
		<< "npPerPole: number of processors used for each pole" << std::endl;
}

int main(int argc, char **argv) 
{
	MPI_Init(&argc, &argv);
	int mpirank, mpisize;
	MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
	MPI_Comm_size( MPI_COMM_WORLD, &mpisize );


	if( argc < 7 || argc%2 == 0 ) {
		if( mpirank == 0 ) Usage();
		MPI_Finalize();
		return 0;
	}
			

	
	try{
		stringstream  ss;
		ss << "logPEXSI" << mpirank;
		statusOFS.open( ss.str().c_str() );


		// *********************************************************************
		// Input parameter
		// *********************************************************************
		std::map<std::string,std::string> options;
		OptionsCreate(argc, argv, options);
		
		Real gap              = 0.0;
		Real temperature      = 300;
		Real numPole          = 80;
		Real numElectronTolerance = 1e-4;
		Real muMaxIter        = 30;
//		// WaterPT
////		pexsiData.mu0              = -0.5;
////		pexsiData.numElectronExact = 1600.0;
////		pexsiData.deltaE           = 15.0;
//
//		// DNA
////		pexsiData.mu0                = 0.00;
////		pexsiData.numElectronExact   = 2442.0;
////		pexsiData.deltaE           = 20.0;
//
//
		Real mu0;
		if( options.find("-mu0") != options.end() ){
			mu0 = std::atof(options["-mu0"].c_str());
		}
		else{
      throw std::logic_error("mu0 must be provided.");
		}

		Real numElectronExact;
    if( options.find("-numel") != options.end() ){
			numElectronExact = std::atof(options["-numel"].c_str());
		}
		else{
      throw std::logic_error("numel must be provided.");
		}

		Real deltaE;
    if( options.find("-deltaE") != options.end() ){
			deltaE = std::atof(options["-deltaE"].c_str());
		}
		else{
      throw std::logic_error("deltaE must be provided.");
		}

		Int npPerPole;
    if( options.find("-npPerPole") != options.end() ){
			npPerPole = std::atoi(options["-npPerPole"].c_str());
		}
		else{
      throw std::logic_error("npPerPole must be provided.");
		}
   

		// *********************************************************************
		// Check the input parameters
		// *********************************************************************
		if( mpisize % npPerPole != 0 ){
			throw std::logic_error( "mpisize cannot be divided evenly by npPerPole." );
		}



//
//
//		// *********************************************************************
//		// Read input matrix
//		// *********************************************************************
//
//		// Test code
//		if(0){
//			Real timeSta, timeEnd;
//			GetTime( timeSta );
//			ReadSparseMatrix( "H.ccs", pexsiData.HMat );
//			GetTime( timeEnd );
//			cout << pexsiData.HMat.size << endl;
//			cout << pexsiData.HMat.nnz  << endl;
//			cout << pexsiData.HMat.colptr.m() << endl;
//			cout << "Time for reading H is " << timeEnd - timeSta << endl;
//		}
// 
//
//		if(1)
//		{
//			ReadSparseMatrix( "H.ccs", pexsiData.HMat );
//			ReadSparseMatrix( "S.ccs", pexsiData.SMat );
//
//			// Make sure that the sparsity of H and S matches.
//			if( pexsiData.HMat.size != pexsiData.SMat.size ||
//					pexsiData.HMat.nnz  != pexsiData.SMat.nnz ){
//					std::ostringstream msg;
//					msg 
//						<< "The dimensions colptr for H and S do not match" << std::endl
//						<< "H.colptr.size = " << pexsiData.HMat.size << std::endl
//						<< "H.colptr.nnz  = " << pexsiData.HMat.nnz  << std::endl
//						<< "S.colptr.size = " << pexsiData.SMat.size << std::endl
//						<< "S.colptr.nnz  = " << pexsiData.SMat.nnz  << std::endl;
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
//
//
//		Print(statusOFS, "mu0                    = ", pexsiData.mu0);
//		Print(statusOFS, "numElectronExact       = ", pexsiData.numElectronExact);
//		Print(statusOFS, "deltaE                 = ", pexsiData.deltaE);
//		Print(statusOFS, "gap                    = ", pexsiData.gap);
//		Print(statusOFS, "temperature            = ", pexsiData.temperature);
//		Print(statusOFS, "numPole                = ", pexsiData.numPole);
//		Print(statusOFS, "numElectronTolerance   = ", pexsiData.numElectronTolerance);
//		Print(statusOFS, "poleTolerance          = ", pexsiData.poleTolerance);
//		Print(statusOFS, "muMaxIter              = ", pexsiData.muMaxIter);
//		Print(statusOFS, "permOrder              = ", pexsiData.permOrder);
//
//
//
//
//
//		// *********************************************************************
//		// Solve
//		// *********************************************************************
//		pexsiData.Setup( );
//		Print( statusOFS, "PEXSI setup finished." );
//
//		pexsiData.Solve( );
//		Print( statusOFS, "PEXSI solve finished." );
//
//		statusOFS.close();
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
