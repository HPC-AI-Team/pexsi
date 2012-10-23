#include "pexsi.hpp"
#include "pluselinv.hpp"

using namespace PEXSI;
using namespace std;

void Usage(){
  std::cout 
		<< "ex3" << std::endl;
// 	"-mu0 [mu0] -numel [numel] -deltaE [deltaE]" << std::endl
//		<< "mu0:    Initial guess for chemical potential" << std::endl
//		<< "numel:  Exact number of electrons (spin-restricted)" << std::endl
//		<< "deltaE: guess for the width of the spectrum of H-mu S" << std::endl;
}


int read_and_dist_csc(SuperMatrix *A, int nrhs, double **rhs,
		int *ldb, double **x, int *ldx,
		FILE *fp, gridinfo_t *grid);


int main(int argc, char **argv) 
{
	MPI_Init(&argc, &argv);
	int mpirank, mpisize;
	MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
	MPI_Comm_size( MPI_COMM_WORLD, &mpisize );


//	if( argc != 7 ) {
//		Usage();
//		MPI_Finalize();
//		return 0;
//	}
			
	try{
		stringstream  ss;
		ss << "logPEXSI";
		statusOFS.open( ss.str().c_str() );

		PEXSIData pexsiData;


		// *********************************************************************
		// Input parameter
		// *********************************************************************
		std::map<std::string,std::string> options;
		OptionsCreate(argc, argv, options);
		
//		pexsiData.gap              = 0.0;
//		pexsiData.temperature      = 300;
//		pexsiData.numPole          = 80;
//		pexsiData.permOrder        = -1;
//		pexsiData.numElectronTolerance = 1e-4;
//		pexsiData.muMaxIter        = 30;
//		pexsiData.poleTolerance    = 1e-4;
		// WaterPT
//		pexsiData.mu0              = -0.5;
//		pexsiData.numElectronExact = 1600.0;
//		pexsiData.deltaE           = 15.0;

		// DNA
//		pexsiData.mu0                = 0.00;
//		pexsiData.numElectronExact   = 2442.0;
//		pexsiData.deltaE           = 20.0;


//		if( options.find("-mu0") != options.end() ){
//			pexsiData.mu0 = std::atof(options["-mu0"].c_str());
//		}
//		else{
//      throw std::logic_error("mu0 must be provided.");
//		}
//
//    if( options.find("-numel") != options.end() ){
//			pexsiData.numElectronExact = std::atof(options["-numel"].c_str());
//		}
//		else{
//      throw std::logic_error("numel must be provided.");
//		}
//
//    if( options.find("-deltaE") != options.end() ){
//			pexsiData.deltaE = std::atof(options["-deltaE"].c_str());
//		}
//		else{
//      throw std::logic_error("deltaE must be provided.");
//		}


		// *********************************************************************
		// Read input matrix
		// *********************************************************************

		// Test code
		if(1){
			Real timeSta, timeEnd;
			GetTime( timeSta );
			ReadSparseMatrix( "H_LU.csc", pexsiData.HMat );
			GetTime( timeEnd );
			cout << pexsiData.HMat.size << endl;
			cout << pexsiData.HMat.nnz  << endl;
			cout << pexsiData.HMat.colptr.m() << endl;
			cout << "Time for reading H is " << timeEnd - timeSta << endl;
		}
 

//		if(0)
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
	
	MPI_Finalize();

	return 0;
}
