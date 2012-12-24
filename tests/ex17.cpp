/// @file ex17.cpp
/// @brief Test for the PEXSI module using SelInv.
///
/// This is an updated version for ex1.cpp.
/// @author Lin Lin
/// @date 2012-12-23
#include "pexsi.hpp"

using namespace PEXSI;
using namespace std;

void Usage(){
  std::cout 
		<< "ex17 -mu0 [mu0] -numel [numel] -numPole [numPole] -deltaE [deltaE] -H [Hfile] -S [Sfile] -muiter [muiter]" << std::endl
		<< "mu0:    Initial guess for chemical potential" << std::endl
		<< "numel:  Exact number of electrons (spin-restricted)" << std::endl
		<< "numPole: Number of poles." << std::endl
		<< "deltaE: guess for the width of the spectrum of H-mu S" << std::endl
		<< "H: Hamiltonian matrix (csc format, both lower triangular and upper triangular)" << std::endl
		<< "S: Overlap     matrix (csc format, both lower triangular and upper triangular)" << std::endl
	  << "muiter:  number of iterations for the chemical potential" << std::endl;
}

int main(int argc, char **argv) 
{
	MPI_Init(&argc, &argv);
	int mpirank, mpisize;
	MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
	MPI_Comm_size( MPI_COMM_WORLD, &mpisize );


	if( argc < 15 || argc%2 == 0 ) {
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
		Real poleTolerance    = 1e-12;
		Real numElectronTolerance = 1e-4;

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
		
		Int numPole;
    if( options.find("-numPole") != options.end() ){
			numPole = std::atof(options["-numPole"].c_str());
		}
		else{
      throw std::logic_error("numPole must be provided.");
		}


		Real deltaE;
    if( options.find("-deltaE") != options.end() ){
			deltaE = std::atof(options["-deltaE"].c_str());
		}
		else{
      throw std::logic_error("deltaE must be provided.");
		}

		std::string Hfile, Sfile;                   
		if( options.find("-H") != options.end() ){ 
			Hfile = options["-H"];
		}
		else{
      throw std::logic_error("Hfile must be provided.");
		}

		if( options.find("-S") != options.end() ){ 
			Sfile = options["-S"];
		}
		else{
      throw std::logic_error("Sfile must be provided.");
		}

		Int muMaxIter;
    if( options.find("-muiter") != options.end() ){
			muMaxIter = std::atoi(options["-muiter"].c_str());
		}
		else{
      throw std::logic_error("muiter must be provided.");
		}


		bool isFreeEnergyDensityMatrix = true;

		bool isEnergyDensityMatrix     = true;

		// *********************************************************************
		// Check the input parameters
		// *********************************************************************
		if( mpisize != 1 ){
			throw std::logic_error( "mpisize must be 1." );
		}

		std::string ColPerm = "MMD_AT_PLUS_A";

		// Initialize
		PEXSIData pexsi;

		// *********************************************************************
		// Read input matrix
		// *********************************************************************


		SparseMatrix<Real> HMat;
		SparseMatrix<Real> SMat;

		Real timeSta, timeEnd;

		// HMat
		ReadSparseMatrix( Hfile.c_str(), HMat );
		// SMat
		ReadSparseMatrix( Sfile.c_str(), SMat );
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

		Print(statusOFS, "mu0                    = ", mu0);
		Print(statusOFS, "numElectronExact       = ", numElectronExact);
		Print(statusOFS, "deltaE                 = ", deltaE);
		Print(statusOFS, "gap                    = ", gap);
		Print(statusOFS, "temperature            = ", temperature);
		Print(statusOFS, "numPole                = ", numPole);
		Print(statusOFS, "poleTolerance          = ", poleTolerance);
		Print(statusOFS, "numElectronTolerance   = ", numElectronTolerance);
		Print(statusOFS, "ColPerm                = ", ColPerm );
		Print(statusOFS, "muMaxIter              = ", muMaxIter);
		Print(statusOFS, "mpisize                = ", mpisize );
		Print(statusOFS, "isFreeEnergyMatrix     = ", isFreeEnergyDensityMatrix );
		Print(statusOFS, "isEnergyMatrix         = ", isEnergyDensityMatrix ); 


		// *********************************************************************
		// Solve
		// *********************************************************************
		std::vector<Real>  muList;
		std::vector<Real>  numElectronList;
		
		Real timeSolveSta, timeSolveEnd;

		GetTime( timeSolveSta );
		pexsi.Solve( 
				numPole,
				temperature,
				numElectronExact,
				gap,
				deltaE,
				mu0,
				HMat,
				SMat,
				muMaxIter,
				poleTolerance,
				numElectronTolerance,
				ColPerm,
				isFreeEnergyDensityMatrix,
				isEnergyDensityMatrix,
				muList,
				numElectronList );

		GetTime( timeSolveEnd );

		statusOFS << "Time for solving the PEXSI is " <<
			timeSolveEnd - timeSolveSta << " [s]" << std::endl;

		statusOFS.close();
	}
	catch( std::exception& e )
	{
		std::cerr << " caught exception with message: "
			<< e.what() << std::endl;
#ifndef _RELEASE_
		DumpCallStack();
#endif
	}
	
	MPI_Finalize();

	return 0;
}
