/// @file run_ppexsi.cpp
/// @brief Test for the PPEXSI module using SuperLU and PSelInv.
///
/// @author Lin Lin
/// @date 2012-12-24
#include "ppexsi.hpp"

using namespace PEXSI;
using namespace std;

void Usage(){
  std::cout 
		<< "run_ppexsi -mu0 [mu0] -muMin [muMin] -muMax [muMax] -numel [numel] -numPole [numPole] -deltaE [deltaE] -gap [gap] -H [Hfile] -S [Sfile] -npPerPole [npPole] -colperm [colperm] -muiter [muiter]" << std::endl << "temp:    Temperature (unit: K)" << std::endl
		<< "mu0:     Initial guess for chemical potential" << std::endl
		<< "muMin:   Lower bound for chemical potential" << std::endl
		<< "muMax:   Upper bound for chemical potential" << std::endl
		<< "numel:   Exact number of electrons (spin-restricted)" << std::endl
		<< "numPole: Number of poles." << std::endl
		<< "deltaE:  guess for the width of the spectrum of H-mu S" << std::endl
		<< "gap:     guess for the distance betweeen the spectrum and mu" << std::endl
		<< "H: Hamiltonian matrix (csc format, both lower triangular and upper triangular)" << std::endl
		<< "S: Overlap     matrix (csc format, both lower triangular and upper triangular)" << std::endl
		<< "npPerPole: number of processors used for each pole" << std::endl
		<< "colperm: permutation method (for SuperLU_DIST)" << std::endl
	  << "muiter:  number of iterations for the chemical potential" << std::endl;
}

int main(int argc, char **argv) 
{
	MPI_Init(&argc, &argv);
	int mpirank, mpisize;
	MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
	MPI_Comm_size( MPI_COMM_WORLD, &mpisize );


	if( argc < 27 || argc%2 == 0 ) {
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
		
		Real poleTolerance    = 1e-30;
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

		Real temperature;
		if( options.find("-temp") != options.end() ){
			temperature = std::atof(options["-temp"].c_str());
		}
		else{
      throw std::logic_error("temp must be provided.");
		}

		Real mu0;
		if( options.find("-mu0") != options.end() ){
			mu0 = std::atof(options["-mu0"].c_str());
		}
		else{
      throw std::logic_error("mu0 must be provided.");
		}

		Real muMin;
		if( options.find("-muMin") != options.end() ){
			muMin = std::atof(options["-muMin"].c_str());
		}
		else{
      throw std::logic_error("muMin must be provided.");
		}

		Real muMax;
		if( options.find("-muMax") != options.end() ){
			muMax = std::atof(options["-muMax"].c_str());
		}
		else{
      throw std::logic_error("muMax must be provided.");
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

		Real gap;
    if( options.find("-gap") != options.end() ){
			gap = std::atof(options["-gap"].c_str());
		}
		else{
      throw std::logic_error("gap must be provided.");
		}


		Int npPerPole;
    if( options.find("-npPerPole") != options.end() ){
			npPerPole = std::atoi(options["-npPerPole"].c_str());
		}
		else{
      throw std::logic_error("npPerPole must be provided.");
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

		std::string ColPerm;
		if( options.find("-colperm") != options.end() ){ 
			ColPerm = options["-colperm"];
		}
		else{
      throw std::logic_error("colperm must be provided.");
		}

		bool isFreeEnergyDensityMatrix = true;

		bool isEnergyDensityMatrix     = true;

		bool isDerivativeTMatrix       = true;

		// *********************************************************************
		// Check the input parameters
		// *********************************************************************
		if( mpisize % npPerPole != 0 ){
			throw std::logic_error( "mpisize cannot be divided evenly by npPerPole." );
		}

	  Int nprow = iround( std::sqrt( (double)npPerPole) );
		Int npcol = npPerPole / nprow;
		if( npPerPole != nprow * npcol || nprow != npcol ){
			throw std::runtime_error( "npPerPole must be a perfect square due to the current implementation of PSelInv." );
		}

		Int muMaxIter;
    if( options.find("-muiter") != options.end() ){
			muMaxIter = std::atoi(options["-muiter"].c_str());
		}
		else{
      throw std::logic_error("muiter must be provided.");
		}


		// Initialize
		Grid gridPole( MPI_COMM_WORLD, mpisize / npPerPole, npPerPole );
		PPEXSIData pexsi( &gridPole, nprow, npcol );

		// *********************************************************************
		// Read input matrix
		// *********************************************************************


		DistSparseMatrix<Real> HMat;
		DistSparseMatrix<Real> SMat;

		Real timeSta, timeEnd;

		// The first processor row gridPole reads first, and then
		// broadcast the information to other row processors in gridPole.
		// This can be improved later by MPI_Bcast directly.

		// HMat
		std::vector<char> sstr;
		Int sizeStm;
		if( MYROW( &gridPole ) == 0 ){
			std::stringstream sstm;
			ReadDistSparseMatrix( Hfile.c_str(), HMat, gridPole.rowComm ); 
			serialize( HMat.size, sstm, NO_MASK );
			serialize( HMat.nnz,  sstm, NO_MASK );
			serialize( HMat.nnzLocal, sstm, NO_MASK );
			serialize( HMat.colptrLocal, sstm, NO_MASK );
			serialize( HMat.rowindLocal, sstm, NO_MASK );
			serialize( HMat.nzvalLocal,  sstm, NO_MASK );
			sstr.resize( Size( sstm ) );
		  sstm.read( &sstr[0], sstr.size() ); 	
			sizeStm = sstr.size();
		}
		MPI_Bcast( &sizeStm, 1, MPI_INT, 0, gridPole.colComm );
		statusOFS << "sizeStm = " << sizeStm << std::endl;

		if( MYROW( &gridPole ) != 0 ) sstr.resize( sizeStm );

		MPI_Bcast( (void*)&sstr[0], sizeStm, MPI_BYTE, 0, gridPole.colComm );

		if( MYROW( &gridPole ) != 0 ){
			std::stringstream sstm;
			sstm.write( &sstr[0], sizeStm );
			deserialize( HMat.size, sstm, NO_MASK );
			deserialize( HMat.nnz,  sstm, NO_MASK );
			deserialize( HMat.nnzLocal, sstm, NO_MASK );
			deserialize( HMat.colptrLocal, sstm, NO_MASK );
			deserialize( HMat.rowindLocal, sstm, NO_MASK );
			deserialize( HMat.nzvalLocal,  sstm, NO_MASK );
		}
		statusOFS << "size     = " << HMat.size     << std::endl;
		statusOFS << "nnzLocal = " << HMat.nnzLocal << std::endl;
		// Communicator
		HMat.comm = gridPole.rowComm;

		sstr.clear();

		// SMat
		if( MYROW( &gridPole ) == 0 ){
			std::stringstream sstm;
			ReadDistSparseMatrix( Sfile.c_str(), SMat, gridPole.rowComm ); 
			serialize( SMat.size, sstm, NO_MASK );
			serialize( SMat.nnz,  sstm, NO_MASK );
			serialize( SMat.nnzLocal, sstm, NO_MASK );
			serialize( SMat.colptrLocal, sstm, NO_MASK );
			serialize( SMat.rowindLocal, sstm, NO_MASK );
			serialize( SMat.nzvalLocal,  sstm, NO_MASK );
			sstr.resize( Size( sstm ) );
		  sstm.read( &sstr[0], sstr.size() ); 	
			sizeStm = sstr.size();
		}

		MPI_Bcast( &sizeStm, 1, MPI_INT, 0, gridPole.colComm );

		if( MYROW( &gridPole ) != 0 ) sstr.resize( sizeStm );

		MPI_Bcast( (void*)&sstr[0], sizeStm, MPI_BYTE, 0, gridPole.colComm );

		if( MYROW( &gridPole ) != 0 ){
			std::stringstream sstm;
			sstm.write( &sstr[0], sizeStm );
			deserialize( SMat.size, sstm, NO_MASK );
			deserialize( SMat.nnz,  sstm, NO_MASK );
			deserialize( SMat.nnzLocal, sstm, NO_MASK );
			deserialize( SMat.colptrLocal, sstm, NO_MASK );
			deserialize( SMat.rowindLocal, sstm, NO_MASK );
			deserialize( SMat.nzvalLocal,  sstm, NO_MASK );
		}
		// Communicator
		SMat.comm = gridPole.rowComm;

		sstr.clear();



		Print(statusOFS, "mu0                    = ", mu0);
		Print(statusOFS, "muMin                  = ", muMin);
		Print(statusOFS, "muMax                  = ", muMax); 
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
		Print(statusOFS, "npPerPole              = ", npPerPole );
		Print(statusOFS, "isFreeEnergyMatrix     = ", isFreeEnergyDensityMatrix );
		Print(statusOFS, "isEnergyMatrix         = ", isEnergyDensityMatrix ); 


		// *********************************************************************
		// Solve
		// *********************************************************************
		std::vector<Real>  muList;
		std::vector<Real>  numElectronList;
		std::vector<Real>  numElectronDrvList;
	  bool isConverged;	
		
		Real timeSolveSta, timeSolveEnd;

		GetTime( timeSolveSta );
		pexsi.Solve( 
				numPole,
				temperature,
				numElectronExact,
				gap,
				deltaE,
				mu0,
				muMin,
				muMax,
				HMat,
				SMat,
				muMaxIter,
				poleTolerance,
				numElectronTolerance,
				ColPerm,
				isFreeEnergyDensityMatrix,
				isEnergyDensityMatrix,
				isDerivativeTMatrix,
				muList,
				numElectronList,
				numElectronDrvList,
			  isConverged	);

		GetTime( timeSolveEnd );

		PrintBlock( statusOFS, "Solve finished." );
		if( isConverged ){
			statusOFS << "PEXSI has converged with " << muList.size() << 
				" iterations" << std::endl;
		}
		else {
			statusOFS << "PEXSI did not converge with " << muList.size() << 
				" iterations" << std::endl;
		}
		Print( statusOFS, "mu                   = ", 
				*muList.rbegin() );
		Print( statusOFS, "Total time for PEXSI = ", 
				timeSolveEnd - timeSolveSta );


		
		// *********************************************************************
		// Solve for the second
		// Using temperature expansion for the chemical potential
		// *********************************************************************

		if(1){
			Real muZeroT = pexsi.EstimateZeroTemperatureChemicalPotential(
					temperature,
					*muList.rbegin(),
					SMat );

			PrintBlock( statusOFS, "Second calculation at low temperature." );

			Print( statusOFS, "mu (T=0)             = ", 
					muZeroT );

			temperature = 300.0;

			mu0   = (muZeroT < *muList.rbegin()) ? muZeroT : *muList.rbegin();
			muMin = mu0 - 0.02;
			muMax = *muList.rbegin();
		
			GetTime( timeSolveSta );

			pexsi.Solve( 
					numPole,
					temperature,
					numElectronExact,
					gap,
					deltaE,
					mu0,
					muMin,
					muMax,
					HMat,
					SMat,
					muMaxIter,
					poleTolerance,
					numElectronTolerance,
					ColPerm,
					isFreeEnergyDensityMatrix,
					isEnergyDensityMatrix,
					isDerivativeTMatrix,
					muList,
					numElectronList,
					numElectronDrvList,
					isConverged	);

			GetTime( timeSolveEnd );

			PrintBlock( statusOFS, "Solve finished." );
			if( isConverged ){
				statusOFS << "PEXSI has converged with " << muList.size() << 
					" iterations" << std::endl;
			}
			else {
				statusOFS << "PEXSI did not converge with " << muList.size() << 
					" iterations" << std::endl;
			}
			Print( statusOFS, "mu                   = ", 
					*muList.rbegin() );
			Print( statusOFS, "Total time for PEXSI = ", 
					timeSolveEnd - timeSolveSta );

		}

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
