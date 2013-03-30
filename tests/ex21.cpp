/// @file ex21.cpp
/// @brief Test for the inertia count.
///
/// @author Lin Lin
/// @date 2013-03-26
#include "ppexsi.hpp"

extern "C" {
  double seekeig_(int *, int *, double *, double *, double *);
}
using namespace PEXSI;
using namespace std;

void Usage(){
  std::cout 
		<< "ex21 -muMin [muMin] -muMax [muMax] -numPole [numPole] -H [Hfile] -S [Sfile] -npPerPole [npPole] -colperm [colperm] -formatted [formatted]" 
		<< "muMin:   Lower bound for shift" << std::endl
		<< "muMax:   Upper bound for shift" << std::endl
		<< "numPole: Number of shifts (poles)." << std::endl
		<< "H: Hamiltonian matrix " << std::endl
		<< "S: Overlap     matrix. if omitted, the overlap matrix is treated as an identity matrix implicitly." << std::endl
		<< "npPerPole: number of processors used for each pole" << std::endl
		<< "colperm: permutation method (for SuperLU_DIST)" << std::endl
		<< "formatted: whether the input of H/S matrices are formatted (1) or unformatted (csc format, 0)" << std::endl;
}

int main(int argc, char **argv) 
{
	MPI_Init(&argc, &argv);
	int mpirank, mpisize;
	MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
	MPI_Comm_size( MPI_COMM_WORLD, &mpisize );


	if( argc < 13 || argc%2 == 0 ) {
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

		
		Int numPole;
    if( options.find("-numPole") != options.end() ){
			numPole = std::atof(options["-numPole"].c_str());
		}
		else{
      throw std::logic_error("numPole must be provided.");
		}


		Int npPerPole;
    if( options.find("-npPerPole") != options.end() ){
			npPerPole = std::atoi(options["-npPerPole"].c_str());
		}
		else{
      throw std::logic_error("npPerPole must be provided.");
		}

		Int isFormatted;
		if( options.find("-formatted") != options.end() ){
			isFormatted = std::atoi(options["-formatted"].c_str());
		}
		else{
			isFormatted = 0; // Binary file
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
			statusOFS << "-S option is not given. " 
				<< "Treat the overlap matrix as an identity matrix." 
				<< std::endl << std::endl;
		}

		std::string ColPerm;
		if( options.find("-colperm") != options.end() ){ 
			ColPerm = options["-colperm"];
		}
		else{
      throw std::logic_error("colperm must be provided.");
		}
		
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

		// Initialize
		Grid gridPole( MPI_COMM_WORLD, mpisize / npPerPole, npPerPole );
		PPEXSIData pexsi( &gridPole, nprow, npcol );

		// *********************************************************************
		// Read input matrix
		// *********************************************************************

		DistSparseMatrix<Real> HMat;

		Real timeSta, timeEnd;

		// The first processor row gridPole reads first, and then
		// broadcast the information to other row processors in gridPole.
		// This can be improved later by MPI_Bcast directly.

		// HMat
		std::vector<char> sstr;
		Int sizeStm;
		if( MYROW( &gridPole ) == 0 ){
			std::stringstream sstm;


			if( isFormatted )
				ReadDistSparseMatrixFormatted( Hfile.c_str(), HMat, gridPole.rowComm ); 
			else
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
		DistSparseMatrix<Real> SMat;
		if( Sfile.empty() ){
			// Set the size to be zero.  This will tell PPEXSI.Solve to treat
			// the overlap matrix as an identity matrix implicitly.
			SMat.size = 0;  
		}
		else{
			// SMat is given directly
			if( MYROW( &gridPole ) == 0 ){
				std::stringstream sstm;
				if( isFormatted )
					ReadDistSparseMatrixFormatted( Sfile.c_str(), SMat, gridPole.rowComm ); 
				else
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
		} // if (Sfile.empty())


		Print(statusOFS, "muMin                  = ", muMin);
		Print(statusOFS, "muMax                  = ", muMax); 
		Print(statusOFS, "numPole                = ", numPole);
		Print(statusOFS, "ColPerm                = ", ColPerm );
		Print(statusOFS, "mpisize                = ", mpisize );
		Print(statusOFS, "npPerPole              = ", npPerPole );


		// *********************************************************************
		// Solve
		// *********************************************************************

		Int  numShift = numPole;
		std::vector<Real>  shiftVec( numShift );
		std::vector<Int>   inertiaVec( numShift );

		for( Int l = 0; l < numShift; l++ ){
			shiftVec[l] = muMin + l * (muMax - muMin) / (numShift-1);
		}

		GetTime( timeSta );
		pexsi.CalculateNegativeInertia( 
				shiftVec,
				inertiaVec,
				HMat,
				SMat,
				ColPerm );

		GetTime( timeEnd );

		PrintBlock( statusOFS, "Inertia count finished." );

		for( Int l = 0; l < numShift; l++ ){
			statusOFS << std::setiosflags(std::ios::left) 
				<< std::setw(LENGTH_VAR_NAME) << "Shift = "
				<< std::setw(LENGTH_VAR_DATA) << shiftVec[l]
				<< std::setw(LENGTH_VAR_NAME) << "Inertia = "
				<< std::setw(LENGTH_VAR_DATA) << inertiaVec[l]
				<< std::endl;
		}

		if( mpirank == 0 ) {
         double mu = 0.0;
         int  nelec = 20;
         double *xs, *ys;

         xs = (double*)malloc(numShift*sizeof(double));
         ys = (double*)malloc(numShift*sizeof(double));

         for (int i = 0; i <numShift; i++) {
            xs[i] = (double)shiftVec[i];
            ys[i] = (double)inertiaVec[i];
            printf("x = %11.3e, y = %11.3e\n", xs[i], ys[i]);
         }

         double mu0 = (muMin  + muMax)/2.0;
         mu = seekeig_(&nelec, &numShift, xs, ys, &mu0); 

         free(xs);
         free(ys); 
      }
		
		Print( statusOFS, "Total time = ", 
				timeEnd - timeSta );

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
