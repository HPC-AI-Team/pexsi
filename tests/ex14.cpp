/// @file ex14.cpp
/// @brief Test for the new input interface for SuperLU together with
/// the new parallel selected inversion.
/// @author Lin Lin
/// @version 0.1
/// @date 2012-11-16
#include  "environment_impl.hpp"
#include  "sparse_matrix.hpp"
#include  "numvec_impl.hpp"
#include  "utility.hpp"
#include  "superlu_dist_interf.hpp"
#include	"pselinv.hpp"

using namespace PEXSI;
using namespace std;

void Usage(){
  std::cout 
		<< "Usage" << std::endl
		<< "ex14 -r [nprow] -c [npcol]" << std::endl;
}

int main(int argc, char **argv) 
{
	if( argc < 5 ) {
		Usage();
		return 0;
	}
	
	MPI_Init(&argc, &argv);
	int mpirank, mpisize;
	MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
	MPI_Comm_size( MPI_COMM_WORLD, &mpisize );

			
	try{
		stringstream  ss;
		ss << "logTest" << mpirank;
		statusOFS.open( ss.str().c_str() );

		// *********************************************************************
		// Input parameter
		// *********************************************************************
		int      nprow, npcol;
		std::map<std::string,std::string> options;
		std::string Hfile, Sfile;                   

		OptionsCreate(argc, argv, options);
		if( options.find("-r") != options.end() ){ 
			nprow = std::atoi(options["-r"].c_str());
		}
		else{
      throw std::logic_error("nprow must be provided.");
		}

		if( options.find("-c") != options.end() ){ 
			npcol = std::atoi(options["-c"].c_str());
		}
		else{
      throw std::logic_error("npcol must be provided.");
		}
		
		if( options.find("-H") != options.end() ){ 
			Hfile = options["-H"];
		}
		else{
			Hfile = "H_LU.csc";
		}

		if( options.find("-S") != options.end() ){ 
			Sfile = options["-S"];
		}
		else{
			Sfile = "S_LU.csc";
		}


		// *********************************************************************
		// Read input matrix
		// *********************************************************************

		int      m, n;
		DistSparseMatrix<Complex>  AMat;

		DistSparseMatrix<Real> HMat;
		DistSparseMatrix<Real> SMat;
		Real timeSta, timeEnd;
		GetTime( timeSta );
		ReadDistSparseMatrix( Hfile.c_str(), HMat, MPI_COMM_WORLD ); 
		ReadDistSparseMatrix( Sfile.c_str(), SMat, MPI_COMM_WORLD ); 
		GetTime( timeEnd );
		if( mpirank == 0 ){
			cout << "Time for reading H and S is " << timeEnd - timeSta << endl;
			cout << "H.size = " << HMat.size << endl;
			cout << "H.nnz  = " << HMat.nnz  << endl;
		}

		GetTime( timeSta );

		AMat.size   = HMat.size;
		AMat.nnz    = HMat.nnz;
		AMat.nnzLocal = HMat.nnzLocal;
		AMat.colptrLocal = HMat.colptrLocal;
		AMat.rowindLocal = HMat.rowindLocal;
		AMat.nzvalLocal.Resize( HMat.nnzLocal );

		Complex *ptr0 = AMat.nzvalLocal.Data();
		Real *ptr1 = HMat.nzvalLocal.Data();
		Real *ptr2 = SMat.nzvalLocal.Data();
		for(Int i = 0; i < HMat.nnzLocal; i++){
			*(ptr0++) = *(ptr1++);// - Z_I * *(ptr2++);
		}
		GetTime( timeEnd );
		if( mpirank == 0 )
			cout << "Time for constructing the matrix A is " << timeEnd - timeSta << endl;

		
		// *********************************************************************
		// Symbolic factorization 
		// *********************************************************************

		// Setup grid.
		SuperLUGrid g( MPI_COMM_WORLD, nprow, npcol );

		GetTime( timeSta );
		SuperLUMatrix luMat( g );
		luMat.DistSparseMatrixToSuperMatrixNRloc( AMat );
		GetTime( timeEnd );
		if( mpirank == 0 )
			cout << "Time for converting to SuperLU format is " << timeEnd - timeSta << endl;

		GetTime( timeSta );
		luMat.SymbolicFactorize();
		luMat.DestroyAOnly();
		GetTime( timeEnd );

		if( mpirank == 0 )
			cout << "Time for performing the symbolic factorization is " << timeEnd - timeSta << endl;


		// *********************************************************************
		// Numerical factorization only 
		// *********************************************************************
		
		// Important: the distribution in pzsymbfact is going to mess up the
		// A matrix.  Recompute the matrix A here.
		luMat.DistSparseMatrixToSuperMatrixNRloc( AMat );
    
		GetTime( timeSta );
		luMat.Distribute();
		GetTime( timeEnd );
		if( mpirank == 0 )
			cout << "Time for distribution is " << timeEnd - timeSta << " sec" << endl; 
		
     		 

		GetTime( timeSta );
		luMat.NumericalFactorize();
		GetTime( timeEnd );

		if( mpirank == 0 )
			cout << "Time for factorization is " << timeEnd - timeSta << " sec" << endl; 


		// *********************************************************************
		// Test the accuracy of factorization by solve
		// *********************************************************************

		if( 1 ) {
			SuperLUMatrix A1( g ), GA( g );
			A1.DistSparseMatrixToSuperMatrixNRloc( AMat );
			A1.ConvertNRlocToNC( GA );

			int n = A1.n();
			int nrhs = 5;
			CpxNumMat xTrueGlobal(n, nrhs), bGlobal(n, nrhs);
			CpxNumMat xTrueLocal, bLocal;
			DblNumVec berr;
			UniformRandom( xTrueGlobal );
			
			GA.MultiplyGlobalMultiVector( xTrueGlobal, bGlobal );

			A1.DistributeGlobalMultiVector( xTrueGlobal, xTrueLocal );
			A1.DistributeGlobalMultiVector( bGlobal,     bLocal );

			luMat.SolveDistMultiVector( bLocal, berr );
			luMat.CheckErrorDistMultiVector( bLocal, xTrueLocal );
		}





		// *********************************************************************
		// Selected inversion
		// *********************************************************************

		if( 1 )
		{
			Grid g1( MPI_COMM_WORLD, nprow, npcol );
			SuperNode super;
			
			GetTime( timeSta );
			luMat.SymbolicToSuperNode( super );
			PMatrix PMloc( &g1, &super );
			luMat.LUstructToPMatrix( PMloc );
			GetTime( timeEnd );

			if( mpirank == 0 )
				cout << "Time for converting LUstruct to PMatrix is " << timeEnd  - timeSta << endl;

			statusOFS << "perm: " << endl << super.perm << endl;
			statusOFS << "superIdx:" << endl << super.superIdx << endl;
			statusOFS << "superPtr:" << endl << super.superPtr << endl; 


			GetTime( timeSta );
			PMloc.ConstructCommunicationPattern();
			GetTime( timeEnd );

			if( mpirank == 0 )
				cout << "Time for constructing the communication pattern is " << timeEnd  - timeSta << endl;

			GetTime( timeSta );
			PMloc.SelInv();
			GetTime( timeEnd );

			if( mpirank == 0 )
				cout << "Time for numerical selected inversion is " << timeEnd  - timeSta << endl;

			NumVec<Scalar> diag;
			
			GetTime( timeSta );
			PMloc.Diagonal( diag );
			GetTime( timeEnd );
			if( mpirank == 0 )
				cout << "Time for getting the diagonal is " << timeEnd  - timeSta << endl;

			if( mpirank == 0 )
				statusOFS << std::endl << "Diagonal of inverse in natural order: " << std::endl << diag << std::endl;

		}
		
		statusOFS.close();
	}
	catch( std::exception& e )
	{
		std::cerr << "Processor " << mpirank << " caught exception with message: "
			<< e.what() << std::endl;
#ifndef _RELEASE_
		DumpCallStack();
#endif
	}
	
	MPI_Finalize();

	return 0;
}
