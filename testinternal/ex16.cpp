/// @file ex16.cpp
/// @brief Check step-by-step the correctness of PSelInv.
/// @author Lin Lin
/// @version 0.1
/// @date 2012-12-19
#include  "environment_impl.hpp"
#include  "sparse_matrix_impl.hpp"
#include  "numvec_impl.hpp"
#include  "utility.hpp"
#include  "superlu_dist_interf.hpp"
#include	"pselinv.hpp"

using namespace PEXSI;
using namespace std;

void Usage(){
  std::cout 
		<< "Usage" << std::endl
		<< "ex16 -H [Hfile] -S [Sfile] -colperm [colperm]" << std::endl;
}

int main(int argc, char **argv) 
{
	if( argc < 5 ) {
		Usage();
		return 0;
	}
	
	MPI_Init( &argc, &argv );
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
		std::map<std::string,std::string> options;
		std::string Hfile, Sfile;                   

		OptionsCreate(argc, argv, options);
	  Int nprow = iround( std::sqrt( (double)mpisize) );
		Int npcol = mpisize / nprow;
		if( mpisize != nprow * npcol || nprow != npcol ){
			throw std::runtime_error( "nprow == npcol is assumed in this test routine." );
		}

		if( mpirank == 0 )
			cout << "nprow = " << nprow << ", npcol = " << npcol << endl;
		
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

		std::string ColPerm;
		if( options.find("-colperm") != options.end() ){ 
			ColPerm = options["-colperm"];
		}
		else{
      throw std::logic_error("colperm must be provided.");
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
		AMat.comm = MPI_COMM_WORLD;

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
		SuperLUOptions luOpt;
		luOpt.ColPerm = ColPerm;
		SuperLUMatrix luMat( g, luOpt );
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
		Real timeTotalFactorizationSta, timeTotalFactorizationEnd; 
	
		
		// Important: the distribution in pzsymbfact is going to mess up the
		// A matrix.  Recompute the matrix A here.
		luMat.DistSparseMatrixToSuperMatrixNRloc( AMat );

		GetTime( timeTotalFactorizationSta );
    
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

		GetTime( timeTotalFactorizationEnd );
		if( mpirank == 0 )
			cout << "Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " sec" << endl; 


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
			Real timeTotalSelInvSta, timeTotalSelInvEnd;
			GetTime( timeTotalSelInvSta );

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
			statusOFS << "permInv: " << endl << super.permInv << endl;
			statusOFS << "superIdx:" << endl << super.superIdx << endl;
			statusOFS << "superPtr:" << endl << super.superPtr << endl; 

			// Convert the LU factors into DistSparseMatrix
			{
				DistSparseMatrix<Scalar> Ainv;
				PMloc.PMatrixToDistSparseMatrix( Ainv );
				ofstream ofs("LU");
				if( !ofs.good() ) 
					throw std::runtime_error("file cannot be opened.");
				serialize( Ainv, ofs, NO_MASK );
				ofs.close();
			}





			// Preparation for the selected inversion
			GetTime( timeSta );
			PMloc.ConstructCommunicationPattern();
			GetTime( timeEnd );

			if( mpirank == 0 )
				cout << "Time for constructing the communication pattern is " << timeEnd  - timeSta << endl;


			GetTime( timeSta );
			PMloc.PreSelInv();
			GetTime( timeEnd );
			if( mpirank == 0 )
				cout << "Time for pre-selected inversion is " << timeEnd  - timeSta << endl;


			// Convert the matrix after PreSelInv into DistSparseMatrix
			{
				DistSparseMatrix<Scalar> Ainv;
				PMloc.PMatrixToDistSparseMatrix( Ainv );
				ofstream ofs("PreInv");
				if( !ofs.good() ) 
					throw std::runtime_error("file cannot be opened.");
				serialize( Ainv, ofs, NO_MASK );
				ofs.close();
			}

			// Main subroutine for selected inversion
			GetTime( timeSta );
			PMloc.SelInv();
			GetTime( timeEnd );

			if( mpirank == 0 )
				cout << "Time for numerical selected inversion is " << timeEnd  - timeSta << endl;


			// Convert the matrix after PreSelInv into DistSparseMatrix
			{
				DistSparseMatrix<Scalar> Ainv;
				PMloc.PMatrixToDistSparseMatrix( Ainv );
				ofstream ofs("Ainv");
				if( !ofs.good() ) 
					throw std::runtime_error("file cannot be opened.");
				serialize( Ainv, ofs, NO_MASK );
				ofs.close();
			}
//
//
//			GetTime( timeTotalSelInvEnd );
//			if( mpirank == 0 )
//				cout << "Time for total selected inversion is " << timeTotalSelInvEnd  - timeTotalSelInvSta << endl;
//
//			NumVec<Scalar> diag;
//			
//			GetTime( timeSta );
//			PMloc.GetDiagonal( diag );
//			GetTime( timeEnd );
//			if( mpirank == 0 )
//				cout << "Time for getting the diagonal is " << timeEnd  - timeSta << endl;
//
//			if( mpirank == 0 ){
//				statusOFS << std::endl << "Diagonal of inverse in natural order: " << std::endl << diag << std::endl;
//				ofstream ofs("diag");
//				if( !ofs.good() ) 
//					throw std::runtime_error("file cannot be opened.");
//				serialize( diag, ofs, NO_MASK );
//				ofs.close();
//			}
//
//
//
//
//			// Convert to DistSparseMatrix and get the diagonal
//			GetTime( timeSta );
//			DistSparseMatrix<Scalar> Ainv;
//			PMloc.PMatrixToDistSparseMatrix( Ainv );
//			GetTime( timeEnd );
//			
//			if( mpirank == 0 )
//				cout << "Time for converting PMatrix to DistSparseMatrix is " << timeEnd  - timeSta << endl;
//
//			NumVec<Scalar> diagDistSparse;
//			GetTime( timeSta );
//			GetDiagonal( Ainv, diagDistSparse );
//			GetTime( timeEnd );
//			if( mpirank == 0 )
//				cout << "Time for getting the diagonal of DistSparseMatrix is " << timeEnd  - timeSta << endl;
//
//			if( mpirank == 0 ){
//				statusOFS << std::endl << "Diagonal of inverse from DistSparseMatrix format : " << std::endl << diagDistSparse << std::endl;
//				Real diffNorm = 0.0;;
//				for( Int i = 0; i < diag.m(); i++ ){
//					diffNorm += pow( std::abs( diag(i) - diagDistSparse(i) ), 2.0 );
//				}
//				diffNorm = std::sqrt( diffNorm );
//				statusOFS << std::endl << "||diag - diagDistSparse||_2 = " << diffNorm << std::endl;
//			}
//
//			// Convert to DistSparseMatrix in the 2nd format and get the diagonal
//			GetTime( timeSta );
//			DistSparseMatrix<Scalar> Ainv2;
//			PMloc.PMatrixToDistSparseMatrix( AMat, Ainv2 );
//			GetTime( timeEnd );
//			
//			if( mpirank == 0 )
//				cout << "Time for converting PMatrix to DistSparseMatrix (2nd format) is " << timeEnd  - timeSta << endl;
//
//			NumVec<Scalar> diagDistSparse2;
//			GetTime( timeSta );
//			GetDiagonal( Ainv2, diagDistSparse2 );
//			GetTime( timeEnd );
//			if( mpirank == 0 )
//				cout << "Time for getting the diagonal of DistSparseMatrix is " << timeEnd  - timeSta << endl;
//
//			if( mpirank == 0 ){
//				statusOFS << std::endl << "Diagonal of inverse from the 2nd conversion into DistSparseMatrix format : " << std::endl << diagDistSparse2 << std::endl;
//				Real diffNorm = 0.0;;
//				for( Int i = 0; i < diag.m(); i++ ){
//					diffNorm += pow( std::abs( diag(i) - diagDistSparse2(i) ), 2.0 );
//				}
//				diffNorm = std::sqrt( diffNorm );
//				statusOFS << std::endl << "||diag - diagDistSparse2||_2 = " << diffNorm << std::endl;
//			}
//
//			Complex traceLocal = blas::Dot( AMat.nnzLocal, AMat.nzvalLocal.Data(), 1, 
//					Ainv2.nzvalLocal.Data(), 1 );
//			Complex trace = Z_ZERO;
//			mpi::Allreduce( &traceLocal, &trace, 1, MPI_SUM, MPI_COMM_WORLD );
//
//			if( mpirank == 0 )
//				statusOFS << std::endl << "Tr[Ainv2 * AMat] = " << std::endl << trace << std::endl;
//

		}
		
		statusOFS.close();
	}
	catch( std::exception& e )
	{
		std::cerr << "Processor " << mpirank << " caught exception with message: "
			<< e.what() << std::endl;
	}
	
	MPI_Finalize();

	return 0;
}
