/// @file ex13.cpp
/// @brief Test for the new input interface for SuperLU with complex
/// arithmetic.
/// @author Lin Lin
/// @version 0.1
/// @date 2012-11-15
#include  "environment_impl.hpp"
#include  "sparse_matrix.hpp"
#include  "numvec_impl.hpp"
#include  "utility.hpp"
#include  "superlu_dist_interf.hpp"

using namespace PEXSI;
using namespace std;

void Usage(){
  std::cout 
		<< "Usage" << std::endl
		<< "ex13 -r [nprow] -c [npcol]" << std::endl;
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
		ss << "logTest";
		statusOFS.open( ss.str().c_str() );


//
//		// For selected inversion
//		PMatrix PMloc;
//		std::vector<std::vector<Int> > localEtree; 
//
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
			*(ptr0++) = *(ptr1++) - Z_I * *(ptr2++);
		}
		GetTime( timeEnd );
		if( mpirank == 0 )
			cout << "Time for constructing the matrix A is " << timeEnd - timeSta << endl;

		GetTime( timeSta );
		SuperLUMatrix luMat;
		luMat.Initialize( MPI_COMM_WORLD, nprow, npcol );
		luMat.DistSparseMatrixToSuperMatrixNRloc( AMat );
		GetTime( timeEnd );
		if( mpirank == 0 )
			cout << "Time for converting to SuperLU format is " << timeEnd - timeSta << endl;
		
		// *********************************************************************
		// Symbolic factorization 
		// *********************************************************************

 
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

		SuperLUGrid g( MPI_COMM_WORLD, nprow, npcol );

		// *********************************************************************
		// Test the accuracy of factorization by solve
		// *********************************************************************

//		if(0){
//			PStatInit(&stat);
//			superlu_options.Fact              = FACTORED;
//
//			int nrhs = 1;
//			PEXSI::CpxNumVec  xTrueGlobal(n), bGlobal(n);
//			{
//				SuperMatrix GA, A1;
//				int needValue = 1;
//				DistSparseMatrixToSuperMatrixNRloc(&A1, AMat, &grid);
//				pzCompRow_loc_to_CompCol_global(needValue, &A1, &grid, &GA);
//				UniformRandom( xTrueGlobal );
//				//				zGenXtrue_dist(n, 1, (doublecomplex*)xTrueGlobal.Data(), n);
//				char trans[1]; *trans='N';
//
//				zFillRHS_dist(trans, 1, (doublecomplex*)xTrueGlobal.Data(), n, &GA,
//						(doublecomplex*)bGlobal.Data(), n);
//				Destroy_CompCol_Matrix_dist(&GA);
//				Destroy_CompRowLoc_Matrix_dist(&A1);
//
//				//				if(mpirank == 0){
//				//					cerr << " x = " << xTrueGlobal << endl;;
//				//					cerr << " b = " << bGlobal << endl;
//				//				}
//			}
//
//			// Get the local solution
//			doublecomplex *bLocal, *xTrueLocal;
//			Int numRowLocal = AMat.colptrLocal.m() - 1;
//			Int numRowLocalFirst = AMat.size / mpisize;
//			Int firstRow = mpirank * numRowLocalFirst;
//			bLocal = doublecomplexMalloc_dist( numRowLocal ); 
//			cout << "Proc " << mpirank << " outputs numRowLocal = " <<
//				numRowLocal << " firstRow = " << firstRow << endl;
//			std::copy( bGlobal.Data()+firstRow, bGlobal.Data()+firstRow+numRowLocal,
//					(Complex*)bLocal );
//			xTrueLocal = doublecomplexMalloc_dist( numRowLocal ); 
//			std::copy( xTrueGlobal.Data()+firstRow, xTrueGlobal.Data()+firstRow+numRowLocal,
//					(Complex*)xTrueLocal );
////			cout << "Energy(bLocal) = " << Energy(PEXSI::CpxNumVec(numRowLocal,false,(Complex*)bLocal)) << endl;
//
//			// Call the linear equation solver. 
//			double *berr=doubleMalloc_dist(nrhs);
//
//			GetTime( timeFactorSta );
//			pzgssvx(&superlu_options, &A, &ScalePermstruct, bLocal, numRowLocal, nrhs, &grid,
//					&LUstruct, &SOLVEstruct, berr, &stat, &info);
//			GetTime( timeFactorEnd );
//
////			cout << "Energy(bLocal) = " << Energy(PEXSI::CpxNumVec(numRowLocal,false,(Complex*)bLocal)) << endl;
////			cout << "Energy(xLocal) = " << Energy(PEXSI::CpxNumVec(numRowLocal,false,(Complex*)xTrueLocal)) << endl;
//
//
//			pzinf_norm_error(mpirank, ((NRformat_loc *)A.Store)->m_loc, nrhs, 
//					bLocal, numRowLocal, xTrueLocal, numRowLocal, &grid);
//
//			PStatPrint(&superlu_options, &stat, &grid);        /* Print the statistics. */
//			PStatFree(&stat);
//			SUPERLU_FREE( bLocal );
//			SUPERLU_FREE( xTrueLocal );
//			SUPERLU_FREE( berr );
//		}

//
//		// *********************************************************************
//		// Selected inversion
//		// *********************************************************************
//		
//		{
//			GetTime(timeSta);
//			SuperLU2SelInv(n, &LUstruct, &grid, PMloc);
//			GetTime(timeEnd);
//			if( mpirank == 0 )
//				cout << "Time for converting the SuperLU to SelInv format is " 
//					<< timeEnd - timeSta << " sec" << endl; 
//
//			GetTime(timeSta);
//			ConstructLocalEtree(n, &grid, PMloc, localEtree);
//			GetTime(timeEnd);
//			if( mpirank == 0 )
//				cout << "Time for converting the SelInv elimination tree is " 
//					<< timeEnd - timeSta << " sec" << endl; 
//			
//			GetTime(timeSta);
//			PLUSelInv(&grid, PMloc, localEtree);
//			GetTime(timeEnd);
//			if( mpirank == 0 )
//				cout << "Time for the selected inversion is " 
//					<< timeEnd - timeSta << " sec" << endl; 
//		}
//
//		// *********************************************************************
//		// Deallocate the storage
//		// *********************************************************************
//		Destroy_CompRowLoc_Matrix_dist(&A);
//		Destroy_LU(n, &grid, &LUstruct);
//		ScalePermstructFree(&ScalePermstruct);
//		LUstructFree(&LUstruct);
//		superlu_gridexit(&grid);

		luMat.Finalize();
		statusOFS.close();
	}
	catch( std::exception& e )
	{
		std::cerr << "Processor " << mpirank << " caught exception with message: "
			<< e.what() << std::endl;
		DumpCallStack();
	}
	
	MPI_Finalize();

	return 0;
}
