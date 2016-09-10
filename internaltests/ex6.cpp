#include  "environment_impl.hpp"
#include  "sparse_matrix_impl.hpp"
#include  "numvec_impl.hpp"
#include  "utility.hpp"
//#include "pluselinv.hpp"
#include "superlu_zdefs.h"
#include "Cnames.h"

using namespace PEXSI;
using namespace std;

void Usage(){
  std::cout 
		<< "Test for factorization using SuperLU for complex arithmetic with parallel input matrices" << std::endl
		<< "ex6 -r [nprow] -c [npcol]" << std::endl;
}

// FIXME: IntNumVec convention.  Assumes a symmetric matrix
	Int mpisize = grid->nprow * grid->npcol;

	Int numRowLocal = A.colptrLocal.m() - 1;
	Int numRowLocalFirst = A.size / mpisize;
	Int firstRow = mpirank * numRowLocalFirst;

  int_t *colindLocal, *rowptrLocal;
	doublecomplex *nzvalLocal;
	rowptrLocal = (int_t*)intMalloc_dist(numRowLocal+1);
	colindLocal = (int_t*)intMalloc_dist(A.nnzLocal); 
	nzvalLocal  = (doublecomplex*)doublecomplexMalloc_dist(A.nnzLocal);
  
	std::copy( A.colptrLocal.Data(), A.colptrLocal.Data() + A.colptrLocal.m(),
			rowptrLocal );
	std::copy( A.rowindLocal.Data(), A.rowindLocal.Data() + A.rowindLocal.m(),
			colindLocal );
	std::copy( A.nzvalLocal.Data(), A.nzvalLocal.Data() + A.nzvalLocal.m(),
			(Complex*)nzvalLocal );

//	std::cout << "Processor " << mpirank << " rowptrLocal[end] = " << 
//		rowptrLocal[numRowLocal] << std::endl;


	// Important to adjust from FORTRAN convention (1 based) to C convention (0 based) indices
	for(Int i = 0; i < A.rowindLocal.m(); i++){
		colindLocal[i]--;
	}

	for(Int i = 0; i < A.colptrLocal.m(); i++){
		rowptrLocal[i]--;
	}


	// Construct the distributed matrix according to the SuperLU_DIST format
	zCreate_CompRowLoc_Matrix_dist(ANRloc, A.size, A.size, A.nnzLocal, 
			numRowLocal, firstRow,
			nzvalLocal, colindLocal, rowptrLocal,
			SLU_NR_loc, SLU_Z, SLU_GE);
}


int main(int argc, char **argv) 
{
	MPI_Init(&argc, &argv);
	int mpirank, mpisize;
	MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
	MPI_Comm_size( MPI_COMM_WORLD, &mpisize );

	if( argc != 5 ) {
		Usage();
		MPI_Finalize();
		return 0;
	}
			
	try{
		stringstream  ss;
		ss << "logPLUSelInv";
		statusOFS.open( ss.str().c_str() );

		superlu_options_t superlu_options;
		SuperLUStat_t stat;
		SuperMatrix A;
		DistSparseMatrix<Complex>  AMat;
		ScalePermstruct_t ScalePermstruct;
		LUstruct_t LUstruct;
		SOLVEstruct_t SOLVEstruct;
		gridinfo_t grid;
		int      nprow, npcol;
		int      m, n;

		// *********************************************************************
		// Input parameter
		// *********************************************************************
		std::map<std::string,std::string> options;
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
		

		// *********************************************************************
		// Read input matrix
		// *********************************************************************

		if(1){
			// Test code
			DistSparseMatrix<Real> HMat;
			DistSparseMatrix<Real> SMat;
			Real timeSta, timeEnd;
			GetTime( timeSta );
			superlu_gridinit(MPI_COMM_WORLD, nprow, npcol, &grid);
			ReadDistSparseMatrix( "H_LU.csc", HMat, MPI_COMM_WORLD ); 
			ReadDistSparseMatrix( "S_LU.csc", SMat, MPI_COMM_WORLD ); 
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

			DistSparseMatrixToSuperMatrixNRloc(&A, AMat, &grid);
			GetTime( timeEnd );
			if( mpirank == 0 )
				cout << "Time for converting the matrix A is " << timeEnd - timeSta << endl;
		}

 
		// Factorization without solve
		if(1){
			set_default_options_dist(&superlu_options);
			superlu_options.Fact              = DOFACT;
			superlu_options.RowPerm           = NOROWPERM;
			superlu_options.IterRefine        = NOREFINE;
			superlu_options.ParSymbFact       = NO;
			superlu_options.Equil             = NO; 
			superlu_options.ReplaceTinyPivot  = NO;
//			superlu_options.ColPerm           = PARMETIS; 
			superlu_options.ColPerm           = MMD_AT_PLUS_A;
//			superlu_options.ColPerm = NATURAL;
			superlu_options.PrintStat         = YES;
			superlu_options.SolveInitialized  = NO;

			m = A.nrow;
			n = A.ncol;

			// Initialize ScalePermstruct and LUstruct.
			ScalePermstructInit(m, n, &ScalePermstruct);
			LUstructInit(m, n, &LUstruct);

			// Initialize the statistics variables.
			PStatInit(&stat);

			// Call the linear equation solver. 
			doublecomplex *b=NULL; 
			double *berr=NULL;
			int nrhs = 0;
			int      info;

			Real timeFactorSta, timeFactorEnd;
			GetTime( timeFactorSta );
			pzgssvx(&superlu_options, &A, &ScalePermstruct, b, n, nrhs, &grid,
					&LUstruct, &SOLVEstruct, berr, &stat, &info);
			GetTime( timeFactorEnd );

			if(mpirank == 0)
				cout << "Time for factorization is " << timeFactorEnd - timeFactorSta << endl;

			PStatPrint(&superlu_options, &stat, &grid);        /* Print the statistics. */
			PStatFree(&stat);
		}

		// *********************************************************************
		// Test the accuracy of factorization by solve
		// *********************************************************************

		if(1){
			Real timeSta, timeEnd;
			// Construct the global matrix A from the distributed matrix A
			SuperMatrix GA;
			SuperMatrix A1;
			GetTime( timeSta );
			DistSparseMatrixToSuperMatrixNRloc(&A1, AMat, &grid);
			if(mpirank == 0) 
				cout << "Sparse matrix generated" << endl;
			int needValue = 1;
			pzCompRow_loc_to_CompCol_global(needValue, &A1, &grid, &GA);
			if(mpirank == 0) 
				cout << "Sparse matrix converted" << endl;
			
			// Generate the true solution and the right hand side
			CpxNumVec  xTrueGlobal(n), bGlobal(n);
			zGenXtrue_dist(n, 1, (doublecomplex*)xTrueGlobal.Data(), n);
			char trans[1]; *trans='N';

			zFillRHS_dist(trans, 1, (doublecomplex*)xTrueGlobal.Data(), n, &GA,
					(doublecomplex*)bGlobal.Data(), n);
			if(mpirank == 0){
				cout << "||xTrue|| = " << Energy(xTrueGlobal) << endl;
				cout << "||b    || = " << Energy(bGlobal) << endl;
//				statusOFS << xTrueGlobal << endl << bGlobal << endl;
			}
		
      // Initialize the solve
	    zSolveInit(&superlu_options, &A, (ScalePermstruct.perm_r), (ScalePermstruct.perm_c), 
					1, &LUstruct, &grid, &SOLVEstruct);

			GetTime( timeEnd );
			if( mpirank == 0 )
				cout << "Time for preparing solve is " << timeEnd - timeSta << endl;

			// Solve
			GetTime( timeSta );
			PStatInit(&stat);
			int info;
			pzgstrs_Bglobal(n, &LUstruct, &grid, (doublecomplex*)bGlobal.Data(), 
					n, 1, &stat, &info);
			PStatFree(&stat);
			GetTime( timeEnd );
			if( mpirank == 0 )
				cout << "Time for solve is " << timeEnd - timeSta << endl;
			
			if(mpirank == 0){
				cout << "||xTrue|| = " << Energy(xTrueGlobal) << endl;
				cout << "||b    || = " << Energy(bGlobal) << endl;
//				statusOFS << bGlobal << endl;
			}

			pzinf_norm_error(mpirank, ((NRformat_loc *)A.Store)->m_loc, 1, 
					(doublecomplex*)bGlobal.Data(), n, 
					(doublecomplex*)xTrueGlobal.Data(), n, &grid);


			// Destroy A and GA
			Destroy_CompRowLoc_Matrix_dist(&A1);
			Destroy_CompCol_Matrix_dist(&GA);
			// Free the solve
			zSolveFinalize(&superlu_options, &SOLVEstruct);
		}


		// *********************************************************************
		// Deallocate the storage
		// *********************************************************************
		if(1){
			Destroy_CompRowLoc_Matrix_dist(&A);
			ScalePermstructFree(&ScalePermstruct);
			Destroy_LU(n, &grid, &LUstruct);
			LUstructFree(&LUstruct);
			superlu_gridexit(&grid);
		}

		statusOFS.close();
	}
	catch( std::exception& e )
	{
		std::cerr << "Processor " << mpirank << " caught exception with message: "
	
	MPI_Finalize();

	return 0;
}
