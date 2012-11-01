#include  "environment_impl.hpp"
#include  "sparse_matrix.hpp"
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
void DistSparseMatrixToSuperMatrixNRloc(SuperMatrix* ANRloc, DistSparseMatrix<Complex>& A, gridinfo_t* grid){
	PushCallStack( "DistSparseMatrixToSuperMatrixNRloc" );
	Int mpirank = grid->iam;
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

	PopCallStack();
	return;
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
			DistSparseMatrix<Complex>  AMat;
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

//			Real timeFactorSta, timeFactorEnd;
//			GetTime( timeFactorSta );
			pzgssvx(&superlu_options, &A, &ScalePermstruct, b, n, nrhs, &grid,
					&LUstruct, &SOLVEstruct, berr, &stat, &info);
//			GetTime( timeFactorEnd );

//			cout << "Time for factorization is " << timeFactorEnd - timeFactorSta << endl;

			PStatPrint(&superlu_options, &stat, &grid);        /* Print the statistics. */
		}



		// *********************************************************************
		// Deallocate the storage
		// *********************************************************************
		if(1){
			PStatFree(&stat);
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
			<< e.what() << std::endl;
		DumpCallStack();
	}
	
	MPI_Finalize();

	return 0;
}
