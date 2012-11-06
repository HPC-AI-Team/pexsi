#include  "environment_impl.hpp"
#include  "sparse_matrix.hpp"
#include  "numvec_impl.hpp"
#include  "utility.hpp"
//#include "pluselinv.hpp"
#include "superlu_ddefs.h"
#include "Cnames.h"

using namespace PEXSI;
using namespace std;

void Usage(){
  std::cout 
		<< "Test for factorization using SuperLU for real arithmetic with parallel input matrices and solve" << std::endl
		<< "ex10 -r [nprow] -c [npcol]" << std::endl;
}

// FIXME: IntNumVec convention.  Assumes a symmetric matrix
void DistSparseMatrixToSuperMatrixNRloc(SuperMatrix* ANRloc, DistSparseMatrix<Real>& A, gridinfo_t* grid){
	PushCallStack( "DistSparseMatrixToSuperMatrixNRloc" );
	Int mpirank = grid->iam;
	Int mpisize = grid->nprow * grid->npcol;

	Int numRowLocal = A.colptrLocal.m() - 1;
	Int numRowLocalFirst = A.size / mpisize;
	Int firstRow = mpirank * numRowLocalFirst;

  int_t *colindLocal, *rowptrLocal;
	double *nzvalLocal;
	rowptrLocal = (int_t*)intMalloc_dist(numRowLocal+1);
	colindLocal = (int_t*)intMalloc_dist(A.nnzLocal); 
	nzvalLocal  = (double*)doubleMalloc_dist(A.nnzLocal);
  
	std::copy( A.colptrLocal.Data(), A.colptrLocal.Data() + A.colptrLocal.m(),
			rowptrLocal );
	std::copy( A.rowindLocal.Data(), A.rowindLocal.Data() + A.rowindLocal.m(),
			colindLocal );
	std::copy( A.nzvalLocal.Data(), A.nzvalLocal.Data() + A.nzvalLocal.m(),
			nzvalLocal );

	std::cout << "Processor " << mpirank << " colindLocal[end] = " 
		<< colindLocal[A.nnzLocal-1] << " nzvalLocal[end] = " 
		<< nzvalLocal[A.nnzLocal-1] << std::endl;


	// Important to adjust from FORTRAN convention (1 based) to C convention (0 based) indices
	for(Int i = 0; i < A.rowindLocal.m(); i++){
		colindLocal[i]--;
	}

	for(Int i = 0; i < A.colptrLocal.m(); i++){
		rowptrLocal[i]--;
	}


	// Construct the distributed matrix according to the SuperLU_DIST format
	dCreate_CompRowLoc_Matrix_dist(ANRloc, A.size, A.size, A.nnzLocal, 
			numRowLocal, firstRow,
			nzvalLocal, colindLocal, rowptrLocal,
			SLU_NR_loc, SLU_D, SLU_GE);

	PopCallStack();
	return;
}


int main(int argc, char **argv) 
{
	MPI_Init(&argc, &argv);
	int mpirank, mpisize;
	MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
	MPI_Comm_size( MPI_COMM_WORLD, &mpisize );
	SetRandomSeed(1);

	if( argc < 5 ) {
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
		DistSparseMatrix<Real>  AMat;
		ScalePermstruct_t ScalePermstruct;
		LUstruct_t LUstruct;
		SOLVEstruct_t SOLVEstruct;
		gridinfo_t grid;
		std::string Hfile, Sfile;
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

		if(1){
			// Test code
			DistSparseMatrix<Real> HMat;
			DistSparseMatrix<Real> SMat;
			Real timeSta, timeEnd;
			GetTime( timeSta );
			superlu_gridinit(MPI_COMM_WORLD, nprow, npcol, &grid);
			ReadDistSparseMatrix( Hfile.c_str(), HMat, MPI_COMM_WORLD ); 
			ReadDistSparseMatrix( Sfile.c_str(), SMat, MPI_COMM_WORLD ); 
			GetTime( timeEnd );
			if( mpirank == 0 ){
				cout << "Time for reading H and S is " << timeEnd - timeSta << endl;
				cout << "H.size = " << HMat.size << endl;
				cout << "H.nnz  = " << HMat.nnz  << endl;
			}
			cout << "Processor " << mpirank << " outputs H.rowindLocal[end] = " 
				<< HMat.rowindLocal[HMat.nnzLocal-1] 
				<< " H.nzvalLocal[end] = " << HMat.nzvalLocal[HMat.nnzLocal-1]
				<< endl;
			cout << "Processor " << mpirank << " outputs S.rowindLocal[end] = " 
				<< SMat.rowindLocal[SMat.nnzLocal-1] 
				<< " S.nzvalLocal[end] = " << SMat.nzvalLocal[HMat.nnzLocal-1]
				<< endl;
			
			GetTime( timeSta );

			AMat.size   = HMat.size;
			AMat.nnz    = HMat.nnz;
			AMat.nnzLocal = HMat.nnzLocal;
			AMat.colptrLocal = HMat.colptrLocal;
			AMat.rowindLocal = HMat.rowindLocal;
			AMat.nzvalLocal.Resize( HMat.nnzLocal );
  
			Real *ptr0 = AMat.nzvalLocal.Data();
			Real *ptr1 = HMat.nzvalLocal.Data();
			Real *ptr2 = SMat.nzvalLocal.Data();
			for(Int i = 0; i < HMat.nnzLocal; i++){
				*(ptr0++) = *(ptr1++) - 3.5 * *(ptr2++);
//				*(ptr0++) = *(ptr1++);
//				*(ptr0++) = *(ptr1++) - Z_I * *(ptr2++);
			}

			DistSparseMatrixToSuperMatrixNRloc(&A, AMat, &grid);
			GetTime( timeEnd );
			if( mpirank == 0 )
				cout << "Time for converting the matrix A is " << timeEnd - timeSta << endl;
		}

 
		// Factorization without solve
		if(1){
			set_default_options_dist(&superlu_options);
//			superlu_options.Fact              = DOFACT;
			superlu_options.RowPerm           = NOROWPERM;
//			superlu_options.IterRefine        = NOREFINE;
//			superlu_options.ParSymbFact       = NO;
//			superlu_options.Equil             = NO; 
//			superlu_options.ReplaceTinyPivot  = NO;
//			superlu_options.ColPerm           = NATURAL; 
//			superlu_options.ColPerm           = METIS_AT_PLUS_A; 
			superlu_options.ColPerm           = MMD_AT_PLUS_A;
//			superlu_options.ColPerm = NATURAL;
//			superlu_options.PrintStat         = YES;
//			superlu_options.SolveInitialized  = NO;

			m = A.nrow;
			n = A.ncol;

			// Initialize ScalePermstruct and LUstruct.
			ScalePermstructInit(m, n, &ScalePermstruct);
			LUstructInit(m, n, &LUstruct);

			// Initialize the statistics variables.
			PStatInit(&stat);

			// Get the right hand side and the true solution
			
			int nrhs = 1;
			DblNumVec  xTrueGlobal(n), bGlobal(n);
			{
				SuperMatrix GA, A1;
				int needValue = 1;
				DistSparseMatrixToSuperMatrixNRloc(&A1, AMat, &grid);
				pdCompRow_loc_to_CompCol_global(needValue, &A1, &grid, &GA);
//				UniformRandom( xTrueGlobal );
				dGenXtrue_dist(n, 1, xTrueGlobal.Data(), n);
				char trans[1]; *trans='N';

				dFillRHS_dist(trans, 1, xTrueGlobal.Data(), n, &GA,
						bGlobal.Data(), n);
				Destroy_CompCol_Matrix_dist(&GA);
				Destroy_CompRowLoc_Matrix_dist(&A1);
			
//				if(mpirank == 0){
//					cerr << " x = " << xTrueGlobal << endl;;
//					cerr << " b = " << bGlobal << endl;
//				}
			}

			// Get the local solution
			double *bLocal, *xTrueLocal;
			Int numRowLocal = AMat.colptrLocal.m() - 1;
			Int numRowLocalFirst = AMat.size / mpisize;
			Int firstRow = mpirank * numRowLocalFirst;
			bLocal = doubleMalloc_dist( numRowLocal ); 
			cout << "Proc " << mpirank << " outputs numRowLocal = " <<
				numRowLocal << " firstRow = " << firstRow << endl;
			std::copy( bGlobal.Data()+firstRow, bGlobal.Data()+firstRow+numRowLocal,
					bLocal );
			xTrueLocal = doubleMalloc_dist( numRowLocal ); 
			std::copy( xTrueGlobal.Data()+firstRow, xTrueGlobal.Data()+firstRow+numRowLocal,
					xTrueLocal );
			cout << "Energy(bLocal) = " << Energy(PEXSI::DblNumVec(numRowLocal,false,bLocal)) << endl;
			
			// Call the linear equation solver. 
			double *berr=doubleMalloc_dist(nrhs);
			int     info;

			Real timeFactorSta, timeFactorEnd;
			GetTime( timeFactorSta );
			pdgssvx(&superlu_options, &A, &ScalePermstruct, bLocal, numRowLocal, nrhs, &grid,
					&LUstruct, &SOLVEstruct, berr, &stat, &info);
			GetTime( timeFactorEnd );

			if(mpirank == 0)
				cout << "Time for factorization is " << timeFactorEnd - timeFactorSta << endl;

			cout << "Energy(bLocal) = " << Energy(PEXSI::DblNumVec(numRowLocal,false,bLocal)) << endl;
			cout << "Energy(xLocal) = " << Energy(PEXSI::DblNumVec(numRowLocal,false,xTrueLocal)) << endl;

			if(mpirank == 0)
				cerr << PEXSI::DblNumVec(numRowLocal,false,bLocal) << endl;

			pdinf_norm_error(mpirank, ((NRformat_loc *)A.Store)->m_loc, nrhs, 
					bLocal, numRowLocal, xTrueLocal, numRowLocal, &grid);

			PStatPrint(&superlu_options, &stat, &grid);        /* Print the statistics. */
			PStatFree(&stat);
			SUPERLU_FREE( bLocal );
			SUPERLU_FREE( xTrueLocal );
			SUPERLU_FREE( berr );
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
			<< e.what() << std::endl;
		DumpCallStack();
	}
	
	MPI_Finalize();

	return 0;
}
