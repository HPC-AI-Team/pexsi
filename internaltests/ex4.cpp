//#include "pluselinv.hpp"
#include "superlu_zdefs.h"
#include "Cnames.h"
//#include "pexsi.hpp" 

using namespace PEXSI;
using namespace std;

void Usage(){
  std::cout 
		<< "Test for factorization using SuperLU for complex arithmetic" << std::endl
		<< "ex4 -r [nprow] -c [npcol]" << std::endl;
}

// FIXME: IntNumVec convention.  Assumes a symmetric matrix
	doublecomplex  *nzval_loc;         /* local */
	int_t    *colind_loc, *rowptr_loc;	 /* local */
	int_t    m_loc, fst_row, nnz_loc;
	int_t    m_loc_fst; /* Record m_loc of the first p-1 processors,
												 when mod(m, p) is not zero. */ 
	int_t    iam;
	int_t    *m_loc_vec;

	iam = grid->iam;
	n = A.size;
	m = n;

	cout << "OK1" << endl;
	cout << grid->nprow << "," << grid->npcol << "," << sizeof(int_t) << endl;

	m_loc_vec = (int_t*)malloc(grid->nprow * grid->npcol*sizeof(int_t));
	cout << "OK2" << endl;
	m_loc_fst = m / (grid->nprow * grid->npcol);
	for (int i = 0; i < grid->nprow * grid->npcol; i++) {
		if (i < grid->nprow * grid->npcol-1 ) {
			m_loc_vec[i] = m_loc_fst;
		}
		else { 
			m_loc_vec[i] = m - m_loc_fst*(grid->nprow * grid->npcol-1);
		} 
	}
	m_loc = m_loc_vec[iam];
	
	
	rowptr_loc = (int_t*)intMalloc_dist((m_loc+1)); 

	/* construct local row pointer. ASSUMING symmetric matrix */
	// Note that -1 cancels
	for (int i = 0; i < m_loc+1; i++)
		rowptr_loc[i] = A.colptr[iam*m_loc_fst+i] - A.colptr[iam*m_loc_fst]; 


	/* calculate nnz_loc on each processor */
	nnz_loc = rowptr_loc[m_loc]-rowptr_loc[0];
	colind_loc = (int_t*)intMalloc_dist(nnz_loc); 
	nzval_loc  = (doublecomplex*)doublecomplexMalloc_dist(nnz_loc);

	cout << "nnz_loc = " << nnz_loc << endl;

	// -1 is VERY IMPORTANT
	int disp = A.colptr[iam*m_loc_fst] - 1;
	for(int i = 0; i < nnz_loc; i++){
		colind_loc[i] = A.rowind[disp+i] - 1 ;
		nzval_loc[i].r  = A.nzval[disp+i].real();
		nzval_loc[i].i  = A.nzval[disp+i].imag();
	}


	fst_row = iam*m_loc_fst;
	zCreate_CompRowLoc_Matrix_dist(ANRloc, m, n, nnz_loc, m_loc, fst_row,
			nzval_loc, colind_loc, rowptr_loc,
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

		superlu_options_t superlu_options;
		SuperLUStat_t stat;
		SuperMatrix A;
		ScalePermstruct_t ScalePermstruct;
		LUstruct_t LUstruct;
		SOLVEstruct_t SOLVEstruct;
		gridinfo_t grid;
		int_t    m, n;
		int_t    nprow, npcol;
		int      iam, info, ldb, ldx, nrhs;
		doublecomplex *b = NULL;
		double *berr = NULL;
//		char     **cpp, c;
		FILE *fp;
//		extern int cpp_defs();

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
		
		// Test code
		if(1){
			Real timeSta, timeEnd;
			GetTime( timeSta );
			nrhs = 0;
			superlu_gridinit(MPI_COMM_WORLD, nprow, npcol, &grid);
			cout << nprow << "," << npcol << endl;
			cout << grid.nprow << "," << grid.npcol << endl;
			if( !(fp = fopen("H_LU.csc", "r")) ){
				throw std::logic_error( "H file not exist." );
			}
			ReadSparseMatrix( "H_LU.csc", pexsiData.HMat );
			ReadSparseMatrix( "S_LU.csc", pexsiData.SMat );
			GetTime( timeEnd );
			cout << "Time for reading H and S is " << timeEnd - timeSta << endl;
		}

		// Generate the distributed SuperMatrix
		if(1){
			Real timeSta, timeEnd;
			GetTime( timeSta );

			SparseMatrix<Complex>  AMat;
			AMat.size   = pexsiData.HMat.size;
			AMat.nnz    = pexsiData.HMat.nnz;
			AMat.colptr = PEXSI::IntNumVec( pexsiData.HMat.colptr.m(), false, pexsiData.HMat.colptr.Data() );
			AMat.rowind = PEXSI::IntNumVec( pexsiData.HMat.rowind.m(), false, pexsiData.HMat.rowind.Data() );
			AMat.nzval.Resize( pexsiData.HMat.nnz );
  
			for(Int i = 0; i < pexsiData.HMat.nnz; i++){
				AMat.nzval(i) = pexsiData.HMat.nzval(i) - Complex(0.0,1.0)*pexsiData.SMat.nzval(i);
			}

			SparseMatrixToSuperMatrixNRloc(&A, AMat, &grid);
			GetTime( timeEnd );
			cout << "Time for converting the matrix A is " << timeEnd - timeSta << endl;

			// Clear the memory
		  pexsiData.HMat.colptr.Resize(0);
			pexsiData.HMat.rowind.Resize(0);
			pexsiData.HMat.nzval.Resize(0);
		  pexsiData.SMat.colptr.Resize(0);
			pexsiData.SMat.rowind.Resize(0);
			pexsiData.SMat.nzval.Resize(0);
//		  AMat.colptr.Resize(0);
//			AMat.rowind.Resize(0);
			AMat.nzval.Resize(0);

		}

 
		// Factorization without solve
		if(1){
			set_default_options_dist(&superlu_options);

			superlu_options.Fact = DOFACT;
			superlu_options.RowPerm = NOROWPERM;
			superlu_options.IterRefine = NOREFINE;
			superlu_options.ParSymbFact       = NO;
			superlu_options.Equil = NO; 
			superlu_options.ReplaceTinyPivot = NO;
			superlu_options.ColPerm = MMD_AT_PLUS_A;
//			superlu_options.ColPerm = NATURAL;
			superlu_options.PrintStat         = YES;
			superlu_options.SolveInitialized  = NO;

			m = A.nrow;
			n = A.ncol;


			/* Initialize ScalePermstruct and LUstruct. */
			ScalePermstructInit(m, n, &ScalePermstruct);
			LUstructInit(m, n, &LUstruct);

			/* Initialize the statistics variables. */
			PStatInit(&stat);

			ldb = n;
			ldx = n;
			nrhs = 0;

			Real timeFactorSta, timeFactorEnd;
			GetTime( timeFactorSta );
			pzgssvx(&superlu_options, &A, &ScalePermstruct, b, ldb, nrhs, &grid,
					&LUstruct, &SOLVEstruct, berr, &stat, &info);
			GetTime( timeFactorEnd );

			cout << "Time for factorization is " << timeFactorEnd - timeFactorSta << endl;
			//#endif

			//		PStatPrint(&superlu_options, &stat, &grid);        /* Print the statistics. */
		
//			SUPERLU_FREE(b);
//			SUPERLU_FREE(xtrue);
//			SUPERLU_FREE(berr);
		}



		// *********************************************************************
		// Deallocate the storage
		// *********************************************************************
		PStatFree(&stat);
		Destroy_CompRowLoc_Matrix_dist(&A);
		ScalePermstructFree(&ScalePermstruct);
		Destroy_LU(n, &grid, &LUstruct);
		LUstructFree(&LUstruct);
		superlu_gridexit(&grid);

		statusOFS.close();
	}
	catch( std::exception& e )
	{
		std::cerr << " caught exception with message: "
	
	MPI_Finalize();

	return 0;
}
