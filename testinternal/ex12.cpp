#include  "environment_impl.hpp"
#include  "sparse_matrix_impl.hpp"
#include  "numvec_impl.hpp"
#include  "utility.hpp"
#include  "pluselinv.hpp"

using namespace PEXSI;
using namespace std;

void Usage(){
  std::cout 
		<< "Test for the factorization and the selected inversion with complex arithmetic." <<std::endl
		<< "The code is divided into different phases: " <<std:: endl
		<< "1. symbolic factorization " << std::endl
		<< "2. data distribution CompRow->SuperLU " << std::endl
		<< "3. numerical factorization" << std::endl
		<< "4. data distribution SuperLU->PMatrix" << std::endl
		<< "5. selected inversion" << std::endl 
		<< std::endl
		<< "ex12 -r [nprow] -c [npcol]" << std::endl;
}

extern "C"{
void
pzsymbfact(superlu_options_t *options, SuperMatrix *A, 
					 ScalePermstruct_t *ScalePermstruct, gridinfo_t *grid,
					 LUstruct_t *LUstruct, SuperLUStat_t *stat, int *info);
}




int main(int argc, char **argv) 
{
	MPI_Init(&argc, &argv);
	int mpirank, mpisize;
	MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
	MPI_Comm_size( MPI_COMM_WORLD, &mpisize );

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
		ScalePermstruct_t ScalePermstruct;
		LUstruct_t LUstruct;
		SOLVEstruct_t SOLVEstruct;
		gridinfo_t grid;
		int      nprow, npcol;
		int      m, n;
		DistSparseMatrix<Complex>  AMat;
		std::string Hfile, Sfile;

		// For selected inversion
		PMatrix PMloc;
		std::vector<std::vector<Int> > localEtree; 

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


		// *********************************************************************
		// Symbolic factorization only
		// *********************************************************************
 
		set_default_options_dist(&superlu_options);
		superlu_options.Fact              = DOFACT;
		superlu_options.RowPerm           = NOROWPERM;
//		superlu_options.IterRefine        = SLU_SINGLE;
		superlu_options.IterRefine        = NOREFINE;
		superlu_options.ParSymbFact       = NO;
		superlu_options.Equil             = NO; 
		superlu_options.ReplaceTinyPivot  = NO;
//		superlu_options.ColPerm           = METIS_AT_PLUS_A;
		superlu_options.ColPerm           = MMD_AT_PLUS_A;
//		superlu_options.ColPerm = NATURAL;
		superlu_options.PrintStat         = YES;
		superlu_options.SolveInitialized  = NO;

		m = AMat.size;
		n = AMat.size;
		int     info;
		Real timeFactorSta, timeFactorEnd;

		ScalePermstructInit(m, n, &ScalePermstruct);
		LUstructInit(m, n, &LUstruct);

		GetTime( timeSta );
		DistSparseMatrixToSuperMatrixNRloc(&A, AMat, &grid);
		GetTime( timeEnd );
		if( mpirank == 0 )
			cout << "Time for converting to SuperLU format is " << timeEnd - timeSta << endl;

		GetTime( timeSta );
		PStatInit(&stat);
		pzsymbfact(&superlu_options, &A, &ScalePermstruct, &grid, 
				&LUstruct, &stat, &info);
		PStatFree(&stat);
		GetTime( timeEnd );

		if( mpirank == 0 )
			cout << "Time for performing the symbolic factorization is " << timeEnd - timeSta << endl;


		// *********************************************************************
		// Numerical factorization only 
		// *********************************************************************
		
		// Important: the distribution in pzsymbfact is going to mess up the
		// A matrix.  Recompute the matrix A here.
		GetTime( timeSta ); 
		Destroy_CompRowLoc_Matrix_dist(&A);
		DistSparseMatrixToSuperMatrixNRloc(&A, AMat, &grid);
		GetTime( timeEnd );
		if( mpirank == 0 )
			cout << "Time for converting to SuperLU format is " << timeEnd - timeSta << endl;

		superlu_options.Fact              = SamePattern_SameRowPerm;

		// Distribution phase
		if(1)
		{
			Int* perm_c = ScalePermstruct.perm_c;
			NRformat_loc *Astore = (NRformat_loc *) A.Store;
			Int* colind = Astore->colind;
			Int  nnzLocal = Astore->nnz_loc;
			if(mpirank == 0) {
				cout << "nnzLocal = " << nnzLocal << endl;
				cout << "colind[end] = " << colind[nnzLocal-1] << endl;
			}
  	  // Apply column permutation to the original distributed A
			for(Int j = 0; j < nnzLocal; j++)
				colind[j] = perm_c[colind[j]];
	    /* Distribute Pc*Pr*diag(R)*A*diag(C)*Pc' into L and U storage. 
	       NOTE: the row permutation Pc*Pr is applied internally in the
  	       distribution routine. */
			GetTime( timeSta );
	    float dist_mem_use = pzdistribute(SamePattern_SameRowPerm, n, &A, &ScalePermstruct,
					NULL, &LUstruct, &grid);
			GetTime( timeEnd );
			if(mpirank == 0)
				cout << "Time for distribution is " << timeEnd - timeSta << " sec" << endl;
		}

		GetTime( timeFactorSta );
		PStatInit(&stat);
		// Estimate the one-norm
		string norm = "1";
		double anorm = pzlangs((char*)norm.c_str(), &A, &grid);
		if(mpirank == 0) cout << "anorm = " << anorm << endl;
		
		// Numerical factorization phase
		pzgstrf(&superlu_options, m, n, anorm, &LUstruct, &grid, &stat, &info); 
		PStatFree(&stat); 
		GetTime( timeFactorEnd );

		if( mpirank == 0 )
			cout << "Time for factorization is " << timeFactorEnd - timeFactorSta << " sec" << endl; 



		// *********************************************************************
		// Test the accuracy of factorization by solve
		// *********************************************************************

		if(0){
			PStatInit(&stat);
			superlu_options.Fact              = FACTORED;

			int nrhs = 1;
			PEXSI::CpxNumVec  xTrueGlobal(n), bGlobal(n);
			{
				SuperMatrix GA, A1;
				int needValue = 1;
				DistSparseMatrixToSuperMatrixNRloc(&A1, AMat, &grid);
				pzCompRow_loc_to_CompCol_global(needValue, &A1, &grid, &GA);
				UniformRandom( xTrueGlobal );
				//				zGenXtrue_dist(n, 1, (doublecomplex*)xTrueGlobal.Data(), n);
				char trans[1]; *trans='N';

				zFillRHS_dist(trans, 1, (doublecomplex*)xTrueGlobal.Data(), n, &GA,
						(doublecomplex*)bGlobal.Data(), n);
				Destroy_CompCol_Matrix_dist(&GA);
				Destroy_CompRowLoc_Matrix_dist(&A1);

				//				if(mpirank == 0){
				//					cerr << " x = " << xTrueGlobal << endl;;
				//					cerr << " b = " << bGlobal << endl;
				//				}
			}

			// Get the local solution
			doublecomplex *bLocal, *xTrueLocal;
			Int numRowLocal = AMat.colptrLocal.m() - 1;
			Int numRowLocalFirst = AMat.size / mpisize;
			Int firstRow = mpirank * numRowLocalFirst;
			bLocal = doublecomplexMalloc_dist( numRowLocal ); 
			cout << "Proc " << mpirank << " outputs numRowLocal = " <<
				numRowLocal << " firstRow = " << firstRow << endl;
			std::copy( bGlobal.Data()+firstRow, bGlobal.Data()+firstRow+numRowLocal,
					(Complex*)bLocal );
			xTrueLocal = doublecomplexMalloc_dist( numRowLocal ); 
			std::copy( xTrueGlobal.Data()+firstRow, xTrueGlobal.Data()+firstRow+numRowLocal,
					(Complex*)xTrueLocal );
//			cout << "Energy(bLocal) = " << Energy(PEXSI::CpxNumVec(numRowLocal,false,(Complex*)bLocal)) << endl;

			// Call the linear equation solver. 
			double *berr=doubleMalloc_dist(nrhs);

			GetTime( timeFactorSta );
			pzgssvx(&superlu_options, &A, &ScalePermstruct, bLocal, numRowLocal, nrhs, &grid,
					&LUstruct, &SOLVEstruct, berr, &stat, &info);
			GetTime( timeFactorEnd );

//			cout << "Energy(bLocal) = " << Energy(PEXSI::CpxNumVec(numRowLocal,false,(Complex*)bLocal)) << endl;
//			cout << "Energy(xLocal) = " << Energy(PEXSI::CpxNumVec(numRowLocal,false,(Complex*)xTrueLocal)) << endl;


			pzinf_norm_error(mpirank, ((NRformat_loc *)A.Store)->m_loc, nrhs, 
					bLocal, numRowLocal, xTrueLocal, numRowLocal, &grid);

			PStatPrint(&superlu_options, &stat, &grid);        /* Print the statistics. */
			PStatFree(&stat);
			SUPERLU_FREE( bLocal );
			SUPERLU_FREE( xTrueLocal );
			SUPERLU_FREE( berr );
		}


		// *********************************************************************
		// Selected inversion
		// *********************************************************************
		
		{
			GetTime(timeSta);
			SuperLU2SelInv(n, &LUstruct, &grid, PMloc);
			GetTime(timeEnd);
			if( mpirank == 0 )
				cout << "Time for converting the SuperLU to SelInv format is " 
					<< timeEnd - timeSta << " sec" << endl; 

			GetTime(timeSta);
			ConstructLocalEtree(n, &grid, PMloc, localEtree);
			GetTime(timeEnd);
			if( mpirank == 0 )
				cout << "Time for converting the SelInv elimination tree is " 
					<< timeEnd - timeSta << " sec" << endl; 
			
			GetTime(timeSta);
			PLUSelInv(&grid, PMloc, localEtree);
			GetTime(timeEnd);
			if( mpirank == 0 )
				cout << "Time for the selected inversion is " 
					<< timeEnd - timeSta << " sec" << endl; 
		}

		// *********************************************************************
		// Deallocate the storage
		// *********************************************************************
		Destroy_CompRowLoc_Matrix_dist(&A);
		Destroy_LU(n, &grid, &LUstruct);
		ScalePermstructFree(&ScalePermstruct);
		LUstructFree(&LUstruct);
		superlu_gridexit(&grid);

		statusOFS.close();
	}
	catch( std::exception& e )
	{
		std::cerr << "Processor " << mpirank << " caught exception with message: "
	
	MPI_Finalize();

	return 0;
}
