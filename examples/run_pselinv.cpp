/// @file run_pselinv.cpp
/// @brief Test for the interface of SuperLU_DIST and SelInv.
/// @author Lin Lin
/// @date 2013-04-15
#include  "ppexsi.hpp"


#ifdef USE_TAU
  #include "TAU.h"
#elif defined (PROFILE) || defined(PMPI)
  #define TAU
  #include "timer.h"
#endif



using namespace PEXSI;
using namespace std;

void Usage(){
  std::cout 
		<< "Usage" << std::endl
		<< "run_pselinv -T [isText] -F [doFacto -E [doTriSolve] -Sinv [doSelInv]]  -H <Hfile> -S [Sfile] -colperm [colperm]" << std::endl;
}

int main(int argc, char **argv) 
{
	if( argc < 3 ) {
		Usage();
		return 0;
	}

#if defined(PROFILE) || defined(USE_TAU) || defined(PMPI)
 TAU_PROFILE_INIT(argc, argv);
#endif
	
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

    OptionsCreate(argc, argv, options);
    Int nprow = iround( std::sqrt( (double)mpisize) );
    Int npcol = mpisize / nprow;
    if( mpisize != nprow * npcol || nprow != npcol ){
      throw std::runtime_error( "nprow == npcol is assumed in this test routine." );
    }

    if( mpirank == 0 )
      cout << "nprow = " << nprow << ", npcol = " << npcol << endl;

    std::string Hfile, Sfile;
    int isCSC = true;
    if( options.find("-T") != options.end() ){ 
      isCSC= ! atoi(options["-T"].c_str());
    }

    int checkAccuracy = true;
    if( options.find("-E") != options.end() ){ 
      checkAccuracy= atoi(options["-E"].c_str());
    }

    int doFacto = true;
    if( options.find("-F") != options.end() ){ 
      doFacto= atoi(options["-F"].c_str());
    }

    int doSelInv = true;
    if( options.find("-Sinv") != options.end() ){ 
      doSelInv= atoi(options["-Sinv"].c_str());
    }

    int doSymbfact = true;
    if( options.find("-Symb") != options.end() ){ 
      doSymbfact= atoi(options["-Symb"].c_str());
    }


    doFacto = doFacto && doSymbfact;


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

    Int maxPipelineDepth = -1;
    if( options.find("-P") != options.end() ){ 
      maxPipelineDepth = atoi(options["-P"].c_str());
    }
    else{
      statusOFS << "-P option is not given. " 
        << "Do not limit SelInv pipelining depth." 
        << std::endl << std::endl;
    }

    Int maxDomains = -1;
    if( options.find("-D") != options.end() ){ 
      maxDomains = atoi(options["-D"].c_str());
    }
    else{
      statusOFS << "-D option is not given. " 
        << "Do not limit partitionning domain sizes." 
        << std::endl << std::endl;
    }



    std::string ColPerm;
    if( options.find("-colperm") != options.end() ){ 
      ColPerm = options["-colperm"];
    }
    else{
      statusOFS << "-colperm option is not given. " 
        << "Use MMD_AT_PLUS_A." 
        << std::endl << std::endl;
      ColPerm = "MMD_AT_PLUS_A";
    }

    // *********************************************************************
    // Read input matrix
    // *********************************************************************

    // Setup grid.
    SuperLUGrid g( MPI_COMM_WORLD, nprow, npcol );

    int      m, n;
    DistSparseMatrix<Complex>  AMat;

    DistSparseMatrix<Real> HMat;
    DistSparseMatrix<Real> SMat;
    Real timeSta, timeEnd;
    GetTime( timeSta );
    if(isCSC)
      ReadDistSparseMatrix( Hfile.c_str(), HMat, MPI_COMM_WORLD ); 
    else
      ReadDistSparseMatrixFormatted( Hfile.c_str(), HMat, MPI_COMM_WORLD ); 
    if( Sfile.empty() ){
      // Set the size to be zero.  This will tell PPEXSI.Solve to treat
      // the overlap matrix as an identity matrix implicitly.
      SMat.size = 0;  
    }
    else{
      if(isCSC)
        ReadDistSparseMatrix( Sfile.c_str(), SMat, MPI_COMM_WORLD ); 
      else
        ReadDistSparseMatrixFormatted( Sfile.c_str(), SMat, MPI_COMM_WORLD ); 
    }

    GetTime( timeEnd );
    if( mpirank == 0 ){
      cout << "Time for reading H and S is " << timeEnd - timeSta << endl;
      cout << "H.size = " << HMat.size << endl;
      cout << "H.nnz  = " << HMat.nnz  << endl;
    }

    // Get the diagonal indices for H and save it n diagIdxLocal_

    std::vector<Int>  diagIdxLocal;
    { 
      Int numColLocal      = HMat.colptrLocal.m() - 1;
      Int numColLocalFirst = HMat.size / mpisize;
      Int firstCol         = mpirank * numColLocalFirst;

      diagIdxLocal.clear();

      for( Int j = 0; j < numColLocal; j++ ){
        Int jcol = firstCol + j + 1;
        for( Int i = HMat.colptrLocal(j)-1; 
            i < HMat.colptrLocal(j+1)-1; i++ ){
          Int irow = HMat.rowindLocal(i);
          if( irow == jcol ){
            diagIdxLocal.push_back( i );
          }
        }
      } // for (j)
    }


    GetTime( timeSta );

    AMat.size          = HMat.size;
    AMat.nnz           = HMat.nnz;
    AMat.nnzLocal      = HMat.nnzLocal;
    AMat.colptrLocal   = HMat.colptrLocal;
    AMat.rowindLocal   = HMat.rowindLocal;
    AMat.nzvalLocal.Resize( HMat.nnzLocal );
    AMat.comm = MPI_COMM_WORLD;

    Complex *ptr0 = AMat.nzvalLocal.Data();
    Real *ptr1 = HMat.nzvalLocal.Data();
    Real *ptr2 = SMat.nzvalLocal.Data();
    Complex zshift = Complex(1.0, 0.0);

    if( SMat.size != 0 ){
      // S is not an identity matrix
      for( Int i = 0; i < HMat.nnzLocal; i++ ){
        AMat.nzvalLocal(i) = HMat.nzvalLocal(i) - zshift * SMat.nzvalLocal(i);
      }
    }
    else{
      // S is an identity matrix
      for( Int i = 0; i < HMat.nnzLocal; i++ ){
        AMat.nzvalLocal(i) = HMat.nzvalLocal(i);
      }

      for( Int i = 0; i < diagIdxLocal.size(); i++ ){
        AMat.nzvalLocal( diagIdxLocal[i] ) -= zshift;
      }
    } // if (SMat.size != 0 )

    GetTime( timeEnd );
    if( mpirank == 0 )
      cout << "Time for constructing the matrix A is " << timeEnd - timeSta << endl;


    // *********************************************************************
    // Symbolic factorization 
    // *********************************************************************

    GetTime( timeSta );
    SuperLUOptions luOpt;
    luOpt.ColPerm = ColPerm;
    luOpt.MaxPipelineDepth = maxPipelineDepth;
    luOpt.MaxDomains = maxDomains;
    SuperLUMatrix luMat( g, luOpt );
    luMat.DistSparseMatrixToSuperMatrixNRloc( AMat );
    GetTime( timeEnd );
    if( mpirank == 0 )
      cout << "Time for converting to SuperLU format is " << timeEnd - timeSta << endl;


    if(doSymbfact){
    GetTime( timeSta );
    luMat.SymbolicFactorize();
    luMat.DestroyAOnly();
    GetTime( timeEnd );

    if( mpirank == 0 )
      cout << "Time for performing the symbolic factorization is " << timeEnd - timeSta << endl;
    }

    // *********************************************************************
    // Numerical factorization only 
    // *********************************************************************

    if(doFacto){
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

      if( checkAccuracy ) {
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

      if( doSelInv )
      {
        Real timeTotalSelInvSta, timeTotalSelInvEnd;
        GetTime( timeTotalSelInvSta );

        Grid g1( MPI_COMM_WORLD, nprow, npcol );
        SuperNode super;

        GetTime( timeSta );
        luMat.SymbolicToSuperNode( super );
        PMatrix PMloc( &g1, &super, &luOpt );
        luMat.LUstructToPMatrix( PMloc );
        GetTime( timeEnd );
#ifdef SANITY_CHECK
        SuperNode superRef;
        luMat.SymbolicToSuperNode( superRef );
        PMatrix PMlocRef( &g1, &superRef, &luOpt  );
        luMat.LUstructToPMatrix( PMlocRef );
#endif
        if( mpirank == 0 )
          cout << "Time for converting LUstruct to PMatrix is " << timeEnd  - timeSta << endl;

        statusOFS << "perm: " << endl << super.perm << endl;
        statusOFS << "permInv: " << endl << super.permInv << endl;
        statusOFS << "superIdx:" << endl << super.superIdx << endl;
        statusOFS << "superPtr:" << endl << super.superPtr << endl; 


        // Preparation for the selected inversion
        GetTime( timeSta );
        PMloc.ConstructCommunicationPattern();
        GetTime( timeEnd );
#ifdef SANITY_CHECK
        PMlocRef.ConstructCommunicationPattern();
#endif
        if( mpirank == 0 )
          cout << "Time for constructing the communication pattern is " << timeEnd  - timeSta << endl;


        GetTime( timeSta );
        PMloc.PreSelInv();
        GetTime( timeEnd );
#ifdef SANITY_CHECK
        PMlocRef.PreSelInv();
#endif
        if( mpirank == 0 )
          cout << "Time for pre-selected inversion is " << timeEnd  - timeSta << endl;

        // Main subroutine for selected inversion
        GetTime( timeSta );
        PMloc.SelInv();
        GetTime( timeEnd );
        if( mpirank == 0 )
          cout << "Time for numerical selected inversion is " << timeEnd  - timeSta << endl;


        GetTime( timeTotalSelInvEnd );
        if( mpirank == 0 )
          cout << "Time for total selected inversion is " << timeTotalSelInvEnd  - timeTotalSelInvSta << endl;

#ifdef SANITY_CHECK
        GetTime( timeSta );
        PMlocRef.SelInvOriginal();
        GetTime( timeEnd );
        GetTime( timeTotalSelInvEnd );
        if( mpirank == 0 )
          cout << "Time for numerical selected inversion (original) is " << timeEnd  - timeSta << endl;
#endif


        NumVec<Scalar> diag;

        GetTime( timeSta );
        PMloc.GetDiagonal( diag );
        GetTime( timeEnd );

#ifdef SANITY_CHECK
        Real maxError = 0.0;
        PMlocRef.CompareDiagonal( PMloc ,maxError);

        if( mpirank == 0 )
          cout << "Max ||diag - diagRef||_2 = " << maxError << std::endl;
#endif

        if( mpirank == 0 )
          cout << "Time for getting the diagonal is " << timeEnd  - timeSta << endl;

        if( mpirank == 0 ){
          statusOFS << std::endl << "Diagonal of inverse in natural order: " << std::endl << diag << std::endl;
          ofstream ofs("diag");
          if( !ofs.good() ) 
            throw std::runtime_error("file cannot be opened.");
          serialize( diag, ofs, NO_MASK );
          ofs.close();
        }




        // Convert to DistSparseMatrix and get the diagonal
        GetTime( timeSta );
        DistSparseMatrix<Scalar> Ainv;
        PMloc.PMatrixToDistSparseMatrix( Ainv );
        GetTime( timeEnd );

        if( mpirank == 0 )
          cout << "Time for converting PMatrix to DistSparseMatrix is " << timeEnd  - timeSta << endl;

        NumVec<Scalar> diagDistSparse;
        GetTime( timeSta );
        GetDiagonal( Ainv, diagDistSparse );
        GetTime( timeEnd );
        if( mpirank == 0 )
          cout << "Time for getting the diagonal of DistSparseMatrix is " << timeEnd  - timeSta << endl;

        if( mpirank == 0 ){
          statusOFS << std::endl << "Diagonal of inverse from DistSparseMatrix format : " << std::endl << diagDistSparse << std::endl;
          Real diffNorm = 0.0;;
          for( Int i = 0; i < diag.m(); i++ ){
            diffNorm += pow( std::abs( diag(i) - diagDistSparse(i) ), 2.0 );
          }
          diffNorm = std::sqrt( diffNorm );
          statusOFS << std::endl << "||diag - diagDistSparse||_2 = " << diffNorm << std::endl;
        }

        // Convert to DistSparseMatrix in the 2nd format and get the diagonal
        GetTime( timeSta );
        DistSparseMatrix<Scalar> Ainv2;
        PMloc.PMatrixToDistSparseMatrix( AMat, Ainv2 );
        GetTime( timeEnd );

        if( mpirank == 0 )
          cout << "Time for converting PMatrix to DistSparseMatrix (2nd format) is " << timeEnd  - timeSta << endl;

        NumVec<Scalar> diagDistSparse2;
        GetTime( timeSta );
        GetDiagonal( Ainv2, diagDistSparse2 );
        GetTime( timeEnd );
        if( mpirank == 0 )
          cout << "Time for getting the diagonal of DistSparseMatrix is " << timeEnd  - timeSta << endl;

        if( mpirank == 0 ){
          statusOFS << std::endl << "Diagonal of inverse from the 2nd conversion into DistSparseMatrix format : " << std::endl << diagDistSparse2 << std::endl;
          Real diffNorm = 0.0;;
          for( Int i = 0; i < diag.m(); i++ ){
            diffNorm += pow( std::abs( diag(i) - diagDistSparse2(i) ), 2.0 );
          }
          diffNorm = std::sqrt( diffNorm );
          statusOFS << std::endl << "||diag - diagDistSparse2||_2 = " << diffNorm << std::endl;
        }

        Complex traceLocal = blas::Dot( AMat.nnzLocal, AMat.nzvalLocal.Data(), 1, 
            Ainv2.nzvalLocal.Data(), 1 );
        Complex trace = Z_ZERO;
        mpi::Allreduce( &traceLocal, &trace, 1, MPI_SUM, MPI_COMM_WORLD );

        if( mpirank == 0 )
          statusOFS << std::endl << "Tr[Ainv2 * AMat] = " << std::endl << trace << std::endl;


      }
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
