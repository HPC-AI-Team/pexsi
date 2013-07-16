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

    int doToDist = true;
    if( options.find("-ToDist") != options.end() ){ 
      doToDist= atoi(options["-ToDist"].c_str());
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

      for(int i=1; i<= doSelInv; ++i )
      {

        Real timeTotalSelInvSta, timeTotalSelInvEnd;






#ifdef USE_BCAST_UL
  NumVec<Scalar> diagBcast;
#endif

#ifdef SANITY_CHECK
  NumVec<Scalar> diagRef;
{
  GetTime( timeTotalSelInvSta );

  Grid g3( MPI_COMM_WORLD, nprow, npcol );
  SuperNode superRef;

  GetTime( timeSta );
  luMat.SymbolicToSuperNode( superRef );
  PMatrix PMlocRef( &g3, &superRef, &luOpt  );
  luMat.LUstructToPMatrix( PMlocRef );
  GetTime( timeEnd );

  if( mpirank == 0 )
    cout << "Time for converting LUstruct to PMatrix (Original) is " << timeEnd  - timeSta << endl;

  GetTime( timeSta );
  PMlocRef.ConstructCommunicationPatternOriginal();
  GetTime( timeEnd );
  if( mpirank == 0 )
    cout << "Time for constructing the communication pattern (Original) is " << timeEnd  - timeSta << endl;


  PMlocRef.PreSelInv();

  GetTime( timeSta );
  PMlocRef.SelInvOriginal();
  GetTime( timeEnd );
  GetTime( timeTotalSelInvEnd );
  if( mpirank == 0 )
    cout << "Time for numerical selected inversion (Original) is " << timeEnd  - timeSta << endl;


  GetTime( timeTotalSelInvEnd );
  if( mpirank == 0 )
    cout << "Time for total selected inversion (Original) is " << timeTotalSelInvEnd  - timeTotalSelInvSta << endl;




#ifdef USE_BCAST_UL
  {
    GetTime( timeTotalSelInvSta );

    Grid g2( MPI_COMM_WORLD, nprow, npcol );
    SuperNode superBcast;

    GetTime( timeSta );
    luMat.SymbolicToSuperNode( superBcast );
    PMatrix PMlocBcast( &g2, &superBcast, &luOpt  );
    luMat.LUstructToPMatrix( PMlocBcast );
    GetTime( timeEnd );

    if( mpirank == 0 )
      cout << "Time for converting LUstruct to PMatrix (Bcast) is " << timeEnd  - timeSta << endl;

    GetTime( timeSta );
    PMlocBcast.ConstructCommunicationPattern_Bcast();
    GetTime( timeEnd );
    if( mpirank == 0 )
      cout << "Time for constructing the communication pattern (Bcast) is " << timeEnd  - timeSta << endl;
    PMlocBcast.PreSelInv();

    GetTime( timeSta );
    PMlocBcast.SelInv_Bcast();
    GetTime( timeEnd );
    if( mpirank == 0 )
      cout << "Time for numerical selected inversion (Bcast) is " << timeEnd  - timeSta << endl;

    GetTime( timeTotalSelInvEnd );
    if( mpirank == 0 )
      cout << "Time for total selected inversion (Bcast) is " << timeTotalSelInvEnd  - timeTotalSelInvSta << endl;

        Real maxError = 0.0;
        PMlocRef.CompareOffDiagonal( PMlocBcast ,maxError);

        if( mpirank == 0 )
          cout << "Max ||OffdiagBcast - OffdiagRef||_2 // || OffdiagRef||_2 = " << maxError << std::endl;
        maxError = 0.0;
        PMlocRef.CompareDiagonal( PMlocBcast ,maxError);

        if( mpirank == 0 )
          cout << "Max ||diagBcast - diagRef||_2 // || diagRef||_2 = " << maxError << std::endl;


    GetTime( timeSta );
    PMlocBcast.GetDiagonal( diagBcast );
    GetTime( timeEnd );


    if( mpirank == 0 ){
      statusOFS << std::endl << "Diagonal (Bcast) of inverse in natural order: " << std::endl << diagBcast << std::endl;
      ofstream ofs("diag_bcast");
      if( !ofs.good() ) 
        throw std::runtime_error("file cannot be opened.");
      serialize( diagBcast, ofs, NO_MASK );
      ofs.close();
    }

//    PMlocBcast.DestructCommunicationPattern_Bcast( );

  }

#endif



  GetTime( timeSta );
  PMlocRef.GetDiagonal( diagRef );
  GetTime( timeEnd );


  if( mpirank == 0 ){
    statusOFS << std::endl << "Diagonal (original algorithm) of inverse in natural order: " << std::endl << diagRef << std::endl;
    ofstream ofs("diag_original");
    if( !ofs.good() ) 
      throw std::runtime_error("file cannot be opened.");
    serialize( diagRef, ofs, NO_MASK );
    ofs.close();
  }




}
#endif


#ifndef SANITY_CHECK
#ifdef USE_BCAST_UL
  {
    GetTime( timeTotalSelInvSta );

    Grid g2( MPI_COMM_WORLD, nprow, npcol );
    SuperNode superBcast;

    GetTime( timeSta );
    luMat.SymbolicToSuperNode( superBcast );
    PMatrix PMlocBcast( &g2, &superBcast, &luOpt  );
    luMat.LUstructToPMatrix( PMlocBcast );
    GetTime( timeEnd );

    if( mpirank == 0 )
      cout << "Time for converting LUstruct to PMatrix (Bcast) is " << timeEnd  - timeSta << endl;

    GetTime( timeSta );
    PMlocBcast.ConstructCommunicationPattern_Bcast();
    GetTime( timeEnd );
    if( mpirank == 0 )
      cout << "Time for constructing the communication pattern (Bcast) is " << timeEnd  - timeSta << endl;
    PMlocBcast.PreSelInv();

    GetTime( timeSta );
    PMlocBcast.SelInv_Bcast();
    GetTime( timeEnd );
    if( mpirank == 0 )
      cout << "Time for numerical selected inversion (Bcast) is " << timeEnd  - timeSta << endl;

    GetTime( timeTotalSelInvEnd );
    if( mpirank == 0 )
      cout << "Time for total selected inversion (Bcast) is " << timeTotalSelInvEnd  - timeTotalSelInvSta << endl;



    GetTime( timeSta );
    PMlocBcast.GetDiagonal( diagBcast );
    GetTime( timeEnd );


    if( mpirank == 0 ){
      statusOFS << std::endl << "Diagonal (Bcast) of inverse in natural order: " << std::endl << diagBcast << std::endl;
      ofstream ofs("diag_bcast");
      if( !ofs.good() ) 
        throw std::runtime_error("file cannot be opened.");
      serialize( diagBcast, ofs, NO_MASK );
      ofs.close();
    }


//    PMlocBcast.DestructCommunicationPattern_Bcast( );

  }
#endif
#endif







        GetTime( timeTotalSelInvSta );

        Grid g1( MPI_COMM_WORLD, nprow, npcol );
        SuperNode super;

        GetTime( timeSta );
        luMat.SymbolicToSuperNode( super );
        PMatrix PMloc( &g1, &super, &luOpt );
        luMat.LUstructToPMatrix( PMloc );
        GetTime( timeEnd );

        if( mpirank == 0 )
          cout << "Time for converting LUstruct to PMatrix is " << timeEnd  - timeSta << endl;

//        statusOFS << "perm: " << endl << super.perm << endl;
//        statusOFS << "permInv: " << endl << super.permInv << endl;
//        statusOFS << "superIdx:" << endl << super.superIdx << endl;
//        statusOFS << "superPtr:" << endl << super.superPtr << endl; 


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

        // Main subroutine for selected inversion
        GetTime( timeSta );
        PMloc.SelInv();
        GetTime( timeEnd );
        if( mpirank == 0 )
          cout << "Time for numerical selected inversion is " << timeEnd  - timeSta << endl;


        GetTime( timeTotalSelInvEnd );
        if( mpirank == 0 )
          cout << "Time for total selected inversion is " << timeTotalSelInvEnd  - timeTotalSelInvSta << endl;




        NumVec<Scalar> diag;

        GetTime( timeSta );
        PMloc.GetDiagonal( diag );
        GetTime( timeEnd );


        if( mpirank == 0 )
          cout << "Time for getting the diagonal is " << timeEnd  - timeSta << endl;



#ifdef SANITY_CHECK
{
  Real maxError = 0.0;
  if( mpirank == 0 ){
    for( Int i = 0; i < diag.m(); i++ ){
      std::stringstream msg;
      Real error = abs(diag(i)-diagRef(i))/abs(diagRef(i));
      msg<< "Row "<<i<<" is wrong : "<< diag(i) << " vs "<<diagRef(i)<< " error is "<< error <<std::endl; 
      if( error > maxError){ maxError = error;}
      if(error > SANITY_PRECISION){
        statusOFS<<msg;
      }
    }
      cout << "Max ||diag - diagRef||_2 // || diagRef||_2 = " << maxError << std::endl;
  }

#ifdef USE_BCAST_UL
  maxError = 0.0;

  if( mpirank == 0 ){
    for( Int i = 0; i < diagBcast.m(); i++ ){
      std::stringstream msg;
      Real error = abs(diagBcast(i)-diagRef(i))/abs(diagRef(i));
      msg<< "Row "<<i<<" is wrong : "<< diagBcast(i) << " vs "<<diagRef(i)<< " error is "<< error <<std::endl; 
      if( error > maxError){ maxError = error;}
      if(error > SANITY_PRECISION){
        statusOFS<<msg;
      }
    }
      cout << "Max ||diagBcast - diagRef||_2 // || diagRef||_2 = " << maxError << std::endl;
  }
#endif

}
#endif


        if( mpirank == 0 ){
          statusOFS << std::endl << "Diagonal of inverse in natural order: " << std::endl << diag << std::endl;
          ofstream ofs("diag");
          if( !ofs.good() ) 
            throw std::runtime_error("file cannot be opened.");
          serialize( diag, ofs, NO_MASK );
          ofs.close();
        }


        if(doToDist){
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
          PMloc.PMatrixToDistSparseMatrix2( AMat, Ainv2 );
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

          if( mpirank == 0 ){

            cout << "H.size = " << std::scientific << HMat.size << endl;
            cout << std::endl << "Tr[Ainv2 * AMat] = " << std::scientific<< trace << std::endl;
            statusOFS << std::endl << "Tr[Ainv2 * AMat] = " << std::endl << trace << std::endl;


          }
        }

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
