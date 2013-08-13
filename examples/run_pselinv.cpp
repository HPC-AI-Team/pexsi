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
  std::cout << "Usage" << std::endl << "run_pselinv -T [isText] -F [doFacto -E [doTriSolve] -Sinv [doSelInv]]  -H <Hfile> -S [Sfile] -colperm [colperm] -npsymbfact [npsymbfact] -P [maxpipelinedepth] -SinvBcast [doSelInvBcast] -SinvPipeline [doSelInvPipeline] -SinvOriginal [doSelInvOriginal] -Shift [imaginary shift] " << std::endl;
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
    MPI_Comm world_comm;

    // *********************************************************************
    // Input parameter
    // *********************************************************************
    std::map<std::string,std::string> options;

    OptionsCreate(argc, argv, options);

    //set the processor grid size at the highest possible square root.
    Int nprow = iround( std::sqrt( (double)mpisize) );
    Int npcol = mpisize / nprow;
    nprow = npcol;

    if( options.find("-r") != options.end() ){
      if( options.find("-c") != options.end() ){
        nprow= atoi(options["-r"].c_str());
        npcol= atoi(options["-c"].c_str());
        if(nprow*npcol > mpisize){
          throw std::runtime_error("The number of used processors cannot be higher than the total number of available processors." );
        } 
      }
      else{
        throw std::runtime_error( "When using -r option, -c also needs to be provided." );
      }
    }
    else if( options.find("-c") != options.end() ){
      if( options.find("-r") != options.end() ){
        nprow= atoi(options["-r"].c_str());
        npcol= atoi(options["-c"].c_str());
        if(nprow*npcol > mpisize){
          throw std::runtime_error("The number of used processors cannot be higher than the total number of available processors." );
        } 
      }
      else{
        throw std::runtime_error( "When using -c option, -r also needs to be provided." );
      }
    }

    //Create a communicator with npcol*nprow processors
    MPI_Comm_split(MPI_COMM_WORLD, mpirank<nprow*npcol, mpirank, &world_comm);

    if (mpirank<nprow*npcol){

      MPI_Comm_rank(world_comm, &mpirank );
      MPI_Comm_size(world_comm, &mpisize );



      stringstream  ss;
      ss << "logTest" << mpirank;
      statusOFS.open( ss.str().c_str() );

      //if( mpisize != nprow * npcol || nprow != npcol ){
      //  throw std::runtime_error( "nprow == npcol is assumed in this test routine." );
      //}

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


      Int numProcSymbFact;
      if( options.find("-npsymbfact") != options.end() ){ 
        numProcSymbFact = atoi( options["-npsymbfact"].c_str() );
      }
      else{
        statusOFS << "-npsymbfact option is not given. " 
          << "Use default value (maximum number of procs)." 
          << std::endl << std::endl;
        numProcSymbFact = 0;
      }







      Int doSinv_Original = 0;
      if( options.find("-SinvOriginal") != options.end() ){ 
        doSinv_Original = atoi(options["-SinvOriginal"].c_str());
      }
      Int doSinv_Bcast = 0;
      if( options.find("-SinvBcast") != options.end() ){ 
        doSinv_Bcast = atoi(options["-SinvBcast"].c_str());
      }

      Int doSinvPipeline = 1;
      if( options.find("-SinvPipeline") != options.end() ){ 
        doSinvPipeline = atoi(options["-SinvPipeline"].c_str());
      }

      Real cshift = 0;
      if( options.find("-Shift") != options.end() ){ 
        cshift = atof(options["-Shift"].c_str());
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
      SuperLUGrid g( world_comm, nprow, npcol );

      int      m, n;
      DistSparseMatrix<Complex>  AMat;

      DistSparseMatrix<Real> HMat;
      DistSparseMatrix<Real> SMat;
      Real timeSta, timeEnd;
      GetTime( timeSta );
      if(isCSC)
        ParaReadDistSparseMatrix( Hfile.c_str(), HMat, world_comm ); 
      else{
        ReadDistSparseMatrixFormatted( Hfile.c_str(), HMat, world_comm ); 
        ParaWriteDistSparseMatrix( "H.csc", HMat, world_comm ); 
      }

      if( Sfile.empty() ){
        // Set the size to be zero.  This will tell PPEXSI.Solve to treat
        // the overlap matrix as an identity matrix implicitly.
        SMat.size = 0;  
      }
      else{
        if(isCSC)
          ParaReadDistSparseMatrix( Sfile.c_str(), SMat, world_comm ); 
        else
          ReadDistSparseMatrixFormatted( Sfile.c_str(), SMat, world_comm ); 
        ParaWriteDistSparseMatrix( "S.csc", SMat, world_comm ); 
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
      AMat.comm = world_comm;

      Complex *ptr0 = AMat.nzvalLocal.Data();
      Real *ptr1 = HMat.nzvalLocal.Data();
      Real *ptr2 = SMat.nzvalLocal.Data();
      Complex zshift = Complex(1.0, cshift);

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
      luOpt.maxPipelineDepth = maxPipelineDepth;
      luOpt.numProcSymbFact = numProcSymbFact;
      SuperLUMatrix * luMatPtr = new SuperLUMatrix( g, luOpt );
      SuperLUMatrix & luMat = *luMatPtr;
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

          A1.DistributeGlobalMultiVector( bGlobal,     bLocal );
          luMat.SolveDistMultiVector( bLocal, berr );
          A1.DistributeGlobalMultiVector( bGlobal,     bLocal );
          luMat.SolveDistMultiVector( bLocal, berr );
          A1.DistributeGlobalMultiVector( bGlobal,     bLocal );
          luMat.SolveDistMultiVector( bLocal, berr );
          //          luMat.CheckErrorDistMultiVector( bLocal, xTrueLocal );
        }

        // *********************************************************************
        // Selected inversion
        // *********************************************************************

        for(int i=1; i<= doSelInv; ++i )
        {

          Real timeTotalSelInvSta, timeTotalSelInvEnd;






          NumVec<Scalar> diagBcast;
          PMatrix * PMlocBcastPtr;
          SuperNode * superBcastPtr;
          Grid * g2Ptr;

          if(doSinv_Bcast)
          {
            GetTime( timeTotalSelInvSta );

            g2Ptr = new Grid( world_comm, nprow, npcol );
            Grid &g2 = *g2Ptr;

            superBcastPtr = new SuperNode();
            SuperNode & superBcast = *superBcastPtr;

            GetTime( timeSta );
            luMat.SymbolicToSuperNode( superBcast );

            PMlocBcastPtr = new PMatrix( &g2, &superBcast, &luOpt  );
            PMatrix & PMlocBcast = *PMlocBcastPtr;

            luMat.LUstructToPMatrix( PMlocBcast );
            GetTime( timeEnd );

            if( mpirank == 0 )
              cout << "Time for converting LUstruct to PMatrix (Bcast) is " << timeEnd  - timeSta << endl;

            GetTime( timeSta );
            PMlocBcast.ConstructCommunicationPattern_Bcast();
            GetTime( timeEnd );
            if( mpirank == 0 )
              cout << "Time for constructing the communication pattern (Bcast) is " << timeEnd  - timeSta << endl;

            GetTime( timeSta );
            PMlocBcast.PreSelInv();
            GetTime( timeEnd );
            if( mpirank == 0 )
              cout << "Time for pre selected inversion (Bcast) is " << timeEnd  - timeSta << endl;

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


            PMlocBcast.DestructCommunicationPattern_Bcast( );

          }


          NumVec<Scalar> diagRef;
          PMatrix * PMlocRefPtr;
          SuperNode * superRefPtr;
          Grid * g3Ptr;

          if(doSinv_Original)
          {
            GetTime( timeTotalSelInvSta );

            g3Ptr = new Grid( world_comm, nprow, npcol );
            Grid &g3 = *g3Ptr;

            superRefPtr = new SuperNode();
            SuperNode & superRef = *superRefPtr;

            GetTime( timeSta );
            luMat.SymbolicToSuperNode( superRef );

            PMlocRefPtr = new PMatrix( &g3, &superRef, &luOpt  );
            PMatrix & PMlocRef = *PMlocRefPtr;

            luMat.LUstructToPMatrix( PMlocRef );
            GetTime( timeEnd );

            if( mpirank == 0 )
              cout << "Time for converting LUstruct to PMatrix (_Original) is " << timeEnd  - timeSta << endl;

            GetTime( timeSta );
            PMlocRef.ConstructCommunicationPattern();
            GetTime( timeEnd );
            if( mpirank == 0 )
              cout << "Time for constructing the communication pattern (_Original) is " << timeEnd  - timeSta << endl;

            GetTime( timeSta );
            PMlocRef.PreSelInv();
            GetTime( timeEnd );
            if( mpirank == 0 )
              cout << "Time for pre selected inversion (_Original) is " << timeEnd  - timeSta << endl;

            GetTime( timeSta );
            PMlocRef.SelInv();
            GetTime( timeEnd );
            GetTime( timeTotalSelInvEnd );
            if( mpirank == 0 )
              cout << "Time for numerical selected inversion (_Original) is " << timeEnd  - timeSta << endl;


            GetTime( timeTotalSelInvEnd );
            if( mpirank == 0 )
              cout << "Time for total selected inversion (_Original) is " << timeTotalSelInvEnd  - timeTotalSelInvSta << endl;





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

          NumVec<Scalar> diag;
          PMatrix * PMlocPtr;
          SuperNode * superPtr;
          Grid * g1Ptr;

          if(doSinvPipeline)
          {
            GetTime( timeTotalSelInvSta );

            g1Ptr = new Grid( world_comm, nprow, npcol );
            Grid &g1 = *g1Ptr;

            superPtr = new SuperNode();
            SuperNode & super = *superPtr;

            GetTime( timeSta );
            luMat.SymbolicToSuperNode( super );

            PMlocPtr = new PMatrix( &g1, &super, &luOpt  );
            PMatrix & PMloc = *PMlocPtr;

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
            PMloc.ConstructCommunicationPattern_Pipeline();
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
            PMloc.SelInv_Pipeline();
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


            if( mpirank == 0 ){
              statusOFS << std::endl << "Diagonal of inverse in natural order: " << std::endl << diag << std::endl;
              ofstream ofs("diag");
              if( !ofs.good() ) 
                throw std::runtime_error("file cannot be opened.");
              serialize( diag, ofs, NO_MASK );
              ofs.close();
            }



          }


#ifdef SANITY_CHECK
          if(doSinv_Original && !checkAccuracy)
          {

            //sanity check for SelInvPipeline
            if(doSinvPipeline)
            {
              SelInvErrors errorsOffDiag;
              SelInvErrors errorsDiag;

              SuperNode & super = *superPtr;
              PMatrix & PMlocRef = *PMlocRefPtr;
              PMatrix & PMloc = *PMlocPtr;
              PMloc.CompareOffDiagonal( PMlocRef ,errorsOffDiag);

              if( mpirank == 0 ){
                cout << errorsOffDiag << std::endl;
              }

              PMloc.CompareDiagonal( PMlocRef , errorsDiag);
              if( mpirank == 0 ){
                cout << errorsDiag << std::endl;
              } 



              //SelInvError & error = (errorsDiag.MaxRelError.Value>errorsOffDiag.MaxRelError.Value)?errorsDiag.MaxRelError:errorsOffDiag.MaxRelError;
              SelInvError & error = (errorsDiag.MaxAbsError.Value>errorsOffDiag.MaxAbsError.Value)?errorsDiag.MaxAbsError:errorsOffDiag.MaxAbsError;

              if(error.Value>0)
              {
                NumVec<Scalar> col,col_Original;

                PMloc.GetColumn( error.j, col );
                PMlocRef.GetColumn( error.j, col_Original );
                NumVec<Scalar> colBcast;
                if(doSinv_Bcast){
                  PMatrix & PMlocBcast = *PMlocBcastPtr;
                  PMlocBcast.GetColumn( error.j, colBcast );
                }

                const IntNumVec& permInv = super.permInv;

                //get the indices in NATURAL ordering
                Int gi = permInv(error.i);
                Int gj = permInv(error.j);


                if( mpirank == 0 ){
                  cout<<"Permuted indices are ("<<error.i<<","<<error.j<<")"<<std::endl;
                  cout<<"Unpermuted indices are ("<<gi<<","<<gj<<")"<<std::endl;
                }
                SuperLUMatrix * AExplicitPtr = new SuperLUMatrix( g );
                SuperLUMatrix & AExplicit = *AExplicitPtr;

                SuperLUMatrix * GAPtr = new SuperLUMatrix( g );
                SuperLUMatrix & GA = *GAPtr;
                AExplicit.DistSparseMatrixToSuperMatrixNRloc( AMat );
                AExplicit.ConvertNRlocToNC( GA );


                int n = AExplicit.n();
                IntNumVec cols(1);
                cols(0)=gj;


                NumMat<Scalar>  * xTrueGlobal,* bGlobal;
                xTrueGlobal= new NumMat<Scalar>(n,1);
                bGlobal= new NumMat<Scalar>(n,1);
                CpxNumMat * xTrueLocal, *bLocal;
                xTrueLocal = new NumMat<Scalar>();
                bLocal = new NumMat<Scalar>();
                DblNumVec *  berr = new DblNumVec();



                IdentityCol(cols, *bGlobal );
                AExplicit.DistributeGlobalMultiVector( *bGlobal,     *bLocal );
                luMat.SolveDistMultiVector( *bLocal, *berr );


                Int displs[g1Ptr->mpisize];
                Int rcounts[g1Ptr->mpisize];
                Int avgSize = n/g1Ptr->mpisize;

                for(Int i=0;i<g1Ptr->mpisize-1;++i){
                  rcounts[i]=avgSize*sizeof(Scalar);
                  displs[i]=i*avgSize*sizeof(Scalar); 
                }
                rcounts[g1Ptr->mpisize-1]=(n-(g1Ptr->mpisize-1)*avgSize)*sizeof(Scalar);
                displs[g1Ptr->mpisize-1]=(g1Ptr->mpisize-1)*avgSize*sizeof(Scalar);


                MPI_Gatherv(bLocal->Data(), bLocal->m()*sizeof(Scalar), MPI_BYTE, xTrueGlobal->Data(),  rcounts,displs, MPI_BYTE, 0, world_comm);

                if( mpirank == 0 ){

                  Real absErrorOriginal, absErrorBcast, absErrorPipeline;
                  Real relErrorOriginal, relErrorBcast, relErrorPipeline;

                  absErrorPipeline = abs(col(gi) - (*xTrueGlobal)(gi,0));
                  relErrorPipeline = abs(absErrorPipeline/ (*xTrueGlobal)(gi,0));

                  absErrorOriginal = abs(col_Original(gi) - (*xTrueGlobal)(gi,0));
                  relErrorOriginal = abs(absErrorOriginal/(*xTrueGlobal)(gi,0));


                  if(doSinv_Bcast){
                    absErrorBcast = abs(colBcast(gi) - (*xTrueGlobal)(gi,0));
                    relErrorBcast = abs(absErrorBcast/(*xTrueGlobal)(gi,0));
                    cout<<std::endl<< "\t\tSolve \t SelInv_Original \t SelInvPipeline \t SelInv_Bcast"<<std::endl;
                    cout<< "Value \t\t"<<(*xTrueGlobal)(gi,0)<<"\t"<< col_Original(gi)<<"\t"<< col(gi)<< "\t" << colBcast(gi)<<std::endl;
                    cout<< "Rel error \t\t"<< "NA" <<"\t"<< relErrorOriginal<<"\t"<< relErrorPipeline<< "\t" << relErrorBcast<<std::endl;
                    cout<< "Abs error \t\t"<< "NA" <<"\t"<< absErrorOriginal<<"\t"<< absErrorPipeline<< "\t" << absErrorBcast<<std::endl<<std::endl;
                  }
                  else{

                    cout<<std::endl<< "\t\tSolve \t SelInv_Original \t SelInvPipeline"<<std::endl;
                    cout<< "Value \t\t"<<(*xTrueGlobal)(gi,0)<<"\t"<< col_Original(gi)<<"\t"<< col(gi)<< std::endl;
                    cout<< "Rel error \t\t"<< "NA" <<"\t"<< relErrorOriginal<<"\t"<< relErrorPipeline<< std::endl;
                    cout<< "Abs error \t\t"<< "NA" <<"\t"<< absErrorOriginal<<"\t"<< absErrorPipeline<< std::endl<<std::endl;
                  }
                }     

                delete bGlobal;
                delete xTrueGlobal;
                delete xTrueLocal;
                delete bLocal;
                delete berr;



                delete GAPtr;
                delete AExplicitPtr;
              }
            }



            if(doSinv_Bcast){
              SelInvErrors errorsOffDiag;
              SelInvErrors errorsDiag;

              PMatrix & PMlocRef = *PMlocRefPtr;
              PMatrix & PMlocBcast = *PMlocBcastPtr;
              PMlocBcast.CompareOffDiagonal( PMlocRef ,errorsOffDiag);

              if( mpirank == 0 ){
                cout << errorsOffDiag << std::endl;
              }

              PMlocBcast.CompareDiagonal( PMlocRef , errorsDiag);
              if( mpirank == 0 ){
                cout << errorsDiag << std::endl;
              } 


              //SelInvError & error = (errorsDiag.MaxRelError.Value>errorsOffDiag.MaxRelError.Value)?errorsDiag.MaxRelError:errorsOffDiag.MaxRelError;
              SelInvError & error = (errorsDiag.MaxAbsError.Value>errorsOffDiag.MaxAbsError.Value)?errorsDiag.MaxAbsError:errorsOffDiag.MaxAbsError;

              if(error.Value>0)
              {
                NumVec<Scalar> col,col_Original,colBcast;
                PMlocBcast.GetColumn( error.j, colBcast );
                PMlocRef.GetColumn( error.j, col_Original );

                if(doSinvPipeline){
                  PMatrix & PMloc = *PMlocPtr;
                  PMloc.GetColumn( error.j, col );
                }

                const IntNumVec& permInv_Bcast = superBcastPtr->permInv;

                //get the indices in NATURAL ordering
                Int gi = permInv_Bcast(error.i);
                Int gj = permInv_Bcast(error.j);

                if( mpirank == 0 ){
                  cout<<"Permuted indices are ("<<error.i<<","<<error.j<<")"<<std::endl;
                  cout<<"Unpermuted indices are ("<<gi<<","<<gj<<")"<<std::endl;
                }

                SuperLUMatrix * AExplicitPtr = new SuperLUMatrix( g );
                SuperLUMatrix & AExplicit = *AExplicitPtr;

                SuperLUMatrix * GAPtr = new SuperLUMatrix( g );
                SuperLUMatrix & GA = *GAPtr;
                AExplicit.DistSparseMatrixToSuperMatrixNRloc( AMat );
                AExplicit.ConvertNRlocToNC( GA );


                int n = AExplicit.n();
                IntNumVec cols(1);
                cols(0)=gj;

                NumMat<Scalar>  * xTrueGlobal,* bGlobal;
                xTrueGlobal= new NumMat<Scalar>(n,1);
                bGlobal= new NumMat<Scalar>(n,1);
                CpxNumMat * xTrueLocal, *bLocal;
                xTrueLocal = new NumMat<Scalar>();
                bLocal = new NumMat<Scalar>();
                DblNumVec *  berr = new DblNumVec();

                IdentityCol(cols, *bGlobal );
                AExplicit.DistributeGlobalMultiVector( *bGlobal,     *bLocal );
                luMat.SolveDistMultiVector( *bLocal, *berr );

                Int displs[g2Ptr->mpisize];
                Int rcounts[g2Ptr->mpisize];
                Int avgSize = n/g2Ptr->mpisize;

                for(Int i=0;i<g2Ptr->mpisize-1;++i){
                  rcounts[i]=avgSize*sizeof(Scalar);
                  displs[i]=i*avgSize*sizeof(Scalar); 
                }
                rcounts[g2Ptr->mpisize-1]=(n-(g2Ptr->mpisize-1)*avgSize)*sizeof(Scalar);
                displs[g2Ptr->mpisize-1]=(g2Ptr->mpisize-1)*avgSize*sizeof(Scalar);


                MPI_Gatherv(bLocal->Data(), bLocal->m()*sizeof(Scalar), MPI_BYTE, xTrueGlobal->Data(),  rcounts,displs, MPI_BYTE, 0, world_comm);

                if( mpirank == 0 ){

                  Real absErrorOriginal, absErrorBcast, absErrorPipeline;
                  Real relErrorOriginal, relErrorBcast, relErrorPipeline;

                  if(doSinvPipeline){
                    absErrorPipeline = abs(col(gi) - (*xTrueGlobal)(gi,0));
                    relErrorPipeline = abs(absErrorPipeline/(*xTrueGlobal)(gi,0));
                  }

                  absErrorOriginal = abs(col_Original(gi) - (*xTrueGlobal)(gi,0));
                  relErrorOriginal = abs(absErrorOriginal/(*xTrueGlobal)(gi,0));


                  absErrorBcast = abs(colBcast(gi) - (*xTrueGlobal)(gi,0));
                  relErrorBcast = abs(absErrorBcast/(*xTrueGlobal)(gi,0));

                  if(doSinvPipeline){

                    cout<<std::endl<< "\t\tSolve \t SelInv_Original \t SelInvPipeline \t SelInv_Bcast"<<std::endl;
                    cout<< "Value \t\t"<<(*xTrueGlobal)(gi,0)<<"\t"<< col_Original(gi)<<"\t"<< col(gi)<< "\t" << colBcast(gi)<<std::endl;
                    cout<< "Rel error \t\t"<< "NA" <<"\t"<< relErrorOriginal<<"\t"<< relErrorPipeline<< "\t" << relErrorBcast<<std::endl;
                    cout<< "Abs error \t\t"<< "NA" <<"\t"<< absErrorOriginal<<"\t"<< absErrorPipeline<< "\t" << absErrorBcast<<std::endl<<std::endl;
                  }
                  else{

                    cout<<std::endl<< "\t\tSolve \t SelInv_Original \t SelInv_Bcast"<<std::endl;
                    cout<< "Value \t\t"<<(*xTrueGlobal)(gi,0)<<"\t"<< col_Original(gi)<<"\t" << colBcast(gi)<<std::endl;
                    cout<< "Rel error \t\t"<< "NA" <<"\t"<< relErrorOriginal<<"\t"<<  relErrorBcast<<std::endl;
                    cout<< "Abs error \t\t"<< "NA" <<"\t"<< absErrorOriginal<<"\t"<<  absErrorBcast<<std::endl<<std::endl;
                  }


                }     

                delete bGlobal;
                delete xTrueGlobal;
                delete xTrueLocal;
                delete bLocal;
                delete berr;



                delete GAPtr;
                delete AExplicitPtr;
              }



            }

          }
#endif

          if(doSinv_Original){
            delete PMlocRefPtr;
            delete superRefPtr;
            delete g3Ptr;
          }




          if(doSinvPipeline || doSinv_Bcast){

            PMatrix * PMloc = doSinvPipeline?PMlocPtr:PMlocBcastPtr;

            if(doToDist){
              // Convert to DistSparseMatrix and get the diagonal
              GetTime( timeSta );
              DistSparseMatrix<Scalar> Ainv;
              PMloc->PMatrixToDistSparseMatrix( Ainv );
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
              PMloc->PMatrixToDistSparseMatrix2( AMat, Ainv2 );
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
              mpi::Allreduce( &traceLocal, &trace, 1, MPI_SUM, world_comm );

              if( mpirank == 0 ){

                cout << "H.size = "  << HMat.size << endl;
                cout << std::endl << "Tr[Ainv2 * AMat] = " <<  trace << std::endl;
                statusOFS << std::endl << "Tr[Ainv2 * AMat] = " << std::endl << trace << std::endl;

                cout << std::endl << "H.size - Tr[Ainv2 * AMat] = " << (Real)HMat.size - (Real)trace.real() << std::endl;
                statusOFS << std::endl << "H.size - Tr[Ainv2 * AMat] = " << (Real)HMat.size- (Real)trace.real() << std::endl;

              }
            }

          }
          if(doSinv_Bcast){
            delete PMlocBcastPtr;
            delete superBcastPtr;
            delete g2Ptr;
          }



          if(doSinvPipeline){
            delete PMlocPtr;
            delete superPtr;
            delete g1Ptr;
          }

        }


      }

      delete luMatPtr;

      statusOFS.close();
    }
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
