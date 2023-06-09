/// @file ex32.cpp
/// @brief Test to use the NGCHOL library in PEXSI.
/// @author Lin Lin
/// @date 2014-07-04

#include <mpi.h>
#include <omp.h>

#include <complex>
#include <string>
#include <sstream>

#include  "ppexsi.hpp"
#include "pexsi/timer.h"

#include "sympack.hpp"
#include "pexsi/sympack_interf.hpp"
#include <memory>

#define _MYCOMPLEX_

typedef double ISCALAR;
#ifdef _MYCOMPLEX_
typedef std::complex<double> SCALAR;
#else
typedef double SCALAR;
#endif




using namespace PEXSI;
using namespace std;

void Usage(){
  std::cout << "Usage" << std::endl << "run_pselinv_sympack -H <symmetric Hfile> -S [symmetric Sfile] -inf [H/S file format: CSC,HARWELL_BOEING] -r [nprow] -c [npcol] -rshift [real shift] -ishift [imaginary shift] -ToDist [doToDist] -Diag [doDiag] -SS [symmetricStorage]" << std::endl;
}


int main(int argc, char **argv) 
{
  int mpisize;
  int mpirank;

  symPACK_Init(&argc,&argv);

  //Create a communicator with npcol*nprow processors
  MPI_Comm world_comm;
  MPI_Comm_dup(MPI_COMM_WORLD, &world_comm);

  MPI_Comm_size(world_comm, &mpisize);
  MPI_Comm_rank(world_comm, &mpirank);

#if defined(PROFILE) || defined(PMPI)
  TAU_PROFILE_INIT(argc, argv);
#endif

#if defined(SPROFILE) || defined(PMPI)
  symPACK::symPACK_set_main_args(argc,argv);
  //  SYMPACK_SPROFILE_INIT(argc, argv);
#endif

  stringstream  ss;
  ss << "logTest" << mpirank;
  statusOFS.open( ss.str().c_str() );

#if defined (PROFILE) || defined(PMPI) || defined(USE_TAU)
  TAU_PROFILE_SET_CONTEXT(world_comm);
#endif

  if(argc<3){
    if( mpirank == 0 ) {
      Usage();
    }
  }


  try{

    //Temporarily required
    //MPI_Comm_size(world_comm, &symPACK::np);
    //MPI_Comm_rank(world_comm, &symPACK::iam);

    std::map<std::string,std::string> options;
    OptionsCreate(argc, argv, options);

    Int nprow = 1, npcol = mpisize;

    if( options.find("-r") != options.end() ){
      if( options.find("-c") != options.end() ){
        nprow= atoi(options["-r"].c_str());
        npcol= atoi(options["-c"].c_str());
      }
      else{
        throw std::runtime_error( "When using -r option, -c also needs to be provided." );
      }
    }
    else if( options.find("-c") != options.end() ){
      if( options.find("-r") != options.end() ){
        nprow= atoi(options["-r"].c_str());
        npcol= atoi(options["-c"].c_str());
      }
      else{
        throw std::runtime_error( "When using -c option, -r also needs to be provided." );
      }
    }
    if(nprow*npcol > mpisize){
      throw std::runtime_error("The number of used processors can't be larger than the total number of available processors." );
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


    std::string Hfile;
    if( options.find("-H") != options.end() ){ 
      Hfile = options["-H"];
    }
    else{
      throw std::logic_error("Hfile must be provided.");
    }

    std::string format("CSC");
    if( options.find("-inf") != options.end() ){ 
      format = options["-inf"];
    }
    int doToDist = 0;
    if( options.find("-ToDist") != options.end() ){ 
      doToDist= atoi(options["-ToDist"].c_str());
    }

    int doDiag = 0;
    if( options.find("-Diag") != options.end() ){ 
      doDiag = atoi(options["-Diag"].c_str());
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

    Int symmetricStorage = 0;
    if( options.find("-SS") != options.end() ){ 
      symmetricStorage = atoi(options["-SS"].c_str());
    }
    else{
      statusOFS << "-SS option is not given. " 
        << "Do not use symmetric storage." 
        << std::endl << std::endl;
    }



    int doFacto = 1;
    if( options.find("-F") != options.end() ){ 
      doFacto= atoi(options["-F"].c_str());
    }

    int doSymbolic = 1;
    if( options.find("-Symb") != options.end() ){ 
      doSymbolic= atoi(options["-Symb"].c_str());
    }

    int doConvert = 1;
    if( options.find("-C") != options.end() ){ 
      doConvert= atoi(options["-C"].c_str());
    }

    int doDistribute = 1;
    if( options.find("-D") != options.end() ){ 
      doDistribute= atoi(options["-D"].c_str());
    }

    Int doSelinv = 1;
    if( options.find("-Sinv") != options.end() ){ 
      doSelinv = atoi(options["-Sinv"].c_str());
    }

    if(doSelinv){
      doConvert=1;
    }

    if(doConvert){
      doFacto=1;
    }

    if(doFacto){
      doDistribute=1;
    }

    if(doDistribute){
      doSymbolic = 1;
    }

    Real rshift = 0.0, ishift = 0.0;
    if( options.find("-rshift") != options.end() ){ 
      rshift = atof(options["-rshift"].c_str());
    }
    if( options.find("-ishift") != options.end() ){ 
      ishift = atof(options["-ishift"].c_str());
    }


    std::string ColPerm;
    if( options.find("-colperm") != options.end() ){ 
      ColPerm = options["-colperm"];
    }
    else{
      statusOFS << "-colperm option is not given. " 
        << "Use MMD" 
        << std::endl << std::endl;
      ColPerm = "MMD";
    }



    std::string refinement;
    if( options.find("-refine") != options.end() ){ 
      refinement = options["-refine"];
    }
    else{
      refinement = "NO";
    }  




#ifdef _MYCOMPLEX_
    int isComplex = 0;
    if( options.find("-Complex") != options.end() ){ 
      isComplex = atoi(options["-Complex"].c_str());
    }
#endif

    //Initialize symPACK logfile
    std::stringstream suffix;
    suffix<<mpirank;
    symPACK::logfileptr = new symPACK::LogFile("status",suffix.str().c_str());
    symPACK::logfileptr->OFS()<<"********* LOGFILE OF P"<<mpirank<<" *********"<<endl;
    symPACK::logfileptr->OFS()<<"**********************************"<<endl;


    mpisize = nprow*npcol;

    MPI_Comm workcomm;
    MPI_Comm_split(world_comm,mpirank<mpisize,mpirank,&workcomm);

    symPACK::symPACKOptions optionsFact;
    optionsFact.NpOrdering = numProcSymbFact;
    optionsFact.decomposition = symPACK::DecompositionType::LDL;
    optionsFact.orderingStr = ColPerm;
    optionsFact.MPIcomm = workcomm;
    optionsFact.verbose=0;
    optionsFact.order_refinement_str = refinement;

    //  optionsFact.MPIcomm = worldcomm;
    optionsFact.factorization = symPACK::FANBOTH;
    optionsFact.factorization = symPACK::FANOUT;
    //symPACK::Multithreading::NumThread = omp_get_max_threads();
    //optionsFact.print_stats=false;

    //    optionsFact.load_balance_str = "NNZ";
    //    optionsFact.maxIsend = 100;
    //    optionsFact.relax.SetMaxSize(0);
    //    optionsFact.relax.SetMaxSize(60);


    //Initialize UPCXX for symPACK
    //upcxx::init(&argc, &argv);

    if(mpirank<mpisize){
      Real timeSta, timeEnd;
      Real timeTotalFactorizationSta, timeTotalFactorizationEnd;
      Real timeTotalSelInvSta, timeTotalSelInvEnd;
      //MPI_Comm_size(workcomm,&symPACK::np);
      //MPI_Comm_rank(workcomm,&symPACK::iam);

#ifdef _MYCOMPLEX_
      SCALAR zshift = SCALAR(rshift, ishift);
#else
      SCALAR zshift = SCALAR(rshift);
#endif

      symPACK::DistSparseMatrix<SCALAR> AMat(workcomm);
      //symPACK::symPACKMatrix<SCALAR>*  symPACKMat = NULL;
#ifdef _MYCOMPLEX_
      if(isComplex){
        symPACK::DistSparseMatrix<Complex> HMatComplex(workcomm);
        symPACK::ReadMatrix<Complex,Complex>(Hfile, format, HMatComplex);

        //Build AMat
        symPACK::DistSparseMatrix<Complex> SMat(workcomm);

        const symPACK::DistSparseMatrixGraph & Hgraph = HMatComplex.GetLocalGraph();
        // Get the diagonal indices for H and save it n diagIdxLocal_
        std::vector<Int>  diagIdxLocal;
        { 

          Int numColLocal      = Hgraph.LocalVertexCount();
          Int firstCol         = Hgraph.LocalFirstVertex();

          diagIdxLocal.clear();
          diagIdxLocal.reserve( HMatComplex.size );
          for( symPACK::Idx j = 0; j < numColLocal; j++ ){
            symPACK::Idx jcol = firstCol + j + 1;
            for( symPACK::Ptr i = Hgraph.colptr[j]-1; 
                i < Hgraph.colptr[j+1]-1; i++ ){
              symPACK::Idx irow = Hgraph.rowind[i];
              if( irow == jcol ){
                diagIdxLocal.push_back( i );
              }
            }
          } // for (j)
        }

        GetTime( timeSta );

        AMat.size          = HMatComplex.size;
        AMat.nnz           = HMatComplex.nnz;
        AMat.SetLocalGraph(HMatComplex.GetLocalGraph());
        AMat.nzvalLocal.resize( HMatComplex.nzvalLocal.size() );

        SCALAR *ptr0 = AMat.nzvalLocal.data();
        Complex *ptr1 = HMatComplex.nzvalLocal.data();
        Complex *ptr2 = SMat.nzvalLocal.data();


        if( SMat.size != 0 ){
          // S is not an identity matrix
          for( Int i = 0; i < HMatComplex.nzvalLocal.size(); i++ ){
            AMat.nzvalLocal[i] = HMatComplex.nzvalLocal[i] - zshift * SMat.nzvalLocal[i];
          }
        }
        else{
          // S is an identity matrix
          for( Int i = 0; i < Hgraph.LocalEdgeCount(); i++ ){
            AMat.nzvalLocal[i] = HMatComplex.nzvalLocal[i];
          }

          for( Int i = 0; i < diagIdxLocal.size(); i++ ){
            AMat.nzvalLocal[ diagIdxLocal[i] ] -= zshift;
          }
        } // if (SMat.size != 0 )

        if( mpirank == 0 ){
          cout << "nonzero in A (DistSparseMatrix format) = " << AMat.nnz << endl;
        }


        GetTime( timeEnd );
        if( mpirank == 0 )
          cout << "Time for constructing the matrix A is " << timeEnd - timeSta << endl;


      }
      else
#endif
      {
        symPACK::DistSparseMatrix<ISCALAR> HMatReal(workcomm);
        symPACK::ReadMatrix<ISCALAR,ISCALAR>(Hfile, format, HMatReal);

        //Build AMat
        symPACK::DistSparseMatrix<ISCALAR> SMat(workcomm);
        if( options.find("-S") != options.end() ){ 
          std::string Sfile;
          Sfile = options["-S"];
          symPACK::ReadMatrix<ISCALAR,ISCALAR>(Sfile, format, SMat);
        }
        else{
          statusOFS << "-S option is not given. " 
            << "Treat the overlap matrix as an identity matrix." 
            << std::endl << std::endl;
          SMat.size=0;
        }



        const symPACK::DistSparseMatrixGraph & Hgraph = HMatReal.GetLocalGraph();
        // Get the diagonal indices for H and save it n diagIdxLocal_
        std::vector<Int>  diagIdxLocal;
        { 

          Int numColLocal      = Hgraph.LocalVertexCount();
          Int firstCol         = Hgraph.LocalFirstVertex();

          diagIdxLocal.clear();
          diagIdxLocal.reserve( HMatReal.size );
          for( symPACK::Idx j = 0; j < numColLocal; j++ ){
            symPACK::Idx jcol = firstCol + j + 1;
            for( symPACK::Ptr i = Hgraph.colptr[j]-1; 
                i < Hgraph.colptr[j+1]-1; i++ ){
              symPACK::Idx irow = Hgraph.rowind[i];
              if( irow == jcol ){
                diagIdxLocal.push_back( i );
              }
            }
          } // for (j)
        }

        GetTime( timeSta );

        AMat.size          = HMatReal.size;
        AMat.nnz           = HMatReal.nnz;
        AMat.SetLocalGraph(HMatReal.GetLocalGraph());
        AMat.nzvalLocal.resize( HMatReal.nzvalLocal.size() );

        SCALAR *ptr0 = AMat.nzvalLocal.data();
        Real *ptr1 = HMatReal.nzvalLocal.data();
        Real *ptr2 = SMat.nzvalLocal.data();


        if( SMat.size != 0 ){
          // S is not an identity matrix
          for( Int i = 0; i < HMatReal.nzvalLocal.size(); i++ ){
            AMat.nzvalLocal[i] = HMatReal.nzvalLocal[i] - zshift * SMat.nzvalLocal[i];
          }
        }
        else{
          // S is an identity matrix
          for( Int i = 0; i < Hgraph.LocalEdgeCount(); i++ ){
            AMat.nzvalLocal[i] = HMatReal.nzvalLocal[i];
          }

          for( Int i = 0; i < diagIdxLocal.size(); i++ ){
            AMat.nzvalLocal[ diagIdxLocal[i] ] -= zshift;
          }
        } // if (SMat.size != 0 )

        if( mpirank == 0 ){
          cout << "nonzero in A (DistSparseMatrix format) = " << AMat.nnz << endl;
        }


        GetTime( timeEnd );
        if( mpirank == 0 )
          cout << "Time for constructing the matrix A is " << timeEnd - timeSta << endl;
      }

      PMatrix<SCALAR> * pMat = NULL; 
      FactorizationOptions factOpt;
      PSelInvOptions selInvOpt;
      GridType *gPtr = new GridType( workcomm, nprow, npcol );
      SuperNodeType *superPtr = new SuperNodeType();
      {
        std::shared_ptr<symPACK::symPACKMatrix<SCALAR> > symPACKMat(new symPACK::symPACKMatrix<SCALAR>());
        if(doSymbolic){
          GetTime( timeSta );
          symPACKMat->Init(optionsFact);
          symPACKMat->SymbolicFactorization(AMat);
          GetTime( timeEnd );
          if( mpirank == 0 )
            cout << "Time for symbolic factorization is " << timeEnd - timeSta << " sec" << endl; 
        }

        if(doDistribute){
          GetTime( timeTotalFactorizationSta );
          GetTime( timeSta );
          symPACKMat->DistributeMatrix(AMat);
          GetTime( timeEnd );
          if( mpirank == 0 )
            cout << "Time for distribution is " << timeEnd - timeSta << " sec" << endl; 
        }

        //MPI_Barrier(workcomm);

        if(doFacto){
          GetTime( timeSta );
          symPACKMat->Factorize();
          GetTime( timeEnd );

          if( mpirank == 0 )
            cout << "Time for factorization is " << timeEnd - timeSta << " sec" << endl; 

          GetTime( timeTotalFactorizationEnd );
          if( mpirank == 0 )
            cout << "Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " sec" << endl; 
          //      symPACKMat->DumpMatlab();

        }

        GetTime( timeTotalSelInvSta );

        factOpt.ColPerm = ColPerm;
        factOpt.Symmetric = 1;

        selInvOpt.maxPipelineDepth = maxPipelineDepth;
        selInvOpt.symmetricStorage = symmetricStorage;

        if (doConvert){
          GetTime( timeSta );
          symPACKMatrixToSuperNode( *symPACKMat, *superPtr );


          pMat = PMatrix<SCALAR>::Create(gPtr,superPtr, &selInvOpt, &factOpt);
          symPACKMatrixToPMatrix( *symPACKMat, *pMat );
          GetTime( timeEnd );

          if( mpirank == 0 )
            cout << "Time for converting symPACK matrix to PMatrix is " << timeEnd  - timeSta << endl;
        }

        //MPI_Barrier(workcomm);
      }

      if (doSelinv>0){
        // Preparation for the selected inversion
        GetTime( timeSta );
        pMat->ConstructCommunicationPattern();
        GetTime( timeEnd );
        if( mpirank == 0 )
          cout << "Time for constructing the communication pattern is " << timeEnd  - timeSta << endl;

        GetTime( timeSta );
        pMat->PreSelInv();
        GetTime( timeEnd );
        if( mpirank == 0 )
          cout << "Time for pre-selected inversion is " << timeEnd  - timeSta << endl;
        statusOFS << "Time for pre-selected inversion is " << timeEnd  - timeSta << endl;

        GetTime( timeSta );
        pMat->SelInv();
        GetTime( timeEnd );
        GetTime( timeTotalSelInvEnd );
        if( mpirank == 0 )
          cout << "Time for numerical selected inversion is " << timeEnd  - timeSta << endl;
        statusOFS << "Time for numerical selected inversion is " << timeEnd  - timeSta << endl;


        if( mpirank == 0 )
          cout << "Time for total selected inversion is " << timeTotalSelInvEnd  - timeTotalSelInvSta << endl;

#ifdef _PRINT_STATS_
        double flops = pMat->GetTotalFlops();
        if( mpirank == 0 )
          cout << "Total FLOPs for selected inversion is " << flops << endl;
#endif

        if( doDiag ){
          NumVec<SCALAR> diag;

          GetTime( timeSta );
          pMat->GetDiagonal( diag );
          GetTime( timeEnd );

          if( mpirank == 0 )
            cout << "Time for getting the diagonal is " << timeEnd  - timeSta << endl;


          if( mpirank == 0 ){
            statusOFS << std::endl << "Diagonal (pipeline) of inverse in natural order: " << std::endl << diag << std::endl;
            ofstream ofs("diag");
            if( !ofs.good() ) 
              ErrorHandling("file cannot be opened.");
            serialize( diag, ofs, NO_MASK );
            ofs.close();
          }
        }

        if(doToDist){
          //Expand AMat to unsymmetric storage
          GetTime( timeSta );
          AMat.ExpandSymmetric();
          AMat.SortGraph();
          AMat.GetLocalGraph().SetBaseval(1);
          GetTime( timeEnd );

          if( mpirank == 0 )
            cout << "Matrix A has been expanded in " << timeEnd  - timeSta << endl;

          DistSparseMatrix<SCALAR> AMatPattern;
          AMatPattern.size = AMat.size;
          AMatPattern.nnz = AMat.nnz;
          AMatPattern.colptrLocal.Resize(AMat.GetLocalGraph().colptr.size());
          std::copy(AMat.GetLocalGraph().colptr.begin(), AMat.GetLocalGraph().colptr.end(), &AMatPattern.colptrLocal[0]);
          AMatPattern.rowindLocal.Resize(AMat.GetLocalGraph().rowind.size());
          std::copy(AMat.GetLocalGraph().rowind.begin(), AMat.GetLocalGraph().rowind.end(), &AMatPattern.rowindLocal[0]);
          AMatPattern.nnzLocal = AMat.GetLocalGraph().rowind.size();

          // Convert to DistSparseMatrix in the 2nd format and get the diagonal
          DistSparseMatrix<SCALAR> Ainv;

          GetTime( timeSta );
          pMat->PMatrixToDistSparseMatrix( AMatPattern, Ainv );
          GetTime( timeEnd );


          if( mpirank == 0 )
            cout << "Time for converting PMatrix to DistSparseMatrix (2nd format) is " << timeEnd  - timeSta << endl;

          SCALAR traceLocal = ZERO<SCALAR>();
          traceLocal = blas::Dotu( Ainv.nnzLocal, Ainv.nzvalLocal.Data(), 1,
              AMat.nzvalLocal.data(), 1 );


          SCALAR trace = ZERO<SCALAR>();
          mpi::Allreduce( &traceLocal, &trace, 1, MPI_SUM, workcomm );

          if( mpirank == 0 ){
            cout << "A.size = "  << AMat.size << endl;
            cout << std::endl << "Tr[Ainv * AMat] = " <<  trace << std::endl;
            statusOFS << std::endl << "Tr[Ainv * AMat] = " << std::endl << trace << std::endl;

            cout << std::endl << "|N - Tr[Ainv * AMat]| = " << std::abs( SCALAR(AMat.size) - trace ) << std::endl;
            statusOFS << std::endl << "|N - Tr[Ainv * AMat]| = " << std::abs( SCALAR(AMat.size) - trace ) << std::endl;
          }
        }
      }

      delete pMat;

      delete superPtr;
      delete gPtr;
    }
    //    delete symPACK::logfileptr;

  }
  catch( std::exception& e )
  {
    std::cerr << "Processor " << mpirank << " caught exception with message: "
      << e.what() << std::endl;
  }

  statusOFS<<"Calling symPACK finalize"<<std::endl;
  symPACK_Finalize();
  statusOFS.close();
  //MPI_Finalize();

  return 0;
}
