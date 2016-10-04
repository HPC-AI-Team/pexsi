/// @file ex32.cpp
/// @brief Test to use the NGCHOL library in PEXSI.
/// @author Lin Lin
/// @date 2014-07-04

#include <mpi.h>

#include <complex>
#include <string>
#include <sstream>

// Libraries from NGCHOL BEGIN
#include  "ppexsi.hpp"

#include "sympack.hpp"
#include "sympack/symPACKMatrix.hpp"
#include "pexsi/sympack_interf.hpp"

using namespace PEXSI;
using namespace std;

void Usage(){
  std::cout 
    << "Test to convert a NGCHOL matrix to a PMatrix" 
    << std::endl << std::endl;
}

//typedef std::complex<double> SCALAR;
typedef double SCALAR;

int main(int argc, char **argv) 
{
  int mpisize;
  int mpirank;

  MPI_Init(&argc,&argv);
  //upcxx::init(&argc, &argv);

  //Create a communicator with npcol*nprow processors
  MPI_Comm world_comm;
  MPI_Comm_dup(MPI_COMM_WORLD, &world_comm);

  MPI_Comm_size(world_comm, &mpisize);
  MPI_Comm_rank(world_comm, &mpirank);

#if defined(SPROFILE) || defined(PMPI)
  symPACK::symPACK_set_main_args(argc,argv);
//  SYMPACK_SPROFILE_INIT(argc, argv);
#endif

  stringstream  ss;
  ss << "logTest" << mpirank;
  statusOFS.open( ss.str().c_str() );


  if( mpirank == 0 ) {
    Usage();
  }

  //try{

    //Temporarily required
    MPI_Comm_size(world_comm, &symPACK::np);
    MPI_Comm_rank(world_comm, &symPACK::iam);

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

    std::string Hfile;
    if( options.find("-H") != options.end() ){ 
      Hfile = options["-H"];
    }
    else{
      throw std::logic_error("Hfile must be provided.");
    }

    std::string format("CSC");
    if( options.find("-Hf") != options.end() ){ 
      format = options["-Hf"];
    }


    std::stringstream suffix;
    suffix<<mpirank;
    symPACK::logfileptr = new symPACK::LogFile("status",suffix.str().c_str());
    symPACK::logfileptr->OFS()<<"********* LOGFILE OF P"<<mpirank<<" *********"<<endl;
    symPACK::logfileptr->OFS()<<"**********************************"<<endl;


    mpisize = nprow*npcol;//mpisize;

    MPI_Comm workcomm;
    MPI_Comm_split(world_comm,mpirank<mpisize,mpirank,&workcomm);

    symPACK::symPACKOptions optionsFact;
    optionsFact.relax.SetMaxSize(200);
//    optionsFact.relax.SetMaxSize(1000);
//    optionsFact.relax.SetNrelax0(60);
//    optionsFact.relax.SetNrelax1(100);
//    optionsFact.relax.SetNrelax2(300);
    optionsFact.factorization = symPACK::FANBOTH;
    optionsFact.decomposition = symPACK::LDL;
    optionsFact.ordering = symPACK::METIS;
    //optionsFact.ordering = symPACK::NATURAL;
    optionsFact.scheduler = symPACK::DL;
    optionsFact.mappingTypeStr = "ROW2D";
    optionsFact.load_balance_str = "SUBCUBE-FO";
    optionsFact.MPIcomm = workcomm;

    //hide it in sympack ?
    upcxx::init(&argc, &argv);

//{
//  int error_code = 671345422;
//  char error_string[500];
//   int length_of_error_string, error_class;
//
//   MPI_Error_class(error_code, &error_class);
//   MPI_Error_string(error_class, error_string, &length_of_error_string);
//   fprintf(stderr, "%3d: %s\n", mpirank, error_string);
//   MPI_Error_string(error_code, error_string, &length_of_error_string);
//   fprintf(stderr, "%3d: %s\n", mpirank, error_string);
//}

    if(mpirank<mpisize){
      Real timeSta, timeEnd;
      Real timeTotalFactorizationSta, timeTotalFactorizationEnd;
      Real timeTotalSelInvSta, timeTotalSelInvEnd;
      MPI_Comm_size(workcomm,&symPACK::np);
      MPI_Comm_rank(workcomm,&symPACK::iam);


      

      symPACK::DistSparseMatrix<SCALAR> HMat(workcomm);
      symPACK::ReadMatrix<SCALAR,SCALAR>(Hfile, format, HMat);

      GetTime( timeTotalFactorizationSta );
      GetTime( timeSta );
      symPACK::symPACKMatrix<SCALAR>*  SMat = new symPACK::symPACKMatrix<SCALAR>();
      optionsFact.commEnv = new symPACK::CommEnvironment(workcomm);
      SMat->Init(optionsFact);
      SMat->SymbolicFactorization(HMat);
      GetTime( timeEnd );
      if( mpirank == 0 )
        cout << "Time for symbolic factorization is " << timeEnd - timeSta << " sec" << endl; 
      GetTime( timeSta );
      SMat->DistributeMatrix(HMat);
      GetTime( timeEnd );
      if( mpirank == 0 )
        cout << "Time for distribution is " << timeEnd - timeSta << " sec" << endl; 
      GetTime( timeSta );
      SMat->Factorize();
      GetTime( timeEnd );

      if( mpirank == 0 )
        cout << "Time for factorization is " << timeEnd - timeSta << " sec" << endl; 

      GetTime( timeTotalFactorizationEnd );
      if( mpirank == 0 )
        cout << "Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " sec" << endl; 

      GetTime( timeTotalSelInvSta );
      GridType *gPtr = new GridType( workcomm, nprow, npcol );
      SuperNodeType *superPtr = new SuperNodeType();

      // deprecated options
      SuperLUOptions luOpt;
      luOpt.ColPerm = "METIS_AT_PLUS_A";
      luOpt.numProcSymbFact = 1;

//#define UNSYM
//#ifdef UNSYM
//      luOpt.Symmetric = 0;
//#else
//      luOpt.Symmetric = 1;
//#endif
      luOpt.Symmetric = 1;
      luOpt.Transpose = 0;//transpose;

      PSelInvOptions selInvOpt;
      selInvOpt.maxPipelineDepth = -1;

      GetTime( timeSta );
      symPACKMatrixToSuperNode( *SMat, *superPtr );


      PMatrix<SCALAR> * pMat = PMatrix<SCALAR>::Create(gPtr,superPtr, &selInvOpt, &luOpt);
      symPACKMatrixToPMatrix( *SMat, *pMat );
      GetTime( timeEnd );

//      statusOFS << "super.perm = " << superPtr->perm << std::endl;
//      statusOFS << "super.permInv = " << superPtr->permInv << std::endl;
//      statusOFS << "super.perm_r = " << superPtr->perm_r << std::endl;
//      statusOFS << "super.permInv_r = " << superPtr->permInv_r << std::endl;
//      statusOFS << "super.superIdx = " << superPtr->superIdx << std::endl;
//      statusOFS << "super.etree = " << superPtr->etree << std::endl;
//      statusOFS << "ColBlockIdx = " << pMat->ColBlockIdx_ << std::endl;
//      statusOFS << "RowBlockIdx = " << pMat->RowBlockIdx_ << std::endl;
      if( mpirank == 0 )
         cout << "Time for converting symPACK matrix to PMatrix is " << timeEnd  - timeSta << endl;

      delete SMat;
      delete optionsFact.commEnv;

      MPI_Barrier(workcomm);

#if 0 //this is not compatible with the symmetric version of pexsi
{
  GridType *g = gPtr;
  Int nprow = g->numProcRow, npcol = g->numProcCol;
  Int mpirow = mpirank / npcol;  
  Int mpicol = mpirank % npcol;
  const SuperNodeType* super = pMat->SuperNode();
  Int numSuper = super->superPtr.m()-1; 
  for( Int ksup = 0; ksup < numSuper; ksup++ ){
    Int pcol = ( ksup % npcol );
    Int prow = ( ksup % nprow );
    if( mpirow == prow ){
      NumMat<SCALAR> DiagBuf;
      DiagBuf.Resize(SuperSize( ksup, super ), SuperSize( ksup, super ));
      if( mpicol == pcol ){
        Int jb = ksup / npcol;
        std::vector<LBlock<SCALAR> >& Lcol = pMat->L(jb);
        LBlock<SCALAR>& LB = Lcol[0];
        for(Int row = 0; row<LB.nzval.m(); row++){
          for(Int col = row+1; col<LB.nzval.n(); col++){
            LB.nzval(row,col) = LB.nzval(col,row) * LB.nzval(row,row);
          }
        } 
        DiagBuf = LB.nzval;
      }
#ifdef UNSYM
      //send the LBlock to the processors in same row
      MPI_Bcast(DiagBuf.Data(),DiagBuf.m()*DiagBuf.n()*sizeof(SCALAR),MPI_BYTE,pcol,g->rowComm);
      Int ib = ksup / nprow;
      std::vector<UBlock<SCALAR> >& Urow = pMat->U(ib);

      for(Int jblk = 0; jblk < Urow.size(); jblk++ ){
        UBlock<SCALAR>& UB = Urow[jblk];
        for(Int row = 0; row<UB.nzval.m(); row++){
          for(Int col = 0; col<UB.nzval.n(); col++){
            UB.nzval(row,col) = UB.nzval(row,col) * DiagBuf(row,row);
          }
        }
      }
#endif
    }
  }
//  pMat->DumpLU2();
}
#endif



#if 1
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


#if 1
        if( options.find("-A") != options.end() ){ 
          std::string Afile;
          if( options.find("-A") != options.end() ){ 
            Afile = options["-A"];
          }
          else{
            throw std::logic_error("Afile must be provided.");
          }





          DistSparseMatrix<SCALAR> AMat;
          ParaReadDistSparseMatrix( Afile.c_str(), AMat, workcomm ); 
          NumVec<SCALAR> diag;


          // Convert to DistSparseMatrix in the 2nd format and get the diagonal
          DistSparseMatrix<SCALAR> Ainv;
          SCALAR traceLocal;


          DistSparseMatrix<SCALAR> * Aptr;
          if(luOpt.Symmetric==0 && luOpt.Transpose==0){
            Aptr = new DistSparseMatrix<SCALAR>();
            //compute the transpose
            CSCToCSR(AMat,*Aptr);
          }
          else{
            Aptr = &AMat;
          }

          GetTime( timeSta );
          pMat->PMatrixToDistSparseMatrix( *Aptr, Ainv );
          GetTime( timeEnd );

          traceLocal = ZERO<SCALAR>();
          traceLocal = blas::Dotu( Aptr->nnzLocal, Ainv.nzvalLocal.Data(), 1,
              Aptr->nzvalLocal.Data(), 1 );

          if(luOpt.Symmetric==0 && luOpt.Transpose==0){
            delete Aptr;
          }


          if( mpirank == 0 )
            cout << "Time for converting PMatrix to DistSparseMatrix (2nd format) is " << timeEnd  - timeSta << endl;

          SCALAR trace = ZERO<SCALAR>();
          mpi::Allreduce( &traceLocal, &trace, 1, MPI_SUM, workcomm );

          if( mpirank == 0 ){

            cout << "A.size = "  << AMat.size << endl;
            cout << std::endl << "Tr[Ainv * AMat] = " <<  trace << std::endl;
            statusOFS << std::endl << "Tr[Ainv * AMat] = " << std::endl << trace << std::endl;
#ifdef _MYCOMPLEX_ 
            cout << std::endl << "|N - Tr[Ainv * AMat]| = " << std::abs( Complex(AMat.size, 0.0) - trace ) << std::endl;
            statusOFS << std::endl << "|N - Tr[Ainv * AMat]| = " << std::abs( Complex(AMat.size, 0.0) - trace ) << std::endl;
#else
            cout << std::endl << "|N - Tr[Ainv * AMat]| = " << std::abs( AMat.size - trace ) << std::endl;
            statusOFS << std::endl << "|N - Tr[Ainv * AMat]| = " << std::abs( AMat.size - trace ) << std::endl;
#endif
          }

          if( true ){
            NumVec<SCALAR> diag;

            GetTime( timeSta );
            pMat->GetDiagonal( diag );
            GetTime( timeEnd );


            if( mpirank == 0 )
              cout << "Time for getting the diagonal is " << timeEnd  - timeSta << endl;


            NumVec<SCALAR> diagDistSparse;
            GetTime( timeSta );
            GetDiagonal( Ainv, diagDistSparse );
            GetTime( timeEnd );
            if( mpirank == 0 )
              cout << "Time for getting the diagonal of DistSparseMatrix is " << timeEnd  - timeSta << endl;

            if( mpirank == 0 ){
              statusOFS << std::endl << "Diagonal of inverse from DistSparseMatrix format: " << std::endl << diagDistSparse << std::endl;
              Real diffNorm = 0.0;;
              for( Int i = 0; i < diag.m(); i++ ){
                diffNorm += pow( std::abs( diag(i) - diagDistSparse(i) ), 2.0 );
              }
              diffNorm = std::sqrt( diffNorm );
              cout << std::endl << "||diag - diagDistSparse||_2 = " << diffNorm << std::endl;
            }


            if( mpirank == 0 ){
              statusOFS << std::endl << "Diagonal of inverse in natural order: " << std::endl << diag << std::endl;
              ofstream ofs("diag");
              if( !ofs.good() ) 
                throw std::runtime_error("file cannot be opened.");
              serialize( diag, ofs, NO_MASK );
              ofs.close();
            }


            {
              NumVec<SCALAR> diagUNP;
              superPtr->perm = superPtr->perm_r;
              superPtr->permInv = superPtr->permInv_r;
              GetTime( timeSta );
              pMat->GetDiagonal( diagUNP );
              GetTime( timeEnd );


              if( mpirank == 0 )
                cout << "Time for getting the permuted diagonal is " << timeEnd  - timeSta << endl;
              statusOFS << std::endl << "Diagonal of inverse in permuted order: " << std::endl << diagUNP << std::endl;

            }



          }
        }
#endif
#endif






      delete pMat;


      delete superPtr;
      delete gPtr;
    }
    delete symPACK::logfileptr;
    statusOFS.close();

//  }
//  catch( std::exception& e )
//  {
//    std::cerr << "Processor " << mpirank << " caught exception with message: "
//      << e.what() << std::endl;
//  }



  //This will also finalize MPI
  upcxx::finalize();

  return 0;
}
