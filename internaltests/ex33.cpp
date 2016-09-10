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
#include "sympack/SupernodalMatrix.hpp"
#include "pexsi/sympack_interf.hpp"
#include <upcxx.h>

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
  upcxx::init(&argc, &argv);

  //Create a communicator with npcol*nprow processors
  MPI_Comm world_comm;
  MPI_Comm_dup(MPI_COMM_WORLD, &world_comm);

  MPI_Comm_size(world_comm, &mpisize);
  MPI_Comm_rank(world_comm, &mpirank);

  stringstream  ss;
  ss << "logTest" << mpirank;
  statusOFS.open( ss.str().c_str() );


  if( mpirank == 0 ) {
    Usage();
  }

  //try{

    //Temporarily required
    MPI_Comm_size(world_comm, &SYMPACK::np);
    MPI_Comm_rank(world_comm, &SYMPACK::iam);

    std::map<std::string,std::string> options;
    OptionsCreate(argc, argv, options);

    Int nprow = 1, npcol = mpisize;

    if( options.find("-r") != options.end() ){
      if( options.find("-c") != options.end() ){
        nprow= atoi(options["-r"].c_str());
        npcol= atoi(options["-c"].c_str());
        if(nprow*npcol != mpisize){
          throw std::runtime_error("The number of used processors can't be larger than the total number of available processors." );
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
          throw std::runtime_error("The number of used processors can't be larger than the total number of available processors." );
        } 
      }
      else{
        throw std::runtime_error( "When using -c option, -r also needs to be provided." );
      }
    }

    std::string Hfile;
    if( options.find("-H") != options.end() ){ 
      Hfile = options["-H"];
    }
    else{
      throw std::logic_error("Hfile must be provided.");
    }

    std::string Afile;
    if( options.find("-A") != options.end() ){ 
      Afile = options["-A"];
    }
    else{
      throw std::logic_error("Afile must be provided.");
    }


    std::stringstream suffix;
    suffix<<mpirank;
    SYMPACK::logfileptr = new SYMPACK::LogFile("status",suffix.str().c_str());
    SYMPACK::logfileptr->OFS()<<"********* LOGFILE OF P"<<mpirank<<" *********"<<endl;
    SYMPACK::logfileptr->OFS()<<"**********************************"<<endl;

    SYMPACK::NGCholOptions optionsFact;
    optionsFact.relax.SetMaxSize(10);
    optionsFact.maxSnode = 10;
    optionsFact.factorization = SYMPACK::FANBOTH;
    optionsFact.decomposition = SYMPACK::LDL;
    optionsFact.ordering = SYMPACK::MMD;
    optionsFact.scheduler = SYMPACK::DL;
    optionsFact.mappingTypeStr = "ROW2D";


    Int all_np = nprow*npcol;//mpisize;
    mpisize = optionsFact.used_procs(all_np);

    MPI_Comm workcomm;
    MPI_Comm_split(world_comm,mpirank<mpisize,mpirank,&workcomm);

    upcxx::team * workteam;
    Int new_rank = (mpirank<mpisize)?mpirank:mpirank-mpisize;
    upcxx::team_all.split(mpirank<mpisize,new_rank, workteam);

    if(mpirank<mpisize){
      Real timeSta, timeEnd;
      Real timeTotalFactorizationSta, timeTotalFactorizationEnd;
      Real timeTotalSelInvSta, timeTotalSelInvEnd;
      MPI_Comm_size(workcomm,&SYMPACK::np);
      MPI_Comm_rank(workcomm,&SYMPACK::iam);


      

      optionsFact.commEnv = new SYMPACK::CommEnvironment(workcomm);
      SYMPACK::DistSparseMatrix<SCALAR> HMat(workcomm);
      std::string format("CSC");
      SYMPACK::ReadMatrix<SCALAR,SCALAR>(Hfile, format, HMat);

      GetTime( timeTotalFactorizationSta );

      SYMPACK::SupernodalMatrix<SCALAR>*  SMat = new SYMPACK::SupernodalMatrix<SCALAR>();
      SMat->team_ = workteam;
      SMat->Init(HMat,optionsFact);

      GetTime( timeSta );
      SMat->Factorize();
      GetTime( timeEnd );

      if( mpirank == 0 )
        cout << "Time for factorization is " << timeEnd - timeSta << " sec" << endl; 

      GetTime( timeTotalFactorizationEnd );
      if( mpirank == 0 )
        cout << "Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " sec" << endl; 

      {
         Int nrhs = 1;
         Int n = HMat.size;
    std::vector<SCALAR> RHS,XTrue;
    if(nrhs>0){
      RHS.resize(n*nrhs);
      XTrue.resize(n*nrhs);

      //Initialize XTrue;
      Int val = 1.0;
      for(Int i = 0; i<n;++i){ 
        for(Int j=0;j<nrhs;++j){
          XTrue[i+j*n] = val;
          val = -val;
        }
      }


      timeSta = get_time();

        //TODO HANDLE MULTIPLE RHS
        SYMPACK::SparseMatrixStructure Local = HMat.GetLocalStructure();
        SYMPACK::SparseMatrixStructure Global;
        Local.ToGlobal(Global,workcomm);
        Global.ExpandSymmetric();

        Int numColFirst = std::max(1,n / mpisize);

        RHS.assign(n*nrhs,0.0);
        for(Int k = 0; k<nrhs; ++k){
          for(Int j = 0; j<n; ++j){
            Int iOwner = std::min(j/numColFirst,mpisize-1);
            if(mpirank == iOwner){
              Int iLocal = (j-(numColFirst)*iOwner);

              //do I own the column ?
              SCALAR t = XTrue[j+k*n];
              //do a dense mat mat mul ?
              for(Int ii = Local.colptr[iLocal]-1; ii< Local.colptr[iLocal+1]-1;++ii){
                Int row = Local.rowind[ii]-1;
                RHS[row+k*n] += t*HMat.nzvalLocal[ii];
                if(row>j){
                  RHS[j+k*n] += XTrue[row+k*n]*HMat.nzvalLocal[ii];
                }
              }
            }
          }
        }
        //Do a reduce of RHS
        mpi::Allreduce((SCALAR*)MPI_IN_PLACE,&RHS[0],RHS.size(),MPI_SUM,workcomm);

      timeEnd = get_time();
      if(mpirank==0){
        cout<<"spGEMM time: "<<timeEnd-timeSta<<endl;
      }

        /**************** SOLVE PHASE ***********/
        std::vector<SCALAR> XFinal;
        XFinal = RHS;

        SMat->Solve(&XFinal[0],nrhs);

        SMat->GetSolution(&XFinal[0],nrhs);




      std::vector<SCALAR> AX(n*nrhs,0.0);

      for(Int k = 0; k<nrhs; ++k){
        for(Int j = 1; j<=n; ++j){
          Int iOwner = std::min((j-1)/numColFirst,mpisize-1);
          if(mpirank == iOwner){
            Int iLocal = (j-(numColFirst)*iOwner);
            //do I own the column ?
            SCALAR t = XFinal[j-1+k*n];
            //do a dense mat mat mul ?
            for(Int ii = Local.colptr[iLocal-1]; ii< Local.colptr[iLocal];++ii){
              Int row = Local.rowind[ii-1];
              AX[row-1+k*n] += t*HMat.nzvalLocal[ii-1];
              if(row>j){
                AX[j-1+k*n] += XFinal[row-1+k*n]*HMat.nzvalLocal[ii-1];
              }
            }
          }
        }
      }

      //Do a reduce of RHS
      mpi::Allreduce((SCALAR*)MPI_IN_PLACE,&AX[0],AX.size(),MPI_SUM,workcomm);

      if(mpirank==0){
        blas::Axpy(AX.size(),-1.0,&RHS[0],1,&AX[0],1);
        double normAX = SYMPACK::lapack::Lange('F',n,nrhs,&AX[0],n);
        double normRHS = SYMPACK::lapack::Lange('F',n,nrhs,&RHS[0],n);
        cout<<"Norm of residual after SPCHOL is "<<normAX/normRHS<<std::endl;
      }


    }


      }


      GetTime( timeTotalSelInvSta );

      GridType *gPtr = new GridType( workcomm, nprow, npcol );
      SuperNodeType *superPtr = new SuperNodeType();

      // deprecated options
      SuperLUOptions luOpt;
      luOpt.ColPerm = "METIS_AT_PLUS_A";
      luOpt.numProcSymbFact = 1;

      PSelInvOptions selInvOpt;
      selInvOpt.maxPipelineDepth = -1;

      GetTime( timeSta );
      symPACKMatrixToSuperNode( *SMat, *superPtr );

      statusOFS << "super.permInv = " << superPtr->permInv << std::endl;
      statusOFS << "super.superIdx = " << superPtr->superIdx << std::endl;
      statusOFS << "super.etree = " << superPtr->etree << std::endl;

      PMatrix<SCALAR> PMat( gPtr, superPtr, &selInvOpt ,&luOpt );
      symPACKMatrixToPMatrix( *SMat, PMat );
      GetTime( timeEnd );

      statusOFS << "ColBlockIdx = " << PMat.ColBlockIdx_ << std::endl;
      statusOFS << "RowBlockIdx = " << PMat.RowBlockIdx_ << std::endl;
      if( mpirank == 0 )
         cout << "Time for converting symPACK matrix to PMatrix is " << timeEnd  - timeSta << endl;

      delete SMat;
      delete optionsFact.commEnv;

        SYMPACK::logfileptr->OFS() << "********************** L ************************" << std::endl;
      // Dump out the L factor
      for( Int jb = 0; jb < PMat.NumLocalBlockCol(); jb++ ){
        SYMPACK::logfileptr->OFS() << "------------ SuperNode " << GBj( jb, PMat.Grid() ) << std::endl;
        std::vector<LBlock<SCALAR> >& Lcol = PMat.L(jb);
        for( Int ib = 0; ib < PMat.NumBlockL(jb); ib++ ){
          LBlock<SCALAR>& LB = Lcol[ib];
          SYMPACK::logfileptr->OFS() << LB << std::endl;
        }
      }

        SYMPACK::logfileptr->OFS() << "********************** U ************************" << std::endl;
      // Dump out the U factor
      for( Int ib = 0; ib < PMat.NumLocalBlockRow(); ib++ ){
        SYMPACK::logfileptr->OFS() << "------------ SuperNode " << GBi( ib, PMat.Grid() ) << std::endl;
        std::vector<UBlock<SCALAR> >& Urow = PMat.U(ib);
        for( Int jb = 0; jb < PMat.NumBlockU(ib); jb++ ){
          UBlock<SCALAR>& UB = Urow[jb];
          SYMPACK::logfileptr->OFS() << UB << std::endl;
        }
      }

      // Preparation for the selected inversion
      GetTime( timeSta );
      PMat.ConstructCommunicationPattern();
      GetTime( timeEnd );
      if( mpirank == 0 )
        cout << "Time for constructing the communication pattern is " << timeEnd  - timeSta << endl;
      PMat.DumpLU();

      GetTime( timeSta );
      PMat.PreSelInv();
      GetTime( timeEnd );
      if( mpirank == 0 )
        cout << "Time for pre-selected inversion is " << timeEnd  - timeSta << endl;
      PMat.DumpLU();

      GetTime( timeSta );
      PMat.SelInv();
      GetTime( timeEnd );
      GetTime( timeTotalSelInvEnd );
      if( mpirank == 0 )
        cout << "Time for numerical selected inversion is " << timeEnd  - timeSta << endl;


      if( mpirank == 0 )
        cout << "Time for total selected inversion is " << timeTotalSelInvEnd  - timeTotalSelInvSta << endl;


      PMat.DumpLU();

      DistSparseMatrix<SCALAR> AMat;
      ParaReadDistSparseMatrix( Afile.c_str(), AMat, workcomm ); 
      NumVec<SCALAR> diag;

//      // Convert to DistSparseMatrix and get the diagonal
//      GetTime( timeSta );
//      DistSparseMatrix<SCALAR> Ainv;
//      PMat.PMatrixToDistSparseMatrix( Ainv );
//      GetTime( timeEnd );
//
//      if( mpirank == 0 )
//        cout << "Time for converting PMatrix to DistSparseMatrix is " << timeEnd  - timeSta << endl;
//
//      NumVec<SCALAR> diagDistSparse;
//      GetTime( timeSta );
//      GetDiagonal( Ainv, diagDistSparse );
//      GetTime( timeEnd );
//      if( mpirank == 0 )
//        cout << "Time for getting the diagonal of DistSparseMatrix is " << timeEnd  - timeSta << endl;
//
//      if( mpirank == 0 ){
//        statusOFS << std::endl << "Diagonal of inverse from DistSparseMatrix format : " << std::endl << diagDistSparse << std::endl;
//        Real diffNorm = 0.0;;
//        for( Int i = 0; i < diag.m(); i++ ){
//          diffNorm += pow( std::abs( diag(i) - diagDistSparse(i) ), 2.0 );
//        }
//        diffNorm = std::sqrt( diffNorm );
//        statusOFS << std::endl << "||diag - diagDistSparse||_2 = " << diffNorm << std::endl;
//      }

      // Convert to DistSparseMatrix in the 2nd format and get the diagonal
      GetTime( timeSta );
      DistSparseMatrix<SCALAR> Ainv2;
      PMat.PMatrixToDistSparseMatrix( AMat, Ainv2 );
      GetTime( timeEnd );

      if( mpirank == 0 )
        cout << "Time for converting PMatrix to DistSparseMatrix (2nd format) is " << timeEnd  - timeSta << endl;

      NumVec<SCALAR> diagDistSparse2;
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

      Complex traceLocal = blas::Dotu( AMat.nnzLocal, AMat.nzvalLocal.Data(), 1,
          Ainv2.nzvalLocal.Data(), 1 );
      Complex trace = Z_ZERO;
      mpi::Allreduce( &traceLocal, &trace, 1, MPI_SUM, world_comm );

      if( mpirank == 0 ){

        cout << "H.size = "  << AMat.size << endl;
        cout << std::endl << "Tr[Ainv2 * AMat] = " <<  trace << std::endl;
        statusOFS << std::endl << "Tr[Ainv2 * AMat] = " << std::endl << trace << std::endl;

        cout << std::endl << "|N - Tr[Ainv2 * AMat]| = " << std::abs( Complex(AMat.size, 0.0) - trace ) << std::endl;
        statusOFS << std::endl << "|N - Tr[Ainv2 * AMat]| = " << std::abs( Complex(AMat.size, 0.0) - trace ) << std::endl;

      }










      delete superPtr;
      delete gPtr;
    }
    delete SYMPACK::logfileptr;
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
