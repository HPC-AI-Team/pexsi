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

#include "ngchol.hpp"
#include "ngchol/sp_blas.hpp"
#include "pexsi/ngchol_interf.hpp"

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

  //Create a communicator with npcol*nprow processors
  MPI_Comm world_comm;
  MPI_Comm_dup(MPI_COMM_WORLD, &world_comm);

  MPI_Comm_size(world_comm, &mpisize);
  MPI_Comm_rank(world_comm, &mpirank);

  stringstream  ss;
  ss << "logTest" << mpirank;
  statusOFS.open( ss.str().c_str() );

  // FIXME
  Int nprow = 1, npcol = mpisize;
  if(argc>2){
    nprow = atoi(argv[2]);
    npcol = atoi(argv[3]);
  }

  if( mpirank == 0 ) {
    Usage();
  }

  try{

    //Temporarily required
    MPI_Comm_size(world_comm, &LIBCHOLESKY::np);
    MPI_Comm_rank(world_comm, &LIBCHOLESKY::iam);

    std::stringstream suffix;
    suffix<<mpirank;
    LIBCHOLESKY::logfileptr = new LIBCHOLESKY::LogFile("status",suffix.str().c_str());
    LIBCHOLESKY::logfileptr->OFS()<<"********* LOGFILE OF P"<<mpirank<<" *********"<<endl;
    LIBCHOLESKY::logfileptr->OFS()<<"**********************************"<<endl;


    LIBCHOLESKY::Modwrap2D * mapping = new LIBCHOLESKY::Modwrap2D(mpisize, mpisize, mpisize, 1);
    //LIBCHOLESKY::Modwrap2D * mapping = new LIBCHOLESKY::Modwrap2D(mpisize, sqrt(mpisize), mpisize, 1);
    LIBCHOLESKY::DistSparseMatrix<SCALAR > HMat(world_comm);
    LIBCHOLESKY::ParaReadDistSparseMatrix( argv[1], HMat, world_comm ); 

//    statusOFS << "size = " << HMat.size << std::endl;
//    statusOFS << "nzvalLocal = " << HMat.nzvalLocal << std::endl;

    LIBCHOLESKY::SupernodalMatrix<SCALAR> SMat(HMat,-1,*mapping,0,0,world_comm);


    GridType *gPtr = new GridType( world_comm, nprow, npcol );;
    SuperNodeType *superPtr = new SuperNodeType();

    // deprecated options
    SuperLUOptions luOpt;
    luOpt.ColPerm = "MMD_AT_PLUS_A";
    luOpt.maxPipelineDepth = -1;
    luOpt.numProcSymbFact = 1;

    NGCHOLMatrixToSuperNode( SMat, *superPtr );

    statusOFS << "super.permInv = " << superPtr->permInv << std::endl;

    statusOFS << "super.superIdx = " << superPtr->superIdx << std::endl;

    PMatrix<SCALAR> PMat( gPtr, superPtr, &luOpt );
    NGCHOLMatrixToPMatrix( SMat, PMat );

    // Dump out the L factor
    for( Int jb = 0; jb < PMat.NumLocalBlockCol(); jb++ ){
      statusOFS << "------------ SuperNode " << GBj( jb, PMat.Grid() ) << std::endl;
      std::vector<LBlock<SCALAR> >& Lcol = PMat.L(jb);
      for( Int ib = 0; ib < PMat.NumBlockL(jb); ib++ ){
        LBlock<SCALAR>& LB = Lcol[ib];
        statusOFS << LB << std::endl;
      }
    }

    // Dump out the U factor
    for( Int ib = 0; ib < PMat.NumLocalBlockRow(); ib++ ){
      statusOFS << "------------ SuperNode " << GBi( ib, PMat.Grid() ) << std::endl;
      std::vector<UBlock<SCALAR> >& Urow = PMat.U(ib);
      for( Int jb = 0; jb < PMat.NumBlockU(ib); jb++ ){
        UBlock<SCALAR>& UB = Urow[jb];
        statusOFS << UB << std::endl;
      }
    }



    delete superPtr;
    delete gPtr;
    delete mapping;
    delete LIBCHOLESKY::logfileptr;
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
