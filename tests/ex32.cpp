/// @file ex32.cpp
/// @brief Test to use the NGCHOL library in PEXSI.
/// @author Lin Lin
/// @date 2014-07-04

#include <mpi.h>

#include <complex>

// Libraries from NGCHOL BEGIN
//#include <time.h>
//#include <random>
//#include <omp.h>

#include "ngchol.hpp"

//#include  "Environment.hpp"
//#include  "DistSparseMatrix.hpp"
////#include  "NumVec.hpp"
////#include  "NumMat.hpp"
//#include  "utility.hpp"
////#include  "ETree.hpp"
//#include  "Mapping.hpp"
//
////#include  "blas.hpp"
////#include  "lapack.hpp"
//////#include  "NZBlock.hpp"
//
//#include  "SuperNode.hpp"
//#include  "SupernodalMatrix.hpp"
//
//#include  "LogFile.hpp"



//using namespace LIBCHOLESKY;
using namespace PEXSI;
using namespace std;

void Usage(){
  std::cout 
    << "Test to use the NGCHOL library in PEXSI." 
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

  if( mpirank == 0 ) {
    Usage();
  }

  LIBCHOLESKY::logfileptr = new LIBCHOLESKY::LogFile(mpirank);
  LIBCHOLESKY::logfileptr->OFS()<<"********* LOGFILE OF P"<<mpirank<<" *********"<<endl;
  LIBCHOLESKY::logfileptr->OFS()<<"**********************************"<<endl;




  //LIBCHOLESKY::Modwrap2D * mapping = new LIBCHOLESKY::Modwrap2D(mpisize, mpisize, mpisize, 1);
  LIBCHOLESKY::Modwrap2D * mapping = new LIBCHOLESKY::Modwrap2D(mpisize, sqrt(mpisize), mpisize, 1);
  LIBCHOLESKY::DistSparseMatrix<SCALAR > HMat(world_comm);
  LIBCHOLESKY::ParaReadDistSparseMatrix( argv[1], HMat, world_comm ); 

  LIBCHOLESKY::SupernodalMatrix<SCALAR> SMat(HMat,-1,*mapping,0,0,world_comm);
  statusOFS<<"Hello"<<endl;
  for(Int i = 0; i<SMat.SupernodeCnt();++i){
    LIBCHOLESKY::SuperNode<SCALAR> curSnode = SMat.GetLocalSupernode(i);
    statusOFS<<curSnode<<endl;
  }

  delete mapping;

  try{
  }
  catch( std::exception& e )
  {
    std::cerr << "Processor " << mpirank << " caught exception with message: "
      << e.what() << std::endl;
#ifndef _RELEASE_
    DumpCallStack();
#endif
  }

  delete LIBCHOLESKY::logfileptr;
  statusOFS.close();


  MPI_Finalize();

  return 0;
}
