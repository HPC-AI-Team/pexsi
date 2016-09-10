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
#include "sympack/sp_blas.hpp"


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



  //Temporarily required
  MPI_Comm_size(world_comm, &SYMPACK::np);
  MPI_Comm_rank(world_comm, &SYMPACK::iam);

  std::stringstream suffix;
  suffix<<mpirank;
  SYMPACK::logfileptr = new SYMPACK::LogFile("status",suffix.str().c_str());
  SYMPACK::logfileptr->OFS()<<"********* LOGFILE OF P"<<mpirank<<" *********"<<endl;
  SYMPACK::logfileptr->OFS()<<"**********************************"<<endl;




  SYMPACK::Modwrap2D * mapping = new SYMPACK::Modwrap2D(mpisize, mpisize, mpisize, 1);
  //SYMPACK::Modwrap2D * mapping = new SYMPACK::Modwrap2D(mpisize, sqrt(mpisize), mpisize, 1);
  SYMPACK::DistSparseMatrix<SCALAR > HMat(world_comm);
  SYMPACK::ParaReadDistSparseMatrix( argv[1], HMat, world_comm ); 


  //Prepare for a solve
  Int n = HMat.size;
  Int nrhs = 5;
  SYMPACK::DblNumMat RHS(n,nrhs);
  SYMPACK::DblNumMat XTrue(n,nrhs);
  SYMPACK::UniformRandom(XTrue);
  SYMPACK::sp_dgemm_dist("N","N", n, XTrue.n(), n, 
  SYMPACK::ONE<SCALAR>(), HMat, XTrue.Data(), XTrue.m(), 
  SYMPACK::ZERO<SCALAR>(), RHS.Data(), RHS.m());


  //Create the supernodal matrix
  SYMPACK::SupernodalMatrix<SCALAR> SMat(HMat,-1,*mapping,0,0,world_comm);

  //sort X the same way (permute rows)
  SYMPACK::DblNumMat X(RHS.m(),RHS.n());
  for(Int row = 1; row<= RHS.m(); ++row){
    for(Int col = 1; col<= RHS.n(); ++col){
      X(row-1,col-1) = RHS(SMat.perm_[row-1]-1,col-1);
    }
  }


  //Factorize the matrix
  SMat.FanOut(world_comm);

  //Do a solve 
  SMat.Solve(&X,world_comm);
  SMat.GetSolution(X,world_comm);

  //Sort back X
  SYMPACK::DblNumMat X2(X.m(),X.n());
  for(Int row = 1; row<= X.m(); ++row){
    for(Int col = 1; col<= X.n(); ++col){
      X2(SMat.perm_[row-1]-1,col-1) = X(row-1,col-1);
    }
  }

  //Check the error
  SYMPACK::blas::Axpy(X2.m()*X2.n(),-1.0,&XTrue(0,0),1,&X2(0,0),1);
  double norm = SYMPACK::lapack::Lange('F',X2.m(),X2.n(),&X2(0,0),X2.m());
  if(mpirank==0){
    cout<<"Norm of residual after SPCHOL is "<<norm<<std::endl;
  }

  //sample loop to convert it to PMatrix
  //for(Int i = 0; i<SMat.SupernodeCnt();++i){
  //  SYMPACK::SuperNode<SCALAR> curSnode = SMat.GetLocalSupernode(i);
  //}


  try{
  }
  catch( std::exception& e )
  {
    std::cerr << "Processor " << mpirank << " caught exception with message: "
      << e.what() << std::endl;
  }

  delete mapping;
  delete SYMPACK::logfileptr;
  statusOFS.close();


  MPI_Finalize();

  return 0;
}
