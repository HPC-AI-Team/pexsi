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
  MPI_Comm_size(world_comm, &symPACK::np);
  MPI_Comm_rank(world_comm, &symPACK::iam);

  std::stringstream suffix;
  suffix<<mpirank;
  symPACK::logfileptr = new symPACK::LogFile("status",suffix.str().c_str());
  symPACK::logfileptr->OFS()<<"********* LOGFILE OF P"<<mpirank<<" *********"<<endl;
  symPACK::logfileptr->OFS()<<"**********************************"<<endl;




  symPACK::Modwrap2D * mapping = new symPACK::Modwrap2D(mpisize, mpisize, mpisize, 1);
  //symPACK::Modwrap2D * mapping = new symPACK::Modwrap2D(mpisize, sqrt(mpisize), mpisize, 1);
  symPACK::DistSparseMatrix<SCALAR > HMat(world_comm);
  symPACK::ParaReadDistSparseMatrix( argv[1], HMat, world_comm ); 


  //Prepare for a solve
  Int n = HMat.size;
  Int nrhs = 5;
  symPACK::DblNumMat RHS(n,nrhs);
  symPACK::DblNumMat XTrue(n,nrhs);
  symPACK::UniformRandom(XTrue);
  symPACK::sp_dgemm_dist("N","N", n, XTrue.n(), n, 
  symPACK::ONE<SCALAR>(), HMat, XTrue.Data(), XTrue.m(), 
  symPACK::ZERO<SCALAR>(), RHS.Data(), RHS.m());


  //Create the supernodal matrix
  symPACK::symPACKMatrix<SCALAR> SMat(HMat,-1,*mapping,0,0,world_comm);

  //sort X the same way (permute rows)
  symPACK::DblNumMat X(RHS.m(),RHS.n());
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
  symPACK::DblNumMat X2(X.m(),X.n());
  for(Int row = 1; row<= X.m(); ++row){
    for(Int col = 1; col<= X.n(); ++col){
      X2(SMat.perm_[row-1]-1,col-1) = X(row-1,col-1);
    }
  }

  //Check the error
  symPACK::blas::Axpy(X2.m()*X2.n(),-1.0,&XTrue(0,0),1,&X2(0,0),1);
  double norm = symPACK::lapack::Lange('F',X2.m(),X2.n(),&X2(0,0),X2.m());
  if(mpirank==0){
    cout<<"Norm of residual after SPCHOL is "<<norm<<std::endl;
  }

  //sample loop to convert it to PMatrix
  //for(Int i = 0; i<SMat.SupernodeCnt();++i){
  //  symPACK::SuperNode<SCALAR> curSnode = SMat.GetLocalSupernode(i);
  //}


  try{
  }
  catch( std::exception& e )
  {
    std::cerr << "Processor " << mpirank << " caught exception with message: "
      << e.what() << std::endl;
  }

  delete mapping;
  delete symPACK::logfileptr;
  statusOFS.close();


  MPI_Finalize();

  return 0;
}
