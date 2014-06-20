/// @file ex30.cpp
/// @brief Test the routine of MPI/OpenMP hybrid ISend/IRecv with the
/// MPI_THREAD_MULTIPLE mode
/// @author Lin Lin
/// @date 2014-06-20

#include  "environment.hpp"
#include  "ppexsi.hpp"

using namespace PEXSI;
using namespace std;

void Usage(){
  std::cout 
		<< "Testing the routine for estimating the spectral radius." << std::endl << std::endl;
}



int main(int argc, char **argv) 
{

  int provided;
  MPI_Init_thread( &argc, &argv, MPI_THREAD_MULTIPLE, &provided );
  int mpirank, mpisize;
  MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
  MPI_Comm_size( MPI_COMM_WORLD, &mpisize );

  if( mpirank == 0 ) {
    Usage();
  }


  try{
    if( provided != MPI_THREAD_MULTIPLE ){
      std::ostringstream msg;
      msg << "MPI thread level is " << provided << 
        ", which does not support MPI_THREAD_MULTIPLE." << std::endl;
      throw std::runtime_error( msg.str().c_str() );  
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
