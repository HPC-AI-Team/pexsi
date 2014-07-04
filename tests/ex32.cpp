/// @file ex32.cpp
/// @brief Test to use the NGCHOL library in PEXSI.
/// @author Lin Lin
/// @date 2014-07-04

#include <mpi.h>

// Libraries from NGCHOL BEGIN
#include <time.h>
#include <random>
#include <omp.h>

#include  "Environment.hpp"
#include  "DistSparseMatrix.hpp"
#include  "NumVec.hpp"
#include  "NumMat.hpp"
#include  "utility.hpp"
#include  "ETree.hpp"
#include  "Mapping.hpp"

#include  "blas.hpp"
#include  "lapack.hpp"
//#include  "NZBlock.hpp"
#include  "SuperNode.hpp"
#include  "SupernodalMatrix.hpp"

//#include <async.h>
#include  "LogFile.hpp"


extern "C" {
#include <bebop/util/config.h>
#include <bebop/smc/sparse_matrix.h>
#include <bebop/smc/csr_matrix.h>
#include <bebop/smc/csc_matrix.h>
#include <bebop/smc/sparse_matrix_ops.h>

#include <bebop/util/get_options.h>
#include <bebop/util/init.h>
#include <bebop/util/malloc.h>
#include <bebop/util/timer.h>
#include <bebop/util/util.h>
}
// Libraries from NGCHOL END



using namespace PEXSI;
using namespace std;

void Usage(){
  std::cout 
		<< "Test to use the NGCHOL library in PEXSI." 
    << std::endl << std::endl;
}



int main(int argc, char **argv) 
{

  if( mpirank == 0 ) {
    Usage();
  }


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

  MPI_Finalize();

  return 0;
}
