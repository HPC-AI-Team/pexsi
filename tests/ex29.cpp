/// @file ex29.cpp
/// @brief Testing the routine for estimating the spectral radius.
/// @author Lin Lin
/// @date 2014-02-20

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

  MPI_Init( &argc, &argv );
  int mpirank, mpisize;
  MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
  MPI_Comm_size( MPI_COMM_WORLD, &mpisize );

  if( mpirank == 0 ) {
    Usage();
  }


  try{
    DistSparseMatrix<Real>   AMat;
    DistSparseMatrix<Real>   BMat;

    ParaReadDistSparseMatrix( "H.csc", AMat, MPI_COMM_WORLD ); 

    BMat.size = 0;  // B is an identity matrix.
    DblNumVec v0;   // Random initial start by having v0.m() == 0
    Real sigma;
    Int numIter;

    GridType gridPole( MPI_COMM_WORLD, 1, mpisize );
    Int nprow = iround(std::sqrt( mpisize ));
    Int npcol = nprow;
    if( mpisize != nprow * npcol ){
      throw std::logic_error( "nprow != npcol" );
    }
    PPEXSIData pexsi( &gridPole, nprow, npcol );

    SetRandomSeed( 10 );
    
    v0.Resize( AMat.size );

    UniformRandom( v0 );

    pexsi.EstimateSpectralRadius(
        0, AMat, BMat, v0, 1e-3, 1000, numIter, sigma );

    if( mpirank == 0 ) {
      std::cout << "iter  = " << numIter << std::endl;
      std::cout << "sigma = " << sigma << std::endl;
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
