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
    std::map<std::string,std::string> options;
    
    OptionsCreate(argc, argv, options);


    std::string Hfile, Sfile;

		if( options.find("-H") != options.end() ){ 
			Hfile = options["-H"];
		}
		else{
      throw std::logic_error("Hfile must be provided.");
		}

		if( options.find("-S") != options.end() ){
			Sfile = options["-S"];
		}
		else{
      if( mpirank == 0 ){
        std::cout << "-S option is not given. " 
          << "Treat the overlap matrix as an identity matrix." 
          << std::endl << std::endl;
      }
		}

    // Read the matrices

    DistSparseMatrix<Real>   HMat;
    DistSparseMatrix<Real>   SMat;


    ParaReadDistSparseMatrix( Hfile.c_str(), HMat, MPI_COMM_WORLD ); 
    if( Sfile.empty() ){
      SMat.size = 0;  // S is an identity matrix.
    }
    else{
      ParaReadDistSparseMatrix( Sfile.c_str(), SMat, MPI_COMM_WORLD ); 
    }

    // Estimate the spectral radius

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
    
    v0.Resize( HMat.size );

    UniformRandom( v0 );

    pexsi.EstimateSpectralRadius(
        0, HMat, SMat, v0, 1e-5, 1000, numIter, sigma );

    if( mpirank == 0 ) {
      std::cout << "iter  = " << numIter << std::endl;
      std::cout << "sigma = " << sigma << std::endl;
    }

  }
  catch( std::exception& e )
  {
    std::cerr << "Processor " << mpirank << " caught exception with message: "
      << e.what() << std::endl;
  }

  MPI_Finalize();

  return 0;
}
