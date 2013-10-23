/// @file ex27.cpp
/// @brief Testing the feature of MPI-3.
/// @author Lin Lin
/// @date 2013-10-22

#include  "environment_impl.hpp"
#include  "ppexsi.hpp"

using namespace PEXSI;
using namespace std;

void Usage(){
  std::cout 
		<< "Test the feature of MPI-3." << std::endl << std::endl;
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
		Real timeSta, timeEnd;
		MPI_Request request;
		MPI_Status  status;

		Int N = 100;
	  DblNumMat A1(N, N), A2(N, N);
		DblNumMat B1(N, N), B2(N, N);
		DblNumMat C1(N, N), C2(N, N);
		SetRandomSeed(1);
		UniformRandom( A1 );
		UniformRandom( A2 );
		SetValue( B1, 0.0 );
		SetValue( B2, 0.0 ); 
		SetValue( C1, 0.0 ); 
		SetValue( C2, 0.0 ); 


		MPI_Barrier( MPI_COMM_WORLD ); 
		GetTime( timeSta );
		blas::Gemm( 'N', 'N', N, N, N, 1.0, A1.Data(), 
				N, A2.Data(), N, 0.0, B1.Data(), N );
		MPI_Barrier( MPI_COMM_WORLD );
		GetTime( timeEnd );

		if( mpirank == 0 ){
			cout << "Wall clock time for Gemm = " << timeEnd - timeSta << endl;
		}
    
		MPI_Barrier( MPI_COMM_WORLD );
		GetTime( timeSta );
		MPI_Iallreduce( B1.Data(), C1.Data(), N*N, MPI_DOUBLE, MPI_SUM,  
				MPI_COMM_WORLD, &request );
		MPI_Wait( &request, &status );
		GetTime( timeEnd );
	
		if( mpirank == 0 ){
			cout << "Wall clock time for Communication = " << timeEnd - timeSta << endl;
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
