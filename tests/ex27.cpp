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
		MPI_Request request;
		MPI_Status  status;

		Int N = 2000;
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

		
		// Blocking collective communcation
		{
			Real timeSta, timeEnd;
			Real timeTotalSta, timeTotalEnd;
			Real timeComp, timeComm;

			timeComp = 0.0;
			timeComm = 0.0;

			if( mpirank == 0 ){
				cout << "Testing the blocking collective communication" << endl;
			}

//			MPI_Barrier( MPI_COMM_WORLD ); 
			GetTime( timeSta );
			GetTime( timeTotalSta );
			blas::Gemm( 'N', 'N', N, N, N, 1.0, A1.Data(), 
					N, A2.Data(), N, 0.0, B1.Data(), N );
			GetTime( timeEnd );
			timeComp  += timeEnd - timeSta;

			GetTime( timeSta );
//			MPI_Allreduce( B1.Data(), C1.Data(), N*N, MPI_DOUBLE, MPI_SUM,  
//					MPI_COMM_WORLD );
			MPI_Reduce( B1.Data(), C1.Data(), N*N, MPI_DOUBLE, MPI_SUM, 0,
					MPI_COMM_WORLD );
			GetTime( timeEnd );
			timeComm  += timeEnd - timeSta;

			GetTime( timeSta );
			blas::Gemm( 'N', 'N', N, N, N, 1.0, A1.Data(), 
					N, A2.Data(), N, 0.0, B2.Data(), N );
			GetTime( timeEnd );
//			MPI_Barrier( MPI_COMM_WORLD );
			GetTime( timeTotalEnd );
			timeComp  += timeEnd - timeSta;

			if( mpirank == 0 ){
				cout << "Wall clock time for all work       = " << timeTotalEnd - timeTotalSta << endl;
				cout << "Wall clock time for computation    = " << timeComp << endl;
				cout << "Wall clock time for communication  = " << timeComm << endl;
				cout << "Wall clock time for communication (total - comp)  = " 
					<< timeTotalEnd - timeTotalSta - timeComp << endl;
			}

		}

		// MPI-3: Non-blocking collective communication
		{
			Real timeSta, timeEnd;
			Real timeTotalSta, timeTotalEnd;
			Real timeComp, timeComm;

			timeComp = 0.0;
			timeComm = 0.0;

			if( mpirank == 0 ){
				cout << endl << "Testing the non-blocking collective communication" << endl;
			}
		
	
//			MPI_Barrier( MPI_COMM_WORLD ); 
			GetTime( timeSta );
			GetTime( timeTotalSta );
			blas::Gemm( 'N', 'N', N, N, N, 1.0, A1.Data(), 
					N, A2.Data(), N, 0.0, B1.Data(), N );
			GetTime( timeEnd );
			timeComp  += timeEnd - timeSta;

//			MPI_Iallreduce( B1.Data(), C1.Data(), N*N, MPI_DOUBLE, MPI_SUM,  
//					MPI_COMM_WORLD, &request );
			MPI_Ireduce( B1.Data(), C1.Data(), N*N, MPI_DOUBLE, MPI_SUM, 0,
					MPI_COMM_WORLD, &request );


			GetTime( timeSta );
			blas::Gemm( 'N', 'N', N, N, N, 1.0, A1.Data(), 
					N, A2.Data(), N, 0.0, B2.Data(), N );
			GetTime( timeEnd );
			timeComp  += timeEnd - timeSta;
			GetTime( timeSta );
			MPI_Wait( &request, &status );
			GetTime( timeEnd );
			timeComm += timeEnd - timeSta;

//			MPI_Barrier( MPI_COMM_WORLD );

			GetTime( timeTotalEnd );


			if( mpirank == 0 ){
				cout << "Wall clock time for all work       = " << timeTotalEnd - timeTotalSta << endl;
				cout << "Wall clock time for computation    = " << timeComp << endl;
				cout << "Wall clock time for extra time of communication   = " << timeComm << endl;
				cout << "Wall clock time for communication (total - comp)  = " 
					<< timeTotalEnd - timeTotalSta - timeComp << endl;
			}
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
