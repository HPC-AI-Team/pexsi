/// @file ex30.cpp
/// @brief Test the routine of MPI/OpenMP hybrid ISend/IRecv with the
/// MPI_THREAD_MULTIPLE mode
/// @author Lin Lin
/// @date 2014-06-20

#include  "environment.hpp"
#include  "ppexsi.hpp"
#include  <omp.h>

using namespace PEXSI;
using namespace std;

void Usage(){
  std::cout 
		<< "Test the routine of MPI/OpenMP hybrid ISend/IRecv with the MPI_THREAD_MULTIPLE mode" 
    << std::endl << std::endl;
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
//      throw std::runtime_error( msg.str().c_str() );  
      if( mpirank == 0 ){
        std::cout << msg.str();
      }
    }

    stringstream  ss;
    ss << "logTest" << mpirank;
    statusOFS.open( ss.str().c_str() );

    Real errorMax = 0.0;
    Int ompsize;
#pragma omp parallel shared(errorMax, statusOFS, mpirank, mpisize)
    {
      ompsize = omp_get_num_threads(); 
    }

    stringstream ofs[ompsize];

#pragma omp parallel shared(errorMax, statusOFS, mpirank, mpisize, ofs)
    {
      Int omprank = omp_get_thread_num();
      Int ompsize = omp_get_num_threads();

      stringstream& msg = ofs[omprank];

      msg << "I have mpirank = " << mpirank << " and thread rank " << omprank 
        << ", thread size " << ompsize << std::endl;

      Int m = 100000, n = 1200/ompsize;
      DblNumMat smat(m, n), rmat(m, n);
      Real timeSta, timeEnd;

      timeSta = omp_get_wtime();
      SetValue( smat, 0.0 );
      SetValue( rmat, 0.0 );
      for( Int j = 0; j < n; j++ ){
        for( Int i = 0; i < m; i++){
          smat(i,j) = j + omprank;
        }
      }
      timeEnd = omp_get_wtime();


      msg << "Thread " << omprank << ": Generate a matrix of size " << 
        m << " x " << n << " takes " << timeEnd - timeSta << " s." << std::endl;

      Int src, dest;
      std::vector<MPI_Request> mpiReq;
      mpiReq.resize( ompsize, MPI_REQUEST_NULL );
      
      timeSta = omp_get_wtime();
      if( mpirank < mpisize / 2 ){
        dest = mpirank + mpisize/2;
        MPI_Isend( smat.Data(), m*n, MPI_DOUBLE, dest, omprank, MPI_COMM_WORLD, &mpiReq[omprank] );
      }
      else{
        src = mpirank - mpisize/2;
        MPI_Irecv( rmat.Data(), m*n, MPI_DOUBLE, src, omprank, MPI_COMM_WORLD, &mpiReq[omprank]);
      }
      mpi::Wait( mpiReq[omprank] );

      timeEnd = omp_get_wtime();


      msg << "Thread " << omprank << ": Send/Recv time for a matrix of size " << 
        m << " x " << n << " takes " << timeEnd - timeSta << " s." << std::endl;


      if( mpirank >= mpisize / 2 ){
        Real tmp, error;
        error = 0.0;
        Real *ptr1 = smat.Data(), *ptr2 = rmat.Data();
        for( Int i = 0; i < m; i++ ){
          tmp = *(ptr2++) - *(ptr1++);
          error += tmp * tmp;
        }
        error = std::sqrt( error );

#pragma critical
        {
          if( errorMax > error ) 
            errorMax = error;
        }
      }
    }

    for( Int i = 0; i < ompsize; i++ ){
      statusOFS << ofs[i].str();
    }

    if( mpirank >= mpisize / 2 ){
      statusOFS << "Max error for the received matrix is " << 
        errorMax << std::endl;
    }

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
