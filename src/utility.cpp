/// @file utility.cpp
/// @brief Implementation of the non-templated utility subroutines.
/// @author Lin Lin
/// @date 2012-09-20
#include "utility.hpp"

using namespace std;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::cerr;

namespace PEXSI{

// *********************************************************************
// IO functions
// TODO Move this to utility.hpp and make them inline functions
// *********************************************************************
//---------------------------------------------------------
Int SeparateRead(std::string name, std::istringstream& is)
{
#ifndef _RELEASE_
	PushCallStack("SeparateRead");
#endif
  MPI_Barrier(MPI_COMM_WORLD);
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  //
  char filename[100];
  sprintf(filename, "%s_%d_%d", name.c_str(), mpirank, mpisize);  //cerr<<filename<<endl;
  ifstream fin(filename);
	if( !fin.good() ){
		throw std::logic_error( "File cannot be openeded!" );
	}
 
 	is.str( std::string(std::istreambuf_iterator<char>(fin), std::istreambuf_iterator<char>()) );
  fin.close();
  //
  MPI_Barrier(MPI_COMM_WORLD);
#ifndef _RELEASE_
	PopCallStack();
#endif
  return 0;
}

//---------------------------------------------------------
Int SeparateWrite(std::string name, std::ostringstream& os)
{
#ifndef _RELEASE_
	PushCallStack("SeparateWrite");
#endif
   MPI_Barrier(MPI_COMM_WORLD);
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  //
  char filename[100];
  sprintf(filename, "%s_%d_%d", name.c_str(), mpirank, mpisize);
  ofstream fout(filename);
	if( !fout.good() ){
		throw std::logic_error( "File cannot be openeded!" );
	}
  fout<<os.str();
  fout.close();
  //
  MPI_Barrier(MPI_COMM_WORLD);
#ifndef _RELEASE_
	PopCallStack();
#endif
  return 0;
}

//---------------------------------------------------------
Int SharedRead(std::string name, std::istringstream& is)
{
#ifndef _RELEASE_
	PushCallStack("SharedRead");
#endif
  MPI_Barrier(MPI_COMM_WORLD);
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  //
  vector<char> tmpstr;
  if(mpirank==0) {
    ifstream fin(name.c_str());
		if( !fin.good() ){
			throw std::logic_error( "File cannot be openeded!" );
		}
    //std::string str(std::istreambuf_iterator<char>(fin), std::istreambuf_iterator<char>());
    //tmpstr.insert(tmpstr.end(), str.begin(), str.end());
    tmpstr.insert(tmpstr.end(), std::istreambuf_iterator<char>(fin), std::istreambuf_iterator<char>());
    fin.close();
    int size = tmpstr.size();	//cerr<<size<<endl;
    MPI_Bcast((void*)&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast((void*)&(tmpstr[0]), size, MPI_BYTE, 0, MPI_COMM_WORLD);
  } else {
    int size;
    MPI_Bcast((void*)&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    tmpstr.resize(size);
    MPI_Bcast((void*)&(tmpstr[0]), size, MPI_BYTE, 0, MPI_COMM_WORLD);
  }
  is.str( std::string(tmpstr.begin(), tmpstr.end()) );
  //
  MPI_Barrier(MPI_COMM_WORLD);
#ifndef _RELEASE_
	PopCallStack();
#endif
  return 0;
}

//---------------------------------------------------------
Int SharedWrite(std::string name, std::ostringstream& os)
{
#ifndef _RELEASE_
	PushCallStack("SharedWrite");
#endif
  MPI_Barrier(MPI_COMM_WORLD);
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  //
  if(mpirank==0) {
    ofstream fout(name.c_str());
		if( !fout.good() ){
			throw std::logic_error( "File cannot be openeded!" );
		}
    fout<<os.str();
    fout.close();
  }
  MPI_Barrier(MPI_COMM_WORLD);
#ifndef _RELEASE_
	PopCallStack();
#endif
  return 0;
}


//---------------------------------------------------------
Int SeparateWriteAscii(std::string name, std::ostringstream& os)
{
#ifndef _RELEASE_
	PushCallStack("SeparateWriteAscii");
#endif
  MPI_Barrier(MPI_COMM_WORLD);
  int mpirank;  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  int mpisize;  MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
  //
  char filename[100];
  sprintf(filename, "%s_%d_%d", name.c_str(), mpirank, mpisize);
  ofstream fout(filename, ios::trunc);
	if( !fout.good() ){
		throw std::logic_error( "File cannot be opened!" );
	}
  fout<<os;
  fout.close();
  //
  MPI_Barrier(MPI_COMM_WORLD);
#ifndef _RELEASE_
	PopCallStack();
#endif
  return 0;
}


// *********************************************************************
// Sparse Matrix
// TODO: Move to sparse_matrix_impl
// *********************************************************************

//---------------------------------------------------------
void ReadSparseMatrix ( const char* filename, SparseMatrix<Real>& spmat )
{
#ifndef _RELEASE_
	PushCallStack("ReadSparseMatrix");
#endif
	
	// FIXME
	// Binary format
	if( 1 ){
		std::istringstream iss;
		SharedRead( filename, iss );
		deserialize( spmat.size, iss, NO_MASK );
		deserialize( spmat.nnz,  iss, NO_MASK );
		deserialize( spmat.colptr, iss, NO_MASK );
		deserialize( spmat.rowind, iss, NO_MASK );
		deserialize( spmat.nzval, iss, NO_MASK );
	}
	
	// Ascii format
  if( 0 ) {	
		ifstream fin(filename);
		fin >> spmat.size >> spmat.nnz;

		spmat.colptr.Resize( spmat.size+1 );
		spmat.rowind.Resize( spmat.nnz );
		spmat.nzval.Resize ( spmat.nnz );

		for( Int i = 0; i < spmat.size + 1; i++ ){
			fin >> spmat.colptr(i);
		}

		for( Int i = 0; i < spmat.nnz; i++ ){
			fin >> spmat.rowind(i);
		}

		for( Int i = 0; i < spmat.nnz; i++ ){
			fin >> spmat.nzval(i);
		}

		fin.close();
	}
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
}		// -----  end of function ReadSparseMatrix  ----- 


//---------------------------------------------------------
void ReadDistSparseMatrix ( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm )
{
#ifndef _RELEASE_
	PushCallStack("ReadDistSparseMatrix");
#endif
	// Get the processor information within the current communicator
  MPI_Barrier( comm );
  Int mpirank;  MPI_Comm_rank(comm, &mpirank);
  Int mpisize;  MPI_Comm_size(comm, &mpisize);
	MPI_Status mpistat;
	std::ifstream fin;

  // Read basic information
	if( mpirank == 0 ){
		fin.open(filename);
		if( !fin.good() ){
			throw std::logic_error( "File cannot be openeded!" );
		}
		fin.read((char*)&pspmat.size, sizeof(Int));
		fin.read((char*)&pspmat.nnz,  sizeof(Int));
	}
	
	MPI_Bcast(&pspmat.size, 1, MPI_INT, 0, comm);
	MPI_Bcast(&pspmat.nnz,  1, MPI_INT, 0, comm);

	// Read colptr

	IntNumVec  colptr(pspmat.size+1);
	if( mpirank == 0 ){
		Int tmp;
		fin.read((char*)&tmp, sizeof(Int));  
		if( tmp != pspmat.size+1 ){
			throw std::logic_error( "colptr is not of the right size." );
		}
		fin.read((char*)colptr.Data(), sizeof(Int)*tmp);
	}

	MPI_Bcast(colptr.Data(), pspmat.size+1, MPI_INT, 0, comm);
//	std::cout << "Proc " << mpirank << " outputs colptr[end]" << colptr[pspmat.size] << endl;

	// Compute the number of columns on each processor
	IntNumVec numColLocalVec(mpisize);
	Int numColLocal, numColFirst;
	numColFirst = pspmat.size / mpisize;
  SetValue( numColLocalVec, numColFirst );
  numColLocalVec[mpisize-1] = pspmat.size - numColFirst * (mpisize-1);  // Modify the last entry	
	numColLocal = numColLocalVec[mpirank];

	pspmat.colptrLocal.Resize( numColLocal + 1 );
	for( Int i = 0; i < numColLocal + 1; i++ ){
		pspmat.colptrLocal[i] = colptr[mpirank * numColFirst+i] - colptr[mpirank * numColFirst] + 1;
	}

	// Calculate nnz_loc on each processor
	pspmat.nnzLocal = pspmat.colptrLocal[numColLocal] - pspmat.colptrLocal[0];

  pspmat.rowindLocal.Resize( pspmat.nnzLocal );
	pspmat.nzvalLocal.Resize ( pspmat.nnzLocal );

	// Read and distribute the row indices
	if( mpirank == 0 ){
		Int tmp;
		fin.read((char*)&tmp, sizeof(Int));  
		if( tmp != pspmat.nnz ){
			std::ostringstream msg;
			msg 
				<< "The number of nonzeros in row indices do not match." << std::endl
				<< "nnz = " << pspmat.nnz << std::endl
				<< "size of row indices = " << tmp << std::endl;
			throw std::logic_error( msg.str().c_str() );
		}
		IntNumVec buf;
		Int numRead;
		for( Int ip = 0; ip < mpisize; ip++ ){
			numRead = colptr[ip*numColFirst + numColLocalVec[ip]] - 
				colptr[ip*numColFirst];
			buf.Resize(numRead);
			fin.read( (char*)buf.Data(), numRead*sizeof(Int) );
			if( ip > 0 ){
				MPI_Send(&numRead, 1, MPI_INT, ip, 0, comm);
				MPI_Send(buf.Data(), numRead, MPI_INT, ip, 1, comm);
			}
			else{
        pspmat.rowindLocal = buf;
			}
		}
	}
	else{
		Int numRead;
		MPI_Recv(&numRead, 1, MPI_INT, 0, 0, comm, &mpistat);
		if( numRead != pspmat.nnzLocal ){
			std::ostringstream msg;
			msg << "The number of columns in row indices do not match." << std::endl
				<< "numRead  = " << numRead << std::endl
				<< "nnzLocal = " << pspmat.nnzLocal << std::endl;
			throw std::logic_error( msg.str().c_str() );
		}

    pspmat.rowindLocal.Resize( numRead );
		MPI_Recv( pspmat.rowindLocal.Data(), numRead, MPI_INT, 0, 1, comm, &mpistat );
	}
		
//	std::cout << "Proc " << mpirank << " outputs rowindLocal.size() = " 
//		<< pspmat.rowindLocal.m() << endl;


	// Read and distribute the nonzero values
	if( mpirank == 0 ){
		Int tmp;
		fin.read((char*)&tmp, sizeof(Int));  
		if( tmp != pspmat.nnz ){
			std::ostringstream msg;
			msg 
				<< "The number of nonzeros in values do not match." << std::endl
				<< "nnz = " << pspmat.nnz << std::endl
				<< "size of values = " << tmp << std::endl;
			throw std::logic_error( msg.str().c_str() );
		}
		NumVec<Real> buf;
		Int numRead;
		for( Int ip = 0; ip < mpisize; ip++ ){
			numRead = colptr[ip*numColFirst + numColLocalVec[ip]] - 
				colptr[ip*numColFirst];
			buf.Resize(numRead);
			fin.read( (char*)buf.Data(), numRead*sizeof(Real) );
			if( ip > 0 ){
				MPI_Send(&numRead, 1, MPI_INT, ip, 0, comm);
				MPI_Send(buf.Data(), numRead, MPI_DOUBLE, ip, 1, comm);
			}
			else{
        pspmat.nzvalLocal = buf;
			}
		}
	}
	else{
		Int numRead;
		MPI_Recv(&numRead, 1, MPI_INT, 0, 0, comm, &mpistat);
		if( numRead != pspmat.nnzLocal ){
			std::ostringstream msg;
			msg << "The number of columns in values do not match." << std::endl
				<< "numRead  = " << numRead << std::endl
				<< "nnzLocal = " << pspmat.nnzLocal << std::endl;
			throw std::logic_error( msg.str().c_str() );
		}

    pspmat.nzvalLocal.Resize( numRead );
		MPI_Recv( pspmat.nzvalLocal.Data(), numRead, MPI_DOUBLE, 0, 1, comm, &mpistat );
	}

	// Close the file
	if( mpirank == 0 ){
    fin.close();
	}



  MPI_Barrier( comm );

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
}		// -----  end of function ReadDistSparseMatrix  ----- 










void ParaReadDistSparseMatrix ( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm )
{
#ifndef _RELEASE_
	PushCallStack("ReadDistSparseMatrix");
#endif
	// Get the processor information within the current communicator
  MPI_Barrier( comm );
  Int mpirank;  MPI_Comm_rank(comm, &mpirank);
  Int mpisize;  MPI_Comm_size(comm, &mpisize);
	MPI_Status mpistat;
	std::ifstream fin;

  // Read basic information
	if( mpirank == 0 ){
		fin.open(filename);
		if( !fin.good() ){
			throw std::logic_error( "File cannot be openeded!" );
		}
		fin.read((char*)&pspmat.size, sizeof(Int));
		fin.read((char*)&pspmat.nnz,  sizeof(Int));
	}
	
	MPI_Bcast(&pspmat.size, 1, MPI_INT, 0, comm);
	MPI_Bcast(&pspmat.nnz,  1, MPI_INT, 0, comm);

	// Read colptr

	IntNumVec  colptr(pspmat.size+1);
	if( mpirank == 0 ){
		Int tmp;
		fin.read((char*)&tmp, sizeof(Int));  
		if( tmp != pspmat.size+1 ){
			throw std::logic_error( "colptr is not of the right size." );
		}
		fin.read((char*)colptr.Data(), sizeof(Int)*tmp);
	}

	MPI_Bcast(colptr.Data(), pspmat.size+1, MPI_INT, 0, comm);
//	std::cout << "Proc " << mpirank << " outputs colptr[end]" << colptr[pspmat.size] << endl;

	// Compute the number of columns on each processor
	IntNumVec numColLocalVec(mpisize);
	Int numColLocal, numColFirst;
	numColFirst = pspmat.size / mpisize;
  SetValue( numColLocalVec, numColFirst );
  numColLocalVec[mpisize-1] = pspmat.size - numColFirst * (mpisize-1);  // Modify the last entry	
	numColLocal = numColLocalVec[mpirank];

	pspmat.colptrLocal.Resize( numColLocal + 1 );
	for( Int i = 0; i < numColLocal + 1; i++ ){
		pspmat.colptrLocal[i] = colptr[mpirank * numColFirst+i] - colptr[mpirank * numColFirst] + 1;
	}

	// Calculate nnz_loc on each processor
	pspmat.nnzLocal = pspmat.colptrLocal[numColLocal] - pspmat.colptrLocal[0];

  pspmat.rowindLocal.Resize( pspmat.nnzLocal );
	pspmat.nzvalLocal.Resize ( pspmat.nnzLocal );

	// Read and distribute the row indices
	if( mpirank == 0 ){
		Int tmp;
		fin.read((char*)&tmp, sizeof(Int));  
		if( tmp != pspmat.nnz ){
			std::ostringstream msg;
			msg 
				<< "The number of nonzeros in row indices do not match." << std::endl
				<< "nnz = " << pspmat.nnz << std::endl
				<< "size of row indices = " << tmp << std::endl;
			throw std::logic_error( msg.str().c_str() );
		}
		IntNumVec buf;
		Int numRead;
		for( Int ip = 0; ip < mpisize; ip++ ){
			numRead = colptr[ip*numColFirst + numColLocalVec[ip]] - 
				colptr[ip*numColFirst];
			buf.Resize(numRead);
			fin.read( (char*)buf.Data(), numRead*sizeof(Int) );
			if( ip > 0 ){
				MPI_Send(&numRead, 1, MPI_INT, ip, 0, comm);
				MPI_Send(buf.Data(), numRead, MPI_INT, ip, 1, comm);
			}
			else{
        pspmat.rowindLocal = buf;
			}
		}
	}
	else{
		Int numRead;
		MPI_Recv(&numRead, 1, MPI_INT, 0, 0, comm, &mpistat);
		if( numRead != pspmat.nnzLocal ){
			std::ostringstream msg;
			msg << "The number of columns in row indices do not match." << std::endl
				<< "numRead  = " << numRead << std::endl
				<< "nnzLocal = " << pspmat.nnzLocal << std::endl;
			throw std::logic_error( msg.str().c_str() );
		}

    pspmat.rowindLocal.Resize( numRead );
		MPI_Recv( pspmat.rowindLocal.Data(), numRead, MPI_INT, 0, 1, comm, &mpistat );
	}
		
//	std::cout << "Proc " << mpirank << " outputs rowindLocal.size() = " 
//		<< pspmat.rowindLocal.m() << endl;


	// Read and distribute the nonzero values
	if( mpirank == 0 ){
		Int tmp;
		fin.read((char*)&tmp, sizeof(Int));  
		if( tmp != pspmat.nnz ){
			std::ostringstream msg;
			msg 
				<< "The number of nonzeros in values do not match." << std::endl
				<< "nnz = " << pspmat.nnz << std::endl
				<< "size of values = " << tmp << std::endl;
			throw std::logic_error( msg.str().c_str() );
		}
		NumVec<Real> buf;
		Int numRead;
		for( Int ip = 0; ip < mpisize; ip++ ){
			numRead = colptr[ip*numColFirst + numColLocalVec[ip]] - 
				colptr[ip*numColFirst];
			buf.Resize(numRead);
			fin.read( (char*)buf.Data(), numRead*sizeof(Real) );
			if( ip > 0 ){
				MPI_Send(&numRead, 1, MPI_INT, ip, 0, comm);
				MPI_Send(buf.Data(), numRead, MPI_DOUBLE, ip, 1, comm);
			}
			else{
        pspmat.nzvalLocal = buf;
			}
		}
	}
	else{
		Int numRead;
		MPI_Recv(&numRead, 1, MPI_INT, 0, 0, comm, &mpistat);
		if( numRead != pspmat.nnzLocal ){
			std::ostringstream msg;
			msg << "The number of columns in values do not match." << std::endl
				<< "numRead  = " << numRead << std::endl
				<< "nnzLocal = " << pspmat.nnzLocal << std::endl;
			throw std::logic_error( msg.str().c_str() );
		}

    pspmat.nzvalLocal.Resize( numRead );
		MPI_Recv( pspmat.nzvalLocal.Data(), numRead, MPI_DOUBLE, 0, 1, comm, &mpistat );
	}

	// Close the file
	if( mpirank == 0 ){
    fin.close();
	}



  MPI_Barrier( comm );

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
}		// -----  end of function ReadDistSparseMatrix  ----- 


















void ReadDistSparseMatrixFormatted ( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm )
{
#ifndef _RELEASE_
	PushCallStack("ReadDistSparseMatrixFormatted");
#endif
	// Get the processor information within the current communicator
  MPI_Barrier( comm );
  Int mpirank;  MPI_Comm_rank(comm, &mpirank);
  Int mpisize;  MPI_Comm_size(comm, &mpisize);
	MPI_Status mpistat;
	std::ifstream fin;

  // Read basic information
	if( mpirank == 0 ){
		fin.open(filename);
		if( !fin.good() ){
			throw std::logic_error( "File cannot be openeded!" );
		}
		Int dummy;
		fin >> pspmat.size >> dummy;
		fin >> pspmat.nnz;
		// FIXME this is temporary and only applies to 4*4 matrix.
//	  fin	>> dummy;
	}
	
	MPI_Bcast(&pspmat.size, 1, MPI_INT, 0, comm);
	MPI_Bcast(&pspmat.nnz,  1, MPI_INT, 0, comm);

	// Read colptr

	IntNumVec  colptr(pspmat.size+1);
	if( mpirank == 0 ){
		Int* ptr = colptr.Data();
		for( Int i = 0; i < pspmat.size+1; i++ )
			fin >> *(ptr++);
	}

	MPI_Bcast(colptr.Data(), pspmat.size+1, MPI_INT, 0, comm);

	// Compute the number of columns on each processor
	IntNumVec numColLocalVec(mpisize);
	Int numColLocal, numColFirst;
	numColFirst = pspmat.size / mpisize;
  SetValue( numColLocalVec, numColFirst );
  numColLocalVec[mpisize-1] = pspmat.size - numColFirst * (mpisize-1);  // Modify the last entry	
	numColLocal = numColLocalVec[mpirank];

	pspmat.colptrLocal.Resize( numColLocal + 1 );
	for( Int i = 0; i < numColLocal + 1; i++ ){
		pspmat.colptrLocal[i] = colptr[mpirank * numColFirst+i] - colptr[mpirank * numColFirst] + 1;
	}

	// Calculate nnz_loc on each processor
	pspmat.nnzLocal = pspmat.colptrLocal[numColLocal] - pspmat.colptrLocal[0];

  pspmat.rowindLocal.Resize( pspmat.nnzLocal );
	pspmat.nzvalLocal.Resize ( pspmat.nnzLocal );

	// Read and distribute the row indices
	if( mpirank == 0 ){
		Int tmp;
		IntNumVec buf;
		Int numRead;
		for( Int ip = 0; ip < mpisize; ip++ ){
			numRead = colptr[ip*numColFirst + numColLocalVec[ip]] - 
				colptr[ip*numColFirst];
			buf.Resize(numRead);
			Int *ptr = buf.Data();
			for( Int i = 0; i < numRead; i++ ){
				fin >> *(ptr++);
			}
			if( ip > 0 ){
				MPI_Send(&numRead, 1, MPI_INT, ip, 0, comm);
				MPI_Send(buf.Data(), numRead, MPI_INT, ip, 1, comm);
			}
			else{
        pspmat.rowindLocal = buf;
			}
		}
	}
	else{
		Int numRead;
		MPI_Recv(&numRead, 1, MPI_INT, 0, 0, comm, &mpistat);
		if( numRead != pspmat.nnzLocal ){
			std::ostringstream msg;
			msg << "The number of columns in row indices do not match." << std::endl
				<< "numRead  = " << numRead << std::endl
				<< "nnzLocal = " << pspmat.nnzLocal << std::endl;
			throw std::logic_error( msg.str().c_str() );
		}

    pspmat.rowindLocal.Resize( numRead );
		MPI_Recv( pspmat.rowindLocal.Data(), numRead, MPI_INT, 0, 1, comm, &mpistat );
	}
		
//	std::cout << "Proc " << mpirank << " outputs rowindLocal.size() = " 
//		<< pspmat.rowindLocal.m() << endl;


	// Read and distribute the nonzero values
	if( mpirank == 0 ){
		Int tmp;
		NumVec<Real> buf;
		Int numRead;
		for( Int ip = 0; ip < mpisize; ip++ ){
			numRead = colptr[ip*numColFirst + numColLocalVec[ip]] - 
				colptr[ip*numColFirst];
			buf.Resize(numRead);
			Real *ptr = buf.Data();
			for( Int i = 0; i < numRead; i++ ){
				fin >> *(ptr++);
			}
			if( ip > 0 ){
				MPI_Send(&numRead, 1, MPI_INT, ip, 0, comm);
				MPI_Send(buf.Data(), numRead, MPI_DOUBLE, ip, 1, comm);
			}
			else{
        pspmat.nzvalLocal = buf;
			}
		}
	}
	else{
		Int numRead;
		MPI_Recv(&numRead, 1, MPI_INT, 0, 0, comm, &mpistat);
		if( numRead != pspmat.nnzLocal ){
			std::ostringstream msg;
			msg << "The number of columns in values do not match." << std::endl
				<< "numRead  = " << numRead << std::endl
				<< "nnzLocal = " << pspmat.nnzLocal << std::endl;
			throw std::logic_error( msg.str().c_str() );
		}

    pspmat.nzvalLocal.Resize( numRead );
		MPI_Recv( pspmat.nzvalLocal.Data(), numRead, MPI_DOUBLE, 0, 1, comm, &mpistat );
	}

	// Close the file
	if( mpirank == 0 ){
    fin.close();
	}



  MPI_Barrier( comm );

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
}		// -----  end of function ReadDistSparseMatrixFormatted  ----- 

void
GetDiagonal ( const DistSparseMatrix<Complex>& A, 
		NumVec<Complex>& diag )
{
#ifndef _RELEASE_
	PushCallStack("GetDiagonal");
#endif
	Int mpirank, mpisize;
	MPI_Comm_rank( A.comm, &mpirank );
	MPI_Comm_size( A.comm, &mpisize );

  NumVec<Complex>	 diagLocal( A.size );
	SetValue( diagLocal, Z_ZERO );
	diag.Resize( A.size );
	SetValue( diag, Z_ZERO );

	Int numColFirst = A.size / mpisize;
	Int firstCol    = mpirank * numColFirst;
	Int numColLocal = A.colptrLocal.m() - 1;

#if ( _DEBUGlevel_ >= 1 )
	statusOFS << "numColFirst = " << numColFirst << std::endl;
	statusOFS << "A.nzvalLocal.size = " << A.nzvalLocal.m() << std::endl;
	statusOFS << "A.nnzLocal = " << A.nnzLocal << std::endl;
#endif

	// Note that the indices in DistSparseMatrix follows the FORTRAN convention
  for( Int j = 0; j < numColLocal; j++ ){
		Int jcol = j + firstCol + 1;
		Int numRow = A.colptrLocal(j+1) - A.colptrLocal(j);
		const Int* rowPtr = &A.rowindLocal( A.colptrLocal(j) - 1 );
		// NOTE: The rows in DistSparseMatrix are not necessarily ordered.
		// So lower_bound cannot be used here for fast searching. find has to be used. 
		const Int* ptr = find( rowPtr, rowPtr + numRow, jcol ); 
		if( ptr == rowPtr + numRow ){
			std::ostringstream msg;
			msg << "Serious problem. Did not find the row corresponding to the column." << std::endl
				<< "This happens when j = " << j << ", jcol = " << jcol << ", and the row indices are " << std::endl
				<< IntNumVec( numRow, false, const_cast<Int*>(rowPtr) ) << std::endl;
			throw std::logic_error( msg.str().c_str() );
		}
		Int diagIdx = ptr - A.rowindLocal.Data();
    diagLocal( jcol - 1 ) = A.nzvalLocal( diagIdx );
	}

	mpi::Allreduce( &diagLocal[0], &diag[0], A.size, MPI_SUM, A.comm );

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
}		// -----  end of function GetDiagonal  ----- 

// *********************************************************************
// Other numerical routines
// *********************************************************************

// Interpolation
/// @brief Linear interpolates from (x,y) to (xx,yy)
///
/// Note: 
///
/// x and xx must be sorted in ascending order.
///
/// if xx[i] < x[0],     yy[i] = y[0]
///    xx[i] > x[end-1], yy[i] = y[end-1]
void
LinearInterpolation ( 
		const std::vector<Real>& x, 
		const std::vector<Real>& y,
		const std::vector<Real>& xx,
		std::vector<Real>& yy )
{
#ifndef _RELEASE_
	PushCallStack("LinearInterpolation");
#endif
	Int numX  = x.size();
	Int numXX = xx.size();

  for( Int i = 1; i < numX; i++ ){
		if( x[i] <= x[i-1] ) 
			throw std::runtime_error("x must be sorted strictly ascendingly.");
	}

	for( Int i = 1; i < numXX; i++){
		if( xx[i] < xx[i-1] )
			throw std::runtime_error("xx must be sorted ascendingly.");
	}


	yy.resize( numXX );
	std::vector<Real>::const_iterator vi;
	Int ix;
	for( Int i = 0; i < numXX; i++ ){
		if( xx[i] <= x[0] ){
			yy[i] = y[0];
		}
		else if( xx[i] >= x[numX-1] ){
			yy[i] = y[numX-1];
		}
		else{
			vi = std::lower_bound( x.begin(), x.end(), xx[i] );
			ix = vi - x.begin();

			yy[i] = y[ix-1] + (y[ix] - y[ix-1]) / (x[ix] - x[ix-1]) 
				* (xx[i] - x[ix-1]);
		}
	} // for (i)

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
}		// -----  end of function LinearInterpolation  ----- 



}  // namespace PEXSI
