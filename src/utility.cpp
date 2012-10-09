#include "utility.hpp"

using namespace std;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::cerr;

namespace PEXSI{

// *********************************************************************
// IO functions
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
		throw std::logic_error( "File cannot be open!" );
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
		throw std::logic_error( "File cannot be open!" );
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
			throw std::logic_error( "File cannot be open!" );
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
			throw std::logic_error( "File cannot be open!" );
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
		throw std::logic_error( "File cannot be open!" );
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
// *********************************************************************

void ReadSparseMatrix ( const char* filename, SparseMatrix<Real>& spmat )
{
#ifndef _RELEASE_
	PushCallStack("ReadSparseMatrix");
#endif
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
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
}		// -----  end of function ReadSparseMatrix  ----- 

}  // namespace PEXSI
