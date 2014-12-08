/*
   Copyright (c) 2012 The Regents of the University of California,
   through Lawrence Berkeley National Laboratory.  
   
   Author: Lin Lin and Mathias Jacquelin
	 
   This file is part of PEXSI. All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

   (1) Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
   (2) Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
   (3) Neither the name of the University of California, Lawrence Berkeley
   National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
   be used to endorse or promote products derived from this software without
   specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
   ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
   ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
   ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   You are under no obligation whatsoever to provide any bug fixes, patches, or
   upgrades to the features, functionality or performance of the source code
   ("Enhancements") to anyone; however, if you choose to make your Enhancements
   available either publicly, or directly to Lawrence Berkeley National
   Laboratory, without imposing a separate written license agreement for such
   Enhancements, then you hereby grant the following license: a non-exclusive,
   royalty-free perpetual license to install, use, modify, prepare derivative
   works, incorporate into other computer software, distribute, and sublicense
   such enhancements or derivative works thereof, in binary and source code form.
*/
/// @file utility_impl.hpp
/// @brief Implementation of utility subroutines.
/// @date 2012-09-27
#ifndef _PEXSI_UTILITY_IMPL_HPP_ 
#define _PEXSI_UTILITY_IMPL_HPP_

using namespace std;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::cerr;

namespace PEXSI{

// *********************************************************************
// Sparse Matrix
// *********************************************************************

//---------------------------------------------------------
inline void ReadSparseMatrix ( const char* filename, SparseMatrix<Real>& spmat )
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
inline void ReadDistSparseMatrix ( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm )
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

  // FIXME Maybe change to MPI_Comm_Dup.
  pspmat.comm = comm;

  // Read basic information
	if( mpirank == 0 ){
		fin.open(filename);
		if( !fin.good() ){
			#ifdef USE_ABORT
abort();
#endif
throw std::logic_error( "File cannot be opened!" );
		}
		fin.read((char*)&pspmat.size, sizeof(Int));
		fin.read((char*)&pspmat.nnz,  sizeof(Int));
	}

  // FIXME Maybe need LongInt format to read the number of nonzeros later

	MPI_Bcast(&pspmat.size, 1, MPI_INT, 0, comm);
	MPI_Bcast(&pspmat.nnz,  1, MPI_INT, 0, comm);

	// Read colptr
	IntNumVec  colptr(pspmat.size+1);
	if( mpirank == 0 ){
		Int tmp;
		fin.read((char*)&tmp, sizeof(Int));  

		if( tmp != pspmat.size+1 ){
			#ifdef USE_ABORT
abort();
#endif
throw std::logic_error( "colptr is not of the right size." );
		}

		fin.read((char*)colptr.Data(), sizeof(Int)*tmp);

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
		fin.read((char*)&tmp, sizeof(Int));  

		if( tmp != pspmat.nnz ){
			std::ostringstream msg;
			msg 
				<< "The number of nonzeros in row indices do not match." << std::endl
				<< "nnz = " << pspmat.nnz << std::endl
				<< "size of row indices = " << tmp << std::endl;
			#ifdef USE_ABORT
abort();
#endif
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
			#ifdef USE_ABORT
abort();
#endif
throw std::logic_error( msg.str().c_str() );
		}

    pspmat.rowindLocal.Resize( numRead );
		MPI_Recv( pspmat.rowindLocal.Data(), numRead, MPI_INT, 0, 1, comm, &mpistat );
	}
		
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
			#ifdef USE_ABORT
abort();
#endif
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
			#ifdef USE_ABORT
abort();
#endif
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

inline void ParaWriteDistSparseMatrix ( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm )
{
#ifndef _RELEASE_
  PushCallStack("ParaWriteDistSparseMatrix");
#endif
  // Get the processor information within the current communicator
  MPI_Barrier( comm );
  Int mpirank;  MPI_Comm_rank(comm, &mpirank);
  Int mpisize;  MPI_Comm_size(comm, &mpisize);
  MPI_Status mpistat;
  Int err = 0;



  int filemode = MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_UNIQUE_OPEN;

  MPI_File fout;
  MPI_Status status;



  err = MPI_File_open(comm,(char*) filename, filemode, MPI_INFO_NULL,  &fout);

  if (err != MPI_SUCCESS) {
    #ifdef USE_ABORT
abort();
#endif
throw std::logic_error( "File cannot be opened!" );
  }

  // FIXME Note that nnz uses the Int data type for consistency of writing / reading
  // Write header
  if( mpirank == 0 ){
    err = MPI_File_write_at(fout, 0,(char*)&pspmat.size, 1, MPI_INT, &status);
    err = MPI_File_write_at(fout, sizeof(Int),(char*)&pspmat.nnz, 1, MPI_INT, &status);
  }


  // Compute the number of columns on each processor
  Int numColLocal = pspmat.colptrLocal.m()-1;
  Int numColFirst = pspmat.size / mpisize;
  IntNumVec  colptrChunk(numColLocal+1);

  Int prev_nz = 0;
  MPI_Exscan(&pspmat.nnzLocal, &prev_nz, 1, MPI_INT, MPI_SUM, comm);

  for( Int i = 0; i < numColLocal + 1; i++ ){
    colptrChunk[i] = pspmat.colptrLocal[i] + prev_nz;
  }


  MPI_Datatype memtype, filetype;
  MPI_Aint disps[6];
  int blklens[6];
  MPI_Datatype types[6] = {MPI_INT,MPI_INT, MPI_INT,MPI_INT, MPI_INT,MPI_DOUBLE};

  /* set block lengths (same for both types) */
  blklens[0] = (mpirank==0)?1:0;
  blklens[1] = numColLocal+1;
  blklens[2] = (mpirank==0)?1:0;
  blklens[3] = pspmat.nnzLocal;
  blklens[4] = (mpirank==0)?1:0;
  blklens[5] = pspmat.nnzLocal;




  //Calculate offsets
  MPI_Offset myColPtrOffset, myRowIdxOffset, myNzValOffset;
  myColPtrOffset = 3*sizeof(int) + (mpirank*numColFirst)*sizeof(Int);
  myRowIdxOffset = 3*sizeof(int) + (pspmat.size +1  +  prev_nz)*sizeof(Int);
  myNzValOffset = 4*sizeof(int) + (pspmat.size +1 +  pspmat.nnz)*sizeof(Int)+ prev_nz*sizeof(Real);
  disps[0] = 2*sizeof(int);
  disps[1] = myColPtrOffset;
  disps[2] = myRowIdxOffset;
  disps[3] = sizeof(int)+myRowIdxOffset;
  disps[4] = myNzValOffset;
  disps[5] = sizeof(int)+myNzValOffset;



#if ( _DEBUGlevel_ >= 1 )
  char msg[200];
  char * tmp = msg;
  tmp += sprintf(tmp,"P%d ",mpirank);
  for(int i = 0; i<6; ++i){
    if(i==5)
      tmp += sprintf(tmp, "%d [%d - %d] | ",i,disps[i],disps[i]+blklens[i]*sizeof(double));
    else
      tmp += sprintf(tmp, "%d [%d - %d] | ",i,disps[i],disps[i]+blklens[i]*sizeof(int));
  }
  tmp += sprintf(tmp,"\n");
  printf("%s",msg);
#endif




  MPI_Type_create_struct(6, blklens, disps, types, &filetype);
  MPI_Type_commit(&filetype);

  /* create memory type */
  Int np1 = pspmat.size+1;
  MPI_Address( (void *)&np1,  &disps[0]);
  MPI_Address(colptrChunk.Data(), &disps[1]);
  MPI_Address( (void *)&pspmat.nnz,  &disps[2]);
  MPI_Address((void *)pspmat.rowindLocal.Data(),  &disps[3]);
  MPI_Address( (void *)&pspmat.nnz,  &disps[4]);
  MPI_Address((void *)pspmat.nzvalLocal.Data(),   &disps[5]);

  MPI_Type_create_struct(6, blklens, disps, types, &memtype);
  MPI_Type_commit(&memtype);



  /* set file view */
  err = MPI_File_set_view(fout, 0, MPI_BYTE, filetype, "native",MPI_INFO_NULL);

  /* everyone writes their own row offsets, columns, and 
   * data with one big noncontiguous write (in memory and 
   * file)
   */
  err = MPI_File_write_all(fout, MPI_BOTTOM, 1, memtype, &status);

  MPI_Type_free(&filetype);
  MPI_Type_free(&memtype);





  MPI_Barrier( comm );

  MPI_File_close(&fout);
#ifndef _RELEASE_
  PopCallStack();
#endif

  return ;
}		// -----  end of function ParaWriteDistSparseMatrix  ----- 

inline void ParaReadDistSparseMatrix ( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm )
{
#ifndef _RELEASE_
  PushCallStack("ParaReadDistSparseMatrix");
#endif
  // Get the processor information within the current communicator
  MPI_Barrier( comm );
  Int mpirank;  MPI_Comm_rank(comm, &mpirank);
  Int mpisize;  MPI_Comm_size(comm, &mpisize);
  MPI_Status mpistat;
  MPI_Datatype type;
  int lens[3];
  MPI_Aint disps[3];
  MPI_Datatype types[3];
  Int err = 0;

  int filemode = MPI_MODE_RDONLY | MPI_MODE_UNIQUE_OPEN;

  MPI_File fin;
  MPI_Status status;

  // FIXME Maybe change to MPI_Comm_Dup.
  pspmat.comm = comm;

  err = MPI_File_open(comm,(char*) filename, filemode, MPI_INFO_NULL,  &fin);

  if (err != MPI_SUCCESS) {
    #ifdef USE_ABORT
abort();
#endif
throw std::logic_error( "File cannot be opened!" );
  }

  // FIXME Note that nnz uses the Int data type for consistency of writing / reading
  // Read header
  if( mpirank == 0 ){
    err = MPI_File_read_at(fin, 0,(char*)&pspmat.size, 1, MPI_INT, &status);
    err = MPI_File_read_at(fin, sizeof(Int),(char*)&pspmat.nnz, 1, MPI_INT, &status);
  }


  /* define a struct that describes all our data */
  lens[0] = 1;
  lens[1] = 1;
  MPI_Address(&pspmat.size, &disps[0]);
  MPI_Address(&pspmat.nnz, &disps[1]);
  types[0] = MPI_INT;
  types[1] = MPI_INT;
  MPI_Type_struct(2, lens, disps, types, &type);
  MPI_Type_commit(&type);


  /* broadcast the header data to everyone */
  MPI_Bcast(MPI_BOTTOM, 1, type, 0, comm);

  MPI_Type_free(&type);

  // Compute the number of columns on each processor
  IntNumVec numColLocalVec(mpisize);
  Int numColLocal, numColFirst;
  numColFirst = pspmat.size / mpisize;
  SetValue( numColLocalVec, numColFirst );
  numColLocalVec[mpisize-1] = pspmat.size - numColFirst * (mpisize-1);  // Modify the last entry	
  numColLocal = numColLocalVec[mpirank];
  pspmat.colptrLocal.Resize( numColLocal + 1 );



  MPI_Offset myColPtrOffset = (2 + ((mpirank==0)?0:1) )*sizeof(int) + (mpirank*numColFirst)*sizeof(Int);

  Int np1 = 0;
  lens[0] = (mpirank==0)?1:0;
  lens[1] = numColLocal + 1;

  MPI_Address(&np1, &disps[0]);
  MPI_Address(pspmat.colptrLocal.Data(), &disps[1]);

  MPI_Type_hindexed(2, lens, disps, MPI_INT, &type);
  MPI_Type_commit(&type);

  err= MPI_File_read_at_all(fin, myColPtrOffset, MPI_BOTTOM, 1, type, &status);

  if (err != MPI_SUCCESS) {
#ifdef USE_ABORT
    abort();
#endif
    throw std::logic_error( "error reading colptr" );
  }
  MPI_Type_free(&type);

  // Calculate nnz_loc on each processor
  pspmat.nnzLocal = pspmat.colptrLocal[numColLocal] - pspmat.colptrLocal[0];


  pspmat.rowindLocal.Resize( pspmat.nnzLocal );
  pspmat.nzvalLocal.Resize ( pspmat.nnzLocal );

  //read rowIdx
  MPI_Offset myRowIdxOffset = (3 + ((mpirank==0)?0:1) )*sizeof(int) + (pspmat.size+1 + (pspmat.colptrLocal[0]-1))*sizeof(int);

  lens[0] = (mpirank==0)?1:0;
  lens[1] = pspmat.nnzLocal;

  MPI_Address(&np1, &disps[0]);
  MPI_Address(pspmat.rowindLocal.Data(), &disps[1]);

  MPI_Type_hindexed(2, lens, disps, MPI_INT, &type);
  MPI_Type_commit(&type);

  err= MPI_File_read_at_all(fin, myRowIdxOffset, MPI_BOTTOM, 1, type,&status);

  if (err != MPI_SUCCESS) {
    #ifdef USE_ABORT
abort();
#endif
throw std::logic_error( "error reading rowind" );
  }
  MPI_Type_free(&type);


  //read nzval
  MPI_Offset myNzValOffset = (4 + ((mpirank==0)?0:1) )*sizeof(int) + (pspmat.size+1 + pspmat.nnz)*sizeof(Int) + (pspmat.colptrLocal[0]-1)*sizeof(double);

  lens[0] = (mpirank==0)?1:0;
  lens[1] = pspmat.nnzLocal;

  MPI_Address(&np1, &disps[0]);
  MPI_Address(pspmat.nzvalLocal.Data(), &disps[1]);

  types[0] = MPI_INT;
  types[1] = MPI_DOUBLE;

  MPI_Type_create_struct(2, lens, disps, types, &type);
  MPI_Type_commit(&type);

  err = MPI_File_read_at_all(fin, myNzValOffset, MPI_BOTTOM, 1, type,&status);

  if (err != MPI_SUCCESS) {
    #ifdef USE_ABORT
abort();
#endif
throw std::logic_error( "error reading nzval" );
  }

  MPI_Type_free(&type);


  //convert to local references
  for( Int i = 1; i < numColLocal + 1; i++ ){
    pspmat.colptrLocal[i] = pspmat.colptrLocal[i] -  pspmat.colptrLocal[0] + 1;
  }
  pspmat.colptrLocal[0]=1;

  MPI_Barrier( comm );

  MPI_File_close(&fin);
#ifndef _RELEASE_
  PopCallStack();
#endif

  return ;
}		// -----  end of function ParaReadDistSparseMatrix  ----- 

inline void ReadDistSparseMatrixFormatted ( const char* filename, DistSparseMatrix<Real>& pspmat, MPI_Comm comm )
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

  // FIXME Maybe change to MPI_Comm_Dup.
  pspmat.comm = comm;

  // FIXME Maybe need LongInt format to read the number of nonzeros later

  // Read basic information
	if( mpirank == 0 ){
		fin.open(filename);
		if( !fin.good() ){
			#ifdef USE_ABORT
abort();
#endif
throw std::logic_error( "File cannot be opened!" );
		}
		Int dummy;

		fin >> pspmat.size >> dummy ;
		fin >> pspmat.nnz >> dummy;
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
			#ifdef USE_ABORT
abort();
#endif
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
			#ifdef USE_ABORT
abort();
#endif
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


inline void ParaReadDistSparseMatrix ( const char* filename, DistSparseMatrix<Complex>& pspmat, MPI_Comm comm )
{
#ifndef _RELEASE_
  PushCallStack("ParaReadDistSparseMatrix");
#endif
  // Get the processor information within the current communicator
  MPI_Barrier( comm );
  Int mpirank;  MPI_Comm_rank(comm, &mpirank);
  Int mpisize;  MPI_Comm_size(comm, &mpisize);
  MPI_Status mpistat;
  MPI_Datatype type;
  int lens[3];
  MPI_Aint disps[3];
  MPI_Datatype types[3];
  Int err = 0;

  int filemode = MPI_MODE_RDONLY | MPI_MODE_UNIQUE_OPEN;

  MPI_File fin;
  MPI_Status status;

  // FIXME Maybe change to MPI_Comm_Dup.
  pspmat.comm = comm;

  err = MPI_File_open(comm,(char*) filename, filemode, MPI_INFO_NULL,  &fin);

  if (err != MPI_SUCCESS) {
#ifdef USE_ABORT
    abort();
#endif
    throw std::logic_error( "File cannot be opened!" );
  }

  // FIXME Note that nnz uses the Int data type for consistency of writing / reading
  // Read header
  if( mpirank == 0 ){
    err = MPI_File_read_at(fin, 0,(char*)&pspmat.size, 1, MPI_INT, &status);
    err = MPI_File_read_at(fin, sizeof(Int),(char*)&pspmat.nnz, 1, MPI_INT, &status);
  }


  /* define a struct that describes all our data */
  lens[0] = 1;
  lens[1] = 1;
  MPI_Address(&pspmat.size, &disps[0]);
  MPI_Address(&pspmat.nnz, &disps[1]);
  types[0] = MPI_INT;
  types[1] = MPI_INT;
  MPI_Type_struct(2, lens, disps, types, &type);
  MPI_Type_commit(&type);


  /* broadcast the header data to everyone */
  MPI_Bcast(MPI_BOTTOM, 1, type, 0, comm);

  MPI_Type_free(&type);

  // Compute the number of columns on each processor
  IntNumVec numColLocalVec(mpisize);
  Int numColLocal, numColFirst;
  numColFirst = pspmat.size / mpisize;
  SetValue( numColLocalVec, numColFirst );
  numColLocalVec[mpisize-1] = pspmat.size - numColFirst * (mpisize-1);  // Modify the last entry	
  numColLocal = numColLocalVec[mpirank];
  pspmat.colptrLocal.Resize( numColLocal + 1 );



  MPI_Offset myColPtrOffset = (2 + ((mpirank==0)?0:1) )*sizeof(int) + (mpirank*numColFirst)*sizeof(Int);

  Int np1 = 0;
  lens[0] = (mpirank==0)?1:0;
  lens[1] = numColLocal + 1;

  MPI_Address(&np1, &disps[0]);
  MPI_Address(pspmat.colptrLocal.Data(), &disps[1]);

  MPI_Type_hindexed(2, lens, disps, MPI_INT, &type);
  MPI_Type_commit(&type);

  err= MPI_File_read_at_all(fin, myColPtrOffset, MPI_BOTTOM, 1, type, &status);

  if (err != MPI_SUCCESS) {
#ifdef USE_ABORT
    abort();
#endif
    throw std::logic_error( "error reading colptr" );
  }
  MPI_Type_free(&type);

  // Calculate nnz_loc on each processor
  pspmat.nnzLocal = pspmat.colptrLocal[numColLocal] - pspmat.colptrLocal[0];


  pspmat.rowindLocal.Resize( pspmat.nnzLocal );
  pspmat.nzvalLocal.Resize ( pspmat.nnzLocal );

  //read rowIdx
  MPI_Offset myRowIdxOffset = (3 + ((mpirank==0)?0:1) )*sizeof(int) + (pspmat.size+1 + pspmat.colptrLocal[0]-1)*sizeof(Int);

  lens[0] = (mpirank==0)?1:0;
  lens[1] = pspmat.nnzLocal;

  MPI_Address(&np1, &disps[0]);
  MPI_Address(pspmat.rowindLocal.Data(), &disps[1]);

  MPI_Type_hindexed(2, lens, disps, MPI_INT, &type);
  MPI_Type_commit(&type);

  err= MPI_File_read_at_all(fin, myRowIdxOffset, MPI_BOTTOM, 1, type,&status);

  if (err != MPI_SUCCESS) {
#ifdef USE_ABORT
    abort();
#endif
    throw std::logic_error( "error reading rowind" );
  }
  MPI_Type_free(&type);


  //read nzval
  MPI_Offset myNzValOffset = (4 + ((mpirank==0)?0:1) )*sizeof(int) + (pspmat.size+1 + pspmat.nnz)*sizeof(int) + (pspmat.colptrLocal[0]-1)*sizeof(Complex);

  lens[0] = (mpirank==0)?1:0;
  lens[1] = pspmat.nnzLocal;

  MPI_Address(&np1, &disps[0]);
  MPI_Address(pspmat.nzvalLocal.Data(), &disps[1]);

  types[0] = MPI_INT;
  types[1] = MPI_DOUBLE_COMPLEX;

  MPI_Type_create_struct(2, lens, disps, types, &type);
  MPI_Type_commit(&type);

  err = MPI_File_read_at_all(fin, myNzValOffset, MPI_BOTTOM, 1, type,&status);

  if (err != MPI_SUCCESS) {
#ifdef USE_ABORT
    abort();
#endif
    throw std::logic_error( "error reading nzval" );
  }

  MPI_Type_free(&type);


  //convert to local references
  for( Int i = 1; i < numColLocal + 1; i++ ){
    pspmat.colptrLocal[i] = pspmat.colptrLocal[i] -  pspmat.colptrLocal[0] + 1;
  }
  pspmat.colptrLocal[0]=1;

  MPI_Barrier( comm );

  MPI_File_close(&fin);
#ifndef _RELEASE_
  PopCallStack();
#endif

  return ;
}		// -----  end of function ParaReadDistSparseMatrix  ----- 

inline void ParaWriteDistSparseMatrix ( const char* filename, DistSparseMatrix<Complex>& pspmat, MPI_Comm comm )
{
#ifndef _RELEASE_
  PushCallStack("ParaWriteDistSparseMatrix");
#endif
  // Get the processor information within the current communicator
  MPI_Barrier( comm );
  Int mpirank;  MPI_Comm_rank(comm, &mpirank);
  Int mpisize;  MPI_Comm_size(comm, &mpisize);
  MPI_Status mpistat;
  Int err = 0;



  int filemode = MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_UNIQUE_OPEN;

  MPI_File fout;
  MPI_Status status;



  err = MPI_File_open(comm,(char*) filename, filemode, MPI_INFO_NULL,  &fout);

  if (err != MPI_SUCCESS) {
    #ifdef USE_ABORT
abort();
#endif
throw std::logic_error( "File cannot be opened!" );
  }

  // FIXME Note that nnz uses the Int data type for consistency of writing / reading
  // Write header
  if( mpirank == 0 ){
    err = MPI_File_write_at(fout, 0,(char*)&pspmat.size, 1, MPI_INT, &status);
    err = MPI_File_write_at(fout, sizeof(Int),(char*)&pspmat.nnz, 1, MPI_INT, &status);
  }


  // Compute the number of columns on each processor
  Int numColLocal = pspmat.colptrLocal.m()-1;
  Int numColFirst = pspmat.size / mpisize;
  IntNumVec  colptrChunk(numColLocal+1);

  Int prev_nz = 0;
  MPI_Exscan(&pspmat.nnzLocal, &prev_nz, 1, MPI_INT, MPI_SUM, comm);

  for( Int i = 0; i < numColLocal + 1; i++ ){
    colptrChunk[i] = pspmat.colptrLocal[i] + prev_nz;
  }


  MPI_Datatype memtype, filetype;
  MPI_Aint disps[6];
  int blklens[6];
  MPI_Datatype types[6] = {MPI_INT,MPI_INT, MPI_INT,MPI_INT, MPI_INT,MPI_DOUBLE_COMPLEX};

  /* set block lengths (same for both types) */
  blklens[0] = (mpirank==0)?1:0;
  blklens[1] = numColLocal+1;
  blklens[2] = (mpirank==0)?1:0;
  blklens[3] = pspmat.nnzLocal;
  blklens[4] = (mpirank==0)?1:0;
  blklens[5] = pspmat.nnzLocal;




  //Calculate offsets
  MPI_Offset myColPtrOffset, myRowIdxOffset, myNzValOffset;
  myColPtrOffset = 3*sizeof(int) + (mpirank*numColFirst)*sizeof(Int);
  myRowIdxOffset = 3*sizeof(int) + (pspmat.size +1  +  prev_nz)*sizeof(Int);
  myNzValOffset = 4*sizeof(int) + (pspmat.size +1 +  pspmat.nnz)*sizeof(Int)+ prev_nz*sizeof(Complex);
  disps[0] = 2*sizeof(int);
  disps[1] = myColPtrOffset;
  disps[2] = myRowIdxOffset;
  disps[3] = sizeof(int)+myRowIdxOffset;
  disps[4] = myNzValOffset;
  disps[5] = sizeof(int)+myNzValOffset;



#if ( _DEBUGlevel_ >= 1 )
  char msg[200];
  char * tmp = msg;
  tmp += sprintf(tmp,"P%d ",mpirank);
  for(int i = 0; i<6; ++i){
    if(i==5)
      tmp += sprintf(tmp, "%d [%d - %d] | ",i,disps[i],disps[i]+blklens[i]*sizeof(Complex));
    else
      tmp += sprintf(tmp, "%d [%d - %d] | ",i,disps[i],disps[i]+blklens[i]*sizeof(int));
  }
  tmp += sprintf(tmp,"\n");
  printf("%s",msg);
#endif




  MPI_Type_create_struct(6, blklens, disps, types, &filetype);
  MPI_Type_commit(&filetype);

  /* create memory type */
  Int np1 = pspmat.size+1;
  MPI_Address( (void *)&np1,  &disps[0]);
  MPI_Address(colptrChunk.Data(), &disps[1]);
  MPI_Address( (void *)&pspmat.nnz,  &disps[2]);
  MPI_Address((void *)pspmat.rowindLocal.Data(),  &disps[3]);
  MPI_Address( (void *)&pspmat.nnz,  &disps[4]);
  MPI_Address((void *)pspmat.nzvalLocal.Data(),   &disps[5]);

  MPI_Type_create_struct(6, blklens, disps, types, &memtype);
  MPI_Type_commit(&memtype);



  /* set file view */
  err = MPI_File_set_view(fout, 0, MPI_BYTE, filetype, "native",MPI_INFO_NULL);

  /* everyone writes their own row offsets, columns, and 
   * data with one big noncontiguous write (in memory and 
   * file)
   */
  err = MPI_File_write_all(fout, MPI_BOTTOM, 1, memtype, &status);

  MPI_Type_free(&filetype);
  MPI_Type_free(&memtype);





  MPI_Barrier( comm );

  MPI_File_close(&fout);
#ifndef _RELEASE_
  PopCallStack();
#endif

  return ;
}		// -----  end of function ParaWriteDistSparseMatrix  ----- 



inline void ReadDistSparseMatrixFormatted ( const char* filename, DistSparseMatrix<Complex>& pspmat, MPI_Comm comm )
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

  // FIXME Maybe change to MPI_Comm_Dup.
  pspmat.comm = comm;

  // FIXME Maybe need LongInt format to read the number of nonzeros later

  // Read basic information
	if( mpirank == 0 ){
		fin.open(filename);
		if( !fin.good() ){
			#ifdef USE_ABORT
abort();
#endif
throw std::logic_error( "File cannot be opened!" );
		}
		Int dummy;

		fin >> pspmat.size >> dummy;
		fin >> pspmat.nnz >> dummy;
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
#ifdef USE_ABORT
      abort();
#endif
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
			numRead = 2* (colptr[ip*numColFirst + numColLocalVec[ip]] - 
				colptr[ip*numColFirst]);
			Real *ptr;
      if( ip > 0 ){
        buf.Resize(numRead);
        ptr = buf.Data();
      }
      else{
        ptr = reinterpret_cast<Real*>(pspmat.nzvalLocal.Data());
      }

      //read data in either buf or pspmat.nzvalLocal
			for( Int i = 0; i < numRead; i++ ){
				fin >> *(ptr++);
			}

			if( ip > 0 ){
				MPI_Send(&numRead, 1, MPI_INT, ip, 0, comm);
				MPI_Send(buf.Data(), numRead, MPI_DOUBLE, ip, 1, comm);
			}
//			else{
//        std::copy( buf.Data(), buf.Data() + numRead,  pspmat.nzvalLocal.Data()  );
//        memcpy( (void*)( pspmat.nzvalLocal.Data() ),
//            (void*)buf.Data(), sizeof(Real)*numRead );
//			}
		}
	}
	else{
		Int numRead;
		MPI_Recv(&numRead, 1, MPI_INT, 0, 0, comm, &mpistat);
		if( numRead != 2*pspmat.nnzLocal ){
			std::ostringstream msg;
			msg << "The number of columns in values do not match." << std::endl
				<< "numRead    = " << numRead << std::endl
				<< "2*nnzLocal = " << 2*pspmat.nnzLocal << std::endl;
#ifdef USE_ABORT
      abort();
#endif
      throw std::logic_error( msg.str().c_str() );
		}

    pspmat.nzvalLocal.Resize( numRead/2 );
		MPI_Recv( reinterpret_cast<Real*>(pspmat.nzvalLocal.Data()), 
        numRead, MPI_DOUBLE, 0, 1, comm, &mpistat );
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


  template<typename T> inline void GetDiagonal ( const DistSparseMatrix<T>& A, NumVec<T>& diag )
{
#ifndef _RELEASE_
	PushCallStack("GetDiagonal");
#endif
	Int mpirank, mpisize;
	MPI_Comm_rank( A.comm, &mpirank );
	MPI_Comm_size( A.comm, &mpisize );

  NumVec<T>	 diagLocal( A.size );
	SetValue( diagLocal, ZERO<T>() );
	diag.Resize( A.size );
	SetValue( diag, ZERO<T>() );

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
			#ifdef USE_ABORT
abort();
#endif
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




}
#endif //_PEXSI_UTILITY_IMPL_HPP_
