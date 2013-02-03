/// @file interface.cpp
/// @brief Interface subroutines of PPEXSI that can be called by both C and FORTRAN.
/// @author Lin Lin
/// @date 2013-02-03
// TODO
//#include "c_pexsi_interface.h"
#include "ppexsi.hpp"
#include "blas.hpp"

// Error handling used in the C interface that is different from the
// throw/catch system.
#define iC(fun)  { int ierr=fun; if(ierr!=0) exit(1); }
#define iA(expr) { if((expr)==0) { std::cerr<<"wrong "<<__LINE__<<" in " <<__FILE__<<std::endl; std::cerr.flush(); exit(1); } }

using namespace PEXSI;

// *********************************************************************
// C interface
// *********************************************************************
/// @brief Dummy interface for the purpose of testing the C interface.
extern "C" 
void DummyInterface( MPI_Comm comm, int a )
{
	int mpirank, mpisize;
	MPI_Comm_rank( comm, &mpirank );
	MPI_Comm_size( comm, &mpisize );
	if( mpirank == 0 ){
		std::cout << "Comm rank = " << mpirank << std::endl;
		std::cout << "Comm size = " << mpisize << std::endl;
		std::cout << "Dummy inteface is working and is outputing an integer " 
			<< a << std::endl;
	}
	return;
}  // -----  end of function DummyInterface  ----- 

/// @brief Read the sizes of a DistSparseMatrix for allocating memory in
/// C.
extern "C"
void ReadDistSparseMatrixFormattedHeadInterface (
		char*    filename,
		int*     size,
		int*     nnz,
		int*     nnzLocal,
		int*     numColLocal,
		MPI_Comm comm )
{
  Int mpirank;  MPI_Comm_rank(comm, &mpirank);
  Int mpisize;  MPI_Comm_size(comm, &mpisize);
	std::ifstream fin;
	if( mpirank == 0 ){
		fin.open(filename);
		if( !fin.good() ){
			throw std::logic_error( "File cannot be openeded!" );
		}
		Int dummy;
		fin >> *size >> dummy;
		fin >> *nnz;
	}
	
	MPI_Bcast( &(*size), 1, MPI_INT, 0, comm);
	MPI_Bcast( &(*nnz),  1, MPI_INT, 0, comm);

	IntNumVec  colptr(*size+1);
	if( mpirank == 0 ){
		Int* ptr = colptr.Data();
		for( Int i = 0; i < *size+1; i++ )
			fin >> *(ptr++);
	}

	MPI_Bcast(colptr.Data(), *size+1, MPI_INT, 0, comm);

	// Compute the number of columns on each processor
	IntNumVec numColLocalVec(mpisize);
	Int numColFirst;
	numColFirst = *size / mpisize;
  SetValue( numColLocalVec, numColFirst );
  numColLocalVec[mpisize-1] = *size - numColFirst * (mpisize-1);  
	// Modify the last entry	

	*numColLocal = numColLocalVec[mpirank];

	*nnzLocal = colptr[mpirank * numColFirst + (*numColLocal)] - 
		colptr[mpirank * numColFirst];
	
	// Close the file
	if( mpirank == 0 ){
    fin.close();
	}

	return;
}  
// -----  end of function ReadDistSparseMatrixFormattedHeadInterface


/// @brief Actual reading the data of a DistSparseMatrix, assuming that
/// the arrays have been allocated outside this subroutine.
extern "C"
void ReadDistSparseMatrixFormattedInterface(
		char*     filename,
		int       size,
		int       nnz,
		int       nnzLocal,
		int       numColLocal,
		int*      colptrLocal,
		int*      rowindLocal,
		double*   nzvalLocal,
		MPI_Comm  comm )
{
	DistSparseMatrix<Real> A;
	ReadDistSparseMatrixFormatted( filename, A, comm );
	iA( size == A.size );
	iA( nnz  == A.nnz  );
	iA( nnzLocal == A.nnzLocal );
	iA( numColLocal + 1 == A.colptrLocal.m() );
	
	memcpy( colptrLocal,  A.colptrLocal.Data(),
			sizeof(int) * (numColLocal + 1) );

	memcpy( rowindLocal,  A.colptrLocal.Data(),
			sizeof(int) * (nnzLocal) );

	blas::Copy( nnzLocal, A.nzvalLocal.Data(), 1,
			nzvalLocal, 1 );

	return;
}  
// -----  end of function ReadDistSparseMatrixFormattedInterface  

/// @brief Main interface between PPEXSI and C.
extern "C" 
void PPEXSIInterface (
		int           nrows,
	  int           nnz,	
		int           nnzLocal,
		int           numColLocal,
		int*          colptrLocal,
		int*          rowindLocal,
		double*       HnzvalLocal,
		double*       SnzvalLocal,
		double*      DMnzvalLocal,
		double*     EDMnzvalLocal,
		double*     FDMnzvalLocal,
		int           numPole,
		double        temperature,
		double        numElectronExact,
		double*       numElectron, 
		double        gap,
		double        deltaE,
		double*       mu,
		double        muMin,
		double        muMax,
		int           muMaxIter,
		double        poleTolerance,
		double        numElectronTolerance,
	  MPI_Comm	    comm,
	  int           npPerPole	)
{
	Int mpirank, mpisize;
	MPI_Comm_rank( comm, &mpirank );
	MPI_Comm_size( comm, &mpisize );

	if( mpisize != npPerPole ){
		throw std::runtime_error( "The current interface only supports mpisize == npPerPole" );
	}
	// log files
	std::stringstream  ss;
	ss << "logPEXSI" << mpirank;
	// append to previous log files
	statusOFS.open( ss.str().c_str(), std::ios_base::app );

	Int nprow = iround( std::sqrt( (double)npPerPole) );
	Int npcol = npPerPole / nprow;

	Grid gridPole( comm, mpisize / npPerPole, npPerPole );
	PPEXSIData pexsi( &gridPole, nprow, npcol );

	// Convert into H and S matrices
	// TODO Simultaneous treatment of multiple poles
	DistSparseMatrix<Real> HMat, SMat;

	HMat.size        = nrows;
	HMat.nnzLocal    = nnzLocal;
	HMat.nnz         = nnz;
	// Pointer without actual copying the data
	HMat.colptrLocal = IntNumVec( numColLocal+1, false, colptrLocal );
	HMat.rowindLocal = IntNumVec( nnzLocal,      false, rowindLocal );
	HMat.comm        = gridPole.rowComm;

	// Copy symbolic data from H to S, including communicator
	CopyPattern( HMat, SMat );                    
	
	// Pointer without actual copying the data
	HMat.nzvalLocal  = DblNumVec( nnzLocal, false, HnzvalLocal );
	SMat.nzvalLocal  = DblNumVec( nnzLocal, false, SnzvalLocal );
  	

	// Parameters
	bool isFreeEnergyDensityMatrix = true;
	bool isEnergyDensityMatrix     = true;
	bool isDerivativeTMatrix       = false;
	bool isConverged               = false; 	
	std::string colPerm            = "PARMETIS";
	// Initial guess of chemical potential
	Real mu0                       = *mu;

	std::vector<Real>  muList;
	std::vector<Real>  numElectronList;
	std::vector<Real>  numElectronDrvList;

	Real timeSolveSta, timeSolveEnd;

	GetTime( timeSolveSta );
	pexsi.Solve( 
			numPole,
			temperature,
			numElectronExact,
			gap,
			deltaE,
			mu0,
			muMin,
			muMax,
			HMat,
			SMat,
			muMaxIter,
			poleTolerance,
			numElectronTolerance,
			colPerm,
			isFreeEnergyDensityMatrix,
			isEnergyDensityMatrix,
			isDerivativeTMatrix,
			muList,
			numElectronList,
			numElectronDrvList,
			isConverged	);
	GetTime( timeSolveEnd );

	PrintBlock( statusOFS, "Solve finished." );
	if( isConverged ){
		statusOFS << "PEXSI has converged with " << muList.size() << 
			" iterations" << std::endl;
	}
	else {
		statusOFS << "PEXSI did not converge with " << muList.size() << 
			" iterations" << std::endl;
	}
	
	// Update chemical potential
	*mu = *(muList.rbegin());
	*numElectron = *(numElectronList.rbegin());

	// Update the density matrices
	blas::Copy( nnzLocal, pexsi.DensityMatrix().nzvalLocal.Data(), 
			1, DMnzvalLocal, 1 );

	blas::Copy( nnzLocal, pexsi.FreeEnergyDensityMatrix().nzvalLocal.Data(), 
			1, FDMnzvalLocal, 1 );

	blas::Copy( nnzLocal, pexsi.EnergyDensityMatrix().nzvalLocal.Data(), 
			1, EDMnzvalLocal, 1 );

	Print( statusOFS, "Total time for PEXSI = ", 
			timeSolveEnd - timeSolveSta );

	statusOFS.close();
	
	return;
}  // -----  end of function PPEXSIInterface ----- 

// *********************************************************************
// FORTRAN interface
// 
// All FORTRAN interfaces call the C-interface subroutines above for
// actual computation.  
//
// NOTE: 
//
// 1. The FORTRAN communicators are converted to C communicators using 
// f2c_comm.
//
// 2. All string arguments should have the "correct input" from FORTRAN
// that is consistent with C-convention. This can be done by either
// using iso_c_binding, or using trim(string)//char(0) instead of string
// in FORTRAN.
//
// *********************************************************************
/// @brief Internal subroutine to convert FORTRAN communicator to C
extern "C" 
MPI_Comm f2c_comm(int *Fcomm)
{
	return MPI_Comm_f2c((MPI_Fint)(*Fcomm));
}  // -----  end of function f2c_comm ----- 



/// @brief Dummy interface for the purpose of testing the FORTRAN
/// interface.
extern "C" 
void FORTRAN(f_dummy_interface)( int* Fcomm, int* a ){
	DummyInterface( f2c_comm(Fcomm), *a );
	return;
}  // -----  end of function f_dummy_interface  ----- 

/// @brief Read the sizes of a DistSparseMatrix for allocating memory in
/// FORTRAN.
extern "C"
void FORTRAN(f_read_distsparsematrix_formatted_head) (
		char*    filename,
		int*     size,
		int*     nnz,
		int*     nnzLocal,
		int*     numColLocal,
		int*     Fcomm )
{
  ReadDistSparseMatrixFormattedHeadInterface(
			filename,
			size,
			nnz,
			nnzLocal,
			numColLocal,
			f2c_comm( Fcomm ) );

	return;
}  // -----  end of function f_read_distsparsematrix_formatted_head  


/// @brief Actual reading the data of a DistSparseMatrix, assuming that
/// the arrays have been allocated in FORTRAN.
extern "C"
void FORTRAN(f_read_distsparsematrix_formatted) (
		char*    filename,
		int*     size,
		int*     nnz,
		int*     nnzLocal,
		int*     numColLocal,
		int*     colptrLocal,
		int*     rowindLocal,
		double*  nzvalLocal,
		int*     Fcomm )
{
	ReadDistSparseMatrixFormattedInterface(
			filename,
			*size,
			*nnz,
			*nnzLocal,
			*numColLocal,
			colptrLocal,
			rowindLocal,
			nzvalLocal,
			f2c_comm( Fcomm ) );
	return;
} // -----  end of function f_read_distsparsematrix_formatted  ----- 

/// @brief Main interface between PPEXSI and FORTRAN.
extern "C" 
void FORTRAN(f_ppexsi_interface)(
		int*          nrows, 
		int*          nnz, 
		int*          nnzLocal,
		int*          numColLocal,
		int*          colptrLocal,
		int*          rowindLocal,
		double*       HnzvalLocal,
		double*       SnzvalLocal,
		double*      DMnzvalLocal,
		double*     EDMnzvalLocal,
		double*     FDMnzvalLocal,
		int*          numPole,
		double*       temperature,
		double*       numElectronExact,
		double*       numElectron,
		double*       gap,
		double*       deltaE,
		double*       mu,
		double*       muMin,
		double*       muMax,
		int*          muMaxIter,
		double*       poleTolerance,
		double*       numElectronTolerance,
		int*    	    Fcomm,
		int*          npPerPole	)
{
	PPEXSIInterface( 
			*nrows,
			*nnz,
			*nnzLocal,
			*numColLocal,
			colptrLocal,
			rowindLocal,
			HnzvalLocal,
			SnzvalLocal,
			DMnzvalLocal,
			EDMnzvalLocal,
			FDMnzvalLocal,
			*numPole,
			*temperature,
			*numElectronExact,
			numElectron,
			*gap,
			*deltaE,
			mu,
			*muMin,
			*muMax,
			*muMaxIter,
			*poleTolerance,
			*numElectronTolerance,
			f2c_comm(Fcomm), 
			*npPerPole );

	return;
} // -----  end of function f_ppexsi_interface  ----- 
