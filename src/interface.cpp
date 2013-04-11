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

// FIXME
#define _DEBUGlevel_ 0

using namespace PEXSI;

// FIXME
extern "C" { 
	double seekeig_(double *, int *, double *, double *, double *);
}


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
	
	blas::Copy( numColLocal+1, A.colptrLocal.Data(), 1,
			colptrLocal, 1 );

	blas::Copy( nnzLocal, A.rowindLocal.Data(), 1,
			rowindLocal, 1 );

	blas::Copy( nnzLocal, A.nzvalLocal.Data(), 1,
			nzvalLocal, 1 );

	return;
}  
// -----  end of function ReadDistSparseMatrixFormattedInterface  

/// @brief Interface between PPEXSI and C for computing the cumulative
/// density of states using inertia counts.
///
/// TODO: Finite temperature dependence of the mu result
extern "C"
void PPEXSIInertiaCountInterface(
		// Input parameters
		int           nrows,                        // Size of the matrix
		int           nnz,                          // Total number of nonzeros in H
		int           nnzLocal,                     // Number of nonzeros in H on this proc
		int           numColLocal,                  // Number of local columns for H
		int*          colptrLocal,                  // Colomn pointer in CSC format
		int*          rowindLocal,                  // Row index pointer in CSC format
		double*       HnzvalLocal,                  // Nonzero value of H in CSC format
		double*       SnzvalLocal,                  // Nonzero falue of S in CSC format
		double        temperature,                  // Temperature, in the same unit as H
		double        numElectronExact,             // Exact number of electrons
		double        muMin0,                       // Initial guess of lower bound for mu
		double        muMax0,                       // Initial guess of upper bound for mu
		int           numPole,                      // Number of shifts in computing the inertia, still called "Pole" for legacy reason
		int           maxIter,                      // Maximum number of iterations for computing the inertia
		double        numElectronTolerance,         // Stopping criterion of inertia count
		int           ordering,                     // SuperLUDIST ordering
		int           npPerPole,                    // Number of processors for each shift, still called "Pole" for legacy reason
		MPI_Comm	    comm,                         // Overall MPI communicator
		// Output parameters
		double*       muMinInertia,                 // Lower bound for mu after inertia count
		double*       muMaxInertia,                 // Upper bound for mu after inertia count
		double*       muLowerEdge,                  // Ne(muLowerEdge) = Ne - eps. For band gapped system
		double*       muUpperEdge,                  // Ne(muUpperEdge) = Ne + eps. For band gapped system
		int*          numIter,                      // Number of actual iterations for inertia count
		double*       muList,                       // The list of shifts
		double*       numElectronList               // The number of electrons corresponding to shifts (0K)
		)
{
	Int mpirank, mpisize;
	MPI_Comm_rank( comm, &mpirank );
	MPI_Comm_size( comm, &mpisize );
	Real timeSta, timeEnd;

	if( mpisize % npPerPole != 0 ){
		std::ostringstream msg;
		msg 
			<< "mpisize    = " << mpisize << std::endl
			<< "npPerPole = " << npPerPole << std::endl
			<< "mpisize is not divisible by npPerPole!" << std::endl;
		throw std::runtime_error( msg.str().c_str() );
	}

	if( numElectronTolerance < 4 ){
		std::ostringstream msg;
		msg 
			<< "numElectronTolerance = " << numElectronTolerance 
			<< ", which is less than 4. " <<
		 	"This is probably too tight for the purpose of inertia count." 
			<< std::endl;
		throw std::runtime_error( msg.str().c_str() );
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
	DistSparseMatrix<Real> HMat, SMat;

	// The first row processors (size: npPerPole) read the matrix, and
	// then distribute the matrix to the rest of the processors.
	//
	// NOTE: The first row processor must have data for H/S.
	std::vector<char> sstr;
	Int sizeStm;
	if( MYROW( &gridPole ) == 0 ){
		std::stringstream sstm;

		HMat.size        = nrows;
		HMat.nnz         = nnz;
		HMat.nnzLocal    = nnzLocal;
		// The first row processor does not need extra copies of the index /
		// value of the matrix. 
		HMat.colptrLocal = IntNumVec( numColLocal+1, false, colptrLocal );
		HMat.rowindLocal = IntNumVec( nnzLocal,      false, rowindLocal );
		// H value
		HMat.nzvalLocal  = DblNumVec( nnzLocal,      false, HnzvalLocal );
		HMat.comm = gridPole.rowComm;

		CopyPattern( HMat, SMat );
		SMat.comm = gridPole.rowComm;

		// S value
		SMat.nzvalLocal  = DblNumVec( nnzLocal,      false, SnzvalLocal );
		
		// Serialization will copy the values regardless of the ownership
		serialize( HMat, sstm, NO_MASK );
		serialize( SMat.nzvalLocal, sstm, NO_MASK );
		
		sstr.resize( Size( sstm ) );
		sstm.read( &sstr[0], sstr.size() ); 	
		sizeStm = sstr.size();
	}
	
	MPI_Bcast( &sizeStm, 1, MPI_INT, 0, gridPole.colComm );
	
#if ( _DEBUGlevel_ >= 0 )
	statusOFS << "sizeStm = " << sizeStm << std::endl;
#endif

	if( MYROW( &gridPole ) != 0 ) sstr.resize( sizeStm );

	MPI_Bcast( (void*)&sstr[0], sizeStm, MPI_BYTE, 0, gridPole.colComm );

	if( MYROW( &gridPole ) != 0 ){
		std::stringstream sstm;
		sstm.write( &sstr[0], sizeStm );
		deserialize( HMat, sstm, NO_MASK );
		// Communicator
		HMat.comm = gridPole.rowComm;
		CopyPattern( HMat, SMat );
		SMat.comm = gridPole.rowComm;
		deserialize( SMat.nzvalLocal, sstm, NO_MASK );
	}
	sstr.clear();


#if ( _DEBUGlevel_ >= 0 )
	statusOFS << "H.size     = " << HMat.size     << std::endl;
	statusOFS << "H.nnzLocal = " << HMat.nnzLocal << std::endl;
	statusOFS << "S.size     = " << SMat.size     << std::endl;
	statusOFS << "S.nnzLocal = " << SMat.nnzLocal << std::endl;
#endif


	// Parameters
	std::string colPerm;
	switch (ordering){
		case 0:
			colPerm = "PARMETIS";
			break;
		case 1:
			colPerm = "METIS_AT_PLUS_A";
			break;
		case 2:
			colPerm = "MMD_AT_PLUS_A";
			break;
		default:
			throw std::logic_error("Unsupported ordering strategy.");
	}

	// In inertia counts, "Pole" is the same as "Shifts" (since there is no actual pole here)
	Int  numShift = numPole;
	std::vector<Real>  shiftVec( numShift );
	std::vector<Real>  inertiaVec( numShift );

	Real muMin = muMin0;
	Real muMax = muMax0;
	const Real EPS = 1e-3;                        // For numerical stability

	bool isConverged = false;
	Int  iter = 0;
	for( iter = 0; iter < maxIter; iter++ ){
		for( Int l = 0; l < numShift; l++ ){
			shiftVec[l] = muMin + l * (muMax - muMin) / (numShift-1);
		}

		GetTime( timeSta );
		pexsi.CalculateNegativeInertia( 
				shiftVec,
				inertiaVec,
				HMat,
				SMat,
				colPerm );

		GetTime( timeEnd );

		for( Int l = 0; l < numShift; l++ ){
			// Inertia is multiplied by 2.0 to reflect the doubly occupied
			// orbitals.
			inertiaVec[l] *= 2.0;

#if ( _DEBUGlevel_ >= 0 )
			statusOFS << std::setiosflags(std::ios::left) 
				<< std::setw(LENGTH_VAR_NAME) << "Shift = "
				<< std::setw(LENGTH_VAR_DATA) << shiftVec[l]
				<< std::setw(LENGTH_VAR_NAME) << "Inertia = "
				<< std::setw(LENGTH_VAR_DATA) << inertiaVec[l]
				<< std::endl << std::endl;
#endif
		}

#if ( _DEBUGlevel_ >= 0 )
		statusOFS << std::endl << "Time for iteration " 
			<< iter << " of the inertia count is " 
			<< timeEnd - timeSta << std::endl;
#endif

		if( inertiaVec[0] >= numElectronExact ||
				inertiaVec[numShift-1] <= numElectronExact ){
			std::ostringstream msg;
			msg 
				<< "The solution is not in the prescribed (muMin, muMax) interval." << std::endl
				<< "(muMin, muMax) ~ (" << shiftVec[0] << " , " << shiftVec[numShift-1] << " ) " << std::endl
				<< "(Ne(muMin), Ne(muMax)) ~ (" << inertiaVec[0] << " , " << inertiaVec[numShift-1] 
				<< " ) " << std::endl
				<< "NeExact = " << numElectronExact << std::endl;
			throw std::runtime_error( msg.str().c_str() );
		}

		if( inertiaVec[numShift-1] - inertiaVec[0] < numElectronTolerance ){
			isConverged = true;
			break;
		}

		// Update muMin, muMax
		std::vector<Real>::iterator vi0, vi1;
		vi0 = std::lower_bound( inertiaVec.begin(), inertiaVec.end(), 
				numElectronExact-2-EPS );
		vi1 = std::lower_bound( inertiaVec.begin(), inertiaVec.end(), 
				numElectronExact+2-EPS );

		Int idx0 = vi0 - inertiaVec.begin();
		Int idx1 = vi1 - inertiaVec.begin();

		muMin = shiftVec[idx0];
		muMax = shiftVec[idx1];
	} // for (iter)

	if( isConverged ){
		statusOFS << std::endl << "Inertia count converged. N(muMax) - N(muMin) = " <<
			inertiaVec[numShift-1] - inertiaVec[0] << std::endl;
	}
	else {
		statusOFS << std::endl << "Inertia count did not converge. N(muMax) - N(muMin) = " <<
			inertiaVec[numShift-1] - inertiaVec[0] << std::endl;
	}

	// TODO Add the finite temperature effect
	Real muLow, muUpp;
	{
		double mu0, netarget;
		mu0 = ( muMin + muMax ) / 2.0;
		netarget = numElectronExact - EPS;
		muLow = seekeig_(&netarget, &numShift, &shiftVec[0],&inertiaVec[0], &mu0);
		mu0 = ( muMin + muMax ) / 2.0;
		netarget = numElectronExact + EPS;
		muUpp = seekeig_(&netarget, &numShift, &shiftVec[0],&inertiaVec[0], &mu0);
	}


	statusOFS << std::endl << "After the inertia count," << std::endl;
	Print( statusOFS, "numIter = ", iter );
	Print( statusOFS, "muLowerEdge   = ", muLow );
	Print( statusOFS, "muUppperEdge  = ", muUpp );
	Print( statusOFS, "muMin         = ", muMin );
	Print( statusOFS, "muMax         = ", muMax );
	statusOFS << std::endl;

	statusOFS.close();

	// Convert the internal variables to output parameters
	*muMinInertia = muMin;
	*muMaxInertia = muMax;
	*muLowerEdge  = muLow;
	*muUpperEdge  = muUpp;
	*numIter      = iter;

	for( Int l = 0; l < numShift; l++ ){
		muList[l]          = shiftVec[l];
		numElectronList[l] = inertiaVec[l];
	}

	return;
}  // -----  end of function PPEXSIInertiaCountInterface ----- 


/// @brief Interface between PPEXSI and C for the PPEXSI solve
/// procedure.
extern "C" 
void PPEXSISolveInterface (
		// Input parameters
		int           nrows,                        // Size of the matrix
	  int           nnz,                          // Total number of nonzeros in H
		int           nnzLocal,                     // Number of nonzeros in H on this proc
		int           numColLocal,                  // Number of local columns for H
		int*          colptrLocal,                  // Colomn pointer in CSC format
		int*          rowindLocal,                  // Row index pointer in CSC format
		double*       HnzvalLocal,                  // Nonzero value of H in CSC format
		double*       SnzvalLocal,                  // Nonzero falue of S in CSC format
		double        temperature,                  // Temperature, in the same unit as H
		double        numElectronExact,             // Exact number of electrons
		double        mu0,                          // Initial guess for mu
		double        muMin0,                       // Initial guess for lower bound of mu
		double        muMax0,                       // Initial guess for upper bound of mu
		double        gap,                          // Energy gap (lower bound)
		double        deltaE,                       // Spectral radius of S^{-1}H
		int           numPole,                      // Number of poles
		int           maxIter,                      // Maximum number of iterations for mu-iteration in PEXSI
		double        numElectronTolerance,         // Stopping criterion of PEXSI mu iteration.
		double        poleTolerance,                // Tolerance for neglecting a pole
		int           ordering,                     // SuperLUDIST ordering
	  int           npPerPole,                    // Number of processors for each pole
	  MPI_Comm	    comm,                         // Overall MPI communicator
		// Output parameters
		double*      DMnzvalLocal,                  // Nonzero value of density matrix in CSC format
		double*     EDMnzvalLocal,                  // Nonzero value of energy density matrix in CSC format
		double*     FDMnzvalLocal,                  // Nonzero value of free energy density matrix in CSC format
		double*       muPEXSI,                      // Final chemical potential
		double*       numElectronPEXSI,             // Computed number of electron at the final chemical potential
		double*       muMinPEXSI,                   // Final lower bound for mu.
		double*       muMaxPEXSI,                   // Final upper bound for mu
		int*          numIter,                      // Number of actual iterations for PEXSI
		double*       muList,                       // The history of mu
		double*       numElectronList,              // The history of number of electrons correspondig to mu
		double*       numElectronDrvList            // The history of dN/dMu
		)
{
	Int mpirank, mpisize;
	MPI_Comm_rank( comm, &mpirank );
	MPI_Comm_size( comm, &mpisize );
	Real timeSta, timeEnd;

	if( mpisize % npPerPole != 0 ){
		std::ostringstream msg;
		msg 
			<< "mpisize    = " << mpisize << std::endl
			<< "npPerPole  = " << npPerPole << std::endl
			<< "mpisize is not divisible by npPerPole!" << std::endl;
		throw std::runtime_error( msg.str().c_str() );
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
	DistSparseMatrix<Real> HMat, SMat;

	// The first row processors (size: npPerPole) read the matrix, and
	// then distribute the matrix to the rest of the processors.
	//
	// NOTE: The first row processor must have data for H/S.
	std::vector<char> sstr;
	Int sizeStm;
	if( MYROW( &gridPole ) == 0 ){
		std::stringstream sstm;

		HMat.size        = nrows;
		HMat.nnz         = nnz;
		HMat.nnzLocal    = nnzLocal;
		// The first row processor does not need extra copies of the index /
		// value of the matrix. 
		HMat.colptrLocal = IntNumVec( numColLocal+1, false, colptrLocal );
		HMat.rowindLocal = IntNumVec( nnzLocal,      false, rowindLocal );
		// H value
		HMat.nzvalLocal  = DblNumVec( nnzLocal,      false, HnzvalLocal );
		HMat.comm = gridPole.rowComm;

		CopyPattern( HMat, SMat );
		SMat.comm = gridPole.rowComm;

		// S value
		SMat.nzvalLocal  = DblNumVec( nnzLocal,      false, SnzvalLocal );
		
		// Serialization will copy the values regardless of the ownership
		serialize( HMat, sstm, NO_MASK );
		serialize( SMat.nzvalLocal, sstm, NO_MASK );
		
		sstr.resize( Size( sstm ) );
		sstm.read( &sstr[0], sstr.size() ); 	
		sizeStm = sstr.size();
	}
	
	MPI_Bcast( &sizeStm, 1, MPI_INT, 0, gridPole.colComm );
	
#if ( _DEBUGlevel_ >= 0 )
	statusOFS << "sizeStm = " << sizeStm << std::endl;
#endif

	if( MYROW( &gridPole ) != 0 ) sstr.resize( sizeStm );

	MPI_Bcast( (void*)&sstr[0], sizeStm, MPI_BYTE, 0, gridPole.colComm );

	if( MYROW( &gridPole ) != 0 ){
		std::stringstream sstm;
		sstm.write( &sstr[0], sizeStm );
		deserialize( HMat, sstm, NO_MASK );
		// Communicator
		HMat.comm = gridPole.rowComm;
		CopyPattern( HMat, SMat );
		SMat.comm = gridPole.rowComm;
		deserialize( SMat.nzvalLocal, sstm, NO_MASK );
	}
	sstr.clear();


#if ( _DEBUGlevel_ >= 0 )
	statusOFS << "H.size     = " << HMat.size     << std::endl;
	statusOFS << "H.nnzLocal = " << HMat.nnzLocal << std::endl;
	statusOFS << "S.size     = " << SMat.size     << std::endl;
	statusOFS << "S.nnzLocal = " << SMat.nnzLocal << std::endl;
#endif


	// Parameters
	bool isFreeEnergyDensityMatrix = true;
	bool isEnergyDensityMatrix     = true;
	bool isDerivativeTMatrix       = false;
	bool isConverged               = false; 	
	std::string colPerm;
	switch (ordering){
		case 0:
			colPerm = "PARMETIS";
			break;
		case 1:
			colPerm = "METIS_AT_PLUS_A";
			break;
		case 2:
			colPerm = "MMD_AT_PLUS_A";
			break;
		default:
			throw std::logic_error("Unsupported ordering strategy.");
	}

	std::vector<Real>  muVec;
	std::vector<Real>  numElectronVec;
	std::vector<Real>  numElectronDrvVec;

	Real timeSolveSta, timeSolveEnd;

	Real muMin = muMin0; 
	Real muMax = muMax0;
	Real mu    = mu0;

	GetTime( timeSolveSta );
	pexsi.Solve( 
			numPole,
			temperature,
			numElectronExact,
			gap,
			deltaE,
			mu,
			muMin,
			muMax,
			HMat,
			SMat,
			maxIter,
			poleTolerance,
			numElectronTolerance,
			colPerm,
			isFreeEnergyDensityMatrix,
			isEnergyDensityMatrix,
			isDerivativeTMatrix,
			muVec,
			numElectronVec,
			numElectronDrvVec,
			isConverged	);
	GetTime( timeSolveEnd );

	Int muIter = muVec.size();

	PrintBlock( statusOFS, "Solve finished." );
	if( isConverged ){
		statusOFS << "PEXSI has converged with " << muIter << 
			" iterations" << std::endl;
	}
	else {
		statusOFS << "PEXSI did not converge with " << muIter << 
			" iterations" << std::endl;
	}



	// Convert the internal variables to output parameters

	// Update the density matrices for the first row processors (size: npPerPole) 
	if( MYROW( &gridPole ) == 0 ){
		blas::Copy( nnzLocal, pexsi.DensityMatrix().nzvalLocal.Data(), 
				1, DMnzvalLocal, 1 );

		blas::Copy( nnzLocal, pexsi.FreeEnergyDensityMatrix().nzvalLocal.Data(), 
				1, FDMnzvalLocal, 1 );

		blas::Copy( nnzLocal, pexsi.EnergyDensityMatrix().nzvalLocal.Data(), 
				1, EDMnzvalLocal, 1 );
	}

	*muPEXSI          = *(muVec.rbegin());
	*numElectronPEXSI = *(numElectronVec.rbegin());
	*muMinPEXSI       = muMin;
	*muMaxPEXSI       = muMax;
	*numIter          = muIter;
	for( Int i = 0; i < muIter; i++ ){
		muList[i]               = muVec[i];
		numElectronList[i]      = numElectronVec[i];
		numElectronDrvList[i]   = numElectronDrvVec[i];
	}


	MPI_Barrier( comm );

	// Compute the guess for T=0 chemical potential. Not used anymore
//	*muZeroT = pexsi.EstimateZeroTemperatureChemicalPotential(
//			temperature,
//			*mu,
//			SMat );
//
//	Print( statusOFS, "guess of mu(T=0) = ", 
//			*muZeroT );

	Print( statusOFS, "Total time for PEXSI = ", 
			timeSolveEnd - timeSolveSta );

	statusOFS.close();
	
	return;
}  // -----  end of function PPEXSISolveInterface ----- 

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

extern "C" 
void FORTRAN(f_ppexsi_inertiacount_interface)(
		// Input parameters
		int*          nrows,                        // Size of the matrix
		int*          nnz,                          // Total number of nonzeros in H
		int*          nnzLocal,                     // Number of nonzeros in H on this proc
		int*          numColLocal,                  // Number of local columns for H
		int*          colptrLocal,                  // Colomn pointer in CSC format
		int*          rowindLocal,                  // Row index pointer in CSC format
		double*       HnzvalLocal,                  // Nonzero value of H in CSC format
		double*       SnzvalLocal,                  // Nonzero falue of S in CSC format
		double*       temperature,                  // Temperature, in the same unit as H
		double*       numElectronExact,             // Exact number of electrons
		double*       muMin0,                       // Initial guess of lower bound for mu
		double*       muMax0,                       // Initial guess of upper bound for mu
		int*          numPole,                      // Number of shifts in computing the inertia, still called "Pole" for legacy reason
		int*          maxIter,                      // Maximum number of iterations for computing the inertia
		double*       numElectronTolerance,         // Stopping criterion of inertia count
		int*          ordering,                     // SuperLUDIST ordering
		int*          npPerPole,                    // Number of processors for each shift, still called "Pole" for legacy reason
		int*    	    Fcomm,                        // Overall MPI communicator
		// Output parameters
		double*       muMinInertia,                 // Lower bound for mu after inertia count
		double*       muMaxInertia,                 // Upper bound for mu after inertia count
		double*       muLowerEdge,                  // Ne(muLowerEdge) = Ne - eps. For band gapped system
		double*       muUpperEdge,                  // Ne(muUpperEdge) = Ne + eps. For band gapped system
		int*          numIter,                      // Number of actual iterations for inertia count
		double*       muList,                       // The list of shifts
		double*       numElectronList               // The number of electrons corresponding to shifts (0K)
		)
{
	PPEXSIInertiaCountInterface( 
			*nrows,
			*nnz,
			*nnzLocal,
			*numColLocal,
			colptrLocal,
			rowindLocal,
			HnzvalLocal,
			SnzvalLocal,
			*temperature,
			*numElectronExact,
			*muMin0,
			*muMax0,
			*numPole,
			*maxIter,
			*numElectronTolerance,
			*ordering,
			*npPerPole,
			f2c_comm(Fcomm),
			muMinInertia,
			muMaxInertia,
			muLowerEdge,
			muUpperEdge,
			numIter,
			muList,
			numElectronList );

	return;
} // -----  end of function f_ppexsi_inertiacount_interface  ----- 


/// @brief Interface between PPEXSI and FORTRAN for the solve procedure.
extern "C" 
void FORTRAN(f_ppexsi_solve_interface)(
		// Input parameters
		int*          nrows,                        // Size of the matrix
	  int*          nnz,                          // Total number of nonzeros in H
		int*          nnzLocal,                     // Number of nonzeros in H on this proc
		int*          numColLocal,                  // Number of local columns for H
		int*          colptrLocal,                  // Colomn pointer in CSC format
		int*          rowindLocal,                  // Row index pointer in CSC format
		double*       HnzvalLocal,                  // Nonzero value of H in CSC format
		double*       SnzvalLocal,                  // Nonzero falue of S in CSC format
		double*       temperature,                  // Temperature, in the same unit as H
		double*       numElectronExact,             // Exact number of electrons
		double*       mu0,                          // Initial guess for mu
		double*       muMin0,                       // Initial guess for lower bound of mu
		double*       muMax0,                       // Initial guess for upper bound of mu
		double*       gap,                          // Energy gap (lower bound)
		double*       deltaE,                       // Spectral radius of S^{-1}H
		int*          numPole,                      // Number of poles
		int*          maxIter,                      // Maximum number of iterations for mu-iteration in PEXSI
		double*       numElectronTolerance,         // Stopping criterion of PEXSI mu iteration.
		double*       poleTolerance,                // Tolerance for neglecting a pole
		int*          ordering,                     // SuperLUDIST ordering
	  int*          npPerPole,                    // Number of processors for each pole
		int*    	    Fcomm,                        // Overall MPI communicator
		// Output parameters
		double*      DMnzvalLocal,                  // Nonzero value of density matrix in CSC format
		double*     EDMnzvalLocal,                  // Nonzero value of energy density matrix in CSC format
		double*     FDMnzvalLocal,                  // Nonzero value of free energy density matrix in CSC format
		double*       muPEXSI,                      // Final chemical potential
		double*       numElectronPEXSI,             // Computed number of electron at the final chemical potential
		double*       muMinPEXSI,                   // Final lower bound for mu.
		double*       muMaxPEXSI,                   // Final upper bound for mu
		int*          numIter,                      // Number of actual iterations for PEXSI
		double*       muList,                       // The history of mu
		double*       numElectronList,              // The history of number of electrons correspondig to mu
		double*       numElectronDrvList            // The history of dN/dMu
		)
{
	PPEXSISolveInterface( 
			*nrows,
			*nnz,
			*nnzLocal,
			*numColLocal,
			colptrLocal,
			rowindLocal,
			HnzvalLocal,
			SnzvalLocal,
			*temperature,
			*numElectronExact,
			*mu0,
			*muMin0,
			*muMax0,
			*gap,
			*deltaE,
			*numPole,
			*maxIter,
			*numElectronTolerance,
			*poleTolerance,
		  *ordering,
			*npPerPole,
			f2c_comm(Fcomm),
			DMnzvalLocal,
			EDMnzvalLocal,
			FDMnzvalLocal,
			muPEXSI,
			numElectronPEXSI,
			muMinPEXSI,
			muMaxPEXSI,
			numIter,
		  muList,
	  	numElectronList,
		  numElectronDrvList);

	return;
} // -----  end of function f_ppexsi_solve_interface  ----- 
