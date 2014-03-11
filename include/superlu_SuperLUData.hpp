#ifndef _SUPERLU_DIST_SUPERLUDATA_HPP_
#define _SUPERLU_DIST_SUPERLUDATA_HPP_

#include "environment.hpp"

#include "superlu_dist_interf.hpp"
#include "pselinv.hpp"

namespace PEXSI{

struct SuperLUData{
	/// @brief SuperLU matrix. 
	SuperMatrix         A;                        

	/// @brief SuperLU options. 
	///
	/// Note
	/// ----
	///
	/// It is important to have 
	///
	/// options.RowPerm           = NOROWPERM;
	/// 
	/// to make sure that symmetric permutation is used.
	///
	superlu_options_t   options;                  

	/// @brief Saves the permutation vectors.  Only perm_c (permutation of
	/// column as well as rows due to the symmetric permutation) will be used.
	ScalePermstruct_t   ScalePermstruct;          

	/// @brief SuperLU grid structure.
	gridinfo_t*         grid;

	/// @brief Saves the supernodal partition as well as the numerical
	/// values and structures of the L and U structure.
	LUstruct_t          LUstruct;

	/// @brief Used for solve for multivectors.
	SOLVEstruct_t       SOLVEstruct;

	/// @brief SuperLU statistics
	SuperLUStat_t       stat;

	/// @brief Number of processors used for parallel symbolic
	/// factorization and PARMETIS/PT-SCOTCH
	Int                 numProcSymbFact;

	/// @brief SuperLU information
	Int                 info;

	// The following are for consistency checks
	bool                isSuperMatrixAllocated;
	bool                isSuperMatrixFactorized;
	bool                isScalePermstructAllocated;
	bool                isLUstructAllocated;

  Int maxDomains;

};


}

#endif
