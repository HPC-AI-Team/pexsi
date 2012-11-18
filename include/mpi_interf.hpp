#ifndef _MPI_INTERF_HPP_
#define _MPI_INTERF_HPP_

#include  "environment_impl.hpp"

namespace PEXSI{

namespace mpi{


// *********************************************************************
// Allgatherv
// *********************************************************************
void Allgatherv( 
		std::vector<Int>& localVec, 
		std::vector<Int>& allVec,
		MPI_Comm          comm );

} // namespace mpi

} // namespace PEXSI

#endif // _MPI_INTERF_HPP_

