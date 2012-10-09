#ifndef _ENVIRONMENT_DECL_HPP_
#define _ENVIRONMENT_DECL_HPP_

// STL libraries
#include <iostream> 
#include <iomanip> 
#include <fstream>
#include <sstream>
#include <unistd.h>

#include <cfloat>
#include <complex>
#include <string>

#include <set>
#include <map>
#include <stack>
#include <vector>

#include <algorithm>
#include <cmath>

#include <cassert>
#include <stdexcept>

// MPI
#include "mpi.h"

/***********************************************************************
 *  Data types and constants
 **********************************************************************/

namespace PEXSI{

// Basic data types

#define BLAS(name) name##_
#define LAPACK(name) name##_

typedef    int                   Int;
typedef    double                Real;
typedef    std::complex<double>  Complex; // Must use elemental form of complex
#ifdef _USE_COMPLEX_
typedef    std::complex<double>  Scalar;  // Must use elemental form of complex
#else
typedef    double                Scalar;
#endif

// IO
extern  std::ofstream  statusOFS;



// *********************************************************************
// Define constants
// *********************************************************************
// Commonly used
const Int I_ZERO = 0;
const Int I_ONE = 1;
const Real D_ZERO = 0.0;
const Real D_ONE  = 1.0;
const Complex Z_ZERO = Complex(0.0, 0.0);
const Complex Z_ONE  = Complex(1.0, 0.0);
const Scalar SCALAR_ZERO    = static_cast<Scalar>(0.0);
const Scalar SCALAR_ONE     = static_cast<Scalar>(1.0);
const char UPPER = 'U';
const char LOWER = 'L';

// Physical constants

const Real au2K = 315774.67;
const Real PI = 3.141592653589793;

} // namespace PEXSI

/***********************************************************************
 *  Error handling
 **********************************************************************/

namespace PEXSI{


#ifndef _RELEASE_
void PushCallStack( std::string s );
void PopCallStack();
void DumpCallStack();
#endif // ifndef _RELEASE_

// We define an output stream that does nothing. This is done so that the 
// root process can be used to print data to a file's ostream while all other 
// processes use a null ostream. 
struct NullStream : std::ostream
{            
	struct NullStreamBuffer : std::streambuf
	{
		Int overflow( Int c ) { return traits_type::not_eof(c); }
	} nullStreamBuffer_;

	NullStream() 
		: std::ios(&nullStreamBuffer_), std::ostream(&nullStreamBuffer_)
		{ }
};  

// *********************************************************************
// Global utility functions 
// These utility functions do not depend on local definitions
// *********************************************************************
// Return the closest integer to a rela number
Int iround( Real a );

// Read the options from command line
Int OptionsCreate(Int argc, char** argv, 
		std::map<std::string,std::string>& options);

} // namespace PEXSI


#endif // _ENVIRONMENT_DECL_HPP_
